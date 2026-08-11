#' Ensure a person identifier column exists for longitudinal LME analyses
#'
#' @param data Data frame containing the longitudinal phenotype-plus-beta data.
#' @param personVar Character. Name of the subject identifier column.
#' @param sidVar Character. Name of the fallback sample identifier column used
#' to
#'   derive `personVar` when it is missing.
#'
#' @return A list containing the updated data frame, whether `personVar` was
#'   created, an optional SID-to-person preview, and counts per person.
#'
#' @description
#' Internal helper that derives a subject identifier from `SID` when the
#' requested `personVar` column is not present.
#'
#' @keywords internal
#' @noRd
ensurePersonColumnMethylationLME <- function(data, personVar = "person",
    sidVar = "SID") {
    person_created <- FALSE
    mapping_preview <- NULL
    if (!(personVar %in% colnames(data))) {
        if (!(sidVar %in% colnames(data))) {
            stop("Column '", personVar,
            "' was not found and cannot be created because '",
                sidVar, "' is missing.", call. = FALSE)
        }
        sid_values <- validateSampleIdentifiersDnaEpico(data[[sidVar]],
            paste0("inputPheno$", sidVar))
        if (any(!grepl("[AB]$", sid_values))) {
            stop("Cannot safely derive '", personVar, "' from '",
                sidVar, "': every SID must end in A or B. ",
            "Supply an explicit subject identifier column.",
                call. = FALSE)
        }
        person_values <- sub("[AB]$", "", sid_values)
        visit_values <- sub("^.*([AB])$", "\\1", sid_values)
        if (any(!nzchar(person_values))) {
            stop("Derived subject identifiers cannot be empty.",
                call. = FALSE)
        }
        if (anyDuplicated(paste(person_values, visit_values, sep = "\r"))) {
            stop("SID contains duplicate subject/visit identifiers.",
                call. = FALSE)
        }
        data[[personVar]] <- person_values
        person_created <- TRUE
        mapping_preview <- utils::head(data[order(data[[personVar]],
            data[[sidVar]]), c(sidVar, personVar), drop = FALSE],
            20L) }
    person_values <- as.character(data[[personVar]])
    if (anyNA(person_values) || any(!nzchar(trimws(person_values)))) {
        stop("personVar contains missing or blank subject identifiers.",
            call. = FALSE) }
    person_counts <- table(data[[personVar]], useNA = "ifany")
    if (length(person_counts) < 2L) {
        stop("At least two subjects are required for mixed-effects modeling.",
            call. = FALSE)
    }
    if (max(person_counts) < 2L) {
        stop("At least one subject must have repeated observations for ",
            "mixed-effects modeling.", call. = FALSE)
    }
    list(data = data, personCreated = person_created, mappingPreview =
        mapping_preview,
        personCounts = person_counts)
}

#' Build a mixed-effects formula for methylationLME helpers
#'
#' @param phenotype Character. Phenotype variable of interest.
#' @param personVar Character. Subject identifier variable used as a random
#'   intercept.
#' @param covariates Character vector of additional fixed-effect covariates.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return Character scalar containing a model formula.
#'
#' @description
#' Internal helper that builds the per-CpG linear mixed-effects formula used by
#' `methylationLME()` and its modular helpers.
#'
#' @keywords internal
#' @noRd
buildFormulaMethylationLME <- function(phenotype,
    personVar, covariates = character(0),
    interactionTerm = NULL, includeRandomTerm = TRUE,
    responseVar = "beta") {
    if (!is.null(interactionTerm) && nzchar(interactionTerm) &&
        identical(interactionTerm, phenotype)) {
        stop("interactionTerm must differ from the phenotype being modelled.",
            call. = FALSE)
    }
    quoted_phenotype <- quoteNamesMethylationGLM(phenotype)
    quoted_person <- quoteNamesMethylationGLM(personVar)
    fixed_terms <- unique(c(covariates))
    if (!is.null(interactionTerm) && nzchar(interactionTerm)) {
        quoted_interaction <- quoteNamesMethylationGLM(interactionTerm)
        interaction_part <- paste(quoted_phenotype,
            quoted_interaction, sep = " * ")
        fixed_terms <- setdiff(fixed_terms,
            c(interactionTerm, phenotype))
        fixed_formula_terms <- c(interaction_part,
            quoteNamesMethylationGLM(fixed_terms))
    }
    else {
        fixed_terms <- setdiff(fixed_terms,
            phenotype)
        fixed_formula_terms <- c(quoted_phenotype,
            quoteNamesMethylationGLM(fixed_terms))
    }
    fixed_formula_terms <- unique(fixed_formula_terms[nzchar(
        fixed_formula_terms)])
    if (length(fixed_formula_terms) == 0L) {
        stop("At least one fixed-effect term is required.",
            call. = FALSE)
    }
    formula_text <- paste(responseVar, "~",
        paste(fixed_formula_terms, collapse = " + "))
    if (isTRUE(includeRandomTerm)) {
        formula_text <- paste(formula_text,
            "+ (1 |", quoted_person, ")")
    }
    formula_text
}

normalizeCorrelationStructureMethylationLME <-
    function(correlationStructure = "none") {
    if (length(correlationStructure) == 0L ||
        all(is.na(correlationStructure))) {
        return("none")
    }
    if (length(correlationStructure) != 1L) {
        stop("correlationStructure must contain one value.",
        call. = FALSE
        )
    }

    structure_value <- tolower(trimws(as.character(correlationStructure[[1L]])))
    if (structure_value %in% c("", "none", "null", "na")) {
        return("none")
    }
    if (identical(structure_value, "ar1")) {
        return("AR1")
    }
    if (identical(structure_value, "car1")) {
        return("CAR1")
    }

    stop("correlationStructure must be one of: none, AR1, CAR1.",
        call. = FALSE
    )
    }

normalizeCorrelationVariableMethylationLME <- function(correlationVar = NULL) {
    if (is.null(correlationVar) || length(correlationVar) ==
    0L || all(is.na(correlationVar))) {
    return(NULL)
    }

    parsed <- splitOptionMinfiEwasWater(correlationVar, sep = ",")
    if (length(parsed) == 0L) {
    return(NULL)
    }
    if (length(parsed) > 1L) {
    stop("correlationVar must contain a single variable name.",
        call. = FALSE
    )
    }

    parsed <- trimws(as.character(parsed[[1L]]))
    if (!nzchar(parsed) || tolower(parsed) %in% c("null", "na")) {
    return(NULL)
    }

    parsed
}

resolveLmeLibrariesMethylationLME <- function(lmeLibs = "lme4,lmerTest") {
    requested <- splitOptionMinfiEwasWater(lmeLibs, sep = ",")
    if (length(requested) == 0L) {
    requested <- c("lme4", "lmerTest")
    }

    requested_lower <- tolower(requested)
    has_nlme <- "nlme" %in% requested_lower
    has_lme4 <- any(requested_lower %in% c("lme4", "lmertest"))

    if (isTRUE(has_nlme) && isTRUE(has_lme4)) {
    stop("lmeLibs must choose either 'lme4,lmerTest' or 'nlme', not both.",
        call. = FALSE
    )
    }

    if (isTRUE(has_nlme)) {
    return(list(
        engine = "nlme", requestedPackages = requested,
        requiredPackages = "nlme"
    ))
    }

    if (isTRUE(has_lme4)) {
    return(list(
        engine = "lme4", requestedPackages = requested,
        requiredPackages = unique(c(requested, "lme4", "lmerTest"))
    ))
    }

    stop("lmeLibs must contain either 'lme4'/'lmerTest' or 'nlme'.",
    call. = FALSE
    )
}

normalizeOmnibusDdfMethylationLME <- function(
    omnibusDdf = "Satterthwaite"
) {
    if (!is.character(omnibusDdf) || length(omnibusDdf) != 1L ||
    is.na(omnibusDdf)) {
    stop(
        "omnibusDdf must be either 'Satterthwaite' or 'Kenward-Roger'.",
        call. = FALSE
    )
    }

    normalized <- tolower(trimws(omnibusDdf))
    if (identical(normalized, "satterthwaite")) {
    return("Satterthwaite")
    }
    if (normalized %in% c("kenward-roger", "kenwardroger")) {
    return("Kenward-Roger")
    }

    stop(
    "omnibusDdf must be either 'Satterthwaite' or 'Kenward-Roger'.",
    call. = FALSE
    )
}

validateOmnibusConfigurationMethylationLME <- function(
    omnibusTest = FALSE, omnibusDdf = "Satterthwaite",
    lmeEngine = "lme4"
) {
    omnibus_test <- validateLogicalScalarDnaEpico(
    omnibusTest,
    "omnibusTest"
    )
    omnibus_ddf <- normalizeOmnibusDdfMethylationLME(omnibusDdf)

    if (isTRUE(omnibus_test) && identical(lmeEngine, "nlme")) {
    stop(
        "omnibusTest is currently available only for the lmerTest/lme4 ",
        "engine, not nlme.",
        call. = FALSE
    )
    }
    if (isTRUE(omnibus_test) &&
    identical(omnibus_ddf, "Kenward-Roger") &&
    !requireNamespace("pbkrtest", quietly = TRUE)) {
    stop(
        "omnibusDdf = 'Kenward-Roger' requires the suggested package ",
        "'pbkrtest'.",
        call. = FALSE
    )
    }

    list(test = omnibus_test, ddf = omnibus_ddf)
}

resolveOmnibusTargetTermMethylationLME <- function(
    formulaText, data, phenotype, interactionTerm = NULL
) {
    fixed_formula <- removeRandomInterceptMethylationModels(formulaText)
    terms_object <- stats::terms(
    stats::as.formula(fixed_formula),
    data = data
    )
    term_labels <- attr(stats::delete.response(terms_object), "term.labels")
    target_variables <- if (!is.null(interactionTerm) &&
    nzchar(interactionTerm)) {
    c(phenotype, interactionTerm)
    } else {
    phenotype
    }

    matches_target <- vapply(term_labels, function(term_label) {
    term_variables <- all.vars(stats::as.formula(paste("~", term_label)))
    length(term_variables) == length(target_variables) &&
        setequal(term_variables, target_variables)
    }, logical(1))
    matched_terms <- term_labels[matches_target]
    if (length(matched_terms) != 1L) {
    target_label <- if (length(target_variables) == 1L) {
        target_variables
    } else {
        paste(target_variables, collapse = ":")
    }
    target_message <- paste(
        "Could not identify one fixed-effect model term",
        "for omnibus testing"
    )
    stop(
        sprintf(
        "%s: %s.",
        target_message, target_label
        ),
        call. = FALSE
    )
    }

    unname(matched_terms[[1L]])
}

emptyOmnibusResultMethylationLME <- function(
    term, method, status = "not_estimable", reason = NA_character_
) {
    list(
    term = term, method = method,
    fValue = NA_real_, numeratorDf = NA_real_,
    denominatorDf = NA_real_, pValue = NA_real_,
    status = status, reason = reason, modelMessage = NA_character_,
    rhs = 0
    )
}

computeOmnibusTestMethylationLME <- function(fit,
    coefficientTerms, omnibusTerm, omnibusDdf = "Satterthwaite") {
    result <- emptyOmnibusResultMethylationLME(term = omnibusTerm,
        method = omnibusDdf)
    captured <- captureModelConditionsDnaEpico({
        fixed_effects <- lme4::fixef(fit)
        fixed_names <- names(fixed_effects)
        mapped_terms <- coefficientTerms[fixed_names]
        selected <- which(!is.na(mapped_terms) &
            mapped_terms == omnibusTerm)
        if (length(selected) == 0L) {
            stop("The omnibus model term has no estimable fixed-effect ",
                "coefficients.", call. = FALSE)
        }
        contrast <- matrix(0, nrow = length(selected),
            ncol = length(fixed_effects), dimnames = list(fixed_names[selected],
                fixed_names))
        contrast[cbind(seq_along(selected), selected)] <- 1
        test <- lmerTest::contestMD(model = fit,
            L = contrast, rhs = 0, ddf = omnibusDdf,
            joint = TRUE)
        required <- c("F value", "NumDF", "DenDF",
            "Pr(>F)")
        if (!all(required %in% colnames(test))) {
            stop("The omnibus test did not return the expected F-test ",
                "statistics.", call. = FALSE)
        }
        values <- unlist(test[1L, required, drop = TRUE],
            use.names = FALSE)
        if (any(!is.finite(values))) {
            stop("The omnibus test returned missing or non-finite statistics.",
                call. = FALSE)
        }
        list(fValue = unname(test[["F value"]][[1L]]),
            numeratorDf = unname(test[["NumDF"]][[1L]]),
            denominatorDf = unname(test[["DenDF"]][[1L]]),
            pValue = unname(test[["Pr(>F)"]][[1L]]))
    })
    result$modelMessage <- captured$modelMessage
    if (!is.null(captured$error)) {
        result$reason <- conditionMessage(captured$error)
        result
    }
    else {
        result[names(captured$value)] <- captured$value
        result$status <- "tested"
        result$reason <- NA_character_
        result
    }
}

collectOmnibusTestsMethylationLME <- function(fits, phenotype) {
    rows <- lapply(names(fits), function(cpg) {
    fit <- fits[[cpg]]
    if (is.null(fit) || inherits(
        fit,
        "dnaEPICO_methylationLME_fit_error"
    ) || is.null(fit$omnibus)) {
        return(NULL)
    }
    omnibus <- fit$omnibus
    data.frame(
        Phenotype = phenotype, CpG = cpg,
        Omnibus.Term = as.character(omnibus$term),
        Omnibus.F.Value = as.numeric(omnibus$fValue),
        Omnibus.Num.DF = as.numeric(omnibus$numeratorDf),
        Omnibus.Den.DF = as.numeric(omnibus$denominatorDf),
        Omnibus.P.Value = as.numeric(omnibus$pValue),
        Omnibus.Method = as.character(omnibus$method),
        Omnibus.Status = as.character(omnibus$status),
        Omnibus.Reason = as.character(omnibus$reason),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) == 0L) {
    return(data.frame(
        Phenotype = character(0), CpG = character(0),
        Omnibus.Term = character(0), Omnibus.F.Value = numeric(0),
        Omnibus.Num.DF = numeric(0), Omnibus.Den.DF = numeric(0),
        Omnibus.P.Value = numeric(0), Omnibus.Method = character(0),
        Omnibus.Status = character(0), Omnibus.Reason = character(0),
        stringsAsFactors = FALSE, check.names = FALSE
    ))
    }
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

coerceCorrelationTimeMethylationLME <- function(x) {
    if (is.numeric(x)) {
    return(as.numeric(x))
    }

    character_x <- trimws(as.character(x))
    non_missing <- !is.na(x)
    numeric_pattern <- "^[+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?$"
    is_numeric_text <- grepl(numeric_pattern, character_x)

    if (all(is_numeric_text[non_missing])) {
    numeric_x <- rep(NA_real_, length(x))
    numeric_x[non_missing] <- as.numeric(character_x[non_missing])
    return(numeric_x)
    }

    stop(
    "correlationVar must be numeric or contain numeric text; ",
    "categorical values cannot define AR1/CAR1 spacing.",
    call. = FALSE
    )
}

addCorrelationTimeVariableMethylationLME <- function(
    modelData,
    correlationVar
) {
    correlation_time_var <- "dnaEPICO_lme_correlation_time"
    while (correlation_time_var %in% colnames(modelData)) {
    correlation_time_var <- paste0(
        correlation_time_var,
        "_"
    )
    }

    modelData[[correlation_time_var]] <-
    coerceCorrelationTimeMethylationLME(modelData[[correlationVar]])

    list(data = modelData, correlationTimeVar = correlation_time_var)
}

validateCorrelationTimeMethylationLME <- function(
    modelData,
    correlationTimeVar, personVar, correlationStructure
) {
    time_values <- modelData[[correlationTimeVar]]
    person_values <- modelData[[personVar]]

    if (anyNA(time_values) || any(!is.finite(time_values))) {
    stop(
        "correlationVar contains missing or non-finite values after ",
        "conversion.",
        call. = FALSE
    )
    }
    person_time <- interaction(person_values, time_values, drop = TRUE)
    if (anyDuplicated(person_time)) {
    stop(
        "correlationVar must be unique within each subject for AR1/CAR1 ",
        "models.",
        call. = FALSE
    )
    }
    if (identical(correlationStructure, "AR1") && any(abs(time_values -
    round(time_values)) > sqrt(.Machine$double.eps))) {
    stop("AR1 requires integer-valued correlationVar observations.",
        call. = FALSE
    )
    }

    invisible(TRUE)
}

buildCorrelationMethylationLME <- function(
    correlationStructure = "none",
    correlationTimeVar = NULL, personVar
) {
    if (identical(correlationStructure, "none")) {
    return(NULL)
    }

    if (is.null(correlationTimeVar) || !nzchar(correlationTimeVar)) {
    stop(
        "An internal correlation time variable is required for AR1/CAR1 ",
        "models.",
        call. = FALSE
    )
    }

    correlation_formula <- stats::as.formula(paste(
    "~", quoteNamesMethylationGLM(correlationTimeVar),
    "|", quoteNamesMethylationGLM(personVar)
    ))

    switch(correlationStructure,
    AR1 = nlme::corAR1(form = correlation_formula),
    CAR1 = nlme::corCAR1(form = correlation_formula),
    stop("Unsupported correlation structure.",
        call. = FALSE
    )
    )
}

coerceCoefficientTableMethylationLME <- function(coefTable, engine = "lme4") {
    coef_table <- as.data.frame(coefTable, check.names = FALSE)

    if (identical(engine, "nlme")) {
    rename_map <- c(
        Value = "Estimate", Std.Error = "Std. Error",
        DF = "df", `t-value` = "t value", `p-value` = "Pr(>|t|)"
    )
    for (old_name in names(rename_map)) {
        if (old_name %in% colnames(coef_table)) {
        colnames(coef_table)[colnames(coef_table) ==
            old_name] <- rename_map[[old_name]]
        }
    }
    }

    as.matrix(coef_table)
}

#' Find coefficient rows for a longitudinal phenotype effect or interaction
#'
#' @param coefNames Character vector of coefficient names.
#' @param phenotype Character. Phenotype variable of interest.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return Character vector of matching coefficient names.
#'
#' @description
#' Internal helper that identifies the longitudinal fixed-effect terms that
#' should be extracted from each fitted mixed-effects model.
#'
#' @keywords internal
#' @noRd
findCoefficientRowsMethylationLME <- function(
    coefNames, phenotype,
    interactionTerm = NULL, coefficientTerms = NULL
) {
    findCoefficientRowsMethylationGLM(
    coefNames = coefNames,
    variable = phenotype, interactionTerm = interactionTerm,
    coefficientTerms = coefficientTerms
    )
}

collectModelMessagesMethylationLME <- function(fits) {
    rows <- list()
    row_index <- 1L
    for (phenotype in names(fits)) {
    fit_group <- fits[[phenotype]]
    for (cpg in names(fit_group)) {
        fit_object <- fit_group[[cpg]]
        rows[[row_index]] <- data.frame(
        Phenotype = phenotype,
        CpG = cpg,
        Model.Message = modelMessageDnaEpico(fit_object),
        P.Value.Available = isTRUE(fit_object$pValueAvailable),
        stringsAsFactors = FALSE, check.names = FALSE
        )
        row_index <- row_index + 1L
    }
    }

    if (length(rows) == 0L) {
    return(data.frame(
        Phenotype = character(0), CpG = character(0),
        Model.Message = character(0), P.Value.Available = logical(0),
        stringsAsFactors = FALSE, check.names = FALSE
    ))
    }

    model_messages <- do.call(rbind, rows)
    rownames(model_messages) <- NULL
    model_messages
}

#' Summarize a single CpG-level mixed-effects model fit
#'
#' @param cpg Character. CpG identifier.
#' @param modelObj List returned by `fitCpGModelMethylationLME()`.
#' @param phenotype Character. Phenotype variable of interest.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return A data frame for the requested CpG, or `NULL`.
#'
#' @description
#' Internal helper that extracts phenotype-specific fixed-effect rows from a
#' single CpG-level mixed-effects model fit.
#'
#' @keywords internal
#' @noRd
summarizeCpGFitMethylationLME <- function(
    cpg, modelObj, phenotype,
    interactionTerm = NULL
) {
    if (is.null(modelObj) || inherits(
    modelObj,
    "dnaEPICO_methylationLME_fit_error"
    )) {
    return(NULL)
    }

    coef_table <- modelObj$coef
    if (is.null(coef_table)) {
    return(NULL)
    }

    matched_terms <- findCoefficientRowsMethylationLME(
    coefNames = rownames(coef_table),
    phenotype = phenotype, interactionTerm = interactionTerm,
    coefficientTerms = modelObj$coefficientTerms
    )
    if (length(matched_terms) == 0L) {
    return(NULL)
    }

    do.call(rbind, lapply(matched_terms, function(term) {
    coef_row <- coef_table[term, ]
    data.frame(
        CpG = cpg, Interaction.Term = term,
        Estimate = unname(coef_row["Estimate"]),
        Std.Error = unname(coef_row["Std. Error"]),
        t.value = unname(coef_row["t value"]),
        P.value = unname(coef_row["Pr(>|t|)"]),
        Model.Message = modelMessageDnaEpico(modelObj),
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    }))
}

filterSummaryByPvalueMethylationLME <- function(summaryDf, pValueFilter) {
    summary_df <- summaryDf
    if (is.null(summary_df) || nrow(summary_df) == 0L) {
    return(data.frame())
    }

    if (nrow(summary_df) > 0L && !is.na(pValueFilter)) {
    keep <- is.finite(summary_df$P.value) & summary_df$P.value <
        pValueFilter
    summary_df <- summary_df[keep, , drop = FALSE]
    }
    rownames(summary_df) <- NULL

    summary_df
}

fitLmeEngineDnaEpico <- function(
    formula, data, personVar, engine, correlationStructure,
    correlationTimeVar
) {
    if (identical(engine, "nlme")) {
    fit <- nlme::lme(
        fixed = formula,
        random = stats::as.formula(paste(
        "~ 1 |", quoteNamesMethylationGLM(personVar)
        )),
        correlation = buildCorrelationMethylationLME(
        correlationStructure = correlationStructure,
        correlationTimeVar = correlationTimeVar,
        personVar = personVar
        ),
        data = data, na.action = stats::na.exclude, method = "REML",
        control = nlme::lmeControl(returnObject = TRUE)
    )
    coefficients <- coerceCoefficientTableMethylationLME(
        summary(fit)$tTable,
        engine = engine
    )
    } else {
    fit <- lmerTest::lmer(
        formula = formula, data = data,
        na.action = stats::na.exclude, REML = TRUE
    )
    coefficients <- coerceCoefficientTableMethylationLME(
        summary(fit)$coefficients,
        engine = engine
    )
    }
    list(fit = fit, coefficients = coefficients)
}

coefficientTermsForLmeDnaEpico <- function(
    supplied, formulaText, data, engine
) {
    if (!is.null(supplied)) {
    return(supplied)
    }
    buildCoefficientTermMapMethylationModels(
    formulaText = formulaText, data = data,
    removeRandomEffects = !identical(engine, "nlme")
    )
}

omnibusForLmeFitDnaEpico <- function(
    fit, coefficientTerms, omnibusTest, omnibusTerm, omnibusDdf, engine
) {
    if (!isTRUE(omnibusTest)) {
    return(NULL)
    }
    if (identical(engine, "nlme")) {
    stop("omnibusTest is not available for nlme fits.", call. = FALSE)
    }
    computeOmnibusTestMethylationLME(
    fit = fit, coefficientTerms = coefficientTerms,
    omnibusTerm = omnibusTerm, omnibusDdf = omnibusDdf
    )
}

retainLmeFitDetailsDnaEpico <- function(result, fit, engine, retainModel) {
    if (!isTRUE(retainModel)) {
    return(result)
    }
    result$residuals <- stats::residuals(fit)
    result$fitted <- stats::fitted(fit)
    result$ranef <- if (identical(engine, "nlme")) {
    nlme::ranef(fit)
    } else {
    lme4::ranef(fit)
    }
    result$fixef <- if (identical(engine, "nlme")) {
    nlme::fixef(fit)
    } else {
    lme4::fixef(fit)
    }
    result$model <- fit
    result
}

completeCapturedLmeFitDnaEpico <- function(captured) {
    if (!is.null(captured$error)) {
    return(newMethylationFitErrorDnaEpico(
        reason = conditionMessage(captured$error),
        errorClass = "dnaEPICO_methylationLME_fit_error",
        modelMessage = captured$modelMessage
    ))
    }
    captured$value$modelMessage <- combineModelMessagesDnaEpico(
    captured$modelMessage,
    if (is.null(captured$value$omnibus)) {
        NULL
    } else {
        captured$value$omnibus$modelMessage
    }
    )
    captured$value
}

fitCpGModelMethylationLME <- function(
    cpg, cpgValues, modelData,
    formulaText, personVar, lmeEngine = "lme4", correlationStructure = "none",
    correlationTimeVar = NULL, responseVar = "beta",
    omnibusTest = FALSE, omnibusDdf = "Satterthwaite",
    omnibusTerm = NULL, retainModel = TRUE,
    formulaObject = NULL, coefficientTerms = NULL
) {
    captured <- captureModelConditionsDnaEpico({
    if (!is.numeric(cpgValues)) {
        stop("The CpG response is not numeric.", call. = FALSE)
    }
    model_data <- modelData
    model_data[[responseVar]] <- cpgValues
    fit_formula <- if (is.null(formulaObject)) {
        stats::as.formula(formulaText)
    } else {
        formulaObject
    }
    engine_fit <- fitLmeEngineDnaEpico(
        fit_formula, model_data, personVar, lmeEngine,
        correlationStructure, correlationTimeVar
    )
    coefficient_terms <- coefficientTermsForLmeDnaEpico(
        coefficientTerms, formulaText, model_data, lmeEngine
    )
    omnibus <- omnibusForLmeFitDnaEpico(
        engine_fit$fit, coefficient_terms, omnibusTest,
        omnibusTerm, omnibusDdf, lmeEngine
    )
    result <- list(
        coef = engine_fit$coefficients,
        coefficientTerms = coefficient_terms,
        omnibus = omnibus, pValueAvailable = FALSE
    )
    retainLmeFitDetailsDnaEpico(
        result, engine_fit$fit, lmeEngine, retainModel
    )
    })
    completeCapturedLmeFitDnaEpico(captured)
}

fitMethylationLMEBatch <- function(cpgBatch, data, modelData,
    formulaText, personVar, lmeEngine = "lme4", correlationStructure = "none",
    correlationTimeVar = NULL, phenotype, interactionTerm = NULL,
    responseVar = "beta", omnibusTest = FALSE, omnibusDdf = "Satterthwaite",
    omnibusTerm = NULL, formulaObject = NULL, coefficientTerms = NULL) {
    fits <- vector("list", length(cpgBatch))
    names(fits) <- cpgBatch
    summaries <- vector("list", length(cpgBatch))
    names(summaries) <- cpgBatch
    for (cpg in cpgBatch) {
        model_obj <- fitCpGModelMethylationLME(cpg = cpg,
            cpgValues = cpgResponseMethylationModels(data,
                cpg), modelData = modelData, formulaText = formulaText,
            personVar = personVar, lmeEngine = lmeEngine,
            correlationStructure = correlationStructure,
            correlationTimeVar = correlationTimeVar, responseVar = responseVar,
            omnibusTest = omnibusTest, omnibusDdf = omnibusDdf,
            omnibusTerm = omnibusTerm, retainModel = FALSE,
            formulaObject = formulaObject, coefficientTerms = coefficientTerms)
        summary_row <- summarizeCpGFitMethylationLME(cpg = cpg,
            modelObj = model_obj, phenotype = phenotype,
            interactionTerm = interactionTerm)
        coefficient_p_available <- !is.null(summary_row) &&
            "P.value" %in% colnames(summary_row) && any(is.finite(
            summary_row$P.value))
        omnibus_p_available <- !is.null(model_obj$omnibus) &&
            is.finite(model_obj$omnibus$pValue)
        p_value_available <- coefficient_p_available || omnibus_p_available
        model_obj$pValueAvailable <- p_value_available
        if (!is.null(summary_row)) {
            summary_row$Model.Message <- modelMessageDnaEpico(model_obj)
        }; fits[[cpg]] <- model_obj
        summaries[[cpg]] <- summary_row
    }
    summaries <- Filter(Negate(is.null), summaries)
    summary_df <- if (length(summaries) == 0L) {
        data.frame() } else {
        out <- do.call(rbind, summaries)
        rownames(out) <- NULL; out }
    list(coefficientResults = compactCoefficientResultsMethylationModels(fits =
        fits,
        cpgOrder = cpgBatch, includeResidualSD = FALSE),
        summaries = summary_df, omnibusTests =
            collectOmnibusTestsMethylationLME(fits = fits,
            phenotype = phenotype), modelMessages =
            collectBatchModelMessagesMethylationModels(fits,
            phenotype), fitFailures = collectBatchFitFailuresMethylationModels(
            fits,
            phenotype, "dnaEPICO_methylationLME_fit_error"))
}

summarizeTimepointGroupDnaEpico <- function(data, timeVar, phenotypes) {
    output <- list()
    output[[timeVar]] <- as.character(data[[timeVar]][[1L]])
    for (phenotype in phenotypes) {
    values <- data[[phenotype]]
    if (is.numeric(values)) {
        finite_values <- values[is.finite(values)]
        output[[paste0(phenotype, "_mean")]] <-
        if (length(finite_values)) mean(finite_values) else NA_real_
        output[[paste0(phenotype, "_sd")]] <-
        if (length(finite_values) > 1L) {
            stats::sd(finite_values)
        } else {
            NA_real_
        }
        output[[paste0(phenotype, "_n")]] <- length(finite_values)
    } else {
        observed <- unique(as.character(values[!is.na(values)]))
        output[[paste0(phenotype, "_n")]] <- sum(!is.na(values))
        output[[paste0(phenotype, "_levels")]] <- paste(
        observed,
        collapse = ","
        )
    }
    }
    as.data.frame(output, stringsAsFactors = FALSE)
}

#' Summarize phenotype values by timepoint for longitudinal methylation analyses
#'
#' @param data Data frame containing the longitudinal phenotype-plus-beta data.
#' @param timeVar Character. Name of the time variable.
#' @param phenotypes Character vector or comma-separated string of phenotype
#'   variables to summarize.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A data frame with one row per timepoint and summary columns for each
#'   requested phenotype.
#'
#' @description
#' Summarize the requested phenotype variables by timepoint. Numeric phenotypes
#' are reported with mean, standard deviation, and non-missing counts;
#' non-numeric
#' phenotypes are reported with non-missing counts and the observed levels.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' timepoint_summary <- summarizeTimepointsMethylationLME(
#'   data = ex$preparedData$data,
#'   timeVar = "Timepoint",
#'   phenotypes = "score",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' nrow(timepoint_summary)
#'
#' @export
summarizeTimepointsMethylationLME <- function(
    data, timeVar = "Timepoint",
    phenotypes, verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir,
    log_file = log_file
    )
    phenotype_list <- unique(splitOptionMinfiEwasWater(phenotypes,
    sep = ","
    ))

    if (!(timeVar %in% colnames(data))) {
    stop("timeVar column not found: ", timeVar, call. = FALSE)
    }

    missing_phenotypes <- setdiff(phenotype_list, colnames(data))
    if (length(missing_phenotypes) > 0L) {
    missing_phenotypes_text <- paste(missing_phenotypes, collapse = ", ")
    stop(
        sprintf(
        "Phenotype columns not found for timepoint summary: %s",
        missing_phenotypes_text
        ),
        call. = FALSE
    )
    }

    split_data <- split(data, as.character(data[[timeVar]]))
    summaries <- lapply(split_data, summarizeTimepointGroupDnaEpico,
    timeVar = timeVar, phenotypes = phenotype_list
    )
    summary_df <- do.call(rbind, summaries)
    rownames(summary_df) <- NULL

    emitLogMinfiEwasWater(
    c(
        "============================================================",
        "Summary statistics for phenotype scores by timepoint:",
        previewLinesMinfiEwasWater(summary_df),
        "============================================================"
    ),
    verbose = verbose, log_path = log_path
    )

    summary_df
}

normalizeLmePreparationConfigDnaEpico <- function(config) {
    config$phenotypeList <- unique(splitOptionMinfiEwasWater(
    config$phenotypes,
    sep = ","
    ))
    config$covariateList <- splitOptionMinfiEwasWater(
    config$covariates,
    sep = ","
    )
    config$factorList <- splitOptionMinfiEwasWater(config$factorVars, sep = ",")
    config$scaleList <- normalizeScaleVariablesDnaEpico(config$scaleVars)
    config$prsMapParsed <- parsePrsMapMethylationGLM(config$prsMap)
    config$cpgLimitValue <- validateCpgLimitMethylationModels(
    normalizeOptionalNumericMethylationGLM(config$cpgLimit)
    )
    config$cpgPrefix <- validateCpgPrefixDnaEpico(config$cpgPrefix)
    config$methylationScaleValue <- normalizeMethylationScaleDnaEpico(
    config$methylationScale
    )
    config$responseLabel <- methylationScaleResponseLabelDnaEpico(
    config$methylationScaleValue
    )
    config$objectPrefix <- methylationScaleObjectPrefixDnaEpico(
    config$methylationScaleValue
    )
    config$responseColumn <- methylationScaleResponseColumnDnaEpico(
    config$methylationScaleValue
    )
    config$logPath <- resolveLogPathMinfiEwasWater(
    logs = config$logs, log_dir = config$log_dir,
    log_file = config$log_file
    )
    config
}

loadMethylationLMEInputDnaEpico <- function(config) {
    input_name <- sub("[.][^.]+$", "", basename(config$inputPheno))
    scale_name <- if (startsWith(input_name, config$objectPrefix)) {
    input_name
    } else {
    character(0)
    }
    legacy_name <- if (identical(config$methylationScaleValue, "beta")) {
    "phenoBT1T2"
    } else {
    character(0)
    }
    data <- loadSavedObjectPreprocessingPheno(
    config$inputPheno,
    preferred_name = c(
        scale_name, paste0(config$objectPrefix, "T1T2"),
        legacy_name, input_name
    )
    )
    if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    }
    data
}

validateLmeColumnsDnaEpico <- function(data, columns, messageTemplate) {
    missing <- setdiff(columns, colnames(data))
    if (length(missing) > 0L) {
    missing_text <- paste(unique(missing), collapse = ", ")
    stop(sprintf(messageTemplate, missing_text),
        call. = FALSE
    )
    }
    invisible(TRUE)
}

validateLmeTimeValuesDnaEpico <- function(data, timeVar) {
    if (!(timeVar %in% colnames(data))) {
    stop("timeVar column not found in inputPheno: ", timeVar,
        call. = FALSE
    )
    }
    values <- data[[timeVar]]
    invalid <- is.na(values) | !nzchar(trimws(as.character(values)))
    if (is.numeric(values)) {
    invalid <- invalid | !is.finite(values)
    }
    if (any(invalid)) {
    stop("timeVar contains missing, blank, or non-finite values.",
        call. = FALSE
    )
    }
    invisible(TRUE)
}

validateLmeModelVariablesDnaEpico <- function(data, config) {
    if (length(config$phenotypeList) == 0L) {
    stop("At least one phenotype must be supplied.", call. = FALSE)
    }
    validateLmeTimeValuesDnaEpico(data, config$timeVar)
    validateLmeColumnsDnaEpico(
    data, config$phenotypeList,
    "Phenotype columns not found in inputPheno: %s"
    )
    validateLmeColumnsDnaEpico(
    data, config$covariateList,
    "Covariate columns not found in inputPheno: %s"
    )
    mapped_prs <- unname(config$prsMapParsed[
    names(config$prsMapParsed) %in% config$phenotypeList
    ])
    validateLmeColumnsDnaEpico(
    data, mapped_prs, "PRS columns not found in inputPheno: %s"
    )
    interaction <- normalizeOptionalColumnMethylationModels(
    config$interactionTerm, "interactionTerm"
    )
    if (!is.null(interaction) && !(interaction %in% colnames(data))) {
    stop("interactionTerm column not found in inputPheno: ", interaction,
        call. = FALSE
    )
    }
    validateLmeColumnsDnaEpico(
    data, config$factorList,
    "Factor columns not found in inputPheno: %s"
    )
    list(mappedPrs = mapped_prs, interaction = interaction)
}

prepareLmeColumnsDnaEpico <- function(data, config, variables) {
    for (var in intersect(config$factorList, colnames(data))) {
    data[[var]] <- as.factor(data[[var]])
    }
    cpg_columns <- grep(
    paste0("^", escapeRegexMethylationGLM(config$cpgPrefix)),
    colnames(data),
    value = TRUE
    )
    if (!is.na(config$cpgLimitValue)) {
    cpg_columns <- utils::head(cpg_columns, config$cpgLimitValue)
    }
    if (length(cpg_columns) == 0L) {
    stop("No CpG columns were found with prefix '", config$cpgPrefix,
        "'.",
        call. = FALSE
    )
    }
    validateMethylationProbeIdentifiersDnaEpico(
    cpg_columns, "CpG columns in the LME input"
    )
    scaled_cpgs <- intersect(config$scaleList, cpg_columns)
    if (length(scaled_cpgs) > 0L) {
    scaled_cpg_text <- paste(scaled_cpgs, collapse = ", ")
    stop(sprintf(
        "%s: %s",
        "CpG methylation response columns cannot be listed in scaleVars",
        scaled_cpg_text
    ), call. = FALSE)
    }
    eligible <- unique(c(
    config$phenotypeList, config$covariateList,
    variables$mappedPrs, variables$interaction
    ))
    eligible <- eligible[!is.na(eligible) & nzchar(eligible)]
    scaling <- scaleModelVariablesDnaEpico(
    data = data[, setdiff(colnames(data), cpg_columns), drop = FALSE],
    scaleVars = config$scaleList, factorVars = config$factorList,
    eligibleVars = eligible, protectedVars = c(config$personVar, cpg_columns)
    )
    list(data = data, cpgColumns = cpg_columns, scaling = scaling)
}

logPreparedLmeDataDnaEpico <- function(config,
    data, personData, columns, timepointSummary) {
    lines <- c("============================================================",
        paste("Loaded longitudinal phenotype + methylation data from:",
            config$inputPheno), paste("Merged modeling object:            ",
            config$objectPrefix, "*"), paste(
            "Data dimensions:                  ",
            paste(dim(data), collapse = " x ")),
        paste("Person variable:                  ",
            config$personVar), paste("Time variable:                    ",
            config$timeVar), paste("Phenotypes:                       ",
            paste(config$phenotypeList, collapse = ", ")),
        paste("Covariates:                       ",
            paste(config$covariateList, collapse = ", ")),
        paste("Factor variables:                 ",
            paste(config$factorList, collapse = ", ")),
        formatScalingMetadataLogDnaEpico(columns$scaling$metadata),
        paste("CpG columns retained:             ",
            length(columns$cpgColumns)),
        if (isTRUE(personData$personCreated)) {
            paste("Created person variable from SID: ",
                config$personVar)
        } else {
            paste("Person variable already present:  ",
                config$personVar)
        }, paste("Values observed in", config$timeVar,
            ":"), previewLinesMinfiEwasWater(table(data[[config$timeVar]],
            useNA = "ifany")))
    if (!is.null(personData$mappingPreview)) {
        lines <- c(lines, "Example mapping of SID to person ID:",
            previewLinesMinfiEwasWater(personData$mappingPreview))
    }
    lines <- c(lines, "Summary statistics for phenotype scores by timepoint:",
        previewLinesMinfiEwasWater(timepointSummary),
        "============================================================")
    emitLogMinfiEwasWater(lines, verbose = config$verbose,
        log_path = config$logPath)
}

newPreparedMethylationLMEDataDnaEpico <- function(
    config, columns, variables, personData, timepointSummary
) {
    structure(list(
    data = columns$data, modelData = columns$scaling$data,
    inputPheno = config$inputPheno,
    inputIdentity = inputIdentityMethylationModels(config$inputPheno),
    personVar = config$personVar, timeVar = config$timeVar,
    phenotypes = config$phenotypeList, covariates = config$covariateList,
    factorVars = config$factorList,
    scaleVars = columns$scaling$scaleVars,
    scalingMetadata = columns$scaling$metadata,
    prsMap = config$prsMapParsed, cpgColumns = columns$cpgColumns,
    cpgPrefix = config$cpgPrefix, cpgLimit = config$cpgLimitValue,
    methylationScale = config$methylationScaleValue,
    responseLabel = config$responseLabel,
    methylationObjectPrefix = config$objectPrefix,
    internalResponseColumn = config$responseColumn,
    interactionTerm = variables$interaction,
    requestedInteractionTerm = config$interactionTerm,
    personCreated = personData$personCreated,
    personCounts = personData$personCounts,
    personMappingPreview = personData$mappingPreview,
    timepointSummary = timepointSummary
    ), class = "dnaEPICO_methylationLME_data")
}

#' Prepare longitudinal phenotype-plus-methylation data for mixed-effects
#' analyses
#'
#' @param inputPheno Character. Path to the merged longitudinal phenotype-plus-
#'   methylation object created by `preprocessingPheno()`.
#' @param personVar Character. Name of the subject identifier column.
#' @param timeVar Character. Name of the time variable.
#' @param phenotypes Character vector or comma-separated string of phenotype
#'   variables to model.
#' @param covariates Character vector or comma-separated string of covariate
#'   variables to adjust for.
#' @param factorVars Character vector or comma-separated string of variables
#' that
#'   should be converted to factors before modeling.
#' @param scaleVars Character vector, comma-separated variable names, or `NULL`.
#'   Numeric fixed-effect variables to standardize before model fitting.
#' @param prsMap Character vector or comma-separated string of phenotype-to-PRS
#'   mappings in the form `'Phenotype:PRS'`.
#' @param cpgPrefix Character. Prefix used to identify methylation columns.
#' @param cpgLimit Integer or `NA`. Maximum number of CpGs to retain. `NA`
#'   keeps all matching CpGs.
#' @param methylationScale Character. Methylation metric represented by the CpG
#'   columns. One of `'Beta'`, `'M'`, or `'CN'`, case-insensitively.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationLME_data'` containing
#'   the prepared analysis data, parsed variable selections, CpG columns,
#'   timepoint summaries, and subject-ID diagnostics.
#'
#' @description
#' Load the merged longitudinal phenotype-plus-methylation object, ensure that
#' a subject
#' identifier column is available, validate the requested modeling variables,
#' convert selected variables to factors, and return a single in-memory object
#' for downstream mixed-effects modeling helpers.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' prepared_data <- prepareMethylationLMEData(
#'   inputPheno = ex$inputPath,
#'   personVar = "person",
#'   timeVar = "Timepoint",
#'   phenotypes = "score",
#'   covariates = "sex",
#'   factorVars = "sex",
#'   cpgLimit = 2,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(prepared_data)
#'
#' @export
prepareMethylationLMEData <- function(
    inputPheno, personVar = "person",
    timeVar = "Timepoint", phenotypes, covariates, factorVars,
    scaleVars = NULL, prsMap = NULL, cpgPrefix = "cg", cpgLimit = NA,
    methylationScale = "beta", interactionTerm = NULL, verbose = FALSE,
    logs = FALSE, log_dir = NULL, log_file = "log_methylationLME.txt"
) {
    config <- normalizeLmePreparationConfigDnaEpico(
    as.list(environment(), all.names = TRUE)
    )
    analysis_data <- loadMethylationLMEInputDnaEpico(config)
    if (length(config$phenotypeList) == 0L) {
    stop("At least one phenotype must be supplied.", call. = FALSE)
    }
    person_data <- ensurePersonColumnMethylationLME(
    data = analysis_data, personVar = config$personVar
    )
    analysis_data <- person_data$data
    variables <- validateLmeModelVariablesDnaEpico(analysis_data, config)
    columns <- prepareLmeColumnsDnaEpico(analysis_data, config, variables)
    timepoint_summary <- summarizeTimepointsMethylationLME(
    data = columns$data, timeVar = config$timeVar,
    phenotypes = config$phenotypeList, verbose = FALSE, logs = FALSE
    )
    logPreparedLmeDataDnaEpico(
    config, columns$data, person_data, columns, timepoint_summary
    )
    newPreparedMethylationLMEDataDnaEpico(
    config, columns, variables, person_data, timepoint_summary
    )
}
validateLmeCorrelationConfigDnaEpico <- function(config, preparedData) {
    if (!identical(config$engine, "nlme") &&
    !identical(config$correlationStructureValue, "none")) {
    stop(
        "correlationStructure can only be AR1 or CAR1 when lmeLibs ",
        "selects 'nlme'.",
        call. = FALSE
    )
    }
    if (!identical(config$correlationStructureValue, "none")) {
    if (is.null(config$correlationVarValue)) {
        stop(
        "correlationVar must be supplied when correlationStructure ",
        "is AR1 or CAR1.",
        call. = FALSE
        )
    }
    if (config$correlationVarValue %in% preparedData$scaleVars) {
        stop(
        "correlationVar cannot also be listed in scaleVars; provide ",
        "a separate scaled fixed-effect column if required.",
        call. = FALSE
        )
    }
    }
    config
}

normalizeLmeFitConfigDnaEpico <- function(config, preparedData) {
    config$resumeFromSummary <- validateLogicalScalarDnaEpico(
    config$resumeFromSummary, "resumeFromSummary"
    )
    config$logPath <- resolveLogPathMinfiEwasWater(
    logs = config$logs, log_dir = config$log_dir,
    log_file = config$log_file
    )
    if (is.null(config$libPath)) {
    config$libPath <- .libPaths()
    }
    libraries <- resolveLmeLibrariesMethylationLME(config$lmeLibs)
    config$requestedPackages <- libraries$requestedPackages
    config$requiredPackages <- libraries$requiredPackages
    config$engine <- libraries$engine
    config$omnibus <- validateOmnibusConfigurationMethylationLME(
    config$omnibusTest, config$omnibusDdf, config$engine
    )
    if (config$omnibus$test && identical(config$omnibus$ddf, "Kenward-Roger")) {
    config$requiredPackages <- unique(c(
        config$requiredPackages, "pbkrtest"
    ))
    }
    config$correlationStructureValue <-
    normalizeCorrelationStructureMethylationLME(
        config$correlationStructure
    )
    config$correlationVarValue <- normalizeCorrelationVariableMethylationLME(
    config$correlationVar
    )
    config$nCoresValue <- validatePositiveIntegerMethylationModels(
    config$nCores, "nCores"
    )
    validateLmeCorrelationConfigDnaEpico(config, preparedData)
}

newLmeFitStateDnaEpico <- function(preparedData, config) {
    state <- new.env(parent = emptyenv())
    state$fits <- list()
    state$summaryCache <- list()
    state$coefficientResults <- list()
    state$phenotypeSummaries <- list()
    state$summaryFiles <- list()
    state$modelMessages <- list()
    state$fitFailures <- list()
    state$resumedPhenotypes <- character(0)
    state$fittedPhenotypes <- character(0)
    state$omnibusTests <- list()
    state$formulas <- stats::setNames(
    character(length(preparedData$phenotypes)), preparedData$phenotypes
    )
    state$failureCounts <- stats::setNames(
    integer(length(preparedData$phenotypes)), preparedData$phenotypes
    )
    state$failureReasons <- list()
    state$omnibusTargets <- stats::setNames(
    character(length(preparedData$phenotypes)), preparedData$phenotypes
    )
    state$parallelPlan <- resolveParallelPlanMethylationModels(
    engine = config$engine, nCores = config$nCoresValue,
    nCpGs = length(preparedData$cpgColumns),
    analysisData = preparedData$data,
    modelData = if (is.null(preparedData$modelData)) {
        preparedData$data
    } else {
        preparedData$modelData
    }
    )
    state$parallelPlan$pilotMemoryMB <- NA_real_
    state$backend <- state$parallelPlan$backend
    state$workerCount <- state$parallelPlan$workerCount
    state$batches <- chunkCpGColumnsMethylationModels(
    preparedData$cpgColumns, state$workerCount,
    batchesPerCore = 8L
    )
    state$cluster <- NULL
    state$clusterUseCount <- 0L
    state$pilotCompleted <- FALSE
    state
}

stopLmeFitClusterDnaEpico <- function(state) {
    if (!is.null(state$cluster)) {
    try(parallel::stopCluster(state$cluster), silent = TRUE)
    state$cluster <- NULL
    }
    invisible(NULL)
}

lmePhenotypeVariablesDnaEpico <- function(
    preparedData, modelData, phenotype, config
) {
    prs <- if (phenotype %in% names(preparedData$prsMap)) {
    unname(preparedData$prsMap[[phenotype]])
    } else {
    character(0)
    }
    covariates <- unique(c(preparedData$covariates, prs))
    variables <- unique(c(
    preparedData$personVar, phenotype, covariates,
    preparedData$interactionTerm,
    if (!identical(config$correlationStructureValue, "none")) {
        config$correlationVarValue
    } else {
        character(0)
    }
    ))
    variables <- variables[!is.na(variables) & nzchar(variables)]
    missing <- setdiff(variables, colnames(modelData))
    if (length(missing)) {
    missing_text <- paste(missing, collapse = ", ")
    stop(sprintf(
        "Model variables not found for phenotype %s: %s",
        phenotype, missing_text
    ), call. = FALSE)
    }
    list(covariates = covariates, variables = variables)
}

lmePhenotypeFormulasDnaEpico <- function(
    preparedData, phenotype, covariates, modelData, config
) {
    display <- buildFormulaMethylationLME(
    phenotype, preparedData$personVar, covariates,
    preparedData$interactionTerm,
    includeRandomTerm = TRUE,
    responseVar = preparedData$internalResponseColumn
    )
    fitted <- buildFormulaMethylationLME(
    phenotype, preparedData$personVar, covariates,
    preparedData$interactionTerm,
    includeRandomTerm = !identical(config$engine, "nlme"),
    responseVar = preparedData$internalResponseColumn
    )
    validateFixedEffectDesignMethylationModels(
    formulaText = fitted, data = modelData,
    removeRandomEffects = !identical(config$engine, "nlme")
    )
    list(
    display = display, fitted = fitted,
    object = stats::as.formula(fitted, env = baseenv()),
    coefficientTerms = buildCoefficientTermMapMethylationModels(
        formulaText = fitted, data = modelData,
        removeRandomEffects = !identical(config$engine, "nlme")
    )
    )
}

addLmeCorrelationDataDnaEpico <- function(data, preparedData, config) {
    if (identical(config$correlationStructureValue, "none")) {
    return(list(data = data, timeVar = NULL))
    }
    correlation <- addCorrelationTimeVariableMethylationLME(
    data, config$correlationVarValue
    )
    validateCorrelationTimeMethylationLME(
    modelData = correlation$data,
    correlationTimeVar = correlation$correlationTimeVar,
    personVar = preparedData$personVar,
    correlationStructure = config$correlationStructureValue
    )
    list(data = correlation$data, timeVar = correlation$correlationTimeVar)
}

prepareLmePhenotypeSpecDnaEpico <- function(
    preparedData, modelData, phenotype, config
) {
    variables <- lmePhenotypeVariablesDnaEpico(
    preparedData, modelData, phenotype, config
    )
    data <- modelData[, variables$variables, drop = FALSE]
    data[[preparedData$personVar]] <- as.factor(
    data[[preparedData$personVar]]
    )
    for (var in intersect(preparedData$factorVars, colnames(data))) {
    data[[var]] <- as.factor(data[[var]])
    }
    formulas <- lmePhenotypeFormulasDnaEpico(
    preparedData, phenotype, variables$covariates, data, config
    )
    omnibus_target <- if (config$omnibus$test) {
    resolveOmnibusTargetTermMethylationLME(
        formulas$fitted, data, phenotype, preparedData$interactionTerm
    )
    } else {
    NULL
    }
    correlation <- addLmeCorrelationDataDnaEpico(data, preparedData, config)
    list(
    phenotype = phenotype, covariates = variables$covariates,
    data = correlation$data, correlationTimeVar = correlation$timeVar,
    formulas = formulas, omnibusTarget = omnibus_target
    )
}

lmePhenotypeSummaryPathDnaEpico <- function(phenotype, config) {
    if (is.null(config$summaryDir)) {
    return(NULL)
    }
    phenotypeSummaryPathMethylationModels(
    outputDir = config$summaryDir, phenotype = phenotype, analysis = "lme"
    )
}

lmePhenotypeSignatureDnaEpico <- function(
    preparedData, spec, config
) {
    buildPhenotypeSignatureMethylationModels(
    analysis = "lme", engine = config$engine,
    phenotype = spec$phenotype, formulaText = spec$formulas$display,
    preparedData = preparedData,
    modelSettings = list(
        fittedFormula = spec$formulas$fitted,
        personVar = preparedData$personVar, timeVar = preparedData$timeVar,
        random = paste0("~ 1 | ", preparedData$personVar), REML = TRUE,
        correlationStructure = config$correlationStructureValue,
        correlationVar = if (identical(
        config$correlationStructureValue, "none"
        )) {
        NULL
        } else {
        config$correlationVarValue
        },
        omnibusTest = config$omnibus$test, omnibusDdf = config$omnibus$ddf,
        omnibusTerm = spec$omnibusTarget,
        omnibusRhs = 0, omnibusJoint = TRUE
    ),
    packages = config$requiredPackages
    )
}

logResumedLmePhenotypeDnaEpico <- function(
    phenotype, formula, path, nCpGs, config
) {
    emitLogMinfiEwasWater(c(
    "============================================================",
    paste("Resumed phenotype:           ", phenotype),
    paste("Formula:                     ", formula),
    paste("Phenotype summary:           ", path),
    paste("CpGs restored:               ", nCpGs),
    "No CpG models were refitted.",
    "============================================================"
    ), verbose = config$verbose, log_path = config$logPath)
}

restoreLmePhenotypeDnaEpico <- function(
    state, phenotype, artifact, summaryPath, config, nCpGs
) {
    state$fits[phenotype] <- list(structure(
    list(),
    class = "dnaEPICO_compact_fit_index"
    ))
    state$summaryCache[[phenotype]] <- artifact$targetSummary
    state$coefficientResults[[phenotype]] <- artifact$coefficientResults
    state$omnibusTests[[phenotype]] <- artifact$omnibusTests
    state$modelMessages[[phenotype]] <- artifact$modelMessages
    state$fitFailures[[phenotype]] <- artifact$fitFailures
    state$formulas[[phenotype]] <- artifact$formula
    state$failureCounts[[phenotype]] <- artifact$failureCount
    state$failureReasons[[phenotype]] <- artifact$failureReasons
    state$phenotypeSummaries[[phenotype]] <- artifact
    state$summaryFiles[[phenotype]] <- summaryPath
    state$resumedPhenotypes <- c(state$resumedPhenotypes, phenotype)
    logResumedLmePhenotypeDnaEpico(
    phenotype, artifact$formula, summaryPath, nCpGs, config
    )
    invisible(TRUE)
}

tryResumeLmePhenotypeDnaEpico <- function(
    state, spec, signature, summaryPath, config, nCpGs
) {
    resumed <- list(object = NULL, reason = "resume was not requested")
    if (config$resumeFromSummary && !is.null(summaryPath)) {
    resumed <- loadPhenotypeSummaryMethylationModels(
        path = summaryPath, expectedSignature = signature
    )
    }
    if (!is.null(resumed$object)) {
    restoreLmePhenotypeDnaEpico(
        state, spec$phenotype, resumed$object, summaryPath,
        config, nCpGs
    )
    return(TRUE)
    }
    if (config$resumeFromSummary && !is.null(summaryPath)) {
    emitLogMinfiEwasWater(paste(
        "Phenotype summary was not reused for", spec$phenotype,
        ":", resumed$reason
    ), verbose = config$verbose, log_path = config$logPath)
    }
    FALSE
}

lmeBatchArgumentsDnaEpico <- function(preparedData, spec, config) {
    list(
    modelData = spec$data, formulaText = spec$formulas$fitted,
    personVar = preparedData$personVar, lmeEngine = config$engine,
    correlationStructure = config$correlationStructureValue,
    correlationTimeVar = spec$correlationTimeVar,
    phenotype = spec$phenotype,
    interactionTerm = preparedData$interactionTerm,
    responseVar = preparedData$internalResponseColumn,
    omnibusTest = config$omnibus$test,
    omnibusDdf = config$omnibus$ddf,
    omnibusTerm = spec$omnibusTarget,
    formulaObject = spec$formulas$object,
    coefficientTerms = spec$formulas$coefficientTerms
    )
}

updateLmePilotPlanDnaEpico <- function(
    state, preparedData, spec, config
) {
    cpgs <- preparedData$cpgColumns
    if (state$pilotCompleted || state$workerCount <= 1L || !length(cpgs)) {
    return(invisible(NULL))
    }
    pilot_cpgs <- utils::head(cpgs, 3L)
    common <- lmeBatchArgumentsDnaEpico(preparedData, spec, config)
    pilot <- measurePilotMemoryMethylationModels(function() {
    do.call(fitMethylationLMEBatch, c(list(
        cpgBatch = pilot_cpgs,
        data = preparedData$data[, pilot_cpgs, drop = FALSE]
    ), common))
    })
    state$parallelPlan <- refineParallelPlanWithPilotMethylationModels(
    state$parallelPlan, pilot$incrementalMB
    )
    state$backend <- state$parallelPlan$backend
    state$workerCount <- state$parallelPlan$workerCount
    state$batches <- chunkCpGColumnsMethylationModels(
    cpgs, state$workerCount,
    batchesPerCore = 8L
    )
    state$pilotCompleted <- TRUE
    invisible(gc(FALSE))
}

lmePsockDependencyNamesDnaEpico <- function() {
    c(
    "validateWorkerPackagesMethylationModels",
    "newMethylationFitErrorDnaEpico",
    "fitMethylationLMEBatch", "fitCpGModelMethylationLME",
    "fitLmeEngineDnaEpico", "coefficientTermsForLmeDnaEpico",
    "omnibusForLmeFitDnaEpico", "retainLmeFitDetailsDnaEpico",
    "completeCapturedLmeFitDnaEpico", "buildCorrelationMethylationLME",
    "coerceCoefficientTableMethylationLME",
    "buildCoefficientTermMapMethylationModels",
    "removeRandomInterceptMethylationModels",
    "computeOmnibusTestMethylationLME", "emptyOmnibusResultMethylationLME",
    "summarizeCpGFitMethylationLME", "captureModelConditionsDnaEpico",
    "combineModelMessagesDnaEpico", "modelMessageDnaEpico",
    "findCoefficientRowsMethylationLME",
    "findCoefficientRowsMethylationGLM", "quoteNamesMethylationGLM",
    "escapeRegexMethylationGLM", "cpgResponseMethylationModels",
    "compactCoefficientResultsMethylationModels",
    "coefficientNamesMethylationModels",
    "newCompactCoefficientStorageMethylationModels",
    "fillCompactCoefficientStorageMethylationModels",
    "collectBatchModelMessagesMethylationModels",
    "collectBatchFitFailuresMethylationModels",
    "emptyModelMessagesDnaEpico", "emptyFitFailuresMethylationModels",
    "collectOmnibusTestsMethylationLME"
    )
}

ensureLmePsockClusterDnaEpico <- function(state, config) {
    if (!identical(state$backend, "psock") || length(state$batches) <= 1L ||
    !is.null(state$cluster)) {
    return(invisible(NULL))
    }
    state$cluster <- makePsockClusterMethylationModels(min(
    state$workerCount, length(state$batches)
    ))
    parallel::clusterExport(
    state$cluster,
    varlist = lmePsockDependencyNamesDnaEpico(),
    envir = environment()
    )
    lib_path <- config$libPath
    packages <- config$requiredPackages
    parallel::clusterExport(
    state$cluster, c("lib_path", "packages"),
    envir = environment()
    )
    parallel::clusterEvalQ(
    state$cluster,
    validateWorkerPackagesMethylationModels(lib_path, packages)
    )
    invisible(NULL)
}

runSerialLmeBatchesDnaEpico <- function(
    batches, preparedData, common, config, worker
) {
    validateWorkerPackagesMethylationModels(
    config$libPath, config$requiredPackages
    )
    lapply(batches, function(batch) {
    do.call(worker, c(list(
        cpgBatch = batch,
        data = preparedData$data[, batch, drop = FALSE]
    ), common))
    })
}

runForkLmeBatchesDnaEpico <- function(
    batches, preparedData, common, config, worker, workers
) {
    parallel::mclapply(batches, function(batch) {
    validateWorkerPackagesMethylationModels(
        config$libPath, config$requiredPackages
    )
    do.call(worker, c(list(
        cpgBatch = batch,
        data = preparedData$data[, batch, drop = FALSE]
    ), common))
    }, mc.cores = min(workers, length(batches)), mc.preschedule = TRUE)
}

runPsockLmeBatchesDnaEpico <- function(
    state, preparedData, common, worker
) {
    state$clusterUseCount <- state$clusterUseCount + 1L
    parallel::clusterExport(
    state$cluster, c("common", "worker"),
    envir = environment()
    )
    results <- vector("list", length(state$batches))
    cluster_size <- min(state$workerCount, length(state$batches))
    waves <- split(
    seq_along(state$batches),
    ceiling(seq_along(state$batches) / cluster_size)
    )
    for (wave in waves) {
    tasks <- lapply(wave, function(index) {
        list(
        cpgBatch = state$batches[[index]],
        responses = as.matrix(preparedData$data[,
            state$batches[[index]],
            drop = FALSE
        ])
        )
    })
    wave_results <- parallel::parLapplyLB(
        state$cluster, tasks, function(task) {
        do.call(worker, c(list(
            cpgBatch = task$cpgBatch, data = task$responses
        ), common))
        }
    )
    results[wave] <- wave_results
    invisible(gc(FALSE))
    }
    results
}

runLmePhenotypeBatchesDnaEpico <- function(
    state, preparedData, spec, config
) {
    ensureLmePsockClusterDnaEpico(state, config)
    worker <- fitMethylationLMEBatch
    common <- lmeBatchArgumentsDnaEpico(preparedData, spec, config)
    if (identical(state$backend, "fork") && length(state$batches) > 1L) {
    return(runForkLmeBatchesDnaEpico(
        state$batches, preparedData, common, config,
        worker, state$workerCount
    ))
    }
    if (identical(state$backend, "psock") && length(state$batches) > 1L) {
    environment(worker) <- .GlobalEnv
    return(runPsockLmeBatchesDnaEpico(
        state, preparedData, common, worker
    ))
    }
    runSerialLmeBatchesDnaEpico(
    state$batches, preparedData, common, config, worker
    )
}

combineLmeBatchResultsDnaEpico <- function(
    batchResults, cpgColumns, phenotype
) {
    coefficients <- combineCompactCoefficientResultsMethylationModels(
    batchResults,
    cpgOrder = cpgColumns
    )
    summaries <- combineBatchTablesMethylationModels(
    batchResults, "summaries", data.frame()
    )
    omnibus <- combineBatchTablesMethylationModels(
    batchResults, "omnibusTests",
    collectOmnibusTestsMethylationLME(list(), phenotype)
    )
    summary_cache <- filterSummaryByPvalueMethylationLME(
    summaries, NA_real_
    )
    messages <- combineBatchTablesMethylationModels(
    batchResults, "modelMessages", emptyModelMessagesDnaEpico()
    )
    failures <- combineBatchTablesMethylationModels(
    batchResults, "fitFailures", emptyFitFailuresMethylationModels()
    )
    messages <- messages[match(cpgColumns, messages$CpG), , drop = FALSE]
    rownames(messages) <- NULL
    errors <- if (!nrow(failures)) {
    integer(0)
    } else {
    sort(table(failures$Error), decreasing = TRUE)
    }
    list(
    coefficients = coefficients, summaryCache = summary_cache,
    omnibus = omnibus, messages = messages, failures = failures,
    errors = errors, pValueAvailable = messages$P.Value.Available,
    failureCount = sum(!messages$P.Value.Available)
    )
}

lmePhenotypeArtifactDnaEpico <- function(
    preparedData, spec, combined, signature, config
) {
    factor_vars <- preparedData$factorVars[
    preparedData$factorVars %in% colnames(preparedData$data)
    ]
    assemblePhenotypeSummaryMethylationModels(
    analysis = "lme", engine = config$engine,
    phenotype = spec$phenotype, signature = signature,
    cpgOrder = preparedData$cpgColumns,
    coefficientResults = combined$coefficients,
    targetSummary = combined$summaryCache,
    omnibusTests = combined$omnibus,
    modelMessages = combined$messages, fitFailures = combined$failures,
    failureCount = combined$failureCount,
    failureReasons = combined$errors,
    formulaText = spec$formulas$display,
    settings = list(
        fittedFormula = spec$formulas$fitted,
        methylationScale = preparedData$methylationScale,
        responseLabel = preparedData$responseLabel,
        interactionTerm = preparedData$interactionTerm,
        covariates = spec$covariates, factorVars = preparedData$factorVars,
        factorLevels = lapply(preparedData$data[factor_vars], levels),
        scaleVars = preparedData$scaleVars,
        scalingMetadata = preparedData$scalingMetadata,
        personVar = preparedData$personVar, timeVar = preparedData$timeVar,
        random = paste0("~ 1 | ", preparedData$personVar), REML = TRUE,
        correlationStructure = config$correlationStructureValue,
        correlationVar = if (identical(
        config$correlationStructureValue, "none"
        )) {
        NULL
        } else {
        config$correlationVarValue
        },
        omnibusTest = config$omnibus$test, omnibusDdf = config$omnibus$ddf,
        omnibusTerm = spec$omnibusTarget,
        omnibusRhs = 0, omnibusJoint = TRUE
    )
    )
}

storeFittedLmePhenotypeDnaEpico <- function(
    state, spec, combined, artifact, summaryPath
) {
    phenotype <- spec$phenotype
    state$fits[phenotype] <- list(structure(
    list(),
    class = "dnaEPICO_compact_fit_index"
    ))
    state$summaryCache[[phenotype]] <- combined$summaryCache
    state$coefficientResults[[phenotype]] <- combined$coefficients
    state$omnibusTests[[phenotype]] <- combined$omnibus
    state$modelMessages[[phenotype]] <- combined$messages
    state$fitFailures[[phenotype]] <- combined$failures
    state$formulas[[phenotype]] <- spec$formulas$display
    state$failureCounts[[phenotype]] <- combined$failureCount
    state$failureReasons[[phenotype]] <- combined$errors
    state$phenotypeSummaries[[phenotype]] <- artifact
    state$fittedPhenotypes <- c(state$fittedPhenotypes, phenotype)
    if (!is.null(summaryPath)) {
    savePhenotypeSummaryMethylationModels(artifact, summaryPath)
    state$summaryFiles[[phenotype]] <- summaryPath
    }
    invisible(NULL)
}

lmeFitResourceLogLinesDnaEpico <- function(state) {
    plan <- state$parallelPlan
    c(
    paste("Parallel backend:            ", state$backend),
    paste("Effective workers:           ", state$workerCount),
    paste("Parallel crossover CpGs:     ", plan$crossoverCpGs),
    paste("Parallel selection:          ", plan$reason),
    paste(
        "Available memory (MB):       ",
        if (is.finite(plan$availableMemoryMB)) {
        round(plan$availableMemoryMB)
        } else {
        "unknown"
        }
    ),
    paste(
        "Estimated memory per worker: ",
        round(plan$estimatedWorkerMemoryMB), "MB"
    ),
    paste(
        "Pilot incremental memory:    ",
        if (is.finite(plan$pilotMemoryMB)) {
        paste(round(plan$pilotMemoryMB), "MB")
        } else {
        "not required"
        }
    ),
    paste("Memory worker cap:           ", plan$memoryWorkerCap),
    paste("Fit batches:                 ", length(state$batches)),
    paste(
        "Fit batch size:              ",
        if (!length(state$batches)) {
        0L
        } else {
        max(vapply(state$batches, length, integer(1)))
        }
    )
    )
}

logFittedLmePhenotypeDnaEpico <- function(
    state, preparedData, spec, combined, config
) {
    emitLogMinfiEwasWater(c(
    "============================================================",
    paste("Fitted phenotype:            ", spec$phenotype),
    paste("Formula:                     ", spec$formulas$display),
    paste(
        "Correlation structure:       ",
        config$correlationStructureValue
    ),
    paste(
        "Correlation variable:        ",
        if (identical(config$correlationStructureValue, "none")) {
        "None"
        } else {
        config$correlationVarValue
        }
    ),
    paste("CpGs attempted:              ", length(preparedData$cpgColumns)),
    paste("CpGs without p-values:       ", combined$failureCount),
    paste("Omnibus tests requested:     ", config$omnibus$test),
    paste(
        "Omnibus target term:         ",
        if (is.null(spec$omnibusTarget)) "None" else spec$omnibusTarget
    ),
    paste(
        "Omnibus denominator DF:      ",
        if (config$omnibus$test) config$omnibus$ddf else "None"
    ),
    paste(
        "Successful omnibus tests:    ",
        sum(combined$omnibus$Omnibus.Status == "tested")
    ),
    paste(
        "Unavailable omnibus tests:   ",
        sum(combined$omnibus$Omnibus.Status != "tested")
    ),
    paste(
        "Top fit errors:              ",
        formatFitErrorsMethylationModels(combined$errors)
    ),
    lmeFitResourceLogLinesDnaEpico(state),
    paste("Fit-time summary rows cached:", nrow(combined$summaryCache)),
    "============================================================"
    ), verbose = config$verbose, log_path = config$logPath)
}

warnUnavailableLmePvaluesDnaEpico <- function(
    preparedData, spec, combined, config
) {
    if (!length(preparedData$cpgColumns) || any(combined$pValueAvailable)) {
    return(invisible(NULL))
    }
    warning(sprintf(
    "%s %s p-values were available for phenotype '%s'. %s: %s. %s",
    "No CpG", toupper(config$engine), spec$phenotype,
    "Top failure reasons", formatFitErrorsMethylationModels(combined$errors),
    "The failure inventory was retained and the analysis continued."
    ), call. = FALSE)
    invisible(NULL)
}

combineLmeFitAuditTablesDnaEpico <- function(state) {
    list(
    failures = combineBatchTablesMethylationModels(
        lapply(state$fitFailures, function(x) list(table = x)),
        "table", emptyFitFailuresMethylationModels()
    ),
    messages = combineBatchTablesMethylationModels(
        lapply(state$modelMessages, function(x) list(table = x)),
        "table", emptyModelMessagesDnaEpico()
    )
    )
}

lmeFitSettingsDnaEpico <- function(state, preparedData, config) {
    plan <- state$parallelPlan
    list(
    nCores = config$nCoresValue, parallelBackend = state$backend,
    workerCount = state$workerCount,
    resourceWorkerCap = plan$resourceWorkerCap,
    parallelCrossoverCpGs = plan$crossoverCpGs,
    parallelSelectionReason = plan$reason,
    availableMemoryMB = plan$availableMemoryMB,
    reservedMemoryMB = plan$reservedMemoryMB,
    estimatedWorkerMemoryMB = plan$estimatedWorkerMemoryMB,
    pilotIncrementalMemoryMB = plan$pilotMemoryMB,
    memoryWorkerCap = plan$memoryWorkerCap,
    clusterReusedAcrossPhenotypes = state$clusterUseCount > 1L,
    fitBatchCount = length(state$batches), libPath = config$libPath,
    lmeLibs = config$requestedPackages, lmeEngine = config$engine,
    correlationStructure = config$correlationStructureValue,
    correlationVar = if (identical(
        config$correlationStructureValue, "none"
    )) {
        NULL
    } else {
        config$correlationVarValue
    },
    methylationScale = preparedData$methylationScale,
    methylationObjectPrefix = preparedData$methylationObjectPrefix,
    responseLabel = preparedData$responseLabel,
    internalResponseColumn = preparedData$internalResponseColumn,
    interactionTerm = preparedData$interactionTerm,
    omnibusTest = config$omnibus$test, omnibusDdf = config$omnibus$ddf,
    omnibusRhs = 0, omnibusJoint = TRUE,
    phenotypes = preparedData$phenotypes,
    covariates = preparedData$covariates,
    factorVars = preparedData$factorVars,
    factorLevels = lapply(stats::setNames(
        preparedData$factorVars, preparedData$factorVars
    ), function(variable) levels(preparedData$data[[variable]])),
    scaleVars = preparedData$scaleVars,
    scalingMetadata = preparedData$scalingMetadata,
    sampleCount = nrow(preparedData$data), personVar = preparedData$personVar,
    timeVar = preparedData$timeVar
    )
}

newMethylationLMEModelsDnaEpico <- function(state, preparedData, config) {
    audit <- combineLmeFitAuditTablesDnaEpico(state)
    structure(list(
    fits = state$fits, summaryCache = state$summaryCache,
    coefficientResults = state$coefficientResults,
    phenotypeSummaries = state$phenotypeSummaries,
    summaryFiles = state$summaryFiles,
    resumedPhenotypes = state$resumedPhenotypes,
    fittedPhenotypes = state$fittedPhenotypes,
    omnibusTests = state$omnibusTests,
    omnibusTargets = state$omnibusTargets,
    formulas = state$formulas, phenotypes = names(state$fits),
    failureCounts = state$failureCounts,
    failureReasons = state$failureReasons,
    fitFailures = audit$failures, modelMessages = audit$messages,
    settings = lmeFitSettingsDnaEpico(state, preparedData, config),
    responseLabel = preparedData$responseLabel
    ), class = "dnaEPICO_methylationLME_models")
}

fitAllLmePhenotypesDnaEpico <- function(
    state, preparedData, modelData, config
) {
    for (phenotype in preparedData$phenotypes) {
    spec <- prepareLmePhenotypeSpecDnaEpico(
        preparedData, modelData, phenotype, config
    )
    if (config$omnibus$test) {
        state$omnibusTargets[[phenotype]] <- spec$omnibusTarget
    }
    signature <- lmePhenotypeSignatureDnaEpico(
        preparedData, spec, config
    )
    summary_path <- lmePhenotypeSummaryPathDnaEpico(phenotype, config)
    if (tryResumeLmePhenotypeDnaEpico(
        state, spec, signature, summary_path, config,
        length(preparedData$cpgColumns)
    )) {
        next
    }
    updateLmePilotPlanDnaEpico(state, preparedData, spec, config)
    batches <- runLmePhenotypeBatchesDnaEpico(
        state, preparedData, spec, config
    )
    combined <- combineLmeBatchResultsDnaEpico(
        batches, preparedData$cpgColumns, phenotype
    )
    artifact <- lmePhenotypeArtifactDnaEpico(
        preparedData, spec, combined, signature, config
    )
    storeFittedLmePhenotypeDnaEpico(
        state, spec, combined, artifact, summary_path
    )
    logFittedLmePhenotypeDnaEpico(
        state, preparedData, spec, combined, config
    )
    warnUnavailableLmePvaluesDnaEpico(
        preparedData, spec, combined, config
    )
    }
    invisible(NULL)
}

#' Fit CpG-wise mixed-effects models for longitudinal methylation analyses
#'
#' @param preparedData Object returned by `prepareMethylationLMEData()`.
#' @param nCores Integer. Maximum number of worker processes to use. Automatic
#'   fitting uses the lme4 or nlme crossover and caps workers by available CpGs,
#'   CPUs, and detected memory.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded
#'   to worker processes.
#' @param lmeLibs Character vector or comma-separated string of package names to
#'   check on worker processes. The default is `'lme4,lmerTest'`.
#' @param correlationStructure Character. Residual correlation structure used
#'   when `lmeLibs` selects `'nlme'`. One of `'none'`, `'AR1'`, or `'CAR1'`.
#' @param correlationVar Character or `NULL`. Variable used to order repeated
#'   observations within `personVar` for `AR1` or `CAR1` residual correlation
#'   structures. Must be supplied explicitly for `AR1` or `CAR1`.
#' @param omnibusTest Logical. If `TRUE`, calculate one joint fixed-effect test
#'   per CpG for the phenotype-by-interaction term, or the phenotype main effect
#'   when no interaction is requested.
#' @param omnibusDdf Character. Denominator degrees-of-freedom method used by
#'   `lmerTest::contestMD()`: `'Satterthwaite'` or `'Kenward-Roger'`.
#' @param summaryDir Character or `NULL`. Directory used for one complete
#'   compact summary per phenotype. `NULL` disables disk persistence.
#' @param resumeFromSummary Logical. If `TRUE`, reuse a complete summary in
#'   `summaryDir` when it was generated from the same input file and model
#'   configuration.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationLME_models'`
#'   containing compact coefficient matrices, unfiltered target summaries,
#'   formulas, model conditions, omnibus results, hard errors, and phenotype
#'   summary artifacts.
#'
#' @description
#' Fit one linear mixed-effects model per CpG for each phenotype requested in
#' the object returned by `prepareMethylationLMEData()`. Each native fit is
#' reduced to compact numerical results and discarded before the next batch is
#' returned.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' model_results <- fitMethylationLMEModels(
#'   preparedData = ex$preparedData,
#'   nCores = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(model_results$fits)
#'
#' @export
fitMethylationLMEModels <- function(
    preparedData, nCores = 1L,
    libPath = NULL, lmeLibs = "lme4,lmerTest", correlationStructure = "none",
    correlationVar = NULL, omnibusTest = FALSE,
    omnibusDdf = "Satterthwaite", summaryDir = NULL,
    resumeFromSummary = TRUE, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationLME.txt"
) {
    config <- normalizeLmeFitConfigDnaEpico(
    as.list(environment(), all.names = TRUE), preparedData
    )
    model_data <- if (is.null(preparedData$modelData)) {
    preparedData$data
    } else {
    preparedData$modelData
    }
    state <- newLmeFitStateDnaEpico(preparedData, config)
    on.exit(stopLmeFitClusterDnaEpico(state), add = TRUE)
    fitAllLmePhenotypesDnaEpico(state, preparedData, model_data, config)
    stopLmeFitClusterDnaEpico(state)
    newMethylationLMEModelsDnaEpico(state, preparedData, config)
}

summarizeOmnibusTestsMethylationLME <- function(
    modelResults, padjmethod = "fdr"
) {
    adjustment_method <- validatePAdjustmentMethodMethylationModels(
    padjmethod
    )
    omnibus_tables <- modelResults$omnibusTests
    if (!is.list(omnibus_tables)) {
    omnibus_tables <- list()
    }

    phenotype_names <- modelResults$phenotypes
    if (is.null(phenotype_names)) {
    phenotype_names <- names(modelResults$fits)
    }
    summaries <- lapply(phenotype_names, function(phenotype) {
    table <- omnibus_tables[[phenotype]]
    if (!is.data.frame(table) || nrow(table) == 0L) {
        return(data.frame())
    }

    table$Omnibus.Adjusted.P.Value <- NA_real_
    valid <- table$Omnibus.Status == "tested" &
        is.finite(table$Omnibus.P.Value)
    table$Omnibus.Adjusted.P.Value[valid] <- stats::p.adjust(
        table$Omnibus.P.Value[valid],
        method = adjustment_method
    )
    table
    })
    names(summaries) <- phenotype_names
    summaries
}

normalizeLmeSummaryConfigDnaEpico <- function(config) {
    config$logPath <- resolveLogPathMinfiEwasWater(
    logs = config$logs, log_dir = config$log_dir,
    log_file = config$log_file
    )
    config$pValueFilter <- validateProbabilityDnaEpico(
    normalizeOptionalNumericMethylationGLM(config$summaryPval),
    "summaryPval",
    allowNA = TRUE
    )
    config$chunkSizeValue <- normalizeChunkSizeMethylationGLM(config$chunkSize)
    config$nCoresValue <- validatePositiveIntegerMethylationModels(
    config$nCores, "nCores"
    )
    config$adjustmentMethod <- validatePAdjustmentMethodMethylationModels(
    config$padjmethod
    )
    config
}

summarizeLmeChunkDnaEpico <- function(
    chunk, fits, phenotype, interactionTerm, worker
) {
    rows <- lapply(chunk, function(cpg) {
    worker(
        cpg = cpg, modelObj = fits[[cpg]], phenotype = phenotype,
        interactionTerm = interactionTerm
    )
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) NULL else do.call(rbind, rows)
}

runLmeSummaryChunksDnaEpico <- function(
    fits, phenotype, interactionTerm, nCores, chunkSize
) {
    cpg_names <- names(fits)
    local_size <- chunkSize
    if (is.null(local_size)) {
    local_size <- max(
        10L, floor(length(cpg_names) / max(nCores * 4L, 1L))
    )
    }
    local_size <- max(1L, as.integer(local_size))
    chunks <- split(
    cpg_names, ceiling(seq_along(cpg_names) / local_size)
    )
    worker <- summarizeCpGFitMethylationLME
    chunk_worker <- summarizeLmeChunkDnaEpico
    run_chunk <- function(chunk) {
    chunk_worker(
        chunk, fits, phenotype, interactionTerm, worker
    )
    }
    if (nCores > 1L && length(chunks) > 1L) {
    environment(worker) <- .GlobalEnv
    environment(chunk_worker) <- .GlobalEnv
    cluster <- parallel::makeCluster(min(nCores, length(chunks)))
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    parallel::clusterExport(cluster, varlist = c(
        "fits", "phenotype", "interactionTerm", "worker", "chunk_worker",
        "summarizeLmeChunkDnaEpico", "findCoefficientRowsMethylationLME",
        "findCoefficientRowsMethylationGLM", "modelMessageDnaEpico"
    ), envir = environment())
    results <- parallel::parLapplyLB(cluster, chunks, run_chunk)
    } else {
    results <- lapply(chunks, run_chunk)
    }
    results <- Filter(Negate(is.null), results)
    summary <- if (!length(results)) data.frame() else do.call(rbind, results)
    rownames(summary) <- NULL
    list(summary = summary, chunkSize = local_size)
}

summarizeLmePhenotypeDnaEpico <- function(
    modelResults, preparedData, phenotype, config
) {
    cached <- modelResults$summaryCache[[phenotype]]
    if (!is.null(cached)) {
    diagnostic <- filterSummaryByPvalueMethylationLME(cached, NA_real_)
    summary <- filterSummaryByPvalueMethylationLME(
        diagnostic, config$pValueFilter
    )
    return(list(
        summary = summary, diagnostic = diagnostic,
        source = "fit-time cache", chunkSize = NULL
    ))
    }
    chunked <- runLmeSummaryChunksDnaEpico(
    fits = modelResults$fits[[phenotype]], phenotype = phenotype,
    interactionTerm = preparedData$interactionTerm,
    nCores = config$nCoresValue, chunkSize = config$chunkSizeValue
    )
    diagnostic <- chunked$summary
    summary <- filterSummaryByPvalueMethylationLME(
    diagnostic, config$pValueFilter
    )
    list(
    summary = summary, diagnostic = diagnostic,
    source = "model fits", chunkSize = chunked$chunkSize
    )
}

logLmePhenotypeSummaryDnaEpico <- function(phenotype, result, config) {
    lines <- c(
    "============================================================",
    paste("Summarized phenotype:        ", phenotype),
    paste("LME summary rows returned:   ", nrow(result$summary)),
    paste("Summary source:              ", result$source)
    )
    if (!is.null(result$chunkSize)) {
    lines <- c(lines, paste(
        "Summary chunk size:          ", result$chunkSize
    ))
    }
    lines <- c(
    lines,
    if (is.na(config$pValueFilter)) {
        "P-value filter:              none"
    } else {
        paste("P-value filter:              ", config$pValueFilter)
    },
    "============================================================"
    )
    emitLogMinfiEwasWater(
    lines,
    verbose = config$verbose, log_path = config$logPath
    )
}

newMethylationLMESummariesDnaEpico <- function(
    summaries, diagnostics, omnibus, modelResults, preparedData, config
) {
    messages <- modelResults$modelMessages
    if (is.null(messages)) {
    messages <- collectModelMessagesMethylationLME(modelResults$fits)
    }
    structure(list(
    summaries = summaries, diagnosticSummaries = diagnostics,
    omnibusTests = omnibus, phenotypes = names(summaries),
    fitFailures = modelResults$fitFailures, modelMessages = messages,
    settings = list(
        summaryPval = config$pValueFilter,
        padjmethod = config$adjustmentMethod,
        chunkSize = config$chunkSizeValue,
        interactionTerm = preparedData$interactionTerm,
        lmeEngine = modelResults$settings$lmeEngine,
        omnibusTest = isTRUE(modelResults$settings$omnibusTest),
        factorLevels = modelResults$settings$factorLevels
    )
    ), class = "dnaEPICO_methylationLME_summaries")
}

#' Summarize CpG-wise mixed-effects results for longitudinal analyses
#'
#' @param modelResults Object returned by `fitMethylationLMEModels()`.
#' @param preparedData Object returned by `prepareMethylationLMEData()`.
#' @param summaryPval Numeric or `NA`. Optional p-value filter applied to the
#'   returned summary tables. `NA` keeps all rows.
#' @param padjmethod Character. Adjustment method passed to `stats::p.adjust()`
#'   for omnibus p-values across CpGs within each phenotype and tested term.
#' @param nCores Integer. Number of worker processes to use while extracting
#'   summary rows.
#' @param chunkSize Integer or `NULL`. Number of CpGs processed per parallel
#'   chunk. `NULL` chooses a value automatically.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationLME_summaries'`
#'   containing the optionally filtered summary tables in `summaries` and the
#'   complete CpG-level tables in `diagnosticSummaries`. Diagnostics,
#'   annotation, and report output use the complete tables so `summaryPval`
#'   does not remove CpGs from those outputs. `modelMessages` retains native
#'   messages, warnings, and errors for every attempted CpG.
#'
#' @description
#' Return phenotype-specific fixed-effect tables from the compact fit-time
#' results produced by `fitMethylationLMEModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' summary_results <- summarizeMethylationLMEModels(
#'   modelResults = ex$modelResults,
#'   preparedData = ex$preparedData,
#'   summaryPval = NA,
#'   nCores = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(summary_results$summaries)
#'
#' @export
summarizeMethylationLMEModels <- function(
    modelResults, preparedData,
    summaryPval = NA, padjmethod = "fdr", nCores = 1L, chunkSize = NULL,
    verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationLME.txt"
) {
    config <- normalizeLmeSummaryConfigDnaEpico(
    as.list(environment(), all.names = TRUE)
    )
    summaries <- list()
    diagnostic_summaries <- list()
    for (phenotype in names(modelResults$fits)) {
    result <- summarizeLmePhenotypeDnaEpico(
        modelResults, preparedData, phenotype, config
    )
    summaries[[phenotype]] <- result$summary
    diagnostic_summaries[[phenotype]] <- result$diagnostic
    logLmePhenotypeSummaryDnaEpico(phenotype, result, config)
    }
    omnibus <- summarizeOmnibusTestsMethylationLME(
    modelResults = modelResults, padjmethod = config$adjustmentMethod
    )
    newMethylationLMESummariesDnaEpico(
    summaries, diagnostic_summaries, omnibus,
    modelResults, preparedData, config
    )
}

lmeCoefficientTableForCpgDnaEpico <- function(modelResults, phenotype, cpg) {
    table <- coefficientTableFromCompactMethylationModels(
    modelResults$coefficientResults[[phenotype]], cpg
    )
    if (!is.null(table)) {
    return(table)
    }
    fit <- modelResults$fits[[phenotype]][[cpg]]
    if (is.null(fit) || inherits(fit, "dnaEPICO_methylationLME_fit_error") ||
    is.null(fit$coef)) {
    return(NULL)
    }
    as.data.frame(fit$coef)
}

collectOmnibusLmeHitsDnaEpico <- function(
    modelResults, phenotype, threshold
) {
    table <- modelResults$omnibusTests[[phenotype]]
    if (is.data.frame(table) && nrow(table) > 0L) {
    hit_cpgs <- table$CpG[
        table$Omnibus.Status == "tested" &
        is.finite(table$Omnibus.P.Value) &
        table$Omnibus.P.Value < threshold
    ]
    } else {
    fits <- modelResults$fits[[phenotype]]
    tested <- vapply(fits, function(fit) {
        !is.null(fit) && !is.null(fit$omnibus) &&
        identical(fit$omnibus$status, "tested") &&
        is.finite(fit$omnibus$pValue) &&
        fit$omnibus$pValue < threshold
    }, logical(1))
    hit_cpgs <- names(fits)[tested]
    }
    hits <- lapply(unique(hit_cpgs), function(cpg) {
    lmeCoefficientTableForCpgDnaEpico(modelResults, phenotype, cpg)
    })
    stats::setNames(Filter(Negate(is.null), hits), unique(hit_cpgs)[
    !vapply(hits, is.null, logical(1))
    ])
}

collectCachedLmeHitsDnaEpico <- function(
    modelResults, phenotype, threshold, interactionTerm
) {
    cached <- modelResults$summaryCache[[phenotype]]
    if (is.null(cached) || !optionalTermMatchesMethylationModels(
    requested = interactionTerm,
    cached = modelResults$settings$interactionTerm
    )) {
    return(NULL)
    }
    hit_cpgs <- character(0)
    if (nrow(cached) > 0L && !is.na(threshold)) {
    hit_cpgs <- unique(cached$CpG[cached$P.value < threshold])
    hit_cpgs <- hit_cpgs[!is.na(hit_cpgs)]
    }
    hits <- lapply(hit_cpgs, function(cpg) {
    lmeCoefficientTableForCpgDnaEpico(modelResults, phenotype, cpg)
    })
    stats::setNames(Filter(Negate(is.null), hits), hit_cpgs[
    !vapply(hits, is.null, logical(1))
    ])
}

collectDirectLmeHitsDnaEpico <- function(
    modelResults, phenotype, threshold, interactionTerm
) {
    fits <- modelResults$fits[[phenotype]]
    hits <- list()
    for (cpg in names(fits)) {
    fit <- fits[[cpg]]
    if (is.null(fit) || inherits(
        fit, "dnaEPICO_methylationLME_fit_error"
    ) || is.null(fit$coef)) {
        next
    }
    matched <- findCoefficientRowsMethylationLME(
        coefNames = rownames(fit$coef), phenotype = phenotype,
        interactionTerm = interactionTerm,
        coefficientTerms = fit$coefficientTerms
    )
    if (length(matched) && any(
        fit$coef[matched, "Pr(>|t|)", drop = TRUE] < threshold,
        na.rm = TRUE
    )) {
        hits[[cpg]] <- as.data.frame(fit$coef)
    }
    }
    hits
}

#' Collect significant longitudinal terms from compact mixed-effects results
#'
#' @param modelResults Object returned by `fitMethylationLMEModels()`.
#' @param pvalThreshold Numeric. Threshold applied to the extracted phenotype or
#'   interaction p-values.
#' @param interactionTerm Character or `NULL`. Optional interaction term. When
#'   `NULL`, phenotype main effects are used.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationLME_significant'`
#'   containing the retained coefficient tables for each phenotype.
#'
#' @description
#' Collect raw coefficient tables for CpGs whose phenotype main effect or
#' requested interaction p-value passes the chosen threshold.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' significant_hits <- collectSignificantInteractionsMethylationLME(
#'   modelResults = ex$modelResults,
#'   pvalThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(significant_hits)
#'
#' @export
collectSignificantInteractionsMethylationLME <- function(
    modelResults,
    pvalThreshold = 0.05, interactionTerm = NULL,
    verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationLME.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    threshold <- validateProbabilityDnaEpico(pvalThreshold, "pvalThreshold")
    phenotypes <- modelResults$phenotypes
    if (is.null(phenotypes)) {
    phenotypes <- names(modelResults$fits)
    }
    use_omnibus <- isTRUE(modelResults$settings$omnibusTest)
    if (use_omnibus && length(phenotypes) &&
        !optionalTermMatchesMethylationModels(
    requested = interactionTerm,
    cached = modelResults$settings$interactionTerm
    )) {
    stop(
        "interactionTerm does not match the term used for the fitted ",
        "omnibus tests.",
        call. = FALSE
    ) }
    retained <- stats::setNames(lapply(phenotypes, function(phenotype) {
    if (use_omnibus) {
        return(collectOmnibusLmeHitsDnaEpico(
        modelResults, phenotype, threshold
        ))
    }
    cached <- collectCachedLmeHitsDnaEpico(
        modelResults, phenotype, threshold, interactionTerm
    )
    if (!is.null(cached)) {
        return(cached)
    }
    collectDirectLmeHitsDnaEpico(
        modelResults, phenotype, threshold, interactionTerm
    )
    }), phenotypes)
    hit_counts <- vapply(retained, length, integer(1))
    emitLogMinfiEwasWater(c(
    "============================================================",
    paste("Significant longitudinal terms retained at p <", threshold, ":"),
    paste(names(hit_counts), hit_counts, sep = ": ", collapse = "; "),
    "============================================================"
    ), verbose = verbose, log_path = log_path)
    structure(retained, class = "dnaEPICO_methylationLME_significant")
}
lmeDiagnosticFilesDnaEpico <- function(outputDir, fileKey, plots) {
    files <- list(
    qqplot = NULL, standardError = NULL,
    standardErrorSignificance = NULL, volcano = NULL,
    effectForest = NULL
    )
    if (is.null(outputDir)) {
    return(files)
    }
    files$qqplot <- file.path(
    outputDir, paste0("qqplot_", fileKey, "_coefficientPvalue.tiff")
    )
    files$standardError <- file.path(
    outputDir,
    paste0("standardError_", fileKey, "_byAverageMethylation.tiff")
    )
    files$standardErrorSignificance <- file.path(
    outputDir,
    paste0("standardErrorSignificance_", fileKey, "_byPvalue.tiff")
    )
    if (!is.null(plots$volcano)) {
    files$volcano <- file.path(
        outputDir, paste0("volcano_", fileKey, "_coefficientEstimate.tiff")
    )
    }
    if (!is.null(plots$effectForest)) {
    files$effectForest <- file.path(
        outputDir,
        paste0("effectForest_", fileKey, "_coefficientEstimate95CI.tiff")
    )
    }
    files
}

drawLmeDiagnosticPlotsDnaEpico <- function(plots, files, config) {
    for (plot_name in names(plots)) {
    plot_object <- plots[[plot_name]]
    if (is.null(plot_object)) {
        next
    }
    square <- identical(plot_name, "qqplot")
    height <- if (square) {
        max(config$plotWidth, config$plotHeight)
    } else if (identical(plot_name, "effectForest")) {
        max(config$plotHeight, 1400L)
    } else {
        config$plotHeight
    }
    runPlotMinfiEwasWater(
        draw_fun = function() drawPlotObjectMinfiEwasWater(plot_object),
        display = config$display, file = files[[plot_name]],
        width = if (square) height else config$plotWidth,
        height = height, res = config$plotDPI
    )
    }
    invisible(NULL)
}

buildLmeTermDiagnosticDnaEpico <- function(
    summaryData, phenotype, term, termIndex, multipleTerms, diagnosticMean,
    config
) {
    diagnostic <- buildMethylationTermDiagnosticsDnaEpico(
    summaryData = summaryData, phenotype = phenotype, term = term,
    termColumn = "Interaction.Term", pValueColumn = "P.value",
    yColumn = "Std.Error", yLabel = "Standard Error",
    diagnosticMean = diagnosticMean,
    fdrThreshold = config$fdrThreshold,
    estimateColumn = "Value", standardErrorColumn = "Std.Error"
    )
    if (is.null(diagnostic)) {
    return(NULL)
    }
    names(diagnostic$plots)[names(diagnostic$plots) == "residualSD"] <-
    "standardError"
    names(diagnostic$plots)[
    names(diagnostic$plots) == "residualSignificance"
    ] <- "standardErrorSignificance"
    file_key <- phenotype
    if (multipleTerms) {
    file_key <- paste0(
        phenotype, "_", sprintf("%02d", termIndex), "_",
        sanitizeDiagnosticTermDnaEpico(term)
    )
    }
    files <- lmeDiagnosticFilesDnaEpico(
    config$outputDir, file_key, diagnostic$plots
    )
    drawLmeDiagnosticPlotsDnaEpico(diagnostic$plots, files, config)
    list(plots = diagnostic$plots, files = files, lambda = diagnostic$lambda)
}

buildLmePhenotypeDiagnosticsDnaEpico <- function(
    summaryData, phenotype, diagnosticMean, config
) {
    if (is.null(summaryData) || !nrow(summaryData)) {
    return(NULL)
    }
    summaryData$FDR <- adjustPvaluesByTermMethylationModels(
    pValues = summaryData$P.value,
    terms = summaryData$Interaction.Term, method = config$padjmethod
    )
    terms <- unique(as.character(summaryData$Interaction.Term))
    terms <- terms[!is.na(terms)]
    multiple <- length(terms) > 1L
    results <- lapply(seq_along(terms), function(index) {
    buildLmeTermDiagnosticDnaEpico(
        summaryData, phenotype, terms[[index]], index, multiple,
        diagnosticMean, config
    )
    })
    names(results) <- terms
    results <- Filter(Negate(is.null), results)
    if (!length(results)) {
    return(NULL)
    }
    if (multiple) {
    return(list(
        plots = lapply(results, `[[`, "plots"),
        files = lapply(results, `[[`, "files"),
        lambda = vapply(results, `[[`, numeric(1), "lambda")
    ))
    }
    list(
    plots = results[[1L]]$plots, files = results[[1L]]$files,
    lambda = results[[1L]]$lambda
    )
}

normalizeLmeDiagnosticConfigDnaEpico <- function(config) {
    config$logPath <- resolveLogPathMinfiEwasWater(
    logs = config$logs, log_dir = config$log_dir,
    log_file = config$log_file
    )
    config$fdrThreshold <- validateProbabilityDnaEpico(
    config$fdrThreshold, "fdrThreshold"
    )
    config$padjmethod <- validatePAdjustmentMethodMethylationModels(
    config$padjmethod
    )
    config
}

#' Plot longitudinal mixed-effects model diagnostics
#'
#' @param modelSummaries Object returned by `summarizeMethylationLMEModels()`.
#' @param preparedData Object returned by `prepareMethylationLMEData()`.
#' @param fdrThreshold Numeric. False-discovery-rate threshold used to highlight
#'   CpGs in the diagnostic plots.
#' @param padjmethod Character. P-value adjustment method passed to
#'   `stats::p.adjust()`.
#' @param outputDir Character or `NULL`. Directory used for TIFF files. When
#'   `NULL`, plots are returned in memory only.
#' @param plotWidth Integer. TIFF width in pixels when plots are written to
#' disk.
#' @param plotHeight Integer. TIFF height in pixels when plots are written to
#'   disk.
#' @param plotDPI Integer. TIFF resolution in DPI when plots are written to
#'   disk.
#' @param display Logical. If `TRUE`, draw plots on the active graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationLME_diagnostic_plots'`
#'   containing the generated `ggplot2` objects, genomic inflation factors, and
#'   any saved TIFF file paths.
#'
#' @description
#' Create Q-Q and standard-error diagnostic plots from the mixed-effects summary
#' tables returned by `summarizeMethylationLMEModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' diagnostic_plots <- plotMethylationLMEDiagnostics(
#'   modelSummaries = ex$modelSummaries,
#'   preparedData = ex$preparedData,
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(diagnostic_plots$plots)
#'
#' @export
plotMethylationLMEDiagnostics <- function(
    modelSummaries, preparedData,
    fdrThreshold = 0.05, padjmethod = "fdr", outputDir = NULL,
    plotWidth = 2000L, plotHeight = 1000L, plotDPI = 150L, display = FALSE,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
    config <- normalizeLmeDiagnosticConfigDnaEpico(
    as.list(environment(), all.names = TRUE)
    )
    summaries <- if (is.null(modelSummaries$diagnosticSummaries)) {
    modelSummaries$summaries
    } else {
    modelSummaries$diagnosticSummaries
    }
    diagnostic_mean <- diagnosticMeanMethylationModels(preparedData)
    results <- lapply(names(summaries), function(phenotype) {
    buildLmePhenotypeDiagnosticsDnaEpico(
        summaries[[phenotype]], phenotype, diagnostic_mean, config
    )
    })
    names(results) <- names(summaries)
    results <- Filter(Negate(is.null), results)
    plots <- lapply(results, `[[`, "plots")
    files <- lapply(results, `[[`, "files")
    inflation <- lapply(results, `[[`, "lambda")
    emitLogMinfiEwasWater(c(
    "============================================================",
    paste("Diagnostic plots generated for phenotypes:", length(plots)),
    if (is.null(config$outputDir)) {
        "Diagnostic plots were returned in memory only."
    } else {
        paste("Diagnostic plots saved to:    ", config$outputDir)
    },
    "============================================================"
    ), verbose = config$verbose, log_path = config$logPath)
    structure(list(
    plots = plots, inflationFactors = inflation, files = files
    ), class = "dnaEPICO_methylationLME_diagnostic_plots")
}

buildAnnotationOmnibusTablesMethylationLME <- function(modelSummaries) {
    omnibus_tables <- modelSummaries$omnibusTests
    if (!is.list(omnibus_tables) || length(omnibus_tables) == 0L) {
    return(list())
    }
    interaction_term <- modelSummaries$settings$interactionTerm

    tables <- lapply(names(omnibus_tables), function(phenotype) {
    table <- omnibus_tables[[phenotype]]
    required <- c(
        "CpG", "Omnibus.F.Value", "Omnibus.Num.DF",
        "Omnibus.Den.DF", "Omnibus.P.Value",
        "Omnibus.Adjusted.P.Value", "Omnibus.Method"
    )
    if (!is.data.frame(table) || nrow(table) == 0L ||
        !all(required %in% colnames(table))) {
        return(NULL)
    }

    prefix_parts <- phenotype
    if (!is.null(interaction_term) && nzchar(interaction_term)) {
        prefix_parts <- c(prefix_parts, interaction_term)
    }
    prefix <- paste(gsub("`", "", prefix_parts, fixed = TRUE),
        collapse = "_"
    )
    result <- table[, required, drop = FALSE]
    colnames(result) <- c(
        "CpG", paste0(prefix, "_Omnibus_F.Value"),
        paste0(prefix, "_Omnibus_Num.DF"),
        paste0(prefix, "_Omnibus_Den.DF"),
        paste0(prefix, "_Omnibus_P.Value"),
        paste0(prefix, "_Omnibus_Adjusted.P.Value"),
        paste0(prefix, "_Omnibus_Method")
    )
    result
    })
    names(tables) <- names(omnibus_tables)
    Filter(Negate(is.null), tables)
}

cleanLmeModelSummaryDnaEpico <- function(data, modelName) {
    required <- c("CpG", "Interaction.Term", "P.value")
    if (is.null(data) || !nrow(data) || !all(required %in% colnames(data))) {
    return(NULL)
    }
    split_data <- split(data, data$Interaction.Term)
    tables <- lapply(names(split_data), function(term) {
    output <- split_data[[term]][, c("CpG", "P.value"), drop = FALSE]
    clean_term <- gsub("`", "", term, fixed = TRUE)
    suffix <- clean_term
    prefix <- paste0(modelName, ".")
    if (startsWith(clean_term, prefix)) {
        suffix <- sub(
        paste0("^", escapeRegexMethylationGLM(modelName), "\\."),
        "", clean_term
        )
    }
    colnames(output)[2L] <- paste0(modelName, "_", suffix, "_P.Value")
    output
    })
    Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), tables)
}

mergeLmeSummaryTablesDnaEpico <- function(summaryList, modelSummaries) {
    cleaned <- lapply(names(summaryList), function(name) {
    cleanLmeModelSummaryDnaEpico(summaryList[[name]], name)
    })
    cleaned <- Filter(Negate(is.null), cleaned)
    merged <- if (!length(cleaned)) {
    data.frame(CpG = character(0))
    } else {
    Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), cleaned)
    }
    omnibus <- buildAnnotationOmnibusTablesMethylationLME(modelSummaries)
    for (table in omnibus) {
    merged <- merge(merged, table, by = "CpG", all = TRUE)
    }
    p_columns <- grep("P\\.Value$|P\\.value$", names(merged), value = TRUE)
    if (length(p_columns) && nrow(merged)) {
    has_p <- apply(merged[, p_columns, drop = FALSE], 1L, function(values) {
        any(is.finite(as.numeric(values)))
    })
    merged <- merged[has_p, , drop = FALSE]
    }
    messages <- buildAnnotationModelMessagesDnaEpico(
    modelSummaries$modelMessages
    )
    if (nrow(messages)) {
    merged <- merge(merged, messages, by = "CpG", all.x = TRUE)
    }
    merged
}

appendLmeGencodeDnaEpico <- function(data, annotationObject, enabled) {
    if (!enabled) {
    return(list(data = data, columns = character(0), result = NULL))
    }
    validateGencodeHubArrayCompatibilityDnaEpico(annotationObject)
    resource <- resolveGencodeHubResourceDnaEpico()
    result <- appendGencodeHubAnnotationDnaEpico(data, resource)
    result$resource <- resource
    result$dictionaryRows <- buildGencodeHubDictionaryDnaEpico(result$release)
    result$metadataRows <- buildGencodeHubMetadataDnaEpico(
    resource, result$counts
    )
    list(
    data = result$data, columns = result$annotationColumns,
    result = result
    )
}

annotateLmeSummaryDataDnaEpico <- function(
    merged, annotationData, annotationCols, annotationObject,
    gencodeHub, modelNames
) {
    requested <- unique(c(
    annotationCols, if (gencodeHub) c("chr", "pos") else character(0)
    ))
    available <- intersect(requested, colnames(annotationData))
    missing <- setdiff(annotationCols, colnames(annotationData))
    annotated <- merge(
    merged, annotationData[, c("CpG", available), drop = FALSE],
    by = "CpG", all.x = TRUE
    )
    if ("CpG" %in% colnames(annotated)) {
    colnames(annotated)[colnames(annotated) == "CpG"] <- "IlmnID"
    }
    gencode <- appendLmeGencodeDnaEpico(
    annotated, annotationObject, gencodeHub
    )
    used <- c(available, gencode$columns)
    data <- orderAnnotatedModelColumnsDnaEpico(
    data = gencode$data, annotationCols = used, modelNames = modelNames
    )
    list(
    data = data, used = used, missing = missing,
    gencode = gencode$result
    )
}

logLmeAnnotationDnaEpico <- function(result, gencodeHub, verbose, logPath) {
    gencode_line <- if (gencodeHub) {
    paste0(
        "GENCODEHub resource:          ",
        result$gencode$resource$annotationHubId,
        " (GENCODE release ", result$gencode$release, "; ",
        "direct=", result$gencode$counts[["annotated"]], ", ",
        "nearest=", result$gencode$counts[["non_annotated"]], ", ",
        "unassigned=", result$gencode$counts[["unassigned"]], ")"
    )
    } else {
    "GENCODEHub resource:          disabled"
    }
    missing_line <- if (!length(result$missing)) {
    "Missing annotation columns:   none"
    } else {
    paste(
        "Missing annotation columns:   ",
        paste(result$missing, collapse = ", ")
    )
    }
    emitLogMinfiEwasWater(c(
    "============================================================",
    paste("Annotated CpG rows:          ", nrow(result$data)),
    paste(
        "Annotation columns used:      ",
        paste(result$used, collapse = ", ")
    ),
    gencode_line, missing_line,
    "============================================================"
    ), verbose = verbose, log_path = logPath)
}

#' Annotate longitudinal mixed-effects summary tables with array annotation
#' metadata
#'
#' @param modelSummaries Object returned by
#'   `summarizeMethylationLMEModels()` or a named list of summary data
#'   frames.
#' @param annotationObject Character package/object name, annotation data frame,
#'   or annotation object understood by `minfi::getAnnotation()`.
#' @param annotationCols Character vector or comma-separated string of
#' annotation
#'   columns to append.
#' @param gencodeHub Logical. If `TRUE`, retrieve the package-managed GENCODE
#'   gene resource from AnnotationHub and append direct gene-body and nearest
#'   transcription-start-site annotations. This requires GRCh38 coordinates.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationLME_annotation'`
#'   containing the annotated summary table and any requested annotation columns
#'   that were unavailable in the chosen annotation object.
#'
#' @description
#' Merge phenotype-specific longitudinal summary tables with probe annotation
#' metadata and return a single annotated result table.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' annotation_data <- annotateMethylationLMESummaries(
#'   modelSummaries = ex$modelSummaries,
#'   annotationObject = ex$annotationData,
#'   annotationCols = "Name,chr,pos",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(annotation_data)
#'
#' @export
annotateMethylationLMESummaries <- function(
    modelSummaries, annotationObject,
    annotationCols = c(
    "Name", "chr", "pos", "UCSC_RefGene_Group",
    "UCSC_RefGene_Name", "Relation_to_Island", "GencodeV41_Group"
    ),
    gencodeHub = FALSE,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
    gencodeHub <- validateLogicalScalarDnaEpico(gencodeHub, "gencodeHub")
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    summary_list <- modelSummaries
    if (!is.null(modelSummaries$diagnosticSummaries)) {
    summary_list <- modelSummaries$diagnosticSummaries
    } else if (!is.null(modelSummaries$summaries)) {
    summary_list <- modelSummaries$summaries
    }
    phenotypes <- modelSummaries$phenotypes
    if (is.null(phenotypes)) {
    phenotypes <- names(summary_list)
    }
    annotation_cols <- splitOptionMinfiEwasWater(annotationCols, sep = ",")
    annotation_data <- coerceAnnotationDataMethylationGLM(annotationObject)
    merged <- mergeLmeSummaryTablesDnaEpico(summary_list, modelSummaries)
    result <- annotateLmeSummaryDataDnaEpico(
    merged, annotation_data, annotation_cols, annotationObject,
    gencodeHub, phenotypes
    )
    logLmeAnnotationDnaEpico(result, gencodeHub, verbose, log_path)
    structure(list(
    data = result$data, fitFailures = modelSummaries$fitFailures,
    modelMessages = modelSummaries$modelMessages,
    annotationColumnsUsed = result$used,
    missingAnnotationCols = result$missing,
    gencodeHub = result$gencode
    ), class = "dnaEPICO_methylationLME_annotation")
}

validateLmeReportAssetsDirDnaEpico <- function(reportAssetsDir) {
    if (!is.null(reportAssetsDir) &&
    (!is.character(reportAssetsDir) || length(reportAssetsDir) != 1L ||
        is.na(reportAssetsDir) || !nzchar(reportAssetsDir))) {
    stop("reportAssetsDir must be NULL or one non-empty directory path.",
        call. = FALSE
    )
    }
    invisible(TRUE)
}

saveLmePhenotypeSummariesDnaEpico <- function(
    modelResults, modelSummaries, outputDir
) {
    files <- stats::setNames(
    character(length(modelSummaries$summaries)),
    names(modelSummaries$summaries)
    )
    phenotypes <- modelResults$phenotypes
    if (is.null(phenotypes)) {
    phenotypes <- names(modelResults$summaryCache)
    }
    for (phenotype in phenotypes) {
    file <- file.path(outputDir, paste0(phenotype, "SummaryLME.rds"))
    artifact <- modelResults$phenotypeSummaries[[phenotype]]
    if (is.null(artifact)) {
        stop(
        "A complete compact phenotype summary is unavailable for ",
        phenotype, ".",
        call. = FALSE
        )
    }
    existing <- loadPhenotypeSummaryMethylationModels(
        path = file, expectedSignature = artifact$signature
    )
    if (is.null(existing$object)) {
        savePhenotypeSummaryMethylationModels(artifact, file)
    }
    files[[phenotype]] <- file
    }
    files
}

writeLmeTextSummariesDnaEpico <- function(
    summaries, outputDir, enabled
) {
    files <- list()
    if (!isTRUE(enabled)) {
    return(files)
    }
    dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
    for (phenotype in names(summaries)) {
    data <- summaries[[phenotype]]
    if (is.null(data) || !nrow(data)) {
        next
    }
    if ("P.value" %in% colnames(data)) {
        data <- data[order(data$P.value), , drop = FALSE]
    }
    file <- file.path(outputDir, paste0(phenotype, "SummaryLME.txt"))
    utils::write.table(
        data,
        file = file, sep = "\t", row.names = FALSE, quote = FALSE
    )
    files[[phenotype]] <- file
    }
    files
}

writeSignificantLmeTablesDnaEpico <- function(
    significant, outputDir, enabled
) {
    files <- list()
    if (!isTRUE(enabled) || is.null(significant)) {
    return(files)
    }
    dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)
    for (phenotype in names(significant)) {
    hits <- significant[[phenotype]]
    if (!length(hits)) {
        next
    }
    phenotype_dir <- file.path(outputDir, phenotype)
    dir.create(phenotype_dir, recursive = TRUE, showWarnings = FALSE)
    files[[phenotype]] <- vapply(names(hits), function(cpg) {
        cpg_dir <- file.path(phenotype_dir, cpg)
        dir.create(cpg_dir, recursive = TRUE, showWarnings = FALSE)
        file <- file.path(cpg_dir, paste0(cpg, ".txt"))
        utils::write.table(hits[[cpg]], file = file, sep = "\t", quote = FALSE)
        file
    }, character(1))
    }
    files
}

addLmeGencodeWorkbookRowsDnaEpico <- function(
    dictionary, metadata, annotatedResults
) {
    gencode <- annotatedResults$gencodeHub
    if (is.null(gencode)) {
    return(list(dictionary = dictionary, metadata = metadata))
    }
    gencode_dictionary <- gencode$dictionaryRows
    if (is.data.frame(gencode_dictionary) && nrow(gencode_dictionary)) {
    dictionary <- dictionary[
        !(dictionary$Column %in% gencode_dictionary$Column), ,
        drop = FALSE
    ]
    dictionary <- rbind(dictionary, gencode_dictionary)
    }
    gencode_metadata <- gencode$metadataRows
    if (is.data.frame(gencode_metadata) && nrow(gencode_metadata)) {
    metadata <- metadata[
        !startsWith(metadata$Key, "gencode."), ,
        drop = FALSE
    ]
    metadata <- rbind(metadata, gencode_metadata)
    }
    list(dictionary = dictionary, metadata = metadata)
}

addLmeVennWorkbookRowsDnaEpico <- function(
    dictionary, metadata, vennDResults
) {
    if (is.null(vennDResults)) {
    return(list(
        dictionary = dictionary, metadata = metadata, sheets = list()
    ))
    }
    if (!inherits(vennDResults, "dnaEPICO_vennD_plots")) {
    stop("vennDResults must be a model-level Venn result.", call. = FALSE)
    }
    if (is.data.frame(vennDResults$metadataRows) &&
    nrow(vennDResults$metadataRows)) {
    metadata <- rbind(metadata, vennDResults$metadataRows)
    }
    if (is.data.frame(vennDResults$dictionaryRows) &&
    nrow(vennDResults$dictionaryRows)) {
    dictionary <- rbind(dictionary, vennDResults$dictionaryRows)
    }
    list(
    dictionary = dictionary, metadata = metadata,
    sheets = vennDResults$sheets
    )
}

writeLmeAnnotatedWorkbookDnaEpico <- function(
    modelResults, modelSummaries, annotatedResults, annotatedLMEOut,
    reportAssetsDir, vennDResults
) {
    annotated <- if (is.null(annotatedResults$data)) {
    annotatedResults
    } else {
    annotatedResults$data
    }
    report_table <- sortReportTableDnaEpico(annotated)
    annotated <- report_table$data
    file <- file.path(annotatedLMEOut, "annotatedLME.xlsx")
    dictionary <- buildAnnotatedWorkbookDictionaryMethylationGLM(
    columns = colnames(annotated), modelDescription = "Pvalue from LME model",
    formulaText = modelResults$formulas, modelLabel = "LME",
    responseLabel = inferMethylationValueLabelMethylationGLM(modelResults)
    )
    metadata <- buildModelWorkbookMetadataDnaEpico(
    modelResults, modelSummaries, annotatedResults,
    analysis = "lme"
    )
    rows <- addLmeGencodeWorkbookRowsDnaEpico(
    dictionary, metadata, annotatedResults
    )
    rows <- addLmeVennWorkbookRowsDnaEpico(
    rows$dictionary, rows$metadata, vennDResults
    )
    writeAnnotatedWorkbookMethylationGLM(
    annotated_df = annotated, file = file, resultSheet = "annotatedLME",
    dictionary = rows$dictionary, metadata = rows$metadata,
    extraSheets = rows$sheets
    )
    sidecar <- list(
    table = NULL, metadata = NULL, dictionary = NULL,
    workbookMetadata = NULL
    )
    if (!is.null(reportAssetsDir)) {
    sidecar <- writeReportTableSidecarDnaEpico(
        tableData = annotated, workbookFile = file,
        sidecarDir = reportAssetsDir, sheet = "annotatedLME",
        idColumn = report_table$idColumn, dictionary = rows$dictionary,
        workbookMetadata = rows$metadata
    )
    }
    list(file = file, sidecar = sidecar, sheets = names(rows$sheets))
}

logLmeOutputFilesDnaEpico <- function(
    summaryFiles, textFiles, significantFiles, workbook,
    fitFailures, saveTxtSummaries, saveSignificantInteractions,
    verbose, logPath
) {
    sidecar_line <- if (is.null(workbook$sidecar$table)) {
    "Report table sidecar:          not requested"
    } else {
    paste("Report table sidecar:         ", workbook$sidecar$table)
    }
    text_count <- if (isTRUE(saveTxtSummaries)) length(textFiles) else 0L
    significant_count <- if (isTRUE(saveSignificantInteractions)) {
    sum(vapply(significantFiles, length, integer(1)))
    } else {
    0L
    }
    emitLogMinfiEwasWater(c(
    "============================================================",
    "Full model files written:       0",
    paste("Compact phenotype summaries: ", length(summaryFiles)),
    paste("Annotated results file:      ", workbook$file),
    sidecar_line,
    paste("Hard model errors retained:  ", nrow(fitFailures)),
    paste("Summary text files written:   ", text_count),
    paste("Significant interaction files:", significant_count),
    "============================================================"
    ), verbose = verbose, log_path = logPath)
}

#' Write optional disk outputs for longitudinal mixed-effects analyses
#'
#' @param modelResults Object returned by `fitMethylationLMEModels()`.
#' @param modelSummaries Object returned by `summarizeMethylationLMEModels()`.
#' @param annotatedResults Object returned by
#'   `annotateMethylationLMESummaries()` or a compatible data frame.
#' @param significantInteractions Object returned by
#'   `collectSignificantInteractionsMethylationLME()` or `NULL`.
#' @param outputRData Character. Directory used for complete compact phenotype
#'   summaries.
#' @param summaryTxtDir Character. Directory used for tab-delimited summary
#'   tables.
#' @param significantInteractionDir Character. Directory used for significant
#'   interaction coefficient tables.
#' @param annotatedLMEOut Character. Directory used for the annotated summary
#'   XLSX workbook.
#' @param reportAssetsDir Character or `NULL`. Report results directory used for
#'   the compressed TSV table and compact metadata sidecars.
#' @param vennDResults Optional model-level Venn result. Its configuration
#'   metadata and threshold tables are added before the workbook dictionary.
#' @param saveTxtSummaries Logical. If `TRUE`, write tab-delimited summary
#' tables.
#' @param saveSignificantInteractions Logical. If `TRUE`, write significant
#'   interaction coefficient tables.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationLME_paths'`
#'   containing the paths of the files written to disk, including the annotated
#'   workbook and any requested report-table sidecars.
#'
#' @description
#' Write compact phenotype summaries, optional text and significant-interaction
#' tables, and annotated results from the longitudinal mixed-effects workflow.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' annotation_data <- annotateMethylationLMESummaries(
#'   modelSummaries = ex$modelSummaries,
#'   annotationObject = ex$annotationData,
#'   annotationCols = "Name,chr,pos",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' significant_hits <- collectSignificantInteractionsMethylationLME(
#'   modelResults = ex$modelResults,
#'   pvalThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' output_paths <- writeMethylationLMEOutputs(
#'   modelResults = ex$modelResults,
#'   modelSummaries = ex$modelSummaries,
#'   annotatedResults = annotation_data,
#'   significantInteractions = significant_hits,
#'   outputRData = file.path(ex$tempDir, "models"),
#'   summaryTxtDir = file.path(ex$tempDir, "summary"),
#'   significantInteractionDir = file.path(ex$tempDir, "significant"),
#'   annotatedLMEOut = file.path(ex$tempDir, "annotated"),
#'   saveTxtSummaries = TRUE,
#'   saveSignificantInteractions = TRUE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(output_paths)
#'
#' @export
writeMethylationLMEOutputs <- function(
    modelResults, modelSummaries,
    annotatedResults, significantInteractions = NULL, outputRData,
    summaryTxtDir, significantInteractionDir, annotatedLMEOut,
    reportAssetsDir = NULL, vennDResults = NULL, saveTxtSummaries = TRUE,
    saveSignificantInteractions = FALSE,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
    validateLmeReportAssetsDirDnaEpico(reportAssetsDir)
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    dir.create(outputRData, recursive = TRUE, showWarnings = FALSE)
    dir.create(annotatedLMEOut, recursive = TRUE, showWarnings = FALSE)
    summary_files <- saveLmePhenotypeSummariesDnaEpico(
    modelResults, modelSummaries, outputRData
    )
    text_files <- writeLmeTextSummariesDnaEpico(
    modelSummaries$summaries, summaryTxtDir, saveTxtSummaries
    )
    significant_files <- writeSignificantLmeTablesDnaEpico(
    significantInteractions, significantInteractionDir,
    saveSignificantInteractions
    )
    workbook <- writeLmeAnnotatedWorkbookDnaEpico(
    modelResults, modelSummaries, annotatedResults, annotatedLMEOut,
    reportAssetsDir, vennDResults
    )
    fit_failures <- modelSummaries$fitFailures
    if (is.null(fit_failures)) {
    fit_failures <- modelResults$fitFailures
    }
    logLmeOutputFilesDnaEpico(
    summary_files, text_files, significant_files, workbook,
    fit_failures, saveTxtSummaries, saveSignificantInteractions,
    verbose, log_path
    )
    structure(list(
    modelFiles = character(0), summaryFiles = summary_files,
    summaryTxtFiles = text_files,
    significantInteractionFiles = significant_files,
    annotatedLME = workbook$file,
    annotatedLMEText = workbook$sidecar$table,
    annotatedLMEReportMetadata = workbook$sidecar$metadata,
    annotatedLMEDictionary = workbook$sidecar$dictionary,
    annotatedLMEMetadata = workbook$sidecar$workbookMetadata,
    vennDSheets = workbook$sheets
    ), class = "dnaEPICO_methylationLME_paths")
}
