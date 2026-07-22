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
ensurePersonColumnMethylationLME <- function(
    data, personVar = "person",
    sidVar = "SID"
) {
    person_created <- FALSE
    mapping_preview <- NULL

    if (!(personVar %in% colnames(data))) {
        if (!(sidVar %in% colnames(data))) {
            stop("Column '", personVar,
                "' was not found and cannot be created because '",
                sidVar, "' is missing.",
                call. = FALSE
            )
        }

        sid_values <- validateSampleIdentifiersDnaEpico(
            data[[sidVar]],
            paste0("inputPheno$", sidVar)
        )
        if (any(!grepl("[AB]$", sid_values))) {
            stop("Cannot safely derive '", personVar, "' from '",
                sidVar, "': every SID must end in A or B. Supply an explicit subject identifier column.",
                call. = FALSE
            )
        }
        person_values <- sub("[AB]$", "", sid_values)
        visit_values <- sub("^.*([AB])$", "\\1", sid_values)
        if (any(!nzchar(person_values))) {
            stop("Derived subject identifiers cannot be empty.",
                call. = FALSE
            )
        }
        if (anyDuplicated(paste(person_values, visit_values,
            sep = "\r"
        ))) {
            stop("SID contains duplicate subject/visit identifiers.",
                call. = FALSE
            )
        }
        data[[personVar]] <- person_values
        person_created <- TRUE
        mapping_preview <- utils::head(
            data[order(
                data[[personVar]],
                data[[sidVar]]
            ), c(sidVar, personVar), drop = FALSE],
            20L
        )
    }

    person_values <- as.character(data[[personVar]])
    if (anyNA(person_values) || any(!nzchar(trimws(person_values)))) {
        stop("personVar contains missing or blank subject identifiers.",
            call. = FALSE
        )
    }
    person_counts <- table(data[[personVar]], useNA = "ifany")
    if (length(person_counts) < 2L) {
        stop("At least two subjects are required for mixed-effects modeling.",
            call. = FALSE
        )
    }
    if (max(person_counts) < 2L) {
        stop("At least one subject must have repeated observations for mixed-effects modeling.",
            call. = FALSE
        )
    }

    list(
        data = data, personCreated = person_created,
            mappingPreview = mapping_preview,
        personCounts = person_counts
    )
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
buildFormulaMethylationLME <- function(
    phenotype, personVar,
    covariates = character(0), interactionTerm = NULL, includeRandomTerm = TRUE,
    responseVar = "beta"
) {
    if (!is.null(interactionTerm) && nzchar(interactionTerm) &&
        identical(interactionTerm, phenotype)) {
        stop("interactionTerm must differ from the phenotype being modelled.",
            call. = FALSE
        )
    }
    quoted_phenotype <- quoteNamesMethylationGLM(phenotype)
    quoted_person <- quoteNamesMethylationGLM(personVar)
    fixed_terms <- unique(c(covariates))

    if (!is.null(interactionTerm) && nzchar(interactionTerm)) {
        quoted_interaction <- quoteNamesMethylationGLM(interactionTerm)
        interaction_part <- paste(quoted_phenotype, quoted_interaction,
            sep = " * "
        )
        fixed_terms <- setdiff(fixed_terms, c(
            interactionTerm,
            phenotype
        ))
        fixed_formula_terms <- c(interaction_part,
            quoteNamesMethylationGLM(fixed_terms))
    } else {
        fixed_terms <- setdiff(fixed_terms, phenotype)
        fixed_formula_terms <- c(quoted_phenotype,
            quoteNamesMethylationGLM(fixed_terms))
    }

    fixed_formula_terms <-
        unique(fixed_formula_terms[nzchar(fixed_formula_terms)])
    if (length(fixed_formula_terms) == 0L) {
        stop("At least one fixed-effect term is required.", call. = FALSE)
    }

    formula_text <- paste(responseVar, "~", paste(fixed_formula_terms,
        collapse = " + "
    ))

    if (isTRUE(includeRandomTerm)) {
        formula_text <- paste(
            formula_text, "+ (1 |", quoted_person,
            ")"
        )
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
            "omnibusTest is currently available only for the lmerTest/lme4 engine, not nlme.",
            call. = FALSE
        )
    }
    if (isTRUE(omnibus_test) &&
        identical(omnibus_ddf, "Kenward-Roger") &&
        !requireNamespace("pbkrtest", quietly = TRUE)) {
        stop(
            "omnibusDdf = 'Kenward-Roger' requires the suggested package 'pbkrtest'.",
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
        stop(
            "Could not identify one fixed-effect model term for omnibus testing: ",
            target_label, ".",
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
        status = status, reason = reason, warning = NA_character_,
        rhs = 0, eps = sqrt(.Machine$double.eps)
    )
}

computeOmnibusTestMethylationLME <- function(
    fit, coefficientTerms, omnibusTerm,
    omnibusDdf = "Satterthwaite"
) {
    result <- emptyOmnibusResultMethylationLME(
        term = omnibusTerm,
        method = omnibusDdf
    )

    tryCatch({
        fixed_effects <- lme4::fixef(fit)
        fixed_names <- names(fixed_effects)
        mapped_terms <- coefficientTerms[fixed_names]
        selected <- which(!is.na(mapped_terms) & mapped_terms == omnibusTerm)
        if (length(selected) == 0L) {
            stop(
                "The omnibus model term has no estimable fixed-effect coefficients.",
                call. = FALSE
            )
        }

        contrast <- matrix(
            0,
            nrow = length(selected), ncol = length(fixed_effects),
            dimnames = list(fixed_names[selected], fixed_names)
        )
        contrast[cbind(seq_along(selected), selected)] <- 1

        warning_state <- new.env(parent = emptyenv())
        warning_state$messages <- character(0)
        test <- withCallingHandlers(
            lmerTest::contestMD(
                model = fit, L = contrast, rhs = 0,
                ddf = omnibusDdf, joint = TRUE,
                eps = sqrt(.Machine$double.eps)
            ),
            warning = function(condition) {
                warning_state$messages <- c(
                    warning_state$messages,
                    conditionMessage(condition)
                )
                invokeRestart("muffleWarning")
            }
        )
        required <- c("F value", "NumDF", "DenDF", "Pr(>F)")
        if (!all(required %in% colnames(test))) {
            stop(
                "The omnibus test did not return the expected F-test statistics.",
                call. = FALSE
            )
        }
        values <- unlist(test[1L, required, drop = TRUE], use.names = FALSE)
        if (any(!is.finite(values))) {
            stop(
                "The omnibus test returned missing or non-finite statistics.",
                call. = FALSE
            )
        }

        result$fValue <- unname(test[["F value"]][[1L]])
        result$numeratorDf <- unname(test[["NumDF"]][[1L]])
        result$denominatorDf <- unname(test[["DenDF"]][[1L]])
        result$pValue <- unname(test[["Pr(>F)"]][[1L]])
        result$status <- "tested"
        result$reason <- NA_character_
        warning_messages <- unique(
            warning_state$messages[nzchar(warning_state$messages)]
        )
        if (length(warning_messages) > 0L) {
            result$warning <- paste(warning_messages, collapse = " | ")
        }
        result
    }, error = function(error) {
        result$reason <- conditionMessage(error)
        result
    })
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

    stop("correlationVar must be numeric or contain numeric text; categorical values cannot define AR1/CAR1 spacing.",
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
        stop("correlationVar contains missing or non-finite values after conversion.",
            call. = FALSE
        )
    }
    person_time <- interaction(person_values, time_values, drop = TRUE)
    if (anyDuplicated(person_time)) {
        stop("correlationVar must be unique within each subject for AR1/CAR1 models.",
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
        stop("An internal correlation time variable is required for AR1/CAR1 models.",
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

fitStatusValuesMethylationLME <- function(modelObj) {
    if (inherits(modelObj, "dnaEPICO_methylationLME_fit_error")) {
        return(list(
            status = if (is.null(modelObj$status)) "failed" else as.character(modelObj$status),
            singular = NA,
                converged = if (is.null(modelObj$converged)) NA else as.logical(modelObj$converged)[[1L]],
            convergenceMessage = if (is.null(modelObj$convergenceMessage)) {
                NA_character_
            } else {
                as.character(modelObj$convergenceMessage)[[1L]]
            }, warning = NA_character_
        ))
    }

    fit_status <- modelObj$fitStatus
    singular_fit <- if (is.null(fit_status$singular)) {
        NA
    } else {
        as.logical(fit_status$singular)[[1L]]
    }
    converged <- if (is.null(fit_status$converged)) {
        length(fit_status$convergenceMessages) == 0L
    } else {
        isTRUE(fit_status$converged)
    }
    convergence_messages <- unique(as.character(fit_status$convergenceMessages))
    convergence_messages <- convergence_messages[!is.na(convergence_messages) &
        nzchar(convergence_messages)]
    fit_warnings <- unique(as.character(fit_status$warnings))
    fit_warnings <- fit_warnings[!is.na(fit_warnings) & nzchar(fit_warnings)]

    list(status = if (!isTRUE(converged)) {
        "fitted_not_converged"
    } else if (isTRUE(singular_fit)) {
        "fitted_singular"
    } else if (length(fit_warnings) > 0L) {
        "fitted_with_warning"
    } else {
        "fitted"
    }, singular = singular_fit, converged = converged,
        convergenceMessage = if (length(convergence_messages) >
        0L) {
        paste(convergence_messages, collapse = " | ")
    } else {
        NA_character_
    }, warning = if (length(fit_warnings) > 0L) {
        paste(fit_warnings, collapse = " | ")
    } else {
        NA_character_
    })
}

collectFitDiagnosticsMethylationLME <- function(fits) {
    rows <- list()
    row_index <- 1L
    for (phenotype in names(fits)) {
        fit_group <- fits[[phenotype]]
        for (cpg in names(fit_group)) {
            fit_object <- fit_group[[cpg]]
            values <- fitStatusValuesMethylationLME(fit_object)
            inference_included <- startsWith(values$status, "fitted")
            omnibus_failure <- !is.null(fit_object$omnibus) &&
                !identical(fit_object$omnibus$status, "tested")
            if (isTRUE(omnibus_failure)) {
                inference_included <- FALSE
            }
            rows[[row_index]] <- data.frame(
                Phenotype = phenotype,
                CpG = cpg, Fit.Status = values$status,
                    Singular.Fit = values$singular,
                Converged = values$converged,
                    Convergence.Message = values$convergenceMessage,
                Fit.Warning = values$warning,
                    Inference.Included = inference_included,
                Exclusion.Reason = if (!inference_included &&
                    !is.null(fit_object$error)) {
                    as.character(fit_object$error)
                } else if (isTRUE(omnibus_failure)) {
                    paste0(
                        "omnibus test not estimable: ",
                        fit_object$omnibus$reason
                    )
                } else {
                    NA_character_
                }, stringsAsFactors = FALSE, check.names = FALSE
            )
            row_index <- row_index + 1L
        }
    }

    if (length(rows) == 0L) {
        return(data.frame(
            Phenotype = character(0), CpG = character(0),
            Fit.Status = character(0), Singular.Fit = logical(0),
            Converged = logical(0), Convergence.Message = character(0),
            Fit.Warning = character(0), Inference.Included = logical(0),
            Exclusion.Reason = character(0), stringsAsFactors = FALSE,
            check.names = FALSE
        ))
    }

    diagnostics <- do.call(rbind, rows)
    rownames(diagnostics) <- NULL
    diagnostics
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
    if (is.null(modelObj) || inherits(modelObj,
        "dnaEPICO_methylationLME_fit_error")) {
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

    fit_status <- fitStatusValuesMethylationLME(modelObj)
    omnibus_failure <- !is.null(modelObj$omnibus) &&
        !identical(modelObj$omnibus$status, "tested")

    do.call(rbind, lapply(matched_terms, function(term) {
        coef_row <- coef_table[term, ]
        data.frame(
            CpG = cpg, Interaction.Term = term,
                Estimate = unname(coef_row["Estimate"]),
            Std.Error = unname(coef_row["Std. Error"]),
                t.value = unname(coef_row["t value"]),
            P.value = unname(coef_row["Pr(>|t|)"]),
                Fit.Status = fit_status$status,
            Singular.Fit = fit_status$singular,
                Converged = fit_status$converged,
            Convergence.Message = fit_status$convergenceMessage,
            Fit.Warning = fit_status$warning,
            Inference.Included = !omnibus_failure,
            Exclusion.Reason = if (isTRUE(omnibus_failure)) {
                paste0(
                    "omnibus test not estimable: ",
                    modelObj$omnibus$reason
                )
            } else {
                NA_character_
            }, stringsAsFactors = FALSE,
            row.names = NULL
        )
    }))
}

applyFitQualityExclusionsMethylationLME <- function(
    summaryDf,
    excludeSingular = FALSE, excludeNonConverged = FALSE
) {
    summary_df <- summaryDf
    if (is.null(summary_df) || nrow(summary_df) == 0L) {
        return(summary_df)
    }
    reasons <- if ("Exclusion.Reason" %in% names(summary_df)) {
        as.character(summary_df$Exclusion.Reason)
    } else {
        rep(NA_character_, nrow(summary_df))
    }
    initially_excluded <- if ("Inference.Included" %in%
        names(summary_df)) {
        included <- as.logical(summary_df$Inference.Included)
        !is.na(included) & !included
    } else {
        rep(FALSE, nrow(summary_df))
    }
    singular <- rep(FALSE, nrow(summary_df))
    if (isTRUE(excludeSingular) && "Singular.Fit" %in% names(summary_df)) {
        singular_values <- as.logical(summary_df[["Singular.Fit"]])
        singular <- !is.na(singular_values) & singular_values
    }
    not_converged <- rep(FALSE, nrow(summary_df))
    if (isTRUE(excludeNonConverged) && "Converged" %in% names(summary_df)) {
        converged <- as.logical(summary_df[["Converged"]])
        not_converged <- !is.na(converged) & !converged
    }
    append_reason <- function(existing, added) {
        ifelse(
            is.na(existing) | !nzchar(existing),
            added,
            paste(existing, added, sep = "; ")
        )
    }
    reasons[singular] <- append_reason(
        reasons[singular],
        "singular random-effects fit"
    )
    reasons[not_converged] <- ifelse(is.na(reasons[not_converged]),
        "model did not converge", paste(reasons[not_converged],
            "model did not converge",
            sep = "; "
        )
    )
    excluded <- initially_excluded | singular | not_converged
    summary_df$Inference.Included <- !excluded
    summary_df$Exclusion.Reason <- reasons
    if (any(excluded) && "P.value" %in% names(summary_df)) {
        summary_df$P.value[excluded] <- NA_real_
    }
    summary_df
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

fitCpGModelMethylationLME <- function(
    cpg, cpgValues, modelData,
    formulaText, personVar, lmeEngine = "lme4", correlationStructure = "none",
    correlationTimeVar = NULL, responseVar = "beta",
    omnibusTest = FALSE, omnibusDdf = "Satterthwaite",
    omnibusTerm = NULL
) {
    tryCatch(
        {
            model_data <- modelData
            model_data[[responseVar]] <- as.numeric(cpgValues)
            observed_response <-
                model_data[[responseVar]][is.finite(model_data[[responseVar]])]
            if (length(unique(observed_response)) < 2L) {
                stop("The CpG response has no observed variation.",
                    call. = FALSE
                )
            }

            warning_state <- new.env(parent = emptyenv())
            warning_state$messages <- character(0)

            if (identical(lmeEngine, "nlme")) {
                fit <- withCallingHandlers(
                    nlme::lme(
                        fixed = stats::as.formula(formulaText),
                        random = stats::as.formula(paste("~ 1 |",
                            quoteNamesMethylationGLM(personVar))),
                        correlation = buildCorrelationMethylationLME(
                            correlationStructure = correlationStructure,
                            correlationTimeVar = correlationTimeVar,
                                personVar = personVar
                        ),
                        data = model_data, na.action = stats::na.exclude,
                        method = "REML",
                            control = nlme::lmeControl(returnObject = TRUE)
                    ),
                    warning = function(condition) {
                        warning_state$messages <- c(
                            warning_state$messages,
                            conditionMessage(condition)
                        )
                        invokeRestart("muffleWarning")
                    }
                )
                coef_table <-
                    coerceCoefficientTableMethylationLME(summary(fit)$tTable,
                    engine = lmeEngine
                )
                ranef_values <- nlme::ranef(fit)
                fixef_values <- nlme::fixef(fit)
                singular_fit <- NA
                fit_warnings <- warning_state$messages
                convergence_messages <- fit_warnings
            } else {
                fit <- withCallingHandlers(lmerTest::lmer(
                    formula = stats::as.formula(formulaText),
                    data = model_data, na.action = stats::na.exclude,
                    REML = TRUE
                ), warning = function(condition) {
                    warning_state$messages <- c(
                        warning_state$messages,
                        conditionMessage(condition)
                    )
                    invokeRestart("muffleWarning")
                })
                fit_warnings <- warning_state$messages
                coef_table <-
                    coerceCoefficientTableMethylationLME(summary(fit)$coefficients,
                    engine = lmeEngine
                )
                ranef_values <- lme4::ranef(fit)
                fixef_values <- lme4::fixef(fit)
                singular_fit <- lme4::isSingular(fit, tol = 1e-04)
                raw_optimizer_codes <- unlist(fit@optinfo$conv$opt,
                    use.names = FALSE
                )
                optimizer_codes <- vapply(raw_optimizer_codes, function(value) {
                    converted <- utils::type.convert(as.character(value),
                        as.is = TRUE, na.strings = character()
                    )
                    if (is.numeric(converted)) {
                        as.numeric(converted)
                    } else {
                        NA_real_
                    }
                }, numeric(1))
                optimizer_codes <- optimizer_codes[is.finite(optimizer_codes) &
                    optimizer_codes != 0]
                optimizer_messages <- if (length(optimizer_codes) >
                    0L) {
                    paste0("lme4 optimizer convergence code: ",
                        paste(unique(optimizer_codes),
                        collapse = ", "
                    ))
                } else {
                    character(0)
                }
                convergence_messages <- c(
                    fit_warnings, fit@optinfo$conv$lme4$messages,
                    optimizer_messages
                )
            }

            fit_warnings <- unique(as.character(convergence_messages))
            fit_warnings <-
                fit_warnings[!is.na(fit_warnings) & nzchar(fit_warnings)]
            convergence_messages <-
                fit_warnings[grepl("converg|unable to evaluate|degenerate|failed to|gradient|iteration limit",
                fit_warnings,
                ignore.case = TRUE
            )]

            required_statistics <- c(
                "Estimate", "Std. Error", "t value",
                "Pr(>|t|)"
            )
            if (!all(required_statistics %in% colnames(coef_table)) ||
                any(!is.finite(as.matrix(coef_table[, required_statistics,
                    drop = FALSE
                ])))) {
                stop("The mixed model returned missing or non-finite coefficient statistics.",
                    call. = FALSE
                )
            }

            coefficient_terms <- buildCoefficientTermMapMethylationModels(
                formulaText = formulaText,
                data = model_data, removeRandomEffects = !identical(
                    lmeEngine,
                    "nlme"
                )
            )
            omnibus_result <- NULL
            if (isTRUE(omnibusTest)) {
                if (identical(lmeEngine, "nlme")) {
                    stop(
                        "omnibusTest is not available for nlme fits.",
                        call. = FALSE
                    )
                }
                omnibus_result <- computeOmnibusTestMethylationLME(
                    fit = fit, coefficientTerms = coefficient_terms,
                    omnibusTerm = omnibusTerm,
                    omnibusDdf = omnibusDdf
                )
            }

            list(
                coef = coef_table, residuals = stats::residuals(fit),
                fitted = stats::fitted(fit), ranef = ranef_values,
                fixef = fixef_values, coefficientTerms = coefficient_terms,
                omnibus = omnibus_result,
                fitStatus = list(singular = singular_fit,
                    converged = length(convergence_messages) ==
                    0L, warnings = fit_warnings,
                        convergenceMessages = convergence_messages)
            )
        },
        error = function(error) {
            reason <- conditionMessage(error)
            convergence_failure <-
                grepl("converg|unable to evaluate|degenerate|failed to|gradient|iteration limit",
                reason,
                ignore.case = TRUE
            )
            newMethylationFitErrorDnaEpico(
                reason = reason,
                    errorClass = "dnaEPICO_methylationLME_fit_error",
                status = if (convergence_failure) {
                    "not_converged"
                } else {
                    "failed"
                }, converged = if (convergence_failure) {
                    FALSE
                } else {
                    NA
                }, convergenceMessage = if (convergence_failure) {
                    reason
                } else {
                    NA_character_
                }
            )
        }
    )
}

fitMethylationLMEBatch <- function(
    cpgBatch, data, modelData,
    formulaText, personVar, lmeEngine = "lme4", correlationStructure = "none",
    correlationTimeVar = NULL, phenotype, interactionTerm = NULL,
    responseVar = "beta", invalidCpgReasons = character(0),
    omnibusTest = FALSE, omnibusDdf = "Satterthwaite",
    omnibusTerm = NULL
) {
    fits <- vector("list", length(cpgBatch))
    names(fits) <- cpgBatch
    summaries <- vector("list", length(cpgBatch))
    names(summaries) <- cpgBatch

    for (cpg in cpgBatch) {
        invalid_reason <- if (cpg %in% names(invalidCpgReasons)) {
            invalidCpgReasons[[cpg]]
        } else {
            NULL
        }
        model_obj <- if (!is.null(invalid_reason)) {
            newMethylationFitErrorDnaEpico(
                reason = invalid_reason,
                errorClass = "dnaEPICO_methylationLME_fit_error",
                status = "invalid"
            )
        } else {
            fitCpGModelMethylationLME(
                cpg = cpg, cpgValues = data[[cpg]],
                modelData = modelData, formulaText = formulaText,
                personVar = personVar, lmeEngine = lmeEngine,
                correlationStructure = correlationStructure,
                correlationTimeVar = correlationTimeVar,
                    responseVar = responseVar,
                omnibusTest = omnibusTest, omnibusDdf = omnibusDdf,
                omnibusTerm = omnibusTerm
            )
        }
        fits[[cpg]] <- model_obj
        summaries[[cpg]] <- summarizeCpGFitMethylationLME(
            cpg = cpg,
            modelObj = model_obj, phenotype = phenotype,
                interactionTerm = interactionTerm
        )
    }

    summaries <- Filter(Negate(is.null), summaries)
    summary_df <- if (length(summaries) == 0L) {
        data.frame()
    } else {
        out <- do.call(rbind, summaries)
        rownames(out) <- NULL
        out
    }

    list(fits = fits, summaries = summary_df)
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
#'     data = ex$preparedData$data,
#'     timeVar = "Timepoint",
#'     phenotypes = "score",
#'     verbose = FALSE,
#'     logs = FALSE
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
        stop("Phenotype columns not found for timepoint summary: ",
            paste(missing_phenotypes, collapse = ", "),
            call. = FALSE
        )
    }

    split_data <- split(data, as.character(data[[timeVar]]))
    summaries <- lapply(split_data, function(df) {
        out <- list()
        out[[timeVar]] <- as.character(df[[timeVar]][[1L]])

        for (phenotype in phenotype_list) {
            values <- df[[phenotype]]
            if (is.numeric(values)) {
                finite_values <- values[is.finite(values)]
                out[[paste0(phenotype, "_mean")]] <-
                    if (length(finite_values)) {
                    mean(finite_values)
                } else {
                    NA_real_
                }
                out[[paste0(phenotype, "_sd")]] <- if (length(finite_values) >
                    1L) {
                    stats::sd(finite_values)
                } else {
                    NA_real_
                }
                out[[paste0(phenotype, "_n")]] <- length(finite_values)
            } else {
                observed_levels <- unique(as.character(values[!is.na(values)]))
                out[[paste0(phenotype, "_n")]] <- sum(!is.na(values))
                out[[paste0(phenotype, "_levels")]] <- paste(observed_levels,
                    collapse = ","
                )
            }
        }

        as.data.frame(out, stringsAsFactors = FALSE)
    })
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
#'     inputPheno = ex$inputPath,
#'     personVar = "person",
#'     timeVar = "Timepoint",
#'     phenotypes = "score",
#'     covariates = "sex",
#'     factorVars = "sex",
#'     cpgLimit = 2,
#'     verbose = FALSE,
#'     logs = FALSE
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
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    phenotype_list <- unique(splitOptionMinfiEwasWater(phenotypes,
        sep = ","
    ))
    covariate_list <- splitOptionMinfiEwasWater(covariates, sep = ",")
    factor_list <- splitOptionMinfiEwasWater(factorVars, sep = ",")
    scale_list <- normalizeScaleVariablesDnaEpico(scaleVars)
    prs_map <- parsePrsMapMethylationGLM(prsMap)
    cpg_limit <-
        validateCpgLimitMethylationModels(normalizeOptionalNumericMethylationGLM(cpgLimit))
    cpgPrefix <- validateCpgPrefixDnaEpico(cpgPrefix)
    methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
    methylation_label <-
        methylationScaleResponseLabelDnaEpico(methylation_scale)
    methylation_prefix <-
        methylationScaleObjectPrefixDnaEpico(methylation_scale)
    response_column <- methylationScaleResponseColumnDnaEpico(methylation_scale)
    input_object_name <- sub("[.][^.]+$", "", basename(inputPheno))
    scale_named_input <- if (startsWith(input_object_name,
        methylation_prefix)) {
        input_object_name
    } else {
        character(0)
    }
    legacy_input_name <- if (identical(methylation_scale, "beta")) {
        "phenoBT1T2"
    } else {
        character(0)
    }
    analysis_data <- loadSavedObjectPreprocessingPheno(inputPheno,
        preferred_name = c(scale_named_input, paste0(
            methylation_prefix,
            "T1T2"
        ), legacy_input_name, input_object_name)
    )

    if (!is.data.frame(analysis_data)) {
        analysis_data <- as.data.frame(analysis_data, stringsAsFactors = FALSE)
    }

    if (length(phenotype_list) == 0L) {
        stop("At least one phenotype must be supplied.", call. = FALSE)
    }

    person_data <- ensurePersonColumnMethylationLME(
        data = analysis_data,
        personVar = personVar
    )
    analysis_data <- person_data$data

    if (!(timeVar %in% colnames(analysis_data))) {
        stop("timeVar column not found in inputPheno: ", timeVar,
            call. = FALSE
        )
    }
    time_values <- analysis_data[[timeVar]]
    invalid_time <-
        is.na(time_values) | !nzchar(trimws(as.character(time_values)))
    if (is.numeric(time_values)) {
        invalid_time <- invalid_time | !is.finite(time_values)
    }
    if (any(invalid_time)) {
        stop("timeVar contains missing, blank, or non-finite values.",
            call. = FALSE
        )
    }

    missing_phenotypes <- setdiff(phenotype_list, colnames(analysis_data))
    if (length(missing_phenotypes) > 0L) {
        stop("Phenotype columns not found in inputPheno: ",
            paste(missing_phenotypes,
            collapse = ", "
        ), call. = FALSE)
    }

    missing_covariates <- setdiff(covariate_list, colnames(analysis_data))
    if (length(missing_covariates) > 0L) {
        stop("Covariate columns not found in inputPheno: ",
            paste(missing_covariates,
            collapse = ", "
        ), call. = FALSE)
    }

    mapped_prs <- unname(prs_map[names(prs_map) %in% phenotype_list])
    missing_prs <- setdiff(mapped_prs, colnames(analysis_data))
    if (length(missing_prs) > 0L) {
        stop("PRS columns not found in inputPheno: ", paste(unique(missing_prs),
            collapse = ", "
        ), call. = FALSE)
    }

    resolved_interaction <- normalizeOptionalColumnMethylationModels(
        interactionTerm,
        "interactionTerm"
    )
    if (!is.null(resolved_interaction) && nzchar(resolved_interaction)) {
        if (!(resolved_interaction %in% colnames(analysis_data))) {
            stop("interactionTerm column not found in inputPheno: ",
                resolved_interaction,
                call. = FALSE
            )
        }
    }

    missing_factor_vars <- setdiff(factor_list, colnames(analysis_data))
    if (length(missing_factor_vars) > 0L) {
        stop("Factor columns not found in inputPheno: ",
            paste(missing_factor_vars,
            collapse = ", "
        ), call. = FALSE)
    }

    for (var in intersect(c(personVar, factor_list), colnames(analysis_data))) {
        if (identical(var, personVar)) {
            next
        }
        analysis_data[[var]] <- as.factor(analysis_data[[var]])
    }

    cpg_columns <- grep(paste0("^", escapeRegexMethylationGLM(cpgPrefix)),
        colnames(analysis_data),
        value = TRUE
    )
    if (!is.na(cpg_limit)) {
        cpg_columns <- utils::head(cpg_columns, as.integer(cpg_limit))
    }
    if (length(cpg_columns) == 0L) {
        stop("No CpG columns were found with prefix '", cpgPrefix,
            "'.",
            call. = FALSE
        )
    }
    validateMethylationProbeIdentifiersDnaEpico(
        cpg_columns,
        "CpG columns in the LME input"
    )
    methylation_validation <- inspectMethylationColumnsDnaEpico(
        data = analysis_data,
        cpgColumns = cpg_columns, methylationScale = methylation_scale
    )
    analysis_data <- methylation_validation$data
    if (methylation_validation$boundaries$NaN.Converted > 0L) {
        warning(methylation_validation$boundaries$NaN.Converted,
            " NaN methylation values were converted to NA. See methylationIssues.",
            call. = FALSE
        )
    }
    if (nrow(methylation_validation$invalidCpGs) > 0L) {
        warning(nrow(methylation_validation$invalidCpGs),
            " invalid CpG columns will be reported as failed fits instead of stopping the analysis.",
            call. = FALSE
        )
    }

    fixed_effect_variables <- unique(c(
        phenotype_list, covariate_list,
        mapped_prs, resolved_interaction
    ))
    fixed_effect_variables <-
        fixed_effect_variables[!is.na(fixed_effect_variables) &
        nzchar(fixed_effect_variables)]
    scaled_cpgs <- intersect(scale_list, cpg_columns)
    if (length(scaled_cpgs) > 0L) {
        stop("CpG methylation response columns cannot be listed in scaleVars: ",
            paste(scaled_cpgs, collapse = ", "),
            call. = FALSE
        )
    }
    model_columns <- setdiff(colnames(analysis_data), cpg_columns)
    scaling <- scaleModelVariablesDnaEpico(
        data = analysis_data[,
            model_columns,
            drop = FALSE
        ], scaleVars = scale_list,
        factorVars = factor_list, eligibleVars = fixed_effect_variables,
        protectedVars = c(personVar, cpg_columns)
    )

    timepoint_summary <- summarizeTimepointsMethylationLME(
        data = analysis_data,
        timeVar = timeVar, phenotypes = phenotype_list, verbose = FALSE,
        logs = FALSE
    )

    log_lines <- c(
        "============================================================",
        paste(
            "Loaded longitudinal phenotype + methylation data from:",
            inputPheno
        ), paste(
            "Merged modeling object:            ",
            methylation_prefix, "*"
        ), paste(
            "Data dimensions:                  ",
            paste(dim(analysis_data), collapse = " x ")
        ), paste(
            "Person variable:                  ",
            personVar
        ), paste(
            "Time variable:                    ",
            timeVar
        ), paste(
            "Phenotypes:                       ",
            paste(phenotype_list, collapse = ", ")
        ), paste(
            "Covariates:                       ",
            paste(covariate_list, collapse = ", ")
        ), paste(
            "Factor variables:                 ",
            paste(factor_list, collapse = ", ")
        ), formatScalingMetadataLogDnaEpico(scaling$metadata),
        paste("CpG columns retained:             ", length(cpg_columns)),
        formatMethylationBoundariesLogDnaEpico(methylation_validation$boundaries),
        if (nrow(methylation_validation$issues) > 0L) {
            c("Methylation value issues:",
                previewLinesMinfiEwasWater(methylation_validation$issues))
        } else {
            "Methylation value issues:            none"
        }, if (isTRUE(person_data$personCreated)) {
            paste("Created person variable from SID: ", personVar)
        } else {
            paste("Person variable already present:  ", personVar)
        }, "Count of records per person ID:",
            previewLinesMinfiEwasWater(person_data$personCounts),
        paste("Values observed in", timeVar, ":"),
            previewLinesMinfiEwasWater(table(analysis_data[[timeVar]],
            useNA = "ifany"
        ))
    )

    if (!is.null(person_data$mappingPreview)) {
        log_lines <- c(
            log_lines, "Example mapping of SID to person ID:",
            previewLinesMinfiEwasWater(person_data$mappingPreview)
        )
    }

    log_lines <- c(
        log_lines, "Summary statistics for phenotype scores by timepoint:",
        previewLinesMinfiEwasWater(timepoint_summary),
            "============================================================"
    )
    emitLogMinfiEwasWater(log_lines, verbose = verbose, log_path = log_path)

    structure(
        list(
            data = analysis_data, modelData = scaling$data,
            personVar = personVar, timeVar = timeVar,
                phenotypes = phenotype_list,
            covariates = covariate_list, factorVars = factor_list,
            scaleVars = scaling$scaleVars, scalingMetadata = scaling$metadata,
            prsMap = prs_map, cpgColumns = cpg_columns, cpgPrefix = cpgPrefix,
            cpgLimit = cpg_limit, methylationScale = methylation_scale,
            responseLabel = methylation_label,
                methylationObjectPrefix = methylation_prefix,
            internalResponseColumn = response_column,
                interactionTerm = resolved_interaction,
            requestedInteractionTerm = interactionTerm,
                methylationBoundaries = methylation_validation$boundaries,
            methylationIssues = methylation_validation$issues,
                invalidCpGs = methylation_validation$invalidCpGs,
            personCreated = person_data$personCreated,
                personCounts = person_data$personCounts,
            personMappingPreview = person_data$mappingPreview,
                timepointSummary = timepoint_summary
        ),
        class = "dnaEPICO_methylationLME_data"
    )
}
#' Fit CpG-wise mixed-effects models for longitudinal methylation analyses
#'
#' @param preparedData Object returned by `prepareMethylationLMEData()`.
#' @param nCores Integer. Maximum number of worker processes to use. Automatic
#'   fitting uses the empirical lme4 or nlme crossover and caps workers by
#'   available CpGs and detected physical cores.
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
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationLME_models'`
#'   containing fitted model lists, model formulas, and counts of failed CpG
#'   fits.
#'
#' @description
#' Fit one linear mixed-effects model per CpG for each phenotype requested in
#' the
#' object returned by `prepareMethylationLMEData()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' model_results <- fitMethylationLMEModels(
#'     preparedData = ex$preparedData,
#'     nCores = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(model_results$fits)
#'
#' @export
fitMethylationLMEModels <- function(
    preparedData, nCores = 1L,
    libPath = NULL, lmeLibs = "lme4,lmerTest", correlationStructure = "none",
    correlationVar = NULL, omnibusTest = FALSE,
    omnibusDdf = "Satterthwaite", verbose = FALSE, logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    if (is.null(libPath)) {
        libPath <- .libPaths()
    }

    lme_config <- resolveLmeLibrariesMethylationLME(lmeLibs)
    lme_lib_list <- lme_config$requestedPackages
    required_lme_lib_list <- lme_config$requiredPackages
    lme_engine <- lme_config$engine
    omnibus_config <- validateOmnibusConfigurationMethylationLME(
        omnibusTest = omnibusTest, omnibusDdf = omnibusDdf,
        lmeEngine = lme_engine
    )
    if (isTRUE(omnibus_config$test) &&
        identical(omnibus_config$ddf, "Kenward-Roger")) {
        required_lme_lib_list <- unique(c(
            required_lme_lib_list,
            "pbkrtest"
        ))
    }
    correlation_structure <-
        normalizeCorrelationStructureMethylationLME(correlationStructure)
    correlation_var <-
        normalizeCorrelationVariableMethylationLME(correlationVar = correlationVar)
    if (!identical(lme_engine, "nlme") && !identical(
        correlation_structure,
        "none"
    )) {
        stop("correlationStructure can only be AR1 or CAR1 when lmeLibs selects 'nlme'.",
            call. = FALSE
        )
    }
    if (!identical(correlation_structure, "none")) {
        if (is.null(correlation_var)) {
            stop("correlationVar must be supplied when correlationStructure is AR1 or CAR1.",
                call. = FALSE
            )
        }
        if (correlation_var %in% preparedData$scaleVars) {
            stop("correlationVar cannot also be listed in scaleVars; provide a separate scaled fixed-effect column if required.",
                call. = FALSE
            )
        }
    }

    analysis_data <- preparedData$data
    model_data <- if (!is.null(preparedData$modelData)) {
        preparedData$modelData
    } else {
        analysis_data
    }
    cpg_columns <- preparedData$cpgColumns
    invalid_cpg_reasons <- invalidMethylationReasonsDnaEpico(preparedData)
    n_cores <- validatePositiveIntegerMethylationModels(
        nCores,
        "nCores"
    )
    fits <- list()
    summary_cache <- list()
    omnibus_tests <- list()
    formulas <- stats::setNames(
        character(length(preparedData$phenotypes)),
        preparedData$phenotypes
    )
    failure_counts <- stats::setNames(
        integer(length(preparedData$phenotypes)),
        preparedData$phenotypes
    )
    failure_reasons <- list()
    singular_counts <- stats::setNames(
        integer(length(preparedData$phenotypes)),
        preparedData$phenotypes
    )
    convergence_issue_counts <- stats::setNames(
        integer(length(preparedData$phenotypes)),
        preparedData$phenotypes
    )
    omnibus_targets <- stats::setNames(
        character(length(preparedData$phenotypes)),
        preparedData$phenotypes
    )
    parallel_plan <- resolveParallelPlanMethylationModels(
        engine = lme_engine,
        nCores = n_cores, nCpGs = length(cpg_columns)
    )
    backend <- parallel_plan$backend
    worker_count <- parallel_plan$workerCount
    cpg_batches <- chunkCpGColumnsMethylationModels(
        cpgColumns = cpg_columns,
        nCores = worker_count, batchesPerCore = 8L
    )
    psock_cluster <- NULL
    if (identical(backend, "psock") && length(cpg_batches) >
        1L) {
        psock_cluster <- makePsockClusterMethylationModels(min(
            worker_count,
            length(cpg_batches)
        ))
        on.exit(
            {
                if (!is.null(psock_cluster)) {
                    try(parallel::stopCluster(psock_cluster), silent = TRUE)
                }
            },
            add = TRUE
        )
        parallel::clusterExport(psock_cluster,
            varlist = c(
                "analysis_data",
                "invalid_cpg_reasons", "libPath", "required_lme_lib_list",
                "validateWorkerPackagesMethylationModels",
                    "newMethylationFitErrorDnaEpico",
                "fitMethylationLMEBatch", "fitCpGModelMethylationLME",
                "buildCorrelationMethylationLME",
                    "coerceCoefficientTableMethylationLME",
                "buildCoefficientTermMapMethylationModels",
                    "removeRandomInterceptMethylationModels",
                "computeOmnibusTestMethylationLME",
                    "emptyOmnibusResultMethylationLME",
                "summarizeCpGFitMethylationLME",
                    "fitStatusValuesMethylationLME",
                "findCoefficientRowsMethylationLME",
                    "findCoefficientRowsMethylationGLM",
                "quoteNamesMethylationGLM", "escapeRegexMethylationGLM"
            ),
            envir = environment()
        )
        parallel::clusterEvalQ(psock_cluster,
            validateWorkerPackagesMethylationModels(
            libPath = libPath,
            packages = required_lme_lib_list
        ))
    }

    for (phenotype in preparedData$phenotypes) {
        prs_var <- character(0)
        if (phenotype %in% names(preparedData$prsMap)) {
            prs_var <- unname(preparedData$prsMap[[phenotype]])
        }
        covariates <- unique(c(preparedData$covariates, prs_var))
        model_vars <- unique(c(
            preparedData$personVar, phenotype,
            covariates, preparedData$interactionTerm, if (!identical(
                correlation_structure,
                "none"
            )) {
                correlation_var
            } else {
                character(0)
            }
        ))
        model_vars <- model_vars[!is.na(model_vars) & nzchar(model_vars)]

        missing_vars <- setdiff(model_vars, colnames(model_data))
        if (length(missing_vars) > 0L) {
            stop("Model variables not found for phenotype ",
                phenotype, ": ", paste(missing_vars, collapse = ", "),
                call. = FALSE
            )
        }

        display_formula_text <- buildFormulaMethylationLME(
            phenotype = phenotype,
            personVar = preparedData$personVar, covariates = covariates,
            interactionTerm = preparedData$interactionTerm,
                includeRandomTerm = TRUE,
            responseVar = preparedData$internalResponseColumn
        )
        formula_text <- buildFormulaMethylationLME(
            phenotype = phenotype,
            personVar = preparedData$personVar, covariates = covariates,
            interactionTerm = preparedData$interactionTerm,
                includeRandomTerm = !identical(
                lme_engine,
                "nlme"
            ), responseVar = preparedData$internalResponseColumn
        )

        base_model_data <- model_data[, model_vars, drop = FALSE]
        factor_vars <- preparedData$factorVars
        person_var <- preparedData$personVar
        base_model_data[[person_var]] <-
            as.factor(base_model_data[[person_var]])
        for (var in intersect(factor_vars, colnames(base_model_data))) {
            base_model_data[[var]] <- as.factor(base_model_data[[var]])
        }
        validateFixedEffectDesignMethylationModels(
            formulaText = formula_text,
            data = base_model_data, removeRandomEffects = !identical(
                lme_engine,
                "nlme"
            )
        )
        omnibus_target <- NULL
        if (isTRUE(omnibus_config$test)) {
            omnibus_target <- resolveOmnibusTargetTermMethylationLME(
                formulaText = formula_text, data = base_model_data,
                phenotype = phenotype,
                interactionTerm = preparedData$interactionTerm
            )
            omnibus_targets[[phenotype]] <- omnibus_target
        }
        correlation_time_var <- NULL
        if (!identical(correlation_structure, "none")) {
            correlation_data <- addCorrelationTimeVariableMethylationLME(
                modelData = base_model_data,
                correlationVar = correlation_var
            )
            base_model_data <- correlation_data$data
            correlation_time_var <- correlation_data$correlationTimeVar
            validateCorrelationTimeMethylationLME(
                modelData = base_model_data,
                correlationTimeVar = correlation_time_var,
                    personVar = person_var,
                correlationStructure = correlation_structure
            )
        }
        batch_worker <- fitMethylationLMEBatch
        resolved_interaction <- preparedData$interactionTerm
        response_var <- preparedData$internalResponseColumn

        if (!identical(backend, "serial") && length(cpg_batches) >
            1L) {
            cluster_size <- min(worker_count, length(cpg_batches))

            if (identical(backend, "fork")) {
                batch_results <- parallel::mclapply(cpg_batches,
                    function(batch) {
                        validateWorkerPackagesMethylationModels(
                            libPath = libPath,
                            packages = required_lme_lib_list
                        )
                        batch_worker(
                            cpgBatch = batch, data = analysis_data,
                            modelData = base_model_data,
                                formulaText = formula_text,
                            personVar = person_var, lmeEngine = lme_engine,
                            correlationStructure = correlation_structure,
                            correlationTimeVar = correlation_time_var,
                            phenotype = phenotype,
                                interactionTerm = resolved_interaction,
                            responseVar = response_var,
                                invalidCpgReasons = invalid_cpg_reasons,
                            omnibusTest = omnibus_config$test,
                            omnibusDdf = omnibus_config$ddf,
                            omnibusTerm = omnibus_target
                        )
                    },
                    mc.cores = cluster_size, mc.preschedule = FALSE
                )
            } else {
                parallel::clusterExport(psock_cluster,
                    varlist = c(
                        "base_model_data",
                        "formula_text", "phenotype", "resolved_interaction",
                        "person_var", "lme_engine", "correlation_structure",
                        "correlation_time_var", "response_var", "batch_worker",
                        "omnibus_config", "omnibus_target"
                    ),
                    envir = environment()
                )
                batch_results <- parallel::parLapplyLB(
                    psock_cluster,
                    cpg_batches, function(batch) {
                        batch_worker(
                            cpgBatch = batch, data = analysis_data,
                            modelData = base_model_data,
                                formulaText = formula_text,
                            personVar = person_var, lmeEngine = lme_engine,
                            correlationStructure = correlation_structure,
                            correlationTimeVar = correlation_time_var,
                            phenotype = phenotype,
                                interactionTerm = resolved_interaction,
                            responseVar = response_var,
                                invalidCpgReasons = invalid_cpg_reasons,
                            omnibusTest = omnibus_config$test,
                            omnibusDdf = omnibus_config$ddf,
                            omnibusTerm = omnibus_target
                        )
                    }
                )
            }
        } else {
            validateWorkerPackagesMethylationModels(
                libPath = libPath,
                packages = required_lme_lib_list
            )
            batch_results <- lapply(cpg_batches, function(batch) {
                batch_worker(
                    cpgBatch = batch, data = analysis_data,
                    modelData = base_model_data, formulaText = formula_text,
                    personVar = person_var, lmeEngine = lme_engine,
                    correlationStructure = correlation_structure,
                    correlationTimeVar = correlation_time_var,
                    phenotype = phenotype,
                        interactionTerm = resolved_interaction,
                    responseVar = response_var,
                        invalidCpgReasons = invalid_cpg_reasons,
                    omnibusTest = omnibus_config$test,
                    omnibusDdf = omnibus_config$ddf,
                    omnibusTerm = omnibus_target
                )
            })
        }

        combined_results <- combineFitBatchResultsMethylationModels(
            batchResults = batch_results,
            cpgColumns = cpg_columns
        )
        fit_list <- combined_results$fits
        phenotype_omnibus <- collectOmnibusTestsMethylationLME(
            fits = fit_list,
            phenotype = phenotype
        )
        phenotype_summary_cache <- filterSummaryByPvalueMethylationLME(
            summaryDf = combined_results$summaries,
            pValueFilter = NA_real_
        )
        failures <- vapply(fit_list, function(x) {
            inherits(
                x,
                "dnaEPICO_methylationLME_fit_error"
            )
        }, logical(1))
        error_counts <- summarizeFitErrorsMethylationModels(
            fitList = fit_list,
            errorClass = "dnaEPICO_methylationLME_fit_error"
        )
        valid_fits <- fit_list[!failures]
        singular_fits <- vapply(
            valid_fits, function(x) isTRUE(x$fitStatus$singular),
            logical(1)
        )
        convergence_issues <- vapply(valid_fits, function(x) {
            identical(
                x$fitStatus$converged,
                FALSE
            )
        }, logical(1))

        fits[[phenotype]] <- fit_list
        summary_cache[[phenotype]] <- phenotype_summary_cache
        omnibus_tests[[phenotype]] <- phenotype_omnibus
        formulas[[phenotype]] <- display_formula_text
        failure_counts[[phenotype]] <- sum(failures)
        failure_reasons[[phenotype]] <- error_counts
        singular_counts[[phenotype]] <- sum(singular_fits)
        convergence_issue_counts[[phenotype]] <- sum(convergence_issues)

        emitLogMinfiEwasWater(
            c(
                "============================================================",
                paste("Fitted phenotype:            ", phenotype),
                paste("Formula:                     ", display_formula_text),
                paste("Correlation structure:       ", correlation_structure),
                paste("Correlation variable:        ", if (identical(
                    correlation_structure,
                    "none"
                )) {
                    "None"
                } else {
                    correlation_var
                }), paste(
                    "CpGs attempted:              ",
                    length(cpg_columns)
                ), paste(
                    "Failed CpG fits:             ",
                    failure_counts[[phenotype]]
                ), paste(
                    "Singular CpG fits:           ",
                    singular_counts[[phenotype]]
                ), paste(
                    "Non-converged CpG fits:      ",
                    convergence_issue_counts[[phenotype]]
                ), paste(
                    "Omnibus tests requested:     ",
                    isTRUE(omnibus_config$test)
                ), paste(
                    "Omnibus target term:         ",
                    if (is.null(omnibus_target)) "None" else omnibus_target
                ), paste(
                    "Omnibus denominator DF:      ",
                    if (isTRUE(omnibus_config$test)) {
                        omnibus_config$ddf
                    } else {
                        "None"
                    }
                ), paste(
                    "Successful omnibus tests:    ",
                    sum(phenotype_omnibus$Omnibus.Status == "tested")
                ), paste(
                    "Unavailable omnibus tests:   ",
                    sum(phenotype_omnibus$Omnibus.Status != "tested")
                ), paste(
                    "Top fit errors:              ",
                    formatFitErrorsMethylationModels(error_counts)
                ),
                paste("Parallel backend:            ", backend),
                paste("Effective workers:           ", worker_count),
                paste("Parallel crossover CpGs:     ",
                    parallel_plan$crossoverCpGs),
                paste("Parallel selection:          ", parallel_plan$reason),
                paste("Fit batches:                 ", length(cpg_batches)),
                paste("Fit batch size:              ",
                    if (length(cpg_batches) ==
                    0L) {
                    0L
                } else {
                    max(vapply(cpg_batches, length, integer(1)))
                }),
                paste("Fit-time summary rows cached:",
                    nrow(phenotype_summary_cache)),
                "============================================================"
            ),
            verbose = verbose, log_path = log_path
        )

        if (length(fit_list) > 0L && all(failures)) {
            warning("All CpG ", toupper(lme_engine),
                " fits failed for phenotype '",
                phenotype, "'. Top failure reasons: ",
                    formatFitErrorsMethylationModels(error_counts),
                ". The failure inventory was retained and the analysis continued.",
                call. = FALSE
            )
        }
    }

    fit_failures <- collectFitFailuresMethylationModels(
        fits = fits,
        errorClass = "dnaEPICO_methylationLME_fit_error"
    )
    fit_diagnostics <- collectFitDiagnosticsMethylationLME(fits)

    if (!is.null(psock_cluster)) {
        parallel::stopCluster(psock_cluster)
        psock_cluster <- NULL
    }

    structure(
        list(
            fits = fits, summaryCache = summary_cache,
            omnibusTests = omnibus_tests, omnibusTargets = omnibus_targets,
            formulas = formulas, phenotypes = names(fits),
                failureCounts = failure_counts,
            failureReasons = failure_reasons, fitFailures = fit_failures,
            fitDiagnostics = fit_diagnostics,
                methylationBoundaries = preparedData$methylationBoundaries,
            methylationIssues = preparedData$methylationIssues,
                invalidCpGs = preparedData$invalidCpGs,
            singularCounts = singular_counts,
                convergenceIssueCounts = convergence_issue_counts,
            settings = list(
                nCores = n_cores, parallelBackend = backend,
                workerCount = worker_count,
                    resourceWorkerCap = parallel_plan$resourceWorkerCap,
                parallelCrossoverCpGs = parallel_plan$crossoverCpGs,
                parallelSelectionReason = parallel_plan$reason,
                    clusterReusedAcrossPhenotypes = identical(
                    backend,
                    "psock"
                ) && length(preparedData$phenotypes) >
                    1L, fitBatchCount = length(cpg_batches), libPath = libPath,
                lmeLibs = lme_lib_list, lmeEngine = lme_engine,
                    correlationStructure = correlation_structure,
                correlationVar = if (identical(
                    correlation_structure,
                    "none"
                )) {
                    NULL
                } else {
                    correlation_var
                }, methylationScale = preparedData$methylationScale,
                methylationObjectPrefix = preparedData$methylationObjectPrefix,
                responseLabel = preparedData$responseLabel,
                    internalResponseColumn = preparedData$internalResponseColumn,
                interactionTerm = preparedData$interactionTerm,
                omnibusTest = omnibus_config$test,
                omnibusDdf = omnibus_config$ddf,
                omnibusRhs = 0,
                omnibusJoint = TRUE,
                omnibusEps = sqrt(.Machine$double.eps),
                    phenotypes = preparedData$phenotypes,
                covariates = preparedData$covariates,
                    factorVars = preparedData$factorVars,
                factorLevels = lapply(stats::setNames(
                    preparedData$factorVars,
                    preparedData$factorVars
                ), function(variable) levels(preparedData$data[[variable]])),
                scaleVars = preparedData$scaleVars,
                    scalingMetadata = preparedData$scalingMetadata,
                sampleCount = nrow(preparedData$data),
                    personVar = preparedData$personVar,
                timeVar = preparedData$timeVar
            ), responseLabel = preparedData$responseLabel
        ),
        class = "dnaEPICO_methylationLME_models"
    )
}

summarizeOmnibusTestsMethylationLME <- function(
    modelResults, padjmethod = "fdr",
    excludeSingular = FALSE, excludeNonConverged = FALSE
) {
    adjustment_method <- validatePAdjustmentMethodMethylationModels(
        padjmethod
    )
    omnibus_tables <- modelResults$omnibusTests
    if (!is.list(omnibus_tables)) {
        omnibus_tables <- list()
    }

    summaries <- lapply(names(modelResults$fits), function(phenotype) {
        table <- omnibus_tables[[phenotype]]
        if (!is.data.frame(table) || nrow(table) == 0L) {
            return(data.frame())
        }

        fits <- modelResults$fits[[phenotype]]
        for (index in seq_len(nrow(table))) {
            cpg <- table$CpG[[index]]
            fit <- fits[[cpg]]
            if (is.null(fit) || inherits(
                fit,
                "dnaEPICO_methylationLME_fit_error"
            )) {
                next
            }
            status <- fitStatusValuesMethylationLME(fit)
            reasons <- character(0)
            if (isTRUE(excludeSingular) && isTRUE(status$singular)) {
                reasons <- c(reasons, "singular random-effects fit")
            }
            if (isTRUE(excludeNonConverged) &&
                identical(status$converged, FALSE)) {
                reasons <- c(reasons, "model did not converge")
            }
            if (length(reasons) > 0L) {
                table$Omnibus.P.Value[[index]] <- NA_real_
                existing <- table$Omnibus.Reason[[index]]
                table$Omnibus.Reason[[index]] <- paste(
                    c(existing[!is.na(existing) & nzchar(existing)], reasons),
                    collapse = "; "
                )
            }
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
    names(summaries) <- names(modelResults$fits)
    summaries
}

#' Summarize CpG-wise mixed-effects model fits for longitudinal analyses
#'
#' @param modelResults Object returned by `fitMethylationLMEModels()`.
#' @param preparedData Object returned by `prepareMethylationLMEData()`.
#' @param summaryPval Numeric or `NA`. Optional p-value filter applied to the
#'   returned summary tables. `NA` keeps all rows.
#' @param excludeSingular Logical. If `TRUE`, retain singular CpG rows and
#'   reasons but set their inferential p-values to `NA`.
#' @param excludeNonConverged Logical. If `TRUE`, retain non-converged CpG rows
#'   and reasons but set their inferential p-values to `NA`.
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
#'   does not remove CpGs from those outputs. `fitDiagnostics` retains the
#'   phenotype-specific fit status, lme4 singular-fit flag, and warning text for
#'   every attempted CpG.
#'
#' @description
#' Extract phenotype-specific fixed-effect tables from the fitted mixed-effects
#' model object returned by `fitMethylationLMEModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' summary_results <- summarizeMethylationLMEModels(
#'     modelResults = ex$modelResults,
#'     preparedData = ex$preparedData,
#'     summaryPval = NA,
#'     nCores = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(summary_results$summaries)
#'
#' @export
summarizeMethylationLMEModels <- function(
    modelResults, preparedData,
    summaryPval = NA, excludeSingular = FALSE, excludeNonConverged = FALSE,
    padjmethod = "fdr", nCores = 1L, chunkSize = NULL,
    verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationLME.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    p_value_filter <- normalizeOptionalNumericMethylationGLM(summaryPval)
    p_value_filter <- validateProbabilityDnaEpico(p_value_filter,
        "summaryPval",
        allowNA = TRUE
    )
    chunk_size <- normalizeChunkSizeMethylationGLM(chunkSize)
    n_cores <- validatePositiveIntegerMethylationModels(
        nCores,
        "nCores"
    )
    adjustment_method <- validatePAdjustmentMethodMethylationModels(
        padjmethod
    )
    summaries <- list()
    diagnostic_summaries <- list()
    lme_engine <- modelResults$settings$lmeEngine
    if (is.null(lme_engine) || !nzchar(lme_engine)) {
        lme_engine <- "lme4"
    }
    if (identical(lme_engine, "nlme") && isTRUE(excludeSingular)) {
        warning("excludeSingular is not applicable to nlme fits and was ignored.",
            call. = FALSE
        )
    }
    remove_inapplicable_singular <- function(data) {
        if (identical(lme_engine, "nlme") && "Singular.Fit" %in%
            colnames(data)) {
            data$Singular.Fit <- NULL
        }
        data
    }

    for (phenotype in names(modelResults$fits)) {
        if (!is.null(modelResults$summaryCache) &&
            !is.null(modelResults$summaryCache[[phenotype]])) {
            diagnostic_df <- filterSummaryByPvalueMethylationLME(
                summaryDf = modelResults$summaryCache[[phenotype]],
                pValueFilter = NA_real_
            )
            diagnostic_df <- applyFitQualityExclusionsMethylationLME(
                summaryDf = diagnostic_df,
                excludeSingular = excludeSingular && identical(
                    lme_engine,
                    "lme4"
                ), excludeNonConverged = excludeNonConverged
            )
            diagnostic_df <- remove_inapplicable_singular(diagnostic_df)
            summary_df <- filterSummaryByPvalueMethylationLME(
                summaryDf = diagnostic_df,
                pValueFilter = p_value_filter
            )
            diagnostic_summaries[[phenotype]] <- diagnostic_df
            summaries[[phenotype]] <- summary_df

            emitLogMinfiEwasWater(
                c(
                    "============================================================",
                    paste("Summarized phenotype:        ", phenotype),
                    paste("LME summary rows returned:   ", nrow(summary_df)),
                    "Summary source:              fit-time cache",
                    if (is.na(p_value_filter)) {
                        "P-value filter:              none"
                    } else {
                        paste("P-value filter:              ", p_value_filter)
                    }, "============================================================"
                ),
                verbose = verbose, log_path = log_path
            )
            next
        }

        fit_list <- modelResults$fits[[phenotype]]
        cpg_names <- names(fit_list)

        local_chunk_size <- chunk_size
        if (is.null(local_chunk_size)) {
            local_chunk_size <- max(10L, floor(length(cpg_names) / max(n_cores *
                4L, 1L)))
        }
        local_chunk_size <- max(1L, as.integer(local_chunk_size))
        cpg_chunks <- split(cpg_names,
            ceiling(seq_along(cpg_names) / local_chunk_size))

        summary_worker <- summarizeCpGFitMethylationLME
        resolved_interaction <- preparedData$interactionTerm

        if (n_cores > 1L && length(cpg_chunks) > 1L) {
            cluster_size <- min(n_cores, length(cpg_chunks))
            cl <- parallel::makeCluster(cluster_size)
            on.exit(parallel::stopCluster(cl), add = TRUE)

            parallel::clusterExport(cl, varlist = c(
                "fit_list",
                "phenotype", "resolved_interaction", "summary_worker",
                "fitStatusValuesMethylationLME"
            ), envir = environment())

            result_chunks <- parallel::parLapplyLB(
                cl, cpg_chunks,
                function(chunk) {
                    rows <- lapply(chunk, function(cpg) {
                        summary_worker(
                            cpg = cpg, modelObj = fit_list[[cpg]],
                            phenotype = phenotype,
                                interactionTerm = resolved_interaction
                        )
                    })
                    rows <- Filter(Negate(is.null), rows)
                    if (length(rows) == 0L) {
                        return(NULL)
                    }
                    do.call(rbind, rows)
                }
            )
            parallel::stopCluster(cl)
            on.exit(NULL, add = FALSE)
        } else {
            result_chunks <- lapply(cpg_chunks, function(chunk) {
                rows <- lapply(chunk, function(cpg) {
                    summary_worker(
                        cpg = cpg, modelObj = fit_list[[cpg]],
                        phenotype = phenotype,
                            interactionTerm = resolved_interaction
                    )
                })
                rows <- Filter(Negate(is.null), rows)
                if (length(rows) == 0L) {
                    return(NULL)
                }
                do.call(rbind, rows)
            })
        }

        result_chunks <- Filter(Negate(is.null), result_chunks)
        if (length(result_chunks) == 0L) {
            summary_df <- data.frame()
        } else {
            summary_df <- do.call(rbind, result_chunks)
            rownames(summary_df) <- NULL
        }

        summary_df <- applyFitQualityExclusionsMethylationLME(
            summaryDf = summary_df,
            excludeSingular = excludeSingular && identical(
                lme_engine,
                "lme4"
            ), excludeNonConverged = excludeNonConverged
        )
        summary_df <- remove_inapplicable_singular(summary_df)
        diagnostic_summaries[[phenotype]] <- summary_df
        if (nrow(summary_df) > 0L && !is.na(p_value_filter)) {
            keep <- is.finite(summary_df$P.value) & summary_df$P.value <
                p_value_filter
            summary_df <- summary_df[keep, , drop = FALSE]
            rownames(summary_df) <- NULL
        }

        summaries[[phenotype]] <- summary_df

        emitLogMinfiEwasWater(
            c(
                "============================================================",
                paste("Summarized phenotype:        ", phenotype),
                paste("LME summary rows returned:   ", nrow(summary_df)),
                paste("Summary chunk size:          ", local_chunk_size),
                if (is.na(p_value_filter)) {
                    "P-value filter:              none"
                } else {
                    paste("P-value filter:              ", p_value_filter)
                }, "============================================================"
            ),
            verbose = verbose, log_path = log_path
        )
    }

    omnibus_summaries <- summarizeOmnibusTestsMethylationLME(
        modelResults = modelResults, padjmethod = adjustment_method,
        excludeSingular = excludeSingular && identical(lme_engine, "lme4"),
        excludeNonConverged = excludeNonConverged
    )

    structure(list(
        summaries = summaries, diagnosticSummaries = diagnostic_summaries,
        omnibusTests = omnibus_summaries,
        phenotypes = names(summaries), fitFailures = modelResults$fitFailures,
        fitDiagnostics = resolveFitDiagnosticsMethylationLME(modelSummaries = list(fitDiagnostics = if (!is.null(modelResults$fitDiagnostics)) {
            modelResults$fitDiagnostics
        } else {
            collectFitDiagnosticsMethylationLME(modelResults$fits)
        }, fitFailures = modelResults$fitFailures),
            summaryList = diagnostic_summaries),
        settings = list(
            summaryPval = p_value_filter,
                excludeSingular = isTRUE(excludeSingular),
            excludeNonConverged = isTRUE(excludeNonConverged),
            padjmethod = adjustment_method,
            chunkSize = chunk_size,
                interactionTerm = preparedData$interactionTerm,
            lmeEngine = modelResults$settings$lmeEngine
        )
    ), class = "dnaEPICO_methylationLME_summaries")
}

#' Collect significant longitudinal model terms from fitted mixed-effects models
#'
#' @param modelResults Object returned by `fitMethylationLMEModels()`.
#' @param pvalThreshold Numeric. Threshold applied to the extracted phenotype or
#'   interaction p-values.
#' @param interactionTerm Character or `NULL`. Optional interaction term. When
#'   `NULL`, phenotype main effects are used.
#' @param excludeSingular Logical. If `TRUE`, do not collect significant CpGs
#'   with singular lme4 fits.
#' @param excludeNonConverged Logical. If `TRUE`, do not collect significant
#'   CpGs whose mixed model did not converge. Excluded CpG rows and reasons
#'   remain available in diagnostic and annotated output.
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
#'     modelResults = ex$modelResults,
#'     pvalThreshold = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(significant_hits)
#'
#' @export
collectSignificantInteractionsMethylationLME <- function(
    modelResults,
    pvalThreshold = 0.05, interactionTerm = NULL, excludeSingular = FALSE,
    excludeNonConverged = FALSE, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationLME.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    threshold <- validateProbabilityDnaEpico(pvalThreshold, "pvalThreshold")
    retained <- list()

    for (phenotype in names(modelResults$fits)) {
        fit_list <- modelResults$fits[[phenotype]]
        phenotype_hits <- list()
        use_omnibus <- isTRUE(modelResults$settings$omnibusTest)
        if (isTRUE(use_omnibus) && !optionalTermMatchesMethylationModels(
            requested = interactionTerm,
            cached = modelResults$settings$interactionTerm
        )) {
            stop(
                "interactionTerm does not match the term used for the fitted omnibus tests.",
                call. = FALSE
            )
        }
        if (isTRUE(use_omnibus)) {
            for (cpg in names(fit_list)) {
                model_obj <- fit_list[[cpg]]
                if (is.null(model_obj) || inherits(
                    model_obj,
                    "dnaEPICO_methylationLME_fit_error"
                ) || is.null(model_obj$omnibus) ||
                    !identical(model_obj$omnibus$status, "tested") ||
                    !is.finite(model_obj$omnibus$pValue) ||
                    model_obj$omnibus$pValue >= threshold) {
                    next
                }
                fit_quality <- fitStatusValuesMethylationLME(model_obj)
                if ((isTRUE(excludeSingular) &&
                    isTRUE(fit_quality$singular)) ||
                    (isTRUE(excludeNonConverged) && identical(
                        fit_quality$converged,
                        FALSE
                    ))) {
                    next
                }
                if (!is.null(model_obj$coef)) {
                    phenotype_hits[[cpg]] <- as.data.frame(model_obj$coef)
                }
            }
            retained[[phenotype]] <- phenotype_hits
            next
        }
        if (!is.null(modelResults$summaryCache) &&
            !is.null(modelResults$summaryCache[[phenotype]]) &&
            optionalTermMatchesMethylationModels(
                requested = interactionTerm,
                cached = modelResults$settings$interactionTerm
            )) {
            cached_summary <- modelResults$summaryCache[[phenotype]]
            if (nrow(cached_summary) > 0L && !is.na(threshold)) {
                hit_cpgs <- unique(cached_summary$CpG[cached_summary$P.value <
                    threshold])
                hit_cpgs <- hit_cpgs[!is.na(hit_cpgs) & hit_cpgs %in%
                    names(fit_list)]
                for (cpg in hit_cpgs) {
                    model_obj <- fit_list[[cpg]]
                    if (is.null(model_obj) || inherits(
                        model_obj,
                        "dnaEPICO_methylationLME_fit_error"
                    )) {
                        next
                    }
                    fit_quality <- fitStatusValuesMethylationLME(model_obj)
                    if ((isTRUE(excludeSingular) &&
                        isTRUE(fit_quality$singular)) ||
                        (isTRUE(excludeNonConverged) && identical(
                            fit_quality$converged,
                            FALSE
                        ))) {
                        next
                    }
                    if (!is.null(model_obj$coef)) {
                        phenotype_hits[[cpg]] <- as.data.frame(model_obj$coef)
                    }
                }
            }

            retained[[phenotype]] <- phenotype_hits
            next
        }

        for (cpg in names(fit_list)) {
            model_obj <- fit_list[[cpg]]
            if (is.null(model_obj) || inherits(model_obj,
                "dnaEPICO_methylationLME_fit_error")) {
                next
            }
            fit_quality <- fitStatusValuesMethylationLME(model_obj)
            if ((isTRUE(excludeSingular) && isTRUE(fit_quality$singular)) ||
                (isTRUE(excludeNonConverged) && identical(
                    fit_quality$converged,
                    FALSE
                ))) {
                next
            }

            coef_table <- model_obj$coef
            if (is.null(coef_table)) {
                next
            }

            matched_rows <- findCoefficientRowsMethylationLME(
                coefNames = rownames(coef_table),
                phenotype = phenotype, interactionTerm = interactionTerm,
                coefficientTerms = model_obj$coefficientTerms
            )
            if (length(matched_rows) == 0L) {
                next
            }

            matched_pvals <- coef_table[matched_rows, "Pr(>|t|)",
                drop = TRUE
            ]
            if (any(matched_pvals < threshold, na.rm = TRUE)) {
                phenotype_hits[[cpg]] <- as.data.frame(coef_table)
            }
        }

        retained[[phenotype]] <- phenotype_hits
    }

    hit_counts <- vapply(retained, length, integer(1))
    emitLogMinfiEwasWater(
        c(
            "============================================================",
            paste(
                "Significant longitudinal terms retained at p <",
                threshold, ":"
            ), paste(names(hit_counts), hit_counts,
                sep = ": ", collapse = "; "
            ), "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(retained, class = "dnaEPICO_methylationLME_significant")
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
#'     modelSummaries = ex$modelSummaries,
#'     preparedData = ex$preparedData,
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
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
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    fdrThreshold <- validateProbabilityDnaEpico(
        fdrThreshold,
        "fdrThreshold"
    )
    padjmethod <- validatePAdjustmentMethodMethylationModels(padjmethod)
    summary_list <- if (!is.null(modelSummaries$diagnosticSummaries)) {
        modelSummaries$diagnosticSummaries
    } else {
        modelSummaries$summaries
    }
    diagnostic_mean <- diagnosticMeanMethylationModels(preparedData)
    plot_list <- list()
    inflation_factors <- list()
    saved_files <- list()

    for (phenotype in names(summary_list)) {
        summary_df <- summary_list[[phenotype]]
        if (is.null(summary_df) || nrow(summary_df) == 0L) {
            next
        }

        summary_df$FDR <- adjustPvaluesByTermMethylationModels(
            pValues = summary_df$P.value,
            terms = summary_df$Interaction.Term, method = padjmethod
        )
        interaction_terms <- unique(as.character(summary_df$Interaction.Term))
        interaction_terms <- interaction_terms[!is.na(interaction_terms)]
        multiple_terms <- length(interaction_terms) > 1L
        phenotype_plots <- phenotype_files <- list()
        phenotype_inflation <- numeric(0)

        for (term_index in seq_along(interaction_terms)) {
            term <- interaction_terms[[term_index]]
            term_diagnostics <- buildMethylationTermDiagnosticsDnaEpico(
                summaryData = summary_df,
                phenotype = phenotype, term = term,
                    termColumn = "Interaction.Term",
                pValueColumn = "P.value", yColumn = "Std.Error",
                yLabel = "Standard Error", diagnosticMean = diagnostic_mean,
                fdrThreshold = fdrThreshold
            )
            if (is.null(term_diagnostics)) {
                next
            }

            file_key <- phenotype
            if (multiple_terms) {
                file_key <- paste0(phenotype, "_", sprintf(
                    "%02d",
                    term_index
                ), "_", sanitizeDiagnosticTermDnaEpico(term))
            }
            term_files <- list(
                qqplot = NULL, residualSD = NULL,
                residualSignificance = NULL
            )
            if (!is.null(outputDir)) {
                term_files <- list(
                    qqplot = file.path(
                        outputDir,
                        paste0("qqplot_", file_key, ".tiff")
                    ), residualSD = file.path(
                        outputDir,
                        paste0("residualSD_", file_key, ".tiff")
                    ),
                    residualSignificance = file.path(
                        outputDir,
                        paste0(
                            "residualSignificance_", file_key,
                            ".tiff"
                        )
                    )
                )
            }

            for (plot_name in names(term_diagnostics$plots)) {
                plot_object <- term_diagnostics$plots[[plot_name]]
                if (is.null(plot_object)) {
                    next
                }
                runPlotMinfiEwasWater(
                    draw_fun = function() drawPlotObjectMinfiEwasWater(plot_object),
                    display = display, file = term_files[[plot_name]],
                    width = plotWidth, height = plotHeight, res = plotDPI
                )
            }

            phenotype_inflation[[term]] <- term_diagnostics$lambda
            if (multiple_terms) {
                phenotype_plots[[term]] <- term_diagnostics$plots
                phenotype_files[[term]] <- term_files
            } else {
                phenotype_plots <- term_diagnostics$plots
                phenotype_files <- term_files
            }
        }

        if (length(phenotype_inflation) > 0L) {
            inflation_factors[[phenotype]] <- if (multiple_terms) {
                phenotype_inflation
            } else {
                unname(phenotype_inflation[[1L]])
            }
            plot_list[[phenotype]] <- phenotype_plots
            saved_files[[phenotype]] <- phenotype_files
        }
    }

    emitLogMinfiEwasWater(
        c(
            "============================================================",
            paste("Diagnostic plots generated for phenotypes:",
                length(plot_list)),
            if (is.null(outputDir)) {
                "Diagnostic plots were returned in memory only."
            } else {
                paste("Diagnostic plots saved to:    ", outputDir)
            }, "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(
        plots = plot_list, inflationFactors = inflation_factors,
        files = saved_files
    ), class = "dnaEPICO_methylationLME_diagnostic_plots")
}

resolveFitDiagnosticsMethylationLME <- function(
    modelSummaries,
    summaryList
) {
    required_columns <- c(
        "Phenotype", "CpG", "Fit.Status", "Singular.Fit",
        "Converged", "Convergence.Message", "Fit.Warning", "Inference.Included",
        "Exclusion.Reason"
    )
    rows <- list()
    row_index <- 1L
    for (phenotype in names(summaryList)) {
        summary_df <- summaryList[[phenotype]]
        if (is.null(summary_df) || nrow(summary_df) == 0L ||
            !("CpG" %in% names(summary_df))) {
            next
        }
        value_or_default <- function(column, default) {
            if (column %in% names(summary_df)) {
                summary_df[[column]]
            } else {
                rep(default, nrow(summary_df))
            }
        }
        phenotype_rows <- unique(data.frame(
            Phenotype = phenotype,
            CpG = as.character(summary_df$CpG),
                Fit.Status = as.character(value_or_default(
                "Fit.Status",
                "fitted"
            )), Singular.Fit = as.logical(value_or_default(
                "Singular.Fit",
                NA
            )), Converged = as.logical(value_or_default(
                "Converged",
                TRUE
            )), Convergence.Message = as.character(value_or_default(
                "Convergence.Message",
                NA_character_
            )), Fit.Warning = as.character(value_or_default(
                "Fit.Warning",
                NA_character_
            )), Inference.Included = as.logical(value_or_default(
                "Inference.Included",
                TRUE
            )), Exclusion.Reason = as.character(value_or_default(
                "Exclusion.Reason",
                NA_character_
            )), stringsAsFactors = FALSE, check.names = FALSE
        ))
        rows[[row_index]] <- phenotype_rows
        row_index <- row_index + 1L
    }

    diagnostics <- modelSummaries$fitDiagnostics
    if (is.data.frame(diagnostics) && all(c(
        "Phenotype", "CpG",
        "Fit.Status"
    ) %in% names(diagnostics))) {
        defaults <- list(
            Singular.Fit = NA, Converged = NA,
                Convergence.Message = NA_character_,
            Fit.Warning = NA_character_, Inference.Included = startsWith(
                as.character(diagnostics$Fit.Status),
                "fitted"
            ), Exclusion.Reason = NA_character_
        )
        for (column in names(defaults)) {
            if (!(column %in% names(diagnostics))) {
                value <- defaults[[column]]
                diagnostics[[column]] <- if (length(value) ==
                    1L) {
                    rep(value, nrow(diagnostics))
                } else {
                    value
                }
            }
        }
        rows[[row_index]] <- diagnostics[, required_columns,
            drop = FALSE
        ]
        row_index <- row_index + 1L
    }

    fit_failures <- modelSummaries$fitFailures
    if (is.data.frame(fit_failures) && nrow(fit_failures) > 0L) {
        rows[[row_index]] <- data.frame(
            Phenotype = as.character(fit_failures$Phenotype),
            CpG = as.character(fit_failures$CpG),
                Fit.Status = as.character(fit_failures$Status),
            Singular.Fit = NA, Converged = ifelse(fit_failures$Status ==
                "not_converged", FALSE, NA),
                    Convergence.Message = ifelse(fit_failures$Status ==
                "not_converged", as.character(fit_failures$Error),
            NA_character_
            ), Fit.Warning = NA_character_,
            Inference.Included = FALSE,
                Exclusion.Reason = as.character(fit_failures$Error),
            stringsAsFactors = FALSE, check.names = FALSE
        )
    }

    if (length(rows) == 0L) {
        return(data.frame(
            Phenotype = character(0), CpG = character(0),
            Fit.Status = character(0), Singular.Fit = logical(0),
            Converged = logical(0), Convergence.Message = character(0),
            Fit.Warning = character(0), Inference.Included = logical(0),
            Exclusion.Reason = character(0), stringsAsFactors = FALSE,
            check.names = FALSE
        ))
    }

    diagnostics <- do.call(rbind, rows)
    key <- paste(diagnostics$Phenotype, diagnostics$CpG, sep = "\r")
    diagnostics <- diagnostics[!duplicated(key), , drop = FALSE]
    rownames(diagnostics) <- NULL
    diagnostics
}

buildAnnotationFitDiagnosticsMethylationLME <- function(fitDiagnostics) {
    if (!is.data.frame(fitDiagnostics) || nrow(fitDiagnostics) ==
        0L) {
        return(list(overall = data.frame(), byPhenotype = data.frame()))
    }

    aggregate_logical <- function(values, falseDominates = TRUE) {
        values <- as.logical(values)
        values <- values[!is.na(values)]
        if (length(values) == 0L) {
            return(NA)
        }
        if (isTRUE(falseDominates) && any(!values)) {
            return(FALSE)
        }
        if (!isTRUE(falseDominates) && any(values)) {
            return(TRUE)
        }
        if (isTRUE(falseDominates)) {
            TRUE
        } else {
            FALSE
        }
    }
    diagnostics_by_cpg <- split(fitDiagnostics, fitDiagnostics$CpG,
        drop = TRUE
    )
    overall_rows <- lapply(names(diagnostics_by_cpg), function(cpg) {
        diagnostic <- diagnostics_by_cpg[[cpg]]
        statuses <- as.character(diagnostic$Fit.Status)
        failed <- statuses %in% c("failed", "invalid", "not_converged")
        status <- if (any(failed)) {
            if (!all(failed)) {
                "partial"
            } else if (all(statuses == "invalid")) {
                "invalid"
            } else if (all(statuses == "not_converged")) {
                "not_converged"
            } else {
                "failed"
            }
        } else if (any(statuses == "fitted_not_converged")) {
            "fitted_not_converged"
        } else if (any(statuses == "fitted_singular")) {
            "fitted_singular"
        } else if (any(statuses == "fitted_with_warning")) {
            "fitted_with_warning"
        } else {
            "fitted"
        }
        reason_index <- !is.na(diagnostic$Exclusion.Reason) &
            nzchar(diagnostic$Exclusion.Reason)
        reasons <- unique(paste0(
            diagnostic$Phenotype[reason_index],
            ": ", diagnostic$Exclusion.Reason[reason_index]
        ))
        data.frame(
            CpG = cpg, Diagnostic.Fit.Status = status,
            Singular.Fit = aggregate_logical(diagnostic$Singular.Fit,
                falseDominates = FALSE
            ), Converged = aggregate_logical(diagnostic$Converged),
            Inference.Included = aggregate_logical(diagnostic$Inference.Included),
            Fit.Quality.Reason = if (length(reasons)) {
                paste(reasons, collapse = " | ")
            } else {
                NA_character_
            }, stringsAsFactors = FALSE, check.names = FALSE
        )
    })
    overall <- do.call(rbind, overall_rows)
    rownames(overall) <- NULL

    diagnostic_columns <- setdiff(names(fitDiagnostics), c(
        "Phenotype",
        "CpG"
    ))
    phenotype_tables <- lapply(
        unique(as.character(fitDiagnostics$Phenotype)),
        function(phenotype) {
            diagnostic <- fitDiagnostics[fitDiagnostics$Phenotype ==
                phenotype, c("CpG", diagnostic_columns), drop = FALSE]
            diagnostic <- unique(diagnostic)
            colnames(diagnostic)[-1L] <- paste0(
                phenotype, "_",
                diagnostic_columns
            )
            diagnostic
        }
    )
    phenotype_tables <- Filter(Negate(is.null), phenotype_tables)
    by_phenotype <- if (length(phenotype_tables) == 0L) {
        data.frame()
    } else if (length(phenotype_tables) == 1L) {
        phenotype_tables[[1L]]
    } else {
        Reduce(
            function(x, y) merge(x, y, by = "CpG", all = TRUE),
            phenotype_tables
        )
    }

    list(overall = overall, byPhenotype = by_phenotype)
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
#'     modelSummaries = ex$modelSummaries,
#'     annotationObject = ex$annotationData,
#'     annotationCols = "Name,chr,pos",
#'     verbose = FALSE,
#'     logs = FALSE
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
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_methylationLME.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    fit_failures <- modelSummaries$fitFailures
    summary_phenotypes <- modelSummaries$phenotypes
    summary_list <- modelSummaries
    if (!is.null(modelSummaries$diagnosticSummaries)) {
        summary_list <- modelSummaries$diagnosticSummaries
    } else if (!is.null(modelSummaries$summaries)) {
        summary_list <- modelSummaries$summaries
    }
    if (is.null(summary_phenotypes)) {
        summary_phenotypes <- names(summary_list)
    }
    fit_diagnostics <- resolveFitDiagnosticsMethylationLME(
        modelSummaries = modelSummaries,
        summaryList = summary_list
    )
    annotation_diagnostics <-
        buildAnnotationFitDiagnosticsMethylationLME(fitDiagnostics = fit_diagnostics)

    annotation_cols <- splitOptionMinfiEwasWater(annotationCols,
        sep = ","
    )
    annotation_df <- coerceAnnotationDataMethylationGLM(annotationObject)

    cleaned_summaries <- lapply(names(summary_list), function(model_name) {
        df <- summary_list[[model_name]]
        if (is.null(df) || nrow(df) == 0L) {
            return(NULL)
        }
        if (!all(c("CpG", "Interaction.Term", "P.value") %in%
            colnames(df))) {
            return(NULL)
        }

        df_split <- split(df, df$Interaction.Term)
        model_tables <- lapply(names(df_split), function(term) {
            sub_df <- df_split[[term]][, c("CpG", "P.value"),
                drop = FALSE
            ]
            clean_term <- gsub("`", "", term, fixed = TRUE)
            interaction_suffix <- clean_term
            if (startsWith(clean_term, paste0(model_name, "."))) {
                interaction_suffix <- sub(paste0(
                    "^", escapeRegexMethylationGLM(model_name),
                    "\\."
                ), "", clean_term)
            }
            p_col <- paste0(
                model_name, "_", interaction_suffix,
                "_P.Value"
            )
            colnames(sub_df)[2] <- p_col
            sub_df
        })

        Reduce(
            function(x, y) merge(x, y, by = "CpG", all = TRUE),
            model_tables
        )
    })
    cleaned_summaries <- Filter(Negate(is.null), cleaned_summaries)

    if (length(cleaned_summaries) == 0L) {
        merged_summary <- data.frame(CpG = character(0))
    } else if (length(cleaned_summaries) == 1L) {
        merged_summary <- cleaned_summaries[[1L]]
    } else {
        merged_summary <- Reduce(function(x, y) {
            merge(x, y,
                by = "CpG",
                all = TRUE
            )
        }, cleaned_summaries)
    }
    omnibus_tables <- buildAnnotationOmnibusTablesMethylationLME(
        modelSummaries
    )
    for (omnibus_table in omnibus_tables) {
        merged_summary <- merge(
            merged_summary, omnibus_table,
            by = "CpG", all = TRUE
        )
    }
    lme_engine <- modelSummaries$settings$lmeEngine
    if (is.null(lme_engine) || !nzchar(lme_engine)) {
        lme_engine <- "lme4"
    }
    if (nrow(annotation_diagnostics$byPhenotype) > 0L) {
        lme_diagnostics <- annotation_diagnostics$byPhenotype
        if (!identical(lme_engine, "lme4")) {
            lme_diagnostics <- lme_diagnostics[, !grepl(
                "_Singular\\.Fit$",
                colnames(lme_diagnostics)
            ), drop = FALSE]
        }
        merged_summary <- merge(merged_summary, lme_diagnostics,
            by = "CpG", all = TRUE
        )
    }

    available_annotation_cols <- intersect(annotation_cols,
        colnames(annotation_df))
    missing_annotation_cols <- setdiff(annotation_cols, colnames(annotation_df))
    annotated_results <- merge(merged_summary, annotation_df[,
        c("CpG", available_annotation_cols),
        drop = FALSE
    ],
    by = "CpG",
    all.x = TRUE
    )
    if ("CpG" %in% colnames(annotated_results)) {
        colnames(annotated_results)[colnames(annotated_results) ==
            "CpG"] <- "IlmnID"
    }
    annotated_results <- orderAnnotatedModelColumnsDnaEpico(
        data = annotated_results,
        annotationCols = available_annotation_cols,
        modelNames = summary_phenotypes,
        includeSingular = identical(
            lme_engine,
            "lme4"
        )
    )

    emitLogMinfiEwasWater(
        c(
            "============================================================",
            paste("Annotated CpG rows:          ", nrow(annotated_results)),
            paste("Annotation columns used:      ",
                paste(available_annotation_cols,
                collapse = ", "
            )), if (length(missing_annotation_cols) ==
                0L) {
                "Missing annotation columns:   none"
            } else {
                paste("Missing annotation columns:   ",
                    paste(missing_annotation_cols,
                    collapse = ", "
                ))
            }, "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(
        data = annotated_results, fitFailures = fit_failures,
        fitDiagnostics = fit_diagnostics,
            annotationColumnsUsed = available_annotation_cols,
        missingAnnotationCols = missing_annotation_cols
    ), class = "dnaEPICO_methylationLME_annotation")
}

#' Write optional disk outputs for longitudinal mixed-effects analyses
#'
#' @param modelResults Object returned by `fitMethylationLMEModels()`.
#' @param modelSummaries Object returned by `summarizeMethylationLMEModels()`.
#' @param annotatedResults Object returned by
#'   `annotateMethylationLMESummaries()` or a compatible data frame.
#' @param significantInteractions Object returned by
#'   `collectSignificantInteractionsMethylationLME()` or `NULL`.
#' @param outputRData Character. Directory used for serialized model and summary
#'   outputs.
#' @param summaryTxtDir Character. Directory used for tab-delimited summary
#'   tables.
#' @param significantInteractionDir Character. Directory used for significant
#'   interaction coefficient tables.
#' @param annotatedLMEOut Character. Directory used for the annotated summary
#'   XLSX workbook.
#' @param reportAssetsDir Character or `NULL`. Report results directory used for
#'   the compressed TSV table and compact metadata sidecars.
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
#' Write optional serialized outputs, summary tables, significant interaction
#' tables, and annotated results from the longitudinal mixed-effects workflow.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' annotation_data <- annotateMethylationLMESummaries(
#'     modelSummaries = ex$modelSummaries,
#'     annotationObject = ex$annotationData,
#'     annotationCols = "Name,chr,pos",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' significant_hits <- collectSignificantInteractionsMethylationLME(
#'     modelResults = ex$modelResults,
#'     pvalThreshold = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' output_paths <- writeMethylationLMEOutputs(
#'     modelResults = ex$modelResults,
#'     modelSummaries = ex$modelSummaries,
#'     annotatedResults = annotation_data,
#'     significantInteractions = significant_hits,
#'     outputRData = file.path(ex$tempDir, "models"),
#'     summaryTxtDir = file.path(ex$tempDir, "summary"),
#'     significantInteractionDir = file.path(ex$tempDir, "significant"),
#'     annotatedLMEOut = file.path(ex$tempDir, "annotated"),
#'     saveTxtSummaries = TRUE,
#'     saveSignificantInteractions = TRUE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(output_paths)
#'
#' @export
writeMethylationLMEOutputs <- function(
    modelResults, modelSummaries,
    annotatedResults, significantInteractions = NULL, outputRData,
    summaryTxtDir, significantInteractionDir, annotatedLMEOut,
    reportAssetsDir = NULL, saveTxtSummaries = TRUE,
    saveSignificantInteractions = FALSE,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_methylationLME.txt"
) {
    if (!is.null(reportAssetsDir) &&
        (!is.character(reportAssetsDir) || length(reportAssetsDir) != 1L ||
            is.na(reportAssetsDir) || !nzchar(reportAssetsDir))) {
        stop("reportAssetsDir must be NULL or one non-empty directory path.",
            call. = FALSE
        )
    }
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    dir.create(outputRData, recursive = TRUE, showWarnings = FALSE)
    dir.create(annotatedLMEOut, recursive = TRUE, showWarnings = FALSE)

    model_files <- stats::setNames(
        character(length(modelResults$fits)),
        names(modelResults$fits)
    )
    summary_files <- stats::setNames(
        character(length(modelSummaries$summaries)),
        names(modelSummaries$summaries)
    )
    summary_txt_files <- list()
    significant_files <- list()
    fit_failure_file <- character(0)

    for (phenotype in names(modelResults$fits)) {
        model_file <- file.path(outputRData, paste0(
            phenotype,
            "LME.rds"
        ))
        summary_file <- file.path(outputRData, paste0(
            phenotype,
            "SummaryLME.rds"
        ))
        saveRDS(modelResults$fits[[phenotype]], file = model_file)
        saveRDS(modelSummaries$summaries[[phenotype]], file = summary_file)
        model_files[[phenotype]] <- model_file
        summary_files[[phenotype]] <- summary_file
    }

    if (isTRUE(saveTxtSummaries)) {
        dir.create(summaryTxtDir, recursive = TRUE, showWarnings = FALSE)
        for (phenotype in names(modelSummaries$summaries)) {
            summary_df <- modelSummaries$summaries[[phenotype]]
            if (is.null(summary_df) || nrow(summary_df) == 0L) {
                next
            }

            if ("P.value" %in% colnames(summary_df)) {
                summary_df <- summary_df[order(summary_df$P.value), ,
                    drop = FALSE
                ]
            }

            output_file <- file.path(summaryTxtDir, paste0(
                phenotype,
                "SummaryLME.txt"
            ))
            utils::write.table(summary_df,
                file = output_file,
                sep = "\t", row.names = FALSE, quote = FALSE
            )
            summary_txt_files[[phenotype]] <- output_file
        }
    }

    if (isTRUE(saveSignificantInteractions) &&
        !is.null(significantInteractions)) {
        dir.create(significantInteractionDir,
            recursive = TRUE,
            showWarnings = FALSE
        )
        for (phenotype in names(significantInteractions)) {
            phenotype_hits <- significantInteractions[[phenotype]]
            if (length(phenotype_hits) == 0L) {
                next
            }

            significant_files[[phenotype]] <- character(0)
            phenotype_dir <- file.path(
                significantInteractionDir,
                phenotype
            )
            dir.create(phenotype_dir, recursive = TRUE, showWarnings = FALSE)

            for (cpg in names(phenotype_hits)) {
                cpg_dir <- file.path(phenotype_dir, cpg)
                dir.create(cpg_dir, recursive = TRUE, showWarnings = FALSE)
                output_file <- file.path(cpg_dir, paste0(
                    cpg,
                    ".txt"
                ))
                utils::write.table(phenotype_hits[[cpg]],
                    file = output_file,
                    sep = "\t", quote = FALSE
                )
                significant_files[[phenotype]] <- c(
                    significant_files[[phenotype]],
                    output_file
                )
            }
        }
    }

    annotated_df <- annotatedResults
    if (!is.null(annotatedResults$data)) {
        annotated_df <- annotatedResults$data
    }
    report_table <- sortReportTableDnaEpico(annotated_df)
    annotated_df <- report_table$data
    annotated_file <- file.path(annotatedLMEOut, "annotatedLME.xlsx")
    fit_failures <- modelSummaries$fitFailures
    if (is.null(fit_failures)) {
        fit_failures <- modelResults$fitFailures
    }
    if (is.data.frame(fit_failures) && nrow(fit_failures) > 0L) {
        fit_failure_file <- file.path(annotatedLMEOut, "CpGFitFailuresLME.txt")
        utils::write.table(fit_failures,
            file = fit_failure_file,
            sep = "\t", row.names = FALSE, quote = FALSE
        )
    }
    dictionary <- buildAnnotatedWorkbookDictionaryMethylationGLM(
        columns = colnames(annotated_df),
        modelDescription = "Pvalue from LME model",
            formulaText = modelResults$formulas,
        modelLabel = "LME",
            responseLabel = inferMethylationValueLabelMethylationGLM(modelResults)
    )
    metadata <- buildModelWorkbookMetadataDnaEpico(
        modelResults = modelResults,
        modelSummaries = modelSummaries, annotatedResults = annotatedResults,
        analysis = "lme"
    )
    writeAnnotatedWorkbookMethylationGLM(
        annotated_df = annotated_df,
        file = annotated_file, resultSheet = "annotatedLME",
        dictionary = dictionary, metadata = metadata
    )
    report_sidecar <- list(
        table = NULL, metadata = NULL, dictionary = NULL,
        workbookMetadata = NULL
    )
    if (!is.null(reportAssetsDir)) {
        report_sidecar <- writeReportTableSidecarDnaEpico(
            tableData = annotated_df,
            workbookFile = annotated_file, sidecarDir = reportAssetsDir,
            sheet = "annotatedLME", idColumn = report_table$idColumn,
            dictionary = dictionary, workbookMetadata = metadata
        )
    }

    emitLogMinfiEwasWater(
        c(
            "============================================================",
            paste("Serialized model outputs:    ", length(model_files)),
            paste("Serialized summary outputs:  ", length(summary_files)),
            paste("Annotated results file:      ", annotated_file),
            if (is.null(report_sidecar$table)) {
                "Report table sidecar:          not requested"
            } else {
                paste("Report table sidecar:         ", report_sidecar$table)
            },
            paste("CpG fit failures reported:   ", nrow(fit_failures)),
            if (isTRUE(saveTxtSummaries)) {
                paste("Summary text files written:   ",
                    length(summary_txt_files))
            } else {
                "Summary text files written:   0"
            }, if (isTRUE(saveSignificantInteractions)) {
                paste("Significant interaction files:", sum(vapply(
                    significant_files,
                    length, integer(1)
                )))
            } else {
                "Significant interaction files: 0"
            }, "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(
        list(
            modelFiles = model_files, summaryFiles = summary_files,
            summaryTxtFiles = summary_txt_files,
                significantInteractionFiles = significant_files,
            fitFailureFile = fit_failure_file, annotatedLME = annotated_file,
            annotatedLMEText = report_sidecar$table,
                annotatedLMEReportMetadata = report_sidecar$metadata,
            annotatedLMEDictionary = report_sidecar$dictionary,
                annotatedLMEMetadata = report_sidecar$workbookMetadata
        ),
        class = "dnaEPICO_methylationLME_paths"
    )
}
