#' Normalize optional numeric inputs for methylationGLM helpers
#'
#' @param value Numeric, character, or `NA` input.
#'
#' @return A numeric scalar or `NA_real_`.
#'
#' @description
#' Internal helper that standardizes optional numeric arguments such as
#' `cpgLimit` and `summaryPval`.
#'
#' @keywords internal
#' @noRd
normalizeOptionalNumericMethylationGLM <- function(value) {
    if (is.null(value) || length(value) == 0L) {
    return(NA_real_)
    }
    if (length(value) != 1L) {
    stop("Optional numeric inputs must contain one value.",
        call. = FALSE
    )
    }

    if (is.character(value)) {
    trimmed <- trimws(value)
    if (is.na(trimmed)) {
        return(NA_real_)
    }
    if (!nzchar(trimmed) || tolower(trimmed) %in% c(
        "na",
        "null"
    )) {
        return(NA_real_)
    }
    numeric_value <- utils::type.convert(trimmed,
        as.is = TRUE,
        na.strings = character()
    )
    if (!is.numeric(numeric_value) || is.na(numeric_value)) {
        stop("Optional numeric inputs must be numeric, NA, or NULL.",
        call. = FALSE
        )
    }
    return(numeric_value)
    }

    if (is.na(value)) {
    return(NA_real_)
    }
    if (!is.numeric(value)) {
    stop("Optional numeric inputs must be numeric, NA, or NULL.",
        call. = FALSE
    )
    }

    as.numeric(value)
}

validateCpgLimitMethylationModels <- function(value) {
    if (is.na(value)) {
    return(NA_integer_)
    }
    if (!is.finite(value) || value < 1 || value != floor(value)) {
    stop("cpgLimit must be NA or a positive whole number.",
        call. = FALSE
    )
    }
    as.integer(value)
}

validatePositiveIntegerMethylationModels <- function(value, name) {
    if (length(value) != 1L) {
    stop(name, " must be a positive whole number.", call. = FALSE)
    }
    numeric_value <- utils::type.convert(as.character(value),
    as.is = TRUE,
    na.strings = character()
    )
    if (!is.numeric(numeric_value) || is.na(numeric_value) ||
    !is.finite(numeric_value) || numeric_value < 1 || numeric_value !=
    floor(numeric_value)) {
    stop(name, " must be a positive whole number.", call. = FALSE)
    }
    as.integer(numeric_value)
}

validatePAdjustmentMethodMethylationModels <- function(method) {
    if (length(method) != 1L || is.na(method) || !(method %in%
    stats::p.adjust.methods)) {
    supported_methods_text <- paste(stats::p.adjust.methods,
        collapse = ", "
    )
    stop(sprintf(
        "padjmethod must be one of: %s", supported_methods_text
    ), call. = FALSE)
    }
    method
}

normalizeOptionalColumnMethylationModels <- function(value, name) {
    if (is.null(value) || length(value) == 0L || all(is.na(value))) {
    return(NULL)
    }
    if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop(name, " must be NULL or one column name.", call. = FALSE)
    }
    value <- trimws(value)
    if (!nzchar(value) || tolower(value) %in% c("null", "na")) {
    return(NULL)
    }
    value
}

adjustPvaluesByTermMethylationModels <- function(
    pValues, terms,
    method
) {
    adjusted <- rep(NA_real_, length(pValues))
    groups <- split(seq_along(pValues), as.character(terms),
    drop = TRUE
    )
    for (indices in groups) {
    adjusted[indices] <- stats::p.adjust(pValues[indices],
        method = method
    )
    }
    adjusted
}

combineModelMessagesDnaEpico <- function(...) {
    messages <- unlist(list(...), use.names = FALSE)
    messages <- trimws(as.character(messages))
    messages <- unique(messages[!is.na(messages) & nzchar(messages)])
    if (length(messages) == 0L) {
    return(NA_character_)
    }
    paste(messages, collapse = " || ")
}

captureModelConditionsDnaEpico <- function(expression) {
    state <- new.env(parent = emptyenv())
    state$conditions <- character(0)
    state$error <- NULL
    value <- tryCatch(
    withCallingHandlers(
        expression,
        message = function(condition) {
        state$conditions <- c(
            state$conditions,
            paste0("MESSAGE: ", conditionMessage(condition))
        )
        invokeRestart("muffleMessage")
        },
        warning = function(condition) {
        state$conditions <- c(
            state$conditions,
            paste0("WARNING: ", conditionMessage(condition))
        )
        invokeRestart("muffleWarning")
        }
    ),
    error = function(condition) {
        state$error <- condition
        state$conditions <- c(
        state$conditions,
        paste0("ERROR: ", conditionMessage(condition))
        )
        NULL
    }
    )
    list(
    value = value,
    error = state$error,
    modelMessage = combineModelMessagesDnaEpico(state$conditions)
    )
}

modelMessageDnaEpico <- function(modelObject) {
    if (is.null(modelObject$modelMessage)) {
    return(NA_character_)
    }
    combineModelMessagesDnaEpico(modelObject$modelMessage)
}

collectFitFailuresMethylationModels <- function(fits, errorClass) {
    rows <- list()
    row_index <- 1L
    for (phenotype in names(fits)) {
    fit_group <- fits[[phenotype]]
    for (cpg in names(fit_group)) {
        fit_object <- fit_group[[cpg]]
        if (!inherits(fit_object, errorClass)) {
        next
        }
        rows[[row_index]] <- data.frame(
        Phenotype = phenotype,
        CpG = cpg, Status = "failed",
        Error = as.character(fit_object$error),
        stringsAsFactors = FALSE, check.names = FALSE
        )
        row_index <- row_index + 1L
    }
    }
    if (length(rows) == 0L) {
    return(data.frame(
        Phenotype = character(0), CpG = character(0),
        Status = character(0), Error = character(0),
        stringsAsFactors = FALSE,
        check.names = FALSE
    ))
    }
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out
}

diagnosticMeanMethylationModels <- function(preparedData) {
    methylation_matrix <- preparedData$data[, preparedData$cpgColumns,
    drop = FALSE
    ]
    means <- vapply(
    methylation_matrix, meanFiniteOrNADnaEpico,
    numeric(1)
    )
    scale <- preparedData$methylationScale
    if (is.null(scale) || identical(scale, "beta")) {
    return(list(values = means, label = "Average Beta"))
    }
    response_label <- preparedData$responseLabel
    if (is.null(response_label) || !nzchar(response_label)) {
    response_label <- toupper(scale)
    }
    list(values = means, label = paste("Average", response_label))
}

sanitizeDiagnosticTermDnaEpico <- function(term) {
    component <- gsub("[^A-Za-z0-9.-]+", "_", as.character(term))
    component <- gsub("^_+|_+$", "", component)
    if (nzchar(component)) {
    component
    } else {
    "term"
    }
}

prepareTermDiagnosticDataDnaEpico <- function(
    summaryData, term, termColumn, pValueColumn, maximumPoints
) {
    term_values <- as.character(summaryData[[termColumn]])
    data <- summaryData[
    !is.na(term_values) & term_values == term, ,
    drop = FALSE
    ]
    valid <- is.finite(data[[pValueColumn]]) & data[[pValueColumn]] >= 0 &
    data[[pValueColumn]] <= 1
    data <- data[valid, , drop = FALSE]
    if (!nrow(data)) {
    return(NULL)
    }
    p_values <- pmax(data[[pValueColumn]], .Machine$double.xmin)
    chi_square <- stats::qchisq(p_values, df = 1, lower.tail = FALSE)
    lambda <- if (length(p_values) >= 100L) {
    round(
        stats::median(chi_square) / stats::qchisq(0.5, df = 1), 3
    )
    } else {
    NA_real_
    }
    index <- seq_along(p_values)
    qq <- data.frame(
    expected = -log10(index / (length(p_values) + 1)),
    observed = -log10(sort(p_values)),
    lower = -log10(stats::qbeta(
        0.975, index, length(p_values) - index + 1L
    )),
    upper = -log10(stats::qbeta(
        0.025, index, length(p_values) - index + 1L
    ))
    )
    rows <- deterministicPlotRowsDnaEpico(
    nrow(qq),
    maximum = maximumPoints,
    priority = utils::head(order(qq$observed, decreasing = TRUE), 500L)
    )
    list(
    data = data, pValues = p_values, lambda = lambda,
    qq = qq[rows, , drop = FALSE]
    )
}

createDiagnosticQqPlotDnaEpico <- function(data) {
    style <- adaptivePointStyleDnaEpico(nrow(data))
    ggplot2::ggplot(data, ggplot2::aes(x = expected, y = observed)) +
    ggplot2::geom_ribbon(
        ggplot2::aes(ymin = lower, ymax = upper),
        fill = "#9ECAE1", alpha = 0.35
    ) +
    ggplot2::geom_point(
        color = "#111827", size = style$size, alpha = style$alpha
    ) +
    ggplot2::geom_abline(
        intercept = 0, slope = 1, color = "#B42318", linewidth = 0.65
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
        title = NULL, x = "Expected -log10(p)",
        y = "Observed -log10(p)"
    ) +
    dnaEpicoModelPlotTheme()
}

createDiagnosticMeanPlotDnaEpico <- function(data, meanInfo, yLabel) {
    style <- adaptivePointStyleDnaEpico(nrow(data))
    plot <- ggplot2::ggplot(
    data, ggplot2::aes(x = meanMethylation, y = diagnosticY)
    )
    if (nrow(data) > 5000L) {
    plot <- plot + ggplot2::geom_bin_2d(bins = 55) +
        ggplot2::scale_fill_gradient(
        low = "#D8EFF3", high = "#176B87", name = "CpGs"
        )
    } else {
    plot <- plot + ggplot2::geom_point(
        alpha = style$alpha, size = style$size, color = "#243746"
    )
    }
    plot <- plot + ggplot2::labs(
    title = NULL, x = meanInfo$label, y = yLabel
    ) + dnaEpicoModelPlotTheme()
    if (nrow(data) >= 3L && nrow(data) <= 10000L) {
    plot <- plot + ggplot2::geom_smooth(
        method = "loess", formula = y ~ x, se = FALSE, color = "red"
    )
    }
    plot
}

createDiagnosticSignificancePlotDnaEpico <- function(
    data, labelData, style, yLabel, fdrThreshold
) {
    ggplot2::ggplot(data, ggplot2::aes(
    x = -log10(diagnosticP), y = diagnosticY,
    color = FDR < fdrThreshold
    )) +
    ggplot2::geom_point(alpha = style$alpha, size = style$size) +
    ggrepel::geom_text_repel(
        data = labelData, ggplot2::aes(label = CpG),
        size = 3, max.overlaps = Inf, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = c(
        `FALSE` = "grey50", `TRUE` = "firebrick"
    )) +
    ggplot2::labs(
        title = NULL, x = "-log10(p-value)", y = yLabel,
        color = paste("FDR <", fdrThreshold)
    ) +
    dnaEpicoModelPlotTheme()
}

buildDiagnosticYPlotsDnaEpico <- function(
    data, yColumn, yLabel, meanInfo, fdrThreshold, maximumPoints
) {
    if (!(yColumn %in% colnames(data))) {
    return(list(mean = NULL, significance = NULL))
    }
    data$diagnosticY <- data[[yColumn]]
    complete <- is.finite(data$meanMethylation) & is.finite(data$diagnosticY)
    diagnostic <- data[complete, , drop = FALSE]
    significant <- !is.na(diagnostic$FDR) & diagnostic$FDR < fdrThreshold
    rows <- deterministicPlotRowsDnaEpico(
    nrow(diagnostic),
    maximum = maximumPoints, priority = which(significant)
    )
    mean_plot <- createDiagnosticMeanPlotDnaEpico(
    diagnostic[rows, , drop = FALSE], meanInfo, yLabel
    )
    significant_data <- data[
    !is.na(data$FDR) & data$FDR < fdrThreshold, ,
    drop = FALSE
    ]
    rows <- deterministicPlotRowsDnaEpico(
    nrow(data),
    maximum = maximumPoints,
    priority = which(!is.na(data$FDR) & data$FDR < fdrThreshold)
    )
    display_data <- data[rows, , drop = FALSE]
    label_data <- utils::head(
    significant_data[order(significant_data$diagnosticP), , drop = FALSE],
    20L
    )
    significance_plot <- createDiagnosticSignificancePlotDnaEpico(
    display_data, label_data,
    adaptivePointStyleDnaEpico(nrow(display_data)),
    yLabel, fdrThreshold
    )
    list(mean = mean_plot, significance = significance_plot)
}

createDiagnosticVolcanoPlotDnaEpico <- function(data, fdrThreshold) {
    style <- adaptivePointStyleDnaEpico(nrow(data))
    ggplot2::ggplot(data, ggplot2::aes(
    x = diagnosticEstimate, y = minusLog10P, colour = significant
    )) +
    ggplot2::geom_point(size = style$size, alpha = style$alpha) +
    ggplot2::geom_vline(
        xintercept = 0, colour = "#64748B", linewidth = 0.4
    ) +
    ggplot2::scale_colour_manual(values = c(
        `FALSE` = "#94A3B8", `TRUE` = "#B42318"
    )) +
    ggplot2::labs(
        title = NULL, x = "Coefficient estimate",
        y = "-log10(p-value)", colour = paste("FDR <", fdrThreshold)
    ) +
    dnaEpicoModelPlotTheme()
}

createDiagnosticForestPlotDnaEpico <- function(
    data, standardErrorColumn, fdrThreshold
) {
    if (is.null(standardErrorColumn) || !(standardErrorColumn %in% names(
        data))) {
    return(NULL)
    }
    data$diagnosticSE <- coerceNumericDnaEpico(data[[standardErrorColumn]])
    data <- data[
    is.finite(data$diagnosticSE) & data$diagnosticSE >= 0, ,
    drop = FALSE
    ]
    data <- utils::head(
    data[order(data$diagnosticP, data$CpG), , drop = FALSE], 20L
    )
    if (!nrow(data)) {
    return(NULL)
    }
    rank_label <- paste0("CpG (top ", nrow(data), " ranked by p-value)")
    data$CpG <- factor(data$CpG, levels = rev(data$CpG))
    data$lower <- data$diagnosticEstimate - 1.96 * data$diagnosticSE
    data$upper <- data$diagnosticEstimate + 1.96 * data$diagnosticSE
    ggplot2::ggplot(data, ggplot2::aes(
    x = diagnosticEstimate, y = CpG, colour = significant
    )) +
    ggplot2::geom_vline(
        xintercept = 0, colour = "#64748B", linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
        ggplot2::aes(xmin = lower, xmax = upper),
        width = 0.18, linewidth = 0.65
    ) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_colour_manual(values = c(
        `FALSE` = "#176B87", `TRUE` = "#B42318"
    )) +
    ggplot2::labs(
        title = NULL, x = "Estimate with 95% CI", y = rank_label,
        colour = paste("FDR <", fdrThreshold)
    ) +
    dnaEpicoModelPlotTheme()
}

buildDiagnosticEffectPlotsDnaEpico <- function(
    data, pValues, estimateColumn, standardErrorColumn,
    fdrThreshold, maximumPoints
) {
    if (is.null(estimateColumn) || !(estimateColumn %in% names(data))) {
    return(list(volcano = NULL, forest = NULL))
    }
    estimates <- coerceNumericDnaEpico(data[[estimateColumn]])
    complete <- is.finite(estimates) & is.finite(pValues)
    effect <- data[complete, , drop = FALSE]
    effect$diagnosticEstimate <- estimates[complete]
    effect$minusLog10P <- -log10(pValues[complete])
    effect$significant <- !is.na(effect$FDR) & effect$FDR < fdrThreshold
    rows <- deterministicPlotRowsDnaEpico(
    nrow(effect),
    maximum = maximumPoints,
    priority = which(effect$significant)
    )
    list(
    volcano = createDiagnosticVolcanoPlotDnaEpico(
        effect[rows, , drop = FALSE], fdrThreshold
    ),
    forest = createDiagnosticForestPlotDnaEpico(
        effect, standardErrorColumn, fdrThreshold
    )
    )
}

buildMethylationTermDiagnosticsDnaEpico <- function(
    summaryData,
    phenotype, term, termColumn, pValueColumn, yColumn, yLabel,
    diagnosticMean, fdrThreshold, estimateColumn = NULL,
    standardErrorColumn = NULL, maximumPoints = 50000L
) {
    prepared <- prepareTermDiagnosticDataDnaEpico(
    summaryData, term, termColumn, pValueColumn, maximumPoints
    )
    if (is.null(prepared)) {
    return(NULL)
    }
    data <- prepared$data
    data$meanMethylation <- diagnosticMean$values[data$CpG]
    data$diagnosticP <- prepared$pValues
    y_plots <- buildDiagnosticYPlotsDnaEpico(
    data, yColumn, yLabel, diagnosticMean, fdrThreshold, maximumPoints
    )
    effect_plots <- buildDiagnosticEffectPlotsDnaEpico(
    data, prepared$pValues, estimateColumn, standardErrorColumn,
    fdrThreshold, maximumPoints
    )
    list(
    lambda = prepared$lambda,
    plots = list(
        qqplot = createDiagnosticQqPlotDnaEpico(prepared$qq),
        residualSD = y_plots$mean,
        residualSignificance = y_plots$significance,
        volcano = effect_plots$volcano,
        effectForest = effect_plots$forest
    )
    )
}

#' Normalize optional chunk-size inputs for methylationGLM helpers
#'
#' @param chunkSize Integer, character, or `NULL` input.
#'
#' @return An integer scalar or `NULL`.
#'
#' @description
#' Internal helper that standardizes optional chunk-size arguments used while
#' extracting CpG-level summaries.
#'
#' @keywords internal
#' @noRd
normalizeChunkSizeMethylationGLM <- function(chunkSize) {
    if (is.null(chunkSize) || length(chunkSize) == 0L) {
    return(NULL)
    }

    if (length(chunkSize) == 1L && is.character(chunkSize)) {
    trimmed <- trimws(chunkSize)
    if (!nzchar(trimmed) || tolower(trimmed) %in% c(
        "na",
        "null"
    )) {
        return(NULL)
    }
    return(validatePositiveIntegerMethylationModels(
        trimmed,
        "chunkSize"
    ))
    }

    validatePositiveIntegerMethylationModels(chunkSize, "chunkSize")
}

#' Escape regular-expression metacharacters for methylationGLM helpers
#'
#' @param x Character scalar to escape.
#'
#' @return Character scalar safe to use in `grep()` or `grepl()`.
#'
#' @description
#' Internal helper that escapes user-supplied model terms before matching
#' coefficient names.
#'
#' @keywords internal
#' @noRd
escapeRegexMethylationGLM <- function(x) {
    gsub("([][{}()+*^$|\\?.])", "\\\\\\1", x)
}

#' Backtick variable names for methylationGLM formulas
#'
#' @param x Character vector of variable names.
#'
#' @return Character vector with names wrapped in backticks.
#'
#' @description
#' Internal helper that keeps formula construction robust when variable names
#' contain punctuation.
#'
#' @keywords internal
#' @noRd
quoteNamesMethylationGLM <- function(x) {
    if (length(x) == 0L) {
    return(character(0))
    }
    paste0("`", gsub("`", "", x, fixed = TRUE), "`")
}

#' Parse phenotype-to-PRS mappings for methylationGLM helpers
#'
#' @param prsMap Character scalar or vector describing mappings in the form
#'   `'Phenotype:PRS'`.
#'
#' @return A named character vector.
#'
#' @description
#' Internal helper that parses the optional `prsMap` argument.
#'
#' @keywords internal
#' @noRd
parsePrsMapMethylationGLM <- function(prsMap = NULL) {
    if (is.null(prsMap) || length(prsMap) == 0L) {
    return(stats::setNames(character(0), character(0)))
    }

    map_entries <- splitOptionMinfiEwasWater(prsMap, sep = ",")
    if (length(map_entries) == 0L) {
    return(stats::setNames(character(0), character(0)))
    }

    pieces <- strsplit(map_entries, split = ":", fixed = TRUE)
    invalid <- !vapply(pieces, function(x) length(x) == 2L, logical(1))
    if (any(invalid)) {
    invalid_entries_text <- paste(map_entries[invalid], collapse = ", ")
    stop(
        sprintf(
        "%s Invalid entries: %s",
        "Each prsMap entry must follow the format 'Phenotype:PRS'.",
        invalid_entries_text
        ),
        call. = FALSE
    )
    }

    keys <- vapply(pieces, function(x) trimws(x[[1L]]), character(1))
    values <- vapply(pieces, function(x) trimws(x[[2L]]), character(1))
    if (any(!nzchar(keys)) || any(!nzchar(values))) {
    stop("prsMap phenotype and PRS column names cannot be blank.",
        call. = FALSE
    )
    }
    duplicated_keys <- unique(keys[duplicated(keys)])
    if (length(duplicated_keys) > 0L) {
    duplicated_keys_text <- paste(duplicated_keys, collapse = ", ")
    stop(
        sprintf(
        "prsMap contains duplicate phenotype mappings: %s",
        duplicated_keys_text
        ),
        call. = FALSE
    )
    }
    stats::setNames(values, keys)
}

#' Find coefficient rows for a phenotype or interaction term
#'
#' @param coefNames Character vector of model coefficient names.
#' @param variable Character. Phenotype variable of interest.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return Character vector of matching coefficient names.
#'
#' @description
#' Internal helper that locates the phenotype main-effect or interaction rows in
#' a CpG-level coefficient table.
#'
#' @keywords internal
#' @noRd
findCoefficientRowsMethylationGLM <- function(
    coefNames, variable,
    interactionTerm = NULL, coefficientTerms = NULL
) {
    if (!is.null(coefficientTerms)) {
    coefficient_names <- intersect(coefNames, names(coefficientTerms))
    term_values <- coefficientTerms[coefficient_names]
    normalize_term <- function(x) {
        gsub("`", "", gsub("[[:space:]]+", "", x), fixed = TRUE)
    }
    normalized_variable <- normalize_term(variable)
    normalized_terms <- normalize_term(term_values)

    if (!is.null(interactionTerm) && nzchar(interactionTerm)) {
        normalized_interaction <- normalize_term(interactionTerm)
        matches <- vapply(strsplit(normalized_terms, ":",
        fixed = TRUE
        ), function(parts) {
        normalized_variable %in% parts && normalized_interaction %in%
            parts
        }, logical(1))
    } else {
        matches <- normalized_terms == normalized_variable
    }

    return(coefficient_names[matches])
    }

    normalized_names <- gsub("`", "", coefNames, fixed = TRUE)
    variable_pattern <- escapeRegexMethylationGLM(variable)

    if (!is.null(interactionTerm) && nzchar(interactionTerm)) {
    interaction_pattern <- escapeRegexMethylationGLM(interactionTerm)
    matches <- grepl(paste0(
        "^", variable_pattern, ".*:",
        interaction_pattern
    ), normalized_names)
    } else {
    matches <- grepl(paste0("^", variable_pattern), normalized_names)
    }

    coefNames[matches]
}

#' Map model-matrix columns to their source formula terms
#'
#' @keywords internal
#' @noRd
removeRandomInterceptMethylationModels <- function(formulaText) {
    fixed_text <- sub("\\s*\\+\\s*\\(\\s*1\\s*\\|\\s*[^)]+\\)\\s*$",
    "", formulaText,
    perl = TRUE
    )
    if (identical(fixed_text, formulaText)) {
    stop("Could not isolate the fixed-effect formula from the mixed model.",
        call. = FALSE
    )
    }
    fixed_text
}

buildCoefficientTermMapMethylationModels <- function(
    formulaText,
    data, removeRandomEffects = FALSE
) {
    if (isTRUE(removeRandomEffects)) {
    formulaText <- removeRandomInterceptMethylationModels(formulaText)
    }
    formula_object <- stats::as.formula(formulaText)
    terms_object <- stats::terms(formula_object, data = data)
    design_terms <- stats::delete.response(terms_object)
    design <- stats::model.matrix(design_terms, data = data)
    assignments <- attr(design, "assign")
    term_labels <- attr(design_terms, "term.labels")
    mapped_terms <- rep("(Intercept)", length(assignments))
    non_intercept <- assignments > 0L
    mapped_terms[non_intercept] <- term_labels[assignments[non_intercept]]
    stats::setNames(mapped_terms, colnames(design))
}

validateOmnibusConfigurationMethylationGLM <- function(omnibusTest = FALSE) {
    omnibus_test <- validateLogicalScalarDnaEpico(
    omnibusTest,
    "omnibusTest"
    )
    if (isTRUE(omnibus_test) &&
    !requireNamespace("car", quietly = TRUE)) {
    stop(
        "omnibusTest = TRUE requires the package 'car'.",
        call. = FALSE
    )
    }
    omnibus_test
}

resolveOmnibusTargetTermMethylationGLM <- function(
    formulaText, data, phenotype, interactionTerm = NULL
) {
    terms_object <- stats::terms(
    stats::as.formula(formulaText),
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

emptyOmnibusResultMethylationGLM <- function(
    term, status = "not_estimable", reason = NA_character_
) {
    list(
    term = term, method = "car::linearHypothesis Wald F",
    fValue = NA_real_, numeratorDf = NA_real_,
    denominatorDf = NA_real_, pValue = NA_real_,
    status = status, reason = reason, modelMessage = NA_character_,
    rhs = 0
    )
}

computeOmnibusTestMethylationGLM <- function(fit,
    coefficientTerms, omnibusTerm) {
    result <- emptyOmnibusResultMethylationGLM(term = omnibusTerm)
    captured <- captureModelConditionsDnaEpico({
        fixed_effects <- stats::coef(fit)
        fixed_names <- names(fixed_effects)
        mapped_terms <- coefficientTerms[fixed_names]
        selected <- which(!is.na(mapped_terms) &
            mapped_terms == omnibusTerm)
        if (length(selected) == 0L) {
            stop("The omnibus model term has no estimable coefficients.",
                call. = FALSE)
        }
        contrast <- matrix(0, nrow = length(selected),
            ncol = length(fixed_effects), dimnames = list(fixed_names[selected],
                fixed_names))
        contrast[cbind(seq_along(selected), selected)] <- 1
        test <- car::linearHypothesis(model = fit,
            hypothesis.matrix = contrast, rhs = 0,
            test = "F", singular.ok = FALSE)
        required <- c("Res.Df", "Df", "F", "Pr(>F)")
        if (!all(required %in% colnames(test)) ||
            nrow(test) < 1L) {
            stop("The omnibus test did not return the expected F-test ",
                "statistics.", call. = FALSE)
        }
        tested_row <- nrow(test)
        values <- unlist(test[tested_row, required,
            drop = TRUE], use.names = FALSE)
        if (any(!is.finite(values))) {
            stop("The omnibus test returned missing or non-finite statistics.",
                call. = FALSE)
        }
        list(fValue = unname(test[["F"]][[tested_row]]),
            numeratorDf = unname(test[["Df"]][[tested_row]]),
            denominatorDf = unname(test[["Res.Df"]][[tested_row]]),
            pValue = unname(test[["Pr(>F)"]][[tested_row]]))
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

collectOmnibusTestsMethylationGLM <- function(fits, phenotype) {
    rows <- lapply(names(fits), function(cpg) {
    fit <- fits[[cpg]]
    if (is.null(fit) || inherits(
        fit,
        "dnaEPICO_methylationGLM_fit_error"
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

#' Validate a fixed-effect model matrix before fitting CpGs
#'
#' @keywords internal
#' @noRd
validateFixedEffectDesignMethylationModels <- function(
    formulaText,
    data, removeRandomEffects = FALSE
) {
    if (isTRUE(removeRandomEffects)) {
    formulaText <- removeRandomInterceptMethylationModels(formulaText)
    }
    formula_object <- stats::as.formula(formulaText)
    design_terms <- stats::delete.response(stats::terms(formula_object,
    data = data
    ))
    design <- stats::model.matrix(design_terms, data = data)

    if (nrow(design) == 0L) {
    stop("No complete observations remain for the requested model.",
        call. = FALSE
    )
    }
    if (any(!is.finite(design))) {
    stop(
        "The fixed-effect design contains non-finite numeric values. ",
        "Check phenotypes, covariates, and interactions.",
        call. = FALSE
    )
    }
    design_rank <- qr(design)$rank
    if (design_rank < ncol(design)) {
    stop(
        sprintf(
        "%s (rank %s of %s). %s",
        "The fixed-effect design matrix is rank deficient",
        design_rank, ncol(design),
        "Check duplicated covariates, factor levels, and interactions."
        ),
        call. = FALSE
    )
    }

    invisible(list(
    rows = nrow(design), columns = ncol(design),
    rank = design_rank
    ))
}

summarizeFitErrorsMethylationModels <- function(fitList, errorClass) {
    failed <- vapply(fitList, inherits, logical(1), what = errorClass)
    if (!any(failed)) {
    return(integer(0))
    }
    messages <- vapply(fitList[failed], function(x) x$error,
    character(1),
    USE.NAMES = FALSE
    )
    sort(table(messages), decreasing = TRUE)
}

formatFitErrorsMethylationModels <- function(errorCounts, maximum = 3L) {
    if (length(errorCounts) == 0L) {
    return("none")
    }
    keep <- seq_len(min(length(errorCounts), as.integer(maximum)))
    paste0(names(errorCounts)[keep], " (n=", as.integer(errorCounts[keep]),
    ")",
    collapse = "; "
    )
}

#' Build a GLM formula string for methylationGLM helpers
#'
#' @param phenotype Character. Phenotype variable of interest.
#' @param covariates Character vector of covariate variables.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return Character scalar containing a formula string.
#'
#' @description
#' Internal helper that builds the one-CpG Gaussian GLM formula.
#'
#' @keywords internal
#' @noRd
buildFormulaMethylationGLM <- function(
    phenotype, covariates = character(0),
    interactionTerm = NULL, responseVar = "beta"
) {
    if (!is.null(interactionTerm) && nzchar(interactionTerm) &&
    identical(interactionTerm, phenotype)) {
    stop("interactionTerm must differ from the phenotype being modelled.",
        call. = FALSE
    )
    }
    quoted_phenotype <- quoteNamesMethylationGLM(phenotype)
    quoted_covariates <- quoteNamesMethylationGLM(covariates)

    if (!is.null(interactionTerm) && nzchar(interactionTerm)) {
    quoted_interaction <- quoteNamesMethylationGLM(interactionTerm)
    interaction_part <- paste(quoted_phenotype, quoted_interaction,
        sep = " * "
    )
    fixed_terms <- setdiff(covariates, interactionTerm)
    quoted_fixed_terms <- quoteNamesMethylationGLM(fixed_terms)
    terms <- c(interaction_part, quoted_fixed_terms)
    } else {
    terms <- c(quoted_phenotype, quoted_covariates)
    }

    terms <- unique(terms[nzchar(terms)])
    if (length(terms) == 0L) {
    stop("At least one phenotype or covariate term is required.",
        call. = FALSE
    )
    }

    paste(responseVar, "~", paste(terms, collapse = " + "))
}

collectModelMessagesMethylationGLM <- function(fits) {
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
        stringsAsFactors = FALSE,
        check.names = FALSE
    ))
    }

    model_messages <- do.call(rbind, rows)
    rownames(model_messages) <- NULL
    model_messages
}

#' Summarize a single CpG-level GLM fit
#'
#' @param cpg Character. CpG identifier.
#' @param modelObj List returned by `fitCpGModelMethylationGLM()`.
#' @param variable Character. Phenotype variable of interest.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#' @param includeResidualSD Logical. If `TRUE`, append residual standard
#'   deviations to the output.
#'
#' @return A data frame for the requested CpG, or `NULL`.
#'
#' @description
#' Internal helper that extracts phenotype-specific coefficient rows from one
#' CpG-level GLM fit.
#'
#' @keywords internal
#' @noRd
summarizeCpGFitMethylationGLM <- function(
    cpg, modelObj, variable,
    interactionTerm = NULL, includeResidualSD = TRUE
) {
    if (is.null(modelObj) || inherits(
    modelObj,
    "dnaEPICO_methylationGLM_fit_error"
    )) {
    return(NULL)
    }

    coef_table <- modelObj$coef
    if (is.null(coef_table)) {
    return(NULL)
    }

    matched_rows <- findCoefficientRowsMethylationGLM(
    coefNames = rownames(coef_table),
    variable = variable, interactionTerm = interactionTerm,
    coefficientTerms = modelObj$coefficientTerms
    )
    if (length(matched_rows) == 0L) {
    return(NULL)
    }

    summary_df <- as.data.frame(coef_table[matched_rows, , drop = FALSE])
    summary_df$CpG <- cpg
    summary_df$Coefficient <- rownames(summary_df)
    summary_df$Model.Message <- modelMessageDnaEpico(modelObj)

    if (isTRUE(includeResidualSD) &&
    (!is.null(modelObj$residualSD) || !is.null(modelObj$residuals))) {
    summary_df$ResidualSD <- if (!is.null(modelObj$residualSD)) {
        modelObj$residualSD
    } else {
        finite_residuals <-
        modelObj$residuals[is.finite(modelObj$residuals)]
        if (length(finite_residuals) > 1L) {
        stats::sd(finite_residuals)
        } else {
        NA_real_
        }
    }
    }

    summary_df
}

resolveParallelBackendMethylationModels <- function(nCores) {
    n_cores <- validatePositiveIntegerMethylationModels(
    nCores,
    "nCores"
    )
    requested_backend <- tolower(Sys.getenv(
    "DNAEPICO_PARALLEL_BACKEND",
    "auto"
    ))
    if (!(requested_backend %in% c("auto", "fork", "psock", "serial"))) {
    requested_backend <- "auto"
    }

    if (n_cores <= 1L || identical(requested_backend, "serial")) {
    return("serial")
    }

    fork_supported <- !identical(.Platform$OS.type, "windows")
    if (identical(requested_backend, "fork")) {
    return(if (fork_supported) "fork" else "psock")
    }
    if (identical(requested_backend, "psock")) {
    return("psock")
    }

    if (fork_supported) {
    return("fork")
    }

    "psock"
}

parallelCrossoverMethylationModels <- function(engine) {
    engine <- tolower(as.character(engine[[1L]]))
    defaults <- c(glm2 = 25000L, lme4 = 1500L, nlme = 5000L)
    if (!(engine %in% names(defaults))) {
    stop("Unsupported methylation model engine: ", engine,
        call. = FALSE
    )
    }

    configured <- getOption("dnaEPICO.parallelCrossover")
    if (!is.null(configured) && !is.null(names(configured)) &&
    engine %in% names(configured)) {
    return(validatePositiveIntegerMethylationModels(
        configured[[engine]],
        paste0("dnaEPICO.parallelCrossover['", engine, "']")
    ))
    }

    env_name <- paste0("DNAEPICO_PARALLEL_CROSSOVER_", toupper(engine))
    env_value <- Sys.getenv(env_name, unset = "")
    if (nzchar(env_value)) {
    return(validatePositiveIntegerMethylationModels(
        env_value,
        env_name
    ))
    }

    unname(defaults[[engine]])
}

availableWorkersMethylationModels <- function() {
    detect_cores <- function(logical) {
    tryCatch(parallel::detectCores(logical = logical),
        warning = function(condition) NA_integer_,
        error = function(condition) NA_integer_
    )
    }
    physical_cores <- detect_cores(FALSE)
    logical_cores <- detect_cores(TRUE)
    available <- if (length(physical_cores) == 1L &&
    is.finite(physical_cores)) {
    physical_cores
    } else if (length(logical_cores) == 1L && is.finite(logical_cores)) {
    logical_cores
    } else {
    1L
    }

    configured_cap <- Sys.getenv("DNAEPICO_MAX_WORKERS", unset = "")
    if (nzchar(configured_cap)) {
    available <- min(available, validatePositiveIntegerMethylationModels(
        configured_cap,
        "DNAEPICO_MAX_WORKERS"
    ))
    }
    scheduler_values <- Sys.getenv(c(
    "SLURM_CPUS_PER_TASK", "NSLOTS", "PBS_NP"
    ), unset = "")
    scheduler_values <- vapply(
    scheduler_values,
    parseFiniteNumericMethylationModels,
    numeric(1)
    )
    scheduler_values <- scheduler_values[
    is.finite(scheduler_values) & scheduler_values >= 1
    ]
    if (length(scheduler_values) > 0L) {
    available <- min(available, scheduler_values)
    }

    max(1L, as.integer(available))
}

resolveParallelPlanMethylationModels <- function(engine, nCores, nCpGs,
    analysisData = NULL, modelData = NULL) {
    requested_cores <- validatePositiveIntegerMethylationModels(nCores,
        "nCores")
    n_cpgs <- max(0L, as.integer(nCpGs))
    requested_backend <- tolower(Sys.getenv("DNAEPICO_PARALLEL_BACKEND",
        "auto"))
    if (!(requested_backend %in% c("auto", "fork", "psock", "serial"))) {
        requested_backend <- "auto" }
    crossover <- parallelCrossoverMethylationModels(engine)
    resource_cap <- availableWorkersMethylationModels()
    worker_count <- min(requested_cores, resource_cap, max(1L, n_cpgs))
    memory_plan <- list(cap = worker_count, availableMB = NA_real_, reserveMB =
        NA_real_,
        estimatedWorkerMB = NA_real_, reason =
            "memory estimation was not requested")
    if (!is.null(analysisData) && !is.null(modelData)) {
        memory_plan <- memoryWorkerCapMethylationModels(engine = engine,
            analysisData = analysisData, modelData = modelData,
            requestedWorkers = worker_count)
        worker_count <- min(worker_count, memory_plan$cap)
    }
    if (requested_cores <= 1L || worker_count <= 1L) {
        backend <- "serial"
        reason <- "one effective worker" }
    else if (identical(requested_backend, "serial")) {
        backend <- "serial"; worker_count <- 1L
        reason <- "serial backend requested"
    }
    else if (identical(requested_backend, "auto") && n_cpgs < crossover) {
        backend <- "serial"
        worker_count <- 1L
        reason <- paste0("workload below the ", engine, " crossover of ",
            crossover, " CpGs") } else {
        backend <- resolveParallelBackendMethylationModels(worker_count)
        reason <- if (identical(requested_backend, "auto")) {
            paste0("workload reached the ", engine, " crossover")
        } else {
            paste0(requested_backend, " backend requested")
        } }
    list(engine = tolower(as.character(engine[[1L]])), backend = backend,
        requestedCores = requested_cores, workerCount = worker_count,
            resourceWorkerCap = resource_cap,
        workloadCpGs = n_cpgs, crossoverCpGs = crossover, reason = reason,
        forced = !identical(requested_backend, "auto"), availableMemoryMB =
            memory_plan$availableMB,
        reservedMemoryMB = memory_plan$reserveMB, estimatedWorkerMemoryMB =
            memory_plan$estimatedWorkerMB,
        memoryWorkerCap = memory_plan$cap, memoryReason = memory_plan$reason)
}

chunkCpGColumnsMethylationModels <- function(
    cpgColumns, nCores = 1L,
    batchesPerCore = 8L
) {
    if (length(cpgColumns) == 0L) {
    return(list())
    }

    n_cores <- validatePositiveIntegerMethylationModels(
    nCores,
    "nCores"
    )
    batches_per_core <- max(1L, as.integer(batchesPerCore))
    target_batches <- min(length(cpgColumns), max(1L, n_cores *
    batches_per_core))
    chunk_size <- max(1L, ceiling(length(cpgColumns) / target_batches))

    split(cpgColumns, ceiling(seq_along(cpgColumns) / chunk_size))
}

validateWorkerPackagesMethylationModels <- function(
    libPath = NULL,
    packages = character(0)
) {
    if (!is.null(libPath)) {
    .libPaths(unique(c(libPath, .libPaths())))
    }

    packages <- packages[!is.na(packages) & nzchar(packages)]
    for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Failed to load package: ", pkg, call. = FALSE)
    }
    }

    invisible(TRUE)
}

makePsockClusterMethylationModels <- function(clusterSize) {
    tryCatch(parallel::makeCluster(clusterSize,
    type = "PSOCK",
    useXDR = FALSE, methods = FALSE
    ), error = function(error) {
    parallel::makeCluster(clusterSize, type = "PSOCK")
    })
}

combineFitBatchResultsMethylationModels <- function(
    batchResults,
    cpgColumns
) {
    fit_chunks <- lapply(batchResults, function(result) result$fits)
    fit_chunks <- Filter(function(x) {
    !is.null(x) && length(x) >
        0L
    }, fit_chunks)
    fit_list <- if (length(fit_chunks) == 0L) {
    list()
    } else {
    out <- list()
    for (chunk in fit_chunks) {
        for (cpg in names(chunk)) {
        out[[cpg]] <- chunk[[cpg]]
        }
    }
    out
    }
    fit_list <- fit_list[cpgColumns[cpgColumns %in% names(fit_list)]]

    summary_chunks <- lapply(batchResults, function(result) result$summaries)
    summary_chunks <- Filter(function(x) {
    !is.null(x) && nrow(x) >
        0L
    }, summary_chunks)
    summary_df <- if (length(summary_chunks) == 0L) {
    data.frame()
    } else {
    out <- do.call(rbind, summary_chunks)
    rownames(out) <- NULL
    out
    }

    list(fits = fit_list, summaries = summary_df)
}

optionalTermMatchesMethylationModels <- function(requested, cached) {
    normalize <- function(value) {
    if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) {
        return("")
    }
    value <- as.character(value[[1L]])
    if (!nzchar(value)) {
        return("")
    }
    value
    }

    identical(normalize(requested), normalize(cached))
}

filterSummaryByPvalueMethylationGLM <- function(
    summaryDf, pValueFilter,
    includeResidualSD = TRUE
) {
    summary_df <- summaryDf
    if (is.null(summary_df) || nrow(summary_df) == 0L) {
    return(data.frame())
    }

    ordered_columns <- c(
    "CpG", "Coefficient", "Estimate", "Std. Error",
    "t value", "Pr(>|t|)", if (isTRUE(includeResidualSD)) "ResidualSD",
    "Model.Message"
    )
    ordered_columns <- intersect(ordered_columns, colnames(summary_df))
    summary_df <- summary_df[, ordered_columns, drop = FALSE]

    if (nrow(summary_df) > 0L && !is.na(pValueFilter)) {
    keep <- is.finite(summary_df[["Pr(>|t|)"]]) & summary_df[["Pr(>|t|)"]] <
        pValueFilter
    summary_df <- summary_df[keep, , drop = FALSE]
    }
    rownames(summary_df) <- NULL

    summary_df
}

fitCpGGlmValueDnaEpico <- function(
    cpgValues, modelData, formulaText, responseVar, retainModel,
    omnibusTest, omnibusTerm, formulaObject, coefficientTerms
) {
    if (!is.numeric(cpgValues)) {
    stop("The CpG response is not numeric.", call. = FALSE)
    }
    data <- modelData
    data[[responseVar]] <- cpgValues
    formula <- if (is.null(formulaObject)) {
    stats::as.formula(formulaText)
    } else {
    formulaObject
    }
    fit <- glm2::glm2(
    formula = formula, data = data, family = stats::gaussian(),
    na.action = stats::na.exclude
    )
    terms <- if (is.null(coefficientTerms)) {
    buildCoefficientTermMapMethylationModels(formulaText, data)
    } else {
    coefficientTerms
    }
    omnibus <- if (isTRUE(omnibusTest)) {
    computeOmnibusTestMethylationGLM(fit, terms, omnibusTerm)
    } else {
    NULL
    }
    residuals <- stats::residuals(fit)
    residual_sd <- if (fit$df.residual > 0L) {
    sqrt(sum(residuals^2, na.rm = TRUE) / fit$df.residual)
    } else {
    NA_real_
    }
    result <- list(
    coef = summary(fit)$coefficients, residualSD = residual_sd,
    coefficientTerms = terms, omnibus = omnibus,
    pValueAvailable = FALSE
    )
    if (isTRUE(retainModel)) {
    result$residuals <- residuals
    result$fitted <- stats::fitted(fit)
    result$model <- fit
    }
    result
}

completeCapturedGlmFitDnaEpico <- function(captured) {
    if (!is.null(captured$error)) {
    return(newMethylationFitErrorDnaEpico(
        reason = conditionMessage(captured$error),
        errorClass = "dnaEPICO_methylationGLM_fit_error",
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

fitCpGModelMethylationGLM <- function(
    cpg, cpgValues, modelData,
    formulaText, responseVar = "beta", retainModel = TRUE,
    omnibusTest = FALSE, omnibusTerm = NULL,
    formulaObject = NULL, coefficientTerms = NULL
) {
    captured <- captureModelConditionsDnaEpico(
    fitCpGGlmValueDnaEpico(
        cpgValues, modelData, formulaText, responseVar, retainModel,
        omnibusTest, omnibusTerm, formulaObject, coefficientTerms
    )
    )
    completeCapturedGlmFitDnaEpico(captured)
}

fitMethylationGLMBatch <- function(cpgBatch, data,
    modelData, formulaText, phenotype, interactionTerm = NULL,
    responseVar = "beta", omnibusTest = FALSE,
    omnibusTerm = NULL, formulaObject = NULL, coefficientTerms = NULL) {
    fits <- vector("list", length(cpgBatch))
    names(fits) <- cpgBatch
    summaries <- vector("list", length(cpgBatch))
    names(summaries) <- cpgBatch
    for (cpg in cpgBatch) {
        model_obj <- fitCpGModelMethylationGLM(cpg = cpg,
            cpgValues = cpgResponseMethylationModels(data,
                cpg), modelData = modelData, formulaText = formulaText,
            responseVar = responseVar, retainModel = FALSE,
            omnibusTest = omnibusTest, omnibusTerm = omnibusTerm,
            formulaObject = formulaObject, coefficientTerms = coefficientTerms)
        summary_row <- summarizeCpGFitMethylationGLM(cpg = cpg,
            modelObj = model_obj, variable = phenotype,
            interactionTerm = interactionTerm,
            includeResidualSD = TRUE)
        p_value_available <- !is.null(summary_row) &&
            "Pr(>|t|)" %in% colnames(summary_row) &&
            any(is.finite(summary_row[["Pr(>|t|)"]]))
        omnibus_p_available <- !is.null(model_obj$omnibus) &&
            is.finite(model_obj$omnibus$pValue)
        p_value_available <- p_value_available ||
            omnibus_p_available
        model_obj$pValueAvailable <- p_value_available
        if (!is.null(summary_row)) {
            summary_row$Model.Message <- modelMessageDnaEpico(model_obj)
        }
        fits[[cpg]] <- model_obj
        summaries[[cpg]] <- summary_row
    }
    summaries <- Filter(Negate(is.null), summaries)
    summary_df <- if (length(summaries) == 0L) {
        data.frame() } else {
        out <- do.call(rbind, summaries)
        rownames(out) <- NULL
        out }
    list(coefficientResults = compactCoefficientResultsMethylationModels(fits =
        fits,
        cpgOrder = cpgBatch, includeResidualSD = TRUE),
        summaries = summary_df, omnibusTests =
            collectOmnibusTestsMethylationGLM(fits = fits,
            phenotype = phenotype), modelMessages =
            collectBatchModelMessagesMethylationModels(fits,
            phenotype), fitFailures = collectBatchFitFailuresMethylationModels(
            fits,
            phenotype, "dnaEPICO_methylationGLM_fit_error"))
}

#' Create a distribution plot used by methylationGLM helpers
#'
#' @param values Vector of phenotype or covariate values to plot.
#' @param variable Character. Variable label.
#' @param type Character. Either `'hist'` or `'bar'`.
#' @param fill Character. Fill colour used for the bars.
#'
#' @return A `ggplot2` object.
#'
#' @description
#' Internal helper that builds exploratory distribution plots.
#'
#' @keywords internal
#' @noRd
createDistributionPlotMethylationGLM <- function(values,
    variable, type = c("hist", "bar"), fill = "steelblue") {
    type <- match.arg(type)
    if (identical(type, "hist")) {
        numeric_values <- coerceNumericDnaEpico(values)
        observed <- numeric_values[is.finite(numeric_values)]
        bins <- max(5L, min(30L, ceiling(sqrt(max(1L, length(observed))))))
        plot_data <- data.frame(value = observed)
        if (length(observed) < 30L) {
            plot_data$row <- 1
            plot_object <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value,
                y = row)) + ggplot2::geom_jitter(height = 0.08,
                width = 0, colour = fill, size = 2.4, alpha = 0.75) +
                ggplot2::geom_boxplot(ggplot2::aes(group = row),
                    width = 0.18, outlier.shape = NA, fill = NA,
                    colour = "#17324D") + ggplot2::scale_y_continuous(NULL,
                breaks = NULL) + ggplot2::labs(title = NULL,
                x = variable) + dnaEpicoModelPlotTheme()
        } else {
            plot_object <- ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
                ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(
            density)),
                    bins = bins, fill = fill, color = "white",
                    alpha = 0.82) + ggplot2::labs(title = NULL,
                x = variable, y = "Density") + dnaEpicoModelPlotTheme()
            if (length(observed) <= 2000L) {
                plot_object <- plot_object + ggplot2::geom_rug(alpha = 0.25)
            }
            if (length(unique(observed)) > 1L) {
                plot_object <- plot_object + ggplot2::geom_density(colour =
            "#111827",
                    linewidth = 0.8, adjust = 1)
            } }; return(plot_object) }
    character_values <- as.character(values)
    character_values[is.na(character_values) | !nzchar(trimws(
        character_values))] <- "Missing"
    counts <- as.data.frame(table(character_values), stringsAsFactors = FALSE)
    names(counts) <- c("value", "count")
    counts$percentage <- 100 * counts$count/sum(counts$count)
    counts$label <- sprintf("%s (%.1f%%)", counts$count,
        counts$percentage)
    counts$value <- stats::reorder(counts$value, counts$count)
    ggplot2::ggplot(counts, ggplot2::aes(x = value, y = count)) +
        ggplot2::geom_col(fill = fill, alpha = 0.85) + ggplot2::geom_text(
            ggplot2::aes(label = label),
        hjust = -0.08, size = 3.5) + ggplot2::coord_flip(clip = "off") +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0,
            0.18))) + ggplot2::labs(title = NULL, x = variable,
        y = "Count") + dnaEpicoModelPlotTheme()
}

#' Resolve an annotation object for methylationGLM helpers
#'
#' @param annotationObject Character package/object name, annotation data frame,
#'   or annotation object understood by `minfi::getAnnotation()`.
#'
#' @return An object accepted by `minfi::getAnnotation()`.
#'
#' @description
#' Internal helper that loads or validates the annotation resource used to
#' annotate CpG-level summary tables.
#'
#' @keywords internal
#' @noRd
resolveAnnotationObjectMethylationGLM <- function(annotationObject) {
    if (!is.character(annotationObject) || length(annotationObject) !=
    1L) {
    return(annotationObject)
    }

    if (requireNamespace(annotationObject, quietly = TRUE)) {
    annotation_namespace <- asNamespace(annotationObject)
    if (exists(annotationObject,
        envir = annotation_namespace,
        inherits = FALSE
    )) {
        return(get(annotationObject,
        envir = annotation_namespace,
        inherits = FALSE
        ))
    }
    }

    loaded_namespaces <- loadedNamespaces()
    for (namespace_name in loaded_namespaces) {
    namespace_env <- asNamespace(namespace_name)
    if (exists(annotationObject, envir = namespace_env, inherits = FALSE)) {
        return(get(annotationObject,
        envir = namespace_env,
        inherits = FALSE
        ))
    }
    }

    if (exists(annotationObject, inherits = TRUE)) {
    return(get(annotationObject, inherits = TRUE))
    }

    stop("Annotation package or object was not found: ", annotationObject,
    call. = FALSE
    )
}

annotationPackageDataMethylationGLM <- function(annotationPackage) {
    if (!requireNamespace(annotationPackage,
        quietly = TRUE)) {
        return(NULL)
    }
    data_index <- utils::data(package = annotationPackage)$results
    if (is.null(data_index) || !(annotationPackage %in%
        data_index[, "Item"])) {
        return(NULL)
    }
    annotation_environment <- new.env(parent = emptyenv())
    utils::data(list = unique(data_index[,
        "Item"]), package = annotationPackage,
        envir = annotation_environment)
    if (!exists(annotationPackage, envir = annotation_environment,
        inherits = FALSE)) {
        return(NULL)
    }
    annotation_resource <- get(annotationPackage,
        envir = annotation_environment, inherits = FALSE)
    annotation_data <- methods::slot(annotation_resource,
        "data")
    annotation_names <- ls(envir = annotation_data,
        all.names = TRUE)
    for (annotation_name in annotation_names) {
        if (exists(annotation_name, envir = annotation_environment,
            inherits = FALSE)) {
            assign(annotation_name, get(annotation_name,
                envir = annotation_environment,
                inherits = FALSE), envir = annotation_data)
        }
    }
    annotation_df <- as.data.frame(minfi::getAnnotation(annotation_resource))
    annotation_df$CpG <- rownames(annotation_df)
    annotation_df
}

coerceAnnotationDataMethylationGLM <- function(annotationObject) {
    if (is.data.frame(annotationObject)) {
    annotation_df <- annotationObject
    if (!("CpG" %in% colnames(annotation_df))) {
        if ("IlmnID" %in% colnames(annotation_df)) {
        annotation_df$CpG <- annotation_df$IlmnID
        } else if (.row_names_info(annotation_df, type = 1L) <
        0L) {
        stop(
            "annotationObject data frames must include a CpG or ",
            "IlmnID column, or explicit probe row names.",
            call. = FALSE
        )
        } else {
        annotation_df$CpG <- rownames(annotation_df)
        }
    }
    } else if (is.character(annotationObject) &&
    length(annotationObject) == 1L) {
    annotation_df <- annotationPackageDataMethylationGLM(annotationObject)
    if (is.null(annotation_df)) {
        annotation_source <-
        resolveAnnotationObjectMethylationGLM(annotationObject)
        annotation_df <- as.data.frame(
        minfi::getAnnotation(annotation_source)
        )
        annotation_df$CpG <- rownames(annotation_df)
    }
    } else {
    annotation_source <-
        resolveAnnotationObjectMethylationGLM(annotationObject)
    annotation_df <- as.data.frame(minfi::getAnnotation(annotation_source))
    annotation_df$CpG <- rownames(annotation_df)
    }

    annotation_df$CpG <- validateMethylationProbeIdentifiersDnaEpico(
    annotation_df$CpG,
    "Annotation probe identifiers"
    )
    annotation_df
}

selectWorkbookFormulaDnaEpico <- function(column, formulaValues,
                                            defaultFormula) {
    if (length(formulaValues) == 0L) {
    return(defaultFormula)
    }
    formula_names <- names(formulaValues)
    if (!is.null(formula_names) && any(nzchar(formula_names))) {
    for (formula_name in formula_names[nzchar(formula_names)]) {
        prefix <- paste0("^", escapeRegexMethylationGLM(formula_name))
        if (grepl(prefix, column)) {
        return(unname(formulaValues[[formula_name]]))
        }
    }
    }
    if (length(formulaValues) == 1L) {
    return(unname(formulaValues[[1L]]))
    }
    paste(unique(unname(formulaValues)), collapse = "; ")
}

displayWorkbookFormulaDnaEpico <- function(column, formulaValues,
                                            defaultFormula, modelLabel,
                                            responseLabel) {
    formula_text <- selectWorkbookFormulaDnaEpico(
    column, formulaValues, defaultFormula
    )
    formula_parts <- trimws(unlist(strsplit(
    formula_text, "\\s*;\\s*",
    perl = TRUE
    ), use.names = FALSE))
    formula_parts <- vapply(formula_parts, function(formula_part) {
    sub(
        "^\\s*[^~]+\\s*~", paste(responseLabel, "~"),
        formula_part,
        ignore.case = TRUE
    )
    }, character(1))
    paste0(modelLabel, ": ", paste(formula_parts, collapse = "; "))
}

workbookColumnDescriptionDnaEpico <- function(column, isPvalue,
                                                modelDescription) {
    descriptions <- c(
    "_Omnibus_F\\.Value$" =
        "F statistic for the joint fixed-effect omnibus test",
    "_Omnibus_Num\\.DF$" = paste0(
        "Numerator degrees of freedom equal to the estimable ",
        "rank of the omnibus contrast"
    ),
    "_Omnibus_Den\\.DF$" =
        "Denominator degrees of freedom for the omnibus F test",
    "_Omnibus_Adjusted\\.P\\.Value$" = paste0(
        "Omnibus p-value adjusted across valid CpGs within ",
        "the phenotype and tested term"
    ),
    "_Omnibus_P\\.Value$" = paste0(
        "Raw p-value for the joint null hypothesis that all ",
        "estimable coefficients in the tested fixed-effect ",
        "term equal zero"
    ),
    "_Omnibus_Method$" = "Method used for the omnibus F test"
    )
    matches <- vapply(names(descriptions), grepl, logical(1), x = column)
    if (any(matches)) {
    return(unname(descriptions[which(matches)[[1L]]]))
    }
    if (isPvalue) {
    return(modelDescription)
    }
    if (column %in% c("IlmnID", "CpG", "Name")) {
    return("CpG probe identifier")
    }
    if (grepl("_Model\\.Message$", column)) {
    return(paste(
        "Messages, warnings, or errors reported while fitting",
        "or summarizing the phenotype-specific CpG model"
    ))
    }
    "Genomic annotation or supporting result column"
}

buildAnnotatedWorkbookDictionaryMethylationGLM <- function(
    columns, modelDescription, formulaText, modelLabel, responseLabel
) {
    pvalue_columns <- grepl("P\\.Value$|P\\.value$", columns)
    omnibus_columns <- grepl("_Omnibus_", columns, fixed = TRUE)
    default_formula <- paste(responseLabel, "~ formula unavailable")
    formula_values <- stats::setNames(
    as.character(formulaText), names(formulaText)
    )
    formula_values <- formula_values[
    !is.na(formula_values) & nzchar(formula_values)
    ]
    descriptions <- vapply(seq_along(columns), function(index) {
    workbookColumnDescriptionDnaEpico(
        columns[[index]], pvalue_columns[[index]], modelDescription
    )
    }, character(1))
    formulas <- rep("", length(columns))
    formula_columns <- pvalue_columns | omnibus_columns
    formulas[formula_columns] <- vapply(
    columns[formula_columns], function(column) {
        displayWorkbookFormulaDnaEpico(
        column, formula_values, default_formula,
        modelLabel, responseLabel
        )
    }, character(1)
    )
    data.frame(
    Column = columns, Description = descriptions,
    Formula = formulas, stringsAsFactors = FALSE, check.names = FALSE
    )
}

orderOmnibusColumnsDnaEpico <- function(values, allColumns) {
    suffix_order <- c(
    "_Omnibus_F.Value", "_Omnibus_Num.DF", "_Omnibus_Den.DF",
    "_Omnibus_P.Value", "_Omnibus_Adjusted.P.Value",
    "_Omnibus_Method"
    )
    ranks <- vapply(values, function(value) {
    matches <- which(endsWith(value, suffix_order))
    if (length(matches) == 0L) {
        length(suffix_order) + 1L
    } else {
        matches[[1L]]
    }
    }, integer(1))
    values[order(ranks, match(values, allColumns))]
}

modelColumnOwnersDnaEpico <- function(columns, modelNames) {
    vapply(columns, function(column) {
    candidates <- modelNames[startsWith(column, paste0(modelNames, "_"))]
    if (length(candidates) == 0L) {
        return(NA_character_)
    }
    candidates[[which.max(nchar(candidates))]]
    }, character(1))
}

groupModelColumnsDnaEpico <- function(columns, owners, modelNames,
                                        primary = NULL, orderOmnibus = FALSE,
                                        allColumns = columns) {
    grouped <- unlist(lapply(modelNames, function(model_name) {
    model_columns <- columns[owners == model_name]
    if (orderOmnibus) {
        return(orderOmnibusColumnsDnaEpico(model_columns, allColumns))
    }
    c(intersect(primary, model_columns), setdiff(model_columns, primary))
    }), use.names = FALSE)
    unique(c(grouped, columns[is.na(owners)]))
}

orderAnnotatedModelColumnsDnaEpico <- function(data,
    annotationCols = character(0), modelNames = NULL) {
    columns <- colnames(data)
    id_columns <- intersect(c("IlmnID", "CpG"),
        columns)
    omnibus_columns <- columns[grepl("_Omnibus_",
        columns, fixed = TRUE)]
    omnibus_pvalue_columns <- omnibus_columns[endsWith(omnibus_columns,
        "_Omnibus_P.Value")]
    omnibus_detail_columns <- setdiff(omnibus_columns,
        omnibus_pvalue_columns)
    pvalue_columns <- setdiff(columns[grepl("P\\.Value$|P\\.value$",
        columns)], omnibus_columns)
    annotation_columns <- intersect(annotationCols,
        columns)
    diagnostic_columns <- columns[grepl("_Model\\.Message$",
        columns)]
    supporting_columns <- setdiff(columns,
        c(id_columns, pvalue_columns, omnibus_columns,
            annotation_columns, diagnostic_columns))
    result_columns <- c(pvalue_columns, omnibus_pvalue_columns)
    model_names <- as.character(modelNames)
    model_names <- model_names[!is.na(model_names) &
        nzchar(model_names)]
    if (length(model_names) > 0L && length(result_columns) >
        0L) {
        owners <- modelColumnOwnersDnaEpico(result_columns,
            model_names)
        result_columns <- groupModelColumnsDnaEpico(result_columns,
            owners, model_names, primary = pvalue_columns)
    }
    ordered_omnibus_details <- orderOmnibusColumnsDnaEpico(
        omnibus_detail_columns,
        columns)
    if (length(model_names) > 0L && length(ordered_omnibus_details) >
        0L) {
        owners <- modelColumnOwnersDnaEpico(ordered_omnibus_details,
            model_names)
        ordered_omnibus_details <- groupModelColumnsDnaEpico(
            ordered_omnibus_details,
            owners, model_names, orderOmnibus = TRUE,
            allColumns = columns)
    }
    ordered_columns <- unique(c(id_columns,
        result_columns, annotation_columns,
        ordered_omnibus_details, supporting_columns,
        diagnostic_columns))
    data[, ordered_columns, drop = FALSE]
}

modelSettingTextDnaEpico <- function(value, empty = "None") {
    if (is.null(value) || length(value) == 0L || all(is.na(value))) {
    return(empty)
    }
    value <- as.character(value)
    value <- value[!is.na(value) & nzchar(value)]
    if (length(value) == 0L) {
    empty
    } else {
    paste(value, collapse = ",")
    }
}

modelMessageMetricsDnaEpico <- function(modelResults) {
    model_messages <- modelResults$modelMessages
    if (!is.data.frame(model_messages)) {
    model_messages <- data.frame()
    }
    message_values <- if ("Model.Message" %in% names(model_messages)) {
    as.character(model_messages$Model.Message)
    } else {
    character(0)
    }
    p_value_available <- if ("P.Value.Available" %in%
    names(model_messages)) {
    as.logical(model_messages$P.Value.Available)
    } else {
    logical(0)
    }
    attempted_cpgs <- if ("CpG" %in% names(model_messages)) {
    unique(as.character(model_messages$CpG))
    } else {
    character(0)
    }
    available_by_cpg <- if (length(attempted_cpgs) > 0L) {
    vapply(attempted_cpgs, function(cpg) {
        any(p_value_available[model_messages$CpG == cpg], na.rm = TRUE)
    }, logical(1))
    } else {
    logical(0)
    }
    has_message <- !is.na(message_values) & nzchar(message_values)
    list(
    messages = model_messages, values = message_values,
    pValueAvailable = p_value_available,
    attemptedCpGs = attempted_cpgs,
    availableByCpG = available_by_cpg,
    hasMessage = has_message
    )
}

modelEngineMetadataDnaEpico <- function(analysis, settings) {
    lme_engine <- if (identical(analysis, "lme")) {
    modelSettingTextDnaEpico(settings$lmeEngine, empty = "lme4")
    } else {
    "glm2"
    }
    engine_package <- if (identical(analysis, "glm")) {
    "glm2"
    } else if (identical(lme_engine, "nlme")) {
    "nlme"
    } else {
    "lmerTest"
    }
    engine_version <-
    tryCatch(as.character(utils::packageVersion(engine_package)),
        error = function(error) "unavailable"
    )
    fitting_function <- if (identical(analysis, "glm")) {
    "glm2::glm2"
    } else if (identical(lme_engine, "nlme")) {
    "nlme::lme"
    } else {
    "lmerTest::lmer"
    }
    list(
    backend = lme_engine, fittingFunction = fitting_function,
    engineVersion = engine_version
    )
}

annotationMetadataDnaEpico <- function(annotatedResults) {
    annotation_used <- if (!is.null(annotatedResults$annotationColumnsUsed)) {
    annotatedResults$annotationColumnsUsed
    } else {
    character(0)
    }
    annotation_missing <-
    if (!is.null(annotatedResults$missingAnnotationCols)) {
        annotatedResults$missingAnnotationCols
    } else {
        character(0)
    }
    list(used = annotation_used, missing = annotation_missing)
}

baseWorkbookMetadataDnaEpico <- function(modelResults,
    analysis, metrics, engine, annotation) {
    settings <- modelResults$settings
    keys <- c("analysis", "backend", "fitting_function",
        "engine_version", "formulas", "response_label",
        "methylation_scale", "sample_count",
        "phenotypes", "covariates", "factor_variables",
        "scaled_variables", "interaction_term",
        "cpg_count", "model_attempt_count",
        "model_with_p_value_count", "model_without_p_value_count",
        "cpg_with_p_value_count", "cpg_without_p_value_count",
        "model_message_count", "model_warning_count",
        "model_error_count", "parallel_backend",
        "worker_count", "annotation_columns",
        "missing_annotation_columns", "created")
    values <- c(if (identical(analysis,
        "glm")) "methylationGLM" else "methylationLME",
        engine$backend, engine$fittingFunction,
        engine$engineVersion, paste(paste(names(modelResults$formulas),
            modelResults$formulas, sep = ": "),
            collapse = " | "), inferMethylationValueLabelMethylationGLM(
            modelResults),
        modelSettingTextDnaEpico(settings$methylationScale),
        modelSettingTextDnaEpico(settings$sampleCount),
        modelSettingTextDnaEpico(modelResults$phenotypes),
        modelSettingTextDnaEpico(settings$covariates),
        modelSettingTextDnaEpico(settings$factorVars),
        modelSettingTextDnaEpico(settings$scaleVars),
        modelSettingTextDnaEpico(settings$interactionTerm),
        as.character(length(metrics$attemptedCpGs)),
        as.character(nrow(metrics$messages)),
        as.character(sum(metrics$pValueAvailable,
            na.rm = TRUE)), as.character(sum(!metrics$pValueAvailable,
            na.rm = TRUE)), as.character(sum(metrics$availableByCpG)),
        as.character(sum(!metrics$availableByCpG)),
        as.character(sum(metrics$hasMessage)),
        as.character(sum(grepl("WARNING:",
            metrics$values, fixed = TRUE),
            na.rm = TRUE)), as.character(sum(grepl("ERROR:",
            metrics$values, fixed = TRUE),
            na.rm = TRUE)), modelSettingTextDnaEpico(settings$parallelBackend),
        modelSettingTextDnaEpico(settings$workerCount),
        modelSettingTextDnaEpico(annotation$used),
        modelSettingTextDnaEpico(annotation$missing),
        format(Sys.time(), tz = "UTC", usetz = TRUE))
    data.frame(Key = keys, Value = values,
        stringsAsFactors = FALSE, check.names = FALSE)
}

factorLevelMetadataDnaEpico <- function(settings) {
    factor_levels <- settings$factorLevels
    if (!is.list(factor_levels) || length(factor_levels) == 0L) {
    return(NULL)
    }
    data.frame(
    Key = paste0("factor.", names(factor_levels), ".levels"),
    Value = vapply(
        factor_levels, modelSettingTextDnaEpico, character(1)
    ),
    stringsAsFactors = FALSE, check.names = FALSE
    )
}

omnibusMetadataCountsDnaEpico <- function(modelSummaries) {
    omnibus_tables <- modelSummaries$omnibusTests
    omnibus_rows <- if (is.list(omnibus_tables)) {
    Filter(function(table) {
        is.data.frame(table) && nrow(table) > 0L
    }, omnibus_tables)
    } else {
    list()
    }
    omnibus_data <- if (length(omnibus_rows) > 0L) {
    do.call(rbind, omnibus_rows)
    } else {
    data.frame()
    }
    if (!"Omnibus.Status" %in% colnames(omnibus_data)) {
    return(c(tested = 0L, unavailable = 0L))
    }
    c(
    tested = sum(omnibus_data$Omnibus.Status == "tested", na.rm = TRUE),
    unavailable = sum(
        omnibus_data$Omnibus.Status != "tested",
        na.rm = TRUE
    )
    )
}

modelAnalysisMetadataDnaEpico <- function(modelResults,
    modelSummaries, analysis) {
    settings <- modelResults$settings
    counts <- omnibusMetadataCountsDnaEpico(modelSummaries)
    scope <- "within each phenotype and omnibus term across valid CpGs"
    if (identical(analysis, "lme")) {
        keys <- c("libraries", "person_variable",
            "time_variable", "correlation_structure",
            "correlation_variable", "random_effect_structure",
            "omnibus_test", "omnibus_target", "omnibus_ddf",
            "omnibus_rhs", "omnibus_joint", "omnibus_tested_count",
            "omnibus_unavailable_count", "p_adjust_method",
            "p_adjust_scope", "lmerTest_version",
            "pbkrtest_version")
        values <- c(modelSettingTextDnaEpico(settings$lmeLibs),
            modelSettingTextDnaEpico(settings$personVar),
            modelSettingTextDnaEpico(settings$timeVar),
            modelSettingTextDnaEpico(settings$correlationStructure,
                empty = "none"), modelSettingTextDnaEpico(
            settings$correlationVar),
            paste0("(1 | ", modelSettingTextDnaEpico(settings$personVar),
                ")"), modelSettingTextDnaEpico(settings$omnibusTest),
            modelSettingTextDnaEpico(modelResults$omnibusTargets),
            modelSettingTextDnaEpico(settings$omnibusDdf),
            modelSettingTextDnaEpico(settings$omnibusRhs),
            modelSettingTextDnaEpico(settings$omnibusJoint),
            as.character(counts[["tested"]]), as.character(counts[[
            "unavailable"]]),
            modelSettingTextDnaEpico(modelSummaries$settings$padjmethod),
            scope, packageVersionTextDnaEpico("lmerTest"),
            packageVersionTextDnaEpico("pbkrtest"))
    } else {
        keys <- c("libraries", "omnibus_test",
            "omnibus_target", "omnibus_method",
            "omnibus_rhs", "omnibus_joint", "omnibus_tested_count",
            "omnibus_unavailable_count", "p_adjust_method",
            "p_adjust_scope", "car_version")
        values <- c(modelSettingTextDnaEpico(settings$glmLibs),
            modelSettingTextDnaEpico(settings$omnibusTest),
            modelSettingTextDnaEpico(modelResults$omnibusTargets),
            modelSettingTextDnaEpico(settings$omnibusMethod),
            modelSettingTextDnaEpico(settings$omnibusRhs),
            modelSettingTextDnaEpico(settings$omnibusJoint),
            as.character(counts[["tested"]]), as.character(counts[[
            "unavailable"]]),
            modelSettingTextDnaEpico(modelSummaries$settings$padjmethod),
            scope, packageVersionTextDnaEpico("car"))
    }
    data.frame(Key = keys, Value = values, stringsAsFactors = FALSE,
        check.names = FALSE) }

packageVersionTextDnaEpico <- function(package) {
    tryCatch(
    as.character(utils::packageVersion(package)),
    error = function(error) "unavailable"
    )
}

scalingMetadataRowsDnaEpico <- function(settings) {
    scaling <- settings$scalingMetadata
    if (!is.data.frame(scaling) || nrow(scaling) == 0L) {
    return(NULL)
    }
    do.call(rbind, lapply(seq_len(nrow(scaling)), function(index) {
    variable <- scaling$Variable[[index]]
    data.frame(
        Key = c(
        paste0("scale.", variable, ".center"),
        paste0("scale.", variable, ".sd"),
        paste0("scale.", variable, ".finite_values"),
        paste0("scale.", variable, ".missing_values")
        ),
        Value = as.character(c(
        scaling$Center[[index]], scaling$Scale[[index]],
        scaling$Finite.Values[[index]],
        scaling$Missing.Values[[index]]
        )),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    }))
}

buildModelWorkbookMetadataDnaEpico <- function(
    modelResults, modelSummaries, annotatedResults,
    analysis = c("glm", "lme")
) {
    analysis <- match.arg(analysis)
    metrics <- modelMessageMetricsDnaEpico(modelResults)
    engine <- modelEngineMetadataDnaEpico(analysis, modelResults$settings)
    annotation <- annotationMetadataDnaEpico(annotatedResults)
    metadata <- baseWorkbookMetadataDnaEpico(
    modelResults, analysis, metrics, engine, annotation
    )
    metadata <- rbind(
    metadata,
    factorLevelMetadataDnaEpico(modelResults$settings),
    modelAnalysisMetadataDnaEpico(modelResults, modelSummaries, analysis),
    scalingMetadataRowsDnaEpico(modelResults$settings)
    )
    rownames(metadata) <- NULL
    metadata
}

writeAnnotatedWorkbookMethylationGLM <- function(
    annotated_df,
    file, resultSheet, dictionary, metadata = NULL, extraSheets = list()
) {
    if (!is.list(extraSheets) || (length(extraSheets) &&
    is.null(names(extraSheets)))) {
    stop("extraSheets must be a named list of data frames.",
        call. = FALSE
    )
    }
    if (length(extraSheets)) {
    invalid <- !vapply(extraSheets, is.data.frame, logical(1))
    if (any(invalid) || any(!nzchar(names(extraSheets))) ||
        anyDuplicated(tolower(names(extraSheets)))) {
        stop("extraSheets must contain uniquely named data frames.",
        call. = FALSE
        )
    }
    reserved <- c(resultSheet, "metadata", "dictionary")
    if (any(tolower(names(extraSheets)) %in% tolower(reserved))) {
        stop("extraSheets cannot use a reserved workbook sheet name.",
        call. = FALSE
        )
    }
    }
    workbook <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(workbook, resultSheet)
    openxlsx::writeData(workbook, sheet = resultSheet, x = annotated_df)
    if (is.data.frame(metadata) && nrow(metadata) > 0L) {
    openxlsx::addWorksheet(workbook, "metadata")
    openxlsx::writeData(workbook, sheet = "metadata", x = metadata)
    }
    for (sheet in names(extraSheets)) {
    openxlsx::addWorksheet(workbook, sheet)
    openxlsx::writeData(workbook, sheet = sheet, x = extraSheets[[sheet]])
    }
    openxlsx::addWorksheet(workbook, "dictionary")
    openxlsx::writeData(workbook, sheet = "dictionary", x = dictionary)
    openxlsx::saveWorkbook(workbook, file = file, overwrite = TRUE)

    invisible(file)
}

sortReportTableDnaEpico <- function(tableData) {
    table_data <- as.data.frame(tableData,
    stringsAsFactors = FALSE,
    check.names = FALSE
    )
    id_column <- if ("IlmnID" %in% names(table_data)) {
    "IlmnID"
    } else if ("CpG" %in% names(table_data)) {
    "CpG"
    } else if (ncol(table_data) > 0L) {
    names(table_data)[[1L]]
    } else {
    ""
    }

    if (nzchar(id_column) && nrow(table_data) > 1L) {
    row_order <- order(as.character(table_data[[id_column]]),
        na.last = TRUE
    )
    if (!identical(row_order, seq_len(nrow(table_data)))) {
        table_data <- table_data[row_order, , drop = FALSE]
        rownames(table_data) <- NULL
    }
    }

    list(data = table_data, idColumn = id_column)
}

reportTableSidecarPathsDnaEpico <- function(sidecarDir, sheet) {
    list(
    table = file.path(sidecarDir, paste0(sheet, ".tsv.gz")),
    metadata = file.path(sidecarDir, paste0(sheet, ".report.tsv")),
    dictionary = file.path(sidecarDir, paste0(sheet, ".dictionary.tsv")),
    workbookMetadata = file.path(sidecarDir, paste0(
        sheet,
        ".metadata.tsv"
    ))
    )
}

writeReportTableSidecarDnaEpico <- function(tableData,
    workbookFile, sidecarDir, sheet, idColumn,
    dictionary, workbookMetadata = NULL) {
    paths <- reportTableSidecarPathsDnaEpico(sidecarDir,
        sheet)
    dir.create(dirname(paths$table), recursive = TRUE,
        showWarnings = FALSE)
    if (!(idColumn %in% names(tableData))) {
        stop("The report table identifier column was not found.",
            call. = FALSE)
    }
    row_order <- order(as.character(tableData[[idColumn]]),
        na.last = TRUE)
    if (!identical(row_order, seq_len(nrow(tableData)))) {
        stop("The report table must be sorted by its identifier column.",
            call. = FALSE)
    }
    data.table::fwrite(tableData, file = paths$table,
        sep = "\t", quote = TRUE, na = "",
        compress = "gzip", showProgress = FALSE)
    sidecar_metadata <- data.frame(Key = c("format_version",
        "sheet", "workbook_md5", "rows",
        "columns", "id_column", "sorted_by_id"),
        Value = c("1", sheet, unname(tools::md5sum(workbookFile)),
            as.character(nrow(tableData)),
            as.character(ncol(tableData)),
            idColumn, "TRUE"), stringsAsFactors = FALSE,
        check.names = FALSE)
    data.table::fwrite(dictionary, file = paths$dictionary,
        sep = "\t", quote = TRUE, na = "")
    if (is.data.frame(workbookMetadata) &&
        nrow(workbookMetadata) > 0L) {
        data.table::fwrite(workbookMetadata,
            file = paths$workbookMetadata,
            sep = "\t", quote = TRUE, na = "")
    }
    else if (file.exists(paths$workbookMetadata)) {
        unlink(paths$workbookMetadata, force = TRUE)
    }
    data.table::fwrite(sidecar_metadata,
        file = paths$metadata, sep = "\t",
        quote = TRUE, na = "")
    paths
}

resolveReportTableSidecarDnaEpico <- function(workbookFile,
    sidecarDir, sheet) {
    paths <- reportTableSidecarPathsDnaEpico(sidecarDir,
        sheet)
    unavailable <- function(reason) {
        list(ok = FALSE, reason = reason,
            table = paths$table, metadata = paths$metadata)
    }
    if (!file.exists(paths$table) || !file.exists(paths$metadata) ||
        !file.exists(paths$dictionary)) {
        return(unavailable("A report table sidecar was not found."))
    }
    metadata <- tryCatch(utils::read.delim(paths$metadata,
        check.names = FALSE, stringsAsFactors = FALSE,
        quote = "\"", comment.char = ""),
        error = function(error) error)
    if (inherits(metadata, "error") || !all(c("Key",
        "Value") %in% names(metadata))) {
        return(unavailable("The report table sidecar metadata is invalid."))
    }
    values <- as.character(metadata$Value)
    names(values) <- as.character(metadata$Key)
    required <- c("format_version", "sheet",
        "workbook_md5", "rows", "columns",
        "id_column", "sorted_by_id")
    if (!all(required %in% names(values))) {
        return(unavailable("The report table sidecar metadata is incomplete."))
    }
    if (!identical(values[["format_version"]],
        "1") || !identical(values[["sheet"]],
        sheet) || !identical(values[["sorted_by_id"]],
        "TRUE")) {
        return(unavailable(
            "The report table sidecar metadata is incompatible."))
    }
    workbook_md5 <- unname(tools::md5sum(workbookFile))
    if (!identical(values[["workbook_md5"]],
        workbook_md5)) {
        return(unavailable(
            "The report table sidecar does not match the workbook."))
    }
    list(ok = TRUE, reason = NULL, table = paths$table,
        metadata = paths$metadata, dictionary = paths$dictionary,
        workbookMetadata = paths$workbookMetadata,
        rows = as.integer(values[["rows"]]),
        columns = as.integer(values[["columns"]]),
        idColumn = values[["id_column"]])
}

streamReportTableDnaEpico <- function(tableFile, chunkSize = 5000L,
    expectedRows = NULL, expectedColumns = NULL, chunkHandler = NULL) {
    chunk_size <- max(1L, as.integer(chunkSize))
    connection <- if (grepl("\\.gz$", tableFile, ignore.case = TRUE)) {
        gzfile(tableFile, open = "rt")
    }
    else {
        file(tableFile, open = "rt")
    }
    on.exit(close(connection), add = TRUE)
    header_line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(header_line) == 0L) {
        stop("The report table sidecar is empty.", call. = FALSE)
    }
    columns <- names(utils::read.delim(text = header_line,
        nrows = 0L, check.names = FALSE, quote = "\"",
        comment.char = ""))
    if (!is.null(expectedColumns) && !identical(length(columns),
        as.integer(expectedColumns))) {
        stop("The report table sidecar columns do not match its metadata.",
            call. = FALSE)
    }
    n_rows <- 0L
    chunk_number <- 1L
    maximum_chunk_rows <- 0L
    repeat {
        chunk <- utils::read.delim(connection, header = FALSE,
            col.names = columns, nrows = chunk_size, sep = "\t",
            quote = "\"", comment.char = "", colClasses = "character",
            check.names = FALSE, na.strings = "")
        if (nrow(chunk) == 0L) {
            break
        }
        if (is.function(chunkHandler)) {
            chunkHandler(chunk, chunk_number)
        }
        n_rows <- n_rows + nrow(chunk)
        maximum_chunk_rows <- max(maximum_chunk_rows, nrow(chunk))
        chunk_number <- chunk_number + 1L
    }
    close(connection)
    on.exit(NULL, add = FALSE)
    if (!is.null(expectedRows) && !identical(as.integer(n_rows),
        as.integer(expectedRows))) {
        stop("The report table sidecar row count does not match its metadata.",
            call. = FALSE)
    }
    list(columns = columns, nRows = n_rows, nChunks = chunk_number -
        1L, maximumChunkRows = maximum_chunk_rows)
}

inferMethylationValueLabelMethylationGLM <- function(modelResults) {
    if (!is.null(modelResults$responseLabel) &&
        nzchar(modelResults$responseLabel)) {
        return(modelResults$responseLabel)
    }
    if (!is.null(modelResults$settings$methylationScale)) {
        return(methylationScaleResponseLabelDnaEpico(
            modelResults$settings$methylationScale))
    }
    if (is.null(modelResults$fits) || length(modelResults$fits) ==
        0L) {
        return("Beta values")
    }
    found_response <- FALSE
    for (fit_group in modelResults$fits) {
        if (!is.list(fit_group)) {
            next
        }
        for (fit_object in fit_group) {
            if (!is.list(fit_object) || inherits(fit_object,
                "dnaEPICO_methylationGLM_fit_error") ||
                inherits(fit_object, "dnaEPICO_methylationLME_fit_error") ||
                is.null(fit_object$fitted) ||
                is.null(fit_object$residuals)) {
                next
            }
            fitted_values <- fit_object$fitted
            residual_values <- fit_object$residuals
            if (!is.numeric(fitted_values) ||
                !is.numeric(residual_values) ||
                length(fitted_values) != length(residual_values)) {
                next
            }
            response_values <- fitted_values +
                residual_values
            response_values <- response_values[is.finite(response_values)]
            if (length(response_values) == 0L) {
                next
            }
            found_response <- TRUE
            if (any(response_values < 0 | response_values >
                1, na.rm = TRUE)) {
                return("M-values")
            }
        } }
    if (isTRUE(found_response)) {
        return("Beta values")
    }
    "Beta values"
}

normalizeGlmDataConfigDnaEpico <- function(config) {
    config$logPath <- resolveLogPathMinfiEwasWater(
    logs = config$logs, log_dir = config$log_dir,
    log_file = config$log_file
    )
    config$phenotypeList <- unique(splitOptionMinfiEwasWater(
    config$phenotypes,
    sep = ","
    ))
    config$covariateList <- splitOptionMinfiEwasWater(
    config$covariates,
    sep = ","
    )
    config$factorList <- splitOptionMinfiEwasWater(
    config$factorVars,
    sep = ","
    )
    config$scaleList <- normalizeScaleVariablesDnaEpico(config$scaleVars)
    config$prsMapValue <- parsePrsMapMethylationGLM(config$prsMap)
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
    config$interactionValue <- normalizeOptionalColumnMethylationModels(
    config$interactionTerm, "interactionTerm"
    )
    config
}

loadGlmAnalysisDataDnaEpico <- function(config) {
    object_name <- sub("[.][^.]+$", "", basename(config$inputPheno))
    scale_name <- if (startsWith(object_name, config$objectPrefix)) {
    object_name
    } else {
    character(0)
    }
    legacy_name <- if (identical(config$methylationScaleValue, "beta")) {
    "phenoBT1"
    } else {
    character(0)
    }
    data <- loadSavedObjectPreprocessingPheno(
    config$inputPheno,
    preferred_name = c(
        scale_name, paste0(config$objectPrefix, "T1"),
        legacy_name, object_name
    )
    )
    if (!is.data.frame(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    }
    data
}

validateGlmDataColumnsDnaEpico <- function(data, config) {
    if (!length(config$phenotypeList)) {
    stop("At least one phenotype must be supplied.", call. = FALSE)
    }
    checks <- list(
    list(
        values = config$phenotypeList,
        message = "Phenotype columns not found in inputPheno: %s"
    ),
    list(
        values = config$covariateList,
        message = "Covariate columns not found in inputPheno: %s"
    ),
    list(
        values = unname(config$prsMapValue[
        names(config$prsMapValue) %in% config$phenotypeList
        ]),
        message = "PRS columns not found in inputPheno: %s"
    ),
    list(
        values = config$factorList,
        message = "Factor columns not found in inputPheno: %s"
    )
    )
    for (check in checks) {
    missing <- setdiff(check$values, colnames(data))
    if (length(missing)) {
        missing_text <- paste(unique(missing), collapse = ", ")
        stop(sprintf(
        check$message, missing_text
        ), call. = FALSE)
    }
    }
    interaction <- config$interactionValue
    if (!is.null(interaction) && !interaction %in% colnames(data)) {
    stop(
        "interactionTerm column not found in inputPheno: ",
        interaction,
        call. = FALSE
    )
    }
    invisible(NULL)
}

prepareGlmColumnsDnaEpico <- function(data, config) {
    for (variable in intersect(config$factorList, colnames(data))) {
    data[[variable]] <- as.factor(data[[variable]])
    }
    cpgs <- grep(
    paste0("^", escapeRegexMethylationGLM(config$cpgPrefix)),
    colnames(data),
    value = TRUE
    )
    if (!is.na(config$cpgLimitValue)) {
    cpgs <- utils::head(cpgs, config$cpgLimitValue)
    }
    if (!length(cpgs)) {
    stop(
        "No CpG columns were found with prefix '",
        config$cpgPrefix, "'.",
        call. = FALSE
    )
    }
    validateMethylationProbeIdentifiersDnaEpico(
    cpgs, "CpG columns in the GLM input"
    )
    list(data = data, cpgs = cpgs)
}

scaleGlmModelDataDnaEpico <- function(data, cpgs, config) {
    mapped_prs <- unname(config$prsMapValue[
    names(config$prsMapValue) %in% config$phenotypeList
    ])
    eligible <- unique(c(
    config$phenotypeList, config$covariateList,
    mapped_prs, config$interactionValue
    ))
    eligible <- eligible[!is.na(eligible) & nzchar(eligible)]
    scaled_cpgs <- intersect(config$scaleList, cpgs)
    if (length(scaled_cpgs)) {
    scaled_cpg_text <- paste(scaled_cpgs, collapse = ", ")
    stop(sprintf(
        "%s: %s",
        "CpG methylation response columns cannot be listed in scaleVars",
        scaled_cpg_text
    ), call. = FALSE)
    }
    scaleModelVariablesDnaEpico(
    data = data[, setdiff(colnames(data), cpgs), drop = FALSE],
    scaleVars = config$scaleList, factorVars = config$factorList,
    eligibleVars = eligible, protectedVars = cpgs
    )
}

summarizeGlmInputDataDnaEpico <- function(data, config) {
    requested <- unique(c(config$phenotypeList, config$covariateList))
    missing <- vapply(
    requested, function(column) sum(is.na(data[[column]])), integer(1)
    )
    variable_summary <- summary(data[, requested, drop = FALSE])
    interaction_table <- NULL
    if (!is.null(config$interactionValue)) {
    interaction_table <- table(
        data[[config$interactionValue]],
        useNA = "ifany"
    )
    }
    list(
    missing = missing, variables = variable_summary,
    interaction = interaction_table
    )
}

logGlmPreparedDataDnaEpico <- function(data, cpgs, scaling, summary, config) {
    lines <- c(
    "============================================================",
    paste("Loaded phenotype + methylation data from:", config$inputPheno),
    paste("Merged modeling object:     ", config$objectPrefix, "*"),
    paste("Data dimensions:             ", paste(dim(data), collapse = " x ")),
    paste(
        "Phenotypes:                  ",
        paste(config$phenotypeList, collapse = ", ")
    ),
    paste(
        "Covariates:                  ",
        paste(config$covariateList, collapse = ", ")
    ),
    paste(
        "Factor variables:            ",
        paste(config$factorList, collapse = ", ")
    ),
    formatScalingMetadataLogDnaEpico(scaling$metadata),
    paste("CpG columns retained:        ", length(cpgs)),
    "Missing summary:",
    paste(names(summary$missing), summary$missing,
        sep = ": ", collapse = "; "
    ),
    "Summary statistics:", previewLinesMinfiEwasWater(summary$variables)
    )
    if (!is.null(summary$interaction)) {
    lines <- c(
        lines,
        paste("Interaction table for", config$interactionValue, ":"),
        previewLinesMinfiEwasWater(summary$interaction)
    )
    }
    emitLogMinfiEwasWater(
    c(lines, "============================================================"),
    verbose = config$verbose, log_path = config$logPath
    )
}

newMethylationGLMDataDnaEpico <- function(
    data, cpgs, scaling, summary, config
) {
    structure(list(
    data = data, modelData = scaling$data,
    inputPheno = config$inputPheno,
    inputIdentity = inputIdentityMethylationModels(config$inputPheno),
    phenotypes = config$phenotypeList,
    covariates = config$covariateList,
    factorVars = config$factorList, scaleVars = scaling$scaleVars,
    scalingMetadata = scaling$metadata, cpgColumns = cpgs,
    cpgPrefix = config$cpgPrefix, cpgLimit = config$cpgLimitValue,
    methylationScale = config$methylationScaleValue,
    responseLabel = config$responseLabel,
    methylationObjectPrefix = config$objectPrefix,
    internalResponseColumn = config$responseColumn,
    prsMap = config$prsMapValue,
    interactionTerm = config$interactionValue,
    requestedInteractionTerm = config$interactionTerm,
    missingCounts = summary$missing,
    variableSummary = summary$variables,
    interactionTable = summary$interaction
    ), class = "dnaEPICO_methylationGLM_data")
}

#' Prepare phenotype-plus-methylation data for one-timepoint GLM analyses
#'
#' @param inputPheno Character. Path to the merged phenotype-plus-methylation
#' object
#'   created by `preprocessingPheno()`.
#' @param phenotypes Character vector or comma-separated string of phenotype
#'   variables to model.
#' @param covariates Character vector or comma-separated string of covariate
#'   variables to adjust for.
#' @param factorVars Character vector or comma-separated string of variables
#' that
#'   should be converted to factors before modeling.
#' @param scaleVars Character vector, comma-separated variable names, or `NULL`.
#'   Numeric fixed-effect variables to standardize before model fitting.
#' @param cpgPrefix Character. Prefix used to identify methylation columns.
#' @param cpgLimit Integer or `NA`. Maximum number of CpGs to retain. `NA`
#'   keeps all matching CpGs.
#' @param methylationScale Character. Methylation metric represented by the CpG
#'   columns. One of `'Beta'`, `'M'`, or `'CN'`, case-insensitively.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#' @param prsMap Character vector or comma-separated string of phenotype-to-PRS
#'   mappings in the form `'Phenotype:PRS'`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationGLM_data'` containing the
#'   prepared analysis data, parsed variable selections, CpG columns, and
#'   exploratory summaries.
#'
#' @description
#' Load the merged phenotype-plus-methylation input object, validate the
#' requested
#' modeling variables, convert selected variables to factors, and return a
#' single in-memory object for downstream helpers.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' prepared_data <- prepareMethylationGLMData(
#'   inputPheno = ex$inputPath,
#'   phenotypes = "status",
#'   covariates = "sex,age",
#'   factorVars = "status,sex",
#'   cpgLimit = 2,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(prepared_data)
#'
#' @export
prepareMethylationGLMData <- function(
    inputPheno, phenotypes, covariates, factorVars, scaleVars = NULL,
    cpgPrefix = "cg", cpgLimit = NA, methylationScale = "beta",
    interactionTerm = NULL, prsMap = NULL,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_methylationGLM.txt"
) {
    config <- normalizeGlmDataConfigDnaEpico(
    as.list(environment(), all.names = TRUE)
    )
    data <- loadGlmAnalysisDataDnaEpico(config)
    validateGlmDataColumnsDnaEpico(data, config)
    columns <- prepareGlmColumnsDnaEpico(data, config)
    data <- columns$data
    scaling <- scaleGlmModelDataDnaEpico(data, columns$cpgs, config)
    summary <- summarizeGlmInputDataDnaEpico(data, config)
    logGlmPreparedDataDnaEpico(
    data, columns$cpgs, scaling, summary, config
    )
    newMethylationGLMDataDnaEpico(
    data, columns$cpgs, scaling, summary, config
    )
}
plotGlmDistributionVariableDnaEpico <- function(
    data, variable, type, fill, outputDir, display,
    width, height, resolution
) {
    if (identical(type, "auto")) {
    type <- if (is.numeric(data[[variable]])) "hist" else "bar"
    }
    plot <- createDistributionPlotMethylationGLM(
    values = data[[variable]], variable = variable,
    type = type, fill = fill
    )
    suffix <- if (identical(type, "hist")) "continuous" else "categorical"
    file <- if (is.null(outputDir)) {
    NULL
    } else {
    file.path(outputDir, paste0(
        "distribution_", safeFigureComponentDnaEpico(variable),
        "_", suffix, ".tiff"
    ))
    }
    runPlotMinfiEwasWater(
    draw_fun = function() drawPlotObjectMinfiEwasWater(plot),
    display = display, file = file, width = width,
    height = height, res = resolution
    )
    list(plot = plot, file = file)
}

plotGlmDistributionVariablesDnaEpico <- function(
    variables, data, type, fill, outputDir, display,
    width, height, resolution
) {
    plots <- list()
    files <- list()
    for (variable in variables) {
    result <- plotGlmDistributionVariableDnaEpico(
        data, variable, type, fill, outputDir, display,
        width, height, resolution
    )
    plots[[variable]] <- result$plot
    if (!is.null(result$file)) {
        files[[variable]] <- result$file
    }
    }
    list(plots = plots, files = files)
}

#' Plot phenotype and covariate distributions for one-timepoint GLM analyses
#'
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
#' @param plotWidth Integer. TIFF width in pixels when plots are written to
#' disk.
#' @param plotHeight Integer. TIFF height in pixels when plots are written to
#'   disk.
#' @param plotDPI Integer. TIFF resolution in DPI when plots are written to
#'   disk.
#' @param outputDir Character or `NULL`. Directory used for TIFF files. When
#'   `NULL`, plots are returned in memory only.
#' @param display Logical. If `TRUE`, draw plots on the active graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationGLM_distribution_plots'`
#'   containing the generated `ggplot2` objects and any saved TIFF file paths.
#'
#' @description
#' Create phenotype, factor-variable, and numeric-covariate distribution plots
#' from the object returned by `prepareMethylationGLMData()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' distribution_plots <- plotMethylationGLMDistributions(
#'   preparedData = ex$preparedData,
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(distribution_plots)
#'
#' @export
plotMethylationGLMDistributions <- function(preparedData,
    plotWidth = 2000L, plotHeight = 1000L, plotDPI = 150L,
    outputDir = NULL, display = FALSE, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationGLM.txt") {
    log_path <- resolveLogPathMinfiEwasWater(logs = logs,
        log_dir = log_dir, log_file = log_file)
    data <- preparedData$data
    factor_vars <- setdiff(intersect(preparedData$factorVars,
        colnames(data)), preparedData$phenotypes)
    numeric_vars <- setdiff(preparedData$covariates, c(preparedData$factorVars,
        preparedData$phenotypes))
    numeric_vars <- intersect(numeric_vars, colnames(data))
    interaction_vars <- preparedData$interactionTerm
    interaction_vars <- interaction_vars[!is.na(interaction_vars) &
        nzchar(interaction_vars)]
    interaction_vars <- setdiff(intersect(interaction_vars,
        colnames(data)), c(preparedData$phenotypes, factor_vars,
        numeric_vars))
    common <- list(data = data, outputDir = outputDir, display = display,
        width = plotWidth, height = plotHeight, resolution = plotDPI)
    phenotype <- do.call(plotGlmDistributionVariablesDnaEpico,
        c(list(variables = intersect(preparedData$phenotypes,
            colnames(data)), type = "auto", fill = "steelblue"),
            common))
    factor <- do.call(plotGlmDistributionVariablesDnaEpico,
        c(list(variables = factor_vars, type = "bar", fill = "darkorange"),
            common))
    numeric <- do.call(plotGlmDistributionVariablesDnaEpico,
        c(list(variables = numeric_vars, type = "hist", fill = "darkgreen"),
            common))
    interaction <- do.call(plotGlmDistributionVariablesDnaEpico,
        c(list(variables = interaction_vars, type = "auto",
            fill = "#7A5195"), common))
    covariate_plots <- c(numeric$plots, interaction$plots)
    covariate_files <- c(numeric$files, interaction$files)
    emitLogMinfiEwasWater(c(
        "============================================================",
        paste("Phenotype distribution plots: ", length(phenotype$plots)),
        paste("Factor distribution plots:    ", length(factor$plots)),
        paste("Numeric covariate plots:      ", length(covariate_plots)),
        if (is.null(outputDir)) {
            "Distribution plots were returned in memory only."
        } else {
            paste("Distribution plots saved to:  ", outputDir)
        }, "============================================================"),
        verbose = verbose, log_path = log_path)
    structure(list(phenotypes = phenotype$plots, factors = factor$plots,
        covariates = covariate_plots, files = list(phenotypes = phenotype$files,
            factors = factor$files, covariates = covariate_files)),
        class = "dnaEPICO_methylationGLM_distribution_plots") }

normalizeGlmFitConfigDnaEpico <- function(config) {
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
    config$glmLibraries <- splitOptionMinfiEwasWater(
    config$glmLibs,
    sep = ","
    )
    if (!length(config$glmLibraries)) {
    config$glmLibraries <- "glm2"
    }
    config$omnibus <- validateOmnibusConfigurationMethylationGLM(
    config$omnibusTest
    )
    config$requiredPackages <- unique(c(
    config$glmLibraries,
    if (config$omnibus) "car" else character(0)
    ))
    config$nCoresValue <- validatePositiveIntegerMethylationModels(
    config$nCores, "nCores"
    )
    config
}

newGlmFitStateDnaEpico <- function(preparedData, modelData, config) {
    state <- new.env(parent = emptyenv())
    list_fields <- c(
    "fits", "summaryCache", "coefficientResults", "phenotypeSummaries",
    "summaryFiles", "modelMessages", "fitFailures", "failureReasons",
    "omnibusTests"
    )
    for (field in list_fields) {
    state[[field]] <- list()
    }
    state$resumedPhenotypes <- character(0)
    state$fittedPhenotypes <- character(0)
    state$formulas <- stats::setNames(
    character(length(preparedData$phenotypes)), preparedData$phenotypes
    )
    state$failureCounts <- stats::setNames(
    integer(length(preparedData$phenotypes)), preparedData$phenotypes
    )
    state$omnibusTargets <- stats::setNames(
    character(length(preparedData$phenotypes)), preparedData$phenotypes
    )
    state$parallelPlan <- resolveParallelPlanMethylationModels(
    engine = "glm2", nCores = config$nCoresValue,
    nCpGs = length(preparedData$cpgColumns),
    analysisData = preparedData$data, modelData = modelData
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

stopGlmFitClusterDnaEpico <- function(state) {
    if (!is.null(state$cluster)) {
    try(parallel::stopCluster(state$cluster), silent = TRUE)
    state$cluster <- NULL
    }
    invisible(NULL)
}

glmPhenotypeVariablesDnaEpico <- function(
    preparedData, modelData, phenotype
) {
    prs <- if (phenotype %in% names(preparedData$prsMap)) {
    unname(preparedData$prsMap[[phenotype]])
    } else {
    character(0)
    }
    covariates <- unique(c(preparedData$covariates, prs))
    variables <- unique(c(
    phenotype, covariates, preparedData$interactionTerm
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

prepareGlmPhenotypeSpecDnaEpico <- function(
    preparedData, modelData, phenotype, config
) {
    variables <- glmPhenotypeVariablesDnaEpico(
    preparedData, modelData, phenotype
    )
    formula_text <- buildFormulaMethylationGLM(
    phenotype = phenotype, covariates = variables$covariates,
    interactionTerm = preparedData$interactionTerm,
    responseVar = preparedData$internalResponseColumn
    )
    data <- modelData[, variables$variables, drop = FALSE]
    for (var in intersect(preparedData$factorVars, colnames(data))) {
    data[[var]] <- as.factor(data[[var]])
    }
    validateFixedEffectDesignMethylationModels(formula_text, data)
    formula_object <- stats::as.formula(formula_text, env = baseenv())
    coefficient_terms <- buildCoefficientTermMapMethylationModels(
    formula_text, data
    )
    omnibus_target <- if (config$omnibus) {
    resolveOmnibusTargetTermMethylationGLM(
        formula_text, data, phenotype, preparedData$interactionTerm
    )
    } else {
    NULL
    }
    list(
    phenotype = phenotype, covariates = variables$covariates,
    data = data, formula = formula_text, formulaObject = formula_object,
    coefficientTerms = coefficient_terms, omnibusTarget = omnibus_target
    )
}

glmPhenotypeSignatureDnaEpico <- function(preparedData, spec, config) {
    buildPhenotypeSignatureMethylationModels(
    analysis = "glm", engine = "glm2", phenotype = spec$phenotype,
    formulaText = spec$formula, preparedData = preparedData,
    modelSettings = list(
        family = "gaussian", link = "identity",
        omnibusTest = config$omnibus,
        omnibusTerm = spec$omnibusTarget,
        omnibusMethod = "car::linearHypothesis Wald F",
        omnibusRhs = 0, omnibusJoint = TRUE
    ),
    packages = config$requiredPackages
    )
}

glmPhenotypeSummaryPathDnaEpico <- function(phenotype, config) {
    if (is.null(config$summaryDir)) {
    return(NULL)
    }
    phenotypeSummaryPathMethylationModels(
    outputDir = config$summaryDir, phenotype = phenotype,
    analysis = "glm"
    )
}

logResumedGlmPhenotypeDnaEpico <- function(
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

restoreGlmPhenotypeDnaEpico <- function(
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
    logResumedGlmPhenotypeDnaEpico(
    phenotype, artifact$formula, summaryPath, nCpGs, config
    )
    invisible(TRUE)
}

tryResumeGlmPhenotypeDnaEpico <- function(
    state, spec, signature, summaryPath, config, nCpGs
) {
    resumed <- list(object = NULL, reason = "resume was not requested")
    if (config$resumeFromSummary && !is.null(summaryPath)) {
    resumed <- loadPhenotypeSummaryMethylationModels(
        summaryPath, signature
    )
    }
    if (!is.null(resumed$object)) {
    restoreGlmPhenotypeDnaEpico(
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

glmBatchArgumentsDnaEpico <- function(preparedData, spec, config) {
    list(
    modelData = spec$data, formulaText = spec$formula,
    phenotype = spec$phenotype,
    interactionTerm = preparedData$interactionTerm,
    responseVar = preparedData$internalResponseColumn,
    omnibusTest = config$omnibus,
    omnibusTerm = spec$omnibusTarget,
    formulaObject = spec$formulaObject,
    coefficientTerms = spec$coefficientTerms
    )
}

updateGlmPilotPlanDnaEpico <- function(
    state, preparedData, spec, config
) {
    cpgs <- preparedData$cpgColumns
    if (state$pilotCompleted || state$workerCount <= 1L || !length(cpgs)) {
    return(invisible(NULL))
    }
    common <- glmBatchArgumentsDnaEpico(preparedData, spec, config)
    pilot_cpgs <- utils::head(cpgs, 3L)
    pilot <- measurePilotMemoryMethylationModels(function() {
    do.call(fitMethylationGLMBatch, c(list(
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

glmPsockDependencyNamesDnaEpico <- function() {
    c(
    "validateWorkerPackagesMethylationModels",
    "newMethylationFitErrorDnaEpico", "captureModelConditionsDnaEpico",
    "combineModelMessagesDnaEpico", "modelMessageDnaEpico",
    "fitMethylationGLMBatch", "fitCpGModelMethylationGLM",
    "fitCpGGlmValueDnaEpico", "completeCapturedGlmFitDnaEpico",
    "buildCoefficientTermMapMethylationModels",
    "removeRandomInterceptMethylationModels",
    "computeOmnibusTestMethylationGLM",
    "emptyOmnibusResultMethylationGLM",
    "summarizeCpGFitMethylationGLM", "findCoefficientRowsMethylationGLM",
    "escapeRegexMethylationGLM", "cpgResponseMethylationModels",
    "compactCoefficientResultsMethylationModels",
    "coefficientNamesMethylationModels",
    "newCompactCoefficientStorageMethylationModels",
    "fillCompactCoefficientStorageMethylationModels",
    "collectBatchModelMessagesMethylationModels",
    "collectBatchFitFailuresMethylationModels",
    "emptyModelMessagesDnaEpico", "emptyFitFailuresMethylationModels",
    "collectOmnibusTestsMethylationGLM"
    )
}

ensureGlmPsockClusterDnaEpico <- function(state, config) {
    if (!identical(state$backend, "psock") || length(state$batches) <= 1L ||
    !is.null(state$cluster)) {
    return(invisible(NULL))
    }
    state$cluster <- makePsockClusterMethylationModels(min(
    state$workerCount, length(state$batches)
    ))
    parallel::clusterExport(
    state$cluster, glmPsockDependencyNamesDnaEpico(),
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

runSerialGlmBatchesDnaEpico <- function(
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

runForkGlmBatchesDnaEpico <- function(
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

runPsockGlmBatchesDnaEpico <- function(
    state, preparedData, common, worker
) {
    state$clusterUseCount <- state$clusterUseCount + 1L
    parallel::clusterExport(
    state$cluster, c("common", "worker"),
    envir = environment()
    )
    results <- vector("list", length(state$batches))
    size <- min(state$workerCount, length(state$batches))
    waves <- split(
    seq_along(state$batches),
    ceiling(seq_along(state$batches) / size)
    )
    for (wave in waves) {
    tasks <- lapply(wave, function(index) {
        batch <- state$batches[[index]]
        list(
        cpgBatch = batch,
        responses = as.matrix(
            preparedData$data[, batch, drop = FALSE]
        )
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

runGlmPhenotypeBatchesDnaEpico <- function(
    state, preparedData, spec, config
) {
    common <- glmBatchArgumentsDnaEpico(preparedData, spec, config)
    worker <- fitMethylationGLMBatch
    ensureGlmPsockClusterDnaEpico(state, config)
    if (identical(state$backend, "fork") && length(state$batches) > 1L) {
    return(runForkGlmBatchesDnaEpico(
        state$batches, preparedData, common, config,
        worker, state$workerCount
    ))
    }
    if (identical(state$backend, "psock") && length(state$batches) > 1L) {
    environment(worker) <- .GlobalEnv
    return(runPsockGlmBatchesDnaEpico(
        state, preparedData, common, worker
    ))
    }
    runSerialGlmBatchesDnaEpico(
    state$batches, preparedData, common, config, worker
    )
}

combineGlmBatchResultsDnaEpico <- function(
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
    collectOmnibusTestsMethylationGLM(list(), phenotype)
    )
    summary_cache <- filterSummaryByPvalueMethylationGLM(
    summaries, NA_real_,
    includeResidualSD = TRUE
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

glmPhenotypeArtifactDnaEpico <- function(
    preparedData, spec, combined, signature, config
) {
    factor_vars <- preparedData$factorVars[
    preparedData$factorVars %in% colnames(preparedData$data)
    ]
    assemblePhenotypeSummaryMethylationModels(
    analysis = "glm", engine = "glm2", phenotype = spec$phenotype,
    signature = signature, cpgOrder = preparedData$cpgColumns,
    coefficientResults = combined$coefficients,
    targetSummary = combined$summaryCache,
    omnibusTests = combined$omnibus,
    modelMessages = combined$messages, fitFailures = combined$failures,
    failureCount = combined$failureCount,
    failureReasons = combined$errors, formulaText = spec$formula,
    settings = list(
        methylationScale = preparedData$methylationScale,
        responseLabel = preparedData$responseLabel,
        interactionTerm = preparedData$interactionTerm,
        covariates = spec$covariates, factorVars = preparedData$factorVars,
        factorLevels = lapply(preparedData$data[factor_vars], levels),
        scaleVars = preparedData$scaleVars,
        scalingMetadata = preparedData$scalingMetadata,
        omnibusTest = config$omnibus,
        omnibusTerm = spec$omnibusTarget,
        omnibusMethod = "car::linearHypothesis Wald F",
        omnibusRhs = 0, omnibusJoint = TRUE
    )
    )
}

storeFittedGlmPhenotypeDnaEpico <- function(
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
    state$formulas[[phenotype]] <- spec$formula
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

glmFitResourceLogLinesDnaEpico <- function(state) {
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

logFittedGlmPhenotypeDnaEpico <- function(
    state, preparedData, spec, combined, config
) {
    emitLogMinfiEwasWater(c(
    "============================================================",
    paste("Fitted phenotype:            ", spec$phenotype),
    paste("Formula:                     ", spec$formula),
    paste("CpGs attempted:              ", length(preparedData$cpgColumns)),
    paste("CpGs without p-values:       ", combined$failureCount),
    paste("Omnibus tests requested:     ", config$omnibus),
    paste(
        "Omnibus target term:         ",
        if (is.null(spec$omnibusTarget)) "None" else spec$omnibusTarget
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
        "Top fit errors:               ",
        formatFitErrorsMethylationModels(combined$errors)
    ),
    glmFitResourceLogLinesDnaEpico(state),
    paste("Fit-time summary rows cached:", nrow(combined$summaryCache)),
    "============================================================"
    ), verbose = config$verbose, log_path = config$logPath)
}

warnUnavailableGlmPvaluesDnaEpico <- function(
    preparedData, spec, combined
) {
    if (!length(preparedData$cpgColumns) || any(combined$pValueAvailable)) {
    return(invisible(NULL))
    }
    warning(sprintf(
    "%s '%s'. %s: %s. %s",
    "No CpG GLM p-values were available for phenotype", spec$phenotype,
    "Top failure reasons",
    formatFitErrorsMethylationModels(combined$errors),
    "The failure inventory was retained and the analysis continued."
    ), call. = FALSE)
    invisible(NULL)
}

combineGlmFitAuditTablesDnaEpico <- function(state) {
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

glmFitSettingsDnaEpico <- function(state, preparedData, config) {
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
    glmLibs = config$glmLibraries,
    methylationScale = preparedData$methylationScale,
    methylationObjectPrefix = preparedData$methylationObjectPrefix,
    responseLabel = preparedData$responseLabel,
    internalResponseColumn = preparedData$internalResponseColumn,
    interactionTerm = preparedData$interactionTerm,
    omnibusTest = config$omnibus,
    omnibusMethod = "car::linearHypothesis Wald F",
    omnibusRhs = 0, omnibusJoint = TRUE,
    phenotypes = preparedData$phenotypes,
    covariates = preparedData$covariates,
    factorVars = preparedData$factorVars,
    factorLevels = lapply(stats::setNames(
        preparedData$factorVars, preparedData$factorVars
    ), function(variable) levels(preparedData$data[[variable]])),
    scaleVars = preparedData$scaleVars,
    scalingMetadata = preparedData$scalingMetadata,
    sampleCount = nrow(preparedData$data)
    )
}

newMethylationGLMModelsDnaEpico <- function(state, preparedData, config) {
    audit <- combineGlmFitAuditTablesDnaEpico(state)
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
    settings = glmFitSettingsDnaEpico(state, preparedData, config),
    responseLabel = preparedData$responseLabel
    ), class = "dnaEPICO_methylationGLM_models")
}

fitAllGlmPhenotypesDnaEpico <- function(
    state, preparedData, modelData, config
) {
    for (phenotype in preparedData$phenotypes) {
    spec <- prepareGlmPhenotypeSpecDnaEpico(
        preparedData, modelData, phenotype, config
    )
    if (config$omnibus) {
        state$omnibusTargets[[phenotype]] <- spec$omnibusTarget
    }
    signature <- glmPhenotypeSignatureDnaEpico(
        preparedData, spec, config
    )
    summary_path <- glmPhenotypeSummaryPathDnaEpico(phenotype, config)
    if (tryResumeGlmPhenotypeDnaEpico(
        state, spec, signature, summary_path, config,
        length(preparedData$cpgColumns)
    )) {
        next
    }
    updateGlmPilotPlanDnaEpico(state, preparedData, spec, config)
    batches <- runGlmPhenotypeBatchesDnaEpico(
        state, preparedData, spec, config
    )
    combined <- combineGlmBatchResultsDnaEpico(
        batches, preparedData$cpgColumns, phenotype
    )
    artifact <- glmPhenotypeArtifactDnaEpico(
        preparedData, spec, combined, signature, config
    )
    storeFittedGlmPhenotypeDnaEpico(
        state, spec, combined, artifact, summary_path
    )
    logFittedGlmPhenotypeDnaEpico(
        state, preparedData, spec, combined, config
    )
    warnUnavailableGlmPvaluesDnaEpico(preparedData, spec, combined)
    }
    invisible(NULL)
}

#' Fit CpG-wise Gaussian GLMs for one-timepoint methylation analyses
#'
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
#' @param nCores Integer. Maximum number of worker processes to use. Automatic
#'   fitting remains serial below the glm2 crossover and caps workers by
#'   available CpGs, CPUs, and detected memory.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded
#'   to worker processes.
#' @param glmLibs Character vector or comma-separated string of package names to
#'   check on worker processes. The default is `'glm2'`.
#' @param omnibusTest Logical. If `TRUE`, use `car::linearHypothesis()` to test
#'   the complete phenotype-by-interaction term, or the phenotype main effect
#'   when no interaction is specified, once per CpG.
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
#' @return A list with class `'dnaEPICO_methylationGLM_models'` containing
#'   compact coefficient matrices, unfiltered target summaries, formulas,
#'   model conditions, hard errors, and phenotype summary artifacts.
#'
#' @description
#' Fit one Gaussian GLM per CpG for each phenotype requested in the object
#' returned by `prepareMethylationGLMData()`. Each native fit is reduced to
#' compact numerical results and discarded before the next batch is returned.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' model_results <- fitMethylationGLMModels(
#'   preparedData = ex$preparedData,
#'   nCores = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(model_results$fits)
#'
#' @export
fitMethylationGLMModels <- function(
    preparedData, nCores = 1L,
    libPath = NULL, glmLibs = "glm2", summaryDir = NULL,
    omnibusTest = FALSE, resumeFromSummary = TRUE,
    verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationGLM.txt"
) {
    config <- normalizeGlmFitConfigDnaEpico(
    as.list(environment(), all.names = TRUE)
    )
    model_data <- if (is.null(preparedData$modelData)) {
    preparedData$data
    } else {
    preparedData$modelData
    }
    state <- newGlmFitStateDnaEpico(preparedData, model_data, config)
    on.exit(stopGlmFitClusterDnaEpico(state), add = TRUE)
    fitAllGlmPhenotypesDnaEpico(
    state, preparedData, model_data, config
    )
    stopGlmFitClusterDnaEpico(state)
    newMethylationGLMModelsDnaEpico(state, preparedData, config)
}

summarizeOmnibusTestsMethylationGLM <- function(
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

normalizeGlmSummaryConfigDnaEpico <- function(config) {
    config$summaryResidualSD <- isTRUE(config$summaryResidualSD)
    config$logPath <- resolveLogPathMinfiEwasWater(
    logs = config$logs, log_dir = config$log_dir,
    log_file = config$log_file
    )
    if (is.null(config$libPath)) {
    config$libPath <- .libPaths()
    }
    config$glmLibraries <- splitOptionMinfiEwasWater(
    config$glmLibs,
    sep = ","
    )
    if (!length(config$glmLibraries)) {
    config$glmLibraries <- "glm2"
    }
    config$pValueFilter <- validateProbabilityDnaEpico(
    normalizeOptionalNumericMethylationGLM(config$summaryPval),
    "summaryPval",
    allowNA = TRUE
    )
    config$chunkSizeValue <- normalizeChunkSizeMethylationGLM(
    config$chunkSize
    )
    config$adjustmentMethod <- validatePAdjustmentMethodMethylationModels(
    config$padjmethod
    )
    config$nCoresValue <- validatePositiveIntegerMethylationModels(
    config$nCores, "nCores"
    )
    config
}

summarizeGlmChunkDnaEpico <- function(
    chunk, fits, phenotype, interactionTerm, includeResidualSD, worker
) {
    rows <- lapply(chunk, function(cpg) {
    worker(
        cpg = cpg, modelObj = fits[[cpg]], variable = phenotype,
        interactionTerm = interactionTerm,
        includeResidualSD = includeResidualSD
    )
    })
    rows <- Filter(Negate(is.null), rows)
    if (!length(rows)) NULL else do.call(rbind, rows)
}

runGlmSummaryChunksDnaEpico <- function(fits, phenotype, interactionTerm,
    includeResidualSD, nCores, chunkSize, libPath, packages) {
    cpg_names <- names(fits)
    local_size <- chunkSize
    if (is.null(local_size)) {
        local_size <- max(10L, floor(length(cpg_names)/max(nCores *
            4L, 1L)))
    }
    local_size <- max(1L, as.integer(local_size))
    chunks <- split(cpg_names, ceiling(seq_along(cpg_names)/local_size))
    worker <- summarizeCpGFitMethylationGLM
    chunk_worker <- summarizeGlmChunkDnaEpico
    run_chunk <- function(chunk) {
        chunk_worker(chunk, fits, phenotype, interactionTerm,
            includeResidualSD, worker)
    }
    if (nCores > 1L && length(chunks) > 1L) {
        environment(worker) <- .GlobalEnv
        environment(chunk_worker) <- .GlobalEnv
        cluster <- parallel::makeCluster(min(nCores, length(chunks)))
        on.exit(parallel::stopCluster(cluster), add = TRUE)
        dependencies <- c("fits", "phenotype", "interactionTerm",
            "includeResidualSD", "worker", "chunk_worker",
            "summarizeGlmChunkDnaEpico", "findCoefficientRowsMethylationGLM",
            "modelMessageDnaEpico")
        parallel::clusterExport(cluster, dependencies, envir = environment())
        parallel::clusterCall(cluster, validateWorkerPackagesMethylationModels,
            libPath, packages)
        results <- parallel::parLapplyLB(cluster, chunks, run_chunk)
    }
    else {
        results <- lapply(chunks, run_chunk)
    }
    results <- Filter(Negate(is.null), results)
    summary <- if (!length(results)) {
        data.frame()
    }
    else {
        do.call(rbind, results)
    }
    if (nrow(summary)) {
        columns <- c("CpG", "Coefficient", "Estimate", "Std. Error",
            "t value", "Pr(>|t|)", if (includeResidualSD) "ResidualSD",
            "Model.Message")
        summary <- summary[, columns, drop = FALSE]
    }
    rownames(summary) <- NULL
    list(summary = summary, chunkSize = local_size)
}

summarizeGlmPhenotypeDnaEpico <- function(
    modelResults, preparedData, phenotype, config
) {
    cached <- modelResults$summaryCache[[phenotype]]
    if (!is.null(cached)) {
    diagnostic <- filterSummaryByPvalueMethylationGLM(
        cached, NA_real_, config$summaryResidualSD
    )
    summary <- filterSummaryByPvalueMethylationGLM(
        diagnostic, config$pValueFilter, config$summaryResidualSD
    )
    return(list(
        summary = summary, diagnostic = diagnostic,
        source = "fit-time cache", chunkSize = NULL
    ))
    }
    chunked <- runGlmSummaryChunksDnaEpico(
    fits = modelResults$fits[[phenotype]], phenotype = phenotype,
    interactionTerm = preparedData$interactionTerm,
    includeResidualSD = config$summaryResidualSD,
    nCores = config$nCoresValue, chunkSize = config$chunkSizeValue,
    libPath = config$libPath, packages = config$glmLibraries
    )
    diagnostic <- chunked$summary
    summary <- filterSummaryByPvalueMethylationGLM(
    diagnostic, config$pValueFilter, config$summaryResidualSD
    )
    list(
    summary = summary, diagnostic = diagnostic,
    source = "model fits", chunkSize = chunked$chunkSize
    )
}

logGlmPhenotypeSummaryDnaEpico <- function(phenotype, result, config) {
    lines <- c(
    "============================================================",
    paste("Summarized phenotype:        ", phenotype),
    paste("CpG summary rows returned:   ", nrow(result$summary)),
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

newMethylationGLMSummariesDnaEpico <- function(
    summaries, diagnostics, omnibus, modelResults, preparedData, config
) {
    messages <- modelResults$modelMessages
    if (is.null(messages)) {
    messages <- collectModelMessagesMethylationGLM(modelResults$fits)
    }
    structure(list(
    summaries = summaries, diagnosticSummaries = diagnostics,
    omnibusTests = omnibus, phenotypes = names(summaries),
    fitFailures = modelResults$fitFailures, modelMessages = messages,
    settings = list(
        summaryResidualSD = config$summaryResidualSD,
        summaryPval = config$pValueFilter,
        padjmethod = config$adjustmentMethod,
        chunkSize = config$chunkSizeValue,
        interactionTerm = preparedData$interactionTerm,
        omnibusTest = isTRUE(modelResults$settings$omnibusTest),
        factorLevels = modelResults$settings$factorLevels
    )
    ), class = "dnaEPICO_methylationGLM_summaries")
}

#' Summarize CpG-wise Gaussian GLM results for one-timepoint analyses
#'
#' @param modelResults Object returned by `fitMethylationGLMModels()`.
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
#' @param summaryResidualSD Logical. If `TRUE`, add residual standard deviations
#'   to each CpG summary row.
#' @param summaryPval Numeric or `NA`. Optional p-value filter applied to the
#'   returned summary tables. `NA` keeps all rows.
#' @param padjmethod Character. Adjustment method passed to `stats::p.adjust()`
#'   for omnibus p-values across CpGs within each phenotype and tested term.
#' @param nCores Integer. Number of worker processes to use while extracting
#'   summary rows.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded
#'   to worker processes.
#' @param glmLibs Character vector or comma-separated string of package names to
#'   check on worker processes. The default is `'glm2'`.
#' @param chunkSize Integer or `NULL`. Number of CpGs to process per parallel
#'   chunk. `NULL` chooses a value automatically.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationGLM_summaries'`
#'   containing the optionally filtered summary tables in `summaries` and the
#'   complete CpG-level tables in `diagnosticSummaries`. Diagnostics,
#'   annotation, and report output use the complete tables so `summaryPval`
#'   does not remove CpGs from those outputs. `modelMessages` retains native
#'   messages, warnings, and errors for every attempted CpG.
#'
#' @description
#' Return phenotype-specific CpG coefficient tables from the compact fit-time
#' results produced by `fitMethylationGLMModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' summary_results <- summarizeMethylationGLMModels(
#'   modelResults = ex$modelResults,
#'   preparedData = ex$preparedData,
#'   summaryResidualSD = TRUE,
#'   summaryPval = NA,
#'   nCores = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(summary_results$summaries)
#'
#' @export
summarizeMethylationGLMModels <- function(
    modelResults, preparedData,
    summaryResidualSD = TRUE, summaryPval = NA, padjmethod = "fdr",
    nCores = 1L, libPath = NULL, glmLibs = "glm2", chunkSize = NULL,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_methylationGLM.txt"
) {
    config <- normalizeGlmSummaryConfigDnaEpico(
    as.list(environment(), all.names = TRUE)
    )
    summaries <- list()
    diagnostics <- list()
    for (phenotype in names(modelResults$fits)) {
    result <- summarizeGlmPhenotypeDnaEpico(
        modelResults, preparedData, phenotype, config
    )
    summaries[[phenotype]] <- result$summary
    diagnostics[[phenotype]] <- result$diagnostic
    logGlmPhenotypeSummaryDnaEpico(phenotype, result, config)
    }
    omnibus <- summarizeOmnibusTestsMethylationGLM(
    modelResults, config$adjustmentMethod
    )
    newMethylationGLMSummariesDnaEpico(
    summaries, diagnostics, omnibus,
    modelResults, preparedData, config
    )
}

glmCoefficientTableForCpgDnaEpico <- function(modelResults, phenotype, cpg) {
    table <- coefficientTableFromCompactMethylationModels(
    modelResults$coefficientResults[[phenotype]], cpg
    )
    if (!is.null(table)) {
    return(table)
    }
    fit <- modelResults$fits[[phenotype]][[cpg]]
    if (is.null(fit) || inherits(
    fit, "dnaEPICO_methylationGLM_fit_error"
    ) || is.null(fit$coef)) {
    return(NULL)
    }
    as.data.frame(fit$coef)
}

collectOmnibusGlmHitsDnaEpico <- function(
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
    glmCoefficientTableForCpgDnaEpico(modelResults, phenotype, cpg)
    })
    keep <- !vapply(hits, is.null, logical(1))
    stats::setNames(Filter(Negate(is.null), hits), unique(hit_cpgs)[keep])
}

collectCachedGlmHitsDnaEpico <- function(
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
    hit_cpgs <- unique(cached$CpG[cached[["Pr(>|t|)"]] < threshold])
    hit_cpgs <- hit_cpgs[!is.na(hit_cpgs)]
    }
    hits <- lapply(hit_cpgs, function(cpg) {
    glmCoefficientTableForCpgDnaEpico(modelResults, phenotype, cpg)
    })
    keep <- !vapply(hits, is.null, logical(1))
    stats::setNames(Filter(Negate(is.null), hits), hit_cpgs[keep])
}

collectDirectGlmHitsDnaEpico <- function(
    modelResults, phenotype, threshold, interactionTerm
) {
    fits <- modelResults$fits[[phenotype]]
    hits <- list()
    for (cpg in names(fits)) {
    fit <- fits[[cpg]]
    if (is.null(fit) || inherits(
        fit, "dnaEPICO_methylationGLM_fit_error"
    ) || is.null(fit$coef)) {
        next
    }
    matched <- findCoefficientRowsMethylationGLM(
        coefNames = rownames(fit$coef), variable = phenotype,
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

#' Collect significant CpG coefficient tables from fitted one-timepoint GLMs
#'
#' @param modelResults Object returned by `fitMethylationGLMModels()`.
#' @param pvalThreshold Numeric. Threshold applied to omnibus p-values when
#'   omnibus testing was enabled during model fitting, or to phenotype
#'   main-effect or interaction coefficient p-values otherwise.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationGLM_significant_cpgs'`.
#'
#' @description
#' Collect coefficient tables for CpGs selected by the configured target test.
#' Selection uses target-term omnibus p-values when omnibus testing was enabled
#' during fitting, and phenotype main-effect or interaction coefficient
#' p-values otherwise.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' significant_cpgs <- collectSignificantCpGsMethylationGLM(
#'   modelResults = ex$modelResults,
#'   pvalThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(significant_cpgs)
#'
#' @export
collectSignificantCpGsMethylationGLM <- function(modelResults,
    pvalThreshold = 0.05, interactionTerm = NULL,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_methylationGLM.txt") {
    log_path <- resolveLogPathMinfiEwasWater(logs = logs,
        log_dir = log_dir, log_file = log_file)
    threshold <- validateProbabilityDnaEpico(pvalThreshold,
        "pvalThreshold")
    phenotypes <- modelResults$phenotypes
    if (is.null(phenotypes)) {
        phenotypes <- names(modelResults$fits)
    }
    use_omnibus <- isTRUE(modelResults$settings$omnibusTest)
    if (use_omnibus && !optionalTermMatchesMethylationModels(requested =
        interactionTerm,
        cached = modelResults$settings$interactionTerm)) {
        stop("interactionTerm does not match the term used for the fitted ",
            "omnibus tests.", call. = FALSE)
    }
    retained <- lapply(phenotypes, function(phenotype) {
        if (use_omnibus) {
            return(collectOmnibusGlmHitsDnaEpico(modelResults,
                phenotype, threshold))
        }
        cached <- collectCachedGlmHitsDnaEpico(modelResults,
            phenotype, threshold, interactionTerm)
        if (!is.null(cached)) {
            return(cached)
        }
        collectDirectGlmHitsDnaEpico(modelResults,
            phenotype, threshold, interactionTerm)
    })
    names(retained) <- phenotypes
    hit_counts <- vapply(retained, length,
        integer(1))
    emitLogMinfiEwasWater(c(
        "============================================================",
        paste("Significant CpGs retained at p <",
            threshold, ":"), paste(names(hit_counts),
            hit_counts, sep = ": ", collapse = "; "),
        "============================================================"),
        verbose = verbose, log_path = log_path)
    structure(retained, class = "dnaEPICO_methylationGLM_significant_cpgs")
}
glmDiagnosticFilesDnaEpico <- function(outputDir, fileKey, plots) {
    files <- list(
    qqplot = NULL, residualSD = NULL,
    residualSignificance = NULL, volcano = NULL,
    effectForest = NULL
    )
    if (is.null(outputDir)) {
    return(files)
    }
    files$qqplot <- file.path(
    outputDir, paste0("qqplot_", fileKey, "_coefficientPvalue.tiff")
    )
    if (!is.null(plots$residualSD)) {
    files$residualSD <- file.path(
        outputDir,
        paste0("residualSD_", fileKey, "_byAverageMethylation.tiff")
    )
    }
    if (!is.null(plots$residualSignificance)) {
    files$residualSignificance <- file.path(
        outputDir,
        paste0("residualSignificance_", fileKey, "_byPvalue.tiff")
    )
    }
    if (!is.null(plots$volcano)) {
    files$volcano <- file.path(
        outputDir,
        paste0("volcano_", fileKey, "_coefficientEstimate.tiff")
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

drawGlmDiagnosticPlotsDnaEpico <- function(plots, files, config) {
    for (plot_name in names(plots)) {
    plot <- plots[[plot_name]]
    if (is.null(plot)) {
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
        draw_fun = function() drawPlotObjectMinfiEwasWater(plot),
        display = config$display, file = files[[plot_name]],
        width = if (square) height else config$plotWidth,
        height = height, res = config$plotDPI
    )
    }
    invisible(NULL)
}

buildGlmTermDiagnosticDnaEpico <- function(
    summaryData, phenotype, term, termIndex, multipleTerms,
    diagnosticMean, config
) {
    diagnostic <- buildMethylationTermDiagnosticsDnaEpico(
    summaryData = summaryData, phenotype = phenotype, term = term,
    termColumn = "Coefficient", pValueColumn = "Pr(>|t|)",
    yColumn = "ResidualSD", yLabel = "Residual SD",
    diagnosticMean = diagnosticMean,
    fdrThreshold = config$fdrThreshold,
    estimateColumn = "Estimate", standardErrorColumn = "Std. Error"
    )
    if (is.null(diagnostic)) {
    return(NULL)
    }
    file_key <- phenotype
    if (multipleTerms) {
    file_key <- paste0(
        phenotype, "_", sprintf("%02d", termIndex), "_",
        sanitizeDiagnosticTermDnaEpico(term)
    )
    }
    files <- glmDiagnosticFilesDnaEpico(
    config$outputDir, file_key, diagnostic$plots
    )
    drawGlmDiagnosticPlotsDnaEpico(diagnostic$plots, files, config)
    list(plots = diagnostic$plots, files = files, lambda = diagnostic$lambda)
}

buildGlmPhenotypeDiagnosticsDnaEpico <- function(
    summaryData, phenotype, diagnosticMean, config
) {
    if (is.null(summaryData) || !nrow(summaryData)) {
    return(NULL)
    }
    summaryData$FDR <- adjustPvaluesByTermMethylationModels(
    pValues = summaryData[["Pr(>|t|)"]],
    terms = summaryData$Coefficient, method = config$padjmethod
    )
    terms <- unique(as.character(summaryData$Coefficient))
    terms <- terms[!is.na(terms)]
    multiple <- length(terms) > 1L
    results <- lapply(seq_along(terms), function(index) {
    buildGlmTermDiagnosticDnaEpico(
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

normalizeGlmDiagnosticConfigDnaEpico <- function(config) {
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

#' Plot diagnostic summaries for one-timepoint methylation GLMs
#'
#' @param modelSummaries Object returned by `summarizeMethylationGLMModels()`.
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
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
#' @return A list with class `'dnaEPICO_methylationGLM_diagnostic_plots'`
#'   containing the generated `ggplot2` objects, genomic inflation factors, and
#'   any saved TIFF file paths.
#'
#' @description
#' Create Q-Q and residual-diagnostic plots from the CpG summary tables returned
#' by `summarizeMethylationGLMModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' diagnostic_plots <- plotMethylationGLMDiagnostics(
#'   modelSummaries = ex$modelSummaries,
#'   preparedData = ex$preparedData,
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(diagnostic_plots$plots)
#'
#' @export
plotMethylationGLMDiagnostics <- function(
    modelSummaries, preparedData,
    fdrThreshold = 0.05, padjmethod = "fdr", outputDir = NULL,
    plotWidth = 2000L, plotHeight = 1000L, plotDPI = 150L,
    display = FALSE, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationGLM.txt"
) {
    config <- normalizeGlmDiagnosticConfigDnaEpico(
    as.list(environment(), all.names = TRUE)
    )
    summaries <- if (is.null(modelSummaries$diagnosticSummaries)) {
    modelSummaries$summaries
    } else {
    modelSummaries$diagnosticSummaries
    }
    diagnostic_mean <- diagnosticMeanMethylationModels(preparedData)
    results <- lapply(names(summaries), function(phenotype) {
    buildGlmPhenotypeDiagnosticsDnaEpico(
        summaries[[phenotype]], phenotype, diagnostic_mean, config
    )
    })
    names(results) <- names(summaries)
    results <- Filter(Negate(is.null), results)
    plots <- lapply(results, function(result) result[["plots"]])
    files <- lapply(results, function(result) result[["files"]])
    inflation <- lapply(results, function(result) result[["lambda"]])
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
    ), class = "dnaEPICO_methylationGLM_diagnostic_plots")
}

buildAnnotationModelMessagesDnaEpico <- function(modelMessages) {
    if (!is.data.frame(modelMessages) || nrow(modelMessages) == 0L) {
    return(data.frame())
    }
    required <- c("Phenotype", "CpG", "Model.Message")
    if (!all(required %in% names(modelMessages))) {
    return(data.frame())
    }
    phenotype_tables <- lapply(
    unique(as.character(modelMessages$Phenotype)),
    function(phenotype) {
        table <- unique(modelMessages[
        modelMessages$Phenotype == phenotype,
        c("CpG", "Model.Message"),
        drop = FALSE
        ])
        colnames(table)[[2L]] <- paste0(phenotype, "_Model.Message")
        table
    }
    )
    if (length(phenotype_tables) == 1L) {
    phenotype_tables[[1L]]
    } else {
    Reduce(
        function(x, y) merge(x, y, by = "CpG", all = TRUE),
        phenotype_tables
    )
    }
}

buildAnnotationOmnibusTablesMethylationGLM <- function(modelSummaries) {
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

duplicatedGlmCoefficientsDnaEpico <- function(summaryList) {
    occurrences <- unlist(lapply(summaryList, function(data) {
    if (is.null(data) || !nrow(data) ||
        !"Coefficient" %in% colnames(data)) {
        return(character(0))
    }
    unique(as.character(data$Coefficient))
    }), use.names = FALSE)
    unique(occurrences[duplicated(occurrences)])
}

cleanGlmModelSummaryDnaEpico <- function(data, phenotype, duplicates) {
    if (is.null(data) || !nrow(data)) {
    return(NULL)
    }
    coefficient_names <- unique(data$Coefficient)
    tables <- lapply(coefficient_names, function(coefficient) {
    output <- data[data$Coefficient == coefficient,
        c("CpG", "Pr(>|t|)"),
        drop = FALSE
    ]
    clean_name <- gsub("`", "", coefficient, fixed = TRUE)
    colnames(output)[[2L]] <- if (coefficient %in% duplicates) {
        paste0(phenotype, "_", clean_name, "_P.Value")
    } else {
        paste0(clean_name, "P.Value")
    }
    output
    })
    if (!length(tables)) {
    return(NULL)
    }
    Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), tables)
}

mergeGlmSummaryTablesDnaEpico <- function(summaryList, modelSummaries) {
    duplicates <- duplicatedGlmCoefficientsDnaEpico(summaryList)
    cleaned <- lapply(names(summaryList), function(phenotype) {
    cleanGlmModelSummaryDnaEpico(
        summaryList[[phenotype]], phenotype, duplicates
    )
    })
    cleaned <- Filter(Negate(is.null), cleaned)
    merged <- if (!length(cleaned)) {
    data.frame(CpG = character(0))
    } else {
    Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), cleaned)
    }
    omnibus <- buildAnnotationOmnibusTablesMethylationGLM(modelSummaries)
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

appendGlmGencodeDnaEpico <- function(data, annotationObject, enabled) {
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

annotateGlmSummaryDataDnaEpico <- function(
    merged, annotationData, annotationCols, annotationObject, gencodeHub
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
    gencode <- appendGlmGencodeDnaEpico(
    annotated, annotationObject, gencodeHub
    )
    used <- c(available, gencode$columns)
    data <- orderAnnotatedModelColumnsDnaEpico(
    data = gencode$data, annotationCols = used
    )
    list(
    data = data, used = used, missing = missing,
    gencode = gencode$result
    )
}

logGlmAnnotationDnaEpico <- function(result, gencodeHub, verbose, logPath) {
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

#' Annotate one-timepoint GLM summary tables with array annotation metadata
#'
#' @param modelSummaries Object returned by `summarizeMethylationGLMModels()`
#'   or a named list of CpG summary data frames.
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
#' @return A list with class `'dnaEPICO_methylationGLM_annotation'`
#'   containing the annotated summary table and any requested annotation columns
#'   that were unavailable in the chosen annotation object.
#'
#' @description
#' Merge phenotype-specific CpG summary tables with probe annotation metadata
#' and
#' return a single annotated result table.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' annotation_data <- annotateMethylationGLMSummaries(
#'   modelSummaries = ex$modelSummaries,
#'   annotationObject = ex$annotationData,
#'   annotationCols = "Name,chr,pos",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(annotation_data)
#'
#' @export
annotateMethylationGLMSummaries <- function(
    modelSummaries, annotationObject,
    annotationCols = c(
    "Name", "chr", "pos", "UCSC_RefGene_Group",
    "UCSC_RefGene_Name", "Relation_to_Island", "GencodeV41_Group"
    ),
    gencodeHub = FALSE, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationGLM.txt"
) {
    gencodeHub <- validateLogicalScalarDnaEpico(gencodeHub, "gencodeHub")
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    summaries <- modelSummaries
    if (!is.null(modelSummaries$diagnosticSummaries)) {
    summaries <- modelSummaries$diagnosticSummaries
    } else if (!is.null(modelSummaries$summaries)) {
    summaries <- modelSummaries$summaries
    }
    merged <- mergeGlmSummaryTablesDnaEpico(summaries, modelSummaries)
    annotation_cols <- splitOptionMinfiEwasWater(annotationCols, sep = ",")
    result <- annotateGlmSummaryDataDnaEpico(
    merged, coerceAnnotationDataMethylationGLM(annotationObject),
    annotation_cols, annotationObject, gencodeHub
    )
    logGlmAnnotationDnaEpico(result, gencodeHub, verbose, log_path)
    structure(list(
    data = result$data, fitFailures = modelSummaries$fitFailures,
    modelMessages = modelSummaries$modelMessages,
    annotationColumnsUsed = result$used,
    missingAnnotationCols = result$missing,
    gencodeHub = result$gencode
    ), class = "dnaEPICO_methylationGLM_annotation")
}

validateGlmReportAssetsDirDnaEpico <- function(reportAssetsDir) {
    if (!is.null(reportAssetsDir) &&
    (!is.character(reportAssetsDir) || length(reportAssetsDir) != 1L ||
        is.na(reportAssetsDir) || !nzchar(reportAssetsDir))) {
    stop("reportAssetsDir must be NULL or one non-empty directory path.",
        call. = FALSE
    )
    }
    invisible(TRUE)
}

saveGlmPhenotypeSummariesDnaEpico <- function(
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
    file <- file.path(outputDir, paste0(phenotype, "SummaryGLM.rds"))
    artifact <- modelResults$phenotypeSummaries[[phenotype]]
    if (is.null(artifact)) {
        stop(
        "A complete compact phenotype summary is unavailable for ",
        phenotype, ".",
        call. = FALSE
        )
    }
    existing <- loadPhenotypeSummaryMethylationModels(
        file, artifact$signature
    )
    if (is.null(existing$object)) {
        savePhenotypeSummaryMethylationModels(artifact, file)
    }
    files[[phenotype]] <- file
    }
    files
}

writeGlmTextSummariesDnaEpico <- function(summaries, outputDir, enabled) {
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
    if ("Pr(>|t|)" %in% colnames(data)) {
        data <- data[order(data[["Pr(>|t|)"]]), , drop = FALSE]
    }
    file <- file.path(outputDir, paste0(phenotype, "SummaryGLM.txt"))
    utils::write.table(
        data,
        file = file, sep = "\t", row.names = FALSE, quote = FALSE
    )
    files[[phenotype]] <- file
    }
    files
}

writeSignificantGlmTablesDnaEpico <- function(
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

addGlmGencodeWorkbookRowsDnaEpico <- function(
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

addGlmVennWorkbookRowsDnaEpico <- function(
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

writeGlmAnnotatedWorkbookDnaEpico <- function(
    modelResults, modelSummaries, annotatedResults, annotatedGLMOut,
    reportAssetsDir, vennDResults
) {
    annotated <- if (is.null(annotatedResults$data)) {
    annotatedResults
    } else {
    annotatedResults$data
    }
    report_table <- sortReportTableDnaEpico(annotated)
    annotated <- report_table$data
    file <- file.path(annotatedGLMOut, "annotatedGLM.xlsx")
    dictionary <- buildAnnotatedWorkbookDictionaryMethylationGLM(
    columns = colnames(annotated),
    modelDescription = "Pvalue from GLM model",
    formulaText = modelResults$formulas, modelLabel = "GLM",
    responseLabel = inferMethylationValueLabelMethylationGLM(modelResults)
    )
    metadata <- buildModelWorkbookMetadataDnaEpico(
    modelResults, modelSummaries, annotatedResults,
    analysis = "glm"
    )
    rows <- addGlmGencodeWorkbookRowsDnaEpico(
    dictionary, metadata, annotatedResults
    )
    rows <- addGlmVennWorkbookRowsDnaEpico(
    rows$dictionary, rows$metadata, vennDResults
    )
    writeAnnotatedWorkbookMethylationGLM(
    annotated_df = annotated, file = file, resultSheet = "annotatedGLM",
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
        sidecarDir = reportAssetsDir, sheet = "annotatedGLM",
        idColumn = report_table$idColumn,
        dictionary = rows$dictionary,
        workbookMetadata = rows$metadata
    )
    }
    list(file = file, sidecar = sidecar, sheets = names(rows$sheets))
}

logGlmOutputFilesDnaEpico <- function(
    summaryFiles, textFiles, significantFiles, workbook,
    fitFailures, saveTxtSummaries, saveSignificantCpGs,
    verbose, logPath
) {
    sidecar_line <- if (is.null(workbook$sidecar$table)) {
    "Report table sidecar:          not requested"
    } else {
    paste("Report table sidecar:         ", workbook$sidecar$table)
    }
    text_count <- if (isTRUE(saveTxtSummaries)) length(textFiles) else 0L
    significant_count <- if (isTRUE(saveSignificantCpGs)) {
    sum(vapply(significantFiles, length, integer(1)))
    } else {
    0L
    }
    emitLogMinfiEwasWater(c(
    "============================================================",
    "Full model files written:       0",
    paste("Compact phenotype summaries: ", length(summaryFiles)),
    paste("Annotated results file:      ", workbook$file), sidecar_line,
    paste("Hard model errors retained:  ", nrow(fitFailures)),
    paste("Summary text files written:   ", text_count),
    paste("Significant CpG text files:  ", significant_count),
    "============================================================"
    ), verbose = verbose, log_path = logPath)
}

#' Write optional disk outputs for one-timepoint GLM analyses
#'
#' @param modelResults Object returned by `fitMethylationGLMModels()`.
#' @param modelSummaries Object returned by `summarizeMethylationGLMModels()`.
#' @param annotatedResults Object returned by
#'   `annotateMethylationGLMSummaries()` or a compatible data frame.
#' @param significantCpGs Object returned by
#'   `collectSignificantCpGsMethylationGLM()` or `NULL`.
#' @param outputRData Character. Directory used for complete compact phenotype
#'   summaries.
#' @param summaryTxtDir Character. Directory used for tab-delimited summary
#'   tables.
#' @param significantCpGDir Character. Directory used for significant-CpG
#'   coefficient tables.
#' @param annotatedGLMOut Character. Directory used for the annotated summary
#'   XLSX workbook.
#' @param reportAssetsDir Character or `NULL`. Report results directory used for
#'   the compressed TSV table and compact metadata sidecars.
#' @param vennDResults Optional model-level Venn result. Its configuration
#'   metadata and threshold tables are added before the workbook dictionary.
#' @param saveTxtSummaries Logical. If `TRUE`, write tab-delimited summary
#'   tables.
#' @param saveSignificantCpGs Logical. If `TRUE`, write significant-CpG
#'   coefficient tables.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationGLM_paths'` containing
#'   the paths of the files written to disk, including the annotated workbook
#'   and any requested report-table sidecars.
#'
#' @description
#' Write compact phenotype summaries, optional text and significant-CpG tables,
#' and annotated results from the one-timepoint GLM workflow.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' annotation_data <- annotateMethylationGLMSummaries(
#'   modelSummaries = ex$modelSummaries,
#'   annotationObject = ex$annotationData,
#'   annotationCols = "Name,chr,pos",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' significant_cpgs <- collectSignificantCpGsMethylationGLM(
#'   modelResults = ex$modelResults,
#'   pvalThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' output_paths <- writeMethylationGLMOutputs(
#'   modelResults = ex$modelResults,
#'   modelSummaries = ex$modelSummaries,
#'   annotatedResults = annotation_data,
#'   significantCpGs = significant_cpgs,
#'   outputRData = file.path(ex$tempDir, "models"),
#'   summaryTxtDir = file.path(ex$tempDir, "summary"),
#'   significantCpGDir = file.path(ex$tempDir, "significant"),
#'   annotatedGLMOut = file.path(ex$tempDir, "annotated"),
#'   saveTxtSummaries = TRUE,
#'   saveSignificantCpGs = TRUE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(output_paths)
#'
#' @export
writeMethylationGLMOutputs <- function(
    modelResults, modelSummaries, annotatedResults,
    significantCpGs = NULL, outputRData, summaryTxtDir,
    significantCpGDir, annotatedGLMOut, reportAssetsDir = NULL,
    vennDResults = NULL, saveTxtSummaries = TRUE,
    saveSignificantCpGs = FALSE, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationGLM.txt"
) {
    validateGlmReportAssetsDirDnaEpico(reportAssetsDir)
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    dir.create(outputRData, recursive = TRUE, showWarnings = FALSE)
    dir.create(annotatedGLMOut, recursive = TRUE, showWarnings = FALSE)
    summary_files <- saveGlmPhenotypeSummariesDnaEpico(
    modelResults, modelSummaries, outputRData
    )
    text_files <- writeGlmTextSummariesDnaEpico(
    modelSummaries$summaries, summaryTxtDir, saveTxtSummaries
    )
    significant_files <- writeSignificantGlmTablesDnaEpico(
    significantCpGs, significantCpGDir, saveSignificantCpGs
    )
    workbook <- writeGlmAnnotatedWorkbookDnaEpico(
    modelResults, modelSummaries, annotatedResults, annotatedGLMOut,
    reportAssetsDir, vennDResults
    )
    fit_failures <- modelSummaries$fitFailures
    if (is.null(fit_failures)) {
    fit_failures <- modelResults$fitFailures
    }
    logGlmOutputFilesDnaEpico(
    summary_files, text_files, significant_files, workbook,
    fit_failures, saveTxtSummaries, saveSignificantCpGs,
    verbose, log_path
    )
    structure(list(
    modelFiles = character(0), summaryFiles = summary_files,
    summaryTxtFiles = text_files,
    significantCpGFiles = significant_files,
    annotatedGLM = workbook$file,
    annotatedGLMText = workbook$sidecar$table,
    annotatedGLMReportMetadata = workbook$sidecar$metadata,
    annotatedGLMDictionary = workbook$sidecar$dictionary,
    annotatedGLMMetadata = workbook$sidecar$workbookMetadata,
    vennDSheets = workbook$sheets
    ), class = "dnaEPICO_methylationGLM_paths")
}
