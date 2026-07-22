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
        stop("padjmethod must be one of: ", paste(stats::p.adjust.methods,
            collapse = ", "
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
                CpG = cpg, Status = if (is.null(fit_object$status)) {
                    "failed"
                } else {
                    as.character(fit_object$status)
                }, Error = as.character(fit_object$error),
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
    invalid_cpgs <- preparedData$invalidCpGs$CpG
    if (length(invalid_cpgs) > 0L) {
        means[names(means) %in% invalid_cpgs] <- NA_real_
    }
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

buildMethylationTermDiagnosticsDnaEpico <- function(
    summaryData,
    phenotype, term, termColumn, pValueColumn, yColumn, yLabel,
    diagnosticMean, fdrThreshold
) {
    term_values <- as.character(summaryData[[termColumn]])
    term_data <- summaryData[!is.na(term_values) & term_values ==
        term, , drop = FALSE]
    valid_p <-
        is.finite(term_data[[pValueColumn]]) & term_data[[pValueColumn]] >=
        0 & term_data[[pValueColumn]] <= 1
    term_data <- term_data[valid_p, , drop = FALSE]
    if (nrow(term_data) == 0L) {
        return(NULL)
    }

    p_values <- pmax(term_data[[pValueColumn]], .Machine$double.xmin)
    chi_square <- stats::qchisq(p_values, df = 1, lower.tail = FALSE)
    lambda <- round(stats::median(chi_square) / stats::qchisq(0.5,
        df = 1
    ), 3)
    diagnostic_label <- paste0(phenotype, " [", term, "]")
    qq_data <- data.frame(
        expected = -log10(stats::ppoints(length(p_values))),
        observed = -log10(sort(p_values))
    )
    qq_plot <- ggplot2::ggplot(qq_data, ggplot2::aes(
        x = expected,
        y = observed
    )) +
        ggplot2::geom_point(color = "black") +
        ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
        ggplot2::labs(title = paste0(
            "Q-Q Plot of p-values for ",
            diagnostic_label, "\nGenomic Inflation Factor = ",
            lambda
        ), x = "Expected -log10(p)", y = "Observed -log10(p)") +
        ggplot2::theme_minimal()

    term_data$meanMethylation <- diagnosticMean$values[term_data$CpG]
    term_data$diagnosticP <- p_values
    residual_plot <- NULL
    residual_significance_plot <- NULL
    if (yColumn %in% colnames(term_data)) {
        term_data$diagnosticY <- term_data[[yColumn]]
        residual_plot <- ggplot2::ggplot(term_data, ggplot2::aes(
            x = meanMethylation,
            y = diagnosticY
        )) +
            ggplot2::geom_point(
                alpha = 0.6,
                color = "black"
            ) +
            ggplot2::labs(
                title = paste(
                    yLabel,
                    "vs", diagnosticMean$label, "for", diagnostic_label
                ),
                x = diagnosticMean$label, y = yLabel
            ) +
            ggplot2::theme_minimal()

        complete_diagnostic <- is.finite(term_data$meanMethylation) &
            is.finite(term_data$diagnosticY)
        if (sum(complete_diagnostic) >= 3L) {
            residual_plot <- residual_plot + ggplot2::geom_smooth(
                method = "loess",
                formula = y ~ x, se = FALSE, color = "red"
            )
        }

        significant_data <- term_data[!is.na(term_data$FDR) &
            term_data$FDR < fdrThreshold, , drop = FALSE]
        residual_significance_plot <- ggplot2::ggplot(
            term_data,
            ggplot2::aes(
                x = -log10(diagnosticP), y = diagnosticY,
                color = FDR < fdrThreshold
            )
        ) +
            ggplot2::geom_point(alpha = 0.6) +
            ggrepel::geom_text_repel(
                data = significant_data,
                ggplot2::aes(label = CpG)
            ) +
            ggplot2::scale_color_manual(values = c(
                `FALSE` = "grey50",
                `TRUE` = "firebrick"
            )) +
            ggplot2::labs(
                title = paste(
                    yLabel,
                    "vs Significance for", diagnostic_label
                ), x = "-log10(p-value)",
                y = yLabel, color = paste("FDR <", fdrThreshold)
            ) +
            ggplot2::theme_minimal()
    }

    list(lambda = lambda, plots = list(
        qqplot = qq_plot, residualSD = residual_plot,
        residualSignificance = residual_significance_plot
    ))
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
        stop("Each prsMap entry must follow the format 'Phenotype:PRS'. Invalid entries: ",
            paste(map_entries[invalid], collapse = ", "),
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
        stop("prsMap contains duplicate phenotype mappings: ",
            paste(duplicated_keys, collapse = ", "),
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
        stop("The fixed-effect design contains non-finite numeric values. Check phenotypes, covariates, and interactions.",
            call. = FALSE
        )
    }
    design_rank <- qr(design)$rank
    if (design_rank < ncol(design)) {
        stop("The fixed-effect design matrix is rank deficient (rank ",
            design_rank, " of ", ncol(design),
                "). Check duplicated covariates, factor levels, and interactions.",
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

fitStatusValuesMethylationGLM <- function(modelObj) {
    if (inherits(modelObj, "dnaEPICO_methylationGLM_fit_error")) {
        converged <- if (is.null(modelObj$converged)) {
            NA
        } else {
            as.logical(modelObj$converged)[[1L]]
        }
        convergence_message <- if (is.null(modelObj$convergenceMessage)) {
            NA_character_
        } else {
            as.character(modelObj$convergenceMessage)[[1L]]
        }
        return(list(
            status = if (is.null(modelObj$status)) "failed" else as.character(modelObj$status),
            singular = NA, converged = converged,
                convergenceMessage = convergence_message,
            warning = NA_character_
        ))
    }

    fit_status <- modelObj$fitStatus
    converged <- if (is.null(fit_status$converged)) {
        TRUE
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
    } else if (length(fit_warnings) > 0L) {
        "fitted_with_warning"
    } else {
        "fitted"
    }, singular = NA, converged = converged,
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

collectFitDiagnosticsMethylationGLM <- function(fits) {
    rows <- list()
    row_index <- 1L
    for (phenotype in names(fits)) {
        fit_group <- fits[[phenotype]]
        for (cpg in names(fit_group)) {
            fit_object <- fit_group[[cpg]]
            values <- fitStatusValuesMethylationGLM(fit_object)
            inference_included <- startsWith(values$status, "fitted")
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
                } else if (identical(values$converged, FALSE)) {
                    "model did not converge"
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

applyFitQualityExclusionsMethylationGLM <- function(
    summaryDf,
    excludeNonConverged = FALSE
) {
    summary_df <- summaryDf
    if (is.null(summary_df) || nrow(summary_df) == 0L) {
        return(summary_df)
    }
    if (!("Inference.Included" %in% names(summary_df))) {
        summary_df$Inference.Included <- TRUE
    }
    if (!("Exclusion.Reason" %in% names(summary_df))) {
        summary_df$Exclusion.Reason <- NA_character_
    }

    excluded <- rep(FALSE, nrow(summary_df))
    if (isTRUE(excludeNonConverged) && "Converged" %in% names(summary_df)) {
        converged <- as.logical(summary_df[["Converged"]])
        excluded <- !is.na(converged) & !converged
    }
    summary_df$Inference.Included[excluded] <- FALSE
    summary_df$Exclusion.Reason[excluded] <- "model did not converge"
    if (any(excluded) && "Pr(>|t|)" %in% names(summary_df)) {
        summary_df[["Pr(>|t|)"]][excluded] <- NA_real_
    }
    summary_df
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
    if (is.null(modelObj) || inherits(modelObj,
        "dnaEPICO_methylationGLM_fit_error")) {
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
    fit_status <- fitStatusValuesMethylationGLM(modelObj)
    summary_df$Fit.Status <- fit_status$status
    summary_df$Converged <- fit_status$converged
    summary_df$Convergence.Message <- fit_status$convergenceMessage
    summary_df$Fit.Warning <- fit_status$warning
    summary_df$Inference.Included <- !identical(
        fit_status$converged,
        FALSE
    )
    summary_df$Exclusion.Reason <- if (identical(
        fit_status$converged,
        FALSE
    )) {
        "model did not converge"
    } else {
        NA_character_
    }

    if (isTRUE(includeResidualSD) && !is.null(modelObj$residuals)) {
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

    max(1L, as.integer(available))
}

resolveParallelPlanMethylationModels <- function(
    engine, nCores,
    nCpGs
) {
    requested_cores <- validatePositiveIntegerMethylationModels(
        nCores,
        "nCores"
    )
    n_cpgs <- max(0L, as.integer(nCpGs))
    requested_backend <- tolower(Sys.getenv(
        "DNAEPICO_PARALLEL_BACKEND",
        "auto"
    ))
    if (!(requested_backend %in% c("auto", "fork", "psock", "serial"))) {
        requested_backend <- "auto"
    }
    crossover <- parallelCrossoverMethylationModels(engine)
    resource_cap <- availableWorkersMethylationModels()
    worker_count <- min(requested_cores, resource_cap, max(
        1L,
        n_cpgs
    ))

    if (requested_cores <= 1L || worker_count <= 1L) {
        backend <- "serial"
        reason <- "one effective worker"
    } else if (identical(requested_backend, "serial")) {
        backend <- "serial"
        worker_count <- 1L
        reason <- "serial backend requested"
    } else if (identical(requested_backend, "auto") && n_cpgs <
        crossover) {
        backend <- "serial"
        worker_count <- 1L
        reason <- paste0(
            "workload below the ", engine, " crossover of ",
            crossover, " CpGs"
        )
    } else {
        backend <- resolveParallelBackendMethylationModels(worker_count)
        reason <- if (identical(requested_backend, "auto")) {
            paste0("workload reached the ", engine, " crossover")
        } else {
            paste0(requested_backend, " backend requested")
        }
    }

    list(
        engine = tolower(as.character(engine[[1L]])), backend = backend,
        requestedCores = requested_cores, workerCount = worker_count,
        resourceWorkerCap = resource_cap, workloadCpGs = n_cpgs,
        crossoverCpGs = crossover, reason = reason, forced = !identical(
            requested_backend,
            "auto"
        )
    )
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
        "Fit.Status", "Converged", "Convergence.Message", "Fit.Warning",
        "Inference.Included", "Exclusion.Reason"
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

fitCpGModelMethylationGLM <- function(
    cpg, cpgValues, modelData,
    formulaText, responseVar = "beta"
) {
    tryCatch(
        {
            model_data <- modelData
            model_data[[responseVar]] <- as.numeric(cpgValues)
            observed_response <-
                model_data[[responseVar]][!is.na(model_data[[responseVar]])]
            if (length(unique(observed_response)) < 2L) {
                stop("The CpG response has no observed variation.",
                    call. = FALSE
                )
            }

            warning_state <- new.env(parent = emptyenv())
            warning_state$messages <- character(0)
            fit <- withCallingHandlers(
                glm2::glm2(
                    formula = stats::as.formula(formulaText),
                    data = model_data, family = stats::gaussian(),
                        na.action = stats::na.exclude
                ),
                warning = function(condition) {
                    warning_state$messages <- c(
                        warning_state$messages,
                        conditionMessage(condition)
                    )
                    invokeRestart("muffleWarning")
                }
            )
            fit_warnings <- warning_state$messages
            if (!isTRUE(fit$converged)) {
                stop("The GLM did not converge.", call. = FALSE)
            }
            if (fit$rank < length(stats::coef(fit)) ||
                anyNA(stats::coef(fit))) {
                stop("The CpG-specific design became rank deficient after excluding missing values.",
                    call. = FALSE
                )
            }
            if (fit$df.residual <= 0L) {
                stop("The GLM has no residual degrees of freedom.",
                    call. = FALSE
                )
            }
            coef_table <- summary(fit)$coefficients
            required_statistics <- c(
                "Estimate", "Std. Error", "t value",
                "Pr(>|t|)"
            )
            if (!all(required_statistics %in% colnames(coef_table)) ||
                any(!is.finite(as.matrix(coef_table[, required_statistics,
                    drop = FALSE
                ])))) {
                stop("The GLM returned missing or non-finite coefficient statistics.",
                    call. = FALSE
                )
            }

            coefficient_terms <- buildCoefficientTermMapMethylationModels(
                formulaText = formulaText,
                data = model_data
            )
            residuals <- stats::residuals(fit)
            residual_sd <- if (fit$df.residual > 0L) {
                sqrt(sum(residuals^2, na.rm = TRUE) / fit$df.residual)
            } else {
                NA_real_
            }

            list(
                coef = coef_table, residuals = residuals,
                    fitted = stats::fitted(fit),
                residualSD = residual_sd, coefficientTerms = coefficient_terms,
                fitStatus = list(
                    converged = TRUE,
                        warnings = unique(as.character(fit_warnings)),
                    convergenceMessages = character(0)
                )
            )
        },
        error = function(error) {
            reason <- conditionMessage(error)
            convergence_failure <- grepl("did not converge|convergence",
                reason,
                ignore.case = TRUE
            )
            newMethylationFitErrorDnaEpico(
                reason = reason,
                    errorClass = "dnaEPICO_methylationGLM_fit_error",
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

fitMethylationGLMBatch <- function(
    cpgBatch, data, modelData,
    formulaText, phenotype, interactionTerm = NULL, responseVar = "beta",
    invalidCpgReasons = character(0)
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
                errorClass = "dnaEPICO_methylationGLM_fit_error",
                status = "invalid"
            )
        } else {
            fitCpGModelMethylationGLM(
                cpg = cpg, cpgValues = data[[cpg]],
                modelData = modelData, formulaText = formulaText,
                responseVar = responseVar
            )
        }
        fits[[cpg]] <- model_obj
        summaries[[cpg]] <- summarizeCpGFitMethylationGLM(
            cpg = cpg,
            modelObj = model_obj, variable = phenotype,
                interactionTerm = interactionTerm,
            includeResidualSD = TRUE
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
createDistributionPlotMethylationGLM <- function(
    values, variable,
    type = c("hist", "bar"), fill = "steelblue"
) {
    type <- match.arg(type)

    if (identical(type, "hist")) {
        plot_data <- data.frame(value = as.numeric(values))
        return(ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
            ggplot2::geom_histogram(bins = 30, fill = fill, color = "white") +
            ggplot2::labs(
                title = paste("Distribution of", variable),
                x = variable, y = "Frequency"
            ) +
            ggplot2::theme_minimal())
    }

    plot_data <- data.frame(value = as.factor(values))
    ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
        ggplot2::geom_bar(fill = fill) +
        ggplot2::labs(
            title = paste("Distribution of", variable),
            x = variable, y = "Count"
        ) +
        ggplot2::theme_minimal()
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
        annotation_lookup <-
            suppressPackageStartupMessages(tryCatch(minfi::getAnnotation(annotationObject),
            error = function(e) NULL
        ))
        if (!is.null(annotation_lookup)) {
            return(annotationObject)
        }

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

coerceAnnotationDataMethylationGLM <- function(annotationObject) {
    if (is.data.frame(annotationObject)) {
        annotation_df <- annotationObject
        if (!("CpG" %in% colnames(annotation_df))) {
            if ("IlmnID" %in% colnames(annotation_df)) {
                annotation_df$CpG <- annotation_df$IlmnID
            } else if (.row_names_info(annotation_df, type = 1L) <
                0L) {
                stop("annotationObject data frames must include a CpG or IlmnID column, or explicit probe row names.",
                    call. = FALSE
                )
            } else {
                annotation_df$CpG <- rownames(annotation_df)
            }
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

buildAnnotatedWorkbookDictionaryMethylationGLM <- function(
    columns,
    modelDescription, formulaText, modelLabel, responseLabel
) {
    pvalue_columns <- grepl("P\\.Value$|P\\.value$", columns)
    omnibus_columns <- grepl("_Omnibus_", columns, fixed = TRUE)
    default_formula <- paste(responseLabel, "~ formula unavailable")
    formula_values <- stats::setNames(
        as.character(formulaText),
        names(formulaText)
    )
    formula_values <- formula_values[!is.na(formula_values) &
        nzchar(formula_values)]

    select_formula <- function(column) {
        if (length(formula_values) == 0L) {
            return(default_formula)
        }

        formula_names <- names(formula_values)
        if (!is.null(formula_names) && any(nzchar(formula_names))) {
            for (formula_name in formula_names[nzchar(formula_names)]) {
                if (grepl(
                    paste0("^", escapeRegexMethylationGLM(formula_name)),
                    column
                )) {
                    return(unname(formula_values[[formula_name]]))
                }
            }
        }

        if (length(formula_values) == 1L) {
            return(unname(formula_values[[1L]]))
        }

        paste(unique(unname(formula_values)), collapse = "; ")
    }

    display_formula <- function(column) {
        formula_text <- select_formula(column)
        formula_parts <- trimws(unlist(strsplit(formula_text,
            "\\s*;\\s*",
            perl = TRUE
        ), use.names = FALSE))
        formula_parts <- vapply(formula_parts, function(formula_part) {
            sub("^\\s*[^~]+\\s*~", paste(responseLabel, "~"),
                formula_part,
                ignore.case = TRUE
            )
        }, character(1))
        paste0(modelLabel, ": ", paste(formula_parts, collapse = "; "))
    }

    descriptions <- vapply(seq_along(columns), function(index) {
        column <- columns[[index]]
        if (grepl("_Omnibus_F\\.Value$", column)) {
            return("F statistic for the joint fixed-effect omnibus test")
        }
        if (grepl("_Omnibus_Num\\.DF$", column)) {
            return("Numerator degrees of freedom equal to the estimable rank of the omnibus contrast")
        }
        if (grepl("_Omnibus_Den\\.DF$", column)) {
            return("Denominator degrees of freedom for the omnibus F test")
        }
        if (grepl("_Omnibus_Adjusted\\.P\\.Value$", column)) {
            return("Omnibus p-value adjusted across valid CpGs within the phenotype and tested term")
        }
        if (grepl("_Omnibus_P\\.Value$", column)) {
            return("Raw p-value for the joint null hypothesis that all estimable coefficients in the tested fixed-effect term equal zero")
        }
        if (grepl("_Omnibus_Method$", column)) {
            return("Denominator degrees-of-freedom method used for the omnibus F test")
        }
        if (pvalue_columns[[index]]) {
            return(modelDescription)
        }
        if (column %in% c("IlmnID", "CpG", "Name")) {
            return("CpG probe identifier")
        }
        if (identical(column, "Fit.Status")) {
            return("Overall CpG model-fitting status across the requested phenotypes")
        }
        if (grepl("_Fit\\.Status$", column)) {
            return("Phenotype-specific CpG model-fitting status")
        }
        if (grepl("_Singular\\.Fit$", column)) {
            return("Whether lme4 detected a singular random-effects fit; NA when not assessed")
        }
        if (identical(column, "Singular.Fit")) {
            return("Whether any requested lme4 model had a singular random-effects fit; NA when not assessed")
        }
        if (identical(column, "Converged")) {
            return("Whether all assessed CpG models converged")
        }
        if (grepl("_Converged$", column)) {
            return("Whether the phenotype-specific CpG model converged")
        }
        if (grepl("_Convergence\\.Message$", column)) {
            return("Phenotype-specific convergence message retained for interpretation")
        }
        if (grepl("_Fit\\.Warning$", column)) {
            return("Phenotype-specific model-fitting warning retained for interpretation")
        }
        if (identical(column, "Inference.Included")) {
            return("Whether all requested phenotype-specific fits were included in inference")
        }
        if (grepl("_Inference\\.Included$", column)) {
            return("Whether the phenotype-specific CpG fit was included in inference")
        }
        if (identical(column, "Fit.Quality.Reason")) {
            return("Reason a CpG fit was excluded from inference while retaining its row")
        }
        if (grepl("_Exclusion\\.Reason$", column)) {
            return("Phenotype-specific reason the CpG fit was excluded from inference")
        }
        if (identical(column, "Failure.Phenotypes")) {
            return("Phenotypes for which the CpG model was invalid or failed")
        }
        if (identical(column, "Failure.Reason")) {
            return("Reported reason for an invalid or failed CpG model")
        }
        "Genomic annotation or supporting result column"
    }, character(1))
    formulas <- rep("", length(columns))
    formula_columns <- pvalue_columns | omnibus_columns
    formulas[formula_columns] <- vapply(
        columns[formula_columns],
        display_formula, character(1)
    )

    data.frame(
        Column = columns, Description = descriptions,
        Formula = formulas, stringsAsFactors = FALSE, check.names = FALSE
    )
}

orderAnnotatedModelColumnsDnaEpico <- function(
    data, annotationCols = character(0),
    includeSingular = FALSE, modelNames = NULL
) {
    columns <- colnames(data)
    id_columns <- intersect(c("IlmnID", "CpG"), columns)
    omnibus_columns <- columns[grepl("_Omnibus_", columns, fixed = TRUE)]
    omnibus_pvalue_columns <- omnibus_columns[endsWith(
        omnibus_columns,
        "_Omnibus_P.Value"
    )]
    omnibus_detail_columns <- setdiff(
        omnibus_columns,
        omnibus_pvalue_columns
    )
    pvalue_columns <- setdiff(columns[grepl(
        "P\\.Value$|P\\.value$",
        columns
    )], omnibus_columns)
    annotation_columns <- intersect(annotationCols, columns)
    diagnostic_suffixes <- c(
        "Fit.Status", if (isTRUE(includeSingular)) "Singular.Fit",
        "Converged", "Convergence.Message", "Fit.Warning", "Inference.Included",
        "Exclusion.Reason"
    )
    diagnostic_pattern <- paste0("_(?:", paste(gsub(
        "\\.", "\\\\.",
        diagnostic_suffixes
    ), collapse = "|"), ")$")
    diagnostic_columns <- columns[grepl(diagnostic_pattern, columns,
        perl = TRUE
    )]
    if (!isTRUE(includeSingular)) {
        diagnostic_columns <- diagnostic_columns[!grepl(
            "_Singular\\.Fit$",
            diagnostic_columns
        )]
    }
    supporting_columns <- setdiff(columns, c(
        id_columns, pvalue_columns,
        omnibus_columns, annotation_columns, diagnostic_columns
    ))
    omnibus_suffix_order <- c(
        "_Omnibus_F.Value", "_Omnibus_Num.DF", "_Omnibus_Den.DF",
        "_Omnibus_P.Value", "_Omnibus_Adjusted.P.Value",
        "_Omnibus_Method"
    )
    order_omnibus <- function(values) {
        ranks <- vapply(values, function(value) {
            matches <- which(endsWith(value, omnibus_suffix_order))
            if (length(matches) == 0L) {
                length(omnibus_suffix_order) + 1L
            } else {
                matches[[1L]]
            }
        }, integer(1))
        values[order(ranks, match(values, columns))]
    }
    result_columns <- c(pvalue_columns, omnibus_pvalue_columns)
    model_names <- as.character(modelNames)
    model_names <- model_names[!is.na(model_names) & nzchar(model_names)]
    if (length(model_names) > 0L && length(result_columns) > 0L) {
        owners <- vapply(result_columns, function(column) {
            candidates <- model_names[startsWith(
                column,
                paste0(model_names, "_")
            )]
            if (length(candidates) == 0L) {
                return(NA_character_)
            }
            candidates[[which.max(nchar(candidates))]]
        }, character(1))
        grouped_results <- unlist(lapply(model_names, function(model_name) {
            model_columns <- result_columns[owners == model_name]
            c(
                intersect(pvalue_columns, model_columns),
                intersect(omnibus_pvalue_columns, model_columns)
            )
        }), use.names = FALSE)
        result_columns <- unique(c(
            grouped_results,
            result_columns[is.na(owners)]
        ))
    }
    ordered_omnibus_details <- order_omnibus(omnibus_detail_columns)
    if (length(model_names) > 0L &&
        length(ordered_omnibus_details) > 0L) {
        detail_owners <- vapply(ordered_omnibus_details, function(column) {
            candidates <- model_names[startsWith(
                column,
                paste0(model_names, "_")
            )]
            if (length(candidates) == 0L) {
                return(NA_character_)
            }
            candidates[[which.max(nchar(candidates))]]
        }, character(1))
        grouped_details <- unlist(lapply(model_names, function(model_name) {
            order_omnibus(ordered_omnibus_details[
                detail_owners == model_name
            ])
        }), use.names = FALSE)
        ordered_omnibus_details <- unique(c(
            grouped_details,
            ordered_omnibus_details[is.na(detail_owners)]
        ))
    }
    ordered_columns <- unique(c(
        id_columns, result_columns, annotation_columns,
        ordered_omnibus_details, supporting_columns,
        diagnostic_columns
    ))
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

buildModelWorkbookMetadataDnaEpico <- function(
    modelResults,
    modelSummaries, annotatedResults, analysis = c("glm", "lme")
) {
    analysis <- match.arg(analysis)
    settings <- modelResults$settings
    diagnostics <- modelResults$fitDiagnostics
    if (!is.data.frame(diagnostics)) {
        diagnostics <- data.frame()
    }
    fit_status <- if ("Fit.Status" %in% colnames(diagnostics)) {
        as.character(diagnostics$Fit.Status)
    } else {
        character(0)
    }
    fitted <- startsWith(fit_status, "fitted")
    converged <- if ("Converged" %in% colnames(diagnostics)) {
        as.logical(diagnostics$Converged)
    } else {
        logical(0)
    }
    included <- if ("Inference.Included" %in% colnames(diagnostics)) {
        as.logical(diagnostics$Inference.Included)
    } else {
        logical(0)
    }
    singular <- if ("Singular.Fit" %in% colnames(diagnostics)) {
        as.logical(diagnostics$Singular.Fit)
    } else {
        logical(0)
    }
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

    keys <- c(
        "analysis", "backend", "fitting_function", "engine_version",
        "formulas", "response_label", "methylation_scale", "sample_count",
        "phenotypes", "covariates", "factor_variables", "scaled_variables",
        "interaction_term", "cpg_count", "fit_count", "fitted_count",
        "failed_or_invalid_count", "non_converged_count", "singular_count",
        "inference_excluded_count", "parallel_backend", "worker_count",
        "annotation_columns", "missing_annotation_columns", "created"
    )
    fitting_function <- if (identical(analysis, "glm")) {
        "glm2::glm2"
    } else if (identical(lme_engine, "nlme")) {
        "nlme::lme"
    } else {
        "lmerTest::lmer"
    }
    values <- c(
        if (identical(analysis, "glm")) "methylationGLM" else "methylationLME",
        lme_engine, fitting_function, engine_version,
            paste(paste(names(modelResults$formulas),
            modelResults$formulas,
            sep = ": "
        ), collapse = " | "),
        inferMethylationValueLabelMethylationGLM(modelResults),
        modelSettingTextDnaEpico(settings$methylationScale),
        modelSettingTextDnaEpico(settings$sampleCount),
            modelSettingTextDnaEpico(modelResults$phenotypes),
        modelSettingTextDnaEpico(settings$covariates),
            modelSettingTextDnaEpico(settings$factorVars),
        modelSettingTextDnaEpico(settings$scaleVars),
            modelSettingTextDnaEpico(settings$interactionTerm),
        as.character(if (!is.null(annotatedResults$data)) nrow(annotatedResults$data) else NA_integer_),
        as.character(nrow(diagnostics)), as.character(sum(fitted,
            na.rm = TRUE
        )), as.character(sum(!fitted, na.rm = TRUE)),
        as.character(sum(!is.na(converged) & !converged)),
            as.character(sum(!is.na(singular) &
            singular)), as.character(sum(!is.na(included) & !included)),
        modelSettingTextDnaEpico(settings$parallelBackend),
            modelSettingTextDnaEpico(settings$workerCount),
        modelSettingTextDnaEpico(annotation_used),
            modelSettingTextDnaEpico(annotation_missing),
        format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
    metadata <- data.frame(
        Key = keys, Value = values, stringsAsFactors = FALSE,
        check.names = FALSE
    )

    factor_levels <- settings$factorLevels
    if (is.list(factor_levels) && length(factor_levels) > 0L) {
        factor_rows <- data.frame(Key = paste0(
            "factor.", names(factor_levels),
            ".levels"
        ), Value = vapply(
            factor_levels, function(levels) modelSettingTextDnaEpico(levels),
            character(1)
        ), stringsAsFactors = FALSE, check.names = FALSE)
        metadata <- rbind(metadata, factor_rows)
    }
    if (identical(analysis, "lme")) {
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
        omnibus_tested <- if ("Omnibus.Status" %in%
            colnames(omnibus_data)) {
            sum(omnibus_data$Omnibus.Status == "tested", na.rm = TRUE)
        } else {
            0L
        }
        omnibus_unavailable <- if ("Omnibus.Status" %in%
            colnames(omnibus_data)) {
            sum(omnibus_data$Omnibus.Status != "tested", na.rm = TRUE)
        } else {
            0L
        }
        lmer_test_version <- tryCatch(
            as.character(utils::packageVersion("lmerTest")),
            error = function(error) "unavailable"
        )
        pbkrtest_version <- tryCatch(
            as.character(utils::packageVersion("pbkrtest")),
            error = function(error) "unavailable"
        )
        lme_rows <- data.frame(
            Key = c(
                "libraries", "person_variable",
                "time_variable", "correlation_structure",
                    "correlation_variable",
                "random_effect_structure", "exclude_singular",
                "exclude_non_converged", "omnibus_test",
                "omnibus_target", "omnibus_ddf", "omnibus_rhs",
                "omnibus_joint", "omnibus_eps", "omnibus_tested_count",
                "omnibus_unavailable_count", "p_adjust_method",
                "p_adjust_scope", "lmerTest_version", "pbkrtest_version"
            ), Value = c(
                modelSettingTextDnaEpico(settings$lmeLibs),
                modelSettingTextDnaEpico(settings$personVar),
                    modelSettingTextDnaEpico(settings$timeVar),
                modelSettingTextDnaEpico(settings$correlationStructure,
                    empty = "none"
                ), modelSettingTextDnaEpico(settings$correlationVar),
                paste0(
                    "(1 | ",
                    modelSettingTextDnaEpico(settings$personVar), ")"
                ),
                modelSettingTextDnaEpico(modelSummaries$settings$excludeSingular),
                modelSettingTextDnaEpico(modelSummaries$settings$excludeNonConverged),
                modelSettingTextDnaEpico(settings$omnibusTest),
                modelSettingTextDnaEpico(modelResults$omnibusTargets),
                modelSettingTextDnaEpico(settings$omnibusDdf),
                modelSettingTextDnaEpico(settings$omnibusRhs),
                modelSettingTextDnaEpico(settings$omnibusJoint),
                modelSettingTextDnaEpico(settings$omnibusEps),
                as.character(omnibus_tested),
                as.character(omnibus_unavailable),
                modelSettingTextDnaEpico(modelSummaries$settings$padjmethod),
                "within each phenotype and omnibus term across valid CpGs",
                lmer_test_version, pbkrtest_version
            ),
            stringsAsFactors = FALSE, check.names = FALSE
        )
        metadata <- rbind(metadata, lme_rows)
    } else {
        glm_rows <- data.frame(
            Key = c("libraries", "exclude_non_converged"),
            Value = c(
                modelSettingTextDnaEpico(settings$glmLibs),
                modelSettingTextDnaEpico(modelSummaries$settings$excludeNonConverged)
            ),
            stringsAsFactors = FALSE, check.names = FALSE
        )
        metadata <- rbind(metadata, glm_rows)
    }

    scaling <- settings$scalingMetadata
    if (is.data.frame(scaling) && nrow(scaling) > 0L) {
        scaling_rows <- do.call(rbind, lapply(
            seq_len(nrow(scaling)),
            function(index) {
                variable <- scaling$Variable[[index]]
                data.frame(
                    Key = c(
                        paste0(
                            "scale.", variable,
                            ".center"
                        ), paste0("scale.", variable, ".sd"),
                        paste0("scale.", variable, ".finite_values"),
                        paste0("scale.", variable, ".missing_values")
                    ),
                    Value = as.character(c(
                        scaling$Center[[index]],
                        scaling$Scale[[index]], scaling$Finite.Values[[index]],
                        scaling$Missing.Values[[index]]
                    )), stringsAsFactors = FALSE,
                    check.names = FALSE
                )
            }
        ))
        metadata <- rbind(metadata, scaling_rows)
    }
    rownames(metadata) <- NULL
    metadata
}

writeAnnotatedWorkbookMethylationGLM <- function(
    annotated_df,
    file, resultSheet, dictionary, metadata = NULL
) {
    workbook <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(workbook, resultSheet)
    openxlsx::writeData(workbook, sheet = resultSheet, x = annotated_df)
    if (is.data.frame(metadata) && nrow(metadata) > 0L) {
        openxlsx::addWorksheet(workbook, "metadata")
        openxlsx::writeData(workbook, sheet = "metadata", x = metadata)
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

writeReportTableSidecarDnaEpico <- function(
    tableData, workbookFile, sidecarDir,
    sheet, idColumn, dictionary, workbookMetadata = NULL
) {
    paths <- reportTableSidecarPathsDnaEpico(sidecarDir, sheet)
    dir.create(dirname(paths$table), recursive = TRUE, showWarnings = FALSE)
    if (!(idColumn %in% names(tableData))) {
        stop("The report table identifier column was not found.",
            call. = FALSE
        )
    }
    row_order <- order(as.character(tableData[[idColumn]]), na.last = TRUE)
    if (!identical(row_order, seq_len(nrow(tableData)))) {
        stop("The report table must be sorted by its identifier column.",
            call. = FALSE
        )
    }
    data.table::fwrite(tableData,
        file = paths$table, sep = "\t",
        quote = TRUE, na = "", compress = "gzip", showProgress = FALSE
    )

    sidecar_metadata <- data.frame(Key = c(
        "format_version",
        "sheet", "workbook_md5", "rows", "columns", "id_column",
        "sorted_by_id"
    ), Value = c(
        "1", sheet, unname(tools::md5sum(workbookFile)),
        as.character(nrow(tableData)), as.character(ncol(tableData)),
        idColumn, "TRUE"
    ), stringsAsFactors = FALSE, check.names = FALSE)
    data.table::fwrite(dictionary,
        file = paths$dictionary, sep = "\t",
        quote = TRUE, na = ""
    )
    if (is.data.frame(workbookMetadata) && nrow(workbookMetadata) >
        0L) {
        data.table::fwrite(workbookMetadata,
            file = paths$workbookMetadata,
            sep = "\t", quote = TRUE, na = ""
        )
    } else if (file.exists(paths$workbookMetadata)) {
        unlink(paths$workbookMetadata, force = TRUE)
    }
    data.table::fwrite(sidecar_metadata,
        file = paths$metadata,
        sep = "\t", quote = TRUE, na = ""
    )

    paths
}

resolveReportTableSidecarDnaEpico <- function(workbookFile, sidecarDir, sheet) {
    paths <- reportTableSidecarPathsDnaEpico(sidecarDir, sheet)
    unavailable <- function(reason) {
        list(
            ok = FALSE, reason = reason, table = paths$table,
            metadata = paths$metadata
        )
    }
    if (!file.exists(paths$table) || !file.exists(paths$metadata) ||
        !file.exists(paths$dictionary)) {
        return(unavailable("A report table sidecar was not found."))
    }

    metadata <- tryCatch(
        utils::read.delim(paths$metadata,
            check.names = FALSE,
            stringsAsFactors = FALSE, quote = "\"", comment.char = ""
        ),
        error = function(error) error
    )
    if (inherits(metadata, "error") || !all(c("Key", "Value") %in%
        names(metadata))) {
        return(unavailable("The report table sidecar metadata is invalid."))
    }
    values <- as.character(metadata$Value)
    names(values) <- as.character(metadata$Key)
    required <- c(
        "format_version", "sheet", "workbook_md5",
        "rows", "columns", "id_column", "sorted_by_id"
    )
    if (!all(required %in% names(values))) {
        return(unavailable("The report table sidecar metadata is incomplete."))
    }
    if (!identical(values[["format_version"]], "1") || !identical(
        values[["sheet"]],
        sheet
    ) || !identical(values[["sorted_by_id"]], "TRUE")) {
        return(unavailable("The report table sidecar metadata is incompatible."))
    }
    workbook_md5 <- unname(tools::md5sum(workbookFile))
    if (!identical(values[["workbook_md5"]], workbook_md5)) {
        return(unavailable("The report table sidecar does not match the workbook."))
    }

    list(
        ok = TRUE, reason = NULL, table = paths$table,
            metadata = paths$metadata,
        dictionary = paths$dictionary,
            workbookMetadata = paths$workbookMetadata,
        rows = as.integer(values[["rows"]]),
            columns = as.integer(values[["columns"]]),
        idColumn = values[["id_column"]]
    )
}

streamReportTableDnaEpico <- function(
    tableFile, chunkSize = 5000L,
    expectedRows = NULL, expectedColumns = NULL, chunkHandler = NULL
) {
    chunk_size <- max(1L, as.integer(chunkSize))
    connection <- if (grepl("\\.gz$", tableFile, ignore.case = TRUE)) {
        gzfile(tableFile, open = "rt")
    } else {
        file(tableFile, open = "rt")
    }
    on.exit(close(connection), add = TRUE)

    header_line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(header_line) == 0L) {
        stop("The report table sidecar is empty.", call. = FALSE)
    }
    columns <- names(utils::read.delim(
        text = header_line, nrows = 0L,
        check.names = FALSE, quote = "\"", comment.char = ""
    ))
    if (!is.null(expectedColumns) && !identical(
        length(columns),
        as.integer(expectedColumns)
    )) {
        stop("The report table sidecar columns do not match its metadata.",
            call. = FALSE
        )
    }

    n_rows <- 0L
    chunk_number <- 1L
    maximum_chunk_rows <- 0L
    repeat {
        chunk <- utils::read.delim(connection,
            header = FALSE,
            col.names = columns, nrows = chunk_size, sep = "\t",
            quote = "\"", comment.char = "", colClasses = "character",
            check.names = FALSE, na.strings = ""
        )
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

    if (!is.null(expectedRows) && !identical(
        as.integer(n_rows),
        as.integer(expectedRows)
    )) {
        stop("The report table sidecar row count does not match its metadata.",
            call. = FALSE
        )
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
        return(methylationScaleResponseLabelDnaEpico(modelResults$settings$methylationScale))
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
            if (!is.list(fit_object) || inherits(
                fit_object,
                "dnaEPICO_methylationGLM_fit_error"
            ) || inherits(
                fit_object,
                "dnaEPICO_methylationLME_fit_error"
            ) || is.null(fit_object$fitted) ||
                is.null(fit_object$residuals)) {
                next
            }

            fitted_values <- fit_object$fitted
            residual_values <- fit_object$residuals
            if (!is.numeric(fitted_values) || !is.numeric(residual_values) ||
                length(fitted_values) != length(residual_values)) {
                next
            }

            response_values <- fitted_values + residual_values
            response_values <- response_values[is.finite(response_values)]
            if (length(response_values) == 0L) {
                next
            }

            found_response <- TRUE
            if (any(response_values < 0 | response_values > 1,
                na.rm = TRUE
            )) {
                return("M-values")
            }
        }
    }

    if (isTRUE(found_response)) {
        return("Beta values")
    }

    "Beta values"
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
#'     inputPheno = ex$inputPath,
#'     phenotypes = "status",
#'     covariates = "sex,age",
#'     factorVars = "status,sex",
#'     cpgLimit = 2,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(prepared_data)
#'
#' @export
prepareMethylationGLMData <- function(
    inputPheno, phenotypes,
    covariates, factorVars, scaleVars = NULL, cpgPrefix = "cg",
    cpgLimit = NA, methylationScale = "beta", interactionTerm = NULL,
    prsMap = NULL, verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_methylationGLM.txt"
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
        "phenoBT1"
    } else {
        character(0)
    }
    analysis_data <- loadSavedObjectPreprocessingPheno(inputPheno,
        preferred_name = c(scale_named_input, paste0(
            methylation_prefix,
            "T1"
        ), legacy_input_name, input_object_name)
    )

    if (!is.data.frame(analysis_data)) {
        analysis_data <- as.data.frame(analysis_data, stringsAsFactors = FALSE)
    }

    if (length(phenotype_list) == 0L) {
        stop("At least one phenotype must be supplied.", call. = FALSE)
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

    for (var in intersect(factor_list, colnames(analysis_data))) {
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
        "CpG columns in the GLM input"
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
        protectedVars = cpg_columns
    )

    requested_columns <- unique(c(phenotype_list, covariate_list))
    missing_counts <- vapply(
        requested_columns,
            function(column_name) sum(is.na(analysis_data[[column_name]])),
        integer(1)
    )
    variable_summary <- summary(analysis_data[, requested_columns,
        drop = FALSE
    ])
    interaction_table <- NULL
    if (!is.null(resolved_interaction) && nzchar(resolved_interaction)) {
        interaction_table <- table(analysis_data[[resolved_interaction]],
            useNA = "ifany"
        )
    }

    log_lines <- c(
        "============================================================",
        paste("Loaded phenotype + methylation data from:", inputPheno),
        paste(
            "Merged modeling object:     ", methylation_prefix,
            "*"
        ), paste("Data dimensions:             ", paste(dim(analysis_data),
            collapse = " x "
        )), paste(
            "Phenotypes:                  ",
            paste(phenotype_list, collapse = ", ")
        ), paste(
            "Covariates:                  ",
            paste(covariate_list, collapse = ", ")
        ), paste(
            "Factor variables:            ",
            paste(factor_list, collapse = ", ")
        ), formatScalingMetadataLogDnaEpico(scaling$metadata),
        paste("CpG columns retained:        ", length(cpg_columns)),
        formatMethylationBoundariesLogDnaEpico(methylation_validation$boundaries),
        if (nrow(methylation_validation$issues) > 0L) {
            c("Methylation value issues:",
                previewLinesMinfiEwasWater(methylation_validation$issues))
        } else {
            "Methylation value issues:     none"
        }, "Missing summary:", paste(names(missing_counts), missing_counts,
            sep = ": ", collapse = "; "
        ), "Summary statistics:",
        previewLinesMinfiEwasWater(variable_summary)
    )
    if (!is.null(interaction_table)) {
        log_lines <- c(log_lines, paste(
            "Interaction table for",
            resolved_interaction, ":"
        ), previewLinesMinfiEwasWater(interaction_table))
    }
    log_lines <- c(log_lines,
        "============================================================")
    emitLogMinfiEwasWater(log_lines, verbose = verbose, log_path = log_path)

    structure(list(
        data = analysis_data, modelData = scaling$data,
        phenotypes = phenotype_list, covariates = covariate_list,
        factorVars = factor_list, scaleVars = scaling$scaleVars,
        scalingMetadata = scaling$metadata, cpgColumns = cpg_columns,
        cpgPrefix = cpgPrefix, cpgLimit = cpg_limit,
            methylationScale = methylation_scale,
        responseLabel = methylation_label,
            methylationObjectPrefix = methylation_prefix,
        internalResponseColumn = response_column, prsMap = prs_map,
        interactionTerm = resolved_interaction,
            requestedInteractionTerm = interactionTerm,
        methylationBoundaries = methylation_validation$boundaries,
        methylationIssues = methylation_validation$issues,
            invalidCpGs = methylation_validation$invalidCpGs,
        missingCounts = missing_counts, variableSummary = variable_summary,
        interactionTable = interaction_table
    ), class = "dnaEPICO_methylationGLM_data")
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
#'     preparedData = ex$preparedData,
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(distribution_plots)
#'
#' @export
plotMethylationGLMDistributions <- function(
    preparedData, plotWidth = 2000L,
    plotHeight = 1000L, plotDPI = 150L, outputDir = NULL, display = FALSE,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_methylationGLM.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    analysis_data <- preparedData$data
    phenotype_plots <- list()
    factor_plots <- list()
    covariate_plots <- list()
    saved_files <- list(
        phenotypes = list(), factors = list(),
        covariates = list()
    )

    for (var in preparedData$phenotypes) {
        if (!var %in% colnames(analysis_data)) {
            next
        }

        plot_type <- if (is.numeric(analysis_data[[var]])) {
            "hist"
        } else {
            "bar"
        }
        plot_object <- createDistributionPlotMethylationGLM(
            values = analysis_data[[var]],
            variable = var, type = plot_type, fill = "steelblue"
        )
        phenotype_plots[[var]] <- plot_object

        file_path <- NULL
        if (!is.null(outputDir)) {
            file_path <- file.path(outputDir, paste0(if (identical(
                plot_type,
                "hist"
            )) {
                "hist_"
            } else {
                "bar_"
            }, var, ".tiff"))
        }

        runPlotMinfiEwasWater(
            draw_fun = function() drawPlotObjectMinfiEwasWater(plot_object),
            display = display, file = file_path, width = plotWidth,
            height = plotHeight, res = plotDPI
        )

        if (!is.null(file_path)) {
            saved_files$phenotypes[[var]] <- file_path
        }
    }

    factor_vars <- intersect(preparedData$factorVars, colnames(analysis_data))
    for (var in factor_vars) {
        plot_object <- createDistributionPlotMethylationGLM(
            values = analysis_data[[var]],
            variable = var, type = "bar", fill = "darkorange"
        )
        factor_plots[[var]] <- plot_object

        file_path <- NULL
        if (!is.null(outputDir)) {
            file_path <- file.path(outputDir, paste0(
                "bar_",
                var, ".tiff"
            ))
        }

        runPlotMinfiEwasWater(
            draw_fun = function() drawPlotObjectMinfiEwasWater(plot_object),
            display = display, file = file_path, width = plotWidth,
            height = plotHeight, res = plotDPI
        )

        if (!is.null(file_path)) {
            saved_files$factors[[var]] <- file_path
        }
    }

    numeric_covariates <- setdiff(preparedData$covariates,
        preparedData$factorVars)
    numeric_covariates <- intersect(numeric_covariates, colnames(analysis_data))
    for (var in numeric_covariates) {
        plot_object <- createDistributionPlotMethylationGLM(
            values = analysis_data[[var]],
            variable = var, type = "hist", fill = "darkgreen"
        )
        covariate_plots[[var]] <- plot_object

        file_path <- NULL
        if (!is.null(outputDir)) {
            file_path <- file.path(outputDir, paste0(
                "hist_",
                var, ".tiff"
            ))
        }

        runPlotMinfiEwasWater(
            draw_fun = function() drawPlotObjectMinfiEwasWater(plot_object),
            display = display, file = file_path, width = plotWidth,
            height = plotHeight, res = plotDPI
        )

        if (!is.null(file_path)) {
            saved_files$covariates[[var]] <- file_path
        }
    }

    emitLogMinfiEwasWater(
        c(
            "============================================================",
            paste("Phenotype distribution plots: ", length(phenotype_plots)),
            paste("Factor distribution plots:    ", length(factor_plots)),
            paste("Numeric covariate plots:      ", length(covariate_plots)),
            if (is.null(outputDir)) {
                "Distribution plots were returned in memory only."
            } else {
                paste("Distribution plots saved to:  ", outputDir)
            }, "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(
        phenotypes = phenotype_plots, factors = factor_plots,
        covariates = covariate_plots, files = saved_files
    ), class = "dnaEPICO_methylationGLM_distribution_plots")
}

#' Fit CpG-wise Gaussian GLMs for one-timepoint methylation analyses
#'
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
#' @param nCores Integer. Maximum number of worker processes to use. Automatic
#'   fitting remains serial below the empirical glm2 crossover and caps workers
#'   by available CpGs and detected physical cores.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded
#'   to worker processes.
#' @param glmLibs Character vector or comma-separated string of package names to
#'   check on worker processes. The default is `'glm2'`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationGLM_models'` containing
#'   fitted model lists, model formulas, and counts of failed CpG fits.
#'
#' @description
#' Fit one Gaussian GLM per CpG for each phenotype requested in the object
#' returned by `prepareMethylationGLMData()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' model_results <- fitMethylationGLMModels(
#'     preparedData = ex$preparedData,
#'     nCores = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(model_results$fits)
#'
#' @export
fitMethylationGLMModels <- function(
    preparedData, nCores = 1L,
    libPath = NULL, glmLibs = "glm2", verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationGLM.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    if (is.null(libPath)) {
        libPath <- .libPaths()
    }

    glm_lib_list <- splitOptionMinfiEwasWater(glmLibs, sep = ",")
    if (length(glm_lib_list) == 0L) {
        glm_lib_list <- "glm2"
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
    formulas <- stats::setNames(
        character(length(preparedData$phenotypes)),
        preparedData$phenotypes
    )
    failure_counts <- stats::setNames(
        integer(length(preparedData$phenotypes)),
        preparedData$phenotypes
    )
    failure_reasons <- list()
    parallel_plan <- resolveParallelPlanMethylationModels(
        engine = "glm2",
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
                "invalid_cpg_reasons", "libPath", "glm_lib_list",
                "validateWorkerPackagesMethylationModels",
                    "newMethylationFitErrorDnaEpico",
                "fitMethylationGLMBatch", "fitCpGModelMethylationGLM",
                "buildCoefficientTermMapMethylationModels",
                    "removeRandomInterceptMethylationModels",
                "summarizeCpGFitMethylationGLM",
                    "fitStatusValuesMethylationGLM",
                "findCoefficientRowsMethylationGLM", "escapeRegexMethylationGLM"
            ),
            envir = environment()
        )
        parallel::clusterEvalQ(psock_cluster,
            validateWorkerPackagesMethylationModels(
            libPath = libPath,
            packages = glm_lib_list
        ))
    }

    for (phenotype in preparedData$phenotypes) {
        prs_var <- character(0)
        if (phenotype %in% names(preparedData$prsMap)) {
            prs_var <- unname(preparedData$prsMap[[phenotype]])
        }
        covariates <- unique(c(preparedData$covariates, prs_var))
        model_vars <- unique(c(phenotype, covariates,
            preparedData$interactionTerm))
        model_vars <- model_vars[!is.na(model_vars) & nzchar(model_vars)]

        missing_vars <- setdiff(model_vars, colnames(model_data))
        if (length(missing_vars) > 0L) {
            stop("Model variables not found for phenotype ",
                phenotype, ": ", paste(missing_vars, collapse = ", "),
                call. = FALSE
            )
        }

        formula_text <- buildFormulaMethylationGLM(
            phenotype = phenotype,
            covariates = covariates,
                interactionTerm = preparedData$interactionTerm,
            responseVar = preparedData$internalResponseColumn
        )

        base_model_data <- model_data[, model_vars, drop = FALSE]
        factor_vars <- preparedData$factorVars
        for (var in intersect(factor_vars, colnames(base_model_data))) {
            base_model_data[[var]] <- as.factor(base_model_data[[var]])
        }
        validateFixedEffectDesignMethylationModels(
            formulaText = formula_text,
            data = base_model_data
        )
        batch_worker <- fitMethylationGLMBatch
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
                            packages = glm_lib_list
                        )
                        batch_worker(
                            cpgBatch = batch, data = analysis_data,
                            modelData = base_model_data,
                                formulaText = formula_text,
                            phenotype = phenotype,
                                interactionTerm = resolved_interaction,
                            responseVar = response_var,
                                invalidCpgReasons = invalid_cpg_reasons
                        )
                    },
                    mc.cores = cluster_size, mc.preschedule = FALSE
                )
            } else {
                parallel::clusterExport(psock_cluster, varlist = c(
                    "base_model_data",
                    "formula_text", "phenotype", "resolved_interaction",
                    "response_var", "batch_worker"
                ), envir = environment())
                batch_results <- parallel::parLapplyLB(
                    psock_cluster,
                    cpg_batches, function(batch) {
                        batch_worker(
                            cpgBatch = batch, data = analysis_data,
                            modelData = base_model_data,
                                formulaText = formula_text,
                            phenotype = phenotype,
                                interactionTerm = resolved_interaction,
                            responseVar = response_var,
                                invalidCpgReasons = invalid_cpg_reasons
                        )
                    }
                )
            }
        } else {
            validateWorkerPackagesMethylationModels(
                libPath = libPath,
                packages = glm_lib_list
            )
            batch_results <- lapply(cpg_batches, function(batch) {
                batch_worker(
                    cpgBatch = batch, data = analysis_data,
                    modelData = base_model_data, formulaText = formula_text,
                    phenotype = phenotype,
                        interactionTerm = resolved_interaction,
                    responseVar = response_var,
                        invalidCpgReasons = invalid_cpg_reasons
                )
            })
        }

        combined_results <- combineFitBatchResultsMethylationModels(
            batchResults = batch_results,
            cpgColumns = cpg_columns
        )
        fit_list <- combined_results$fits
        phenotype_summary_cache <- filterSummaryByPvalueMethylationGLM(
            summaryDf = combined_results$summaries,
            pValueFilter = NA_real_, includeResidualSD = TRUE
        )
        failures <- vapply(fit_list, function(x) {
            inherits(
                x,
                "dnaEPICO_methylationGLM_fit_error"
            )
        }, logical(1))
        error_counts <- summarizeFitErrorsMethylationModels(
            fitList = fit_list,
            errorClass = "dnaEPICO_methylationGLM_fit_error"
        )
        fits[[phenotype]] <- fit_list
        summary_cache[[phenotype]] <- phenotype_summary_cache
        formulas[[phenotype]] <- formula_text
        failure_counts[[phenotype]] <- sum(failures)
        failure_reasons[[phenotype]] <- error_counts

        emitLogMinfiEwasWater(
            c(
                "============================================================",
                paste("Fitted phenotype:            ", phenotype),
                paste("Formula:                     ", formula_text),
                paste("CpGs attempted:              ", length(cpg_columns)),
                paste("Failed CpG fits:             ",
                    failure_counts[[phenotype]]),
                paste("Top fit errors:               ",
                    formatFitErrorsMethylationModels(error_counts)),
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
            warning("All CpG GLM fits failed for phenotype '",
                phenotype, "'. Top failure reasons: ",
                    formatFitErrorsMethylationModels(error_counts),
                ". The failure inventory was retained and the analysis continued.",
                call. = FALSE
            )
        }
    }

    fit_failures <- collectFitFailuresMethylationModels(
        fits = fits,
        errorClass = "dnaEPICO_methylationGLM_fit_error"
    )
    fit_diagnostics <- collectFitDiagnosticsMethylationGLM(fits)

    if (!is.null(psock_cluster)) {
        parallel::stopCluster(psock_cluster)
        psock_cluster <- NULL
    }

    structure(
        list(
            fits = fits, summaryCache = summary_cache,
            formulas = formulas, phenotypes = names(fits),
                failureCounts = failure_counts,
            failureReasons = failure_reasons, fitFailures = fit_failures,
            fitDiagnostics = fit_diagnostics,
                methylationBoundaries = preparedData$methylationBoundaries,
            methylationIssues = preparedData$methylationIssues,
                invalidCpGs = preparedData$invalidCpGs,
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
                glmLibs = glm_lib_list,
                    methylationScale = preparedData$methylationScale,
                methylationObjectPrefix = preparedData$methylationObjectPrefix,
                responseLabel = preparedData$responseLabel,
                    internalResponseColumn = preparedData$internalResponseColumn,
                interactionTerm = preparedData$interactionTerm,
                    phenotypes = preparedData$phenotypes,
                covariates = preparedData$covariates,
                    factorVars = preparedData$factorVars,
                factorLevels = lapply(stats::setNames(
                    preparedData$factorVars,
                    preparedData$factorVars
                ), function(variable) levels(preparedData$data[[variable]])),
                scaleVars = preparedData$scaleVars,
                    scalingMetadata = preparedData$scalingMetadata,
                sampleCount = nrow(preparedData$data)
            ), responseLabel = preparedData$responseLabel
        ),
        class = "dnaEPICO_methylationGLM_models"
    )
}

#' Summarize CpG-wise Gaussian GLM fits for one-timepoint analyses
#'
#' @param modelResults Object returned by `fitMethylationGLMModels()`.
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
#' @param summaryResidualSD Logical. If `TRUE`, add residual standard deviations
#'   to each CpG summary row.
#' @param summaryPval Numeric or `NA`. Optional p-value filter applied to the
#'   returned summary tables. `NA` keeps all rows.
#' @param excludeNonConverged Logical. If `TRUE`, retain non-converged CpG rows
#'   and their reasons but set their inferential p-values to `NA`.
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
#'   does not remove CpGs from those outputs.
#'
#' @description
#' Extract phenotype-specific CpG coefficient tables from the fitted model
#' object returned by `fitMethylationGLMModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' summary_results <- summarizeMethylationGLMModels(
#'     modelResults = ex$modelResults,
#'     preparedData = ex$preparedData,
#'     summaryResidualSD = TRUE,
#'     summaryPval = NA,
#'     nCores = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(summary_results$summaries)
#'
#' @export
summarizeMethylationGLMModels <- function(
    modelResults, preparedData,
    summaryResidualSD = TRUE, summaryPval = NA, excludeNonConverged = FALSE,
    nCores = 1L, libPath = NULL, glmLibs = "glm2", chunkSize = NULL,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_methylationGLM.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    if (is.null(libPath)) {
        libPath <- .libPaths()
    }

    glm_lib_list <- splitOptionMinfiEwasWater(glmLibs, sep = ",")
    if (length(glm_lib_list) == 0L) {
        glm_lib_list <- "glm2"
    }

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
    summaries <- list()
    diagnostic_summaries <- list()

    for (phenotype in names(modelResults$fits)) {
        if (!is.null(modelResults$summaryCache) &&
            !is.null(modelResults$summaryCache[[phenotype]])) {
            diagnostic_df <- filterSummaryByPvalueMethylationGLM(
                summaryDf = modelResults$summaryCache[[phenotype]],
                pValueFilter = NA_real_,
                    includeResidualSD = isTRUE(summaryResidualSD)
            )
            diagnostic_df <- applyFitQualityExclusionsMethylationGLM(
                summaryDf = diagnostic_df,
                excludeNonConverged = excludeNonConverged
            )
            summary_df <- filterSummaryByPvalueMethylationGLM(
                summaryDf = diagnostic_df,
                pValueFilter = p_value_filter,
                    includeResidualSD = isTRUE(summaryResidualSD)
            )
            diagnostic_summaries[[phenotype]] <- diagnostic_df
            summaries[[phenotype]] <- summary_df

            emitLogMinfiEwasWater(
                c(
                    "============================================================",
                    paste("Summarized phenotype:        ", phenotype),
                    paste("CpG summary rows returned:   ", nrow(summary_df)),
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

        summary_worker <- summarizeCpGFitMethylationGLM
        resolved_interaction <- preparedData$interactionTerm
        include_residual_sd <- isTRUE(summaryResidualSD)

        if (n_cores > 1L && length(cpg_chunks) > 1L) {
            cluster_size <- min(n_cores, length(cpg_chunks))
            cl <- parallel::makeCluster(cluster_size)
            on.exit(parallel::stopCluster(cl), add = TRUE)

            parallel::clusterExport(cl,
                varlist = c(
                    "fit_list",
                    "phenotype", "resolved_interaction", "include_residual_sd",
                    "summary_worker", "libPath", "glm_lib_list"
                ),
                envir = environment()
            )

            parallel::clusterEvalQ(cl, {
                if (!is.null(libPath)) {
                    .libPaths(unique(c(libPath, .libPaths())))
                }

                for (pkg in glm_lib_list) {
                    if (!requireNamespace(pkg, quietly = TRUE)) {
                        stop("Failed to load package: ", pkg, call. = FALSE)
                    }
                }

                NULL
            })

            result_chunks <- parallel::parLapplyLB(
                cl, cpg_chunks,
                function(chunk) {
                    rows <- lapply(chunk, function(cpg) {
                        summary_worker(
                            cpg = cpg, modelObj = fit_list[[cpg]],
                            variable = phenotype,
                                interactionTerm = resolved_interaction,
                            includeResidualSD = include_residual_sd
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
                        variable = phenotype,
                            interactionTerm = resolved_interaction,
                        includeResidualSD = include_residual_sd
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
            summary_df <- summary_df[, c(
                "CpG", "Coefficient",
                "Estimate", "Std. Error", "t value", "Pr(>|t|)",
                if (isTRUE(summaryResidualSD)) "ResidualSD",
                "Fit.Status", "Singular.Fit", "Converged",
                    "Convergence.Message",
                "Fit.Warning", "Inference.Included", "Exclusion.Reason"
            ),
            drop = FALSE
            ]
            rownames(summary_df) <- NULL
        }

        summary_df <- applyFitQualityExclusionsMethylationGLM(
            summaryDf = summary_df,
            excludeNonConverged = excludeNonConverged
        )
        diagnostic_summaries[[phenotype]] <- summary_df
        if (nrow(summary_df) > 0L && !is.na(p_value_filter)) {
            keep <- is.finite(summary_df[["Pr(>|t|)"]]) & summary_df[["Pr(>|t|)"]] <
                p_value_filter
            summary_df <- summary_df[keep, , drop = FALSE]
            rownames(summary_df) <- NULL
        }

        summaries[[phenotype]] <- summary_df

        emitLogMinfiEwasWater(
            c(
                "============================================================",
                paste("Summarized phenotype:        ", phenotype),
                paste("CpG summary rows returned:   ", nrow(summary_df)),
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

    structure(list(
        summaries = summaries, diagnosticSummaries = diagnostic_summaries,
        phenotypes = names(summaries), fitFailures = modelResults$fitFailures,
        fitDiagnostics = resolveFitDiagnosticsMethylationGLM(modelSummaries = list(fitDiagnostics = if (!is.null(modelResults$fitDiagnostics)) {
            modelResults$fitDiagnostics
        } else {
            collectFitDiagnosticsMethylationGLM(modelResults$fits)
        }, fitFailures = modelResults$fitFailures),
            summaryList = diagnostic_summaries),
        settings = list(
            summaryResidualSD = isTRUE(summaryResidualSD),
            summaryPval = p_value_filter,
                excludeNonConverged = isTRUE(excludeNonConverged),
            chunkSize = chunk_size
        )
    ), class = "dnaEPICO_methylationGLM_summaries")
}

#' Collect significant CpG coefficient tables from fitted one-timepoint GLMs
#'
#' @param modelResults Object returned by `fitMethylationGLMModels()`.
#' @param pvalThreshold Numeric. Threshold applied to phenotype main-effect or
#'   interaction p-values.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#' @param excludeNonConverged Logical. If `TRUE`, do not collect CpGs whose GLM
#'   did not converge; their rows and reasons remain in diagnostic and annotated
#'   output.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_methylationGLM_significant_cpgs'`.
#'
#' @description
#' Collect the raw coefficient tables for CpGs whose phenotype main effect or
#' interaction p-value passes the requested threshold.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' significant_cpgs <- collectSignificantCpGsMethylationGLM(
#'     modelResults = ex$modelResults,
#'     pvalThreshold = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(significant_cpgs)
#'
#' @export
collectSignificantCpGsMethylationGLM <- function(
    modelResults,
    pvalThreshold = 0.05, interactionTerm = NULL, excludeNonConverged = FALSE,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_methylationGLM.txt"
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
        if (!is.null(modelResults$summaryCache) &&
            !is.null(modelResults$summaryCache[[phenotype]]) &&
            optionalTermMatchesMethylationModels(
                requested = interactionTerm,
                cached = modelResults$settings$interactionTerm
            )) {
            cached_summary <- modelResults$summaryCache[[phenotype]]
            if (nrow(cached_summary) > 0L && !is.na(threshold)) {
                hit_cpgs <-
                    unique(cached_summary$CpG[cached_summary[["Pr(>|t|)"]] <
                    threshold])
                hit_cpgs <- hit_cpgs[!is.na(hit_cpgs) & hit_cpgs %in%
                    names(fit_list)]
                for (cpg in hit_cpgs) {
                    model_obj <- fit_list[[cpg]]
                    if (is.null(model_obj) || inherits(
                        model_obj,
                        "dnaEPICO_methylationGLM_fit_error"
                    )) {
                        next
                    }
                    if (isTRUE(excludeNonConverged) && identical(
                        fitStatusValuesMethylationGLM(model_obj)$converged,
                        FALSE
                    )) {
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
                "dnaEPICO_methylationGLM_fit_error")) {
                next
            }
            if (isTRUE(excludeNonConverged) && identical(
                fitStatusValuesMethylationGLM(model_obj)$converged,
                FALSE
            )) {
                next
            }

            coef_table <- model_obj$coef
            if (is.null(coef_table)) {
                next
            }

            matched_rows <- findCoefficientRowsMethylationGLM(
                coefNames = rownames(coef_table),
                variable = phenotype, interactionTerm = interactionTerm,
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
                "Significant CpGs retained at p <", threshold,
                ":"
            ), paste(names(hit_counts), hit_counts,
                sep = ": ",
                collapse = "; "
            ), "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(retained, class = "dnaEPICO_methylationGLM_significant_cpgs")
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
#'     modelSummaries = ex$modelSummaries,
#'     preparedData = ex$preparedData,
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(diagnostic_plots$plots)
#'
#' @export
plotMethylationGLMDiagnostics <- function(
    modelSummaries, preparedData,
    fdrThreshold = 0.05, padjmethod = "fdr", outputDir = NULL,
    plotWidth = 2000L, plotHeight = 1000L, plotDPI = 150L, display = FALSE,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_methylationGLM.txt"
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
            pValues = summary_df[["Pr(>|t|)"]],
            terms = summary_df$Coefficient, method = padjmethod
        )
        coefficient_terms <- unique(as.character(summary_df$Coefficient))
        coefficient_terms <- coefficient_terms[!is.na(coefficient_terms)]
        multiple_terms <- length(coefficient_terms) > 1L
        phenotype_plots <- phenotype_files <- list()
        phenotype_inflation <- numeric(0)

        for (term_index in seq_along(coefficient_terms)) {
            term <- coefficient_terms[[term_index]]
            term_diagnostics <- buildMethylationTermDiagnosticsDnaEpico(
                summaryData = summary_df,
                phenotype = phenotype, term = term, termColumn = "Coefficient",
                pValueColumn = "Pr(>|t|)", yColumn = "ResidualSD",
                yLabel = "Residual SD", diagnosticMean = diagnostic_mean,
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
                term_files <- list(qqplot = file.path(
                    outputDir,
                    paste0("qqplot_", file_key, ".tiff")
                ), residualSD = if (!is.null(term_diagnostics$plots$residualSD)) {
                    file.path(outputDir, paste0(
                        "residualSD_",
                        file_key, ".tiff"
                    ))
                } else {
                    NULL
                }, residualSignificance = if (!is.null(term_diagnostics$plots$residualSignificance)) {
                    file.path(outputDir, paste0(
                        "residualSignificance_",
                        file_key, ".tiff"
                    ))
                } else {
                    NULL
                })
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
    ), class = "dnaEPICO_methylationGLM_diagnostic_plots")
}

resolveFitDiagnosticsMethylationGLM <- function(
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
        rows[[row_index]] <- unique(data.frame(
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
        row_index <- row_index + 1L
    }

    diagnostics <- modelSummaries$fitDiagnostics
    if (is.data.frame(diagnostics) && all(required_columns %in%
        names(diagnostics))) {
        rows[[row_index]] <- diagnostics[, required_columns,
            drop = FALSE
        ]
    }
    if (length(rows) == 0L) {
        empty <- lapply(required_columns, function(x) character(0))
        names(empty) <- required_columns
        empty$Singular.Fit <- logical(0)
        empty$Converged <- logical(0)
        empty$Inference.Included <- logical(0)
        return(as.data.frame(empty,
            stringsAsFactors = FALSE,
            check.names = FALSE
        ))
    }

    out <- do.call(rbind, rows)
    key <- paste(out$Phenotype, out$CpG, sep = "\r")
    out <- out[!duplicated(key), , drop = FALSE]
    rownames(out) <- NULL
    out
}

buildAnnotationFitDiagnosticsMethylationGLM <- function(fitDiagnostics) {
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
        } else if (any(statuses == "fitted_with_warning")) {
            "fitted_with_warning"
        } else {
            "fitted"
        }
        reasons <-
            unique(paste0(diagnostic$Phenotype[!is.na(diagnostic$Exclusion.Reason) &
            nzchar(diagnostic$Exclusion.Reason)], ": ",
                diagnostic$Exclusion.Reason[!is.na(diagnostic$Exclusion.Reason) &
            nzchar(diagnostic$Exclusion.Reason)]))
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
    by_phenotype <- if (length(phenotype_tables) == 1L) {
        phenotype_tables[[1L]]
    } else {
        Reduce(
            function(x, y) merge(x, y, by = "CpG", all = TRUE),
            phenotype_tables
        )
    }
    list(overall = overall, byPhenotype = by_phenotype)
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
#'     modelSummaries = ex$modelSummaries,
#'     annotationObject = ex$annotationData,
#'     annotationCols = "Name,chr,pos",
#'     verbose = FALSE,
#'     logs = FALSE
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
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_methylationGLM.txt"
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
    fit_diagnostics <- resolveFitDiagnosticsMethylationGLM(
        modelSummaries = modelSummaries,
        summaryList = summary_list
    )
    annotation_diagnostics <-
        buildAnnotationFitDiagnosticsMethylationGLM(fitDiagnostics = fit_diagnostics)

    annotation_cols <- splitOptionMinfiEwasWater(annotationCols,
        sep = ","
    )
    annotation_df <- coerceAnnotationDataMethylationGLM(annotationObject)
    coefficient_occurrences <- unlist(lapply(summary_list, function(x) {
        if (is.null(x) || nrow(x) == 0L || !("Coefficient" %in%
            colnames(x))) {
            return(character(0))
        }
        unique(as.character(x$Coefficient))
    }), use.names = FALSE)
    duplicated_coefficients <-
        unique(coefficient_occurrences[duplicated(coefficient_occurrences)])

    merged_summary_list <- lapply(names(summary_list), function(phenotype) {
        summary_df <- summary_list[[phenotype]]
        if (is.null(summary_df) || nrow(summary_df) == 0L) {
            return(NULL)
        }

        coefficient_names <- unique(summary_df$Coefficient)
        coefficient_tables <- lapply(coefficient_names,
            function(coefficient_name) {
            sub_df <- summary_df[summary_df$Coefficient == coefficient_name,
                c("CpG", "Pr(>|t|)"),
                drop = FALSE
            ]
            clean_name <- gsub("`", "", coefficient_name, fixed = TRUE)
            colnames(sub_df)[2] <- if (coefficient_name %in%
                duplicated_coefficients) {
                paste0(phenotype, "_", clean_name, "_P.Value")
            } else {
                paste0(clean_name, "P.Value")
            }
            sub_df
        })

        if (length(coefficient_tables) == 0L) {
            return(NULL)
        }

        Reduce(
            function(x, y) merge(x, y, by = "CpG", all = TRUE),
            coefficient_tables
        )
    })
    merged_summary_list <- Filter(Negate(is.null), merged_summary_list)

    if (length(merged_summary_list) == 0L) {
        merged_summary <- data.frame(CpG = character(0))
    } else if (length(merged_summary_list) == 1L) {
        merged_summary <- merged_summary_list[[1L]]
    } else {
        merged_summary <- Reduce(function(x, y) {
            merge(x, y,
                by = "CpG",
                all = TRUE
            )
        }, merged_summary_list)
    }
    if (nrow(annotation_diagnostics$byPhenotype) > 0L) {
        glm_diagnostics <- annotation_diagnostics$byPhenotype
        glm_diagnostics <- glm_diagnostics[, !grepl(
            "_Singular\\.Fit$",
            colnames(glm_diagnostics)
        ), drop = FALSE]
        merged_summary <- merge(merged_summary, glm_diagnostics,
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
        annotationCols = available_annotation_cols, includeSingular = FALSE
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
    ), class = "dnaEPICO_methylationGLM_annotation")
}

#' Write optional disk outputs for one-timepoint GLM analyses
#'
#' @param modelResults Object returned by `fitMethylationGLMModels()`.
#' @param modelSummaries Object returned by `summarizeMethylationGLMModels()`.
#' @param annotatedResults Object returned by
#'   `annotateMethylationGLMSummaries()` or a compatible data frame.
#' @param significantCpGs Object returned by
#'   `collectSignificantCpGsMethylationGLM()` or `NULL`.
#' @param outputRData Character. Directory used for serialized model and summary
#'   outputs.
#' @param summaryTxtDir Character. Directory used for tab-delimited summary
#'   tables.
#' @param significantCpGDir Character. Directory used for significant-CpG
#'   coefficient tables.
#' @param annotatedGLMOut Character. Directory used for the annotated summary
#'   XLSX workbook.
#' @param reportAssetsDir Character or `NULL`. Report results directory used for
#'   the compressed TSV table and compact metadata sidecars.
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
#' Write optional serialized outputs, summary tables, significant-CpG tables,
#' and annotated results from the one-timepoint GLM workflow.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' annotation_data <- annotateMethylationGLMSummaries(
#'     modelSummaries = ex$modelSummaries,
#'     annotationObject = ex$annotationData,
#'     annotationCols = "Name,chr,pos",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' significant_cpgs <- collectSignificantCpGsMethylationGLM(
#'     modelResults = ex$modelResults,
#'     pvalThreshold = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' output_paths <- writeMethylationGLMOutputs(
#'     modelResults = ex$modelResults,
#'     modelSummaries = ex$modelSummaries,
#'     annotatedResults = annotation_data,
#'     significantCpGs = significant_cpgs,
#'     outputRData = file.path(ex$tempDir, "models"),
#'     summaryTxtDir = file.path(ex$tempDir, "summary"),
#'     significantCpGDir = file.path(ex$tempDir, "significant"),
#'     annotatedGLMOut = file.path(ex$tempDir, "annotated"),
#'     saveTxtSummaries = TRUE,
#'     saveSignificantCpGs = TRUE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(output_paths)
#'
#' @export
writeMethylationGLMOutputs <- function(
    modelResults, modelSummaries,
    annotatedResults, significantCpGs = NULL, outputRData, summaryTxtDir,
    significantCpGDir, annotatedGLMOut, reportAssetsDir = NULL,
    saveTxtSummaries = TRUE,
    saveSignificantCpGs = FALSE, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_methylationGLM.txt"
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
    dir.create(annotatedGLMOut, recursive = TRUE, showWarnings = FALSE)

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
            "GLM.rds"
        ))
        summary_file <- file.path(outputRData, paste0(
            phenotype,
            "SummaryGLM.rds"
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

            if ("Pr(>|t|)" %in% colnames(summary_df)) {
                summary_df <- summary_df[order(summary_df[["Pr(>|t|)"]]), ,
                    drop = FALSE
                ]
            }

            output_file <- file.path(summaryTxtDir, paste0(
                phenotype,
                "SummaryGLM.txt"
            ))
            utils::write.table(summary_df,
                file = output_file,
                sep = "\t", row.names = FALSE, quote = FALSE
            )
            summary_txt_files[[phenotype]] <- output_file
        }
    }

    if (isTRUE(saveSignificantCpGs) && !is.null(significantCpGs)) {
        dir.create(significantCpGDir, recursive = TRUE, showWarnings = FALSE)
        for (phenotype in names(significantCpGs)) {
            phenotype_hits <- significantCpGs[[phenotype]]
            if (length(phenotype_hits) == 0L) {
                next
            }

            significant_files[[phenotype]] <- character(0)
            phenotype_dir <- file.path(significantCpGDir, phenotype)
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
    annotated_file <- file.path(annotatedGLMOut, "annotatedGLM.xlsx")
    fit_failures <- modelSummaries$fitFailures
    if (is.null(fit_failures)) {
        fit_failures <- modelResults$fitFailures
    }
    if (is.data.frame(fit_failures) && nrow(fit_failures) > 0L) {
        fit_failure_file <- file.path(annotatedGLMOut, "CpGFitFailuresGLM.txt")
        utils::write.table(fit_failures,
            file = fit_failure_file,
            sep = "\t", row.names = FALSE, quote = FALSE
        )
    }
    dictionary <- buildAnnotatedWorkbookDictionaryMethylationGLM(
        columns = colnames(annotated_df),
        modelDescription = "Pvalue from GLM model",
            formulaText = modelResults$formulas,
        modelLabel = "GLM",
            responseLabel = inferMethylationValueLabelMethylationGLM(modelResults)
    )
    metadata <- buildModelWorkbookMetadataDnaEpico(
        modelResults = modelResults,
        modelSummaries = modelSummaries, annotatedResults = annotatedResults,
        analysis = "glm"
    )
    writeAnnotatedWorkbookMethylationGLM(
        annotated_df = annotated_df,
        file = annotated_file, resultSheet = "annotatedGLM",
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
            sheet = "annotatedGLM", idColumn = report_table$idColumn,
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
            }, if (isTRUE(saveSignificantCpGs)) {
                paste("Significant CpG text files:  ", sum(vapply(
                    significant_files,
                    length, integer(1)
                )))
            } else {
                "Significant CpG text files:  0"
            }, "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(
        list(
            modelFiles = model_files, summaryFiles = summary_files,
            summaryTxtFiles = summary_txt_files,
                significantCpGFiles = significant_files,
            fitFailureFile = fit_failure_file, annotatedGLM = annotated_file,
            annotatedGLMText = report_sidecar$table,
                annotatedGLMReportMetadata = report_sidecar$metadata,
            annotatedGLMDictionary = report_sidecar$dictionary,
                annotatedGLMMetadata = report_sidecar$workbookMetadata
        ),
        class = "dnaEPICO_methylationGLM_paths"
    )
}
