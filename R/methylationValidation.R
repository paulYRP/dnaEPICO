# Internal methylation-scale definitions and value
# validation

validateProbabilityDnaEpico <- function(value, name, allowNA = FALSE) {
    if (length(value) != 1L) {
        stop(name, " must be a single number between 0 and 1.",
            call. = FALSE
        )
    }
    numeric_value <- utils::type.convert(as.character(value),
        as.is = TRUE,
        na.strings = character()
    )
    if (isTRUE(allowNA) && is.na(numeric_value)) {
        return(NA_real_)
    }
    if (!is.numeric(numeric_value) || is.na(numeric_value) ||
        !is.finite(numeric_value) || numeric_value < 0 || numeric_value >
        1) {
        stop(name, " must be a single number between 0 and 1.",
            call. = FALSE
        )
    }
    numeric_value
}

validateCpgPrefixDnaEpico <- function(cpgPrefix) {
    if (!is.character(cpgPrefix) || length(cpgPrefix) != 1L ||
        is.na(cpgPrefix) || !nzchar(cpgPrefix)) {
        stop("cpgPrefix must be one non-empty character string.",
            call. = FALSE
        )
    }
    cpgPrefix
}

normalizeScaleVariablesDnaEpico <- function(scaleVars = NULL) {
    variables <- unique(splitOptionMinfiEwasWater(scaleVars,
        sep = ","
    ))
    variables <- variables[!is.na(variables) & nzchar(variables) &
        !(tolower(variables) %in% c("null", "none", "na"))]
    variables
}

scaleModelVariablesDnaEpico <- function(
    data, scaleVars = NULL,
    factorVars = NULL, eligibleVars, protectedVars = NULL
) {
    scale_vars <- normalizeScaleVariablesDnaEpico(scaleVars)
    factor_vars <- unique(splitOptionMinfiEwasWater(factorVars,
        sep = ","
    ))
    eligible_vars <- unique(as.character(eligibleVars))
    protected_vars <- unique(as.character(protectedVars))

    empty_metadata <- data.frame(
        Variable = character(0), Center = numeric(0),
        Scale = numeric(0), Finite.Values = integer(0),
            Missing.Values = integer(0),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    if (length(scale_vars) == 0L) {
        return(list(data = data, scaleVars = scale_vars,
            metadata = empty_metadata))
    }

    missing_vars <- setdiff(scale_vars, colnames(data))
    if (length(missing_vars) > 0L) {
        stop("Scale columns not found in inputPheno: ", paste(missing_vars,
            collapse = ", "
        ), call. = FALSE)
    }
    factor_overlap <- intersect(scale_vars, factor_vars)
    if (length(factor_overlap) > 0L) {
        stop("Variables cannot be listed in both factorVars and scaleVars: ",
            paste(factor_overlap, collapse = ", "),
            call. = FALSE
        )
    }
    protected_overlap <- intersect(scale_vars, protected_vars)
    if (length(protected_overlap) > 0L) {
        stop("Identifier or protected variables cannot be scaled: ",
            paste(protected_overlap, collapse = ", "),
            call. = FALSE
        )
    }
    unused_vars <- setdiff(scale_vars, eligible_vars)
    if (length(unused_vars) > 0L) {
        stop("scaleVars must contain fixed-effect variables used by at least one model: ",
            paste(unused_vars, collapse = ", "),
            call. = FALSE
        )
    }

    model_data <- data
    metadata_rows <- vector("list", length(scale_vars))
    for (index in seq_along(scale_vars)) {
        variable <- scale_vars[[index]]
        values <- model_data[[variable]]
        if (!is.numeric(values)) {
            stop("scaleVars must be numeric: ", variable, call. = FALSE)
        }
        if (any(is.infinite(values))) {
            stop("scaleVars contains infinite values: ", variable,
                call. = FALSE
            )
        }
        finite_values <- values[is.finite(values)]
        if (length(finite_values) < 2L) {
            stop("scaleVars requires at least two finite observations: ",
                variable,
                call. = FALSE
            )
        }
        center <- mean(finite_values)
        scale_value <- stats::sd(finite_values)
        if (!is.finite(scale_value) || scale_value <= 0) {
            stop("scaleVars has zero or undefined variance: ",
                variable,
                call. = FALSE
            )
        }
        transformed <- values
        finite_index <- is.finite(values)
        transformed[finite_index] <- (values[finite_index] -
            center) / scale_value
        model_data[[variable]] <- transformed
        metadata_rows[[index]] <- data.frame(
            Variable = variable,
            Center = center, Scale = scale_value,
                Finite.Values = length(finite_values),
            Missing.Values = sum(is.na(values)), stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }

    list(data = model_data, scaleVars = scale_vars, metadata = do.call(
        rbind,
        metadata_rows
    ))
}

formatScalingMetadataLogDnaEpico <- function(metadata) {
    if (!is.data.frame(metadata) || nrow(metadata) == 0L) {
        return("Scaled variables:             none")
    }
    c(
        paste("Scaled variables:            ", paste(metadata$Variable,
            collapse = ", "
        )), "Scaling parameters (sample SD):",
        previewLinesMinfiEwasWater(metadata)
    )
}

validateMethylationProbeIdentifiersDnaEpico <- function(
    probeIDs,
    label = "Methylation probe identifiers"
) {
    probe_ids <- as.character(probeIDs)
    if (length(probe_ids) == 0L) {
        stop(label, " must contain at least one CpG identifier.",
            call. = FALSE
        )
    }
    if (anyNA(probe_ids) || any(!nzchar(trimws(probe_ids)))) {
        stop(label, " contains missing or blank CpG identifiers.",
            call. = FALSE
        )
    }
    duplicated_ids <- unique(probe_ids[duplicated(probe_ids)])
    if (length(duplicated_ids) > 0L) {
        stop(label, " contains duplicate CpG identifiers: ",
            paste(duplicated_ids, collapse = ", "),
            call. = FALSE
        )
    }
    probe_ids
}

normalizeMethylationScaleDnaEpico <- function(methylationScale = "beta") {
    choices <- c("beta", "m", "cn")
    if (!is.character(methylationScale) || length(methylationScale) !=
        1L || is.na(methylationScale)) {
        stop("methylationScale must be one of: Beta, M, CN (case-insensitive).",
            call. = FALSE
        )
    }

    methylation_scale <- tolower(trimws(methylationScale))
    if (!(methylation_scale %in% choices)) {
        stop("methylationScale must be one of: Beta, M, CN (case-insensitive).",
            call. = FALSE
        )
    }

    methylation_scale
}

methylationScaleDefinitionDnaEpico <- function(methylationScale = "beta") {
    definitions <- list(
        beta = list(
            objectPrefix = "phenoB",
            responseLabel = "Beta values", responseColumn = "beta"
        ), m = list(
            objectPrefix = "phenoM",
            responseLabel = "M-values", responseColumn = "m"
        ),
        cn = list(
            objectPrefix = "phenoCN", responseLabel = "Copy number values",
            responseColumn = "cn"
        )
    )
    definitions[[normalizeMethylationScaleDnaEpico(methylationScale)]]
}

summarizeMethylationRangeDnaEpico <- function(
    values, methylationScale = "beta"
) {
    methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
    if (!is.matrix(values) || !is.numeric(values)) {
        stop("Methylation values must be supplied as a numeric matrix.",
            call. = FALSE
        )
    }

    finite_minimum <- Inf
    finite_maximum <- -Inf
    for (column in seq_len(ncol(values))) {
        column_values <- values[, column]
        finite_values <- column_values[is.finite(column_values)]
        if (length(finite_values) == 0L) {
            next
        }
        finite_minimum <- min(finite_minimum, min(finite_values))
        finite_maximum <- max(finite_maximum, max(finite_values))
    }

    if (!is.finite(finite_minimum)) {
        finite_minimum <- NA_real_
        finite_maximum <- NA_real_
    }

    data.frame(
        Scale = methylation_scale,
        Observed.Minimum = finite_minimum,
        Observed.Maximum = finite_maximum,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
}

formatMethylationRangeLogDnaEpico <- function(rangeSummary) {
    if (!is.data.frame(rangeSummary) || nrow(rangeSummary) == 0L) {
        return(c(
            "Methylation scale:            unavailable",
            "Observed finite minimum:      NA",
            "Observed finite maximum:      NA"
        ))
    }
    row <- rangeSummary[1L, , drop = FALSE]
    c(
        paste("Methylation scale:           ", row$Scale),
        paste("Observed finite minimum:     ", row$Observed.Minimum),
        paste("Observed finite maximum:     ", row$Observed.Maximum)
    )
}

methylationScaleObjectPrefixDnaEpico <- function(methylationScale = "beta") {
    methylationScaleDefinitionDnaEpico(methylationScale)$objectPrefix
}

methylationScaleResponseLabelDnaEpico <- function(methylationScale = "beta") {
    methylationScaleDefinitionDnaEpico(methylationScale)$responseLabel
}

methylationScaleResponseColumnDnaEpico <- function(methylationScale = "beta") {
    methylationScaleDefinitionDnaEpico(methylationScale)$responseColumn
}

newMethylationFitErrorDnaEpico <- function(
    reason, errorClass, modelMessage = NULL
) {
    resolved_message <- if (is.null(modelMessage)) {
        paste0("ERROR: ", as.character(reason))
    } else {
        as.character(modelMessage)
    }
    structure(
        list(
            error = as.character(reason),
            modelMessage = resolved_message,
            pValueAvailable = FALSE
        ),
        class = errorClass
    )
}
