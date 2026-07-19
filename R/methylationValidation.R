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
            responseLabel = "Beta values", responseColumn = "beta",
            domain = "[0, 1]"
        ), m = list(
            objectPrefix = "phenoM",
            responseLabel = "M-values", responseColumn = "m",
                domain = "Unbounded; finite values required for modeling"
        ),
        cn = list(
            objectPrefix = "phenoCN", responseLabel = "Copy number values",
            responseColumn = "cn",
                domain = "Unbounded; finite values required for modeling"
        )
    )
    definitions[[normalizeMethylationScaleDnaEpico(methylationScale)]]
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

methylationScaleDomainDnaEpico <- function(methylationScale = "beta") {
    methylationScaleDefinitionDnaEpico(methylationScale)$domain
}

buildMethylationValidationDnaEpico <- function(
    cpgIDs, methylationScale,
    numericCpGs, finiteCounts, naCounts, nanCounts, positiveInfCounts,
    negativeInfCounts, outOfRangeCounts, lowerBoundaryCounts,
    upperBoundaryCounts, observedMinimum, observedMaximum,
        overallMinimum = NA_real_,
    overallMaximum = NA_real_
) {
    methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
    n_cpgs <- length(cpgIDs)
    invalid_cpgs <- !numericCpGs | positiveInfCounts > 0L | negativeInfCounts >
        0L | outOfRangeCounts > 0L | (numericCpGs & finiteCounts ==
        0L)
    issue_cpgs <- invalid_cpgs | nanCounts > 0L
    issue_indices <- which(issue_cpgs)

    reason_rules <- list(
        `CpG column is not numeric` = !numericCpGs,
        `NaN converted to NA` = nanCounts > 0L,
            `infinite methylation value` = positiveInfCounts >
            0L | negativeInfCounts > 0L,
                `beta value outside [0, 1]` = outOfRangeCounts >
            0L, `no finite methylation values` = numericCpGs &
            finiteCounts == 0L
    )
    reasons <- vapply(issue_indices, function(index) {
        active <- vapply(reason_rules, `[[`, logical(1), index)
        paste(names(reason_rules)[active], collapse = "; ")
    }, character(1))

    issues <- data.frame(
        CpG = cpgIDs[issue_indices], Scale = rep(
            methylation_scale,
            length(issue_indices)
        ), Status = ifelse(invalid_cpgs[issue_indices],
            "invalid", "NaN converted to NA"
        ), Observed.Minimum = observedMinimum[issue_indices],
        Observed.Maximum = observedMaximum[issue_indices],
            Finite.Values = finiteCounts[issue_indices],
        NA.Values = naCounts[issue_indices],
            NaN.Converted = nanCounts[issue_indices],
        Positive.Inf = positiveInfCounts[issue_indices],
            Negative.Inf = negativeInfCounts[issue_indices],
        Out.Of.Range = outOfRangeCounts[issue_indices], Reason = reasons,
        stringsAsFactors = FALSE, check.names = FALSE
    )

    beta_scale <- identical(methylation_scale, "beta")
    boundaries <- data.frame(
        Scale = methylation_scale,
            Defined.Domain = methylationScaleDomainDnaEpico(methylation_scale),
        Observed.Minimum = overallMinimum, Observed.Maximum = overallMaximum,
        Finite.Values = sum(finiteCounts), NA.Values = sum(naCounts),
        NaN.Converted = sum(nanCounts), Positive.Inf = sum(positiveInfCounts),
        Negative.Inf = sum(negativeInfCounts),
            Out.Of.Range = sum(outOfRangeCounts),
        At.Lower.Boundary = if (beta_scale) {
            sum(lowerBoundaryCounts)
        } else {
            NA_real_
        }, At.Upper.Boundary = if (beta_scale) {
            sum(upperBoundaryCounts)
        } else {
            NA_real_
        }, CpGs = n_cpgs, Numeric.CpGs = sum(numericCpGs),
        Invalid.CpGs = sum(invalid_cpgs), stringsAsFactors = FALSE,
        check.names = FALSE
    )

    list(boundaries = boundaries, issues = issues,
        invalidCpGs = issues[issues$Status ==
        "invalid", , drop = FALSE])
}

inspectMethylationMatrixDnaEpico <- function(values,
    methylationScale = "beta") {
    methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
    if (!is.matrix(values) || !is.numeric(values)) {
        stop("Methylation values must be supplied as a numeric matrix.",
            call. = FALSE
        )
    }

    cpg_ids <- rownames(values)
    if (is.null(cpg_ids)) {
        cpg_ids <- as.character(seq_len(nrow(values)))
    }

    nan_by_cpg <- rowSums(is.nan(values))
    values[is.nan(values)] <- NA_real_
    finite_index <- is.finite(values)
    finite_by_cpg <- rowSums(finite_index)
    na_by_cpg <- rowSums(is.na(values))
    positive_inf_by_cpg <- rowSums(is.infinite(values) & values >
        0)
    negative_inf_by_cpg <- rowSums(is.infinite(values) & values <
        0)
    beta_scale <- identical(methylation_scale, "beta")
    out_of_range_by_cpg <- if (beta_scale) {
        rowSums(finite_index & (values < 0 | values > 1))
    } else {
        integer(nrow(values))
    }
    lower_boundary_by_cpg <- if (beta_scale) {
        rowSums(values == 0, na.rm = TRUE)
    } else {
        integer(nrow(values))
    }
    upper_boundary_by_cpg <- if (beta_scale) {
        rowSums(values == 1, na.rm = TRUE)
    } else {
        integer(nrow(values))
    }

    finite_values <- values[finite_index]
    overall_minimum <- if (length(finite_values)) {
        min(finite_values)
    } else {
        NA_real_
    }
    overall_maximum <- if (length(finite_values)) {
        max(finite_values)
    } else {
        NA_real_
    }
    observed_minimum <- rep(NA_real_, nrow(values))
    observed_maximum <- rep(NA_real_, nrow(values))
    issue_indices <- which(nan_by_cpg > 0L | positive_inf_by_cpg >
        0L | negative_inf_by_cpg > 0L | out_of_range_by_cpg >
        0L | finite_by_cpg == 0L)
    for (index in issue_indices) {
        finite_row <- values[index, is.finite(values[index, ]),
            drop = TRUE
        ]
        if (length(finite_row)) {
            observed_minimum[[index]] <- min(finite_row)
            observed_maximum[[index]] <- max(finite_row)
        }
    }

    validation <- buildMethylationValidationDnaEpico(
        cpgIDs = cpg_ids,
        methylationScale = methylation_scale, numericCpGs = rep(
            TRUE,
            nrow(values)
        ), finiteCounts = finite_by_cpg, naCounts = na_by_cpg,
        nanCounts = nan_by_cpg, positiveInfCounts = positive_inf_by_cpg,
        negativeInfCounts = negative_inf_by_cpg,
            outOfRangeCounts = out_of_range_by_cpg,
        lowerBoundaryCounts = lower_boundary_by_cpg,
            upperBoundaryCounts = upper_boundary_by_cpg,
        observedMinimum = observed_minimum, observedMaximum = observed_maximum,
        overallMinimum = overall_minimum, overallMaximum = overall_maximum
    )
    c(list(values = values), validation)
}

inspectMethylationColumnsDnaEpico <- function(
    data, cpgColumns,
    methylationScale = "beta"
) {
    methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
    n_cpgs <- length(cpgColumns)
    numeric_cpgs <- vapply(data[cpgColumns], is.numeric, logical(1))
    finite_counts <- na_counts <- nan_counts <- positive_inf_counts <-
        negative_inf_counts <- out_of_range_counts <- lower_boundary_counts <-
        upper_boundary_counts <- integer(n_cpgs)
    observed_minimum <- observed_maximum <- rep(NA_real_, n_cpgs)
    beta_scale <- identical(methylation_scale, "beta")

    for (index in which(numeric_cpgs)) {
        cpg <- cpgColumns[[index]]
        values <- data[[cpg]]
        nan_index <- is.nan(values)
        nan_counts[[index]] <- sum(nan_index)
        values[nan_index] <- NA_real_
        data[[cpg]] <- values

        finite_values <- values[is.finite(values)]
        finite_counts[[index]] <- length(finite_values)
        na_counts[[index]] <- sum(is.na(values))
        positive_inf_counts[[index]] <- sum(is.infinite(values) &
            values > 0)
        negative_inf_counts[[index]] <- sum(is.infinite(values) &
            values < 0)
        if (length(finite_values)) {
            observed_minimum[[index]] <- min(finite_values)
            observed_maximum[[index]] <- max(finite_values)
        }
        if (beta_scale) {
            out_of_range_counts[[index]] <- sum(finite_values <
                0 | finite_values > 1)
            lower_boundary_counts[[index]] <- sum(finite_values ==
                0)
            upper_boundary_counts[[index]] <- sum(finite_values ==
                1)
        }
    }
    na_counts[!numeric_cpgs] <- vapply(
        data[cpgColumns[!numeric_cpgs]],
        function(values) sum(is.na(values)), integer(1)
    )

    finite_minimum <- observed_minimum[is.finite(observed_minimum)]
    finite_maximum <- observed_maximum[is.finite(observed_maximum)]
    validation <- buildMethylationValidationDnaEpico(
        cpgIDs = cpgColumns,
        methylationScale = methylation_scale, numericCpGs = numeric_cpgs,
        finiteCounts = finite_counts, naCounts = na_counts,
            nanCounts = nan_counts,
        positiveInfCounts = positive_inf_counts,
            negativeInfCounts = negative_inf_counts,
        outOfRangeCounts = out_of_range_counts,
            lowerBoundaryCounts = lower_boundary_counts,
        upperBoundaryCounts = upper_boundary_counts,
            observedMinimum = observed_minimum,
        observedMaximum = observed_maximum,
            overallMinimum = if (length(finite_minimum)) {
            min(finite_minimum)
        } else {
            NA_real_
        }, overallMaximum = if (length(finite_maximum)) {
            max(finite_maximum)
        } else {
            NA_real_
        }
    )
    c(list(data = data), validation)
}

formatMethylationBoundariesLogDnaEpico <- function(boundaries) {
    if (!is.data.frame(boundaries) || nrow(boundaries) == 0L) {
        return("Methylation boundaries: unavailable")
    }
    row <- boundaries[1L, , drop = FALSE]
    c(
        paste("Methylation scale:           ", row$Scale), paste(
            "Defined domain:              ",
            row$Defined.Domain
        ), paste(
            "Observed finite minimum:     ",
            row$Observed.Minimum
        ), paste(
            "Observed finite maximum:     ",
            row$Observed.Maximum
        ), paste(
            "Values at lower boundary:    ",
            row$At.Lower.Boundary
        ), paste(
            "Values at upper boundary:    ",
            row$At.Upper.Boundary
        ), paste(
            "NaN values converted to NA:  ",
            row$NaN.Converted
        ), paste(
            "Positive/negative infinity:  ",
            paste(row$Positive.Inf, row$Negative.Inf, sep = "/")
        ),
        paste("Out-of-range beta values:    ", row$Out.Of.Range),
        paste("Invalid CpGs:                ", row$Invalid.CpGs)
    )
}

invalidMethylationReasonsDnaEpico <- function(validation) {
    invalid <- validation$invalidCpGs
    if (!is.data.frame(invalid) || nrow(invalid) == 0L) {
        return(stats::setNames(character(0), character(0)))
    }
    stats::setNames(as.character(invalid$Reason), as.character(invalid$CpG))
}

newMethylationFitErrorDnaEpico <- function(
    reason, errorClass,
    status = "failed", converged = NA, convergenceMessage = NA_character_
) {
    structure(
        list(
            error = as.character(reason), status = as.character(status),
            converged = as.logical(converged),
                convergenceMessage = as.character(convergenceMessage)
        ),
        class = errorClass
    )
}
