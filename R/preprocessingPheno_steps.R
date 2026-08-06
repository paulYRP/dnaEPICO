#' Parse timepoint selections for preprocessingPheno helpers
#'
#' @param values Character vector or comma-separated string of requested
#'   timepoints.
#' @param label Character. Label used in validation messages.
#'
#' @return Character vector of parsed timepoints.
#'
#' @description
#' Internal helper that standardizes timepoint parsing for
#' `preprocessingPheno()` and related helper functions.
#'
#' @keywords internal
#' @noRd
parseTimepointsPreprocessingPheno <- function(values, label = "timepoints") {
    parsed <- splitOptionMinfiEwasWater(values, sep = ",")

    if (length(parsed) == 0L) {
        stop(label, " must contain at least one timepoint.",
            call. = FALSE
        )
    }

    unique(parsed)
}

#' Extract a preferred named element from a loaded object
#'
#' @param object Loaded object.
#' @param preferred_name Character vector or `NULL`. Preferred element names,
#'   in selection order, when `object` is a named list.
#'
#' @return The original object or its preferred named element.
#'
#' @description
#' Internal helper that unwraps list-style saved objects used by
#' `preprocessingPheno()` and downstream wrappers.
#'
#' @keywords internal
#' @noRd
extractPreferredObjectPreprocessingPheno <- function(
    object,
    preferred_name = NULL
) {
    if (is.null(preferred_name)) {
        return(object)
    }

    preferred_names <- unique(as.character(preferred_name))
    preferred_names <- preferred_names[!is.na(preferred_names) &
        nzchar(preferred_names)]
    matching_names <- preferred_names[preferred_names %in% names(object)]
    if (is.list(object) && !is.data.frame(object) && length(matching_names) >
        0L) {
        return(object[[matching_names[[1L]]]])
    }

    object
}

#' Load a saved object used by preprocessingPheno helpers
#'
#' @param path Character. Path to an `.RData` or `.rds` file.
#' @param preferred_name Character vector or `NULL`. Preferred object names,
#'   in selection order, when `path` points to an `.RData` file containing
#'   multiple objects.
#'
#' @return The loaded object.
#'
#' @description
#' Internal helper that loads a saved object from disk while keeping
#' `preprocessingPheno()` compatible with the package's existing `.RData`
#' outputs.
#'
#' @keywords internal
#' @noRd
loadSavedObjectPreprocessingPheno <- function(path, preferred_name = NULL) {
    if (!file.exists(path)) {
        stop("Input file does not exist: ", path, call. = FALSE)
    }

    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        return(extractPreferredObjectPreprocessingPheno(
            object = readRDS(path),
            preferred_name = preferred_name
        ))
    }

    load_env <- new.env(parent = emptyenv())
    loaded_names <- load(path, envir = load_env)

    preferred_names <- unique(as.character(preferred_name))
    preferred_names <- preferred_names[!is.na(preferred_names) &
        nzchar(preferred_names)]
    matching_names <- preferred_names[preferred_names %in% loaded_names]
    if (length(matching_names) > 0L) {
        selected_name <- matching_names[[1L]]
        return(extractPreferredObjectPreprocessingPheno(
            object = load_env[[selected_name]],
            preferred_name = preferred_names
        ))
    }

    if (length(loaded_names) == 1L) {
        return(extractPreferredObjectPreprocessingPheno(
            object = load_env[[loaded_names[[1L]]]],
            preferred_name = preferred_name
        ))
    }

    stop("Could not determine which object to load from ", path,
        ". Available objects: ", paste(loaded_names, collapse = ", "),
        call. = FALSE
    )
}

#' Merge phenotype rows with methylation values for preprocessingPheno helpers
#'
#' @param phenoFrame Data frame containing phenotype rows.
#' @param methylationMatrix Numeric matrix with probes in rows and samples in
#'   columns.
#' @param id Character. Name of the sample identifier column in `phenoFrame`.
#' @param methylationScale Character. Selected methylation scale: `'Beta'`,
#'   `'M'`, or `'CN'`, case-insensitively.
#'
#' @return A data frame with phenotype columns followed by transposed
#' methylation
#'   values.
#'
#' @description
#' Internal helper that aligns phenotype rows with methylation-matrix columns
#' before
#' building merged objects used by downstream modeling functions.
#'
#' @keywords internal
#' @noRd
mergePhenoMethylationPreprocessingPheno <- function(
    phenoFrame,
    methylationMatrix, id = "Sample_Name", methylationScale = "beta"
) {
    methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
    methylation_label <-
        methylationScaleResponseLabelDnaEpico(methylation_scale)

    if (!(id %in% colnames(phenoFrame))) {
        stop("Sample identifier column not found in phenotype data: ",
            id,
            call. = FALSE
        )
    }

    if (!is.matrix(methylationMatrix) || !is.numeric(methylationMatrix) ||
        is.null(rownames(methylationMatrix)) ||
            is.null(colnames(methylationMatrix))) {
        stop(methylation_label,
            " must be a numeric matrix with sample identifiers in column names.",
            call. = FALSE
        )
    }
    validateMethylationProbeIdentifiersDnaEpico(
        rownames(methylationMatrix),
        paste0(methylation_label, " row names")
    )

    sample_ids <- validateSampleIdentifiersDnaEpico(
        phenoFrame[[id]],
        paste0("phenotype column '", id, "'")
    )
    matrix_ids <- validateSampleIdentifiersDnaEpico(
        colnames(methylationMatrix),
        paste0(methylation_label, " column names")
    )
    missing_samples <- setdiff(sample_ids, matrix_ids)
    if (length(missing_samples) > 0L) {
        stop("Phenotype samples are missing from ", methylation_label,
            ": ", paste(utils::head(missing_samples, 10L), collapse = ", "),
            call. = FALSE
        )
    }
    matched_ids <- sample_ids

    if (length(matched_ids) == 0L) {
        stop("No matching sample identifiers were found between phenotype data and ",
            methylation_label, ".",
            call. = FALSE
        )
    }

    if (anyDuplicated(matched_ids) > 0L) {
        stop("Duplicate sample identifiers were found while merging phenotype and methylation data: ",
            paste(unique(matched_ids[duplicated(matched_ids)]),
                collapse = ", "
            ),
            call. = FALSE
        )
    }

    phenoFrame <- phenoFrame[match(matched_ids, sample_ids), ,
        drop = FALSE
    ]
    methylationMatrix <- methylationMatrix[, matched_ids, drop = FALSE]

    cbind(phenoFrame, as.data.frame(t(methylationMatrix),
        stringsAsFactors = FALSE))
}

#' Load methylation metric matrices for preprocessingPheno
#'
#' @param betaPath Character. Path to the saved beta-value object.
#' @param mPath Character. Path to the saved M-value object.
#' @param cnPath Character. Path to the saved copy-number object.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_preprocessingPheno_metrics'`
#'   containing `beta`, `m`, and `cn`.
#'
#' @examples
#' ex <- dnaEPICO:::examplePreprocessingPhenoStateDnaEpico()
#' metrics_data <- loadMetricsPreprocessingPheno(
#'     betaPath = ex$betaPath,
#'     mPath = ex$mPath,
#'     cnPath = ex$cnPath,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(metrics_data)
#'
#' @description
#' Load the metric matrices generated by `preprocessingMinfiEwasWater()` and
#' return them as a single in-memory object for downstream phenotype alignment.
#'
#' @export
loadMetricsPreprocessingPheno <- function(
    betaPath, mPath, cnPath,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_loadMetricsPreprocessingPheno.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    beta <- loadSavedObjectPreprocessingPheno(betaPath, preferred_name = "beta")
    m <- loadSavedObjectPreprocessingPheno(mPath, preferred_name = "m")
    cn <- loadSavedObjectPreprocessingPheno(cnPath, preferred_name = "cn")

    metric_objects <- list(beta = beta, m = m, cn = cn)
    invalid_metrics <- names(metric_objects)[!vapply(
        metric_objects,
        function(x) is.matrix(x) && is.numeric(x), logical(1)
    )]
    if (length(invalid_metrics) > 0L) {
        stop("Metric objects must be numeric matrices: ", paste(invalid_metrics,
            collapse = ", "
        ), call. = FALSE)
    }

    if (is.null(rownames(beta)) || is.null(colnames(beta))) {
        stop("The beta matrix must have probe and sample names.",
            call. = FALSE
        )
    }
    validateMethylationProbeIdentifiersDnaEpico(
        rownames(beta),
        "Beta-matrix row names"
    )
    validateSampleIdentifiersDnaEpico(colnames(beta),
        "Beta-matrix sample identifiers")
    for (metric_name in c("m", "cn")) {
        metric <- metric_objects[[metric_name]]
        if (!identical(dim(metric), dim(beta)) || !identical(
            rownames(metric),
            rownames(beta)
        ) || !identical(colnames(metric), colnames(beta))) {
            stop(metric_name,
                " must have the same dimensions, probe order, and sample order as beta.",
                call. = FALSE
            )
        }
    }

    range_summaries <- lapply(names(metric_objects), function(metric_name) {
        summarizeMethylationRangeDnaEpico(
            values = metric_objects[[metric_name]],
            methylationScale = metric_name
        )
    })
    names(range_summaries) <- names(metric_objects)

    preview_rows <- seq_len(min(nrow(beta), 5L))
    preview_cols <- seq_len(min(ncol(beta), 5L))

    emitLogMinfiEwasWater(
        c(
            "============================================================",
            paste("Beta path:                ", betaPath), paste(
                "M-values path:            ",
                mPath
            ), paste("CN path:                  ", cnPath),
            paste("Beta dimensions:          ", paste(dim(beta),
                collapse = " x "
            )), paste(
                "M dimensions:             ",
                paste(dim(m), collapse = " x ")
            ), paste(
                "CN dimensions:            ",
                paste(dim(cn), collapse = " x ")
            ), unlist(lapply(
                range_summaries,
                function(range_summary) {
                    formatMethylationRangeLogDnaEpico(range_summary)
                }
            ), use.names = FALSE), "Preview of beta values:",
            previewLinesMinfiEwasWater(beta[preview_rows, preview_cols,
                drop = FALSE
            ]), "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(
        beta = beta, m = m, cn = cn,
        methylationRanges = range_summaries
    ),
        class = "dnaEPICO_preprocessingPheno_metrics"
    )
}

#' Split phenotype and methylation data by timepoint
#'
#' @param pheno Data frame containing phenotype information.
#' @param metricsData Object returned by `loadMetricsPreprocessingPheno()`.
#' @param SampleID Character. Name of the sample identifier column in `pheno`.
#' @param timeVar Character. Name of the timepoint column in `pheno`.
#' @param timepoints Character vector or comma-separated string of timepoints to
#'   retain.
#' @param methylationScale Character. Methylation metric to use for the merged
#'   modeling table. One of `'Beta'`, `'M'`, or `'CN'`, case-insensitively.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_preprocessingPheno_timepoints'`
#'   containing the parsed timepoints and aligned per-timepoint subsets.
#'
#' @examples
#' ex <- dnaEPICO:::examplePreprocessingPhenoStateDnaEpico()
#' timepoint_data <- splitTimepointsPreprocessingPheno(
#'     pheno = ex$pheno,
#'     metricsData = ex$metricsData,
#'     SampleID = "Sample_Name",
#'     timeVar = "Timepoint",
#'     timepoints = "1,2",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' timepoint_data$timepoints
#'
#' @description
#' Align phenotype rows and metric matrices for each requested timepoint, and
#' precompute the per-timepoint phenotype-plus-beta objects used by downstream
#' modeling functions.
#'
#' @export
splitTimepointsPreprocessingPheno <- function(
    pheno, metricsData,
    SampleID = "Sample_Name", timeVar = "Timepoint", timepoints = "1,2",
    methylationScale = "beta", verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_splitTimepointsPreprocessingPheno.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
    methylation_prefix <-
        methylationScaleObjectPrefixDnaEpico(methylation_scale)

    if (!(SampleID %in% colnames(pheno))) {
        stop("SampleID column not found in phenotype data: ",
            SampleID,
            call. = FALSE
        )
    }

    if (!(timeVar %in% colnames(pheno))) {
        stop("timeVar column not found in phenotype data: ",
            timeVar,
            call. = FALSE
        )
    }

    required_metrics <- c("beta", "m", "cn")
    if (!all(required_metrics %in% names(metricsData))) {
        stop("metricsData must contain beta, m, and cn matrices.",
            call. = FALSE
        )
    }

    requested_timepoints <- parseTimepointsPreprocessingPheno(
        values = timepoints,
        label = "timepoints"
    )
    available_timepoints <- table(pheno[[timeVar]], useNA = "ifany")
    subsets <- list()
    summary_lines <- c(paste("Requested timepoints:      ",
        paste(requested_timepoints,
        collapse = ", "
    )), paste(
        "Available values in", timeVar,
        "column:"
    ))
    summary_lines <- c(summary_lines,
        previewLinesMinfiEwasWater(available_timepoints))

    time_values <- as.character(pheno[[timeVar]])
    for (tp in requested_timepoints) {
        selected_timepoint <- !is.na(time_values) & time_values ==
            tp
        pheno_subset <- pheno[selected_timepoint, , drop = FALSE]

        if (nrow(pheno_subset) == 0L) {
            stop("No phenotype rows were found for timepoint: ",
                tp,
                call. = FALSE
            )
        }

        sample_ids <- validateSampleIdentifiersDnaEpico(
            pheno_subset[[SampleID]],
            paste0(
                "Phenotype sample identifiers at timepoint ",
                tp
            )
        )
        missing_metric_samples <- setdiff(sample_ids,
            colnames(metricsData$beta))
        if (length(missing_metric_samples) > 0L) {
            stop("Samples from timepoint ", tp,
                " are missing from the metric matrices: ",
                paste(utils::head(missing_metric_samples, 10L),
                    collapse = ", "
                ),
                call. = FALSE
            )
        }
        common_ids <- sample_ids

        pheno_subset <- pheno_subset[match(common_ids,
            as.character(pheno_subset[[SampleID]])), ,
            drop = FALSE
        ]

        beta_subset <- metricsData$beta[, common_ids, drop = FALSE]
        m_subset <- metricsData$m[, common_ids, drop = FALSE]
        cn_subset <- metricsData$cn[, common_ids, drop = FALSE]
        selected_subset <- metricsData[[methylation_scale]][,
            common_ids,
            drop = FALSE
        ]
        pheno_methylation <- mergePhenoMethylationPreprocessingPheno(
            phenoFrame = pheno_subset,
            methylationMatrix = selected_subset, id = SampleID,
            methylationScale = methylation_scale
        )

        subset_data <- list(
            pheno = pheno_subset, beta = beta_subset,
            m = m_subset, cn = cn_subset, phenoMethylation = pheno_methylation,
            methylationScale = methylation_scale,
                methylationObjectPrefix = methylation_prefix
        )
        subset_data[[methylation_prefix]] <- pheno_methylation
        subsets[[tp]] <- subset_data

        summary_lines <- c(
            summary_lines, paste("Timepoint ",
                tp, " rows:          ", nrow(pheno_subset),
                sep = ""
            ),
            paste("Timepoint ", tp, " matched samples:", length(common_ids),
                sep = " "
            ), paste("Timepoint ", tp, " merged object:  ",
                methylation_prefix, "T", tp,
                sep = ""
            )
        )
    }

    emitLogMinfiEwasWater(c(summary_lines,
        "============================================================"),
        verbose = verbose, log_path = log_path
    )

    structure(
        list(
            timepoints = requested_timepoints, data = subsets,
            methylationScale = methylation_scale,
                methylationObjectPrefix = methylation_prefix
        ),
        class = "dnaEPICO_preprocessingPheno_timepoints"
    )
}

#' Combine selected timepoints for downstream longitudinal modeling
#'
#' @param timepointData Object returned by
#' `splitTimepointsPreprocessingPheno()`.
#' @param combineTimepoints Character vector or comma-separated string of
#'   timepoints to combine.
#' @param methylationScale Character. Methylation metric used in the merged
#'   modeling table. One of `'Beta'`, `'M'`, or `'CN'`, case-insensitively.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_preprocessingPheno_combined'`
#'   containing the combined phenotype table, merged phenotype-plus-beta table,
#'   selected timepoints, and output suffix.
#'
#' @examples
#' ex <- dnaEPICO:::examplePreprocessingPhenoStateDnaEpico()
#' combined_data <- combineTimepointsPreprocessingPheno(
#'     timepointData = ex$timepointData,
#'     combineTimepoints = "1,2",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' combined_data$suffix
#'
#' @description
#' Combine selected timepoints that were already aligned by
#' `splitTimepointsPreprocessingPheno()` into the wide phenotype-plus-beta
#' objects used by downstream longitudinal models.
#'
#' @export
combineTimepointsPreprocessingPheno <- function(
    timepointData,
    combineTimepoints = "1,2", methylationScale = "beta", verbose = FALSE,
    logs = FALSE, log_dir = NULL,
        log_file = "log_combineTimepointsPreprocessingPheno.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    methylation_scale <- if (!is.null(timepointData$methylationScale)) {
        normalizeMethylationScaleDnaEpico(timepointData$methylationScale)
    } else {
        normalizeMethylationScaleDnaEpico(methylationScale)
    }
    methylation_prefix <-
        methylationScaleObjectPrefixDnaEpico(methylation_scale)

    requested_timepoints <- parseTimepointsPreprocessingPheno(
        values = combineTimepoints,
        label = "combineTimepoints"
    )

    available_timepoints <- names(timepointData$data)
    missing_timepoints <- setdiff(requested_timepoints, available_timepoints)

    if (length(missing_timepoints) > 0L) {
        stop("Requested combined timepoints are missing from timepointData: ",
            paste(missing_timepoints, collapse = ", "),
            call. = FALSE
        )
    }

    combined_pheno <- do.call(rbind, lapply(
        requested_timepoints,
        function(tp) timepointData$data[[tp]]$pheno
    ))
    combined_pheno_beta <- do.call(rbind, lapply(
        requested_timepoints,
        function(tp) {
            tp_data <- timepointData$data[[tp]]
            if (!is.null(tp_data$phenoMethylation)) {
                return(tp_data$phenoMethylation)
            }
            tp_data[[methylation_prefix]]
        }
    ))
    combine_suffix <- paste0("T", paste(requested_timepoints,
        collapse = "T"
    ))

    emitLogMinfiEwasWater(
        c(paste(
            "Combining timepoints:      ",
            paste(requested_timepoints, collapse = ", ")
        ), paste(
            "Combined phenotype rows:   ",
            nrow(combined_pheno)
        ), paste("Combined ", methylation_prefix,
            " rows:   ", nrow(combined_pheno_beta),
            sep = ""
        ), paste(
            "Combined suffix:           ",
            combine_suffix
        ), "============================================================"),
        verbose = verbose, log_path = log_path
    )

    combined_data <- list(
        timepoints = requested_timepoints,
        suffix = combine_suffix, pheno = combined_pheno,
            phenoMethylation = combined_pheno_beta,
        methylationScale = methylation_scale,
            methylationObjectPrefix = methylation_prefix
    )
    combined_data[[methylation_prefix]] <- combined_pheno_beta

    structure(combined_data, class = "dnaEPICO_preprocessingPheno_combined")
}

#' Build Clock Foundation input tables from preprocessingPheno data
#'
#' @param beta Numeric matrix of beta values with probes in rows and samples in
#'   columns.
#' @param pheno Phenotype data frame aligned with the beta matrix columns.
#' @param SampleID Character. Name of the phenotype sample identifier column.
#' @param sexColumn Character. Name of the phenotype sex column.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_preprocessingPheno_clock'`
#'   containing `betaCSV` and `phenoCF`.
#'
#' @examples
#' ex <- dnaEPICO:::examplePreprocessingPhenoStateDnaEpico()
#' clock_inputs <- buildClockFoundationInputsPreprocessingPheno(
#'     beta = ex$timepointData$data[["1"]]$beta,
#'     pheno = ex$timepointData$data[["1"]]$pheno,
#'     SampleID = "Sample_Name",
#'     sexColumn = "Sex",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(clock_inputs)
#'
#' @description
#' Prepare the beta and phenotype tables commonly exported for Clock Foundation
#' style downstream workflows, without writing them to disk.
#'
#' @export
buildClockFoundationInputsPreprocessingPheno <- function(
    beta,
    pheno, SampleID = "Sample_Name", sexColumn = "Sex", verbose = FALSE,
    logs = FALSE, log_dir = NULL,
        log_file = "log_buildClockFoundationInputsPreprocessingPheno.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    if (!(SampleID %in% colnames(pheno))) {
        stop("SampleID column not found in phenotype data: ",
            SampleID,
            call. = FALSE
        )
    }

    if (!(sexColumn %in% colnames(pheno))) {
        stop("sexColumn not found in phenotype data: ", sexColumn,
            call. = FALSE
        )
    }

    if ("id" %in% colnames(pheno) && !identical(SampleID, "id")) {
        stop("The phenotype data already contains a column named 'id'. ",
            "Rename that column or SampleID before building Clock Foundation inputs.",
            call. = FALSE
        )
    }

    if (!is.matrix(beta) || !is.numeric(beta) || is.null(rownames(beta)) ||
        is.null(colnames(beta))) {
        stop("beta must be a numeric matrix with probe and sample names.",
            call. = FALSE
        )
    }
    validateMethylationProbeIdentifiersDnaEpico(
        rownames(beta),
        "Clock Foundation beta row names"
    )
    beta_range <- summarizeMethylationRangeDnaEpico(beta, "beta")
    beta_sample_ids <- validateSampleIdentifiersDnaEpico(
        colnames(beta),
        "Beta-matrix sample identifiers"
    )
    pheno_sample_ids <- validateSampleIdentifiersDnaEpico(
        pheno[[SampleID]],
        paste0("Phenotype column '", SampleID, "'")
    )
    pheno_match <- matchSampleIdentifiersDnaEpico(
        query = beta_sample_ids,
        reference = pheno_sample_ids,
            queryLabel = "Beta-matrix sample identifiers",
        referenceLabel = paste0(
            "phenotype column '", SampleID,
            "'"
        ), requireSameSet = TRUE
    )
    pheno <- pheno[pheno_match, , drop = FALSE]
    rownames(pheno) <- NULL

    beta_csv <- as.data.frame(beta, stringsAsFactors = FALSE)
    beta_csv <- cbind(
        ProbeID = rownames(beta_csv), beta_csv,
        row.names = NULL, stringsAsFactors = FALSE
    )

    pheno_cf <- pheno
    colnames(pheno_cf)[colnames(pheno_cf) == SampleID] <- "id"

    sex_info <- canonicalizeSexDnaEpico(pheno_cf[[sexColumn]])
    if (anyNA(sex_info$code)) {
        stop("The Clock Foundation sex column contains missing or unsupported values",
            if (length(sex_info$unknown) > 0L) {
                paste0(": ", paste(sex_info$unknown, collapse = ", "))
            } else {
                "."
            },
            call. = FALSE
        )
    }
    pheno_cf[[sexColumn]] <- ifelse(sex_info$code == 0L, "Female",
        "Male"
    )
    sex_line <- "Sex values were standardized to Female and Male."

    preview_rows <- seq_len(min(nrow(beta_csv), 5L))
    preview_cols <- seq_len(min(ncol(beta_csv), 5L))

    emitLogMinfiEwasWater(
        c(
            paste(
                "Clock Foundation beta rows:",
                nrow(beta_csv)
            ), paste(
                "Clock Foundation beta cols:",
                ncol(beta_csv)
            ), paste(
                "Clock Foundation pheno rows:",
                nrow(pheno_cf)
            ), formatMethylationRangeLogDnaEpico(beta_range),
            sex_line, "Preview of Clock Foundation beta table:",
            previewLinesMinfiEwasWater(beta_csv[preview_rows, preview_cols,
                drop = FALSE
            ]), "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(
        list(
            betaCSV = beta_csv, phenoCF = pheno_cf,
            methylationRange = beta_range
        ),
        class = "dnaEPICO_preprocessingPheno_clock"
    )
}

#' Write preprocessingPheno outputs to disk
#'
#' @param preprocessingData Object returned by `preprocessingPheno()` or a list
#'   with the same components.
#' @param outputPheno Character. Directory used for saved phenotype CSV files.
#' @param outputRData Character. Directory used for saved metric `.RData` files.
#' @param outputRDataMerge Character. Directory used for saved merged
#'   phenotype-plus-beta `.RData` files.
#' @param outputDir Character. Directory used for the Clock Foundation export
#'   files.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_preprocessingPheno_paths'`
#'   containing the paths written to disk.
#'
#' @examples
#' ex <- dnaEPICO:::examplePreprocessingPhenoStateDnaEpico()
#' output_paths <- writePreprocessingPhenoOutputs(
#'     preprocessingData = ex$preprocessingData,
#'     outputPheno = file.path(ex$tempDir, "pheno"),
#'     outputRData = file.path(ex$tempDir, "metrics"),
#'     outputRDataMerge = file.path(ex$tempDir, "merge"),
#'     outputDir = file.path(ex$tempDir, "clock"),
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(output_paths)
#'
#' @description
#' Write the CSV, ZIP, and `.RData` outputs produced by
#' `preprocessingPheno()`. This helper keeps file writing separate from the
#' in-memory preprocessing steps.
#'
#' @export
writePreprocessingPhenoOutputs <- function(
    preprocessingData,
    outputPheno = "data/preprocessingPheno",
        outputRData = "rData/preprocessingPheno/metrics",
    outputRDataMerge = "rData/preprocessingPheno/mergeData",
    outputDir = "data/preprocessingPheno", verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_writePreprocessingPhenoOutputs.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    dir.create(outputPheno, recursive = TRUE, showWarnings = FALSE)
    dir.create(outputRData, recursive = TRUE, showWarnings = FALSE)
    dir.create(outputRDataMerge, recursive = TRUE, showWarnings = FALSE)
    dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

    methylation_scale <- if (!is.null(preprocessingData$methylationScale)) {
        normalizeMethylationScaleDnaEpico(preprocessingData$methylationScale)
    } else if (!is.null(preprocessingData$combinedData$methylationScale)) {
        normalizeMethylationScaleDnaEpico(preprocessingData$combinedData$methylationScale)
    } else {
        normalizeMethylationScaleDnaEpico("beta")
    }
    methylation_prefix <-
        methylationScaleObjectPrefixDnaEpico(methylation_scale)

    timepoint_paths <- list()

    for (tp in preprocessingData$timepointData$timepoints) {
        tp_data <- preprocessingData$timepointData$data[[tp]]

        pheno_path <- file.path(outputPheno, paste0(
            "phenoT",
            tp, ".csv"
        ))
        beta_path <- file.path(outputRData, paste0(
            "betaT", tp,
            ".RData"
        ))
        m_path <- file.path(outputRData, paste0("mT", tp, ".RData"))
        cn_path <- file.path(outputRData, paste0("cnT", tp, ".RData"))
        merged_path <- file.path(outputRDataMerge, paste0(
            methylation_prefix,
            "T", tp, ".RData"
        ))
        merged_object_name <- paste0(
            methylation_prefix, "T",
            tp
        )
        merged_data <- if (!is.null(tp_data$phenoMethylation)) {
            tp_data$phenoMethylation
        } else {
            tp_data[[methylation_prefix]]
        }

        utils::write.csv(tp_data$pheno, file = pheno_path, row.names = FALSE)
        saveNamedObjectMinfiEwasWater(tp_data$beta, paste0(
            "betaT",
            tp
        ), beta_path)
        saveNamedObjectMinfiEwasWater(tp_data$m, paste0(
            "mT",
            tp
        ), m_path)
        saveNamedObjectMinfiEwasWater(tp_data$cn, paste0(
            "cnT",
            tp
        ), cn_path)
        saveNamedObjectMinfiEwasWater(
            merged_data, merged_object_name,
            merged_path
        )

        tp_paths <- list(
            pheno = pheno_path, beta = beta_path,
            m = m_path, cn = cn_path, phenoMethylation = merged_path
        )
        tp_paths[[methylation_prefix]] <- merged_path
        timepoint_paths[[tp]] <- tp_paths
    }

    combined_pheno_path <- file.path(outputPheno, paste0(
        "pheno",
        preprocessingData$combinedData$suffix, ".csv"
    ))
    combined_merge_path <- file.path(outputRDataMerge, paste0(
        methylation_prefix,
        preprocessingData$combinedData$suffix, ".RData"
    ))
    combined_merge_object_name <- paste0(
        methylation_prefix,
        preprocessingData$combinedData$suffix
    )
    combined_merge_data <-
        if (!is.null(preprocessingData$combinedData$phenoMethylation)) {
        preprocessingData$combinedData$phenoMethylation
    } else {
        preprocessingData$combinedData[[methylation_prefix]]
    }
    beta_csv_path <- file.path(outputDir, "beta.csv")
    zip_path <- file.path(outputDir, "beta.zip")
    pheno_cf_path <- file.path(outputDir, "phenoCF.csv")

    utils::write.csv(preprocessingData$combinedData$pheno,
        file = combined_pheno_path,
        row.names = FALSE
    )
    saveNamedObjectMinfiEwasWater(
        combined_merge_data, combined_merge_object_name,
        combined_merge_path
    )
    utils::write.csv(preprocessingData$clockFoundation$betaCSV,
        file = beta_csv_path, row.names = FALSE
    )
    utils::write.csv(preprocessingData$clockFoundation$phenoCF,
        file = pheno_cf_path, row.names = FALSE
    )
    zip_status <- utils::zip(
        zipfile = zip_path, files = beta_csv_path,
        flags = "-j"
    )
    beta_zip_path <- if (identical(zip_status, 0L) && file.exists(zip_path)) {
        zip_path
    } else {
        NULL
    }

    emitLogMinfiEwasWater(
        c(
            paste(
                "Output phenotype dir:     ",
                outputPheno
            ), paste("RData metrics dir:        ", outputRData),
            paste("RData merge dir:          ", outputRDataMerge),
            paste("Clock Foundation dir:     ", outputDir), paste(
                "Saved combined phenotype: ",
                combined_pheno_path
            ), paste("Saved combined ", methylation_prefix,
                ": ", combined_merge_path,
                sep = ""
            ), paste(
                "Saved beta CSV:           ",
                beta_csv_path
            ), if (is.null(beta_zip_path)) {
                "Beta ZIP file was not created."
            } else {
                paste("Saved beta ZIP:           ", beta_zip_path)
            }, paste("Saved phenoCF:            ", pheno_cf_path),
            paste("ZIP status code:          ", zip_status),
                "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    output_paths <- list(
        timepoints = timepoint_paths, combinedPheno = combined_pheno_path,
        combinedPhenoMethylation = combined_merge_path,
            methylationScale = methylation_scale,
        methylationObjectPrefix = methylation_prefix, betaCSV = beta_csv_path,
        betaZIP = beta_zip_path, phenoCF = pheno_cf_path
    )
    combined_key <- switch(methylation_scale,
        beta = "combinedPhenoB",
        m = "combinedPhenoM",
        cn = "combinedPhenoCN"
    )
    output_paths[[combined_key]] <- combined_merge_path

    structure(output_paths, class = "dnaEPICO_preprocessingPheno_paths")
}
