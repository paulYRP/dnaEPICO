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

    loaded_names_text <- paste(loaded_names, collapse = ", ")
    stop(
    sprintf(
        "%s %s. Available objects: %s",
        "Could not determine which object to load from",
        path, loaded_names_text
    ),
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
mergePhenoMethylationPreprocessingPheno <- function(phenoFrame,
    methylationMatrix, id = "Sample_Name", methylationScale = "beta") {
    methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
    methylation_label <- methylationScaleResponseLabelDnaEpico(
        methylation_scale)
    if (!(id %in% colnames(phenoFrame))) {
        stop("Sample identifier column not found in phenotype data: ",
            id, call. = FALSE) }
    if (!is.matrix(methylationMatrix) || !is.numeric(methylationMatrix) ||
        is.null(rownames(methylationMatrix)) ||
        is.null(colnames(methylationMatrix))) {
        stop(methylation_label,
            " must be a numeric matrix with sample identifiers in column ",
            "names.", call. = FALSE)
    }
    validateMethylationProbeIdentifiersDnaEpico(rownames(methylationMatrix),
        paste0(methylation_label, " row names"))
    sample_ids <- validateSampleIdentifiersDnaEpico(phenoFrame[[id]],
        paste0("phenotype column '", id, "'"))
    matrix_ids <- validateSampleIdentifiersDnaEpico(colnames(methylationMatrix),
        paste0(methylation_label, " column names"))
    missing_samples <- setdiff(sample_ids, matrix_ids)
    if (length(missing_samples) > 0L) {
        missing_samples_text <- paste(utils::head(missing_samples,
            10L), collapse = ", ")
        stop(sprintf("Phenotype samples are missing from %s: %s",
            methylation_label, missing_samples_text),
            call. = FALSE) }
    matched_ids <- sample_ids
    if (length(matched_ids) == 0L) {
        stop(sprintf("%s %s %s.",
            "No matching sample identifiers were found between",
            "phenotype data and", methylation_label),
            call. = FALSE) }
    if (anyDuplicated(matched_ids) > 0L) {
        duplicated_ids_text <- paste(unique(matched_ids[duplicated(
            matched_ids)]),
            collapse = ", ")
        message_prefix <- paste(
            "Duplicate sample identifiers were found while merging",
            "phenotype and methylation data")
        stop(sprintf("%s: %s", message_prefix,
            duplicated_ids_text), call. = FALSE)
    }
    phenoFrame <- phenoFrame[match(matched_ids,
        sample_ids), , drop = FALSE]
    methylationMatrix <- methylationMatrix[,
        matched_ids, drop = FALSE]
    cbind(phenoFrame, as.data.frame(t(methylationMatrix),
        stringsAsFactors = FALSE)) }

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
#'   betaPath = ex$betaPath,
#'   mPath = ex$mPath,
#'   cnPath = ex$cnPath,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(metrics_data)
#'
#' @description
#' Load the metric matrices generated by `preprocessingMinfiEwasWater()` and
#' return them as a single in-memory object for downstream phenotype alignment.
#'
#' @export
loadMetricsPreprocessingPheno <- function(betaPath, mPath, cnPath,
    verbose = FALSE, logs = FALSE, log_dir = NULL, log_file =
        "log_loadMetricsPreprocessingPheno.txt") {
    log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir,
        log_file = log_file)
    beta <- loadSavedObjectPreprocessingPheno(betaPath, preferred_name = "beta")
    m <- loadSavedObjectPreprocessingPheno(mPath, preferred_name = "m")
    cn <- loadSavedObjectPreprocessingPheno(cnPath, preferred_name = "cn")
    metric_objects <- list(beta = beta, m = m, cn = cn)
    invalid_metrics <- names(metric_objects)[!vapply(metric_objects,
        function(x) is.matrix(x) && is.numeric(x), logical(1))]
    if (length(invalid_metrics) > 0L) {
        invalid_metrics_text <- paste(invalid_metrics, collapse = ", ")
        stop(sprintf("Metric objects must be numeric matrices: %s",
            invalid_metrics_text), call. = FALSE) }
    if (is.null(rownames(beta)) || is.null(colnames(beta))) {
        stop("The beta matrix must have probe and sample names.",
            call. = FALSE) }
    validateMethylationProbeIdentifiersDnaEpico(rownames(beta),
        "Beta-matrix row names")
    validateSampleIdentifiersDnaEpico(colnames(beta),
        "Beta-matrix sample identifiers"); for (metric_name in c("m", "cn")) {
        metric <- metric_objects[[metric_name]]
        if (!identical(dim(metric), dim(beta)) || !identical(rownames(metric),
            rownames(beta)) || !identical(colnames(metric), colnames(beta))) {
            stop(metric_name,
            " must have the same dimensions, probe order, and sample ",
                "order as beta.", call. = FALSE) } }
    range_summaries <- lapply(names(metric_objects), function(metric_name) {
        summarizeMethylationRangeDnaEpico(values = metric_objects[[
            metric_name]], methylationScale = metric_name) })
    names(range_summaries) <- names(metric_objects)
    preview_rows <- seq_len(min(nrow(beta), 5L))
    preview_cols <- seq_len(min(ncol(beta), 5L)); emitLogMinfiEwasWater(c(
        "============================================================",
        paste("Beta path:                ", betaPath), paste(
            "M-values path:            ",
            mPath), paste("CN path:                  ", cnPath),
        paste("Beta dimensions:          ", paste(dim(beta), collapse = " x ")),
        paste("M dimensions:             ", paste(dim(m), collapse = " x ")),
        paste("CN dimensions:            ", paste(dim(cn), collapse = " x ")),
        unlist(lapply(range_summaries, function(range_summary) {
            formatMethylationRangeLogDnaEpico(range_summary)
        }), use.names = FALSE), "Preview of beta values:",
            previewLinesMinfiEwasWater(beta[preview_rows,
            preview_cols, drop = FALSE]),
            "============================================================"),
        verbose = verbose, log_path = log_path)
    structure(list(beta = beta, m = m, cn = cn, methylationRanges =
        range_summaries), class = "dnaEPICO_preprocessingPheno_metrics") }

validateTimepointSplitInputsPreprocessingPheno <- function(
    pheno, metricsData, SampleID, timeVar
) {
    if (!(SampleID %in% colnames(pheno))) {
    stop("SampleID column not found in phenotype data: ", SampleID,
        call. = FALSE
    )
    }
    if (!(timeVar %in% colnames(pheno))) {
    stop("timeVar column not found in phenotype data: ", timeVar,
        call. = FALSE
    )
    }
    if (!all(c("beta", "m", "cn") %in% names(metricsData))) {
    stop("metricsData must contain beta, m, and cn matrices.",
        call. = FALSE
    )
    }
    invisible(NULL)
}

oneTimepointPreprocessingPheno <- function(
    pheno, metricsData, SampleID, timeVar, timepoint, scale, prefix
) {
    selected <- !is.na(pheno[[timeVar]]) &
    as.character(pheno[[timeVar]]) == timepoint
    pheno_subset <- pheno[selected, , drop = FALSE]
    if (!nrow(pheno_subset)) {
    stop("No phenotype rows were found for timepoint: ", timepoint,
        call. = FALSE
    )
    }
    sample_ids <- validateSampleIdentifiersDnaEpico(
    pheno_subset[[SampleID]],
    paste0("Phenotype sample identifiers at timepoint ", timepoint)
    )
    missing_samples <- setdiff(sample_ids, colnames(metricsData$beta))
    if (length(missing_samples)) {
    missing_text <- paste(utils::head(missing_samples, 10L), collapse = ", ")
    stop(sprintf(
        "%s %s are missing from the metric matrices: %s",
        "Samples from timepoint", timepoint, missing_text
    ), call. = FALSE)
    }
    pheno_subset <- pheno_subset[
    match(sample_ids, as.character(pheno_subset[[SampleID]])), ,
    drop = FALSE
    ]
    selected_metric <- metricsData[[scale]][, sample_ids, drop = FALSE]
    merged <- mergePhenoMethylationPreprocessingPheno(
    phenoFrame = pheno_subset, methylationMatrix = selected_metric,
    id = SampleID, methylationScale = scale
    )
    result <- list(
    pheno = pheno_subset,
    beta = metricsData$beta[, sample_ids, drop = FALSE],
    m = metricsData$m[, sample_ids, drop = FALSE],
    cn = metricsData$cn[, sample_ids, drop = FALSE],
    phenoMethylation = merged, methylationScale = scale,
    methylationObjectPrefix = prefix
    )
    result[[prefix]] <- merged
    result
}

timepointSplitLogLinesPreprocessingPheno <- function(
    requested, pheno, timeVar, subsets, prefix
) {
    lines <- c(
    paste("Requested timepoints:      ", paste(requested, collapse = ", ")),
    paste("Available values in", timeVar, "column:"),
    previewLinesMinfiEwasWater(table(pheno[[timeVar]], useNA = "ifany"))
    )
    for (timepoint in requested) {
    subset <- subsets[[timepoint]]
    lines <- c(
        lines,
        paste0(
        "Timepoint ", timepoint, " rows:          ",
        nrow(subset$pheno)
        ),
        paste("Timepoint ", timepoint, " matched samples:",
        nrow(subset$pheno),
        sep = " "
        ),
        paste0(
        "Timepoint ", timepoint, " merged object:  ",
        prefix, "T", timepoint
        )
    )
    }
    c(lines, "============================================================")
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
#'   pheno = ex$pheno,
#'   metricsData = ex$metricsData,
#'   SampleID = "Sample_Name",
#'   timeVar = "Timepoint",
#'   timepoints = "1,2",
#'   verbose = FALSE,
#'   logs = FALSE
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
    pheno, metricsData, SampleID = "Sample_Name", timeVar = "Timepoint",
    timepoints = "1,2", methylationScale = "beta", verbose = FALSE,
    logs = FALSE, log_dir = NULL,
    log_file = "log_splitTimepointsPreprocessingPheno.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    scale <- normalizeMethylationScaleDnaEpico(methylationScale)
    prefix <- methylationScaleObjectPrefixDnaEpico(scale)
    validateTimepointSplitInputsPreprocessingPheno(
    pheno, metricsData, SampleID, timeVar
    )
    requested <- parseTimepointsPreprocessingPheno(
    values = timepoints, label = "timepoints"
    )
    subsets <- stats::setNames(vector("list", length(requested)), requested)
    for (timepoint in requested) {
    subsets[[timepoint]] <- oneTimepointPreprocessingPheno(
        pheno, metricsData, SampleID, timeVar, timepoint, scale, prefix
    )
    }
    emitLogMinfiEwasWater(
    timepointSplitLogLinesPreprocessingPheno(
        requested, pheno, timeVar, subsets, prefix
    ),
    verbose = verbose, log_path = log_path
    )
    structure(list(
    timepoints = requested, data = subsets, methylationScale = scale,
    methylationObjectPrefix = prefix
    ), class = "dnaEPICO_preprocessingPheno_timepoints")
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
#'   timepointData = ex$timepointData,
#'   combineTimepoints = "1,2",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' combined_data$suffix
#'
#' @description
#' Combine selected timepoints that were already aligned by
#' `splitTimepointsPreprocessingPheno()` into the wide phenotype-plus-beta
#' objects used by downstream longitudinal models.
#'
#' @export
combineTimepointsPreprocessingPheno <- function(timepointData,
    combineTimepoints = "1,2", methylationScale = "beta",
    verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_combineTimepointsPreprocessingPheno.txt") {
    log_path <- resolveLogPathMinfiEwasWater(logs = logs,
        log_dir = log_dir, log_file = log_file)
    methylation_scale <- if (!is.null(timepointData$methylationScale)) {
        normalizeMethylationScaleDnaEpico(timepointData$methylationScale)
    } else {
        normalizeMethylationScaleDnaEpico(methylationScale)
    }
    methylation_prefix <- methylationScaleObjectPrefixDnaEpico(
        methylation_scale)
    requested_timepoints <- parseTimepointsPreprocessingPheno(values =
        combineTimepoints, label = "combineTimepoints")
    available_timepoints <- names(timepointData$data)
    missing_timepoints <- setdiff(requested_timepoints,
        available_timepoints)
    if (length(missing_timepoints) > 0L) {
        missing_timepoints_text <- paste(missing_timepoints,
            collapse = ", "); stop(sprintf(
            "Requested combined timepoints are missing from timepointData: %s",
            missing_timepoints_text), call. = FALSE) }
    combined_pheno <- do.call(rbind, lapply(requested_timepoints,
        function(tp) timepointData$data[[tp]]$pheno))
    combined_pheno_beta <- do.call(rbind, lapply(requested_timepoints,
        function(tp) {
            tp_data <- timepointData$data[[tp]]
            if (!is.null(tp_data$phenoMethylation)) {
                return(tp_data$phenoMethylation)
            }; tp_data[[methylation_prefix]] }))
    combine_suffix <- paste0("T", paste(requested_timepoints,
        collapse = "T"))
    emitLogMinfiEwasWater(c(paste("Combining timepoints:      ",
        paste(requested_timepoints, collapse = ", ")),
        paste("Combined phenotype rows:   ", nrow(combined_pheno)),
        paste("Combined ", methylation_prefix, " rows:   ",
            nrow(combined_pheno_beta), sep = ""), paste(
            "Combined suffix:           ",
            combine_suffix),
            "============================================================"),
        verbose = verbose, log_path = log_path)
    combined_data <- list(timepoints = requested_timepoints,
        suffix = combine_suffix, pheno = combined_pheno,
        phenoMethylation = combined_pheno_beta, methylationScale =
            methylation_scale,
        methylationObjectPrefix = methylation_prefix)
    combined_data[[methylation_prefix]] <- combined_pheno_beta
    structure(combined_data, class = "dnaEPICO_preprocessingPheno_combined")
}

validateClockInputsPreprocessingPheno <- function(
    beta, pheno, SampleID, sexColumn
) {
    if (!(SampleID %in% colnames(pheno))) {
    stop("SampleID column not found in phenotype data: ", SampleID,
        call. = FALSE
    )
    }
    if (!(sexColumn %in% colnames(pheno))) {
    stop("sexColumn not found in phenotype data: ", sexColumn,
        call. = FALSE
    )
    }
    if ("id" %in% colnames(pheno) && !identical(SampleID, "id")) {
    stop(
        "The phenotype data already contains a column named 'id'. ",
        "Rename that column or SampleID before building Clock Foundation ",
        "inputs.",
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
    rownames(beta), "Clock Foundation beta row names"
    )
    invisible(NULL)
}

alignClockPhenotypePreprocessingPheno <- function(beta, pheno, SampleID) {
    beta_ids <- validateSampleIdentifiersDnaEpico(
    colnames(beta), "Beta-matrix sample identifiers"
    )
    pheno_ids <- validateSampleIdentifiersDnaEpico(
    pheno[[SampleID]], paste0("Phenotype column '", SampleID, "'")
    )
    matched <- matchSampleIdentifiersDnaEpico(
    query = beta_ids, reference = pheno_ids,
    queryLabel = "Beta-matrix sample identifiers",
    referenceLabel = paste0("phenotype column '", SampleID, "'"),
    requireSameSet = TRUE
    )
    pheno <- pheno[matched, , drop = FALSE]
    rownames(pheno) <- NULL
    pheno
}

clockPhenotypePreprocessingPheno <- function(pheno, SampleID, sexColumn) {
    colnames(pheno)[colnames(pheno) == SampleID] <- "id"
    sex_info <- canonicalizeSexDnaEpico(pheno[[sexColumn]])
    if (anyNA(sex_info$code)) {
    if (length(sex_info$unknown)) {
        unknown_text <- paste(sex_info$unknown, collapse = ", ")
        message_template <- paste0(
            "The Clock Foundation sex column contains missing or ",
            "unsupported values: %s"
        )
        stop(sprintf(
        message_template, unknown_text
        ), call. = FALSE)
    }
    stop("The Clock Foundation sex column contains missing or ",
        "unsupported values.",
        call. = FALSE
    )
    }
    pheno[[sexColumn]] <- ifelse(sex_info$code == 0L, "Female", "Male")
    pheno
}

clockInputLogLinesPreprocessingPheno <- function(betaCSV, phenoCF, range) {
    preview_rows <- seq_len(min(nrow(betaCSV), 5L))
    preview_cols <- seq_len(min(ncol(betaCSV), 5L))
    c(
    paste("Clock Foundation beta rows:", nrow(betaCSV)),
    paste("Clock Foundation beta cols:", ncol(betaCSV)),
    paste("Clock Foundation pheno rows:", nrow(phenoCF)),
    formatMethylationRangeLogDnaEpico(range),
    "Sex values were standardized to Female and Male.",
    "Preview of Clock Foundation beta table:",
    previewLinesMinfiEwasWater(
        betaCSV[preview_rows, preview_cols, drop = FALSE]
    ),
    "============================================================"
    )
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
#'   beta = ex$timepointData$data[["1"]]$beta,
#'   pheno = ex$timepointData$data[["1"]]$pheno,
#'   SampleID = "Sample_Name",
#'   sexColumn = "Sex",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(clock_inputs)
#'
#' @description
#' Prepare the beta and phenotype tables commonly exported for Clock Foundation
#' style downstream workflows, without writing them to disk.
#'
#' @export
buildClockFoundationInputsPreprocessingPheno <- function(
    beta, pheno, SampleID = "Sample_Name", sexColumn = "Sex",
    verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_buildClockFoundationInputsPreprocessingPheno.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    validateClockInputsPreprocessingPheno(beta, pheno, SampleID, sexColumn)
    methylation_range <- summarizeMethylationRangeDnaEpico(beta, "beta")
    pheno <- alignClockPhenotypePreprocessingPheno(beta, pheno, SampleID)
    beta_csv <- as.data.frame(beta, stringsAsFactors = FALSE)
    beta_csv <- cbind(
    ProbeID = rownames(beta_csv), beta_csv,
    row.names = NULL, stringsAsFactors = FALSE
    )
    pheno_cf <- clockPhenotypePreprocessingPheno(pheno, SampleID, sexColumn)
    emitLogMinfiEwasWater(
    clockInputLogLinesPreprocessingPheno(
        beta_csv, pheno_cf, methylation_range
    ),
    verbose = verbose, log_path = log_path
    )
    structure(list(
    betaCSV = beta_csv, phenoCF = pheno_cf,
    methylationRange = methylation_range
    ), class = "dnaEPICO_preprocessingPheno_clock")
}

outputMethylationScalePreprocessingPheno <- function(data) {
    if (!is.null(data$methylationScale)) {
    return(normalizeMethylationScaleDnaEpico(data$methylationScale))
    }
    if (!is.null(data$combinedData$methylationScale)) {
    return(normalizeMethylationScaleDnaEpico(
        data$combinedData$methylationScale
    ))
    }
    normalizeMethylationScaleDnaEpico("beta")
}

writeTimepointOutputPreprocessingPheno <- function(
    data, timepoint, prefix, outputPheno, outputRData, outputRDataMerge
) {
    pheno_path <- file.path(outputPheno, paste0("phenoT", timepoint, ".csv"))
    paths <- list(
    pheno = pheno_path,
    beta = file.path(outputRData, paste0("betaT", timepoint, ".RData")),
    m = file.path(outputRData, paste0("mT", timepoint, ".RData")),
    cn = file.path(outputRData, paste0("cnT", timepoint, ".RData")),
    phenoMethylation = file.path(
        outputRDataMerge, paste0(prefix, "T", timepoint, ".RData")
    )
    )
    merged <- if (!is.null(data$phenoMethylation)) {
    data$phenoMethylation
    } else {
    data[[prefix]]
    }
    utils::write.csv(data$pheno, file = paths$pheno, row.names = FALSE)
    for (metric in c("beta", "m", "cn")) {
    saveNamedObjectMinfiEwasWater(
        data[[metric]], paste0(metric, "T", timepoint), paths[[metric]]
    )
    }
    saveNamedObjectMinfiEwasWater(
    merged, paste0(prefix, "T", timepoint), paths$phenoMethylation
    )
    paths[[prefix]] <- paths$phenoMethylation
    paths
}

writeCombinedOutputPreprocessingPheno <- function(
    combined, prefix, outputPheno, outputRDataMerge
) {
    if (is.null(combined)) {
    return(list(pheno = NULL, merged = NULL))
    }
    pheno_path <- file.path(
    outputPheno, paste0("pheno", combined$suffix, ".csv")
    )
    merged_path <- file.path(
    outputRDataMerge, paste0(prefix, combined$suffix, ".RData")
    )
    merged <- if (!is.null(combined$phenoMethylation)) {
    combined$phenoMethylation
    } else {
    combined[[prefix]]
    }
    utils::write.csv(combined$pheno, file = pheno_path, row.names = FALSE)
    saveNamedObjectMinfiEwasWater(
    merged, paste0(prefix, combined$suffix), merged_path
    )
    list(pheno = pheno_path, merged = merged_path)
}

writeClockOutputPreprocessingPheno <- function(clock, outputDir) {
    beta_path <- file.path(outputDir, "beta.csv")
    zip_path <- file.path(outputDir, "beta.zip")
    pheno_path <- file.path(outputDir, "phenoCF.csv")
    utils::write.csv(clock$betaCSV, file = beta_path, row.names = FALSE)
    utils::write.csv(clock$phenoCF, file = pheno_path, row.names = FALSE)
    zip_status <- utils::zip(zipfile = zip_path, files = beta_path, flags =
        "-j")
    list(
    beta = beta_path,
    zip = if (identical(zip_status, 0L) && file.exists(zip_path)) {
        zip_path
    } else {
        NULL
    },
    pheno = pheno_path, zipStatus = zip_status
    )
}

outputLogLinesPreprocessingPheno <- function(
    outputPheno, outputRData, outputRDataMerge, outputDir,
    combined, clock, prefix
) {
    c(
    paste("Output phenotype dir:     ", outputPheno),
    paste("RData metrics dir:        ", outputRData),
    paste("RData merge dir:          ", outputRDataMerge),
    paste("Clock Foundation dir:     ", outputDir),
    if (is.null(combined$pheno)) {
        "Saved combined phenotype: disabled"
    } else {
        paste("Saved combined phenotype: ", combined$pheno)
    },
    if (is.null(combined$merged)) {
        paste0("Saved combined ", prefix, ": disabled")
    } else {
        paste0("Saved combined ", prefix, ": ", combined$merged)
    },
    paste("Saved beta CSV:           ", clock$beta),
    if (is.null(clock$zip)) {
        "Beta ZIP file was not created."
    } else {
        paste("Saved beta ZIP:           ", clock$zip)
    },
    paste("Saved phenoCF:            ", clock$pheno),
    paste("ZIP status code:          ", clock$zipStatus),
    "============================================================"
    )
}

asOutputPathsPreprocessingPheno <- function(
    timepoints, combined, clock, scale, prefix
) {
    output <- list(
    timepoints = timepoints, combinedPheno = combined$pheno,
    combinedPhenoMethylation = combined$merged,
    methylationScale = scale, methylationObjectPrefix = prefix,
    betaCSV = clock$beta, betaZIP = clock$zip, phenoCF = clock$pheno
    )
    key <- switch(scale,
    beta = "combinedPhenoB",
    m = "combinedPhenoM",
    cn = "combinedPhenoCN"
    )
    output[key] <- list(combined$merged)
    structure(output, class = "dnaEPICO_preprocessingPheno_paths")
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
#'   preprocessingData = ex$preprocessingData,
#'   outputPheno = file.path(ex$tempDir, "pheno"),
#'   outputRData = file.path(ex$tempDir, "metrics"),
#'   outputRDataMerge = file.path(ex$tempDir, "merge"),
#'   outputDir = file.path(ex$tempDir, "clock"),
#'   verbose = FALSE,
#'   logs = FALSE
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
    preprocessingData, outputPheno = "data/preprocessingPheno",
    outputRData = "rData/preprocessingPheno/metrics",
    outputRDataMerge = "rData/preprocessingPheno/mergeData",
    outputDir = "data/preprocessingPheno", verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_writePreprocessingPhenoOutputs.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    directories <- c(outputPheno, outputRData, outputRDataMerge, outputDir)
    for (directory in directories) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    }
    scale <- outputMethylationScalePreprocessingPheno(preprocessingData)
    prefix <- methylationScaleObjectPrefixDnaEpico(scale)
    timepoint_paths <- list()
    for (timepoint in preprocessingData$timepointData$timepoints) {
    timepoint_paths[[timepoint]] <- writeTimepointOutputPreprocessingPheno(
        preprocessingData$timepointData$data[[timepoint]], timepoint,
        prefix, outputPheno, outputRData, outputRDataMerge
    )
    }
    combined <- writeCombinedOutputPreprocessingPheno(
    preprocessingData$combinedData, prefix, outputPheno, outputRDataMerge
    )
    clock <- writeClockOutputPreprocessingPheno(
    preprocessingData$clockFoundation, outputDir
    )
    emitLogMinfiEwasWater(
    outputLogLinesPreprocessingPheno(
        outputPheno, outputRData, outputRDataMerge, outputDir,
        combined, clock, prefix
    ),
    verbose = verbose, log_path = log_path
    )
    asOutputPathsPreprocessingPheno(
    timepoint_paths, combined, clock, scale, prefix
    )
}
