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
    stop(label, " must contain at least one timepoint.", call. = FALSE)
  }

  parsed
}

normalizeMethylationScaleDnaEpico <- function(methylationScale = "beta") {
  if (length(methylationScale) != 1L || is.na(methylationScale)) {
    stop("methylationScale must be one of: beta, m, cn.", call. = FALSE)
  }

  methylation_scale <- tolower(trimws(as.character(methylationScale)))
  if (!(methylation_scale %in% c("beta", "m", "cn"))) {
    stop("methylationScale must be one of: beta, m, cn.", call. = FALSE)
  }

  methylation_scale
}

methylationScaleObjectPrefixDnaEpico <- function(methylationScale = "beta") {
  switch(
    normalizeMethylationScaleDnaEpico(methylationScale),
    beta = "phenoBeta",
    m = "phenoM",
    cn = "phenoCN"
  )
}

methylationScaleResponseLabelDnaEpico <- function(methylationScale = "beta") {
  switch(
    normalizeMethylationScaleDnaEpico(methylationScale),
    beta = "Beta values",
    m = "M-values",
    cn = "Copy number values"
  )
}

#' Extract a preferred named element from a loaded object
#'
#' @param object Loaded object.
#' @param preferred_name Character or `NULL`. Name of the preferred element to
#'   extract when `object` is a named list.
#'
#' @return The original object or its preferred named element.
#'
#' @description
#' Internal helper that unwraps legacy list-style saved objects used by
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

  if (is.list(object) && preferred_name %in% names(object)) {
    return(object[[preferred_name]])
  }

  object
}

#' Load a saved object used by preprocessingPheno helpers
#'
#' @param path Character. Path to an `.RData` or `.rds` file.
#' @param preferred_name Character or `NULL`. Preferred object name when `path`
#'   points to an `.RData` file containing multiple objects.
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
    return(
      extractPreferredObjectPreprocessingPheno(
        object = readRDS(path),
        preferred_name = preferred_name
      )
    )
  }

  load_env <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = load_env)

  if (!is.null(preferred_name) && preferred_name %in% loaded_names) {
    return(
      extractPreferredObjectPreprocessingPheno(
        object = load_env[[preferred_name]],
        preferred_name = preferred_name
      )
    )
  }

  if (length(loaded_names) == 1L) {
    return(
      extractPreferredObjectPreprocessingPheno(
        object = load_env[[loaded_names[[1L]]]],
        preferred_name = preferred_name
      )
    )
  }

  stop(
    "Could not determine which object to load from ",
    path,
    ". Available objects: ",
    paste(loaded_names, collapse = ", "),
    call. = FALSE
  )
}

#' Merge phenotype rows with methylation values for preprocessingPheno helpers
#'
#' @param phenoFrame Data frame containing phenotype rows.
#' @param methylationMatrix Numeric matrix with probes in rows and samples in
#'   columns.
#' @param id Character. Name of the sample identifier column in `phenoFrame`.
#' @param methylationScale Character. Selected methylation scale.
#'
#' @return A data frame with phenotype columns followed by transposed methylation
#'   values.
#'
#' @description
#' Internal helper that aligns phenotype rows with methylation-matrix columns before
#' building merged objects used by downstream modeling functions.
#'
#' @keywords internal
#' @noRd
mergePhenoMethylationPreprocessingPheno <- function(
    phenoFrame,
    methylationMatrix,
    id = "Sample_Name",
    methylationScale = "beta"
) {
  methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
  methylation_label <- methylationScaleResponseLabelDnaEpico(methylation_scale)

  if (!(id %in% colnames(phenoFrame))) {
    stop("Sample identifier column not found in phenotype data: ", id, call. = FALSE)
  }

  sample_ids <- as.character(phenoFrame[[id]])
  matched_ids <- sample_ids[sample_ids %in% colnames(methylationMatrix)]

  if (length(matched_ids) == 0L) {
    stop(
      "No matching sample identifiers were found between phenotype data and ",
      methylation_label,
      ".",
      call. = FALSE
    )
  }

  if (anyDuplicated(matched_ids) > 0L) {
    stop(
      "Duplicate sample identifiers were found while merging phenotype and methylation data: ",
      paste(unique(matched_ids[duplicated(matched_ids)]), collapse = ", "),
      call. = FALSE
    )
  }

  phenoFrame <- phenoFrame[match(matched_ids, sample_ids), , drop = FALSE]
  methylationMatrix <- methylationMatrix[, matched_ids, drop = FALSE]

  cbind(
    phenoFrame,
    as.data.frame(t(methylationMatrix), stringsAsFactors = FALSE)
  )
}

mergePhenoBetaPreprocessingPheno <- function(
    phenoFrame,
    betaMatrix,
    id = "Sample_Name"
) {
  mergePhenoMethylationPreprocessingPheno(
    phenoFrame = phenoFrame,
    methylationMatrix = betaMatrix,
    id = id,
    methylationScale = "beta"
  )
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
#' @return A list with class `"dnaEPICO_preprocessingPheno_metrics"`
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
loadMetricsPreprocessingPheno <- function(
    betaPath,
    mPath,
    cnPath,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_loadMetricsPreprocessingPheno.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  beta <- loadSavedObjectPreprocessingPheno(betaPath, preferred_name = "beta")
  m <- loadSavedObjectPreprocessingPheno(mPath, preferred_name = "m")
  cn <- loadSavedObjectPreprocessingPheno(cnPath, preferred_name = "cn")

  preview_rows <- seq_len(min(nrow(beta), 5L))
  preview_cols <- seq_len(min(ncol(beta), 5L))

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      paste("Beta path:                ", betaPath),
      paste("M-values path:            ", mPath),
      paste("CN path:                  ", cnPath),
      paste("Beta dimensions:          ", paste(dim(beta), collapse = " x ")),
      paste("M dimensions:             ", paste(dim(m), collapse = " x ")),
      paste("CN dimensions:            ", paste(dim(cn), collapse = " x ")),
      "Preview of beta values:",
      previewLinesMinfiEwasWater(beta[preview_rows, preview_cols, drop = FALSE]),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      beta = beta,
      m = m,
      cn = cn
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
#'   modeling table. One of `"beta"`, `"m"`, or `"cn"`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_preprocessingPheno_timepoints"`
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
    pheno,
    metricsData,
    SampleID = "Sample_Name",
    timeVar = "Timepoint",
    timepoints = "1,2",
    methylationScale = "beta",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_splitTimepointsPreprocessingPheno.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )
  methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
  methylation_prefix <- methylationScaleObjectPrefixDnaEpico(methylation_scale)
  methylation_label <- methylationScaleResponseLabelDnaEpico(methylation_scale)

  if (!(SampleID %in% colnames(pheno))) {
    stop("SampleID column not found in phenotype data: ", SampleID, call. = FALSE)
  }

  if (!(timeVar %in% colnames(pheno))) {
    stop("timeVar column not found in phenotype data: ", timeVar, call. = FALSE)
  }

  required_metrics <- c("beta", "m", "cn")
  if (!all(required_metrics %in% names(metricsData))) {
    stop(
      "metricsData must contain beta, m, and cn matrices.",
      call. = FALSE
    )
  }

  requested_timepoints <- parseTimepointsPreprocessingPheno(
    values = timepoints,
    label = "timepoints"
  )
  available_timepoints <- table(pheno[[timeVar]], useNA = "ifany")
  subsets <- list()
  summary_lines <- c(
    paste(
      "Requested timepoints:      ",
      paste(requested_timepoints, collapse = ", ")
    ),
    paste("Modeling methylation scale:", methylation_label),
    paste("Available values in", timeVar, "column:")
  )
  summary_lines <- c(
    summary_lines,
    previewLinesMinfiEwasWater(available_timepoints)
  )

  for (tp in requested_timepoints) {
    pheno_subset <- pheno[as.character(pheno[[timeVar]]) == tp, , drop = FALSE]

    if (nrow(pheno_subset) == 0L) {
      stop("No phenotype rows were found for timepoint: ", tp, call. = FALSE)
    }

    sample_ids <- as.character(pheno_subset[[SampleID]])
    common_ids <- sample_ids[
      sample_ids %in% colnames(metricsData$beta) &
      sample_ids %in% colnames(metricsData$m) &
      sample_ids %in% colnames(metricsData$cn)
    ]

    if (length(common_ids) == 0L) {
      stop(
        "No samples from timepoint ",
        tp,
        " were found across beta, m, and cn matrices.",
        call. = FALSE
      )
    }

    if (anyDuplicated(common_ids) > 0L) {
      stop(
        "Duplicate sample identifiers detected for timepoint ",
        tp,
        ": ",
        paste(unique(common_ids[duplicated(common_ids)]), collapse = ", "),
        call. = FALSE
      )
    }

    pheno_subset <- pheno_subset[
      match(common_ids, as.character(pheno_subset[[SampleID]])),
      ,
      drop = FALSE
    ]

    beta_subset <- metricsData$beta[, common_ids, drop = FALSE]
    m_subset <- metricsData$m[, common_ids, drop = FALSE]
    cn_subset <- metricsData$cn[, common_ids, drop = FALSE]
    selected_subset <- metricsData[[methylation_scale]][, common_ids, drop = FALSE]
    pheno_methylation <- mergePhenoMethylationPreprocessingPheno(
      phenoFrame = pheno_subset,
      methylationMatrix = selected_subset,
      id = SampleID,
      methylationScale = methylation_scale
    )

    subset_data <- list(
      pheno = pheno_subset,
      beta = beta_subset,
      m = m_subset,
      cn = cn_subset,
      phenoMethylation = pheno_methylation,
      methylationScale = methylation_scale,
      methylationObjectPrefix = methylation_prefix
    )
    subset_data[[methylation_prefix]] <- pheno_methylation
    subsets[[tp]] <- subset_data

    summary_lines <- c(
      summary_lines,
      paste("Timepoint ", tp, " rows:          ", nrow(pheno_subset), sep = ""),
      paste("Timepoint ", tp, " matched samples:", length(common_ids), sep = " "),
      paste("Timepoint ", tp, " merged object:  ", methylation_prefix, "T", tp, sep = "")
    )
  }

  emitLogMinfiEwasWater(
    c(summary_lines, "======================================================================="),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      timepoints = requested_timepoints,
      data = subsets,
      methylationScale = methylation_scale,
      methylationObjectPrefix = methylation_prefix
    ),
    class = "dnaEPICO_preprocessingPheno_timepoints"
  )
}

#' Combine selected timepoints for downstream longitudinal modeling
#'
#' @param timepointData Object returned by `splitTimepointsPreprocessingPheno()`.
#' @param combineTimepoints Character vector or comma-separated string of
#'   timepoints to combine.
#' @param methylationScale Character. Methylation metric used in the merged
#'   modeling table. One of `"beta"`, `"m"`, or `"cn"`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_preprocessingPheno_combined"`
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
combineTimepointsPreprocessingPheno <- function(
    timepointData,
    combineTimepoints = "1,2",
    methylationScale = "beta",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_combineTimepointsPreprocessingPheno.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )
  methylation_scale <- if (!is.null(timepointData$methylationScale)) {
    normalizeMethylationScaleDnaEpico(timepointData$methylationScale)
  } else {
    normalizeMethylationScaleDnaEpico(methylationScale)
  }
  methylation_prefix <- methylationScaleObjectPrefixDnaEpico(methylation_scale)
  methylation_label <- methylationScaleResponseLabelDnaEpico(methylation_scale)

  requested_timepoints <- parseTimepointsPreprocessingPheno(
    values = combineTimepoints,
    label = "combineTimepoints"
  )

  available_timepoints <- names(timepointData$data)
  missing_timepoints <- setdiff(requested_timepoints, available_timepoints)

  if (length(missing_timepoints) > 0L) {
    stop(
      "Requested combined timepoints are missing from timepointData: ",
      paste(missing_timepoints, collapse = ", "),
      call. = FALSE
    )
  }

  combined_pheno <- do.call(
    rbind,
    lapply(requested_timepoints, function(tp) timepointData$data[[tp]]$pheno)
  )
  combined_pheno_beta <- do.call(
    rbind,
    lapply(
      requested_timepoints,
      function(tp) {
        tp_data <- timepointData$data[[tp]]
        if (!is.null(tp_data$phenoMethylation)) {
          return(tp_data$phenoMethylation)
        }
        tp_data[[methylation_prefix]]
      }
    )
  )
  combine_suffix <- paste0("T", paste(requested_timepoints, collapse = "T"))

  emitLogMinfiEwasWater(
    c(
      paste(
        "Combining timepoints:      ",
        paste(requested_timepoints, collapse = ", ")
      ),
      paste("Combined phenotype rows:   ", nrow(combined_pheno)),
      paste("Modeling methylation scale:", methylation_label),
      paste("Combined ", methylation_prefix, " rows:   ", nrow(combined_pheno_beta), sep = ""),
      paste("Combined suffix:           ", combine_suffix),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  combined_data <- list(
    timepoints = requested_timepoints,
    suffix = combine_suffix,
    pheno = combined_pheno,
    phenoMethylation = combined_pheno_beta,
    methylationScale = methylation_scale,
    methylationObjectPrefix = methylation_prefix
  )
  combined_data[[methylation_prefix]] <- combined_pheno_beta

  structure(
    combined_data,
    class = "dnaEPICO_preprocessingPheno_combined"
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
#' @return A list with class `"dnaEPICO_preprocessingPheno_clock"`
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
    beta,
    pheno,
    SampleID = "Sample_Name",
    sexColumn = "Sex",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_buildClockFoundationInputsPreprocessingPheno.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  if (!(SampleID %in% colnames(pheno))) {
    stop("SampleID column not found in phenotype data: ", SampleID, call. = FALSE)
  }

  if (!(sexColumn %in% colnames(pheno))) {
    stop("sexColumn not found in phenotype data: ", sexColumn, call. = FALSE)
  }

  if ("id" %in% colnames(pheno) && !identical(SampleID, "id")) {
    stop(
      "The phenotype data already contains a column named 'id'. ",
      "Please rename either that column or SampleID before building Clock Foundation inputs.",
      call. = FALSE
    )
  }

  beta_csv <- as.data.frame(beta, stringsAsFactors = FALSE)
  beta_csv <- cbind(
    ProbeID = rownames(beta_csv),
    beta_csv,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  pheno_cf <- pheno
  colnames(pheno_cf)[colnames(pheno_cf) == SampleID] <- "id"

  unique_sex <- unique(stats::na.omit(as.character(pheno_cf[[sexColumn]])))
  sex_line <- "Sex column values were left unchanged."

  if (all(unique_sex %in% c("Male", "Female"))) {
    sex_line <- "Sex column already contains 'Male' and 'Female'. Skipping recoding."
  } else if (all(unique_sex %in% c("0", "1"))) {
    pheno_cf[[sexColumn]] <- ifelse(
      as.character(pheno_cf[[sexColumn]]) == "0",
      "Female",
      "Male"
    )
    sex_line <- "Re-encoding Sex: 0 = Female, 1 = Male"
  } else {
    pheno_cf[[sexColumn]] <- as.character(pheno_cf[[sexColumn]])
    sex_line <- paste(
      "Sex column values are not limited to 0/1 or Male/Female; keeping original values in",
      sexColumn
    )
  }

  preview_rows <- seq_len(min(nrow(beta_csv), 5L))
  preview_cols <- seq_len(min(ncol(beta_csv), 5L))

  emitLogMinfiEwasWater(
    c(
      paste("Clock Foundation beta rows:", nrow(beta_csv)),
      paste("Clock Foundation beta cols:", ncol(beta_csv)),
      paste("Clock Foundation pheno rows:", nrow(pheno_cf)),
      sex_line,
      "Preview of Clock Foundation beta table:",
      previewLinesMinfiEwasWater(beta_csv[preview_rows, preview_cols, drop = FALSE]),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      betaCSV = beta_csv,
      phenoCF = pheno_cf
    ),
    class = "dnaEPICO_preprocessingPheno_clock"
  )
}

#' Write legacy preprocessingPheno outputs to disk
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
#' @return A list with class `"dnaEPICO_preprocessingPheno_paths"`
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
#' Write the legacy CSV, ZIP, and `.RData` outputs produced by
#' `preprocessingPheno()`. This helper keeps file writing separate from the
#' in-memory preprocessing steps.
#'
#' @export
writePreprocessingPhenoOutputs <- function(
    preprocessingData,
    outputPheno = "data/preprocessingPheno",
    outputRData = "rData/preprocessingPheno/metrics",
    outputRDataMerge = "rData/preprocessingPheno/mergeData",
    outputDir = "data/preprocessingPheno",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_writePreprocessingPhenoOutputs.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
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
  methylation_prefix <- methylationScaleObjectPrefixDnaEpico(methylation_scale)

  timepoint_paths <- list()

  for (tp in preprocessingData$timepointData$timepoints) {
    tp_data <- preprocessingData$timepointData$data[[tp]]

    pheno_path <- file.path(outputPheno, paste0("phenoT", tp, ".csv"))
    beta_path <- file.path(outputRData, paste0("betaT", tp, ".RData"))
    m_path <- file.path(outputRData, paste0("mT", tp, ".RData"))
    cn_path <- file.path(outputRData, paste0("cnT", tp, ".RData"))
    merged_path <- file.path(outputRDataMerge, paste0(methylation_prefix, "T", tp, ".RData"))
    merged_object_name <- paste0(methylation_prefix, "T", tp)
    merged_data <- if (!is.null(tp_data$phenoMethylation)) {
      tp_data$phenoMethylation
    } else {
      tp_data[[methylation_prefix]]
    }

    utils::write.csv(tp_data$pheno, file = pheno_path, row.names = FALSE)
    saveNamedObjectMinfiEwasWater(tp_data$beta, paste0("betaT", tp), beta_path)
    saveNamedObjectMinfiEwasWater(tp_data$m, paste0("mT", tp), m_path)
    saveNamedObjectMinfiEwasWater(tp_data$cn, paste0("cnT", tp), cn_path)
    saveNamedObjectMinfiEwasWater(
      merged_data,
      merged_object_name,
      merged_path
    )

    tp_paths <- list(
      pheno = pheno_path,
      beta = beta_path,
      m = m_path,
      cn = cn_path,
      phenoMethylation = merged_path
    )
    tp_paths[[methylation_prefix]] <- merged_path
    timepoint_paths[[tp]] <- tp_paths
  }

  combined_pheno_path <- file.path(
    outputPheno,
    paste0("pheno", preprocessingData$combinedData$suffix, ".csv")
  )
  combined_merge_path <- file.path(
    outputRDataMerge,
    paste0(methylation_prefix, preprocessingData$combinedData$suffix, ".RData")
  )
  combined_merge_object_name <- paste0(methylation_prefix, preprocessingData$combinedData$suffix)
  combined_merge_data <- if (!is.null(preprocessingData$combinedData$phenoMethylation)) {
    preprocessingData$combinedData$phenoMethylation
  } else {
    preprocessingData$combinedData[[methylation_prefix]]
  }
  beta_csv_path <- file.path(outputDir, "beta.csv")
  zip_path <- file.path(outputDir, "beta.zip")
  pheno_cf_path <- file.path(outputDir, "phenoCF.csv")

  utils::write.csv(
    preprocessingData$combinedData$pheno,
    file = combined_pheno_path,
    row.names = FALSE
  )
  saveNamedObjectMinfiEwasWater(
    combined_merge_data,
    combined_merge_object_name,
    combined_merge_path
  )
  utils::write.csv(
    preprocessingData$clockFoundation$betaCSV,
    file = beta_csv_path,
    row.names = FALSE
  )
  utils::write.csv(
    preprocessingData$clockFoundation$phenoCF,
    file = pheno_cf_path,
    row.names = FALSE
  )
  zip_status <- utils::zip(zipfile = zip_path, files = beta_csv_path, flags = "-j")
  beta_zip_path <- if (identical(zip_status, 0L) && file.exists(zip_path)) {
    zip_path
  } else {
    NULL
  }

  emitLogMinfiEwasWater(
    c(
      paste("Output phenotype dir:     ", outputPheno),
      paste("RData metrics dir:        ", outputRData),
      paste("RData merge dir:          ", outputRDataMerge),
      paste("Clock Foundation dir:     ", outputDir),
      paste("Saved combined phenotype: ", combined_pheno_path),
      paste("Saved combined ", methylation_prefix, ": ", combined_merge_path, sep = ""),
      paste("Saved beta CSV:           ", beta_csv_path),
      if (is.null(beta_zip_path)) {
        "Beta ZIP file was not created."
      } else {
        paste("Saved beta ZIP:           ", beta_zip_path)
      },
      paste("Saved phenoCF:            ", pheno_cf_path),
      paste("ZIP status code:          ", zip_status),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  output_paths <- list(
    timepoints = timepoint_paths,
    combinedPheno = combined_pheno_path,
    combinedPhenoMethylation = combined_merge_path,
    methylationScale = methylation_scale,
    methylationObjectPrefix = methylation_prefix,
    betaCSV = beta_csv_path,
    betaZIP = beta_zip_path,
    phenoCF = pheno_cf_path
  )
  combined_key <- switch(
    methylation_scale,
    beta = "combinedPhenoBeta",
    m = "combinedPhenoM",
    cn = "combinedPhenoCN"
  )
  output_paths[[combined_key]] <- combined_merge_path

  structure(
    output_paths,
    class = "dnaEPICO_preprocessingPheno_paths"
  )
}
