logPreprocessingPhenoStartDnaEpico <- function(config) {
    separator <- resolveSeparatorMinfiEwasWater(config$sepType)
    emitLogMinfiEwasWater(c(
        "==== Starting preprocessingPheno ====",
        paste("Start Time:               ", format(Sys.time())),
        paste("Log file path:            ",
            if (is.null(config$log_path)) "disabled" else config$log_path
        ),
        paste("Phenotype file:           ", config$phenoFile),
        paste("Separator type:           ",
            if (is.null(separator)) "default (',')" else config$sepType
        ),
        paste("Beta path:                ", config$betaPath),
        paste("M-values path:            ", config$mPath),
        paste("CN path:                  ", config$cnPath),
        paste("Identifier column:        ", config$SampleID),
        paste("Timepoint column:         ", config$timeVar),
        paste("Timepoints:               ", config$timepoints),
        paste("Combine timepoints:       ",
            if (is.null(config$combineTimepoints)) "disabled" else {
                paste(config$combineTimepoints, collapse = ", ")
            }
        ),
        paste("Merged modeling object:   ",
            config$methylationObjectPrefix, "*"
        ),
        "Clock Foundation scale:     Beta values",
        paste("Sex column:               ", config$sexColumn),
        paste("Output phenotype dir:     ", config$outputPheno),
        paste("RData metrics dir:        ", config$outputRData),
        paste("RData merge dir:          ", config$outputRDataMerge),
        paste("Clock Foundation dir:     ", config$outputDir),
        paste("Save outputs:             ", config$saveOutputs),
        "============================================================"
    ), verbose = config$verbose, log_path = config$log_path)
}

preparePreprocessingPhenoStagesDnaEpico <- function(config) {
    pheno <- readPhenotypeTargets(
        phenoFile = config$phenoFile, sepType = config$sepType,
        SampleID = config$SampleID, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$log_file
    )
    metrics <- loadMetricsPreprocessingPheno(
        betaPath = config$betaPath, mPath = config$mPath,
        cnPath = config$cnPath, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$log_file
    )
    timepoints <- splitTimepointsPreprocessingPheno(
        pheno = pheno, metricsData = metrics, SampleID = config$SampleID,
        timeVar = config$timeVar, timepoints = config$timepoints,
        methylationScale = config$methylationScale,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
    combined <- if (is.null(config$combineTimepoints)) NULL else {
        combineTimepointsPreprocessingPheno(
            timepointData = timepoints,
            combineTimepoints = config$combineTimepoints,
            methylationScale = config$methylationScale,
            verbose = config$verbose, logs = config$logs,
            log_dir = config$outputLogs, log_file = config$log_file
        )
    }
    list(
        pheno = pheno, metricsData = metrics,
        timepointData = timepoints, combinedData = combined
    )
}

buildPreprocessingPhenoResultDnaEpico <- function(stages, config) {
    result <- c(stages, list(
        methylationScale = config$methylationScale,
        methylationLabel = config$methylationLabel,
        methylationObjectPrefix = config$methylationObjectPrefix,
        clockFoundation = buildClockFoundationInputsPreprocessingPheno(
            beta = stages$metricsData$beta, pheno = stages$pheno,
            SampleID = config$SampleID, sexColumn = config$sexColumn,
            verbose = config$verbose, logs = config$logs,
            log_dir = config$outputLogs, log_file = config$log_file
        ),
        savedFiles = NULL, logFile = config$log_path
    ))
    if (isTRUE(config$saveOutputs)) {
        result$savedFiles <- writePreprocessingPhenoOutputs(
            preprocessingData = result, outputPheno = config$outputPheno,
            outputRData = config$outputRData,
            outputRDataMerge = config$outputRDataMerge,
            outputDir = config$outputDir, verbose = config$verbose,
            logs = config$logs, log_dir = config$outputLogs,
            log_file = config$log_file
        )
    }
    result
}

runPreprocessingPhenoDnaEpico <- function(config) {
    stages <- preparePreprocessingPhenoStagesDnaEpico(config)
    result <- buildPreprocessingPhenoResultDnaEpico(stages, config)
    emitLogMinfiEwasWater(c(
        "==== Finished preprocessingPheno ====",
        paste("End Time:                 ", format(Sys.time())),
        "============================================================"
    ), verbose = config$verbose, log_path = config$log_path)
    structure(result, class = "dnaEPICO_preprocessingPheno")
}

#' Prepare phenotype and methylation matrices for downstream modeling
#'
#' Align the phenotype table with preprocessed beta, M-value, and copy-number
#' matrices, split the data by timepoint, optionally prepare longitudinal
#' objects for the selected modeling scale, and build Clock Foundation export
#' tables. The function writes files only when `saveOutputs = TRUE`.
#'
#' @param phenoFile Character. Path to the phenotype CSV file.
#' @param sepType Character or `NULL`. Field separator used in `phenoFile`. Use
#'   `NULL` for a comma-separated file, `'\\t'` for a tab-delimited file, or
#'   another separator accepted by `utils::read.csv()`.
#' @param betaPath Character. Path to the saved beta-value object. Both `.RData`
#'   and `.rds` files are supported.
#' @param mPath Character. Path to the saved M-value object. Both `.RData` and
#'   `.rds` files are supported.
#' @param cnPath Character. Path to the saved copy-number object. Both `.RData`
#'   and `.rds` files are supported.
#' @param SampleID Character. Name of the phenotype column containing sample
#'   identifiers used to align phenotype and methylation data.
#' @param timeVar Character. Name of the phenotype column containing timepoint
#'   labels.
#' @param timepoints Character vector or comma-separated string of timepoints to
#'   retain and split into separate in-memory subsets.
#' @param combineTimepoints Character vector or comma-separated string of
#'   timepoints to combine into the longitudinal phenotype-plus-methylation
#'   object, or `NULL` to skip the combined object.
#' @param methylationScale Character. Methylation metric to use in merged
#'   modeling tables. One of `'Beta'`, `'M'`, or `'CN'`, in any combination
#'   of upper- and lower-case letters. The default is `'beta'`. Beta values
#'   are always used for Clock Foundation exports.
#' @param outputPheno Character. Directory used for saved phenotype CSV files
#'   when `saveOutputs = TRUE`.
#' @param outputRData Character. Directory used for saved metric `.RData` files
#'   when `saveOutputs = TRUE`.
#' @param outputRDataMerge Character. Directory used for saved merged
#'   phenotype-plus-methylation `.RData` files when `saveOutputs = TRUE`.
#' @param sexColumn Character. Name of the phenotype sex column used when
#'   building Clock Foundation exports. Values are preserved as supplied,
#'   including missing, blank, unknown, or other character values; `PredSex` is
#'   not substituted.
#' @param outputLogs Character. Directory used for log files when `logs = TRUE`.
#' @param outputDir Character. Directory used for Clock Foundation export files
#'   when `saveOutputs = TRUE`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#'   The default is `FALSE`.
#' @param logs Logical. If `TRUE`, write the same progress messages to
#'   `outputLogs`. The default is `FALSE`.
#' @param saveOutputs Logical. If `TRUE`, write the CSV, ZIP, and `.RData`
#'   outputs to disk. The default is `FALSE`.
#'
#' @return A list with class `'dnaEPICO_preprocessingPheno'`.
#' \describe{
#'   \item{pheno}{Phenotype table read from `phenoFile`.}
#'   \item{metricsData}{Object returned by [loadMetricsPreprocessingPheno()]
#'   containing the beta-value, M-value, and copy-number matrices loaded from
#'   `betaPath`, `mPath`, and `cnPath`.}
#' \item{timepointData}{Object returned by
#' [splitTimepointsPreprocessingPheno()]
#'   containing per-timepoint phenotype tables and methylation matrices.}
#'   \item{combinedData}{Object returned by
#'   [combineTimepointsPreprocessingPheno()] containing the merged longitudinal
#'   phenotype-plus-methylation object and the timepoint combination metadata,
#'   or `NULL` when `combineTimepoints = NULL`.}
#'   \item{clockFoundation}{Object returned by
#'   [buildClockFoundationInputsPreprocessingPheno()] containing the beta table
#'   and phenotype table prepared for Clock Foundation export, with the sex
#'   column preserved as supplied.}
#' \item{savedFiles}{Object returned by [writePreprocessingPhenoOutputs()] when
#'   `saveOutputs = TRUE`, otherwise `NULL`.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#' See [dnaEPICO_preprocessingPheno-class] for a class-level overview.
#'
#' @examples
#' tmp <- tempdir()
#' pheno <- data.frame(
#'     Sample_Name = c("S1", "S2", "S3"),
#'     Timepoint = c("1", "1", "2"),
#'     Sex = c(0, 1, 0),
#'     stringsAsFactors = FALSE
#' )
#' beta <- matrix(
#'     c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60),
#'     nrow = 2,
#'     dimnames = list(c("cg1", "cg2"), pheno$Sample_Name)
#' )
#' m <- beta * 10
#' cn <- beta * 100
#' pheno_file <- file.path(tmp, "pheno.csv")
#' beta_path <- file.path(tmp, "beta.RData")
#' m_path <- file.path(tmp, "m.RData")
#' cn_path <- file.path(tmp, "cn.RData")
#' utils::write.csv(pheno, pheno_file, row.names = FALSE)
#' save(beta, file = beta_path)
#' save(m, file = m_path)
#' save(cn, file = cn_path)
#' result <- preprocessingPheno(
#'     phenoFile = pheno_file,
#'     betaPath = beta_path,
#'     mPath = m_path,
#'     cnPath = cn_path,
#'     SampleID = "Sample_Name",
#'     timeVar = "Timepoint",
#'     timepoints = "1,2",
#'     combineTimepoints = "1,2",
#'     outputPheno = file.path(tmp, "data", "preprocessingPheno"),
#'     outputRData = file.path(tmp, "rData", "preprocessingPheno", "metrics"),
#'     outputRDataMerge = file.path(
#'         tmp, "rData", "preprocessingPheno", "mergeData"
#'     ),
#'     sexColumn = "Sex",
#'     outputLogs = file.path(tmp, "logs"),
#'     outputDir = file.path(tmp, "clockFoundation"),
#'     saveOutputs = FALSE
#' )
#' stopifnot(inherits(result, "dnaEPICO_preprocessingPheno"))
#'
#' @seealso [dnaEPICO_preprocessingPheno-class]
#'
#' @export
preprocessingPheno <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
    sepType = NULL,
    betaPath = paste0(
        "rData/preprocessingMinfiEwasWater/metrics/",
        "beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData"
    ),
    mPath = paste0(
        "rData/preprocessingMinfiEwasWater/metrics/",
        "m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData"
    ),
    cnPath = paste0(
        "rData/preprocessingMinfiEwasWater/metrics/",
        "cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData"
    ),
    SampleID = "Sample_Name", timeVar = "Timepoint", timepoints = "1,2",
    combineTimepoints = "1,2", methylationScale = "beta",
    outputPheno = "data/preprocessingPheno",
    outputRData = "rData/preprocessingPheno/metrics",
    outputRDataMerge = "rData/preprocessingPheno/mergeData",
    sexColumn = "Sex", outputLogs = "logs",
    outputDir = "data/preprocessingPheno",
    verbose = FALSE, logs = FALSE, saveOutputs = FALSE
) {
    methylationScale <- normalizeMethylationScaleDnaEpico(methylationScale)
    methylationLabel <- methylationScaleResponseLabelDnaEpico(methylationScale)
    methylationObjectPrefix <-
        methylationScaleObjectPrefixDnaEpico(methylationScale)
    log_file <- "log_preprocessingPheno.txt"
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = outputLogs, log_file = log_file
    )
    config <- as.list(environment(), all.names = TRUE)
    logPreprocessingPhenoStartDnaEpico(config)
    withLoggedErrorsMinfiEwasWater(
        expr = runPreprocessingPhenoDnaEpico(config),
        log_path = log_path, verbose = verbose, context = "preprocessingPheno"
    )
}
