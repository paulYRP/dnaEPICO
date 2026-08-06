#' Prepare phenotype and methylation matrices for downstream modeling
#'
#' Align the phenotype table with preprocessed beta, M-value, and copy-number
#' matrices, split the data by timepoint, prepare longitudinal objects for the
#' selected modeling scale, and build Clock Foundation export tables. The
#' function writes files only when `saveOutputs = TRUE`.
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
#'   object.
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
#'   building Clock Foundation exports.
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
#'   phenotype-plus-methylation object and the timepoint combination metadata.}
#'   \item{clockFoundation}{Object returned by
#'   [buildClockFoundationInputsPreprocessingPheno()] containing the beta table
#'   and phenotype table prepared for Clock Foundation export.}
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
#'     outputRDataMerge = file.path(tmp, "rData", "preprocessingPheno", "mergeData"),
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
        betaPath = "rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    mPath = "rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    cnPath = "rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
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
        logs = logs, log_dir = outputLogs,
        log_file = log_file
    )

    emitLogMinfiEwasWater(
        c(
            "==== Starting preprocessingPheno ====",
            paste("Start Time:               ", format(Sys.time())),
            paste("Log file path:            ",
                if (is.null(log_path)) "disabled" else log_path),
            paste("Phenotype file:           ", phenoFile), paste(
                "Separator type:           ",
                if (is.null(resolveSeparatorMinfiEwasWater(sepType))) {
                    "default (',')"
                } else {
                    sepType
                }
            ), paste("Beta path:                ", betaPath),
            paste("M-values path:            ", mPath), paste(
                "CN path:                  ",
                cnPath
            ), paste("Identifier column:        ", SampleID),
            paste("Timepoint column:         ", timeVar), paste(
                "Timepoints:               ",
                timepoints
            ), paste(
                "Combine timepoints:       ",
                combineTimepoints
            ), paste(
                "Merged modeling object:   ",
                methylationObjectPrefix, "*"
            ), "Clock Foundation scale:     Beta values",
            paste("Sex column:               ", sexColumn), paste(
                "Output phenotype dir:     ",
                outputPheno
            ), paste(
                "RData metrics dir:        ",
                outputRData
            ), paste(
                "RData merge dir:          ",
                outputRDataMerge
            ), paste(
                "Clock Foundation dir:     ",
                outputDir
            ), paste("Save outputs:             ", saveOutputs),
            "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    withLoggedErrorsMinfiEwasWater(expr = {
        pheno <- readPhenotypeTargets(
            phenoFile = phenoFile,
            sepType = sepType, SampleID = SampleID, verbose = verbose,
            logs = logs, log_dir = outputLogs, log_file = log_file
        )
        metricsData <- loadMetricsPreprocessingPheno(
            betaPath = betaPath,
            mPath = mPath, cnPath = cnPath, verbose = verbose,
            logs = logs, log_dir = outputLogs, log_file = log_file
        )
        timepointData <- splitTimepointsPreprocessingPheno(
            pheno = pheno,
            metricsData = metricsData, SampleID = SampleID, timeVar = timeVar,
            timepoints = timepoints, methylationScale = methylationScale,
            verbose = verbose, logs = logs, log_dir = outputLogs,
            log_file = log_file
        )
        combinedData <- combineTimepointsPreprocessingPheno(
            timepointData = timepointData,
            combineTimepoints = combineTimepoints,
                methylationScale = methylationScale,
            verbose = verbose, logs = logs, log_dir = outputLogs,
            log_file = log_file
        )
        clockFoundation <- buildClockFoundationInputsPreprocessingPheno(
            beta = metricsData$beta,
            pheno = pheno, SampleID = SampleID, sexColumn = sexColumn,
            verbose = verbose, logs = logs, log_dir = outputLogs,
            log_file = log_file
        )

        result <- list(
            pheno = pheno, metricsData = metricsData,
            timepointData = timepointData, combinedData = combinedData,
            methylationScale = methylationScale,
                methylationLabel = methylationLabel,
            methylationObjectPrefix = methylationObjectPrefix,
            clockFoundation = clockFoundation, savedFiles = NULL,
            logFile = log_path
        )

        if (isTRUE(saveOutputs)) {
            result$savedFiles <- writePreprocessingPhenoOutputs(
                preprocessingData = result,
                outputPheno = outputPheno, outputRData = outputRData,
                outputRDataMerge = outputRDataMerge, outputDir = outputDir,
                verbose = verbose, logs = logs, log_dir = outputLogs,
                log_file = log_file
            )
        }

        emitLogMinfiEwasWater(
            c(
                "==== Finished preprocessingPheno ====",
                paste("End Time:                 ", format(Sys.time())),
                "============================================================"
            ),
            verbose = verbose, log_path = log_path
        )

        structure(result, class = "dnaEPICO_preprocessingPheno")
    }, log_path = log_path, verbose = verbose, context = "preprocessingPheno")
}
