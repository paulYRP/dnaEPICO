#' Validate a loaded RGSet for svaEnmix
#'
#' @param RGSet Object loaded from `rgsetData`.
#' @param rgsetData Character. Source file path used in validation messages.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return The validated `RGSet`.
#'
#' @description
#' Internal helper that logs the class and dimensions of the loaded object and
#' stops with a clearer error when the saved file does not contain a usable
#' `RGChannelSet`-like object.
#'
#' @keywords internal
#' @noRd
svaEnmixGetRGSetSampleNames <- function(RGSet) {
    sample_names <- tryCatch(colnames(RGSet), error = function(...) NULL)

    if (length(sample_names) > 0L) {
        return(sample_names)
    }

    if (!isS4(RGSet)) {
        return(character(0))
    }

    slots <- methods::slotNames(RGSet)

    if ("colData" %in% slots) {
        col_data <- methods::slot(RGSet, "colData")
        sample_names <- rownames(as.data.frame(col_data))

        if (length(sample_names) > 0L) {
            return(sample_names)
        }
    }

    if ("phenoData" %in% slots) {
        pheno_data <- Biobase::pData(Biobase::phenoData(RGSet))
        sample_names <- rownames(pheno_data)

        if (length(sample_names) > 0L) {
            return(sample_names)
        }
    }

    character(0)
}

#' @keywords internal
#' @noRd
svaEnmixGetRGSetSampleCount <- function(RGSet) {
    rgset_dim <- tryCatch(dim(RGSet), error = function(...) NULL)

    if (length(rgset_dim) == 2L) {
        return(rgset_dim[[2L]])
    }

    sample_names <- svaEnmixGetRGSetSampleNames(RGSet)

    if (length(sample_names) > 0L) {
        return(length(sample_names))
    }

    if (!isS4(RGSet)) {
        return(NULL)
    }

    slots <- methods::slotNames(RGSet)

    if ("colData" %in% slots) {
        return(nrow(as.data.frame(methods::slot(RGSet, "colData"))))
    }

    if ("phenoData" %in% slots) {
        return(nrow(Biobase::pData(Biobase::phenoData(RGSet))))
    }

    NULL
}

#' @keywords internal
#' @noRd
svaEnmixSetRGSetSampleNames <- function(RGSet, sampleNames) {
    updated_rgset <- tryCatch(
        {
            colnames(RGSet) <- sampleNames
            RGSet
        },
        error = function(...) NULL
    )

    if (!is.null(updated_rgset)) {
        return(updated_rgset)
    }

    if (!isS4(RGSet)) {
        stop("Could not set sample names on the loaded RGSet object.",
            call. = FALSE
        )
    }

    slots <- methods::slotNames(RGSet)

    if ("colData" %in% slots) {
        col_data <- methods::slot(RGSet, "colData")
        rownames(col_data) <- sampleNames
        methods::slot(RGSet, "colData") <- col_data
        return(RGSet)
    }

    if ("phenoData" %in% slots) {
        pheno_data <- Biobase::phenoData(RGSet)
        rownames(Biobase::pData(pheno_data)) <- sampleNames
        Biobase::phenoData(RGSet) <- pheno_data
        return(RGSet)
    }

    stop("Could not set sample names on the loaded RGSet object.",
        call. = FALSE
    )
}

#' @keywords internal
#' @noRd
svaEnmixGetRGSetColData <- function(RGSet) {
    col_data <- tryCatch(SummarizedExperiment::colData(RGSet),
        error = function(...) NULL
    )

    if (!is.null(col_data)) {
        return(col_data)
    }

    if (!isS4(RGSet)) {
        stop("Could not retrieve column metadata from the loaded RGSet.",
            call. = FALSE
        )
    }

    slots <- methods::slotNames(RGSet)

    if ("colData" %in% slots) {
        return(methods::slot(RGSet, "colData"))
    }

    if ("phenoData" %in% slots) {
        return(Biobase::pData(Biobase::phenoData(RGSet)))
    }

    stop("Could not retrieve column metadata from the loaded RGSet.",
        call. = FALSE
    )
}

#' @keywords internal
#' @noRd
svaEnmixValidateRGSet <- function(
    RGSet, rgsetData, verbose = FALSE,
    logs = FALSE, log_dir = NULL, log_file = "log_svaEnmix.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    rgset_class <- paste(class(RGSet), collapse = ", ")
    rgset_dims <- dim(RGSet)
    rgset_ncol <- svaEnmixGetRGSetSampleCount(RGSet)
    slots <- if (isS4(RGSet)) {
        methods::slotNames(RGSet)
    } else {
        character(0)
    }
    log_lines <- c(paste("Loaded RGSet object class: ", rgset_class))

    if (length(rgset_dims) == 2L) {
        log_lines <- c(log_lines, paste(
            "Loaded RGSet dimensions:   ",
            paste(rgset_dims, collapse = " x ")
        ))
    }

    log_lines <- c(log_lines, paste(
        "Loaded RGSet sample count:  ",
        if (length(rgset_ncol) == 1L && !is.na(rgset_ncol)) {
            rgset_ncol
        } else {
            "unavailable"
        }
    ), paste("Loaded RGSet slots:         ", if (length(slots) >
        0L) {
        paste(slots, collapse = ", ")
    } else {
        "unavailable"
    }))

    emitLogMinfiEwasWater(log_lines, verbose = verbose, log_path = log_path)

    if (!(length(rgset_ncol) == 1L && !is.na(rgset_ncol))) {
        stop("The object loaded from ", rgsetData, " has class ",
            rgset_class,
                " and does not expose a usable sample count after loading.",
            call. = FALSE
        )
    }

    list(RGSet = RGSet, sampleCount = rgset_ncol)
}

#' Estimate surrogate variables from ENmix control probes
#'
#' Read the phenotype table and a saved `RGChannelSet`, estimate surrogate
#' variables from ENmix control probes, analyze their association with Sentrix
#' chip and position factors, and return a structured in-memory result. CSV,
#' `.RData`, text-summary, and figure outputs are written only when
#' `saveOutputs = TRUE`.
#'
#' @param phenoFile Character. Path to the phenotype file with cell-composition
#'   data. When `saveOutputs = TRUE`, the validated PC columns are appended to
#'   this file using a rollback-safe replacement.
#' @param rgsetData Character. Path to a saved `RGChannelSet` object. Both
#'   `.RData` and `.rds` files are supported.
#' @param sepType Character or `NULL`. Field separator used in `phenoFile`. Use
#'   `NULL` for a comma-separated file, `'\\t'` for a tab-delimited file, or
#'   another separator accepted by `utils::read.csv()`.
#' @param outputLogs Character. Directory used for log files when `logs = TRUE`.
#' @param nSamples Integer or `NA`. Number of rows to keep from the phenotype
#'   table. Use `NA` to keep all samples.
#' @param SampleID Character. Name of the phenotype column containing sample
#'   identifiers.
#' @param arrayType Character. Illumina array identifier assigned to
#'   `Biobase::annotation(RGSet)`.
#' @param annotationVersion Character. Annotation build assigned to
#'   `Biobase::annotation(RGSet)`.
#' @param SentrixIDColumn Character. Name of the chip identifier column in the
#'   phenotype data.
#' @param SentrixPositionColumn Character. Name of the chip position column in
#'   the phenotype data.
#' @param ctrlSvaPercVar Numeric. Proportion of control-probe variance explained
#'   when running `ENmix::ctrlsva()`.
#' @param ctrlSvaFlag Integer. Control-probe flag passed to `ENmix::ctrlsva()`.
#' @param scriptLabel Character. Label used to name output folders when
#'   `saveOutputs = TRUE`.
#' @param tiffWidth Integer. Width of saved TIFF plots in pixels.
#' @param tiffHeight Integer. Height of saved TIFF plots in pixels.
#' @param tiffRes Integer. Resolution in DPI for saved TIFF plots.
#' @param figureBaseDir Character. Base directory used for saved figure outputs
#'   when `saveOutputs = TRUE`.
#' @param dataBaseDir Character. Base directory used for saved CSV and text
#'   outputs when `saveOutputs = TRUE`.
#' @param rBaseDir Character. Base directory used for saved `.RData` outputs
#'   when `saveOutputs = TRUE`.
#' @param display Logical. If `TRUE`, draw plots on the active graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#'   The default is `FALSE`.
#' @param logs Logical. If `TRUE`, write the same progress messages to
#'   `outputLogs`. The default is `FALSE`.
#' @param saveOutputs Logical. If `TRUE`, write the CSV, `.RData`, text,
#'   and TIFF outputs to disk. The default is `FALSE`.
#'
#' @return A list with class `'dnaEPICO_svaEnmix'`.
#' \describe{
#'   \item{targets}{Phenotype table read from `phenoFile` after any optional row
#'   subsetting.}
#'   \item{RGSet}{Loaded `RGChannelSet` with sample names realigned to
#'   `targets[[SampleID]]`.}
#'   \item{svaData}{Object returned by [estimateSvaEnmixControls()] containing
#'   the surrogate-variable matrix and the control-probe settings used to
#'   estimate it.}
#'   \item{mergedPheno}{Phenotype table returned by [mergeSvaTargetsEnmix()]
#'   after the surrogate variables were appended as additional columns.}
#'   \item{analysisData}{Object returned by [analyzeSvaEnmix()] containing the
#'   surrogate-variable association models, ANOVA tables, and Sentrix metadata.}
#' \item{plotFiles}{Named list describing the plot file paths requested for the
#'   SVA figures. When `saveOutputs = FALSE`, the entries are typically `NULL`.}
#'   \item{savedFiles}{Object returned by [writeSvaEnmixOutputs()] when
#'   `saveOutputs = TRUE`, otherwise `NULL`.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#' See [dnaEPICO_svaEnmix-class] for a class-level overview.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' if (requireNamespace("minfiData", quietly = TRUE)) {
#'     ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#'     pheno_file <- file.path(tmp, "pheno.csv")
#'     rgset_path <- file.path(tmp, "RGSet.RData")
#'     RGSet <- ex$RGSet
#'     utils::write.csv(ex$targets, pheno_file, row.names = FALSE)
#'     save(RGSet, file = rgset_path)
#'     sva_result <- svaEnmix(
#'         phenoFile = pheno_file,
#'         rgsetData = rgset_path,
#'         SampleID = "Sample_Name",
#'         arrayType = "IlluminaHumanMethylation450k",
#'         annotationVersion = "ilmn12.hg19",
#'         SentrixIDColumn = "Sentrix_ID",
#'         SentrixPositionColumn = "Sentrix_Position",
#'         outputLogs = file.path(tmp, "logs"),
#'         figureBaseDir = file.path(tmp, "figures"),
#'         dataBaseDir = file.path(tmp, "data"),
#'         rBaseDir = file.path(tmp, "rData"),
#'         saveOutputs = FALSE
#'     )
#'     stopifnot(inherits(sva_result, "dnaEPICO_svaEnmix"))
#' }
#'
#' @seealso [dnaEPICO_svaEnmix-class]
#'
#' @export
svaEnmix <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
    rgsetData = "rData/preprocessingMinfiEwasWater/objects/RGSet.RData",
    sepType = NULL, outputLogs = "logs", nSamples = NA,
        SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
        annotationVersion = "20a1.hg38",
    SentrixIDColumn = "Sentrix_ID", SentrixPositionColumn = "Sentrix_Position",
    ctrlSvaPercVar = 0.9, ctrlSvaFlag = 1, scriptLabel = "svaEnmix",
    tiffWidth = 2000, tiffHeight = 1000, tiffRes = 150,
        figureBaseDir = "figures",
    dataBaseDir = "data", rBaseDir = "rData", display = FALSE,
    verbose = FALSE, logs = FALSE, saveOutputs = FALSE
) {
    log_file <- "log_svaEnmix.txt"
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = outputLogs,
        log_file = log_file
    )
    figure_dir <- file.path(figureBaseDir, scriptLabel)

    emitLogMinfiEwasWater(
        c(
            "==== Starting SVA Estimation with Enmix ====",
            paste("Start time:               ", format(Sys.time())),
            paste("Log file path:            ",
                if (is.null(log_path)) "disabled" else log_path),
            paste("Pheno file:               ", phenoFile), paste(
                "RGSet path:               ",
                rgsetData
            ), paste("Separator type:           ",
                if (is.null(resolveSeparatorMinfiEwasWater(sepType))) {
                "default (',')"
            } else {
                sepType
            }), paste("Sample limit:             ",
                if (is.na(nSamples)) "all" else nSamples),
            paste("SampleID column:          ", SampleID), paste(
                "Array type:               ",
                arrayType
            ), paste("Annotation version:       ", annotationVersion),
            paste("Sentrix ID column:        ", SentrixIDColumn),
            paste("Sentrix position column:  ", SentrixPositionColumn),
            paste("ctrlSva percvar:          ", ctrlSvaPercVar),
            paste("ctrlSva flag:             ", ctrlSvaFlag), paste(
                "Script label:             ",
                scriptLabel
            ), paste(
                "TIFF dimensions (WxH):    ",
                tiffWidth, " x ", tiffHeight, " @ ", tiffRes
            ), paste(
                "Display plots:            ",
                display
            ), paste("Save outputs:             ", saveOutputs),
            "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    withLoggedErrorsMinfiEwasWater(expr = {
        targets <- readPhenotypeTargets(
            phenoFile = phenoFile,
            sepType = sepType, nSamples = nSamples, SampleID = SampleID,
            verbose = verbose, logs = logs, log_dir = outputLogs,
            log_file = log_file
        )

        RGSet <- loadSavedObjectPreprocessingPheno(
            path = rgsetData,
            preferred_name = "RGSet"
        )
        rgset_validation <- svaEnmixValidateRGSet(
            RGSet = RGSet,
            rgsetData = rgsetData, verbose = verbose, logs = logs,
            log_dir = outputLogs, log_file = log_file
        )
        RGSet <- rgset_validation$RGSet
        rgset_ncol <- rgset_validation$sampleCount

        if (rgset_ncol != nrow(targets)) {
            stop("The saved RGSet contains ", rgset_ncol,
                " samples but the phenotype table contains ",
                nrow(targets), ".",
                call. = FALSE
            )
        }

        target_sample_ids <- validateSampleIdentifiersDnaEpico(
            targets[[SampleID]],
            paste0("Phenotype column '", SampleID, "'")
        )
        rgset_sample_ids <- svaEnmixGetRGSetSampleNames(RGSet)
        rgset_match <- matchSampleIdentifiersDnaEpico(
            query = target_sample_ids,
            reference = rgset_sample_ids, queryLabel = paste0(
                "phenotype column '",
                SampleID, "'"
            ), referenceLabel = "saved RGSet sample identifiers",
            requireSameSet = TRUE
        )
        RGSet <- RGSet[, rgset_match]

        RGSet <- svaEnmixSetRGSetSampleNames(RGSet = RGSet,
            sampleNames = target_sample_ids)
        Biobase::annotation(RGSet) <- c(array = arrayType,
            annotation = annotationVersion)

        emitLogMinfiEwasWater(
            c(
                paste(
                    "RGSet loaded with          ",
                    rgset_ncol, " samples."
                ), paste(
                    "Applied annotation:        ",
                    paste(Biobase::annotation(RGSet), collapse = ", ")
                ),
                "============================================================"
            ),
            verbose = verbose, log_path = log_path
        )

        svaData <- estimateSvaEnmixControls(
            RGSet = RGSet, ctrlSvaPercVar = ctrlSvaPercVar,
            ctrlSvaFlag = ctrlSvaFlag, verbose = verbose, logs = logs,
            log_dir = outputLogs, log_file = log_file
        )
        mergedPheno <- mergeSvaTargetsEnmix(
            targets = targets,
            sva = svaData$sva, SampleID = SampleID, verbose = verbose,
            logs = logs, log_dir = outputLogs, log_file = log_file
        )
        analysisData <- analyzeSvaEnmix(
            sva = svaData$sva, RGSet = RGSet,
            SentrixIDColumn = SentrixIDColumn,
                SentrixPositionColumn = SentrixPositionColumn,
            verbose = verbose, logs = logs, log_dir = outputLogs,
            log_file = log_file
        )

        plot_files <- list(
            sentrixID = plotSvaEnmix(
                analysisData = analysisData,
                plot = "sentrix_id", display = display,
                    file = if (isTRUE(saveOutputs)) {
                    file.path(
                        figure_dir,
                        "sva_SentrixID.tiff"
                    )
                } else {
                    NULL
                }, width = tiffWidth,
                height = tiffHeight, res = tiffRes, verbose = verbose,
                logs = logs, log_dir = outputLogs, log_file = log_file
            ),
            sentrixPosition = plotSvaEnmix(
                analysisData = analysisData,
                plot = "sentrix_position", display = display,
                file = if (isTRUE(saveOutputs)) {
                    file.path(
                        figure_dir,
                        "sva_SentrixPosition.tiff"
                    )
                } else {
                    NULL
                }, width = tiffWidth,
                height = tiffHeight, res = tiffRes, verbose = verbose,
                logs = logs, log_dir = outputLogs, log_file = log_file
            ),
            matrix = plotSvaEnmix(
                analysisData = analysisData,
                plot = "matrix", display = display,
                    file = if (isTRUE(saveOutputs)) {
                    file.path(
                        figure_dir,
                        "sva_SentrixIDPosition.tiff"
                    )
                } else {
                    NULL
                }, width = tiffWidth,
                height = tiffHeight, res = tiffRes, verbose = verbose,
                logs = logs, log_dir = outputLogs, log_file = log_file
            )
        )

        savedFiles <- NULL
        if (isTRUE(saveOutputs)) {
            savedFiles <- writeSvaEnmixOutputs(
                svaData = svaData,
                mergedPheno = mergedPheno, analysisData = analysisData,
                phenoFile = phenoFile, SampleID = SampleID, sepType = sepType,
                dataBaseDir = dataBaseDir, rBaseDir = rBaseDir,
                scriptLabel = scriptLabel, verbose = verbose,
                logs = logs, log_dir = outputLogs, log_file = log_file
            )
        }

        emitLogMinfiEwasWater(
            c(
                "==== Finished SVA Estimation with Enmix ====",
                paste("End time:                 ", format(Sys.time())),
                "============================================================"
            ),
            verbose = verbose, log_path = log_path
        )

        structure(list(
            targets = targets, RGSet = RGSet, svaData = svaData,
            mergedPheno = mergedPheno, analysisData = analysisData,
            plotFiles = plot_files, savedFiles = savedFiles,
            logFile = log_path
        ), class = "dnaEPICO_svaEnmix")
    }, log_path = log_path, verbose = verbose, context = "svaEnmix")
}
