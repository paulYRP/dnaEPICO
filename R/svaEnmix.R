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
svaEnmixRGSetLogLines <- function(RGSet) {
    rgset_dims <- dim(RGSet)
    rgset_ncol <- svaEnmixGetRGSetSampleCount(RGSet)
    slots <- if (isS4(RGSet)) methods::slotNames(RGSet) else character(0)
    lines <- paste("Loaded RGSet object class: ",
        paste(class(RGSet), collapse = ", ")
    )
    if (length(rgset_dims) == 2L) {
        lines <- c(lines, paste(
            "Loaded RGSet dimensions:   ",
            paste(rgset_dims, collapse = " x ")
        ))
    }
    c(lines,
        paste("Loaded RGSet sample count:  ",
            if (length(rgset_ncol) == 1L && !is.na(rgset_ncol)) {
                rgset_ncol
            } else {
                "unavailable"
            }
        ),
        paste("Loaded RGSet slots:         ",
            if (length(slots) > 0L) paste(slots, collapse = ", ") else {
                "unavailable"
            }
        )
    )
}

svaEnmixValidateRGSet <- function(
    RGSet, rgsetData, verbose = FALSE,
    logs = FALSE, log_dir = NULL, log_file = "log_svaEnmix.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    rgset_ncol <- svaEnmixGetRGSetSampleCount(RGSet)
    emitLogMinfiEwasWater(
        svaEnmixRGSetLogLines(RGSet), verbose = verbose, log_path = log_path
    )

    if (!(length(rgset_ncol) == 1L && !is.na(rgset_ncol))) {
        rgset_class <- paste(class(RGSet), collapse = ", ")
        stop("The object loaded from ", rgsetData, " has class ",
            rgset_class,
                " and does not expose a usable sample count after loading.",
            call. = FALSE
        )
    }

    list(RGSet = RGSet, sampleCount = rgset_ncol)
}

logSvaEnmixStartDnaEpico <- function(config, logPath) {
    separator <- resolveSeparatorMinfiEwasWater(config$sepType)
    emitLogMinfiEwasWater(c(
        "==== Starting SVA Estimation with Enmix ====",
        paste("Start time:               ", format(Sys.time())),
        paste("Log file path:            ",
            if (is.null(logPath)) "disabled" else logPath
        ),
        paste("Pheno file:               ", config$phenoFile),
        paste("RGSet path:               ", config$rgsetData),
        paste("Separator type:           ",
            if (is.null(separator)) "default (',')" else config$sepType
        ),
        paste("Sample limit:             ",
            if (is.na(config$nSamples)) "all" else config$nSamples
        ),
        paste("SampleID column:          ", config$SampleID),
        paste("Array type:               ", config$arrayType),
        paste("Annotation version:       ", config$annotationVersion),
        paste("Sentrix ID column:        ", config$SentrixIDColumn),
        paste("Sentrix position column:  ", config$SentrixPositionColumn),
        paste("ctrlSva percvar:          ", config$ctrlSvaPercVar),
        paste("ctrlSva flag:             ", config$ctrlSvaFlag),
        paste("Script label:             ", config$scriptLabel),
        paste("TIFF dimensions (WxH):    ", config$tiffWidth, " x ",
            config$tiffHeight, " @ ", config$tiffRes
        ),
        paste("Display plots:            ", config$display),
        paste("Save outputs:             ", config$saveOutputs),
        "============================================================"
    ), verbose = config$verbose, log_path = logPath)
}

prepareSvaEnmixInputsDnaEpico <- function(config, logPath) {
    targets <- readPhenotypeTargets(
        phenoFile = config$phenoFile, sepType = config$sepType,
        nSamples = config$nSamples, SampleID = config$SampleID,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$logFile
    )
    RGSet <- loadSavedObjectPreprocessingPheno(
        path = config$rgsetData, preferred_name = "RGSet"
    )
    validation <- svaEnmixValidateRGSet(
        RGSet = RGSet, rgsetData = config$rgsetData,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$logFile
    )
    if (validation$sampleCount != nrow(targets)) {
        stop("The saved RGSet contains ", validation$sampleCount,
            " samples but the phenotype table contains ", nrow(targets), ".",
            call. = FALSE
        )
    }
    target_ids <- validateSampleIdentifiersDnaEpico(
        targets[[config$SampleID]],
        paste0("Phenotype column '", config$SampleID, "'")
    )
    rgset_match <- matchSampleIdentifiersDnaEpico(
        query = target_ids,
        reference = svaEnmixGetRGSetSampleNames(validation$RGSet),
        queryLabel = paste0("phenotype column '", config$SampleID, "'"),
        referenceLabel = "saved RGSet sample identifiers",
        requireSameSet = TRUE
    )
    RGSet <- validation$RGSet[, rgset_match]
    RGSet <- svaEnmixSetRGSetSampleNames(RGSet, target_ids)
    Biobase::annotation(RGSet) <- c(
        array = config$arrayType, annotation = config$annotationVersion
    )
    emitLogMinfiEwasWater(c(
        paste("RGSet loaded with          ", validation$sampleCount,
            " samples."
        ),
        paste("Applied annotation:        ",
            paste(Biobase::annotation(RGSet), collapse = ", ")
        ), "============================================================"
    ), verbose = config$verbose, log_path = logPath)
    list(targets = targets, RGSet = RGSet)
}

analyzeSvaEnmixWorkflowDnaEpico <- function(inputs, config) {
    sva_data <- estimateSvaEnmixControls(
        RGSet = inputs$RGSet, ctrlSvaPercVar = config$ctrlSvaPercVar,
        ctrlSvaFlag = config$ctrlSvaFlag, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$logFile
    )
    merged_pheno <- mergeSvaTargetsEnmix(
        targets = inputs$targets, sva = sva_data$sva,
        SampleID = config$SampleID, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$logFile
    )
    analysis_data <- analyzeSvaEnmix(
        sva = sva_data$sva, RGSet = inputs$RGSet,
        SentrixIDColumn = config$SentrixIDColumn,
        SentrixPositionColumn = config$SentrixPositionColumn,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$logFile
    )
    list(
        svaData = sva_data, mergedPheno = merged_pheno,
        analysisData = analysis_data
    )
}

plotSvaEnmixWorkflowFigureDnaEpico <- function(
    analysisData, plot, filename, config
) {
    file <- if (isTRUE(config$saveOutputs)) {
        file.path(config$figureDir, filename)
    } else {
        NULL
    }
    plotSvaEnmix(
        analysisData = analysisData, plot = plot,
        display = config$display, file = file,
        width = config$tiffWidth, height = config$tiffHeight,
        res = config$tiffRes, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$logFile
    )
}

plotSvaEnmixWorkflowDnaEpico <- function(analysisData, config) {
    plots <- c(
        sentrixID = "sentrix_id",
        sentrixPosition = "sentrix_position",
        matrix = "matrix", association = "association"
    )
    files <- c(
        sentrixID = "surrogateVariables_by_SentrixID.tiff",
        sentrixPosition = "surrogateVariables_by_SentrixPosition.tiff",
        matrix = "surrogateVariableMatrix_by_SentrixIDPosition.tiff",
        association = "surrogateVariable_technicalFactorAssociations.tiff"
    )
    lapply(names(plots), function(name) {
        plotSvaEnmixWorkflowFigureDnaEpico(
            analysisData, plots[[name]], files[[name]], config
        )
    }) |> stats::setNames(names(plots))
}

saveSvaEnmixWorkflowDnaEpico <- function(analysis, config) {
    if (!isTRUE(config$saveOutputs)) {
        return(NULL)
    }
    writeSvaEnmixOutputs(
        svaData = analysis$svaData, mergedPheno = analysis$mergedPheno,
        analysisData = analysis$analysisData, phenoFile = config$phenoFile,
        SampleID = config$SampleID, sepType = config$sepType,
        dataBaseDir = config$dataBaseDir, rBaseDir = config$rBaseDir,
        scriptLabel = config$scriptLabel, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$logFile
    )
}

runSvaEnmixWorkflowDnaEpico <- function(config, logPath) {
    inputs <- prepareSvaEnmixInputsDnaEpico(config, logPath)
    analysis <- analyzeSvaEnmixWorkflowDnaEpico(inputs, config)
    plot_files <- plotSvaEnmixWorkflowDnaEpico(
        analysis$analysisData, config
    )
    saved_files <- saveSvaEnmixWorkflowDnaEpico(analysis, config)
    emitLogMinfiEwasWater(c(
        "==== Finished SVA Estimation with Enmix ====",
        paste("End time:                 ", format(Sys.time())),
        "============================================================"
    ), verbose = config$verbose, log_path = logPath)
    structure(list(
        targets = inputs$targets, RGSet = inputs$RGSet,
        svaData = analysis$svaData, mergedPheno = analysis$mergedPheno,
        analysisData = analysis$analysisData, plotFiles = plot_files,
        savedFiles = saved_files, logFile = logPath
    ), class = "dnaEPICO_svaEnmix")
}

#' Estimate surrogate variables from ENmix control probes
#'
#' Estimate surrogate variables from ENmix control probes and analyze their
#' association with Sentrix chip and position factors. The function returns a
#' structured result and writes files only when `saveOutputs = TRUE`.
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
    config <- list(
        phenoFile = phenoFile, rgsetData = rgsetData, sepType = sepType,
        outputLogs = outputLogs, nSamples = nSamples, SampleID = SampleID,
        arrayType = arrayType, annotationVersion = annotationVersion,
        SentrixIDColumn = SentrixIDColumn,
        SentrixPositionColumn = SentrixPositionColumn,
        ctrlSvaPercVar = ctrlSvaPercVar, ctrlSvaFlag = ctrlSvaFlag,
        scriptLabel = scriptLabel, tiffWidth = tiffWidth,
        tiffHeight = tiffHeight, tiffRes = tiffRes,
        figureDir = file.path(figureBaseDir, scriptLabel),
        dataBaseDir = dataBaseDir, rBaseDir = rBaseDir, display = display,
        verbose = verbose, logs = logs, saveOutputs = saveOutputs,
        logFile = "log_svaEnmix.txt"
    )
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = outputLogs, log_file = config$logFile
    )
    logSvaEnmixStartDnaEpico(config, log_path)
    withLoggedErrorsMinfiEwasWater(
        expr = runSvaEnmixWorkflowDnaEpico(config, log_path),
        log_path = log_path, verbose = verbose, context = "svaEnmix"
    )
}
