normalizePreprocessingConfigDnaEpico <- function(config) {
    config$saveOutputs <- isTRUE(config$saveOutputs)
    config$removeSexMismatch <- validateLogicalScalarDnaEpico(
    config$removeSexMismatch, "removeSexMismatch"
    )
    if (!is.null(config$crossReactivePath)) {
    config$probeExclusionPath <- config$crossReactivePath
    }
    if (!is.null(config$crossReactiveIdColumn)) {
    config$probeExclusionIdColumn <- config$crossReactiveIdColumn
    }
    config$logFile <- "log_preprocessingMinfiEwasWater.txt"
    config$logPath <- resolveLogPathMinfiEwasWater(
    logs = config$logs, log_dir = config$outputLogs,
    log_file = config$logFile
    )
    config$normMethodList <- splitOptionMinfiEwasWater(
    config$normMethods,
    sep = ";"
    )
    config$manifestFlags <- normalizeEpicV2ManifestFlagsMinfiEwasWater(
    config$epicV2ManifestFlags
    )
    config
}

preprocessingProbeIdColumnTextDnaEpico <- function(value) {
    if (!length(value) || is.na(value[[1L]]) ||
    !nzchar(trimws(as.character(value[[1L]]))) ||
    identical(toupper(trimws(as.character(value[[1L]]))), "NULL")) {
    "auto"
    } else {
    value
    }
}

preprocessingInputLogLinesDnaEpico <- function(config) {
    c(
    paste("==== Starting", config$scriptLabel, "===="),
    paste("Start Time:               ", format(Sys.time())),
    paste(
        "Log file path:            ",
        if (is.null(config$logPath)) "disabled" else config$logPath
    ),
    paste("Phenotype file:           ", config$phenoFile),
    paste(
        "Separator type:           ",
        if (is.null(resolveSeparatorMinfiEwasWater(config$sepType))) {
        "default (',')"
        } else {
        config$sepType
        }
    ),
    paste("IDAT folder:              ", config$idatFolder),
    paste(
        "nSamples limit:           ",
        if (is.na(config$nSamples)) "all" else config$nSamples
    ),
    paste("SampleID column:          ", config$SampleID),
    paste("Array type:               ", config$arrayType),
    paste("Annotation version:       ", config$annotationVersion),
    paste("Force IDAT read:          ", config$idatForce),
    paste("Base RData folder:        ", config$baseDataFolder),
    paste("Base Figure folder:       ", config$figureBaseDir),
    paste(
        "TIFF size (w x h @ dpi):  ", config$tiffWidth, " x ",
        config$tiffHeight, " @ ", config$tiffRes
    ),
    paste("QC cutoff (median):       ", config$qcCutoff),
    paste("Detection P-value type:   ", config$detPtype),
    paste("Detection p-value threshold:", config$detPThreshold),
    paste(
        "Normalization methods:    ",
        paste(config$normMethodList, collapse = ", ")
    ),
    paste("Sex column:               ", config$sexColumn),
    paste("Remove sex mismatches:    ", config$removeSexMismatch),
    paste("Plot grouping variable:   ", config$plotGroupVar)
    )
}

preprocessingFilterLogLinesDnaEpico <- function(config) {
    c(
    "Probe filtering:",
    paste("  P-value threshold:      ", config$pvalThreshold),
    paste("  Chromosomes to remove:  ", config$chrToRemove),
    paste("  SNP positions filter:   ", config$snpsToRemove),
    paste("  MAF threshold:          ", config$mafThreshold),
    paste("  Probe-exclusion file(s):", config$probeExclusionPath),
    paste(
        "  Probe-exclusion ID col: ",
        preprocessingProbeIdColumnTextDnaEpico(
        config$probeExclusionIdColumn
        )
    ),
    paste("  Use EPICv2 manifest:    ", config$useEpicV2Manifest),
    paste("  EPICv2 manifest flags:  ", paste(paste(
        names(config$manifestFlags), config$manifestFlags,
        sep = "="
    ), collapse = ";")),
    "Cell composition (estimateLC):",
    paste("  Reference:              ", config$lcRef),
    paste("  Leading pheno order:    ", config$phenoOrder),
    paste("Display plots:            ", config$display),
    paste("Save outputs:             ", config$saveOutputs),
    "============================================================"
    )
}

logPreprocessingStartDnaEpico <- function(config) {
    emitLogMinfiEwasWater(c(
    preprocessingInputLogLinesDnaEpico(config),
    preprocessingFilterLogLinesDnaEpico(config)
    ), verbose = config$verbose, log_path = config$logPath)
}

preprocessingPathsDnaEpico <- function(config) {
    paths <- list(
    object = file.path(
        config$baseDataFolder, config$scriptLabel, "objects"
    ),
    norm = file.path(
        config$baseDataFolder, config$scriptLabel, "normObjects"
    ),
    metrics = file.path(
        config$baseDataFolder, config$scriptLabel, "metrics"
    ),
    filter = file.path(
        config$baseDataFolder, config$scriptLabel, "filterObjects"
    ),
    qc = file.path(config$baseDataFolder, config$scriptLabel, "qc"),
    metricsFigures = file.path(
        config$figureBaseDir, config$scriptLabel, "metrics"
    ),
    qcFigures = file.path(
        config$figureBaseDir, config$scriptLabel, "qc"
    ),
    enmix = file.path(
        config$figureBaseDir, config$scriptLabel, "enmix"
    )
    )
    if (config$saveOutputs) {
    for (path in c(unlist(paths, use.names = FALSE), config$lcPhenoDir)) {
        dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    }
    paths
}

preprocessingFigureFileDnaEpico <- function(config, directory, filename) {
    if (config$saveOutputs) file.path(directory, filename) else NULL
}

preprocessingPlotArgumentsDnaEpico <- function(config, file) {
    list(
    display = config$display, file = file,
    width = config$tiffWidth, height = config$tiffHeight,
    res = config$tiffRes, verbose = config$verbose,
    logs = config$logs, log_dir = config$outputLogs,
    log_file = config$logFile
    )
}

loadAndAssessPreprocessingSamplesDnaEpico <- function(config, paths) {
    targets <- readPhenotypeTargets(
    phenoFile = config$phenoFile, sepType = config$sepType,
    nSamples = config$nSamples, SampleID = config$SampleID,
    verbose = config$verbose, logs = config$logs,
    log_dir = config$outputLogs, log_file = config$logFile
    )
    rgset <- readRGSetMinfiEwasWater(
    idatFolder = config$idatFolder, targets = targets,
    SampleID = config$SampleID, arrayType = config$arrayType,
    annotationVersion = config$annotationVersion,
    force = config$idatForce, verbose = config$verbose,
    logs = config$logs, log_dir = config$outputLogs,
    log_file = config$logFile
    )
    plotCtrlMinfiEwasWater(
    RGSet = rgset,
    output_dir = if (config$saveOutputs) paths$enmix else NULL,
    verbose = config$verbose, logs = config$logs,
    log_dir = config$outputLogs, log_file = config$logFile
    )
    raw <- buildRawMinfiEwasWater(
    RGSet = rgset, verbose = config$verbose, logs = config$logs,
    log_dir = config$outputLogs, log_file = config$logFile
    )
    assessment <- assessSamplesMinfiEwasWater(
    rawData = raw, RGSet = rgset, qcCutoff = config$qcCutoff,
    detPtype = config$detPtype, detPThreshold = config$detPThreshold,
    verbose = config$verbose, logs = config$logs,
    log_dir = config$outputLogs, log_file = config$logFile
    )
    list(targets = targets, RGSet = rgset, assessment = assessment)
}

plotPreprocessingAssessmentDnaEpico <- function(state, config, paths) {
    qc_args <- preprocessingPlotArgumentsDnaEpico(
    config,
    preprocessingFigureFileDnaEpico(
        config, paths$qcFigures, "quality_control(MSet).tiff"
    )
    )
    do.call(plotAssessmentMinfiEwasWater, c(list(
    assessment = state$assessment, plot = "qc"
    ), qc_args))
    detection_args <- preprocessingPlotArgumentsDnaEpico(
    config,
    preprocessingFigureFileDnaEpico(
        config, paths$qcFigures,
        "detectionPvalue_sampleMean_byRank.tiff"
    )
    )
    do.call(plotAssessmentMinfiEwasWater, c(list(
    assessment = state$assessment, plot = "detection"
    ), detection_args))
    invisible(NULL)
}

filterInitialPreprocessingSamplesDnaEpico <- function(state, config, paths) {
    sample_data <- filterSamplesMinfiEwasWater(
    RGSet = state$RGSet, targets = state$targets,
    failedSamples = state$assessment$failedSamples,
    SampleID = config$SampleID, verbose = config$verbose,
    logs = config$logs, log_dir = config$outputLogs,
    log_file = config$logFile
    )
    raw <- buildRawMinfiEwasWater(
    RGSet = sample_data$RGSet, verbose = config$verbose,
    logs = config$logs, log_dir = config$outputLogs,
    log_file = config$logFile
    )
    args <- preprocessingPlotArgumentsDnaEpico(
    config,
    preprocessingFigureFileDnaEpico(
        config, paths$qcFigures, "densityBeta(MSet).tiff"
    )
    )
    do.call(plotRawDensityMinfiEwasWater, c(list(
    rawData = raw, targets = sample_data$targets,
    plotGroupVar = config$plotGroupVar
    ), args))
    list(sampleData = sample_data, rawData = raw)
}

applySexConcordanceDnaEpico <- function(
    sampleData, rawData, initialRGSet, assessment, config
) {
    sex <- predictSexMinfiEwasWater(
    rawData = rawData, targets = sampleData$targets,
    SampleID = config$SampleID, sexColumn = config$sexColumn,
    verbose = config$verbose, logs = config$logs,
    log_dir = config$outputLogs, log_file = config$logFile
    )
    mismatch_ids <- as.character(sex$mismatches[[config$SampleID]])
    mismatch_ids <- mismatch_ids[!is.na(mismatch_ids) & nzchar(mismatch_ids)]
    sex$removeSexMismatch <- config$removeSexMismatch
    sex$removedSampleIDs <- character(0)
    if (config$removeSexMismatch && length(mismatch_ids)) {
    if (length(mismatch_ids) >= ncol(sampleData$RGSet)) {
        stop(
        "Removing sex mismatches would leave no samples for ",
        "normalization.",
        call. = FALSE
        )
    }
    sampleData <- filterSamplesMinfiEwasWater(
        RGSet = sampleData$RGSet, targets = sex$targets,
        failedSamples = mismatch_ids, SampleID = config$SampleID,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$logFile
    )
    rawData <- buildRawMinfiEwasWater(
        RGSet = sampleData$RGSet, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$logFile
    )
    sex$removedSampleIDs <- mismatch_ids
    } else {
    sampleData$targets <- sex$targets
    }
    sex$retainedTargets <- sampleData$targets
    list(sampleData = sampleData, rawData = rawData, sexData = sex)
}

logSexConcordanceDnaEpico <- function(sexData, config) {
    emitLogMinfiEwasWater(c(
    paste("Remove sex mismatches:    ", config$removeSexMismatch),
    paste("Sex mismatches removed:   ", length(sexData$removedSampleIDs)),
    if (length(sexData$removedSampleIDs)) {
        paste(
        "Removed sample IDs:        ",
        paste(sexData$removedSampleIDs, collapse = ", ")
        )
    } else {
        "Removed sample IDs:         none"
    },
    "============================================================"
    ), verbose = config$verbose, log_path = config$logPath)
}

updatePreprocessingColDataDnaEpico <- function(sampleData, sexData, config) {
    col_data <- SummarizedExperiment::colData(sampleData$RGSet)
    if (config$sexColumn %in% colnames(col_data)) {
    col_data[[config$sexColumn]] <- sampleData$targets[[config$sexColumn]]
    }
    col_data$PredSex <- sampleData$targets$PredSex
    col_data$SexMismatch <- sampleData$targets$SexMismatch
    SummarizedExperiment::colData(sampleData$RGSet) <- col_data
    sampleData
}

plotPreprocessingSexDnaEpico <- function(
    sampleData, sexData, initialRGSet, assessment, config, paths
) {
    plotRetentionMinfiEwasWater(
    counts = c(
        `Read from IDAT` = ncol(initialRGSet),
        `Detection p-value filter` = ncol(initialRGSet) -
        length(unique(assessment$failedSamples)),
        `Sex-concordance filter` = ncol(sampleData$RGSet)
    ), unit = "samples", display = config$display,
    file = preprocessingFigureFileDnaEpico(
        config, paths$qcFigures, "sampleRetention_processingStages.tiff"
    ), width = config$tiffWidth, height = config$tiffHeight,
    res = config$tiffRes
    )
    for (type in c("predicted", "clinical")) {
    filename <- if (identical(type, "predicted")) {
        "sexPrediction(GSet).tiff"
    } else {
        "sexClinical(GSet).tiff"
    }
    args <- preprocessingPlotArgumentsDnaEpico(
        config,
        preprocessingFigureFileDnaEpico(config, paths$qcFigures, filename)
    )
    do.call(plotSexMinfiEwasWater, c(list(
        sexData = sexData, type = type
    ), args))
    }
    invisible(NULL)
}

normalizeFilterPreprocessingDnaEpico <- function(sampleData, config, paths) {
    norm <- normalizeMinfiEwasWater(
    sampleData = sampleData, sexColumn = config$sexColumn,
    normMethods = config$normMethods, verbose = config$verbose,
    logs = config$logs, log_dir = config$outputLogs,
    log_file = config$logFile
    )
    args <- preprocessingPlotArgumentsDnaEpico(
    config,
    preprocessingFigureFileDnaEpico(
        config, paths$qcFigures, "sexComparison_RawNorm(MSetF).tiff"
    )
    )
    do.call(plotNormalizationMinfiEwasWater, c(list(
    RGSet = sampleData$RGSet, normData = norm,
    targets = sampleData$targets, sexColumn = config$sexColumn
    ), args))
    filtered <- filterProbesMinfiEwasWater(
    normData = norm, RGSet = sampleData$RGSet,
    pvalThreshold = config$pvalThreshold,
    chrToRemove = config$chrToRemove,
    snpsToRemove = config$snpsToRemove,
    mafThreshold = config$mafThreshold,
    probeExclusionPath = config$probeExclusionPath,
    probeExclusionIdColumn = config$probeExclusionIdColumn,
    useEpicV2Manifest = config$useEpicV2Manifest,
    epicV2ManifestFlags = config$epicV2ManifestFlags,
    detPtype = config$detPtype, verbose = config$verbose,
    logs = config$logs, log_dir = config$outputLogs,
    log_file = config$logFile
    )
    list(normData = norm, filterData = filtered)
}

plotProbeRetentionDnaEpico <- function(filterData, config, paths) {
    plotRetentionMinfiEwasWater(
    counts = stats::setNames(as.numeric(filterData$counts), c(
        "Normalised", "Detection p-value filter", "Chromosome filter",
        "SNP filter", "Probe-exclusion filter"
    )), unit = "CpGs", display = config$display,
    file = preprocessingFigureFileDnaEpico(
        config, paths$qcFigures, "probeRetention_processingStages.tiff"
    ), width = config$tiffWidth, height = config$tiffHeight,
    res = config$tiffRes
    )
    invisible(NULL)
}

extractPreprocessingMetricsDnaEpico <- function(sampleData, filterData,
                                                config, paths) {
    metrics <- extractMetricsMinfiEwasWater(
    filteredData = filterData, verbose = config$verbose,
    logs = config$logs, log_dir = config$outputLogs,
    log_file = config$logFile
    )
    filenames <- c(
    mds = paste0(
        "examineMDS_PostFilteringCrossRect(",
        "MSetF_Flt_Rxy_Ds_Rc).tiff"
    ),
    density = "densityBeta&M(MSetF_Flt_Rxy_Ds_Rc).tiff"
    )
    for (plot_type in names(filenames)) {
    args <- preprocessingPlotArgumentsDnaEpico(
        config,
        preprocessingFigureFileDnaEpico(
        config, paths$metricsFigures, filenames[[plot_type]]
        )
    )
    do.call(plotMetricsMinfiEwasWater, c(list(
        metricsData = metrics, targets = sampleData$targets,
        plot = plot_type, plotGroupVar = config$plotGroupVar,
        sexColumn = config$sexColumn
    ), args))
    }
    metrics
}

estimatePreprocessingCellsDnaEpico <- function(
    sampleData, metricsData, config, paths
) {
    cells <- estimateLCMinfiEwasWater(
    beta = metricsData$beta, targets = sampleData$targets,
    SampleID = config$SampleID, lcRef = config$lcRef,
    phenoOrder = config$phenoOrder, constrained = FALSE,
    verbose = config$verbose, logs = config$logs,
    log_dir = config$outputLogs, log_file = config$logFile
    )
    plotCellCompositionMinfiEwasWater(
    lcData = cells, display = config$display,
    file = preprocessingFigureFileDnaEpico(
        config, paths$metricsFigures,
        "cellComposition_estimatedDistributions.tiff"
    ), width = config$tiffWidth, height = config$tiffHeight,
    res = config$tiffRes
    )
    cells
}

savePreprocessingCoreObjectsDnaEpico <- function(state, paths) {
    objects <- list(
    RGSet = state$sampleData$RGSet,
    MSet = state$rawData$MSet,
    RatioSet = state$rawData$RatioSet,
    GSet = state$rawData$GSet
    )
    for (name in names(objects)) {
    saveNamedObjectMinfiEwasWater(
        objects[[name]], name, file.path(paths$object, paste0(name, ".RData"))
    )
    }
    saveNamedObjectMinfiEwasWater(
    state$assessment$detP, "detP",
    file.path(paths$qc, "detP_RGSet.RData")
    )
    invisible(NULL)
}

savePreprocessingFilteredObjectsDnaEpico <- function(state, paths, config) {
    filtered <- list(
    MSetF_Flt = state$filterData$detPFiltered,
    MSetF_Flt_Rxy = state$filterData$chrFiltered,
    MSetF_Flt_Rxy_Ds = state$filterData$snpFiltered,
    MSetF_Flt_Rxy_Ds_Rc = state$filterData$filtered
    )
    files <- c(
    "removProbes_MSetF_Flt.RData", "removChrXY_MSetF_Flt_Rxy.RData",
    paste0(
        "removSNPs_MAF", config$mafThreshold,
        "_MSetF_Flt_Rxy_Ds.RData"
    ), "removCrossReactive_MSetF_Flt_Rxy_Ds_Rc.RData"
    )
    for (index in seq_along(filtered)) {
    saveNamedObjectMinfiEwasWater(
        filtered[[index]], names(filtered)[[index]],
        file.path(paths$filter, files[[index]])
    )
    }
    invisible(NULL)
}

savePreprocessingMetricsDnaEpico <- function(state, paths) {
    filenames <- c(
    m = "m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    beta = "beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    cn = "cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData"
    )
    for (name in names(filenames)) {
    saveNamedObjectMinfiEwasWater(
        state$metricsData[[name]], name,
        file.path(paths$metrics, filenames[[name]])
    )
    }
    invisible(NULL)
}

savePreprocessingOutputsDnaEpico <- function(state, paths, config) {
    if (!config$saveOutputs) {
    return(invisible(NULL))
    }
    savePreprocessingCoreObjectsDnaEpico(state, paths)
    for (index in seq_along(state$normData$methods)) {
    saveNamedObjectMinfiEwasWater(
        state$normData$normalized[[index]], "normObj",
        file.path(paths$norm, paste0(
        "norm_", state$normData$methods[[index]], "_RGSet.RData"
        ))
    )
    }
    savePreprocessingFilteredObjectsDnaEpico(state, paths, config)
    savePreprocessingMetricsDnaEpico(state, paths)
    writePhenoLCMinfiEwasWater(
    lcData = state$lcData,
    file = file.path(config$lcPhenoDir, "phenoLC.csv"),
    verbose = config$verbose, logs = config$logs,
    log_dir = config$outputLogs, log_file = config$logFile
    )
    invisible(NULL)
}

runPreprocessingWorkflowDnaEpico <- function(config) {
    paths <- preprocessingPathsDnaEpico(config)
    state <- loadAndAssessPreprocessingSamplesDnaEpico(config, paths)
    plotPreprocessingAssessmentDnaEpico(state, config, paths)
    initial <- filterInitialPreprocessingSamplesDnaEpico(state, config, paths)
    sex <- applySexConcordanceDnaEpico(
    initial$sampleData, initial$rawData,
    state$RGSet, state$assessment, config
    )
    logSexConcordanceDnaEpico(sex$sexData, config)
    sex$sampleData <- updatePreprocessingColDataDnaEpico(
    sex$sampleData, sex$sexData, config
    )
    plotPreprocessingSexDnaEpico(
    sex$sampleData, sex$sexData, state$RGSet,
    state$assessment, config, paths
    )
    normalized <- normalizeFilterPreprocessingDnaEpico(
    sex$sampleData, config, paths
    )
    plotProbeRetentionDnaEpico(normalized$filterData, config, paths)
    metrics <- extractPreprocessingMetricsDnaEpico(
    sex$sampleData, normalized$filterData, config, paths
    )
    cells <- estimatePreprocessingCellsDnaEpico(
    sex$sampleData, metrics, config, paths
    )
    result <- list(
    sampleData = sex$sampleData, rawData = sex$rawData,
    assessment = state$assessment, sexData = sex$sexData,
    normData = normalized$normData, filterData = normalized$filterData,
    metricsData = metrics, lcData = cells
    )
    savePreprocessingOutputsDnaEpico(result, paths, config)
    emitLogMinfiEwasWater(c(
    paste("==== Finished", config$scriptLabel, "===="),
    paste("End Time:                 ", format(Sys.time())),
    "============================================================"
    ), verbose = config$verbose, log_path = config$logPath)
    structure(list(
    targets = result$sampleData$targets,
    RGSet = result$sampleData$RGSet, rawData = result$rawData,
    assessment = result$assessment, sexData = result$sexData,
    normData = result$normData, filterData = result$filterData,
    metricsData = result$metricsData, lcData = result$lcData,
    logFile = config$logPath
    ), class = "dnaEPICO_preprocessingMinfiEwasWater")
}

#' Convenience preprocessing pipeline for Illumina methylation arrays
#'
#' Run the `dnaEPICO` preprocessing workflow with the package's
#' minfi/ENmix/wateRmelon helper functions. The function returns the output from
#' each stage and writes files only when `saveOutputs = TRUE`.
#'
#' @param phenoFile Character. Path to the phenotype CSV file.
#' @param idatFolder Character. Directory containing the IDAT files.
#' @param outputLogs Character. Directory used for log files when `logs = TRUE`.
#' @param nSamples Integer or `NA`. Number of rows to keep from the phenotype
#'   table. Use `NA` to keep all samples.
#' @param SampleID Character. Name of the phenotype column containing sample
#'   identifiers.
#' @param arrayType Character. Illumina array identifier passed to
#'   `Biobase::annotation()`, for example `'IlluminaHumanMethylationEPICv2'`.
#' @param annotationVersion Character. Annotation build passed to
#'   `Biobase::annotation()`, for example `'20a1.hg38'` or `'ilmn12.hg19'`.
#' @param idatForce Logical. Passed to `minfi::read.metharray.exp()` when
#'   reading IDAT files. Use `TRUE` only after confirming that the selected IDAT
#'   files should be read together.
#' @param scriptLabel Character. Label used to name output folders when
#'   `saveOutputs = TRUE`.
#' @param baseDataFolder Character. Base directory used for saved `.RData`
#'   outputs when `saveOutputs = TRUE`.
#' @param figureBaseDir Character. Base directory used for saved figure outputs
#'   when `saveOutputs = TRUE`.
#' @param sepType Character or `NULL`. Field separator used in `phenoFile`. Use
#'   `NULL` for a comma-separated file, `'\\t'` for a tab-delimited file, or
#'   another separator accepted by `utils::read.csv()`.
#' @param tiffWidth Integer. Width of saved TIFF plots in pixels.
#' @param tiffHeight Integer. Height of saved TIFF plots in pixels.
#' @param tiffRes Integer. Resolution in DPI for saved TIFF plots.
#' @param qcCutoff Numeric. QC cutoff passed to `minfi::plotQC()`.
#' @param detPtype Character. Detection P-value mode passed to
#'   `minfi::detectionP()`. Common values in minfi workflows are `'m+u'` and
#'   `'negative'`. The default here is `'m+u'`.
#' @param detPThreshold Numeric. Samples with mean detection P value above this
#'   threshold are removed.
#' @param normMethods Character vector or semicolon-separated string of
#'   normalization methods. Supported values are `'adjustedfunnorm'`,
#'   `'funnorm'`, `'illumina'`, `'quantile'`, and `'swan'`.
#' @param sexColumn Character. Name of the phenotype column containing reported
#'   sex. Sex-aware normalization methods use `PredSex` when this column is
#'   missing, blank, unknown, unsupported, or otherwise not coded as female or
#'   male. Every substitution is recorded without overwriting this column.
#' @param removeSexMismatch Logical. If `TRUE`, remove samples whose reported
#'   sex and methylation-predicted sex are both known and disagree. Missing or
#'   unknown sex values are retained. The default is `FALSE`.
#' @param pvalThreshold Numeric. Probe-level detection P-value threshold used in
#'   the probe filter.
#' @param chrToRemove Character vector or comma-separated string of chromosome
#'   names to remove, for example `'chrX,chrY'`.
#' @param snpsToRemove Character vector or comma-separated string of SNP probe
#'   types to remove, for example `'SBE,CpG'`.
#' @param mafThreshold Numeric. Minor allele frequency threshold passed to
#'   `minfi::dropLociWithSnps()`.
#' @param probeExclusionPath Character vector or semicolon-separated string of
#'   CSV files containing probe IDs to remove.
#' @param probeExclusionIdColumn Character or `NULL`. Column containing probe
#'   IDs. When `NULL` or `''`, each file is auto-detected using `ProbeID`,
#'   `TargetID`, `IlmnID`, or `Name`, then falling back to an unlabeled or
#'   probe-like first column.
#' @param useEpicV2Manifest Logical. If `TRUE`, also remove EPICv2 probes
#'   flagged in the Peters et al. expanded manifest from AnnotationHub resource
#'   `AH116484`.
#' @param epicV2ManifestFlags Named logical vector controlling which EPICv2
#'   manifest flags are removed. Defaults remove `CH_WGBS_evidence`, `CH_BLAT`,
#'   and `MissingPos`, but not `MismatchPos`.
#' @param crossReactivePath Deprecated alias for `probeExclusionPath`.
#' @param crossReactiveIdColumn Deprecated alias for `probeExclusionIdColumn`.
#' @param plotGroupVar Character. Phenotype column used for density and MDS
#'   grouping plots.
#' @param lcRef Character. Reference panel used for cell composition estimation.
#'   `'saliva'` and `'salivaEPIC'` use `estimateLC()`. Other values are passed
#'   to `ENmix::estimateCellProp()`.
#' @param phenoOrder Character vector or semicolon-separated phenotype columns
#'   to place first in the merged `phenoLC` table.
#' @param lcPhenoDir Character. Directory used for the saved `phenoLC.csv` file
#'   when `saveOutputs = TRUE`.
#' @param display Logical. If `TRUE`, draw plots on the active graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#'   The default is `FALSE`.
#' @param logs Logical. If `TRUE`, write log messages to `outputLogs`. The
#'   default is `FALSE`.
#' @param saveOutputs Logical. If `TRUE`, write the `.RData`, figure, and
#'   `phenoLC.csv` outputs to disk. The default is `FALSE`.
#'
#' @return A list with class `'dnaEPICO_preprocessingMinfiEwasWater'`.
#' \describe{
#'   \item{targets}{Filtered phenotype table aligned to the retained samples.}
#'   \item{RGSet}{Filtered `RGChannelSet` used in downstream preprocessing and
#'   available for direct interactive inspection.}
#'   \item{rawData}{Object returned by [buildRawMinfiEwasWater()] containing the
#'   raw `MSet`, `RatioSet`, and genome-mapped object derived from `RGSet`.}
#'   \item{assessment}{Object returned by [assessSamplesMinfiEwasWater()]
#'   containing detection P values, QC summaries, and failed-sample tracking.}
#'   \item{sexData}{Object returned by [predictSexMinfiEwasWater()] containing
#'   predicted sex labels, mismatch summaries, normalization-sex provenance,
#'   plotting data, and the sample IDs removed when
#'   `removeSexMismatch = TRUE`.}
#'   \item{normData}{Object returned by [normalizeMinfiEwasWater()] containing
#'   the selected normalized objects, method metadata, and the reported-sex or
#'   `PredSex` value used for each sample.}
#'   \item{filterData}{Object returned by [filterProbesMinfiEwasWater()]
#'   containing the probe-filtered methylation objects at each filtering stage.}
#'   \item{metricsData}{Object returned by [extractMetricsMinfiEwasWater()]
#'   containing the beta-value, M-value, and copy-number matrices used by later
#'   workflow steps.}
#'   \item{lcData}{Object returned by [estimateLCMinfiEwasWater()] containing
#'   the estimated cell-type proportions and the phenotype table augmented with
#'   those proportions.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#' See [dnaEPICO_preprocessingMinfiEwasWater-class] for a class-level overview.
#'
#' @examples
#' if (requireNamespace("minfiData", quietly = TRUE) &&
#'   requireNamespace(
#'     "IlluminaHumanMethylation450kmanifest",
#'     quietly = TRUE
#'   ) &&
#'   requireNamespace(
#'     "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#'     quietly = TRUE
#'   )) {
#'   ex <- dnaEPICO:::exampleMinfiIdatInputsDnaEpico(n = 4)
#'   result <- preprocessingMinfiEwasWater(
#'     phenoFile = ex$phenoFile,
#'     idatFolder = ex$idatFolder,
#'     outputLogs = file.path(ex$tempDir, "logs"),
#'     nSamples = 4,
#'     SampleID = "Sample_Name",
#'     arrayType = ex$arrayType,
#'     annotationVersion = ex$annotationVersion,
#'     scriptLabel = "preprocessingMinfiEwasWater",
#'     baseDataFolder = file.path(ex$tempDir, "rData"),
#'     figureBaseDir = file.path(ex$tempDir, "figures"),
#'     detPThreshold = 1,
#'     normMethods = "quantile",
#'     sexColumn = "Sex",
#'     pvalThreshold = 1,
#'     chrToRemove = "",
#'     snpsToRemove = "SBE",
#'     mafThreshold = 1,
#'     probeExclusionPath = ex$probeExclusionPath,
#'     plotGroupVar = "Sex",
#'     lcRef = "saliva",
#'     phenoOrder = "Sample_Name;Sex;Basename;Sentrix_ID;Sentrix_Position",
#'     lcPhenoDir = ex$tempDir,
#'     saveOutputs = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
#'   )
#'   inherits(result, "dnaEPICO_preprocessingMinfiEwasWater")
#' }
#'
#' @seealso [dnaEPICO_preprocessingMinfiEwasWater-class]
#'
#' @export
preprocessingMinfiEwasWater <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/pheno.csv",
    idatFolder = "data/preprocessingMinfiEwasWater/idats",
    outputLogs = "logs", nSamples = NA, SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
    annotationVersion = "20a1.hg38", idatForce = FALSE,
    scriptLabel = "preprocessingMinfiEwasWater",
    baseDataFolder = "rData", figureBaseDir = "figures", sepType = NULL,
    tiffWidth = 2000, tiffHeight = 1000, tiffRes = 150,
    qcCutoff = 10.5, detPtype = "m+u", detPThreshold = 0.05,
    normMethods = "adjustedfunnorm", sexColumn = "Sex",
    removeSexMismatch = FALSE, pvalThreshold = 0.01,
    chrToRemove = "chrX,chrY", snpsToRemove = "SBE,CpG",
    mafThreshold = 0.1,
    probeExclusionPath = paste0(
    "data/preprocessingMinfiEwasWater/",
    "12864_2024_10027_MOESM8_ESM.csv"
    ),
    probeExclusionIdColumn = NULL, useEpicV2Manifest = FALSE,
    epicV2ManifestFlags = c(
    CH_WGBS_evidence = TRUE, CH_BLAT = TRUE,
    MissingPos = TRUE, MismatchPos = FALSE
    ),
    plotGroupVar = "Sex", lcRef = "salivaEPIC",
    phenoOrder = paste0(
    "Sample_Name;Timepoint;Sex;PredSex;Basename;",
    "Sentrix_ID;Sentrix_Position"
    ),
    lcPhenoDir = "data/preprocessingMinfiEwasWater",
    display = FALSE, verbose = FALSE, logs = FALSE,
    saveOutputs = FALSE, crossReactivePath = NULL,
    crossReactiveIdColumn = NULL
) {
    config <- normalizePreprocessingConfigDnaEpico(
    as.list(environment(), all.names = TRUE)
    )
    logPreprocessingStartDnaEpico(config)
    withLoggedErrorsMinfiEwasWater(
    expr = runPreprocessingWorkflowDnaEpico(config),
    log_path = config$logPath, verbose = config$verbose,
    context = "preprocessingMinfiEwasWater"
    )
}
