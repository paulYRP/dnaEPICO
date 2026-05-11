#' Convenience preprocessing pipeline for Illumina methylation arrays
#'
#' Run the `dnaEPICO` preprocessing workflow as a convenience wrapper around the
#' smaller minfi/ENmix/wateRmelon helper functions in this package. The wrapper
#' now returns a structured result object containing the in-memory outputs from
#' each stage. Legacy files are written only when `saveOutputs = TRUE`.
#'
#' @param phenoFile Character. Path to the phenotype CSV file.
#' @param idatFolder Character. Directory containing the IDAT files.
#' @param outputLogs Character. Directory used for log files when `logs = TRUE`.
#' @param nSamples Integer or `NA`. Number of rows to keep from the phenotype
#'   table. Use `NA` to keep all samples.
#' @param SampleID Character. Name of the phenotype column containing sample
#'   identifiers.
#' @param arrayType Character. Illumina array identifier passed to
#'   `Biobase::annotation()`, for example `"IlluminaHumanMethylationEPICv2"`.
#' @param annotationVersion Character. Annotation build passed to
#'   `Biobase::annotation()`, for example `"20a1.hg38"` or `"ilmn12.hg19"`.
#' @param scriptLabel Character. Label used to name output folders when
#'   `saveOutputs = TRUE`.
#' @param baseDataFolder Character. Base directory used for saved `.RData`
#'   outputs when `saveOutputs = TRUE`.
#' @param figureBaseDir Character. Base directory used for saved figure outputs
#'   when `saveOutputs = TRUE`.
#' @param sepType Character. Field separator used in `phenoFile`. Use `""` for
#'   a comma-separated file, `"\\t"` for a tab-delimited file, or another
#'   separator accepted by `utils::read.csv()`.
#' @param tiffWidth Integer. Width of saved TIFF plots in pixels.
#' @param tiffHeight Integer. Height of saved TIFF plots in pixels.
#' @param tiffRes Integer. Resolution in DPI for saved TIFF plots.
#' @param qcCutoff Numeric. QC cutoff passed to `minfi::plotQC()`.
#' @param detPtype Character. Detection P-value mode passed to
#'   `minfi::detectionP()`. Common values in minfi workflows are `"m+u"` and
#'   `"negative"`. The default here is `"m+u"`.
#' @param detPThreshold Numeric. Samples with mean detection P value above this
#'   threshold are removed.
#' @param normMethods Character vector or semicolon-separated string of
#'   normalization methods. Supported values are `"adjustedfunnorm"`,
#'   `"funnorm"`, `"illumina"`, `"quantile"`, and `"swan"`.
#' @param sexColumn Character. Name of the phenotype column containing reported
#'   sex.
#' @param pvalThreshold Numeric. Probe-level detection P-value threshold used in
#'   the probe filter.
#' @param chrToRemove Character vector or comma-separated string of chromosome
#'   names to remove, for example `"chrX,chrY"`.
#' @param snpsToRemove Character vector or comma-separated string of SNP probe
#'   types to remove, for example `"SBE,CpG"`.
#' @param mafThreshold Numeric. Minor allele frequency threshold passed to
#'   `minfi::dropLociWithSnps()`.
#' @param crossReactivePath Character. Path to a CSV file containing a `ProbeID`
#'   column of cross-reactive probes to remove.
#' @param plotGroupVar Character. Phenotype column used for density and MDS
#'   grouping plots.
#' @param lcRef Character. Reference panel used for cell composition estimation.
#'   `"saliva"` and `"salivaEPIC"` use `estimateLC()`. Other values are passed
#'   to `ENmix::estimateCellProp()`.
#' @param phenoOrder Character vector or semicolon-separated string describing
#'   which phenotype columns should appear first in the merged `phenoLC` table.
#' @param lcPhenoDir Character. Directory used for the saved `phenoLC.csv` file
#'   when `saveOutputs = TRUE`.
#' @param display Logical. If `TRUE`, draw plots on the active graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#'   The default is `FALSE`.
#' @param logs Logical. If `TRUE`, write log messages to `outputLogs`. The
#'   default is `FALSE`.
#' @param saveOutputs Logical. If `TRUE`, write the legacy `.RData`, figure, and
#'   `phenoLC.csv` outputs to disk. The default is `FALSE`, so the function can
#'   be used in the more traditional in-memory Bioconductor style.
#'
#' @return A list with class `"dnaEPICO_preprocessingMinfiEwasWater"`.
#' \describe{
#'   \item{targets}{Filtered phenotype table aligned to the retained samples.}
#'   \item{RGSet}{Filtered `RGChannelSet` used in downstream preprocessing and
#'   available for direct interactive inspection.}
#'   \item{rawData}{Object returned by [buildRawMinfiEwasWater()] containing the
#'   raw `MSet`, `RatioSet`, and genome-mapped object derived from `RGSet`.}
#'   \item{assessment}{Object returned by [assessSamplesMinfiEwasWater()]
#'   containing detection P values, QC summaries, and failed-sample tracking.}
#'   \item{sexData}{Object returned by [predictSexMinfiEwasWater()] containing
#'   predicted sex labels, mismatch summaries, and plotting data.}
#'   \item{normData}{Object returned by [normalizeMinfiEwasWater()] containing
#'   the requested normalized objects and metadata on the methods that were run.}
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
#'     requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE) &&
#'     requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
#'   ex <- dnaEPICO:::exampleMinfiIdatInputsDnaEpico()
#'   result <- preprocessingMinfiEwasWater(
#'     phenoFile = ex$phenoFile,
#'     idatFolder = ex$idatFolder,
#'     outputLogs = file.path(ex$tempDir, "logs"),
#'     nSamples = 6,
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
#'     crossReactivePath = ex$crossReactivePath,
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
    outputLogs = "logs",
    nSamples = NA,
    SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
    annotationVersion = "20a1.hg38",
    scriptLabel = "preprocessingMinfiEwasWater",
    baseDataFolder = "rData",
    figureBaseDir = "figures",
    sepType = "",
    tiffWidth = 2000,
    tiffHeight = 1000,
    tiffRes = 150,
    qcCutoff = 10.5,
    detPtype = "m+u",
    detPThreshold = 0.05,
    normMethods = "adjustedfunnorm",
    sexColumn = "Sex",
    pvalThreshold = 0.01,
    chrToRemove = "chrX,chrY",
    snpsToRemove = "SBE,CpG",
    mafThreshold = 0.1,
    crossReactivePath =
      "data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv",
    plotGroupVar = "Sex",
    lcRef = "salivaEPIC",
    phenoOrder = "Sample_Name;Timepoint;Sex;PredSex;Basename;Sentrix_ID;Sentrix_Position",
    lcPhenoDir = "data/preprocessingMinfiEwasWater",
    display = FALSE,
    verbose = FALSE,
    logs = FALSE,
    saveOutputs = FALSE
) {
  log_file <- "log_preprocessingMinfiEwasWater.txt"
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  norm_method_list <- splitOptionMinfiEwasWater(normMethods, sep = ";")

  emitLogMinfiEwasWater(
    c(
      paste("==== Starting", scriptLabel, "===="),
      paste("Start Time:               ", format(Sys.time())),
      paste(
        "Log file path:            ",
        if (is.null(log_path)) "disabled" else log_path
      ),
      paste("Phenotype file:           ", phenoFile),
      paste(
        "Separator type:           ",
        if (is.null(resolveSeparatorMinfiEwasWater(sepType))) {
          "default (',')"
        } else {
          sepType
        }
      ),
      paste("IDAT folder:              ", idatFolder),
      paste("nSamples limit:           ", if (is.na(nSamples)) "all" else nSamples),
      paste("SampleID column:          ", SampleID),
      paste("Array type:               ", arrayType),
      paste("Annotation version:       ", annotationVersion),
      paste("Base RData folder:        ", baseDataFolder),
      paste("Base Figure folder:       ", figureBaseDir),
      paste(
        "TIFF size (w x h @ dpi):  ",
        tiffWidth,
        " x ",
        tiffHeight,
        " @ ",
        tiffRes
      ),
      paste("QC cutoff (median):       ", qcCutoff),
      paste("Detection P-value type:   ", detPtype),
      paste("Detection p-value threshold:", detPThreshold),
      paste(
        "Normalization methods:    ",
        paste(norm_method_list, collapse = ", ")
      ),
      paste("Sex column:               ", sexColumn),
      paste("Plot grouping variable:   ", plotGroupVar),
      "Probe filtering:",
      paste("  P-value threshold:      ", pvalThreshold),
      paste("  Chromosomes to remove:  ", chrToRemove),
      paste("  SNP positions filter:   ", snpsToRemove),
      paste("  MAF threshold:          ", mafThreshold),
      paste("  Cross-reactive file:    ", crossReactivePath),
      "Cell composition (estimateLC):",
      paste("  Reference:              ", lcRef),
      paste("  Leading pheno order:    ", phenoOrder),
      paste("Display plots:            ", display),
      paste("Verbose messages:         ", verbose),
      paste("Write logs:               ", logs),
      paste("Save outputs:             ", saveOutputs),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  withLoggedErrorsMinfiEwasWater(
    expr = {
  objectDir <- file.path(baseDataFolder, scriptLabel, "objects")
  normDir <- file.path(baseDataFolder, scriptLabel, "normObjects")
  metricsDir <- file.path(baseDataFolder, scriptLabel, "metrics")
  filterDir <- file.path(baseDataFolder, scriptLabel, "filterObjects")
  qcDir <- file.path(baseDataFolder, scriptLabel, "qc")
  metricsFigDir <- file.path(figureBaseDir, scriptLabel, "metrics")
  qcFigDir <- file.path(figureBaseDir, scriptLabel, "qc")
  enmixDir <- file.path(figureBaseDir, scriptLabel, "enmix")

  if (isTRUE(saveOutputs)) {
    dirs_to_create <- c(
      objectDir,
      normDir,
      metricsDir,
      filterDir,
      qcDir,
      metricsFigDir,
      qcFigDir,
      enmixDir,
      lcPhenoDir
    )
    for (dir_path in dirs_to_create) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
    }
  }

  targets <- readPhenotypeTargets(
    phenoFile = phenoFile,
    sepType = sepType,
    nSamples = nSamples,
    SampleID = SampleID,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  RGSet <- readRGSetMinfiEwasWater(
    idatFolder = idatFolder,
    targets = targets,
    SampleID = SampleID,
    arrayType = arrayType,
    annotationVersion = annotationVersion,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  plotCtrlMinfiEwasWater(
    RGSet = RGSet,
    output_dir = if (isTRUE(saveOutputs)) enmixDir else NULL,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  rawDataInitial <- buildRawMinfiEwasWater(
    RGSet = RGSet,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  assessment <- assessSamplesMinfiEwasWater(
    rawData = rawDataInitial,
    RGSet = RGSet,
    qcCutoff = qcCutoff,
    detPtype = detPtype,
    detPThreshold = detPThreshold,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  plotAssessmentMinfiEwasWater(
    assessment = assessment,
    plot = "qc",
    display = display,
    file = if (isTRUE(saveOutputs)) {
      file.path(qcFigDir, "quality_control(MSet).tiff")
    } else {
      NULL
    },
    width = tiffWidth,
    height = tiffHeight,
    res = tiffRes,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  plotAssessmentMinfiEwasWater(
    assessment = assessment,
    plot = "detection",
    display = display,
    file = if (isTRUE(saveOutputs)) {
      file.path(qcFigDir, "detection_pvalues(RGSet).tiff")
    } else {
      NULL
    },
    width = tiffWidth,
    height = tiffHeight,
    res = tiffRes,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  sampleData <- filterSamplesMinfiEwasWater(
    RGSet = RGSet,
    targets = targets,
    failedSamples = assessment$failedSamples,
    SampleID = SampleID,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  rawData <- buildRawMinfiEwasWater(
    RGSet = sampleData$RGSet,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  plotRawDensityMinfiEwasWater(
    rawData = rawData,
    targets = sampleData$targets,
    plotGroupVar = plotGroupVar,
    display = display,
    file = if (isTRUE(saveOutputs)) {
      file.path(qcFigDir, "densityBeta(MSet).tiff")
    } else {
      NULL
    },
    width = tiffWidth,
    height = tiffHeight,
    res = tiffRes,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  sexData <- predictSexMinfiEwasWater(
    rawData = rawData,
    targets = sampleData$targets,
    SampleID = SampleID,
    sexColumn = sexColumn,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  sampleData$targets <- sexData$targets
  rgset_col_data <- SummarizedExperiment::colData(sampleData$RGSet)
  if (sexColumn %in% colnames(rgset_col_data)) {
    rgset_col_data[[sexColumn]] <- sampleData$targets[[sexColumn]]
  }
  rgset_col_data$PredSex <- sampleData$targets$PredSex
  SummarizedExperiment::colData(sampleData$RGSet) <- rgset_col_data

  plotSexMinfiEwasWater(
    sexData = sexData,
    type = "predicted",
    display = display,
    file = if (isTRUE(saveOutputs)) {
      file.path(qcFigDir, "sexPrediction(GSet).tiff")
    } else {
      NULL
    },
    width = tiffWidth,
    height = tiffHeight,
    res = 70L,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  plotSexMinfiEwasWater(
    sexData = sexData,
    type = "clinical",
    display = display,
    file = if (isTRUE(saveOutputs)) {
      file.path(qcFigDir, "sexClinical(GSet).tiff")
    } else {
      NULL
    },
    width = tiffWidth,
    height = tiffHeight,
    res = 70L,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  normData <- normalizeMinfiEwasWater(
    sampleData = sampleData,
    sexColumn = sexColumn,
    normMethods = normMethods,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  plotNormalizationMinfiEwasWater(
    RGSet = sampleData$RGSet,
    normData = normData,
    targets = sampleData$targets,
    sexColumn = sexColumn,
    display = display,
    file = if (isTRUE(saveOutputs)) {
      file.path(qcFigDir, "sexComparison_RawNorm(MSetF).tiff")
    } else {
      NULL
    },
    width = tiffWidth,
    height = tiffHeight,
    res = tiffRes,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  filterData <- filterProbesMinfiEwasWater(
    normData = normData,
    RGSet = sampleData$RGSet,
    pvalThreshold = pvalThreshold,
    chrToRemove = chrToRemove,
    snpsToRemove = snpsToRemove,
    mafThreshold = mafThreshold,
    crossReactivePath = crossReactivePath,
    detPtype = detPtype,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  metricsData <- extractMetricsMinfiEwasWater(
    filteredData = filterData,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  plotMetricsMinfiEwasWater(
    metricsData = metricsData,
    targets = sampleData$targets,
    plot = "mds",
    plotGroupVar = plotGroupVar,
    sexColumn = sexColumn,
    display = display,
    file = if (isTRUE(saveOutputs)) {
      file.path(
        metricsFigDir,
        "examineMDS_PostFilteringCrossRect(MSetF_Flt_Rxy_Ds_Rc).tiff"
      )
    } else {
      NULL
    },
    width = tiffWidth,
    height = tiffHeight,
    res = tiffRes,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  plotMetricsMinfiEwasWater(
    metricsData = metricsData,
    targets = sampleData$targets,
    plot = "density",
    plotGroupVar = plotGroupVar,
    sexColumn = sexColumn,
    display = display,
    file = if (isTRUE(saveOutputs)) {
      file.path(metricsFigDir, "densityBeta&M(MSetF_Flt_Rxy_Ds_Rc).tiff")
    } else {
      NULL
    },
    width = tiffWidth,
    height = tiffHeight,
    res = tiffRes,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  lcData <- estimateLCMinfiEwasWater(
    beta = metricsData$beta,
    targets = sampleData$targets,
    lcRef = lcRef,
    phenoOrder = phenoOrder,
    constrained = FALSE,
    verbose = verbose,
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )

  if (isTRUE(saveOutputs)) {
    saveNamedObjectMinfiEwasWater(
      sampleData$RGSet,
      "RGSet",
      file.path(objectDir, "RGSet.RData")
    )
    saveNamedObjectMinfiEwasWater(
      rawData$MSet,
      "MSet",
      file.path(objectDir, "MSet.RData")
    )
    saveNamedObjectMinfiEwasWater(
      rawData$RatioSet,
      "RatioSet",
      file.path(objectDir, "RatioSet.RData")
    )
    saveNamedObjectMinfiEwasWater(
      rawData$GSet,
      "GSet",
      file.path(objectDir, "GSet.RData")
    )
    saveNamedObjectMinfiEwasWater(
      assessment$detP,
      "detP",
      file.path(qcDir, "detP_RGSet.RData")
    )

    for (i in seq_along(normData$methods)) {
      saveNamedObjectMinfiEwasWater(
        normData$normalized[[i]],
        "normObj",
        file.path(
          normDir,
          paste0("norm_", normData$methods[[i]], "_RGSet.RData")
        )
      )
    }

    saveNamedObjectMinfiEwasWater(
      filterData$detPFiltered,
      "MSetF_Flt",
      file.path(filterDir, "removProbes_MSetF_Flt.RData")
    )
    saveNamedObjectMinfiEwasWater(
      filterData$chrFiltered,
      "MSetF_Flt_Rxy",
      file.path(filterDir, "removChrXY_MSetF_Flt_Rxy.RData")
    )
    saveNamedObjectMinfiEwasWater(
      filterData$snpFiltered,
      "MSetF_Flt_Rxy_Ds",
      file.path(
        filterDir,
        paste0("removSNPs_MAF", mafThreshold, "_MSetF_Flt_Rxy_Ds.RData")
      )
    )
    saveNamedObjectMinfiEwasWater(
      filterData$filtered,
      "MSetF_Flt_Rxy_Ds_Rc",
      file.path(filterDir, "removCrossReactive_MSetF_Flt_Rxy_Ds_Rc.RData")
    )
    saveNamedObjectMinfiEwasWater(
      metricsData$m,
      "m",
      file.path(metricsDir, "m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
    )
    saveNamedObjectMinfiEwasWater(
      metricsData$beta,
      "beta",
      file.path(metricsDir, "beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
    )
    saveNamedObjectMinfiEwasWater(
      metricsData$cn,
      "cn",
      file.path(metricsDir, "cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
    )

    writePhenoLCMinfiEwasWater(
      lcData = lcData,
      file = file.path(lcPhenoDir, "phenoLC.csv"),
      verbose = verbose,
      logs = logs,
      log_dir = outputLogs,
      log_file = log_file
    )
  }

  emitLogMinfiEwasWater(
    c(
      paste("==== Finished", scriptLabel, "===="),
      paste("End Time:                 ", format(Sys.time())),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      targets = sampleData$targets,
      RGSet = sampleData$RGSet,
      rawData = rawData,
      assessment = assessment,
      sexData = sexData,
      normData = normData,
      filterData = filterData,
      metricsData = metricsData,
      lcData = lcData,
      logFile = log_path
    ),
    class = "dnaEPICO_preprocessingMinfiEwasWater"
  )
    },
    log_path = log_path,
    verbose = verbose,
    context = "preprocessingMinfiEwasWater"
  )
}
