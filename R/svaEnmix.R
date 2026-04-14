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
svaEnmixValidateRGSet <- function(
    RGSet,
    rgsetData,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_svaEnmix.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )
  rgset_class <- paste(class(RGSet), collapse = ", ")
  rgset_dims <- dim(RGSet)

  emitLogMinfiEwasWater(
    c(
      paste("Loaded RGSet object class: ", rgset_class),
      paste(
        "Loaded RGSet dimensions:   ",
        if (length(rgset_dims) == 2L) {
          paste(rgset_dims, collapse = " x ")
        } else {
          "unavailable"
        }
      )
    ),
    verbose = verbose,
    log_path = log_path
  )

  if (length(rgset_dims) != 2L) {
    stop(
      "The object loaded from ",
      rgsetData,
      " has class ",
      rgset_class,
      " and does not contain a usable RGSet with two dimensions.",
      call. = FALSE
    )
  }

  if (!methods::is(RGSet, "SummarizedExperiment")) {
    stop(
      "The object loaded from ",
      rgsetData,
      " has class ",
      rgset_class,
      " and is not compatible with SummarizedExperiment-based RGSet processing.",
      call. = FALSE
    )
  }

  RGSet
}

#' Estimate surrogate variables from ENmix control probes
#'
#' Read the phenotype table and a saved `RGChannelSet`, estimate surrogate
#' variables from ENmix control probes, analyze their association with Sentrix
#' chip and position factors, and return a structured in-memory result. Legacy
#' CSV, `.RData`, text-summary, and figure outputs are written only when
#' `saveOutputs = TRUE`.
#'
#' @param phenoFile Character. Path to the phenotype file with cell-composition
#'   data.
#' @param rgsetData Character. Path to a saved `RGChannelSet` object. Both
#'   `.RData` and `.rds` files are supported.
#' @param sepType Character. Field separator used in `phenoFile`. Use `""` for
#'   a comma-separated file, `"\\t"` for a tab-delimited file, or another
#'   separator accepted by `utils::read.csv()`.
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
#' @param saveOutputs Logical. If `TRUE`, write the legacy CSV, `.RData`, text,
#'   and TIFF outputs to disk. The default is `FALSE`.
#'
#' @return A list with class `"dnaEPICO_svaEnmix"` containing the phenotype
#'   data, loaded `RGChannelSet`, surrogate-variable matrix, merged phenotype,
#'   association-analysis objects, optional saved-file paths, and the resolved
#'   log file path.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' if (requireNamespace("minfiData", quietly = TRUE)) {
#'   ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#'   pheno_file <- file.path(tmp, "pheno.csv")
#'   rgset_path <- file.path(tmp, "RGSet.RData")
#'   RGSet <- ex$RGSet
#'   utils::write.csv(ex$targets, pheno_file, row.names = FALSE)
#'   save(RGSet, file = rgset_path)
#'   sva_result <- svaEnmix(
#'     phenoFile = pheno_file,
#'     rgsetData = rgset_path,
#'     SampleID = "Sample_Name",
#'     arrayType = "IlluminaHumanMethylation450k",
#'     annotationVersion = "ilmn12.hg19",
#'     SentrixIDColumn = "Sentrix_ID",
#'     SentrixPositionColumn = "Sentrix_Position",
#'     outputLogs = file.path(tmp, "logs"),
#'     figureBaseDir = file.path(tmp, "figures"),
#'     dataBaseDir = file.path(tmp, "data"),
#'     rBaseDir = file.path(tmp, "rData"),
#'     saveOutputs = FALSE
#'   )
#'   stopifnot(inherits(sva_result, "dnaEPICO_svaEnmix"))
#' }
#'
#' @export
svaEnmix <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
    rgsetData = "rData/preprocessingMinfiEwasWater/objects/RGSet.RData",
    sepType = "",
    outputLogs = "logs",
    nSamples = NA,
    SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
    annotationVersion = "20a1.hg38",
    SentrixIDColumn = "Sentrix_ID",
    SentrixPositionColumn = "Sentrix_Position",
    ctrlSvaPercVar = 0.90,
    ctrlSvaFlag = 1,
    scriptLabel = "svaEnmix",
    tiffWidth = 2000,
    tiffHeight = 1000,
    tiffRes = 150,
    figureBaseDir = "figures",
    dataBaseDir = "data",
    rBaseDir = "rData",
    display = FALSE,
    verbose = FALSE,
    logs = FALSE,
    saveOutputs = FALSE
) {
  log_file <- "log_svaEnmix.txt"
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = outputLogs,
    log_file = log_file
  )
  figure_dir <- file.path(figureBaseDir, scriptLabel)

  emitLogMinfiEwasWater(
    c(
      "==== Starting SVA Estimation with Enmix ====",
      paste("Start time:               ", format(Sys.time())),
      paste(
        "Log file path:            ",
        if (is.null(log_path)) "disabled" else log_path
      ),
      paste("Pheno file:               ", phenoFile),
      paste("RGSet path:               ", rgsetData),
      paste(
        "Separator type:           ",
        if (is.null(resolveSeparatorMinfiEwasWater(sepType))) {
          "default (',')"
        } else {
          sepType
        }
      ),
      paste("Sample limit:             ", if (is.na(nSamples)) "all" else nSamples),
      paste("SampleID column:          ", SampleID),
      paste("Array type:               ", arrayType),
      paste("Annotation version:       ", annotationVersion),
      paste("Sentrix ID column:        ", SentrixIDColumn),
      paste("Sentrix position column:  ", SentrixPositionColumn),
      paste("ctrlSva percvar:          ", ctrlSvaPercVar),
      paste("ctrlSva flag:             ", ctrlSvaFlag),
      paste("Script label:             ", scriptLabel),
      paste("TIFF dimensions (WxH):    ", tiffWidth, " x ", tiffHeight, " @ ", tiffRes),
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

      RGSet <- loadSavedObjectPreprocessingPheno(
        path = rgsetData,
        preferred_name = "RGSet"
      )
      RGSet <- svaEnmixValidateRGSet(
        RGSet = RGSet,
        rgsetData = rgsetData,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = log_file
      )
      rgset_ncol <- ncol(RGSet)

      if (rgset_ncol != nrow(targets)) {
        stop(
          "The saved RGSet contains ",
          rgset_ncol,
          " samples but the phenotype table contains ",
          nrow(targets),
          ".",
          call. = FALSE
        )
      }

      colnames(RGSet) <- targets[[SampleID]]
      Biobase::annotation(RGSet) <- c(
        array = arrayType,
        annotation = annotationVersion
      )

      emitLogMinfiEwasWater(
        c(
          paste("RGSet loaded with          ", rgset_ncol, " samples."),
          paste("Applied annotation:        ", paste(Biobase::annotation(RGSet), collapse = ", ")),
          "======================================================================="
        ),
        verbose = verbose,
        log_path = log_path
      )

      svaData <- estimateSvaEnmixControls(
        RGSet = RGSet,
        ctrlSvaPercVar = ctrlSvaPercVar,
        ctrlSvaFlag = ctrlSvaFlag,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = log_file
      )
      mergedPheno <- mergeSvaTargetsEnmix(
        targets = targets,
        sva = svaData$sva,
        SampleID = SampleID,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = log_file
      )
      analysisData <- analyzeSvaEnmix(
        sva = svaData$sva,
        RGSet = RGSet,
        SentrixIDColumn = SentrixIDColumn,
        SentrixPositionColumn = SentrixPositionColumn,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = log_file
      )

      plot_files <- list(
        sentrixID = plotSvaEnmix(
          analysisData = analysisData,
          plot = "sentrix_id",
          display = display,
          file = if (isTRUE(saveOutputs)) file.path(figure_dir, "sva_SentrixID.tiff") else NULL,
          width = tiffWidth,
          height = tiffHeight,
          res = tiffRes,
          verbose = verbose,
          logs = logs,
          log_dir = outputLogs,
          log_file = log_file
        ),
        sentrixPosition = plotSvaEnmix(
          analysisData = analysisData,
          plot = "sentrix_position",
          display = display,
          file = if (isTRUE(saveOutputs)) file.path(figure_dir, "sva_SentrixPosition.tiff") else NULL,
          width = tiffWidth,
          height = tiffHeight,
          res = tiffRes,
          verbose = verbose,
          logs = logs,
          log_dir = outputLogs,
          log_file = log_file
        ),
        matrix = plotSvaEnmix(
          analysisData = analysisData,
          plot = "matrix",
          display = display,
          file = if (isTRUE(saveOutputs)) file.path(figure_dir, "sva_SentrixIDPosition.tiff") else NULL,
          width = tiffWidth,
          height = tiffHeight,
          res = tiffRes,
          verbose = verbose,
          logs = logs,
          log_dir = outputLogs,
          log_file = log_file
        )
      )

      savedFiles <- NULL
      if (isTRUE(saveOutputs)) {
        savedFiles <- writeSvaEnmixOutputs(
          svaData = svaData,
          mergedPheno = mergedPheno,
          analysisData = analysisData,
          phenoFile = phenoFile,
          dataBaseDir = dataBaseDir,
          rBaseDir = rBaseDir,
          scriptLabel = scriptLabel,
          verbose = verbose,
          logs = logs,
          log_dir = outputLogs,
          log_file = log_file
        )
      }

      emitLogMinfiEwasWater(
        c(
          "==== Finished SVA Estimation with Enmix ====",
          paste("End time:                 ", format(Sys.time())),
          "======================================================================="
        ),
        verbose = verbose,
        log_path = log_path
      )

      structure(
        list(
          targets = targets,
          RGSet = RGSet,
          svaData = svaData,
          mergedPheno = mergedPheno,
          analysisData = analysisData,
          plotFiles = plot_files,
          savedFiles = savedFiles,
          logFile = log_path
        ),
        class = "dnaEPICO_svaEnmix"
      )
    },
    log_path = log_path,
    verbose = verbose,
    context = "svaEnmix"
  )
}
