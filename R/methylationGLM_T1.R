#' Fit CpG-wise GLMs for one-timepoint methylation analyses
#'
#' @param inputPheno Character. Path to the merged phenotype-plus-beta `.RData`
#'   or `.rds` object created by `preprocessingPheno()`. The default points to
#'   the timepoint-1 object produced by the package workflow.
#' @param outputLogs Character. Directory used for optional log files.
#' @param outputRData Character. Directory used for optional serialized model and
#'   summary outputs.
#' @param outputPlots Character. Directory used for optional TIFF plots.
#' @param phenotypes Character vector or comma-separated phenotype variables to
#'   model.
#' @param covariates Character. Comma-separated covariate variables included in
#'   each GLM.
#' @param factorVars Character. Comma-separated variables that should be treated
#'   as factors before modeling.
#' @param cpgPrefix Character. Prefix used to identify methylation columns in the
#'   merged phenotype-plus-beta input object. The default is `"cg"`.
#' @param cpgLimit Integer or `NA`. Maximum number of CpGs to analyse. Use `NA`
#'   to keep all CpGs matching `cpgPrefix`.
#' @param nCores Integer. Number of worker processes to use while fitting models
#'   and extracting summaries.
#' @param plotWidth Integer. TIFF width in pixels when plots are written to disk.
#' @param plotHeight Integer. TIFF height in pixels when plots are written to
#'   disk.
#' @param plotDPI Integer. TIFF resolution in DPI when plots are written to
#'   disk.
#' @param interactionTerm Character or `NULL`. Optional interaction term. When
#'   supplied and present in the input data, the phenotype is modeled together
#'   with its interaction against this variable.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded to
#'   worker processes. By default, the current `.libPaths()` are used.
#' @param glmLibs Character. Comma-separated package names to validate on worker
#'   processes. The default is `"glm2"`.
#' @param prsMap Character or `NULL`. Optional phenotype-to-PRS mapping in the
#'   form `"Phenotype1:PRS_1,Phenotype2:PRS_2"`.
#' @param summaryPval Numeric or `NA`. Optional p-value threshold applied to the
#'   returned CpG summary tables. Use `NA` to keep all summary rows.
#' @param summaryResidualSD Logical. If `TRUE`, append residual standard
#'   deviations to the CpG summary tables and residual diagnostic plots.
#' @param saveSignificantCpGs Logical. If `TRUE`, collect significant CpG
#'   coefficient tables in the returned object and optionally write them to disk
#'   when `saveOutputs = TRUE`.
#' @param significantCpGDir Character. Directory used for optional significant
#'   CpG coefficient tables.
#' @param significantCpGPval Numeric. P-value threshold used to collect or write
#'   significant CpG coefficient tables.
#' @param saveTxtSummaries Logical. If `TRUE` and `saveOutputs = TRUE`, write
#'   tab-delimited summary tables to `summaryTxtDir`.
#' @param chunkSize Integer or `NULL`. Number of CpGs processed per summary
#'   extraction chunk. `NULL` chooses a value automatically.
#' @param summaryTxtDir Character. Directory used for optional tab-delimited GLM
#'   summary tables.
#' @param fdrThreshold Numeric. False-discovery-rate threshold used to highlight
#'   CpGs in the residual-significance diagnostic plots.
#' @param padjmethod Character. P-value adjustment method passed to
#'   `stats::p.adjust()`. The default is `"fdr"`.
#' @param annotationPackage Character. Annotation package or object name passed
#'   to `minfi::getAnnotation()`, for example
#'   `"IlluminaHumanMethylationEPICv2anno.20a1.hg38"`.
#' @param annotationCols Character vector or comma-separated annotation columns
#'   to append to the combined GLM summary table. Available columns depend on
#'   the selected annotation package.
#' @param annotatedGLMOut Character. Directory used for the optional annotated
#'   GLM summary CSV file.
#' @param display Logical. If `TRUE`, draw exploratory and diagnostic plots on
#'   the active graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#'   The default is `FALSE`, so the function is quiet unless requested.
#' @param logs Logical. If `TRUE`, write the same progress messages to
#'   `file.path(outputLogs, "log_methylationGLM_T1.txt")`.
#' @param saveOutputs Logical. If `TRUE`, write optional serialized model files,
#'   summary tables, significant-CpG tables, annotated results, and TIFF plots to
#'   the requested output directories. The default is `FALSE`, so the function
#'   returns in-memory results without writing files.
#'
#' @return A list with class `"dnaEPICO_methylationGLM_T1"`.
#' \describe{
#'   \item{preparedData}{Object returned by [prepareMethylationGLM_T1Data()]
#'   containing the merged phenotype-plus-beta analysis table and modeling
#'   metadata.}
#'   \item{distributionPlots}{Object returned by
#'   [plotMethylationGLM_T1Distributions()] describing any exploratory plots that
#'   were generated or written.}
#'   \item{modelFits}{Object returned by [fitMethylationGLM_T1Models()]
#'   containing the per-phenotype CpG model fits.}
#'   \item{modelSummaries}{Object returned by
#'   [summarizeMethylationGLM_T1Models()] containing the combined CpG summary
#'   tables used for reporting and annotation.}
#'   \item{significantCpGs}{Object returned by
#'   [collectSignificantCpGsMethylationGLM_T1()] containing optional
#'   phenotype-specific significant-CpG tables.}
#'   \item{diagnosticPlots}{Object returned by
#'   [plotMethylationGLM_T1Diagnostics()] describing the diagnostic plot objects
#'   and any written TIFF files.}
#'   \item{annotation}{Object returned by
#'   [annotateMethylationGLM_T1Summaries()] containing the annotated combined
#'   summary table.}
#'   \item{savedFiles}{Object returned by [writeMethylationGLM_T1Outputs()] when
#'   `saveOutputs = TRUE`, otherwise `NULL`.}
#' }
#' See [dnaEPICO_methylationGLM_T1-class] for a class-level overview.
#'
#' @description
#' `methylationGLM_T1()` is the high-level coordinator for the one-timepoint GLM
#' stage of the `dnaEPICO` workflow. It prepares the merged phenotype-plus-beta
#' input, optionally creates exploratory plots, fits one Gaussian GLM per CpG for
#' each requested phenotype, extracts CpG-level summaries, optionally collects
#' significant CpG coefficient tables, generates diagnostic plots, annotates the
#' combined summary table, and optionally writes legacy-style outputs to disk.
#' The default behavior is now in-memory and quiet, which makes the function
#' easier to compose with other package functions and more aligned with typical
#' Bioconductor usage.
#'
#' @examples
#' if (requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
#'   tmp <- tempdir()
#'   toy_path <- file.path(tmp, "phenoBetaT1.RData")
#'   phenoBT1 <- data.frame(
#'     Sample_Name = c("S1", "S2", "S3", "S4"),
#'     status = factor(c("Case", "Case", "Control", "Control")),
#'     sex = factor(c("F", "M", "F", "M")),
#'     cg00000029 = c(0.20, 0.25, 0.22, 0.27),
#'     cg00000108 = c(0.60, 0.55, 0.52, 0.58),
#'     check.names = FALSE
#'   )
#'   save(phenoBT1, file = toy_path)
#'
#'   result <- methylationGLM_T1(
#'     inputPheno = toy_path,
#'     phenotypes = "status",
#'     covariates = "sex",
#'     factorVars = "status,sex",
#'     cpgLimit = 2,
#'     nCores = 1,
#'     summaryPval = 1,
#'     annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#'     annotationCols = "Name,chr,pos",
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE,
#'     saveOutputs = FALSE
#'   )
#'
#'   class(result)
#' }
#'
#' @seealso [dnaEPICO_methylationGLM_T1-class]
#'
#' @export
methylationGLM_T1 <- function(
    inputPheno = "rData/preprocessingPheno/mergeData/phenoBetaT1.RData",
    outputLogs = "logs",
    outputRData = "rData/methylationGLM_T1/models",
    outputPlots = "figures/methylationGLM_T1",
    phenotypes = c(
      "DASS_Depression",
      "DASS_Anxiety",
      "DASS_Stress",
      "PCL5_TotalScore",
      "MHCSF_TotalScore",
      "BRS_TotalScore"
    ),
    covariates = "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells",
    factorVars = "Sex,Ethnicity,TraumaDefinition",
    cpgPrefix = "cg",
    cpgLimit = NA,
    nCores = 32,
    plotWidth = 2000,
    plotHeight = 1000,
    plotDPI = 150,
    interactionTerm = NULL,
    libPath = NULL,
    glmLibs = "glm2",
    prsMap = NULL,
    summaryPval = NA,
    summaryResidualSD = TRUE,
    saveSignificantCpGs = FALSE,
    significantCpGDir = "preliminaryResults/cpgs/methylationGLM_T1",
    significantCpGPval = 0.05,
    saveTxtSummaries = TRUE,
    chunkSize = NULL,
    summaryTxtDir = "preliminaryResults/summary/methylationGLM_T1/glm",
    fdrThreshold = 0.05,
    padjmethod = "fdr",
    annotationPackage = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    annotationCols = c(
      "Name",
      "chr",
      "pos",
      "UCSC_RefGene_Group",
      "UCSC_RefGene_Name",
      "Relation_to_Island",
      "GencodeV41_Group"
    ),
    annotatedGLMOut = "data/methylationGLM_T1",
    display = FALSE,
    verbose = FALSE,
    logs = FALSE,
    saveOutputs = FALSE
) {
  cpgLimit <- normalizeOptionalNumericMethylationGLM_T1(cpgLimit)
  summaryPval <- normalizeOptionalNumericMethylationGLM_T1(summaryPval)
  chunkSize <- normalizeChunkSizeMethylationGLM_T1(chunkSize)

  if (is.null(libPath)) {
    libPath <- .libPaths()
  }

  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = outputLogs,
    log_file = "log_methylationGLM_T1.txt"
  )

  emitLogMinfiEwasWater(
    c(
      "==== Starting DNAm GLM Analysis (Timepoint 1) ====",
      paste("Start time:                ", format(Sys.time())),
      paste("Input phenotype + beta:    ", inputPheno),
      paste("Output RData folder:       ", outputRData),
      paste("Output logs folder:        ", outputLogs),
      paste("Output plots folder:       ", outputPlots),
      paste("Phenotypes:                ", phenotypes),
      paste("Covariates:                ", covariates),
      paste("Factor variables:          ", factorVars),
      paste("CpG column prefix:         ", cpgPrefix),
      paste("CpG limit:                 ", if (is.na(cpgLimit)) "All" else cpgLimit),
      paste("Number of cores:           ", as.integer(nCores)),
      paste("Interaction term:          ", if (is.null(interactionTerm)) "None" else interactionTerm),
      paste("GLM libraries:             ", glmLibs),
      paste("PRS mapping:               ", if (is.null(prsMap)) "None" else prsMap),
      paste("Summary p-value filter:    ", if (is.na(summaryPval)) "None" else summaryPval),
      paste("Include Residual SD:       ", isTRUE(summaryResidualSD)),
      paste("Chunk size:                ", if (is.null(chunkSize)) "Auto" else chunkSize),
      paste("FDR threshold:             ", fdrThreshold),
      paste("P-value adjustment method: ", padjmethod),
      paste("Display plots:             ", isTRUE(display)),
      paste("Verbose messages:          ", isTRUE(verbose)),
      paste("Write logs:                ", isTRUE(logs)),
      paste("Save outputs:              ", isTRUE(saveOutputs)),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  withLoggedErrorsMinfiEwasWater(
    expr = {
      preparedData <- prepareMethylationGLM_T1Data(
        inputPheno = inputPheno,
        phenotypes = phenotypes,
        covariates = covariates,
        factorVars = factorVars,
        cpgPrefix = cpgPrefix,
        cpgLimit = cpgLimit,
        interactionTerm = interactionTerm,
        prsMap = prsMap,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = "log_methylationGLM_T1.txt"
      )

      distributionPlots <- plotMethylationGLM_T1Distributions(
        preparedData = preparedData,
        plotWidth = plotWidth,
        plotHeight = plotHeight,
        plotDPI = plotDPI,
        outputDir = if (isTRUE(saveOutputs)) outputPlots else NULL,
        display = display,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = "log_methylationGLM_T1.txt"
      )

      modelFits <- fitMethylationGLM_T1Models(
        preparedData = preparedData,
        nCores = nCores,
        libPath = libPath,
        glmLibs = glmLibs,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = "log_methylationGLM_T1.txt"
      )

      modelSummaries <- summarizeMethylationGLM_T1Models(
        modelResults = modelFits,
        preparedData = preparedData,
        summaryResidualSD = summaryResidualSD,
        summaryPval = summaryPval,
        nCores = nCores,
        libPath = libPath,
        glmLibs = glmLibs,
        chunkSize = chunkSize,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = "log_methylationGLM_T1.txt"
      )

      significantCpGs <- NULL
      if (isTRUE(saveSignificantCpGs)) {
        significantCpGs <- collectSignificantCpGsMethylationGLM_T1(
          modelResults = modelFits,
          pvalThreshold = significantCpGPval,
          interactionTerm = preparedData$interactionTerm,
          verbose = verbose,
          logs = logs,
          log_dir = outputLogs,
          log_file = "log_methylationGLM_T1.txt"
        )
      }

      diagnosticPlots <- plotMethylationGLM_T1Diagnostics(
        modelSummaries = modelSummaries,
        preparedData = preparedData,
        fdrThreshold = fdrThreshold,
        padjmethod = padjmethod,
        outputDir = if (isTRUE(saveOutputs)) outputPlots else NULL,
        plotWidth = plotWidth,
        plotHeight = plotHeight,
        plotDPI = plotDPI,
        display = display,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = "log_methylationGLM_T1.txt"
      )

      annotation <- annotateMethylationGLM_T1Summaries(
        modelSummaries = modelSummaries,
        annotationObject = annotationPackage,
        annotationCols = annotationCols,
        verbose = verbose,
        logs = logs,
        log_dir = outputLogs,
        log_file = "log_methylationGLM_T1.txt"
      )

      savedFiles <- NULL
      if (isTRUE(saveOutputs)) {
        savedFiles <- writeMethylationGLM_T1Outputs(
          modelResults = modelFits,
          modelSummaries = modelSummaries,
          annotatedResults = annotation,
          significantCpGs = significantCpGs,
          outputRData = outputRData,
          summaryTxtDir = summaryTxtDir,
          significantCpGDir = significantCpGDir,
          annotatedGLMOut = annotatedGLMOut,
          saveTxtSummaries = saveTxtSummaries,
          saveSignificantCpGs = saveSignificantCpGs,
          verbose = verbose,
          logs = logs,
          log_dir = outputLogs,
          log_file = "log_methylationGLM_T1.txt"
        )
      }

      emitLogMinfiEwasWater(
        c(
          "=======================================================================",
          paste("Finished DNAm GLM Analysis (Timepoint 1):", format(Sys.time())),
          if (isTRUE(saveOutputs)) {
            paste("Annotated GLM output:         ", savedFiles$annotatedGLM)
          } else {
            "Outputs were returned in memory only."
          },
          "======================================================================="
        ),
        verbose = verbose,
        log_path = log_path
      )

      structure(
        list(
          preparedData = preparedData,
          distributionPlots = distributionPlots,
          modelFits = modelFits,
          modelSummaries = modelSummaries,
          significantCpGs = significantCpGs,
          diagnosticPlots = diagnosticPlots,
          annotation = annotation,
          savedFiles = savedFiles
        ),
        class = "dnaEPICO_methylationGLM_T1"
      )
    },
    log_path = log_path,
    verbose = verbose,
    context = "methylationGLM_T1"
  )
}
