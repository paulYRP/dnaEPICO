glmInputLogLinesDnaEpico <- function(config) {
    c(
        "==== Starting DNAm GLM Analysis ====",
        paste("Start time:                ", format(Sys.time())),
        paste("Input phenotype + methylation:", config$inputPheno),
        paste("Merged modeling object:    ",
            config$methylationObjectPrefix, "*"
        ),
        paste("Output RData folder:       ", config$outputRData),
        paste("Output logs folder:        ", config$outputLogs),
        paste("Output plots folder:       ", config$outputPlots),
        paste("Report assets folder:      ",
            if (is.null(config$reportAssetsDir)) {
                "None"
            } else {
                config$reportAssetsDir
            }
        ),
        paste("Phenotypes:                ", config$phenotypes),
        paste("Covariates:                ", config$covariates),
        paste("Factor variables:          ", config$factorVars),
        paste("Scale variables:           ",
            if (length(normalizeScaleVariablesDnaEpico(config$scaleVars))) {
                paste(normalizeScaleVariablesDnaEpico(config$scaleVars),
                    collapse = ","
                )
            } else {
                "None"
            }
        )
    )
}

glmModelLogLinesDnaEpico <- function(config) {
    c(
        paste("CpG column prefix:         ", config$cpgPrefix),
        paste("CpG limit:                 ",
            if (is.na(config$cpgLimit)) "All" else config$cpgLimit
        ),
        paste("Number of cores:           ", as.integer(config$nCores)),
        paste("Interaction term:          ",
            if (is.null(config$interactionTerm)) {
                "None"
            } else {
                config$interactionTerm
            }
        ),
        paste("Omnibus test:             ", config$omnibusTest),
        paste("Venn coefficient phenotypes:",
            if (is.null(config$vennDPhenotypes)) "None" else {
                paste(config$vennDPhenotypes, collapse = ",")
            }
        ),
        paste("Venn omnibus phenotypes:  ",
            if (is.null(config$vennDOmnibusPhenotypes)) "None" else {
                paste(config$vennDOmnibusPhenotypes, collapse = ",")
            }
        ),
        paste("GLM libraries:             ", config$glmLibs),
        paste("PRS mapping:               ",
            if (is.null(config$prsMap)) "None" else config$prsMap
        ),
        paste("Summary p-value filter:    ",
            if (is.na(config$summaryPval)) "None" else config$summaryPval
        ),
        paste("Include Residual SD:       ",
            isTRUE(config$summaryResidualSD)
        ),
        paste("Chunk size:                ",
            if (is.null(config$chunkSize)) "Auto" else config$chunkSize
        ),
        paste("FDR threshold:             ", config$fdrThreshold),
        paste("P-value adjustment method: ", config$padjmethod),
        paste("Display plots:             ", isTRUE(config$display)),
        paste("Save outputs:              ", isTRUE(config$saveOutputs)),
        "============================================================"
    )
}

logMethylationGLMStartDnaEpico <- function(config) {
    emitLogMinfiEwasWater(
        c(glmInputLogLinesDnaEpico(config), glmModelLogLinesDnaEpico(config)),
        verbose = config$verbose, log_path = config$log_path
    )
}

prepareMethylationGLMWorkflowDnaEpico <- function(config) {
    prepared <- prepareMethylationGLMData(
        inputPheno = config$inputPheno, phenotypes = config$phenotypes,
        covariates = config$covariates, factorVars = config$factorVars,
        scaleVars = config$scaleVars, cpgPrefix = config$cpgPrefix,
        cpgLimit = config$cpgLimit,
        methylationScale = config$methylationScale,
        interactionTerm = config$interactionTerm, prsMap = config$prsMap,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
    output_dir <- if (isTRUE(config$saveOutputs)) config$outputPlots else NULL
    distributions <- plotMethylationGLMDistributions(
        preparedData = prepared, plotWidth = config$plotWidth,
        plotHeight = config$plotHeight, plotDPI = config$plotDPI,
        outputDir = output_dir, display = config$display,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
    design <- plotModelDesignDnaEpico(
        preparedData = prepared, analysis = "GLM", outputDir = output_dir,
        plotWidth = config$plotWidth, plotHeight = config$plotHeight,
        plotDPI = config$plotDPI, display = config$display
    )
    list(
        preparedData = prepared, distributionPlots = distributions,
        designPlots = design
    )
}

fitMethylationGLMWorkflowDnaEpico <- function(prepared, config) {
    summary_dir <- if (isTRUE(config$saveOutputs)) config$outputRData else NULL
    fits <- fitMethylationGLMModels(
        preparedData = prepared, nCores = config$nCores,
        libPath = config$libPath, glmLibs = config$glmLibs,
        omnibusTest = config$omnibusTest, summaryDir = summary_dir,
        resumeFromSummary = config$resumeFromSummary,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
    summaries <- summarizeMethylationGLMModels(
        modelResults = fits, preparedData = prepared,
        summaryResidualSD = config$summaryResidualSD,
        summaryPval = config$summaryPval, padjmethod = config$padjmethod,
        nCores = config$nCores, libPath = config$libPath,
        glmLibs = config$glmLibs, chunkSize = config$chunkSize,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
    significant <- if (!isTRUE(config$saveSignificantCpGs)) NULL else {
        collectSignificantCpGsMethylationGLM(
            modelResults = fits, pvalThreshold = config$significantCpGPval,
            interactionTerm = prepared$interactionTerm,
            verbose = config$verbose, logs = config$logs,
            log_dir = config$outputLogs, log_file = config$log_file
        )
    }
    list(
        modelFits = fits, modelSummaries = summaries,
        significantCpGs = significant
    )
}

diagnoseMethylationGLMWorkflowDnaEpico <- function(stages, config) {
    output_dir <- if (isTRUE(config$saveOutputs)) config$outputPlots else NULL
    diagnostics <- plotMethylationGLMDiagnostics(
        modelSummaries = stages$modelSummaries,
        preparedData = stages$preparedData,
        fdrThreshold = config$fdrThreshold, padjmethod = config$padjmethod,
        outputDir = output_dir, plotWidth = config$plotWidth,
        plotHeight = config$plotHeight, plotDPI = config$plotDPI,
        display = config$display, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$log_file
    )
    annotation <- annotateMethylationGLMSummaries(
        modelSummaries = stages$modelSummaries,
        annotationObject = config$annotationPackage,
        annotationCols = config$annotationCols, gencodeHub = config$gencodeHub,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
    list(diagnosticPlots = diagnostics, annotation = annotation)
}

plotMethylationGLMResultsDnaEpico <- function(stages, config) {
    output_dir <- if (isTRUE(config$saveOutputs)) config$outputPlots else NULL
    manhattan <- plotAnnotatedManhattanDnaEpico(
        annotatedResults = stages$annotation, analysis = "GLM",
        outputDir = output_dir, plotWidth = config$plotWidth,
        plotHeight = config$plotHeight, plotDPI = config$plotDPI,
        display = config$display
    )
    venn <- generateModelVennDDnaEpico(
        annotatedResults = stages$annotation,
        modelSummaries = stages$modelSummaries, analysis = "GLM",
        vennDPhenotypes = config$vennDPhenotypes,
        vennDLabels = config$vennDLabels,
        vennDOmnibusPhenotypes = config$vennDOmnibusPhenotypes,
        vennDOmnibusLabels = config$vennDOmnibusLabels,
        outputDir = output_dir, plotWidth = config$plotWidth,
        plotHeight = config$plotHeight, plotDPI = config$plotDPI,
        display = config$display, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$log_file
    )
    list(manhattanPlots = manhattan, vennDPlots = venn)
}

saveMethylationGLMWorkflowDnaEpico <- function(stages, config) {
    if (!isTRUE(config$saveOutputs)) {
        return(NULL)
    }
    writeMethylationGLMOutputs(
        modelResults = stages$modelFits,
        modelSummaries = stages$modelSummaries,
        annotatedResults = stages$annotation,
        significantCpGs = stages$significantCpGs,
        outputRData = config$outputRData,
        summaryTxtDir = config$summaryTxtDir,
        significantCpGDir = config$significantCpGDir,
        annotatedGLMOut = config$annotatedGLMOut,
        reportAssetsDir = config$reportAssetsDir,
        vennDResults = stages$vennDPlots,
        saveTxtSummaries = config$saveTxtSummaries,
        saveSignificantCpGs = config$saveSignificantCpGs,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
}

newMethylationGLMResultDnaEpico <- function(stages, config, savedFiles) {
    structure(c(stages, list(
        savedFiles = savedFiles,
        runSettings = list(
            analysisLabel = "methylationGLM",
            methylationScale = config$methylationScale,
            methylationLabel = config$methylationLabel,
            methylationObjectPrefix = config$methylationObjectPrefix,
            internalResponseColumn = config$responseColumn,
            scaleVars = stages$preparedData$scaleVars,
            scalingMetadata = stages$preparedData$scalingMetadata,
            omnibusTest = config$omnibusTest,
            vennDPhenotypes = config$vennDPhenotypes,
            vennDOmnibusPhenotypes = config$vennDOmnibusPhenotypes,
            gencodeHub = config$gencodeHub,
            reportAssetsDir = config$reportAssetsDir
        )
    )), class = "dnaEPICO_methylationGLM")
}

runMethylationGLMWorkflowDnaEpico <- function(config) {
    prepared <- prepareMethylationGLMWorkflowDnaEpico(config)
    fitted <- fitMethylationGLMWorkflowDnaEpico(prepared$preparedData, config)
    stages <- c(prepared, fitted)
    stages <- c(stages, diagnoseMethylationGLMWorkflowDnaEpico(stages, config))
    stages <- c(stages, plotMethylationGLMResultsDnaEpico(stages, config))
    saved_files <- saveMethylationGLMWorkflowDnaEpico(stages, config)
    emitLogMinfiEwasWater(c(
        "============================================================",
        paste("Finished DNAm GLM Analysis:", format(Sys.time())),
        if (isTRUE(config$saveOutputs)) {
            paste("Annotated GLM output:         ", saved_files$annotatedGLM)
        } else {
            "Outputs were returned in memory only."
        }, "============================================================"
    ), verbose = config$verbose, log_path = config$log_path)
    newMethylationGLMResultDnaEpico(stages, config, saved_files)
}

normalizeMethylationGLMConfigDnaEpico <- function(config) {
    config$resumeFromSummary <- validateLogicalScalarDnaEpico(
        config$resumeFromSummary, "resumeFromSummary"
    )
    config$gencodeHub <- validateLogicalScalarDnaEpico(
        config$gencodeHub, "gencodeHub"
    )
    config$omnibusTest <-
        validateOmnibusConfigurationMethylationGLM(config$omnibusTest)
    config$cpgLimit <- normalizeOptionalNumericMethylationGLM(config$cpgLimit)
    config$summaryPval <-
        normalizeOptionalNumericMethylationGLM(config$summaryPval)
    config$chunkSize <- normalizeChunkSizeMethylationGLM(config$chunkSize)
    config$methylationScale <-
        normalizeMethylationScaleDnaEpico(config$methylationScale)
    config$methylationLabel <-
        methylationScaleResponseLabelDnaEpico(config$methylationScale)
    config$methylationObjectPrefix <-
        methylationScaleObjectPrefixDnaEpico(config$methylationScale)
    config$responseColumn <-
        methylationScaleResponseColumnDnaEpico(config$methylationScale)
    if (is.null(config$libPath)) {
        config$libPath <- .libPaths()
    }
    config$log_file <- "log_methylationGLM.txt"
    config$log_path <- resolveLogPathMinfiEwasWater(
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$log_file
    )
    config
}

#' Fit CpG-wise GLMs for one-timepoint methylation analyses
#'
#' @param inputPheno Character. Path to the merged phenotype-plus-methylation
#' `.RData`
#'   or `.rds` object created by `preprocessingPheno()`. The default points to
#'   the timepoint-1 object produced by the package workflow.
#' @param outputLogs Character. Directory used for optional log files.
#' @param outputRData Character. Directory used for compact, resumable
#'   phenotype summaries.
#' @param outputPlots Character. Directory used for optional TIFF plots.
#' @param phenotypes Character vector or comma-separated phenotype variables to
#'   model.
#' @param covariates Character. Comma-separated covariate variables included in
#'   each GLM.
#' @param factorVars Character. Comma-separated variables to convert to factors
#'   before modeling.
#' @param scaleVars Character vector, comma-separated variable names, or `NULL`.
#' Numeric fixed-effect variables to centre and divide by their sample standard
#'   deviations before model fitting.
#' @param cpgPrefix Character. Prefix used to identify methylation columns in
#' the
#'   merged phenotype-plus-methylation input object. The default is `'cg'`.
#' @param cpgLimit Integer or `NA`. Maximum number of CpGs to analyse. Use `NA`
#'   to keep all CpGs matching `cpgPrefix`.
#' @param methylationScale Character. Methylation metric represented by the CpG
#'   columns. One of `'Beta'`, `'M'`, or `'CN'`, in any combination of
#'   upper- and lower-case letters. The default is `'beta'`.
#' @param nCores Integer. Maximum number of worker processes to use while
#'   fitting models. Automatic fitting remains serial below the glm2 crossover
#'   and caps workers by the CpG workload, available CPUs, and detected memory.
#' @param plotWidth Integer. TIFF width in pixels when plots are written to
#' disk.
#' @param plotHeight Integer. TIFF height in pixels when plots are written to
#'   disk.
#' @param plotDPI Integer. TIFF resolution in DPI when plots are written to
#'   disk.
#' @param interactionTerm Character or `NULL`. Optional interaction term. When
#'   supplied and present in the input data, the phenotype is modeled together
#'   with its interaction against this variable.
#' @param omnibusTest Logical. If `TRUE`, use `car::linearHypothesis()` to test
#'   the complete phenotype-by-interaction term, or the phenotype main effect
#'   when `interactionTerm = NULL`, once per CpG. One-degree-of-freedom terms
#'   are tested and therefore reproduce the corresponding coefficient p-value.
#' @param vennDPhenotypes Character vector, comma-separated phenotype names, or
#'   `NULL`. Selected phenotypes are expanded to all coefficient p-value
#'   columns for model-level Venn diagrams and workbook tables.
#' @param vennDLabels Character vector, comma-separated display labels, or
#'   `NULL`. Labels follow the resolved coefficient order and preserve case.
#' @param vennDOmnibusPhenotypes Character vector, comma-separated phenotype
#'   names, or `NULL`. These use only omnibus p-value columns and require
#'   `omnibusTest = TRUE`.
#' @param vennDOmnibusLabels Character vector, comma-separated display labels,
#'   or `NULL`, supplied in the resolved omnibus phenotype order.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded
#' to
#'   worker processes. By default, the current `.libPaths()` are used.
#' @param glmLibs Character. Comma-separated package names to validate on worker
#'   processes. The default is `'glm2'`.
#' @param prsMap Character or `NULL`. Optional phenotype-to-PRS mapping in the
#'   form `'Phenotype1:PRS_1,Phenotype2:PRS_2'`.
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
#'   significant CpG coefficient tables. The threshold is applied to omnibus
#'   p-values when `omnibusTest = TRUE`, and to target coefficient p-values
#'   otherwise.
#' @param saveTxtSummaries Logical. If `TRUE` and `saveOutputs = TRUE`, write
#'   tab-delimited summary tables to `summaryTxtDir`.
#' @param chunkSize Integer or `NULL`. Number of CpGs processed per summary
#'   extraction chunk. `NULL` chooses a value automatically.
#' @param summaryTxtDir Character. Directory used for optional tab-delimited GLM
#'   summary tables.
#' @param fdrThreshold Numeric. False-discovery-rate threshold used to highlight
#'   CpGs in the residual-significance diagnostic plots.
#' @param padjmethod Character. P-value adjustment method passed to
#'   `stats::p.adjust()`. The default is `'fdr'`.
#' @param annotationPackage Character. Annotation package or object name passed
#'   to `minfi::getAnnotation()`, for example
#'   `'IlluminaHumanMethylationEPICv2anno.20a1.hg38'`.
#' @param annotationCols Character vector or comma-separated annotation columns
#'   to append to the combined GLM summary table. Available columns depend on
#'   the selected annotation package.
#' @param gencodeHub Logical. If `TRUE`, append release-aware GENCODE gene-body
#'   and nearest-TSS annotations obtained through AnnotationHub. The selected
#'   array annotation must use GRCh38 coordinates.
#' @param annotatedGLMOut Character. Directory used for the optional annotated
#'   GLM summary XLSX workbook.
#' @param reportAssetsDir Character or `NULL`. Report results directory used for
#'   the compressed TSV table and its compact metadata sidecars. `NULL` writes
#'   only the model outputs and annotated workbook.
#' @param display Logical. If `TRUE`, draw exploratory and diagnostic plots on
#'   the active graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#'   The default is `FALSE`.
#' @param logs Logical. If `TRUE`, write the same progress messages to
#'   `file.path(outputLogs, 'log_methylationGLM.txt')`.
#' @param saveOutputs Logical. If `TRUE`, write compact phenotype summaries,
#'   text summaries, significant-CpG tables, annotated results, and TIFF plots.
#'   The default is `FALSE`.
#' @param resumeFromSummary Logical. If `TRUE` and `saveOutputs = TRUE`, reuse a
#'   complete phenotype summary when its input file and model configuration
#'   match the current analysis. If processing stops before a phenotype summary
#'   is complete, that phenotype is fitted again from its first CpG.
#'
#' @return A list with class `'dnaEPICO_methylationGLM'`.
#' \describe{
#'   \item{preparedData}{Object returned by [prepareMethylationGLMData()]
#' containing the merged phenotype-plus-methylation analysis table and modeling
#'   metadata.}
#'   \item{distributionPlots}{Object returned by
#'   [plotMethylationGLMDistributions()] describing any exploratory plots that
#'   were generated or written.}
#'   \item{designPlots}{Missingness and numeric-correlation plots for the
#'   variables used in the model.}
#'   \item{modelFits}{Object returned by [fitMethylationGLMModels()]
#'   containing compact per-phenotype coefficient, omnibus, and condition
#'   results.}
#'   \item{modelSummaries}{Object returned by
#'   [summarizeMethylationGLMModels()] containing the combined CpG summary
#'   tables used for reporting and annotation.}
#'   \item{significantCpGs}{Object returned by
#'   [collectSignificantCpGsMethylationGLM()] containing optional
#'   phenotype-specific significant-CpG tables.}
#'   \item{diagnosticPlots}{Object returned by
#'   [plotMethylationGLMDiagnostics()] describing the diagnostic plot objects
#'   and any written TIFF files.}
#'   \item{annotation}{Object returned by
#'   [annotateMethylationGLMSummaries()] containing the annotated combined
#'   summary table.}
#'   \item{manhattanPlots}{Versioned circular and rectangular Manhattan plots
#'   for every raw p-value column in the annotated results.}
#'   \item{vennDPlots}{Requested coefficient and omnibus model-level Venn plots,
#'   worksheet tables, and label mappings.}
#'   \item{savedFiles}{Object returned by [writeMethylationGLMOutputs()] when
#'   `saveOutputs = TRUE`, otherwise `NULL`.}
#'   \item{runSettings}{High-level run metadata including the generic analysis
#'   label, methylation scale, display label, selected merged-object prefix, and
#'   internal response-column name.}
#' }
#' See [dnaEPICO_methylationGLM-class] for a class-level overview.
#'
#' @description
#' `methylationGLM()` prepares phenotype-plus-methylation data, fits one
#' Gaussian GLM per CpG and phenotype, summarizes and annotates the results, and
#' creates optional coefficient tables and diagnostic plots. It writes outputs
#' only when `saveOutputs = TRUE`.
#' Numeric CpG columns are passed to `glm2::glm2()` without a separate
#' methylation-domain filter. Native model messages, warnings, and errors are
#' recorded in one phenotype-specific `Model.Message` field. Annotated outputs
#' contain CpGs with at least one returned coefficient or omnibus p-value;
#' aggregate availability and condition counts are recorded in workbook
#' metadata.
#'
#' @examples
#' if (requireNamespace(
#'     "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#'     quietly = TRUE
#' )) {
#'     tmp <- tempdir()
#'     toy_path <- file.path(tmp, "phenoBT1.RData")
#'     phenoBT1 <- data.frame(
#'         Sample_Name = c("S1", "S2", "S3", "S4"),
#'         status = factor(c("Case", "Case", "Control", "Control")),
#'         sex = factor(c("F", "M", "F", "M")),
#'         cg00000029 = c(0.20, 0.25, 0.22, 0.27),
#'         cg00000108 = c(0.60, 0.55, 0.52, 0.58),
#'         check.names = FALSE
#'     )
#'     save(phenoBT1, file = toy_path)
#'
#'     result <- methylationGLM(
#'         inputPheno = toy_path,
#'         phenotypes = "status",
#'         covariates = "sex",
#'         factorVars = "status,sex",
#'         cpgLimit = 2,
#'         nCores = 1,
#'         summaryPval = 1,
#'         annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#'         annotationCols = "Name,chr,pos",
#'         display = FALSE,
#'         verbose = FALSE,
#'         logs = FALSE,
#'         saveOutputs = FALSE
#'     )
#'
#'     class(result)
#' }
#'
#' @seealso [dnaEPICO_methylationGLM-class]
#'
#' @export
methylationGLM <- function(
    inputPheno = "rData/preprocessingPheno/mergeData/phenoBT1.RData",
    outputLogs = "logs", outputRData = "rData/methylationGLM/models",
    outputPlots = "figures/methylationGLM",
    phenotypes = c(
        "DASS_Depression", "DASS_Anxiety", "DASS_Stress",
        "PCL5_TotalScore", "MHCSF_TotalScore", "BRS_TotalScore"
    ),
    covariates = paste0(
        "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,",
        "Epithelial.cells"
    ),
    factorVars = "Sex,Ethnicity,TraumaDefinition", scaleVars = NULL,
    cpgPrefix = "cg", cpgLimit = NA, methylationScale = "beta",
    nCores = 32, plotWidth = 2000, plotHeight = 1000, plotDPI = 150,
    interactionTerm = NULL, omnibusTest = FALSE,
    vennDPhenotypes = NULL, vennDLabels = NULL,
    vennDOmnibusPhenotypes = NULL, vennDOmnibusLabels = NULL,
    libPath = NULL, glmLibs = "glm2", prsMap = NULL,
    summaryPval = NA, summaryResidualSD = TRUE,
    saveSignificantCpGs = FALSE,
    significantCpGDir = "preliminaryResults/cpgs/methylationGLM",
    significantCpGPval = 0.05, saveTxtSummaries = TRUE, chunkSize = NULL,
    summaryTxtDir = "preliminaryResults/summary/methylationGLM",
    fdrThreshold = 0.05, padjmethod = "fdr",
    annotationPackage = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    annotationCols = c(
        "Name", "chr", "pos", "UCSC_RefGene_Group",
        "UCSC_RefGene_Name", "Relation_to_Island", "GencodeV41_Group"
    ),
    gencodeHub = FALSE,
    annotatedGLMOut = "data/methylationGLM", reportAssetsDir = NULL,
    display = FALSE, verbose = FALSE, logs = FALSE, saveOutputs = FALSE,
    resumeFromSummary = TRUE
) {
    config <- normalizeMethylationGLMConfigDnaEpico(
        as.list(environment(), all.names = TRUE)
    )
    logMethylationGLMStartDnaEpico(config)
    withLoggedErrorsMinfiEwasWater(
        expr = runMethylationGLMWorkflowDnaEpico(config),
        log_path = config$log_path, verbose = config$verbose,
        context = "methylationGLM"
    )
}
