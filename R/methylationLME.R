lmeInputLogLinesDnaEpico <- function(config) {
    c(
        "==== Starting DNAm LME Analysis ====",
        paste("Start time:                     ", format(Sys.time())),
        paste("Input phenotype + methylation:  ", config$inputPheno),
        paste("Merged modeling object:         ",
            config$methylationObjectPrefix, "*"
        ),
        paste("Output RData folder:            ", config$outputRData),
        paste("Output logs folder:             ", config$outputLogs),
        paste("Output plots folder:            ", config$outputPlots),
        paste("Report assets folder:           ",
            if (is.null(config$reportAssetsDir)) {
                "None"
            } else {
                config$reportAssetsDir
            }
        ),
        paste("Person ID variable:             ", config$personVar),
        paste("Time variable:                  ", config$timeVar),
        paste("Phenotypes:                     ", config$phenotypes),
        paste("Covariates:                     ", config$covariates),
        paste("Factor variables:               ", config$factorVars),
        paste("Scale variables:                ",
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

lmeModelLogLinesDnaEpico <- function(config) {
    c(
        paste("LME libraries:                  ", config$lmeLibs),
        paste("Correlation structure:          ", config$correlationStructure),
        paste("Correlation variable:           ",
            if (identical(config$correlationStructure, "none")) {
                "None"
            } else {
                config$correlationVar
            }
        ),
        paste("PRS mapping:                    ",
            if (is.null(config$prsMap)) "None" else config$prsMap
        ),
        paste("CpG column prefix:              ", config$cpgPrefix),
        paste("CpG limit:                      ",
            if (is.na(config$cpgLimit)) "All" else config$cpgLimit
        ),
        paste("Number of cores:                ", as.integer(config$nCores)),
        paste("Summary p-value filter:         ",
            if (is.na(config$summaryPval)) "None" else config$summaryPval
        ),
        paste("Interaction term:               ",
            if (is.null(config$interactionTerm)) {
                "None"
            } else {
                config$interactionTerm
            }
        ),
        paste("Omnibus test:                  ", config$omnibusTest),
        paste("Omnibus denominator DF:        ",
            if (isTRUE(config$omnibusTest)) config$omnibusDdf else "None"
        )
    )
}

lmeOutputLogLinesDnaEpico <- function(config) {
    c(
        paste("Venn coefficient phenotypes:   ",
            if (is.null(config$vennDPhenotypes)) "None" else {
                paste(config$vennDPhenotypes, collapse = ",")
            }
        ),
        paste("Venn omnibus phenotypes:       ",
            if (is.null(config$vennDOmnibusPhenotypes)) "None" else {
                paste(config$vennDOmnibusPhenotypes, collapse = ",")
            }
        ),
        paste("Save significant interactions:  ",
            isTRUE(config$saveSignificantInteractions)
        ),
        paste("Significant interaction p-value:",
            config$significantInteractionPval
        ),
        paste("Save text summaries:            ",
            isTRUE(config$saveTxtSummaries)
        ),
        paste("Chunk size:                     ",
            if (is.null(config$chunkSize)) "Auto" else config$chunkSize
        ),
        paste("FDR threshold:                  ", config$fdrThreshold),
        paste("P-value adjustment method:      ", config$padjmethod),
        paste("Display plots:                  ", isTRUE(config$display)),
        paste("Save outputs:                   ", isTRUE(config$saveOutputs)),
        "============================================================"
    )
}

logMethylationLMEStartDnaEpico <- function(config) {
    emitLogMinfiEwasWater(c(
        lmeInputLogLinesDnaEpico(config), lmeModelLogLinesDnaEpico(config),
        lmeOutputLogLinesDnaEpico(config)
    ), verbose = config$verbose, log_path = config$log_path)
}

prepareMethylationLMEWorkflowDnaEpico <- function(config) {
    prepareMethylationLMEData(
        inputPheno = config$inputPheno, personVar = config$personVar,
        timeVar = config$timeVar, phenotypes = config$phenotypes,
        covariates = config$covariates, factorVars = config$factorVars,
        scaleVars = config$scaleVars, prsMap = config$prsMap,
        cpgPrefix = config$cpgPrefix, cpgLimit = config$cpgLimit,
        methylationScale = config$methylationScale,
        interactionTerm = config$interactionTerm,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
}

fitMethylationLMEWorkflowDnaEpico <- function(prepared, config) {
    summary_dir <- if (isTRUE(config$saveOutputs)) config$outputRData else NULL
    fitMethylationLMEModels(
        preparedData = prepared, nCores = config$nCores,
        libPath = config$libPath, lmeLibs = config$lmeLibs,
        summaryDir = summary_dir,
        resumeFromSummary = config$resumeFromSummary,
        correlationStructure = config$correlationStructure,
        correlationVar = config$correlationVar,
        omnibusTest = config$omnibusTest, omnibusDdf = config$omnibusDdf,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
}

plotMethylationLMEInputsDnaEpico <- function(prepared, config) {
    output_dir <- if (isTRUE(config$saveOutputs)) config$outputPlots else NULL
    distributions <- plotMethylationLMEDistributions(
        preparedData = prepared, outputDir = output_dir,
        plotWidth = config$plotWidth, plotHeight = config$plotHeight,
        plotDPI = config$plotDPI, display = config$display
    )
    design <- plotModelDesignDnaEpico(
        preparedData = prepared, analysis = "LME", outputDir = output_dir,
        plotWidth = config$plotWidth, plotHeight = config$plotHeight,
        plotDPI = config$plotDPI, display = config$display
    )
    list(distributionPlots = distributions, designPlots = design)
}

summarizeMethylationLMEWorkflowDnaEpico <- function(fits, prepared, config) {
    summaries <- summarizeMethylationLMEModels(
        modelResults = fits, preparedData = prepared,
        summaryPval = config$summaryPval, padjmethod = config$padjmethod,
        nCores = config$nCores, chunkSize = config$chunkSize,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
    significant <- if (!isTRUE(config$saveSignificantInteractions)) NULL else {
        collectSignificantInteractionsMethylationLME(
            modelResults = fits,
            pvalThreshold = config$significantInteractionPval,
            interactionTerm = prepared$interactionTerm,
            verbose = config$verbose, logs = config$logs,
            log_dir = config$outputLogs, log_file = config$log_file
        )
    }
    list(
        modelSummaries = summaries,
        significantInteractions = significant
    )
}

diagnoseMethylationLMEWorkflowDnaEpico <- function(stages, config) {
    output_dir <- if (isTRUE(config$saveOutputs)) config$outputPlots else NULL
    diagnostics <- plotMethylationLMEDiagnostics(
        modelSummaries = stages$modelSummaries,
        preparedData = stages$preparedData,
        fdrThreshold = config$fdrThreshold, padjmethod = config$padjmethod,
        outputDir = output_dir, plotWidth = config$plotWidth,
        plotHeight = config$plotHeight, plotDPI = config$plotDPI,
        display = config$display, verbose = config$verbose,
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$log_file
    )
    annotation <- annotateMethylationLMESummaries(
        modelSummaries = stages$modelSummaries,
        annotationObject = config$annotationPackage,
        annotationCols = config$annotationCols, gencodeHub = config$gencodeHub,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
    list(diagnosticPlots = diagnostics, annotation = annotation)
}

plotMethylationLMEResultsDnaEpico <- function(stages, config) {
    output_dir <- if (isTRUE(config$saveOutputs)) config$outputPlots else NULL
    manhattan <- plotAnnotatedManhattanDnaEpico(
        annotatedResults = stages$annotation, analysis = "LME",
        outputDir = output_dir, plotWidth = config$plotWidth,
        plotHeight = config$plotHeight, plotDPI = config$plotDPI,
        display = config$display
    )
    venn <- generateModelVennDDnaEpico(
        annotatedResults = stages$annotation,
        modelSummaries = stages$modelSummaries, analysis = "LME",
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

saveMethylationLMEWorkflowDnaEpico <- function(stages, config) {
    if (!isTRUE(config$saveOutputs)) {
        return(NULL)
    }
    writeMethylationLMEOutputs(
        modelResults = stages$modelFits,
        modelSummaries = stages$modelSummaries,
        annotatedResults = stages$annotation,
        significantInteractions = stages$significantInteractions,
        outputRData = config$outputRData,
        summaryTxtDir = config$summaryTxtDir,
        significantInteractionDir = config$significantInteractionDir,
        annotatedLMEOut = config$annotatedLMEOut,
        reportAssetsDir = config$reportAssetsDir,
        vennDResults = stages$vennDPlots,
        saveTxtSummaries = config$saveTxtSummaries,
        saveSignificantInteractions = config$saveSignificantInteractions,
        verbose = config$verbose, logs = config$logs,
        log_dir = config$outputLogs, log_file = config$log_file
    )
}

newMethylationLMEResultDnaEpico <- function(stages, config, savedFiles) {
    structure(c(stages, list(
        savedFiles = savedFiles,
        runSettings = list(
            analysisLabel = "methylationLME",
            methylationScale = config$methylationScale,
            methylationLabel = config$methylationLabel,
            methylationObjectPrefix = config$methylationObjectPrefix,
            internalResponseColumn = config$responseColumn,
            scaleVars = stages$preparedData$scaleVars,
            scalingMetadata = stages$preparedData$scalingMetadata,
            timeVar = config$timeVar, omnibusTest = config$omnibusTest,
            omnibusDdf = config$omnibusDdf,
            vennDPhenotypes = config$vennDPhenotypes,
            vennDOmnibusPhenotypes = config$vennDOmnibusPhenotypes,
            gencodeHub = config$gencodeHub,
            reportAssetsDir = config$reportAssetsDir,
            correlationStructure = config$correlationStructure,
            correlationVar = if (identical(
                config$correlationStructure, "none"
            )) NULL else config$correlationVar
        )
    )), class = "dnaEPICO_methylationLME")
}

runMethylationLMEWorkflowDnaEpico <- function(config) {
    prepared <- prepareMethylationLMEWorkflowDnaEpico(config)
    fits <- fitMethylationLMEWorkflowDnaEpico(prepared, config)
    plots <- plotMethylationLMEInputsDnaEpico(prepared, config)
    summarized <- summarizeMethylationLMEWorkflowDnaEpico(
        fits, prepared, config
    )
    stages <- c(list(preparedData = prepared), plots,
        list(modelFits = fits), summarized
    )
    stages <- c(stages, diagnoseMethylationLMEWorkflowDnaEpico(stages, config))
    stages <- c(stages, plotMethylationLMEResultsDnaEpico(stages, config))
    saved_files <- saveMethylationLMEWorkflowDnaEpico(stages, config)
    emitLogMinfiEwasWater(c(
        "============================================================",
        paste("Finished DNAm LME Analysis:", format(Sys.time())),
        if (isTRUE(config$saveOutputs)) {
            paste("Annotated LME output:          ", saved_files$annotatedLME)
        } else {
            "Outputs were returned in memory only."
        }, "============================================================"
    ), verbose = config$verbose, log_path = config$log_path)
    newMethylationLMEResultDnaEpico(stages, config, saved_files)
}

normalizeMethylationLMEConfigDnaEpico <- function(config) {
    config$resumeFromSummary <- validateLogicalScalarDnaEpico(
        config$resumeFromSummary, "resumeFromSummary"
    )
    config$gencodeHub <- validateLogicalScalarDnaEpico(
        config$gencodeHub, "gencodeHub"
    )
    config$cpgLimit <- normalizeOptionalNumericMethylationGLM(config$cpgLimit)
    config$summaryPval <-
        normalizeOptionalNumericMethylationGLM(config$summaryPval)
    config$chunkSize <- normalizeChunkSizeMethylationGLM(config$chunkSize)
    config$padjmethod <-
        validatePAdjustmentMethodMethylationModels(config$padjmethod)
    config$correlationStructure <-
        normalizeCorrelationStructureMethylationLME(
            config$correlationStructure
        )
    config$correlationVar <- normalizeCorrelationVariableMethylationLME(
        correlationVar = config$correlationVar
    )
    lme_config <- resolveLmeLibrariesMethylationLME(config$lmeLibs)
    omnibus <- validateOmnibusConfigurationMethylationLME(
        omnibusTest = config$omnibusTest, omnibusDdf = config$omnibusDdf,
        lmeEngine = lme_config$engine
    )
    config$omnibusTest <- omnibus$test
    config$omnibusDdf <- omnibus$ddf
    if (!identical(config$correlationStructure, "none") &&
        is.null(config$correlationVar)) {
        stop("correlationVar must be supplied when correlationStructure ",
            "is AR1 or CAR1.", call. = FALSE
        )
    }
    config$methylationScale <-
        normalizeMethylationScaleDnaEpico(config$methylationScale)
    config$methylationLabel <-
        methylationScaleResponseLabelDnaEpico(config$methylationScale)
    config$methylationObjectPrefix <-
        methylationScaleObjectPrefixDnaEpico(config$methylationScale)
    config$responseColumn <-
        methylationScaleResponseColumnDnaEpico(config$methylationScale)
    if (is.null(config$libPath)) config$libPath <- .libPaths()
    config$log_file <- "log_methylationLME.txt"
    config$log_path <- resolveLogPathMinfiEwasWater(
        logs = config$logs, log_dir = config$outputLogs,
        log_file = config$log_file
    )
    config
}

#' Fit CpG-wise linear mixed-effects models for longitudinal methylation
#' analyses
#'
#' @param inputPheno Character. Path to the merged longitudinal
#'   phenotype-plus-methylation `.RData` or `.rds` object created by
#'   `preprocessingPheno()`. The default points to the combined timepoint object
#'   produced by the package workflow.
#' @param outputLogs Character. Directory used for optional log files.
#' @param outputRData Character. Directory used for compact, resumable
#'   phenotype summaries.
#' @param outputPlots Character. Directory used for optional TIFF diagnostic
#' plots.
#' @param personVar Character. Subject identifier variable used for the random
#'   intercept. When this column is missing, it is derived from `SID` using the
#'   package's existing sample naming convention.
#' @param timeVar Character. Name of the longitudinal time variable used for
#'   timepoint summaries and preprocessing checks.
#' @param phenotypes Character vector or comma-separated phenotype variables to
#'   model.
#' @param covariates Character. Comma-separated fixed-effect covariates included
#'   in every mixed model.
#' @param factorVars Character. Comma-separated variables to convert to factors
#'   before modeling, including categorical phenotypes, covariates, or
#'   interaction variables.
#' @param scaleVars Character vector, comma-separated variable names, or `NULL`.
#' Numeric fixed-effect variables to centre and divide by their sample standard
#'   deviations before fitting.
#' @param lmeLibs Character. Comma-separated package names to validate on worker
#'   processes and select the LME backend. Use `'lme4,lmerTest'` or `'lme4'`
#'   for the `lmerTest`/`lme4` path, or `'nlme'` for the `nlme::lme()` path.
#' @param correlationStructure Character. Residual correlation structure used
#'   when `lmeLibs = 'nlme'`. One of `'none'`, `'AR1'`, or `'CAR1'`. The
#'   default is `'none'`.
#' @param correlationVar Character or `NULL`. Variable used to order repeated
#'   observations within `personVar` for `AR1` or `CAR1` residual correlation
#'   structures. Must be supplied explicitly for `AR1` or `CAR1`.
#' @param prsMap Character or `NULL`. Optional phenotype-to-PRS mapping in the
#'   form `'Phenotype1:PRS_1,Phenotype2:PRS_2'`.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded
#' to
#'   worker processes. By default, the current `.libPaths()` are used.
#' @param cpgPrefix Character. Prefix used to identify methylation columns in
#' the
#'   merged phenotype-plus-methylation input object. The default is `'cg'`.
#' @param cpgLimit Integer or `NA`. Maximum number of CpGs to analyse. Use `NA`
#'   to keep all CpGs matching `cpgPrefix`.
#' @param methylationScale Character. Methylation metric represented by the CpG
#'   columns. One of `'Beta'`, `'M'`, or `'CN'`, in any combination of
#'   upper- and lower-case letters. The default is `'beta'`.
#' @param nCores Integer. Maximum number of worker processes to use while
#'   fitting models. Automatic fitting uses an engine-specific crossover and
#'   caps workers by the CpG workload, available CPUs, and detected memory.
#' @param summaryPval Numeric or `NA`. Optional p-value threshold applied to the
#' returned longitudinal CpG summary tables. Use `NA` to keep all summary rows.
#' @param plotWidth Integer. TIFF width in pixels when plots are written to
#' disk.
#' @param plotHeight Integer. TIFF height in pixels when plots are written to
#' disk.
#' @param plotDPI Integer. TIFF resolution in DPI when plots are written to
#' disk.
#' @param interactionTerm Character or `NULL`. Optional interaction term. When
#'   supplied and present in the input data, the phenotype is modeled together
#'   with its interaction against this variable.
#' @param omnibusTest Logical. If `TRUE`, use `lmerTest::contestMD()` to test
#'   the complete phenotype-by-interaction term, or the phenotype main effect
#'   when `interactionTerm = NULL`, once per CpG. This is available only for
#'   the lmerTest/lme4 engine.
#' @param omnibusDdf Character. Denominator degrees-of-freedom method for the
#'   omnibus F test. One of `'Satterthwaite'` or `'Kenward-Roger'`.
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
#' @param saveSignificantInteractions Logical. If `TRUE`, collect coefficient
#'   tables for CpGs passing `significantInteractionPval` in the returned object
#'   and optionally write them to disk when `saveOutputs = TRUE`.
#' @param significantInteractionDir Character. Directory used for optional
#'   significant-interaction coefficient tables.
#' @param significantInteractionPval Numeric. P-value threshold used to collect
#' or
#'   write significant interaction coefficient tables.
#' @param saveTxtSummaries Logical. If `TRUE` and `saveOutputs = TRUE`, write
#'   tab-delimited summary tables to `summaryTxtDir`.
#' @param chunkSize Integer or `NULL`. Number of CpGs processed per summary
#'   extraction chunk. `NULL` chooses a value automatically.
#' @param summaryTxtDir Character. Directory used for optional tab-delimited LME
#'   summary tables.
#' @param fdrThreshold Numeric. False-discovery-rate threshold used to highlight
#'   CpGs in the residual-significance diagnostic plots.
#' @param padjmethod Character. P-value adjustment method passed to
#'   `stats::p.adjust()`. The default is `'fdr'`.
#' @param annotationPackage Character. Annotation package or object name passed
#' to
#'   `minfi::getAnnotation()`, for example
#'   `'IlluminaHumanMethylationEPICv2anno.20a1.hg38'`.
#' @param annotationCols Character vector or comma-separated annotation columns
#'   to append to the combined LME summary table. Available columns depend on
#'   the selected annotation package.
#' @param gencodeHub Logical. If `TRUE`, append release-aware GENCODE gene-body
#'   and nearest-TSS annotations obtained through AnnotationHub. The selected
#'   array annotation must use GRCh38 coordinates.
#' @param annotatedLMEOut Character. Directory used for the optional annotated
#' LME
#'   summary XLSX workbook.
#' @param reportAssetsDir Character or `NULL`. Report results directory used for
#'   the compressed TSV table and its compact metadata sidecars. `NULL` writes
#'   only the model outputs and annotated workbook.
#' @param display Logical. If `TRUE`, draw diagnostic plots on the active
#'   graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#'   The default is `FALSE`.
#' @param logs Logical. If `TRUE`, write the same progress messages to
#'   `file.path(outputLogs, 'log_methylationLME.txt')`.
#' @param saveOutputs Logical. If `TRUE`, write compact phenotype summaries,
#'   text summaries, significant-interaction tables, annotated results, and
#'   TIFF plots. The default is `FALSE`.
#' @param resumeFromSummary Logical. If `TRUE` and `saveOutputs = TRUE`, reuse a
#'   complete phenotype summary when its input file and model configuration
#'   match the current analysis. If processing stops before a phenotype summary
#'   is complete, that phenotype is fitted again from its first CpG.
#'
#' @return A list with class `'dnaEPICO_methylationLME'`.
#' \describe{
#'   \item{preparedData}{Object returned by [prepareMethylationLMEData()]
#' containing the merged longitudinal phenotype-plus-methylation analysis table
#' and
#'   modeling metadata.}
#'   \item{distributionPlots}{Distributions of the model variables,
#'   longitudinal observation counts, time values, and numeric phenotype
#'   trajectories.}
#'   \item{designPlots}{Missingness and numeric-correlation plots for the
#'   variables used in the model.}
#'   \item{modelFits}{Object returned by [fitMethylationLMEModels()]
#'   containing compact per-phenotype coefficient, omnibus, and condition
#'   results.}
#'   \item{modelSummaries}{Object returned by
#'   [summarizeMethylationLMEModels()] containing the combined CpG summary
#'   tables used for reporting and annotation.}
#'   \item{significantInteractions}{Object returned by
#'   [collectSignificantInteractionsMethylationLME()] containing optional
#'   phenotype-specific significant-interaction tables.}
#'   \item{diagnosticPlots}{Object returned by
#'   [plotMethylationLMEDiagnostics()] describing the diagnostic plot
#'   objects and any written TIFF files.}
#'   \item{annotation}{Object returned by
#'   [annotateMethylationLMESummaries()] containing the annotated combined
#'   summary table.}
#'   \item{manhattanPlots}{Versioned circular and rectangular Manhattan plots
#'   for every raw p-value column in the annotated results.}
#'   \item{vennDPlots}{Requested coefficient and omnibus model-level Venn plots,
#'   worksheet tables, and label mappings.}
#'   \item{savedFiles}{Object returned by
#'   [writeMethylationLMEOutputs()] when `saveOutputs = TRUE`, otherwise
#'   `NULL`.}
#'   \item{runSettings}{High-level run metadata including the generic analysis
#'   label, methylation scale, display label, selected merged-object prefix,
#'   internal response-column name, and longitudinal time variable.}
#' }
#' See [dnaEPICO_methylationLME-class] for a class-level overview.
#'
#' @description
#' `methylationLME()` prepares longitudinal phenotype-plus-methylation data,
#' fits one mixed-effects model per CpG and phenotype, summarizes and annotates
#' the results, and creates optional interaction tables and diagnostic plots. It
#' writes outputs only when `saveOutputs = TRUE`.
#' Numeric CpG columns are passed to the selected lmerTest/lme4 or nlme engine
#' without a separate methylation-domain filter. Native model messages,
#' warnings, and errors are recorded in one phenotype-specific `Model.Message`
#' field. Annotated outputs contain CpGs with at least one returned coefficient
#' or omnibus p-value; aggregate availability and condition counts are recorded
#' in workbook metadata.
#'
#' @examples
#' if (
#'     requireNamespace(
#'         "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#'         quietly = TRUE
#'     ) &&
#'         requireNamespace("lmerTest", quietly = TRUE)
#' ) {
#'     tmp <- tempdir()
#'     toy_path <- file.path(tmp, "phenoBT1T2.RData")
#'     phenoBT1T2 <- data.frame(
#'         SID = c("P1A", "P1B", "P2A", "P2B", "P3A", "P3B", "P4A", "P4B"),
#'         person = c(1, 1, 2, 2, 3, 3, 4, 4),
#'         Timepoint = factor(c("1", "2", "1", "2", "1", "2", "1", "2")),
#'         score = c(10, 12, 9, 11, 13, 14, 8, 9),
#'         sex = factor(c("F", "F", "M", "M", "F", "F", "M", "M")),
#'         cg00000029 = c(0.25, 0.27, 0.20, 0.22, 0.30, 0.31, 0.18, 0.20),
#'         cg00000108 = c(0.50, 0.53, 0.55, 0.57, 0.48, 0.49, 0.60, 0.61),
#'         check.names = FALSE
#'     )
#'     save(phenoBT1T2, file = toy_path)
#'
#'     result <- methylationLME(
#'         inputPheno = toy_path,
#'         phenotypes = "score",
#'         covariates = "sex",
#'         factorVars = "sex",
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
#' @seealso [dnaEPICO_methylationLME-class]
#'
#' @export
methylationLME <- function(
    inputPheno = "rData/preprocessingPheno/mergeData/phenoBT1T2.RData",
    outputLogs = "logs", outputRData = "rData/methylationLME/models",
    outputPlots = "figures/methylationLME", personVar = "person",
    timeVar = "Timepoint", phenotypes = c(
        "DASS_Depression", "DASS_Anxiety", "DASS_Stress",
        "PCL5_TotalScore", "MHCSF_TotalScore", "BRS_TotalScore"
    ),
    covariates = paste0(
        "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,",
        "Epithelial.cells"
    ),
    factorVars = "Sex,Ethnicity,TraumaDefinition", scaleVars = NULL,
    lmeLibs = "lme4,lmerTest", correlationStructure = "none",
    correlationVar = NULL, prsMap = NULL, libPath = NULL,
    cpgPrefix = "cg", cpgLimit = NA, methylationScale = "beta",
    nCores = 32, summaryPval = NA, plotWidth = 2000,
    plotHeight = 1000, plotDPI = 150, interactionTerm = NULL,
    omnibusTest = FALSE, omnibusDdf = "Satterthwaite",
    vennDPhenotypes = NULL, vennDLabels = NULL,
    vennDOmnibusPhenotypes = NULL, vennDOmnibusLabels = NULL,
    saveSignificantInteractions = TRUE,
    significantInteractionDir =
        "preliminaryResults/cpgs/methylationLME",
    significantInteractionPval = 0.05, saveTxtSummaries = TRUE,
    chunkSize = NULL,
    summaryTxtDir = "preliminaryResults/summary/methylationLME",
    fdrThreshold = 0.05, padjmethod = "fdr",
    annotationPackage =
        "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    annotationCols = c(
        "Name", "chr", "pos", "UCSC_RefGene_Group",
        "UCSC_RefGene_Name", "Relation_to_Island",
        "GencodeV41_Group"
    ),
    gencodeHub = FALSE, annotatedLMEOut = "data/methylationLME",
    reportAssetsDir = NULL, display = FALSE, verbose = FALSE,
    logs = FALSE, saveOutputs = FALSE, resumeFromSummary = TRUE
) {
    config <- normalizeMethylationLMEConfigDnaEpico(
        as.list(environment(), all.names = TRUE)
    )
    logMethylationLMEStartDnaEpico(config)
    withLoggedErrorsMinfiEwasWater(
        expr = runMethylationLMEWorkflowDnaEpico(config),
        log_path = config$log_path,
        verbose = config$verbose,
        context = "methylationLME"
    )
}
