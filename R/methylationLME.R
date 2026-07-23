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
#' @param factorVars Character. Comma-separated variables that should be coerced
#'   to factors before modeling. This usually includes categorical phenotypes,
#'   covariates, or interaction variables.
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
#' @param annotatedLMEOut Character. Directory used for the optional annotated
#' LME
#'   summary XLSX workbook.
#' @param reportAssetsDir Character or `NULL`. Report results directory used for
#'   the compressed TSV table and its compact metadata sidecars. `NULL` writes
#'   only the model outputs and annotated workbook.
#' @param display Logical. If `TRUE`, draw diagnostic plots on the active
#'   graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' The
#'   default is `FALSE`, so the function is quiet unless requested.
#' @param logs Logical. If `TRUE`, write the same progress messages to
#'   `file.path(outputLogs, 'log_methylationLME.txt')`.
#' @param saveOutputs Logical. If `TRUE`, write compact phenotype summaries,
#'   text summaries, significant-interaction tables, annotated results, and TIFF
#'   plots to the requested output directories. The default is `FALSE`, so the
#'   function returns in-memory results without writing files.
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
#' `methylationLME()` is the high-level coordinator for the longitudinal
#' linear mixed-effects stage of the `dnaEPICO` workflow. It prepares the merged
#' phenotype-plus-methylation input, fits one mixed-effects model per CpG for
#' each
#' requested phenotype, extracts phenotype-specific coefficient summaries,
#' optionally collects significant interaction tables, generates diagnostic
#' plots,
#' annotates the combined summary table, and optionally writes analysis outputs
#' to disk. By default, the function runs quietly and returns its results in
#' memory.
#' Numeric CpG columns are passed to the selected lmerTest/lme4 or nlme engine
#' without a separate methylation-domain filter. Native model messages,
#' warnings, and errors are recorded in one phenotype-specific `Model.Message`
#' field. Annotated outputs contain CpGs with at least one returned coefficient
#' or omnibus p-value; aggregate availability and condition counts are recorded
#' in workbook metadata.
#'
#' @examples
#' if (
#'     requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE) &&
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
        "DASS_Depression",
        "DASS_Anxiety", "DASS_Stress", "PCL5_TotalScore", "MHCSF_TotalScore",
        "BRS_TotalScore"
    ), covariates = "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells",
    factorVars = "Sex,Ethnicity,TraumaDefinition", scaleVars = NULL,
    lmeLibs = "lme4,lmerTest", correlationStructure = "none",
    correlationVar = NULL, prsMap = NULL, libPath = NULL, cpgPrefix = "cg",
    cpgLimit = NA, methylationScale = "beta", nCores = 32, summaryPval = NA,
    plotWidth = 2000,
    plotHeight = 1000, plotDPI = 150, interactionTerm = NULL,
    omnibusTest = FALSE, omnibusDdf = "Satterthwaite",
    saveSignificantInteractions = TRUE,
        significantInteractionDir = "preliminaryResults/cpgs/methylationLME",
    significantInteractionPval = 0.05, saveTxtSummaries = TRUE,
    chunkSize = NULL,
        summaryTxtDir = "preliminaryResults/summary/methylationLME",
    fdrThreshold = 0.05, padjmethod = "fdr",
        annotationPackage = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    annotationCols = c(
        "Name", "chr", "pos", "UCSC_RefGene_Group",
        "UCSC_RefGene_Name", "Relation_to_Island", "GencodeV41_Group"
    ),
    annotatedLMEOut = "data/methylationLME", reportAssetsDir = NULL,
    display = FALSE,
    verbose = FALSE, logs = FALSE, saveOutputs = FALSE,
    resumeFromSummary = TRUE
) {
    resumeFromSummary <- validateLogicalScalarDnaEpico(
        resumeFromSummary,
        "resumeFromSummary"
    )
    cpgLimit <- normalizeOptionalNumericMethylationGLM(cpgLimit)
    summaryPval <- normalizeOptionalNumericMethylationGLM(summaryPval)
    chunkSize <- normalizeChunkSizeMethylationGLM(chunkSize)
    padjmethod <- validatePAdjustmentMethodMethylationModels(padjmethod)
    correlationStructure <-
        normalizeCorrelationStructureMethylationLME(correlationStructure)
    correlationVar <-
        normalizeCorrelationVariableMethylationLME(correlationVar = correlationVar)
    lme_config <- resolveLmeLibrariesMethylationLME(lmeLibs)
    omnibus_config <- validateOmnibusConfigurationMethylationLME(
        omnibusTest = omnibusTest, omnibusDdf = omnibusDdf,
        lmeEngine = lme_config$engine
    )
    omnibusTest <- omnibus_config$test
    omnibusDdf <- omnibus_config$ddf
    if (!identical(correlationStructure, "none")) {
        if (is.null(correlationVar)) {
            stop("correlationVar must be supplied when correlationStructure is AR1 or CAR1.",
                call. = FALSE
            )
        }
    }
    methylationScale <- normalizeMethylationScaleDnaEpico(methylationScale)
    methylationLabel <- methylationScaleResponseLabelDnaEpico(methylationScale)
    methylationObjectPrefix <-
        methylationScaleObjectPrefixDnaEpico(methylationScale)
    responseColumn <- methylationScaleResponseColumnDnaEpico(methylationScale)
    log_file <- "log_methylationLME.txt"

    if (is.null(libPath)) {
        libPath <- .libPaths()
    }

    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = outputLogs,
        log_file = log_file
    )

    emitLogMinfiEwasWater(
        c(
            "==== Starting DNAm LME Analysis ====",
            paste("Start time:                     ", format(Sys.time())),
            paste("Input phenotype + methylation:  ", inputPheno),
            paste(
                "Merged modeling object:         ", methylationObjectPrefix,
                "*"
            ), paste("Output RData folder:            ", outputRData),
            paste("Output logs folder:             ", outputLogs),
            paste("Output plots folder:            ", outputPlots),
            paste("Report assets folder:           ",
                if (is.null(reportAssetsDir)) "None" else reportAssetsDir),
            paste("Person ID variable:             ", personVar),
            paste("Time variable:                  ", timeVar), paste(
                "Phenotypes:                     ",
                phenotypes
            ), paste(
                "Covariates:                     ",
                covariates
            ), paste(
                "Factor variables:               ",
                factorVars
            ), paste(
                "Scale variables:                ",
                if (length(normalizeScaleVariablesDnaEpico(scaleVars))) {
                    paste(normalizeScaleVariablesDnaEpico(scaleVars),
                        collapse = ","
                    )
                } else {
                    "None"
                }
            ), paste("LME libraries:                  ", lmeLibs),
            paste("Correlation structure:          ", correlationStructure),
            paste("Correlation variable:           ", if (identical(
                correlationStructure,
                "none"
            )) {
                "None"
            } else {
                correlationVar
            }), paste(
                "PRS mapping:                    ",
                if (is.null(prsMap)) "None" else prsMap
            ), paste(
                "CpG column prefix:              ",
                cpgPrefix
            ), paste(
                "CpG limit:                      ",
                if (is.na(cpgLimit)) "All" else cpgLimit
            ), paste(
                "Number of cores:                ",
                as.integer(nCores)
            ), paste(
                "Summary p-value filter:         ",
                if (is.na(summaryPval)) "None" else summaryPval
            ),
            paste("Interaction term:               ",
                if (is.null(interactionTerm)) "None" else interactionTerm),
            paste("Omnibus test:                  ", omnibusTest),
            paste("Omnibus denominator DF:        ",
                if (isTRUE(omnibusTest)) omnibusDdf else "None"),
            paste("Save significant interactions:  ",
                isTRUE(saveSignificantInteractions)),
            paste("Significant interaction p-value:",
                significantInteractionPval),
            paste("Save text summaries:            ", isTRUE(saveTxtSummaries)),
            paste("Chunk size:                     ",
                if (is.null(chunkSize)) "Auto" else chunkSize),
            paste("FDR threshold:                  ", fdrThreshold),
            paste("P-value adjustment method:      ", padjmethod),
            paste("Display plots:                  ", isTRUE(display)),
            paste("Save outputs:                   ", isTRUE(saveOutputs)),
            "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    withLoggedErrorsMinfiEwasWater(expr = {
        preparedData <- prepareMethylationLMEData(
            inputPheno = inputPheno,
            personVar = personVar, timeVar = timeVar, phenotypes = phenotypes,
            covariates = covariates, factorVars = factorVars,
            scaleVars = scaleVars, prsMap = prsMap, cpgPrefix = cpgPrefix,
            cpgLimit = cpgLimit, methylationScale = methylationScale,
            interactionTerm = interactionTerm, verbose = verbose,
            logs = logs, log_dir = outputLogs, log_file = log_file
        )

        modelFits <- fitMethylationLMEModels(
            preparedData = preparedData,
            nCores = nCores, libPath = libPath, lmeLibs = lmeLibs,
            summaryDir = if (isTRUE(saveOutputs)) outputRData else NULL,
            resumeFromSummary = resumeFromSummary,
            correlationStructure = correlationStructure,
                correlationVar = correlationVar,
            omnibusTest = omnibusTest, omnibusDdf = omnibusDdf,
            verbose = verbose, logs = logs, log_dir = outputLogs,
            log_file = log_file
        )

        modelSummaries <- summarizeMethylationLMEModels(
            modelResults = modelFits,
            preparedData = preparedData, summaryPval = summaryPval,
            padjmethod = padjmethod,
            nCores = nCores, chunkSize = chunkSize, verbose = verbose,
            logs = logs, log_dir = outputLogs, log_file = log_file
        )

        significantInteractions <- NULL
        if (isTRUE(saveSignificantInteractions)) {
            significantInteractions <-
                collectSignificantInteractionsMethylationLME(
                modelResults = modelFits,
                pvalThreshold = significantInteractionPval,
                    interactionTerm = preparedData$interactionTerm,
                verbose = verbose, logs = logs, log_dir = outputLogs,
                log_file = log_file
            )
        }

        diagnosticPlots <- plotMethylationLMEDiagnostics(
            modelSummaries = modelSummaries,
            preparedData = preparedData, fdrThreshold = fdrThreshold,
            padjmethod = padjmethod, outputDir = if (isTRUE(saveOutputs)) {
                outputPlots
            } else {
                NULL
            }, plotWidth = plotWidth, plotHeight = plotHeight,
            plotDPI = plotDPI, display = display, verbose = verbose,
            logs = logs, log_dir = outputLogs, log_file = log_file
        )

        annotation <- annotateMethylationLMESummaries(
            modelSummaries = modelSummaries,
            annotationObject = annotationPackage,
                annotationCols = annotationCols,
            verbose = verbose, logs = logs, log_dir = outputLogs,
            log_file = log_file
        )

        savedFiles <- NULL
        if (isTRUE(saveOutputs)) {
            savedFiles <- writeMethylationLMEOutputs(
                modelResults = modelFits,
                modelSummaries = modelSummaries, annotatedResults = annotation,
                significantInteractions = significantInteractions,
                outputRData = outputRData, summaryTxtDir = summaryTxtDir,
                significantInteractionDir = significantInteractionDir,
                annotatedLMEOut = annotatedLMEOut,
                    reportAssetsDir = reportAssetsDir,
                    saveTxtSummaries = saveTxtSummaries,
                saveSignificantInteractions = saveSignificantInteractions,
                verbose = verbose, logs = logs, log_dir = outputLogs,
                log_file = log_file
            )
        }

        emitLogMinfiEwasWater(
            c(
                "============================================================",
                paste("Finished DNAm LME Analysis:", format(Sys.time())),
                if (isTRUE(saveOutputs)) {
                    paste("Annotated LME output:          ",
                        savedFiles$annotatedLME)
                } else {
                    "Outputs were returned in memory only."
                }, "============================================================"
            ),
            verbose = verbose, log_path = log_path
        )

        structure(list(
            preparedData = preparedData, modelFits = modelFits,
            modelSummaries = modelSummaries,
                significantInteractions = significantInteractions,
            diagnosticPlots = diagnosticPlots, annotation = annotation,
            savedFiles = savedFiles, runSettings = list(
                analysisLabel = "methylationLME",
                methylationScale = methylationScale,
                    methylationLabel = methylationLabel,
                methylationObjectPrefix = methylationObjectPrefix,
                internalResponseColumn = responseColumn,
                    scaleVars = preparedData$scaleVars,
                scalingMetadata = preparedData$scalingMetadata,
                timeVar = timeVar,
                omnibusTest = omnibusTest, omnibusDdf = omnibusDdf,
                reportAssetsDir = reportAssetsDir,
                correlationStructure = correlationStructure,
                correlationVar = if (identical(
                    correlationStructure,
                    "none"
                )) {
                    NULL
                } else {
                    correlationVar
                }
            )
        ), class = "dnaEPICO_methylationLME")
    }, log_path = log_path, verbose = verbose, context = "methylationLME")
}
