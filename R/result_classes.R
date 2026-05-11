#' Result class returned by preprocessingMinfiEwasWater
#'
#' Objects of class `"dnaEPICO_preprocessingMinfiEwasWater"` are list-based
#' results returned by [preprocessingMinfiEwasWater()]. They are lightweight
#' S3-style containers rather than formal S4 classes.
#'
#' @section Structure:
#' \describe{
#'   \item{targets}{Filtered phenotype table aligned to the retained samples.}
#'   \item{RGSet}{Filtered `RGChannelSet` used in downstream preprocessing.}
#'   \item{rawData}{Object returned by [buildRawMinfiEwasWater()].}
#'   \item{assessment}{Object returned by [assessSamplesMinfiEwasWater()].}
#'   \item{sexData}{Object returned by [predictSexMinfiEwasWater()].}
#'   \item{normData}{Object returned by [normalizeMinfiEwasWater()].}
#'   \item{filterData}{Object returned by [filterProbesMinfiEwasWater()].}
#'   \item{metricsData}{Object returned by [extractMetricsMinfiEwasWater()].}
#'   \item{lcData}{Object returned by [estimateLCMinfiEwasWater()].}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#'
#' @seealso [preprocessingMinfiEwasWater()]
#' @name dnaEPICO_preprocessingMinfiEwasWater-class
#' @aliases dnaEPICO_preprocessingMinfiEwasWater
#' @docType class
NULL

#' Result class returned by svaEnmix
#'
#' Objects of class `"dnaEPICO_svaEnmix"` are list-based results returned by
#' [svaEnmix()]. They collect the loaded inputs, surrogate-variable results,
#' association-analysis summaries, and optional file outputs.
#'
#' @section Structure:
#' \describe{
#'   \item{targets}{Phenotype table read from `phenoFile` after any optional row
#'   subsetting.}
#'   \item{RGSet}{Loaded `RGChannelSet` with sample names reset to match
#'   `targets`.}
#'   \item{svaData}{Object returned by [estimateSvaEnmixControls()].}
#'   \item{mergedPheno}{Phenotype table returned by [mergeSvaTargetsEnmix()]
#'   after surrogate variables were appended.}
#'   \item{analysisData}{Object returned by [analyzeSvaEnmix()].}
#'   \item{plotFiles}{Named list of TIFF output paths for the SVA figures when
#'   `saveOutputs = TRUE`, otherwise `NULL` entries for plots not written.}
#'   \item{savedFiles}{Object returned by [writeSvaEnmixOutputs()] when
#'   `saveOutputs = TRUE`, otherwise `NULL`.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#'
#' @seealso [svaEnmix()]
#' @name dnaEPICO_svaEnmix-class
#' @aliases dnaEPICO_svaEnmix
#' @docType class
NULL

#' Result class returned by preprocessingPheno
#'
#' Objects of class `"dnaEPICO_preprocessingPheno"` are list-based results
#' returned by [preprocessingPheno()]. They describe the phenotype data,
#' methylation matrices, timepoint splits, longitudinal merges, and optional
#' exported files.
#'
#' @section Structure:
#' \describe{
#'   \item{pheno}{Phenotype table read from `phenoFile`.}
#'   \item{metricsData}{Object returned by [loadMetricsPreprocessingPheno()].}
#'   \item{timepointData}{Object returned by [splitTimepointsPreprocessingPheno()].}
#'   \item{combinedData}{Object returned by [combineTimepointsPreprocessingPheno()].}
#'   \item{clockFoundation}{Object returned by
#'   [buildClockFoundationInputsPreprocessingPheno()].}
#'   \item{savedFiles}{Object returned by [writePreprocessingPhenoOutputs()]
#'   when `saveOutputs = TRUE`, otherwise `NULL`.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#'
#' @seealso [preprocessingPheno()]
#' @name dnaEPICO_preprocessingPheno-class
#' @aliases dnaEPICO_preprocessingPheno
#' @docType class
NULL

#' Result class returned by methylationGLM_T1
#'
#' Objects of class `"dnaEPICO_methylationGLM_T1"` are list-based results
#' returned by [methylationGLM_T1()]. They collect the prepared analysis table,
#' fitted models, summaries, diagnostics, annotations, and optional saved files.
#'
#' @section Structure:
#' \describe{
#'   \item{preparedData}{Object returned by [prepareMethylationGLM_T1Data()].}
#'   \item{distributionPlots}{Object returned by
#'   [plotMethylationGLM_T1Distributions()].}
#'   \item{modelFits}{Object returned by [fitMethylationGLM_T1Models()].}
#'   \item{modelSummaries}{Object returned by
#'   [summarizeMethylationGLM_T1Models()].}
#'   \item{significantCpGs}{Object returned by
#'   [collectSignificantCpGsMethylationGLM_T1()].}
#'   \item{diagnosticPlots}{Object returned by
#'   [plotMethylationGLM_T1Diagnostics()].}
#'   \item{annotation}{Object returned by
#'   [annotateMethylationGLM_T1Summaries()].}
#'   \item{savedFiles}{Object returned by [writeMethylationGLM_T1Outputs()] when
#'   `saveOutputs = TRUE`, otherwise `NULL`.}
#' }
#'
#' @seealso [methylationGLM_T1()]
#' @name dnaEPICO_methylationGLM_T1-class
#' @aliases dnaEPICO_methylationGLM_T1
#' @docType class
NULL

#' Result class returned by methylationGLMM_T1T2
#'
#' Objects of class `"dnaEPICO_methylationGLMM_T1T2"` are list-based results
#' returned by [methylationGLMM_T1T2()]. They collect the prepared longitudinal
#' analysis table, fitted mixed models, summaries, diagnostics, annotations, and
#' optional saved files.
#'
#' @section Structure:
#' \describe{
#'   \item{preparedData}{Object returned by [prepareMethylationGLMM_T1T2Data()].}
#'   \item{modelFits}{Object returned by [fitMethylationGLMM_T1T2Models()].}
#'   \item{modelSummaries}{Object returned by
#'   [summarizeMethylationGLMM_T1T2Models()].}
#'   \item{significantInteractions}{Object returned by
#'   [collectSignificantInteractionsMethylationGLMM_T1T2()].}
#'   \item{diagnosticPlots}{Object returned by
#'   [plotMethylationGLMM_T1T2Diagnostics()].}
#'   \item{annotation}{Object returned by
#'   [annotateMethylationGLMM_T1T2Summaries()].}
#'   \item{savedFiles}{Object returned by
#'   [writeMethylationGLMM_T1T2Outputs()] when `saveOutputs = TRUE`, otherwise
#'   `NULL`.}
#' }
#'
#' @seealso [methylationGLMM_T1T2()]
#' @name dnaEPICO_methylationGLMM_T1T2-class
#' @aliases dnaEPICO_methylationGLMM_T1T2
#' @docType class
NULL

#' Result class returned by prepareDnamReportInputs
#'
#' Objects of class `"dnaEPICO_dnamReport_prepared"` are list-based results
#' returned by [prepareDnamReportInputs()]. They capture normalized report paths,
#' available figures, and logging metadata before rendering.
#'
#' @section Structure:
#' \describe{
#'   \item{output}{Requested output file name.}
#'   \item{outputDir}{Normalized report output directory.}
#'   \item{outputFile}{Normalized full path to the intended report output file.}
#'   \item{figDir}{Normalized directory used by the report template for copied
#'   figures.}
#'   \item{figureInventory}{Named list describing the available figures for each
#'   report section.}
#'   \item{missingFigureDirectories}{Character vector of expected figure
#'   directories that were not present at preparation time.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#'
#' @seealso [prepareDnamReportInputs()]
#' @name dnaEPICO_dnamReport_prepared-class
#' @aliases dnaEPICO_dnamReport_prepared
#' @docType class
NULL

#' Result class returned by renderDnamReport
#'
#' Objects of class `"dnaEPICO_dnamReport_render"` are list-based results
#' returned by [renderDnamReport()]. They describe whether a prepared report was
#' rendered, skipped, or failed.
#'
#' @section Structure:
#' \describe{
#'   \item{preparedReport}{The input object supplied to [renderDnamReport()].}
#'   \item{status}{Render status string such as `"rendered"`, `"skipped"`, or
#'   `"failed"`.}
#'   \item{renderedFile}{Normalized path to the rendered PDF file when rendering
#'   succeeded, otherwise `NULL`.}
#'   \item{errorMessage}{Render error or skip message when available, otherwise
#'   `NULL`.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#'
#' @seealso [renderDnamReport()]
#' @name dnaEPICO_dnamReport_render-class
#' @aliases dnaEPICO_dnamReport_render
#' @docType class
NULL

#' Result class returned by dnamReport
#'
#' Objects of class `"dnaEPICO_dnamReport"` are list-based results returned by
#' [dnamReport()]. They combine the prepared report inputs, render result, and
#' final status metadata into one convenience object.
#'
#' @section Structure:
#' \describe{
#'   \item{preparedReport}{Object returned by [prepareDnamReportInputs()].}
#'   \item{renderResult}{Structured render metadata created by [dnamReport()].}
#'   \item{status}{Final status string such as `"rendered"`, `"skipped"`, or
#'   `"failed"`.}
#'   \item{outputFile}{Path to `docs/index.html`.}
#'   \item{errorMessage}{Final render error or skip message when available,
#'   otherwise `NULL`.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#'
#' @seealso [dnamReport()]
#' @name dnaEPICO_dnamReport-class
#' @aliases dnaEPICO_dnamReport
#' @docType class
NULL
