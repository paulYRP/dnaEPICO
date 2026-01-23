#' Run svaEnmix as an external pipeline script
#'
#' Wrapper function that executes the svaEnmix.R script located in
#' inst/scripts using a system call.
#'
#' @param phenoFile Character. Path to phenotype file with cell composition data.
#' @param rgsetData Character. Path to RGSet RData file.
#' @param sepType Character. Field separator for phenotype file.
#' @param outputLogs Character. Directory for log files.
#' @param nSamples Integer or NA. Number of samples to subset for testing.
#' @param SampleID Character. Sample identifier column name.
#' @param arrayType Character. Illumina array type.
#' @param annotationVersion Character. Annotation version.
#' @param SentrixIDColumn Character. Sentrix ID column name.
#' @param SentrixPositionColumn Character. Sentrix position column name.
#' @param ctrlSvaPercVar Numeric. Proportion of variance explained by control probes.
#' @param ctrlSvaFlag Integer. Flag indicating use of control probes.
#' @param scriptLabel Character. Label used in output naming.
#' @param tiffWidth Integer. Width of TIFF plots in pixels.
#' @param tiffHeight Integer. Height of TIFF plots in pixels.
#' @param tiffRes Integer. Resolution (DPI) for TIFF plots.
#' @param figureBaseDir Character. Base directory for Figures outputs.
#' @param dataBaseDir Character. Base directory for Data outputs.
#'
#' @return
#' Invisibly returns \code{NULL}. This function is called for its side effects,
#' executing the external \code{svaEnmix.R} script and writing results, figures,
#' and logs to disk.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' svaEnmix(
#'   phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
#'   rgsetData = "rData/preprocessingMinfiEwasWater/objects/RGSet.RData",
#'   sepType = "",
#'   outputLogs = "logs",
#'   nSamples = 5,
#'   SampleID = "Sample_Name",
#'   arrayType = "IlluminaHumanMethylationEPICv2",
#'   annotationVersion = "20a1.hg38",
#'   SentrixIDColumn = "Sentrix_ID",
#'   SentrixPositionColumn = "Sentrix_Position",
#'   ctrlSvaPercVar = 0.90,
#'   ctrlSvaFlag = 1,
#'   scriptLabel = "svaEnmix",
#'   tiffWidth = 2000,
#'   tiffHeight = 1000,
#'   tiffRes = 150,
#'   figureBaseDir = "figures",
#'   dataBaseDir = "data"
#' )
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
    dataBaseDir = "data"
) {

  # Locate script inside installed package
  script <- system.file("scripts", "svaEnmix.R", package = "dnaEPICO")
  if (script == "")
    stop("Script svaEnmix.R not found in package.")

  # Build argument list
  arg_list <- c(
    "--phenoFile", shQuote(phenoFile),
    "--rgsetData", shQuote(rgsetData),
    "--sepType", shQuote(sepType),
    "--outputLogs", shQuote(outputLogs),
    "--nSamples", nSamples,
    "--SampleID", shQuote(SampleID),
    "--arrayType", shQuote(arrayType),
    "--annotationVersion", shQuote(annotationVersion),
    "--SentrixIDColumn", shQuote(SentrixIDColumn),
    "--SentrixPositionColumn", shQuote(SentrixPositionColumn),
    "--ctrlSvaPercVar", ctrlSvaPercVar,
    "--ctrlSvaFlag", ctrlSvaFlag,
    "--scriptLabel", shQuote(scriptLabel),
    "--tiffWidth", tiffWidth,
    "--tiffHeight", tiffHeight,
    "--tiffRes", tiffRes,
    "--figureBaseDir", shQuote(figureBaseDir),
    "--dataBaseDir", shQuote(dataBaseDir)
  )

  # Build system command
  cmd <- paste("Rscript", shQuote(script), paste(arg_list, collapse = " "))

  message("Running svaEnmix:")
  message(cmd)

  # Platform specific execution
  if (.Platform$OS.type == "windows") {
    shell(cmd)
  } else {
    system(cmd)
  }
}
