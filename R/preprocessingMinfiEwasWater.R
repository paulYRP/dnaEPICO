#' Run preprocessingMinfiEwasWater.R as an external script from inst/script/.R
#' Run preprocessingMinfiEwasWater as an external pipeline script
#'
#' Wrapper function that executes the preprocessingMinfiEwasWater.R script
#' shipped in inst/scripts using system calls.
#'
#' @param phenoFile Character. Path to phenotype CSV file.
#' @param idatFolder Character. Path to IDAT files directory.
#' @param outputLogs Character. Directory for log files.
#' @param nSamples Integer or NA. Number of samples to subset for testing.
#' @param SampleID Character. Sample identifier column name.
#' @param arrayType Character. Illumina array type.
#' @param annotationVersion Character. Annotation version.
#' @param scriptLabel Character. Label used in output naming.
#' @param baseDataFolder Character. Base directory for RData outputs.
#' @param figureBaseDir Character. Base directory for Figures outputs.
#' @param sepType Character. Field separator for phenotype file.
#' @param tiffWidth Integer. Width of TIFF plots in pixels.
#' @param tiffHeight Integer. Height of TIFF plots in pixels.
#' @param tiffRes Integer. Resolution (DPI) for TIFF plots.
#' @param qcCutoff Numeric. Quality-control cutoff threshold.
#' @param detPtype Character. Detection p-value calculation type.
#' @param detPThreshold Numeric. Detection p-value threshold.
#' @param funnormSeed Integer. Random seed for normalization.
#' @param normMethods Character. Normalization method(s).
#' @param sexColumn Character. Sex column name in phenotype data.
#' @param pvalThreshold Numeric. Probe-level p-value threshold.
#' @param chrToRemove Character. Chromosomes to remove.
#' @param snpsToRemove Character. SNP probe types to remove.
#' @param mafThreshold Numeric. Minor allele frequency threshold.
#' @param crossReactivePath Character. Path to cross-reactive probe file.
#' @param plotGroupVar Character. Variable used for grouping plots.
#' @param lcRef Character. Reference panel for cell composition.
#' @param phenoOrder Character. Semicolon-separated phenotype column order.
#' @param lcPhenoDir Character. Output directory for cell composition phenotype.
#'
#' @return
#' Invisibly returns \code{NULL}. This function is called for its side effects,
#' performing Illumina EPICv2 preprocessing, quality control, normalisation,
#' probe filtering, cell composition estimation, and writing plots, logs,
#' and RData objects to disk.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' preprocessingMinfiEwasWater(
#'   phenoFile = "data/preprocessingMinfiEwasWater/pheno.csv",
#'   idatFolder = "data/preprocessingMinfiEwasWater/idats",
#'   outputLogs = "logs",
#'   nSamples = NA,
#'   SampleID = "Sample_Name",
#'   arrayType = "IlluminaHumanMethylationEPICv2",
#'   annotationVersion = "20a1.hg38",
#'   scriptLabel = "preprocessingMinfiEwasWater",
#'   baseDataFolder = "rData",
#'   figureBaseDir = "figures",
#'   sepType = "",
#'   tiffWidth = 2000,
#'   tiffHeight = 1000,
#'   tiffRes = 150,
#'   qcCutoff = 10.5,
#'   detPtype = "m+u",
#'   detPThreshold = 0.05,
#'   funnormSeed = 123,
#'   normMethods = "adjustedfunnorm",
#'   sexColumn = "Sex",
#'   pvalThreshold = 0.01,
#'   chrToRemove = "chrX,chrY",
#'   snpsToRemove = "SBE,CpG",
#'   mafThreshold = 0.1,
#'   crossReactivePath = "data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv",
#'   plotGroupVar = "Sex",
#'   lcRef = "salivaEPIC",
#'   phenoOrder = "Sample_Name;Timepoint;Sex;PredSex;Basename;Sentrix_ID;Sentrix_Position",
#'   lcPhenoDir = "data/preprocessingMinfiEwasWater"
#' )
#' }
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
    funnormSeed = 123,
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
    lcPhenoDir = "data/preprocessingMinfiEwasWater"
) {

  # Resolve script inside the package
  script <- system.file("scripts", "preprocessingMinfiEwasWater.R", package = "dnaEPICO")

  if (script == "") stop("Script preprocessingMinfiEwasWater.R not found in package.")

  # Build command-line arguments
  arg_list <- c(
    "--phenoFile", shQuote(phenoFile),
    "--sepType", shQuote(sepType),
    "--idatFolder", shQuote(idatFolder),
    "--outputLogs", shQuote(outputLogs),
    "--nSamples", nSamples,
    "--SampleID", shQuote(SampleID),
    "--arrayType", shQuote(arrayType),
    "--annotationVersion", shQuote(annotationVersion),
    "--scriptLabel", shQuote(scriptLabel),
    "--baseDataFolder", shQuote(baseDataFolder),
    "--figureBaseDir", shQuote(figureBaseDir),
    "--tiffWidth", tiffWidth,
    "--tiffHeight", tiffHeight,
    "--tiffRes", tiffRes,
    "--qcCutoff", qcCutoff,
    "--detPtype", shQuote(detPtype),
    "--detPThreshold", detPThreshold,
    "--funnormSeed", funnormSeed,
    "--normMethods", shQuote(normMethods),
    "--sexColumn", shQuote(sexColumn),
    "--pvalThreshold", pvalThreshold,
    "--chrToRemove", shQuote(chrToRemove),
    "--snpsToRemove", shQuote(snpsToRemove),
    "--mafThreshold", mafThreshold,
    "--crossReactivePath", shQuote(crossReactivePath),
    "--plotGroupVar", shQuote(plotGroupVar),
    "--lcRef", shQuote(lcRef),
    "--phenoOrder", shQuote(phenoOrder),
    "--lcPhenoDir", shQuote(lcPhenoDir)
  )

  # Build system command
  cmd <- paste("Rscript", shQuote(script), paste(arg_list, collapse = " "))

  message("Running preprocessingMinfiEwasWater:")
  message(cmd)

  # Platform specific execution
  if (.Platform$OS.type == "windows") {
    shell(cmd)
  } else {
    system(cmd)
  }
}
