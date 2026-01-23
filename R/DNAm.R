#' Generate a DNA methylation PDF report
#'
#' Wrapper function that runs DNAm.R and renders DNAm.Rmd from inst/scripts
#' to generate a PDF report summarising preprocessing, QC, SVA, GLM, and GLMM
#' results from the dnaEPICO pipeline.
#'
#' @param output Character. Name of the output PDF file.
#' @param outputDir Character. Directory where the report will be saved.
#' @param qcDir Character. Directory containing ENmix QC figures.
#' @param preprocessingDir Character. Directory containing preprocessing QC figures.
#' @param postprocessingDir Character. Directory containing postprocessing metric figures.
#' @param svaDir Character. Directory containing SVA figures.
#' @param glmDir Character. Directory containing GLM figures.
#' @param glmmDir Character. Directory containing GLMM figures.
#' @param reportTitle Character. Title of the report.
#' @param author Character. Author name displayed in the report.
#' @param date Character. Report date.
#'
#' @return
#' Invisibly returns \code{NULL}. This function is called for its side effect
#' of generating a PDF report and associated output files.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' dnamReport(
#'   output = "DNAm_Report.pdf",
#'   outputDir = "reports",
#'   qcDir = "figures/preprocessingMinfiEwasWater/enMix",
#'   preprocessingDir = "figures/preprocessingMinfiEwasWater/qc",
#'   postprocessingDir = "figures/preprocessingMinfiEwasWater/metrics",
#'   svaDir = "figures/svaEnmix/sva",
#'   glmDir = "figures/methylationGLM_T1",
#'   glmmDir = "figures/methylationGLMM_T1T2",
#'   reportTitle = "DNA methylation analysis",
#'   author = "School of Biomedical Sciences",
#'   date = format(Sys.Date(), "%B %d, %Y")
#' )
#' }
#'
#' @export
dnamReport <- function(
    output = "DNAm_Report.pdf",
    outputDir = "reports",

    qcDir = "figures/preprocessingMinfiEwasWater/enMix",
    preprocessingDir = "figures/preprocessingMinfiEwasWater/qc",
    postprocessingDir = "figures/preprocessingMinfiEwasWater/metrics",
    svaDir = "figures/svaEnmix/sva",
    glmDir = "figures/methylationGLM_T1",
    glmmDir = "figures/methylationGLMM_T1T2",
    reportTitle = "DNA methylation",
    author = "School of Biomedical Sciences",
    date = format(Sys.Date(), "%B %d, %Y")
) {

  # Helper: Convert to absolute path
  makeAbs <- function(path) {
    if (grepl("^([A-Za-z]:|/)", path)) return(path)
    return(file.path(getwd(), path))
  }

  # Locate Rmd + Runner Script
  rmd <- system.file("scripts", "DNAm.Rmd", package = "dnaEPICO")
  script <- system.file("scripts", "DNAm.R", package = "dnaEPICO")

  if (rmd == "" || script == "")
    stop("DNAm.Rmd or DNAm.R not found in package.")

  outDirAbs <- makeAbs(outputDir)
  dir.create(outDirAbs, recursive = TRUE, showWarnings = FALSE)

  figDir <- file.path(outputDir, "figures")
  dir.create(figDir, recursive = TRUE, showWarnings = FALSE)


  # Build argument list for external Rscript
  arg_list <- c(
    "--rmd",        shQuote(rmd),
    "--output",     shQuote(output),
    "--outputDir",  shQuote(outDirAbs),
    "--qcDir",      shQuote(makeAbs(qcDir)),
    "--preDir",     shQuote(makeAbs(preprocessingDir)),
    "--postDir",    shQuote(makeAbs(postprocessingDir)),
    "--svaDir",     shQuote(makeAbs(svaDir)),
    "--glmDir",     shQuote(makeAbs(glmDir)),
    "--glmmDir",    shQuote(makeAbs(glmmDir)),
    "--title",      shQuote(reportTitle),
    "--author",     shQuote(author),
    "--date",       shQuote(date),
    "--figDir", shQuote(figDir)
  )

  # Construct Rscript command
  cmd <- paste("Rscript", shQuote(script), paste(arg_list, collapse = " "))

  # Run command
  message("Running dnamReport():")
  message(cmd)

  if (.Platform$OS.type == "windows") {
    shell(cmd)
  } else {
    system(cmd)
  }
}
