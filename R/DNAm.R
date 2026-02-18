#' Generate a DNA methylation PDF report
#' @import rmarkdown
#' @import tiff
#'
#' @param output Character. Name of the output PDF file.
#' @param outputDir Character. Directory where the report will be saved.
#' @param qcDir Character. Directory containing ENmix QC figures.
#' @param preprocessingDir Character. Directory containing preprocessing QC figures.
#' @param postprocessingDir Character. Directory containing postprocessing metric figures.
#' @param svaDir Character. Directory containing SVA figures.
#' @param glmDir Character. Directory containing GLM figures.
#' @param glmmDir Character. Directory containing GLMM figures.
#' @param figDir Character. Directory where figures will be copied for the report.
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
#'   figDir = "reports/figures",
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
    figDir = "reports/figures",
    reportTitle = "DNA methylation",
    author = "School of Biomedical Sciences",
    date = format(Sys.Date(), "%B %d, %Y")
) {

# Fix Windows path handling for Pandoc/LaTeX:
normalize_path <- function(x) gsub("\\\\", "/", normalizePath(x, winslash = "/", mustWork = FALSE))

rmd <- system.file("extdata", "DNAm.Rmd", package = "dnaEPICO")

if (rmd == "") {
    stop("DNAm.Rmd not found inside package (inst/extdata/).")
  }

if (!requireNamespace("tinytex", quietly = TRUE) || !tinytex::is_tinytex()) {
    message("TinyTeX/LaTeX not available, skipping PDF report generation.")
    return(invisible(NULL))
  }

rmd        <- normalize_path(rmd)
outputDir  <- normalize_path(outputDir)
qcDir      <- normalize_path(qcDir)
preDir     <- normalize_path(preprocessingDir)
postDir    <- normalize_path(postprocessingDir)
svaDir     <- normalize_path(svaDir)
glmDir     <- normalize_path(glmDir)
glmmDir    <- normalize_path(glmmDir)
figDir    <- normalize_path(figDir)

if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
if (!dir.exists(figDir)) dir.create(figDir, recursive = TRUE)

rmarkdown::render(
  input = rmd,
  output_file = output,
  output_dir = outputDir,
  params = list(
    reportTitle = reportTitle,
    author = author,
    date = date,
    qcDir = qcDir,
    preprocessingDir = preDir,
    postprocessingDir = postDir,
    svaDir = svaDir,
    glmDir = glmDir,
    glmmDir = glmmDir,
    figDir = figDir
  ),
  knit_root_dir = dirname(rmd)
)
}