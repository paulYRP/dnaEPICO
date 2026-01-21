#' Copy dnaEPICO Makefile to a user directory
#'
#' Copies the example Makefile pipeline shipped with dnaEPICO to a
#' user-specified directory for local execution or modification.
#'
#' @param destDir Character. Destination directory where the Makefile
#'   will be copied.
#' @param overwrite Logical. Whether to overwrite an existing file.
#'   Default is FALSE.
#'
#' @return Invisibly returns the path to the copied Makefile.
#' @export
extractMake <- function(destDir, overwrite = FALSE) {

  if (!dir.exists(destDir)) {
    stop("Destination directory does not exist: ", destDir)
  }

  makefileSrc <- system.file(
    "extdata",
    "make",
    "Makefile.model.pipeline",
    package = "dnaEPICO",
    mustWork = TRUE
  )

  makefileDest <- file.path(destDir, "Makefile")

  file.copy(
    makefileSrc,
    makefileDest,
    overwrite
  )

  invisible(makefileDest)
}
