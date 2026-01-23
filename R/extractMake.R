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
#' @return
#' Invisibly returns \code{NULL}.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' ## Copy the dnaEPICO Makefile pipeline into a project directory
#' extractMake(
#'   destDir = getwd(),
#'   overwrite = FALSE
#' )
#' }
#'
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
