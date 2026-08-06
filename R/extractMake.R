#' Copy the dnaEPICO Makefile to a user directory
#'
#' Copy the example Makefile supplied with `dnaEPICO` to a user directory.
#'
#' @details
#' The exported workflow requires GNU Make 4.3 or later and `Rscript` on the
#' system path. Quarto is also required for report targets. The rules use
#' relative project paths and support project directories containing spaces.
#'
#' @param destDir Character. Destination directory for the Makefile.
#' @param overwrite Logical. Whether to overwrite an existing `Makefile` in
#'   `destDir`. The default is `FALSE`.
#'
#' @return Character scalar containing the path to the copied `Makefile`.
#'
#' @examples
#' tmp <- file.path(tempdir(), "dnaEPICO-make-example")
#' dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
#' makefile_path <- extractMake(
#'     destDir = tmp,
#'     overwrite = TRUE
#' )
#' stopifnot(file.exists(makefile_path))
#'
#' @export
extractMake <- function(destDir, overwrite = FALSE) {
    if (!is.character(destDir) || length(destDir) != 1L || !nzchar(destDir)) {
        stop("destDir must be a single, non-empty character path.",
            call. = FALSE
        )
    }

    if (!is.logical(overwrite) || length(overwrite) != 1L ||
        is.na(overwrite)) {
        stop("overwrite must be a single TRUE or FALSE value.",
            call. = FALSE
        )
    }

    if (!dir.exists(destDir)) {
        stop("Destination directory does not exist: ", destDir,
            call. = FALSE
        )
    }

    makefileSrc <- system.file("extdata", "make", "Makefile.model.pipeline",
        package = "dnaEPICO", mustWork = TRUE
    )

    makefileDest <- file.path(destDir, "Makefile")

    copied <- file.copy(makefileSrc, makefileDest, overwrite = overwrite)

    if (!isTRUE(copied)) {
        stop("Failed to copy the packaged Makefile to: ", makefileDest,
            call. = FALSE
        )
    }

    normalizePath(makefileDest, winslash = "/", mustWork = FALSE)
}
