#' Export the dnaEPICO Makefile to a user directory
#'
#' Export the example Makefile supplied with `dnaEPICO` to a user directory.
#'
#' @details
#' The exported workflow requires GNU Make 4.3 or later. When the Makefile is
#' created, `extractMake()` records the absolute path to the `Rscript`
#' executable from the current R installation in the `RSCRIPT` Make variable.
#' This value can be overridden when running Make, for example
#' `make f4 RSCRIPT=/path/to/Rscript`, when another R installation or an HPC
#' module should be used.
#'
#' Quarto is required only for targets that render reports. The rules use
#' relative project paths and support project directories containing spaces.
#' GNU Make and Quarto are not required to install or load `dnaEPICO`.
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

    rscriptPath <- file.path(
        R.home("bin"),
        if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    )
    rscriptPath <- normalizePath(rscriptPath, winslash = "/", mustWork = TRUE)
    rscriptPath <- gsub("$", "$$", rscriptPath, fixed = TRUE)
    rscriptPath <- gsub("#", "\\#", rscriptPath, fixed = TRUE)

    makefile <- readLines(makefileDest, warn = FALSE)
    placeholderLine <- which(makefile == "RSCRIPT ?= @DNAEPICO_RSCRIPT@")
    if (length(placeholderLine) != 1L) {
        stop("The packaged Makefile has an invalid RSCRIPT placeholder.",
            call. = FALSE)
    }
    makefile[[placeholderLine]] <- paste("RSCRIPT ?=", rscriptPath)
    writeLines(makefile, makefileDest, useBytes = TRUE)

    normalizePath(makefileDest, winslash = "/", mustWork = FALSE)
}
