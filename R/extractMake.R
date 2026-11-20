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
#' `MODEL` selects one analysis configuration: its input files, settings and
#' output directory. `MODELS` lists configurations for a grouped run.
#' Each analysis can fit many CpG-wise statistical models.
#'
#' `N_CORES` limits the total CPU cores used by one Make workflow. The rules
#' share the available cores among the analyses in `MODELS`. An R worker is a
#' R process that fits a subset of CpGs. With 8 available cores and
#' 3 analyses, the maximum R workers per analysis are 3, 3 and 2. Analyses
#' wait when all available workers are in use. Steps within each analysis
#' run in order.
#'
#' The shared rules set the OpenMP, OpenBLAS, MKL, Accelerate and NumExpr
#' thread limits to one before starting R. This limits the threads used by
#' numerical libraries within each R process. `make status` displays the
#' available cores and maximum R workers per analysis. Fitting can use fewer
#' workers for small CpG sets or limited memory. Memory requirements apply
#' separately to each analysis.
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
