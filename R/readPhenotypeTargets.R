validateReadPhenotypeInputsDnaEpico <- function(phenoFile, SampleID, nSamples) {
    if (!file.exists(phenoFile)) {
        stop("phenoFile does not exist: ", phenoFile, call. = FALSE)
    }
    if (!is.character(SampleID) || length(SampleID) != 1L ||
        is.na(SampleID) || !nzchar(trimws(SampleID))) {
        stop("SampleID must be a single, non-empty character string.",
            call. = FALSE
        )
    }
    if (length(nSamples) != 1L || (!is.na(nSamples) &&
        (!is.numeric(nSamples) || !is.finite(nSamples) || nSamples < 1L ||
            nSamples != floor(nSamples)))) {
        stop("nSamples must be a single positive integer or NA.",
            call. = FALSE
        )
    }
    if (is.na(nSamples)) NA_integer_ else as.integer(nSamples)
}

logPhenotypeReadStartDnaEpico <- function(
    phenoFile, sepType, SampleID, nSamples, verbose, logPath
) {
    resolved_separator <- resolveSeparatorMinfiEwasWater(sepType)
    emitLogMinfiEwasWater(lines = c(
        "============================================================",
        paste("Phenotype file:           ", phenoFile),
        paste("Separator type:           ", if (is.null(resolved_separator)) {
            "default (',')"
        } else {
            sepType
        }),
        paste("SampleID column:          ", SampleID),
        paste("nSamples limit:           ", if (is.na(nSamples)) {
            "all"
        } else {
            nSamples
        })
    ), verbose = verbose, log_path = logPath)
}

readPhenotypeFileDnaEpico <- function(phenoFile, sepType) {
    separator <- resolveSeparatorMinfiEwasWater(sepType)
    if (is.null(separator)) {
        return(utils::read.csv(phenoFile, stringsAsFactors = FALSE))
    }
    utils::read.csv(phenoFile, sep = separator, stringsAsFactors = FALSE)
}

preparePhenotypeTargetsDnaEpico <- function(targets, SampleID, nSamples) {
    if (!(SampleID %in% colnames(targets))) {
        stop("SampleID column not found in phenotype data: ", SampleID,
            call. = FALSE
        )
    }
    if (nrow(targets) == 0L) {
        stop("The phenotype file contains no sample rows.", call. = FALSE)
    }
    targets[[SampleID]] <- validateSampleIdentifiersDnaEpico(
        targets[[SampleID]], paste0("Phenotype column '", SampleID, "'")
    )
    if (!is.na(nSamples) && nSamples < nrow(targets)) {
        return(list(
            targets = targets[seq_len(nSamples), , drop = FALSE],
            subsetLine = paste("Using the first", nSamples, "samples.")
        ))
    }
    list(
        targets = targets,
        subsetLine = paste("Using all", nrow(targets), "samples.")
    )
}

logPhenotypeReadResultDnaEpico <- function(prepared, verbose, logPath) {
    targets <- prepared$targets
    preview_cols <- seq_len(min(ncol(targets), 5L))
    preview_lines <- utils::capture.output(utils::head(
        targets[, preview_cols, drop = FALSE]
    ))
    emitLogMinfiEwasWater(lines = c(
        prepared$subsetLine,
        paste("Phenotype file loaded with", nrow(targets), "samples and",
            ncol(targets), "columns."
        ),
        "Preview of targets:", preview_lines,
        "============================================================"
    ), verbose = verbose, log_path = logPath)
}

#' Read phenotype targets for shared dnaEPICO workflows
#'
#' Read the phenotype table used by shared `dnaEPICO` workflows, validate the
#' sample identifier column, optionally subset the first `nSamples`, and return
#' the targets as a base `data.frame`.
#'
#' @param phenoFile Character. Path to the phenotype table on disk.
#' @param sepType Character or `NULL`. Field separator used in `phenoFile`. Use
#'   `NULL` for a standard comma-separated file, `'\\t'` for a tab-delimited
#' file, or another single-character separator accepted by `utils::read.csv()`.
#' @param nSamples Integer or `NA`. Number of rows to keep from the start of the
#'   phenotype table. The default `NA` reads and returns all rows.
#' @param SampleID Character. Name of the column containing sample identifiers
#'   that will later be used to name methylation-array samples.
#' @param verbose Logical. If `TRUE`, emit progress and preview messages with
#'   `message()`. The default is `FALSE`.
#' @param logs Logical. If `TRUE`, write the same progress messages to a log
#'   file. The default is `FALSE`.
#' @param log_dir Character or `NULL`. Directory where the log file should be
#'   written when `logs = TRUE`. If `NULL`, the current working directory is
#'   used.
#' @param log_file Character. File name used when `logs = TRUE`. The default is
#'   `'log_readPhenotypeTargets.txt'`.
#'
#' @return A `data.frame` containing the phenotype targets.
#'
#' @examples
#' tmp <- tempdir()
#' pheno <- data.frame(
#'     Sample_Name = c("S1", "S2"),
#'     Sex = c("F", "M"),
#'     stringsAsFactors = FALSE
#' )
#' pheno_file <- file.path(tmp, "pheno.csv")
#' utils::write.csv(pheno, pheno_file, row.names = FALSE)
#' targets <- readPhenotypeTargets(
#'     phenoFile = pheno_file,
#'     SampleID = "Sample_Name"
#' )
#' stopifnot(is.data.frame(targets))
#' stopifnot(nrow(targets) == 2L)
#'
#' @export
readPhenotypeTargets <- function(
    phenoFile, sepType = NULL, nSamples = NA,
    SampleID = "Sample_Name", verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_readPhenotypeTargets.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    nSamples <- validateReadPhenotypeInputsDnaEpico(
        phenoFile, SampleID, nSamples
    )
    logPhenotypeReadStartDnaEpico(
        phenoFile, sepType, SampleID, nSamples, verbose, log_path
    )
    prepared <- preparePhenotypeTargetsDnaEpico(
        readPhenotypeFileDnaEpico(phenoFile, sepType), SampleID, nSamples
    )
    logPhenotypeReadResultDnaEpico(prepared, verbose, log_path)
    prepared$targets
}
