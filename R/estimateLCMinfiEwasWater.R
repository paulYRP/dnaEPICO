#' Estimate cell composition for preprocessingMinfiEwasWater
#'
#' Estimate cell proportions from beta values using `estimateLC()` for saliva
#' reference panels or `ENmix::estimateCellProp()` for other supported
#' references, then merge the estimates into the phenotype table.
#'
#' @param beta Numeric matrix of beta values with probes in rows and samples in
#'   columns.
#' @param targets Phenotype data frame aligned with the columns of `beta`.
#' @param SampleID Character. Phenotype column containing sample identifiers.
#' @param lcRef Character. Cell-composition reference. Internal saliva-based
#'   references supported through `estimateLC()` are `'saliva'` and
#'   `'salivaEPIC'`. Other references are passed to
#'   `ENmix::estimateCellProp()`.
#' @param phenoOrder Character vector or semicolon-separated phenotype columns
#'   to place first in the merged `phenoLC` output.
#' @param constrained Logical. Passed to `estimateLC()` when an internal saliva
#'   reference is used. If `TRUE`, estimated proportions are constrained to sum
#'   to one.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_minfiEwasWater_lc'` containing the cell
#'   proportion matrix, merged phenotype table, reference name, and method used.
#'
#' @examples
#' ref_file <- system.file("extdata", "saliva.txt", package = "dnaEPICO")
#' beta <- as.matrix(utils::read.table(ref_file))[1:20, , drop = FALSE]
#' colnames(beta) <- c("sample1", "sample2")
#' targets <- data.frame(
#'     Sample_Name = colnames(beta),
#'     Timepoint = c("T1", "T2"),
#'     stringsAsFactors = FALSE
#' )
#' lc_data <- estimateLCMinfiEwasWater(
#'     beta = beta,
#'     targets = targets,
#'     lcRef = "saliva",
#'     phenoOrder = "Sample_Name;Timepoint"
#' )
#' stopifnot(is.data.frame(lc_data$phenoLC))
#'
#' @export
estimateLCMinfiEwasWater <- function(
    beta, targets, SampleID = "Sample_Name",
    lcRef = "salivaEPIC",
        phenoOrder = "Sample_Name;Timepoint;Sex;PredSex;Basename;Sentrix_ID;Sentrix_Position",
    constrained = FALSE, verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_estimateLCMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    constrained <- validateLogicalScalarDnaEpico(
        constrained,
        "constrained"
    )
    if (!is.character(lcRef) || length(lcRef) != 1L || is.na(lcRef) ||
        !nzchar(trimws(lcRef))) {
        stop("lcRef must be one non-empty reference name.", call. = FALSE)
    }
    ewas_refs <- c("saliva", "salivaEPIC")
    use_internal <- lcRef %in% ewas_refs

    validateCellCompositionBetaStructureDnaEpico(
        beta = beta,
        requireSampleNames = TRUE, objectName = "beta"
    )
    beta_sample_ids <- as.character(colnames(beta))
    beta_range <- summarizeMethylationRangeDnaEpico(beta, "beta")
    if (!(SampleID %in% colnames(targets))) {
        stop("SampleID column not found in targets: ", SampleID,
            call. = FALSE
        )
    }
    target_sample_ids <- validateSampleIdentifiersDnaEpico(
        targets[[SampleID]],
        paste0("Phenotype column '", SampleID, "'")
    )
    target_match <- matchSampleIdentifiersDnaEpico(
        query = beta_sample_ids,
        reference = target_sample_ids,
            queryLabel = "Beta-matrix sample identifiers",
        referenceLabel = paste0(
            "phenotype column '", SampleID,
            "'"
        ), requireSameSet = TRUE
    )
    targets <- targets[target_match, , drop = FALSE]
    rownames(targets) <- NULL

    emitLogMinfiEwasWater(c(
        paste(
            "Cell composition reference:",
            lcRef
        ), if (use_internal) {
            "Using internal Houseman implementation (estimateLC)."
        } else {
            "Using ENmix Houseman-based cell composition."
        }, formatMethylationRangeLogDnaEpico(beta_range)
    ), verbose = verbose, log_path = log_path)

    if (use_internal) {
        lc <- estimateLCFromBetaDnaEpico(
            meth = beta,
            ref = lcRef, constrained = constrained
        )
        method <- "estimateLC"
    } else {
        lc <- ENmix::estimateCellProp(
            userdata = beta, refdata = lcRef,
            nonnegative = TRUE, normalize = FALSE, nProbes = 50,
            refplot = FALSE
        )
        method <- "ENmix::estimateCellProp"
    }

    if (nrow(lc) != ncol(beta)) {
        stop("Cell-composition output has ", nrow(lc), " rows for ",
            ncol(beta), " beta-matrix samples.",
            call. = FALSE
        )
    }

    phenoLC <- cbind(targets, as.data.frame(lc, stringsAsFactors = FALSE))
    phenoLC <- phenoLC[, !duplicated(colnames(phenoLC)), drop = FALSE]

    lead_cols <- splitOptionMinfiEwasWater(phenoOrder, sep = ";")
    lead_cols <- lead_cols[lead_cols %in% colnames(phenoLC)]
    remaining_cols <- setdiff(colnames(phenoLC), lead_cols)
    phenoLC <- phenoLC[, c(lead_cols, remaining_cols), drop = FALSE]

    emitLogMinfiEwasWater(
        c(
            paste(
                "Cell composition method:   ",
                method
            ), paste("phenoLC columns:           ", ncol(phenoLC)),
            paste("phenoLC rows:              ", nrow(phenoLC)),
            "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(
        list(
            lc = lc, phenoLC = phenoLC, sampleIDs = beta_sample_ids,
            SampleID = SampleID, ref = lcRef, method = method,
            methylationRange = beta_range
        ),
        class = "dnaEPICO_minfiEwasWater_lc"
    )
}

#' Write the merged phenotype plus cell-composition table
#'
#' @param lcData Object returned by `estimateLCMinfiEwasWater()`.
#' @param file Character. Path to the CSV file to write.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return Invisibly returns `file`.
#'
#' @examples
#' ref_file <- system.file("extdata", "saliva.txt", package = "dnaEPICO")
#' beta <- as.matrix(utils::read.table(ref_file))[1:20, , drop = FALSE]
#' colnames(beta) <- c("sample1", "sample2")
#' targets <- data.frame(
#'     Sample_Name = colnames(beta),
#'     Timepoint = c("T1", "T2"),
#'     stringsAsFactors = FALSE
#' )
#' lc_data <- estimateLCMinfiEwasWater(
#'     beta = beta,
#'     targets = targets,
#'     lcRef = "saliva",
#'     phenoOrder = "Sample_Name;Timepoint"
#' )
#' output_file <- file.path(tempdir(), "phenoLC.csv")
#' writePhenoLCMinfiEwasWater(lcData = lc_data, file = output_file)
#' file.exists(output_file)
#'
#' @export
writePhenoLCMinfiEwasWater <- function(
    lcData, file, verbose = FALSE,
    logs = FALSE, log_dir = NULL,
        log_file = "log_writePhenoLCMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(lcData$phenoLC, file = file, row.names = FALSE)

    emitLogMinfiEwasWater(
        c(paste(
            "Saved phenoLC:             ",
            file
        ), "============================================================"),
        verbose = verbose, log_path = log_path
    )

    invisible(file)
}
