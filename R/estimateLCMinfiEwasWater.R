prepareCellCompositionInputsDnaEpico <- function(
    beta, targets, SampleID, lcRef, constrained
) {
    constrained <- validateLogicalScalarDnaEpico(constrained, "constrained")
    if (!is.character(lcRef) || length(lcRef) != 1L || is.na(lcRef) ||
        !nzchar(trimws(lcRef))) {
        stop("lcRef must be one non-empty reference name.", call. = FALSE)
    }
    validateCellCompositionBetaStructureDnaEpico(
        beta = beta, requireSampleNames = TRUE, objectName = "beta"
    )
    if (!(SampleID %in% colnames(targets))) {
        stop("SampleID column not found in targets: ", SampleID,
            call. = FALSE
        )
    }
    beta_sample_ids <- as.character(colnames(beta))
    target_sample_ids <- validateSampleIdentifiersDnaEpico(
        targets[[SampleID]], paste0("Phenotype column '", SampleID, "'")
    )
    target_match <- matchSampleIdentifiersDnaEpico(
        query = beta_sample_ids, reference = target_sample_ids,
        queryLabel = "Beta-matrix sample identifiers",
        referenceLabel = paste0("phenotype column '", SampleID, "'"),
        requireSameSet = TRUE
    )
    list(
        constrained = constrained,
        internal = lcRef %in% c("saliva", "salivaEPIC"),
        sampleIDs = beta_sample_ids,
        range = summarizeMethylationRangeDnaEpico(beta, "beta"),
        targets = targets[target_match, , drop = FALSE]
    )
}

estimateCellCompositionDnaEpico <- function(
    beta, lcRef, constrained, internal
) {
    if (isTRUE(internal)) {
        return(list(
            values = estimateLCFromBetaDnaEpico(
                meth = beta, ref = lcRef, constrained = constrained
            ),
            method = "estimateLC"
        ))
    }
    list(
        values = ENmix::estimateCellProp(
            userdata = beta, refdata = lcRef,
            nonnegative = TRUE, normalize = FALSE, nProbes = 50,
            refplot = FALSE
        ),
        method = "ENmix::estimateCellProp"
    )
}

mergeCellCompositionDnaEpico <- function(targets, lc, phenoOrder) {
    rownames(targets) <- NULL
    pheno_lc <- cbind(targets, as.data.frame(lc, stringsAsFactors = FALSE))
    pheno_lc <- pheno_lc[, !duplicated(colnames(pheno_lc)), drop = FALSE]
    lead_cols <- splitOptionMinfiEwasWater(phenoOrder, sep = ";")
    lead_cols <- lead_cols[lead_cols %in% colnames(pheno_lc)]
    remaining_cols <- setdiff(colnames(pheno_lc), lead_cols)
    pheno_lc[, c(lead_cols, remaining_cols), drop = FALSE]
}

logCellCompositionStartDnaEpico <- function(inputs, lcRef, verbose, logPath) {
    emitLogMinfiEwasWater(c(
        paste("Cell composition reference:", lcRef),
        if (inputs$internal) {
            "Using internal Houseman implementation (estimateLC)."
        } else {
            "Using ENmix Houseman-based cell composition."
        },
        formatMethylationRangeLogDnaEpico(inputs$range)
    ), verbose = verbose, log_path = logPath)
}

logCellCompositionResultDnaEpico <- function(
    phenoLC, method, verbose, logPath
) {
    emitLogMinfiEwasWater(c(
        paste("Cell composition method:   ", method),
        paste("phenoLC columns:           ", ncol(phenoLC)),
        paste("phenoLC rows:              ", nrow(phenoLC)),
        "============================================================"
    ), verbose = verbose, log_path = logPath)
}

newCellCompositionResultDnaEpico <- function(
    estimate, phenoLC, inputs, SampleID, lcRef
) {
    structure(list(
        lc = estimate$values, phenoLC = phenoLC,
        sampleIDs = inputs$sampleIDs, SampleID = SampleID, ref = lcRef,
        method = estimate$method, methylationRange = inputs$range
    ), class = "dnaEPICO_minfiEwasWater_lc")
}

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
    phenoOrder = paste0(
        "Sample_Name;Timepoint;Sex;PredSex;Basename;",
        "Sentrix_ID;Sentrix_Position"
    ),
    constrained = FALSE, verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_estimateLCMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    inputs <- prepareCellCompositionInputsDnaEpico(
        beta, targets, SampleID, lcRef, constrained
    )

    logCellCompositionStartDnaEpico(inputs, lcRef, verbose, log_path)

    estimate <- estimateCellCompositionDnaEpico(
        beta, lcRef, inputs$constrained, inputs$internal
    )
    if (nrow(estimate$values) != ncol(beta)) {
        stop("Cell-composition output has ", nrow(estimate$values),
            " rows for ",
            ncol(beta), " beta-matrix samples.",
            call. = FALSE
        )
    }
    phenoLC <- mergeCellCompositionDnaEpico(
        inputs$targets, estimate$values, phenoOrder
    )

    logCellCompositionResultDnaEpico(
        phenoLC, estimate$method, verbose, log_path
    )
    newCellCompositionResultDnaEpico(
        estimate, phenoLC, inputs, SampleID, lcRef
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
