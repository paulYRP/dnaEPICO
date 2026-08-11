#' Estimate surrogate variables from ENmix control probes
#'
#' @param RGSet An `RGChannelSet`.
#' @param ctrlSvaPercVar Numeric. Proportion of variance explained by control
#'   probes, passed to `ENmix::ctrlsva()`.
#' @param ctrlSvaFlag Integer. Control-probe flag passed to
#'   `ENmix::ctrlsva()`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_svaEnmix_sva'` containing the surrogate
#'   variable matrix and the parameters used to estimate it.
#'
#' @description
#' Run `ENmix::ctrlsva()` on an `RGChannelSet` and return the surrogate variable
#' matrix as an in-memory object.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' sva_data <- estimateSvaEnmixControls(
#'   RGSet = ex$RGSet,
#'   ctrlSvaPercVar = 0.5,
#'   ctrlSvaFlag = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' sva_data$K
#'
#' @export
estimateSvaEnmixControls <- function(RGSet, ctrlSvaPercVar = 0.9,
    ctrlSvaFlag = 1, verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_estimateSvaEnmixControls.txt") {
    log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir,
        log_file = log_file)
    ctrlSvaPercVar <- validateProbabilityDnaEpico(ctrlSvaPercVar,
        "ctrlSvaPercVar")
    if (length(ctrlSvaFlag) != 1L || !is.numeric(ctrlSvaFlag) ||
        is.na(ctrlSvaFlag) || !is.finite(ctrlSvaFlag) || !(ctrlSvaFlag %in%
        c(1, 2))) {
        stop("ctrlSvaFlag must be 1 (percvar) or 2 (one requested component).",
            call. = FALSE) }
    sva <- ENmix::ctrlsva(rgSet = RGSet, percvar = ctrlSvaPercVar,
        flag = ctrlSvaFlag)
    sva <- as.matrix(sva)
    nan_sva_converted <- sum(is.nan(sva))
    if (nan_sva_converted > 0L) {
        sva[is.nan(sva)] <- NA_real_ }
    sample_names <- validateSampleIdentifiersDnaEpico(
        svaEnmixGetRGSetSampleNames(RGSet),
        "RGSet sample identifiers")
    if (!is.numeric(sva) || nrow(sva) != length(sample_names) ||
        ncol(sva) == 0L || anyNA(sva) || any(!is.finite(sva)) ||
        is.null(rownames(sva))) {
        if (nan_sva_converted > 0L) {
            ctrlsva_requirement <- paste(
            "ENmix::ctrlsva() must return a finite numeric matrix",
                "with one named row per RGSet sample")
            stop(sprintf("%s; %s upstream NaN value(s) were converted to NA.",
                ctrlsva_requirement, nan_sva_converted), call. = FALSE)
        }
        stop("ENmix::ctrlsva() must return a finite numeric matrix with one ",
            "named row per RGSet sample.", call. = FALSE)
    }
    sample_match <- matchSampleIdentifiersDnaEpico(query = sample_names,
        reference = rownames(sva), queryLabel = "RGSet sample identifiers",
        referenceLabel = "ENmix SVA row names", requireSameSet = TRUE)
    sva <- sva[sample_match, , drop = FALSE]
    colnames(sva) <- paste0("PC", seq_len(ncol(sva)))
    emitLogMinfiEwasWater(c(paste("ctrlSva percvar:          ",
        ctrlSvaPercVar), paste("ctrlSva flag:             ",
        ctrlSvaFlag), paste("ctrlSva NaN to NA:        ", nan_sva_converted),
        paste("Number of surrogate variables:", ncol(sva)),
        "Surrogate variables matrix (first few rows):",
            previewLinesMinfiEwasWater(sva),
        "============================================================"),
        verbose = verbose, log_path = log_path)
    structure(list(sva = sva, K = ncol(sva), ctrlSvaPercVar = ctrlSvaPercVar,
        ctrlSvaFlag = ctrlSvaFlag, nanConverted = nan_sva_converted),
        class = "dnaEPICO_svaEnmix_sva") }

#' Merge surrogate variables into the phenotype table
#'
#' @param targets Phenotype data frame aligned with the samples in `sva`.
#' @param sva Numeric matrix of surrogate variables with samples in rows.
#' @param SampleID Character. Name of the phenotype sample identifier column.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A phenotype data frame with the surrogate variables appended.
#'
#' @description
#' Merge the surrogate variable matrix back into the phenotype table while
#' preserving the original row order of `targets`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleSvaAnalysisStateDnaEpico()
#' merged_pheno <- mergeSvaTargetsEnmix(
#'   targets = ex$targets,
#'   sva = ex$sva,
#'   SampleID = "Sample_Name",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' colnames(merged_pheno)[seq_len(4)]
#'
#' @export
mergeSvaTargetsEnmix <- function(targets, sva, SampleID = "Sample_Name",
    verbose = FALSE, logs = FALSE, log_dir = NULL, log_file =
        "log_mergeSvaTargetsEnmix.txt") {
    log_path <- resolveLogPathMinfiEwasWater(logs = logs,
        log_dir = log_dir, log_file = log_file)
    if (!(SampleID %in% colnames(targets))) {
        stop("SampleID column not found in targets: ",
            SampleID, call. = FALSE)
    }
    if (!is.matrix(sva) || !is.numeric(sva) || ncol(sva) ==
        0L) {
        stop("sva must be a numeric matrix with at least one column.",
            call. = FALSE) }
    if (is.null(rownames(sva))) {
        stop(
            "sva must have row names that match the phenotype SampleID column.",
            call. = FALSE) }
    if (anyNA(sva) || any(!is.finite(sva))) {
        stop("sva contains missing or non-finite values.",
            call. = FALSE) }
    sva_data <- data.frame(sample_id = rownames(sva), sva,
        row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(sva_data)[1] <- SampleID
    surrogate_columns <- setdiff(colnames(sva_data), SampleID)
    duplicate_columns <- intersect(colnames(targets), surrogate_columns)
    if (length(duplicate_columns) > 0L) {
        duplicate_columns_text <- paste(duplicate_columns,
            collapse = ", ")
        duplicate_message <- paste(
            "The phenotype table already contains output columns that would",
            "be duplicated")
        duplicate_advice <- paste(
            "Start from the pre-SVA phenotype file or remove the existing",
            "columns explicitly.")
        stop(sprintf("%s: %s. %s", duplicate_message, duplicate_columns_text,
            duplicate_advice), call. = FALSE)
    }
    match_idx <- matchSampleIdentifiersDnaEpico(query = targets[[SampleID]],
        reference = sva_data[[SampleID]], queryLabel = paste0(
            "Phenotype column '",
            SampleID, "'"), referenceLabel = "SVA row names",
        requireSameSet = TRUE)
    merged_pheno <- cbind(targets, sva_data[match_idx,
        surrogate_columns, drop = FALSE])
    emitLogMinfiEwasWater(c(paste("Merged phenotype rows:     ",
        nrow(merged_pheno)), paste("Merged phenotype columns:  ",
        ncol(merged_pheno)),
            "============================================================"),
        verbose = verbose, log_path = log_path)
    merged_pheno }

alignSvaAnalysisInputDnaEpico <- function(sva, RGSet) {
    sample_names <- svaEnmixGetRGSetSampleNames(RGSet)
    if (!length(sample_names)) {
    stop("Could not determine sample names from the loaded RGSet.",
        call. = FALSE
    )
    }
    if (!is.matrix(sva) || !is.numeric(sva) || !ncol(sva)) {
    stop("sva must be a numeric matrix with at least one surrogate ",
        "variable.",
        call. = FALSE
    )
    }
    if (anyNA(sva) || any(!is.finite(sva))) {
    stop("sva contains missing or non-finite values.", call. = FALSE)
    }
    if (is.null(rownames(sva))) {
    stop("sva must have sample identifiers in its row names.",
        call. = FALSE
    )
    }
    matched <- matchSampleIdentifiersDnaEpico(
    query = sample_names, reference = rownames(sva),
    queryLabel = "RGSet sample identifiers", referenceLabel = "SVA row names",
    requireSameSet = TRUE
    )
    list(sva = sva[matched, , drop = FALSE], sampleNames = sample_names)
}

technicalFactorsSvaEnmixDnaEpico <- function(
    RGSet, SentrixIDColumn, SentrixPositionColumn
) {
    col_data <- svaEnmixGetRGSetColData(RGSet)
    required <- c(SentrixIDColumn, SentrixPositionColumn)
    missing <- setdiff(required, colnames(col_data))
    if (length(missing)) {
    missing_text <- paste(missing, collapse = ", ")
    stop(sprintf(
        "RGSet column data is missing: %s",
        missing_text
    ), call. = FALSE)
    }
    sentrix_id <- as.factor(col_data[[SentrixIDColumn]])
    sentrix_position <- as.factor(col_data[[SentrixPositionColumn]])
    if (anyNA(sentrix_id) || anyNA(sentrix_position)) {
    stop("Sentrix chip and position columns cannot contain missing values.",
        call. = FALSE
    )
    }
    data <- data.frame(
    SentrixID = sentrix_id, SentrixPosition = sentrix_position
    )
    terms <- colnames(data)[vapply(
    data, function(values) nlevels(droplevels(values)) > 1L, logical(1)
    )]
    formula <- if (length(terms)) {
    stats::reformulate(termlabels = terms, response = "sva_value")
    } else {
    stats::as.formula("sva_value ~ 1")
    }
    list(
    data = data, terms = terms, formula = formula,
    sentrixID = sentrix_id, sentrixPosition = sentrix_position
    )
}

fullModelsSvaEnmixDnaEpico <- function(sva, technical) {
    lapply(seq_len(ncol(sva)), function(index) {
    model_data <- technical$data
    model_data$sva_value <- sva[, index]
    stats::lm(technical$formula, data = model_data)
    })
}

reduceOneSvaModelDnaEpico <- function(model, captureWarnings) {
    steps <- list()
    repeat {
    drop_table <- captureWarnings(stats::drop1(model, test = "F"))
    probabilities <- drop_table$`Pr(F)`
    valid <- which(!is.na(probabilities))
    if (!length(valid)) break
    max_index <- valid[which.max(probabilities[valid])]
    max_probability <- probabilities[[max_index]]
    term <- rownames(drop_table)[[max_index]]
    if (!is.finite(max_probability) || max_probability <= 0.05 ||
        identical(term, "<none>")) {
        break
    }
    model <- stats::update(model, paste(". ~ . -", term))
    steps[[length(steps) + 1L]] <- list(
        dropterm = drop_table, summary = summary(model)
    )
    }
    list(model = model, steps = steps)
}

reduceSvaModelsDnaEpico <- function(fullModels) {
    warning_state <- new.env(parent = emptyenv())
    warning_state$messages <- character(0)
    capture_warnings <- function(expression) {
    withCallingHandlers(expression, warning = function(condition) {
        warning_state$messages <- unique(c(
        warning_state$messages, conditionMessage(condition)
        ))
        invokeRestart("muffleWarning")
    })
    }
    reduced <- lapply(fullModels, reduceOneSvaModelDnaEpico,
    captureWarnings = capture_warnings
    )
    reduced_models <- lapply(reduced, `[[`, "model")
    list(
    models = reduced_models, steps = lapply(reduced, `[[`, "steps"),
    anovaFull = lapply(
        fullModels,
        function(model) capture_warnings(stats::anova(model))
    ),
    anovaReduced = lapply(
        reduced_models,
        function(model) capture_warnings(stats::anova(model))
    ),
    warnings = warning_state$messages
    )
}

svaAnalysisLogLinesDnaEpico <- function(aligned, technical, reduced) {
    terms <- if (length(technical$terms)) {
    paste(technical$terms, collapse = ", ")
    } else {
    "intercept only"
    }
    c(
    paste("Number of surrogate variables (K):", ncol(aligned$sva)),
    paste("SentrixID class:            ", class(technical$sentrixID)[1]),
    paste(
        "SentrixID unique levels:    ",
        length(unique(technical$sentrixID))
    ),
    previewLinesMinfiEwasWater(table(technical$sentrixID)),
    paste(
        "SentrixPosition class:      ",
        class(technical$sentrixPosition)[1]
    ),
    paste(
        "SentrixPosition unique levels:",
        length(unique(technical$sentrixPosition))
    ),
    previewLinesMinfiEwasWater(table(technical$sentrixPosition)),
    paste("Technical terms modelled:     ", terms),
    "First row of SVA matrix:",
    previewLinesMinfiEwasWater(aligned$sva[1, , drop = FALSE]),
    paste("Sample names in SVA matrix:", paste(
        rownames(aligned$sva)[seq_len(min(5L, nrow(aligned$sva)))],
        collapse = ", "
    )),
    paste("Sample names in colData(RGSet):", paste(
        aligned$sampleNames[seq_len(min(5L, length(aligned$sampleNames)))],
        collapse = ", "
    )),
    paste("Model-selection warnings:    ", length(reduced$warnings)),
    if (length(reduced$warnings)) reduced$warnings else character(0),
    "============================================================"
    )
}

#' Analyze surrogate variables against Sentrix chip and position factors
#'
#' @param sva Numeric matrix of surrogate variables with samples in rows.
#' @param RGSet An `RGChannelSet` aligned with `sva`.
#' @param SentrixIDColumn Character. Name of the chip identifier column in
#'   `SummarizedExperiment::colData(RGSet)`.
#' @param SentrixPositionColumn Character. Name of the chip position column in
#'   `SummarizedExperiment::colData(RGSet)`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_svaEnmix_analysis'` containing the
#'   aligned Sentrix factors, full and reduced linear models, and ANOVA tables.
#'
#' @description
#' Fit linear models for each surrogate variable against informative Sentrix
#' chip and position factors, perform backward elimination with
#' `stats::drop1()`, and return the in-memory analysis objects. A constant
#' technical factor is omitted; if both are constant, an intercept-only model
#' is retained.
#'
#' @examples
#' ex <- dnaEPICO:::exampleSvaAnalysisStateDnaEpico()
#' analysis_data <- analyzeSvaEnmix(
#'   sva = ex$sva,
#'   RGSet = ex$RGSet,
#'   SentrixIDColumn = "Sentrix_ID",
#'   SentrixPositionColumn = "Sentrix_Position",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' analysis_data$K
#'
#' @export
analyzeSvaEnmix <- function(
    sva, RGSet, SentrixIDColumn = "Sentrix_ID",
    SentrixPositionColumn = "Sentrix_Position", verbose = FALSE,
    logs = FALSE, log_dir = NULL, log_file = "log_analyzeSvaEnmix.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    aligned <- alignSvaAnalysisInputDnaEpico(sva, RGSet)
    technical <- technicalFactorsSvaEnmixDnaEpico(
    RGSet, SentrixIDColumn, SentrixPositionColumn
    )
    full_models <- fullModelsSvaEnmixDnaEpico(aligned$sva, technical)
    reduced <- reduceSvaModelsDnaEpico(full_models)
    emitLogMinfiEwasWater(
    svaAnalysisLogLinesDnaEpico(aligned, technical, reduced),
    verbose = verbose, log_path = log_path
    )
    structure(list(
    sva = aligned$sva, K = ncol(aligned$sva),
    sentrixID = technical$sentrixID,
    sentrixPosition = technical$sentrixPosition,
    technicalTerms = technical$terms, fullModels = full_models,
    reducedModels = reduced$models, droptermSteps = reduced$steps,
    anovaFull = reduced$anovaFull, anovaReduced = reduced$anovaReduced,
    warnings = reduced$warnings
    ), class = "dnaEPICO_svaEnmix_analysis")
}

#' @keywords internal
#' @noRd
svaEnmixPositiveNumber <- function(value, fallback) {
    if (!is.numeric(value) || length(value) == 0L) {
    return(fallback)
    }

    value <- value[1]
    if (!is.finite(value) || value <= 0) {
    return(fallback)
    }

    value
}

#' @keywords internal
#' @noRd
svaEnmixSplitIndices <- function(K, page_limit = 6L) {
    K <- as.integer(svaEnmixPositiveNumber(K, 1))
    page_limit <- as.integer(svaEnmixPositiveNumber(
    page_limit,
    6
    ))
    indices <- seq_len(K)
    split(indices, ceiling(indices / page_limit))
}

#' @keywords internal
#' @noRd
svaEnmixMatrixPages <- function(K, page_limit = 6L) {
    row_blocks <- svaEnmixSplitIndices(K, page_limit)
    col_blocks <- svaEnmixSplitIndices(K, page_limit)
    pages <- vector("list", length(row_blocks) * length(col_blocks))
    page <- 1L

    for (row_block in row_blocks) {
    for (col_block in col_blocks) {
        pages[[page]] <- list(
        page = page, rowIndices = row_block,
        colIndices = col_block
        )
        page <- page + 1L
    }
    }

    pages
}

#' @keywords internal
#' @noRd
svaEnmixMatrixPlotDimensions <- function(
    width, height, res,
    panel_count
) {
    fallback_width <- 2000
    fallback_height <- 1000
    width <- svaEnmixPositiveNumber(width, fallback_width)
    height <- svaEnmixPositiveNumber(height, fallback_height)
    res <- svaEnmixPositiveNumber(res, NA_real_)
    panel_count <- svaEnmixPositiveNumber(panel_count, NA_real_)

    if (!is.finite(res) || res <= 0 || !is.finite(panel_count) ||
    panel_count <= 0L) {
    return(c(
        width = as.integer(ceiling(width)),
        height = as.integer(ceiling(height))
    ))
    }

    minimum_pixels <- ceiling(round(((panel_count * 1.6) + 1) *
    res, digits = 6L))
    c(width = as.integer(max(width, minimum_pixels)), height = as.integer(max(
    height,
    minimum_pixels
    )))
}

#' @keywords internal
#' @noRd
svaEnmixMatrixScale <- function(K, multiplier, minimum, maximum) {
    K <- svaEnmixPositiveNumber(K, 1)

    max(minimum, min(maximum, multiplier / K))
}

#' @keywords internal
#' @noRd
svaEnmixMatrixPanelMargin <- function(K) {
    if (K <= 5L) {
    return(c(2.4, 2.4, 0.5, 0.2))
    }
    if (K <= 8L) {
    return(c(1.9, 1.9, 0.3, 0.1))
    }

    c(1.5, 1.5, 0.2, 0.1)
}

#' @keywords internal
#' @noRd
svaEnmixMatrixPlotState <- function(
    analysisData, width, height,
    res
) {
    sentrix_id_levels <- levels(analysisData$sentrixID)
    sentrix_id_legend_limit <- 20L
    page_panel_limit <- 6L
    pages <- svaEnmixMatrixPages(analysisData$K, page_panel_limit)
    page_panel_count <- min(analysisData$K, page_panel_limit)
    dimensions <- svaEnmixMatrixPlotDimensions(
    width = width,
    height = height, res = res, panel_count = page_panel_count
    )

    list(
    width = dimensions[["width"]], height = dimensions[["height"]],
    pages = pages, pageCount = length(pages),
    pagePanelLimit = page_panel_limit,
    sentrixIDLevels = sentrix_id_levels,
    sentrixIDLegendLimit = sentrix_id_legend_limit,
    showSentrixIDLegend = length(sentrix_id_levels) <=
        sentrix_id_legend_limit
    )
}

#' @keywords internal
#' @noRd
svaEnmixLogMatrixPlotState <- function(state, res, verbose, log_path) {
    log_lines <- c(paste(
    "SVA matrix plot dimensions:", state$width,
    "x", state$height, "@", res
    ))

    if (state$pageCount > 1L) {
    log_lines <- c(log_lines, paste(
        "SVA matrix plot pages:",
        state$pageCount, "pages with up to", state$pagePanelLimit,
        "x", state$pagePanelLimit, "panels per page"
    ))
    }

    if (!isTRUE(state$showSentrixIDLegend)) {
    log_lines <- c(log_lines, paste(
        "SentrixID legend suppressed for matrix plot:",
        length(state$sentrixIDLevels), "levels exceed limit of",
        state$sentrixIDLegendLimit
    ))
    }

    emitLogMinfiEwasWater(log_lines, verbose = verbose, log_path = log_path)
}

#' @keywords internal
#' @noRd
svaEnmixMatrixPageFile <- function(file, page, page_count) {
    if (is.null(file) || page_count <= 1L) {
    return(file)
    }

    extension <- tools::file_ext(file)
    stem <- tools::file_path_sans_ext(file)
    if (!nzchar(extension)) {
    return(sprintf("%s_page%02d", stem, page))
    }

    sprintf("%s_page%02d.%s", stem, page, extension)
}

#' @keywords internal
#' @noRd
svaEnmixPlotSentrixID <- function(analysisData) {
    color_map <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(8, "Dark2")
    )(length(levels(analysisData$sentrixID)))
    sva <- analysisData$sva
    plot_data <- data.frame(
    SV1 = sva[, 1L], SV2 = sva[, 2L],
    group = analysisData$sentrixID
    )
    style <- adaptivePointStyleDnaEpico(nrow(plot_data))
    plot_object <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = SV1, y = SV2, colour = group
    )) +
    ggplot2::geom_point(size = style$size, alpha = max(style$alpha, 0.6)) +
    ggplot2::scale_colour_manual(values = color_map) +
    ggplot2::labs(
        title = NULL, x = "Surrogate Variable 1",
        y = "Surrogate Variable 2", colour = "Sentrix ID"
    ) +
    dnaEpicoModelPlotTheme()
    if (length(levels(analysisData$sentrixID)) > 20L) {
    plot_object <- plot_object + ggplot2::theme(legend.position = "none")
    }
    drawPlotObjectMinfiEwasWater(plot_object)
}

#' @keywords internal
#' @noRd
svaEnmixPlotSentrixPosition <- function(analysisData) {
    color_map <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(8, "Dark2")
    )(length(levels(analysisData$sentrixPosition)))
    sva <- analysisData$sva
    plot_data <- data.frame(
    SV1 = sva[, 1L], SV2 = sva[, 2L],
    group = analysisData$sentrixPosition
    )
    style <- adaptivePointStyleDnaEpico(nrow(plot_data))
    plot_object <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = SV1, y = SV2, colour = group
    )) +
    ggplot2::geom_point(size = style$size, alpha = max(style$alpha, 0.6)) +
    ggplot2::scale_colour_manual(values = color_map) +
    ggplot2::labs(
        title = NULL, x = "Surrogate Variable 1",
        y = "Surrogate Variable 2", colour = "Sentrix position"
    ) +
    dnaEpicoModelPlotTheme()
    drawPlotObjectMinfiEwasWater(plot_object)
}

#' @keywords internal
#' @noRd
svaEnmixPlotAssociations <- function(analysisData) {
    rows <- list()
    for (index in seq_along(analysisData$anovaFull)) {
    table <- as.data.frame(analysisData$anovaFull[[index]])
    p_column <- names(table)[grepl("Pr\\(", names(table))][1L]
    if (is.na(p_column) || !nzchar(p_column)) {
        next
    }
    terms <- rownames(table)
    keep <- !grepl("Residual", terms, ignore.case = TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
        surrogateVariable = paste0("SV", index),
        technicalFactor = terms[keep],
        pvalue = coerceNumericDnaEpico(table[[p_column]][keep]),
        stringsAsFactors = FALSE
    )
    }
    association_data <- if (length(rows)) {
    do.call(rbind, rows)
    } else {
    data.frame()
    }
    association_data <- association_data[
    is.finite(association_data$pvalue) & association_data$pvalue >= 0 &
        association_data$pvalue <= 1, ,
    drop = FALSE
    ]
    if (!nrow(association_data)) {
    return(NULL)
    }
    association_data$minusLog10P <- -log10(pmax(
    association_data$pvalue, .Machine$double.xmin
    ))
    association_data$label <- ifelse(
    association_data$pvalue < 0.001, "<0.001",
    sprintf("%.3f", association_data$pvalue)
    )
    ggplot2::ggplot(association_data, ggplot2::aes(
    x = surrogateVariable, y = technicalFactor, fill = minusLog10P
    )) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.2) +
    ggplot2::scale_fill_gradient(
        low = "#EAF4F6", high = "#B42318", name = "-log10(p)"
    ) +
    ggplot2::labs(
        title = NULL, x = "Surrogate variable", y = "Technical factor"
    ) +
    dnaEpicoModelPlotTheme()
}

#' @keywords internal
#' @noRd
svaEnmixMatrixPlotScales <- function(K) {
    list(axis = svaEnmixMatrixScale(K,
    multiplier = 3.4, minimum = 0.35,
    maximum = 0.7
    ), label = svaEnmixMatrixScale(K,
    multiplier = 3.8,
    minimum = 0.45, maximum = 0.8
    ), point = svaEnmixMatrixScale(K,
    multiplier = 4, minimum = 0.45, maximum = 0.8
    ), legend = svaEnmixMatrixScale(K,
    multiplier = 3.2, minimum = 0.35, maximum = 0.6
    ))
}

#' @keywords internal
#' @noRd
svaEnmixPlotMatrixLegends <- function(
    state, sentrix_position_levels,
    color_map, pch_map, legend_cex
) {
    if (isTRUE(state$showSentrixIDLegend)) {
    graphics::legend("topright",
        legend = state$sentrixIDLevels,
        col = color_map, pch = 15, title = "SentrixID", cex = legend_cex,
        bty = "n"
    )
    }

    graphics::legend(
    if (isTRUE(state$showSentrixIDLegend)) {
        "bottomright"
    } else {
        "topright"
    },
    legend = sentrix_position_levels, pch = pch_map,
    title = "SentrixPosition", cex = legend_cex, bty = "n"
    )
}

#' @keywords internal
#' @noRd
svaEnmixPlotMatrixPageTitle <- function(page, page_count) {
    title <- "Effects of Sentrix ID (color) and Sentrix Position (shape)"
    if (page_count <= 1L) {
    return(title)
    }

    paste(title, sprintf(
    "(page %d of %d: SV%d-SV%d by SV%d-SV%d)",
    page$page, page_count, min(page$rowIndices), max(page$rowIndices),
    min(page$colIndices), max(page$colIndices)
    ))
}

#' @keywords internal
#' @noRd
svaEnmixPlotMatrixPage <- function(analysisData, state, page) {
    sva <- analysisData$sva
    row_indices <- page$rowIndices
    col_indices <- page$colIndices
    page_panel_count <- max(length(row_indices), length(col_indices))
    sentrix_position_levels <- levels(analysisData$sentrixPosition)
    color_map <- (grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,
        "Dark2")))(length(state$sentrixIDLevels))
    pch_map <- seq_along(sentrix_position_levels)
    scales <- svaEnmixMatrixPlotScales(page_panel_count)
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(length(row_indices), length(col_indices)), family =
        "Times",
        las = 1, mar = svaEnmixMatrixPanelMargin(page_panel_count), oma = c(3.2,
            3.2, 0.5, 0.5), mgp = c(1.1, 0.35, 0), tcl = -0.18, cex.axis =
            scales$axis,
        cex.lab = scales$label)
    for (row_pos in seq_along(row_indices)) {
        for (col_pos in seq_along(col_indices)) {
            row_index <- row_indices[[row_pos]]
            col_index <- col_indices[[col_pos]]
            graphics::plot(sva[, col_index], sva[, row_index], col = color_map[
            analysisData$sentrixID],
                pch = pch_map[analysisData$sentrixPosition], cex = scales$point,
                xlab = if (row_pos == length(row_indices)) {
                    paste("SV", col_index)
                } else {
                    ""
                }, ylab = if (col_pos == 1L) {
                    paste("SV", row_index)
                } else {
                    ""
                }, main = "", xaxt = if (row_pos == length(row_indices)) {
                    "s"
                }
                else {
                    "n"
                }, yaxt = if (col_pos == 1L) {
                    "s"
                }
                else {
                    "n"
                })
            if (row_pos == 1L && col_pos == 1L) {
                svaEnmixPlotMatrixLegends(state, sentrix_position_levels,
            color_map,
                    pch_map, scales$legend)
            } } }
    invisible(NULL) }

#' @keywords internal
#' @noRd
svaEnmixSaveMatrixPage <- function(
    file, draw_fun, width, height,
    res
) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::tiff(
    filename = file, width = width, height = height,
    res = res, type = "cairo", compression = "lzw", bg = "white"
    )
    tryCatch(draw_fun(), finally = grDevices::dev.off())

    invisible(file)
}

#' @keywords internal
#' @noRd
svaEnmixRunMatrixPlot <- function(
    analysisData, state, display,
    file, res
) {
    if (is.null(file) && !isTRUE(display)) {
    return(invisible(NULL))
    }

    output_files <- NULL
    if (!is.null(file)) {
    output_files <- vapply(state$pages, function(page) {
        svaEnmixMatrixPageFile(
        file,
        page$page, state$pageCount
        )
    }, character(1))
    output_files <- unname(output_files)
    }

    if (!is.null(file)) {
    for (page in state$pages) {
        page_file <- svaEnmixMatrixPageFile(
        file, page$page,
        state$pageCount
        )
        svaEnmixSaveMatrixPage(
        file = page_file, draw_fun = function() {
            svaEnmixPlotMatrixPage(
            analysisData,
            state, page
            )
        }, width = state$width, height = state$height,
        res = res
        )
    }
    }

    if (isTRUE(display)) {
    for (page in state$pages) {
        svaEnmixPlotMatrixPage(analysisData, state, page)
    }
    }

    invisible(if (is.null(file)) NULL else output_files)
}

#' @keywords internal
#' @noRd
svaEnmixSavedPlotLogLines <- function(output) {
    if (is.null(output)) {
    return("Saved plot path:           not saved")
    }
    if (length(output) <= 1L) {
    return(paste("Saved plot path:          ", output))
    }

    c(
    paste("Saved plot paths:         ", length(output), "files"),
    paste("  ", output)
    )
}

prepareSvaPlotDnaEpico <- function(
    analysisData, plot, file, display, width, height, res,
    verbose, logPath
) {
    sva <- analysisData$sva
    if (plot %in% c("sentrix_id", "sentrix_position") && ncol(sva) < 2L) {
    stop("At least two surrogate variables are required for this plot.",
        call. = FALSE
    )
    }
    matrix_state <- list(width = width, height = height)
    if (identical(plot, "matrix")) {
    matrix_state <- svaEnmixMatrixPlotState(analysisData, width, height, res)
    if (!is.null(file) || isTRUE(display)) {
        svaEnmixLogMatrixPlotState(
        matrix_state, res, verbose, logPath
        )
    }
    }
    association <- if (identical(plot, "association")) {
    svaEnmixPlotAssociations(analysisData)
    } else {
    NULL
    }
    list(matrix = matrix_state, association = association)
}

omitEmptySvaAssociationDnaEpico <- function(
    plot, association, file, verbose, logPath
) {
    if (!identical(plot, "association") || !is.null(association)) {
    return(FALSE)
    }
    if (!is.null(file) && file.exists(file)) {
    status <- unlink(file, force = TRUE)
    if (status != 0L || file.exists(file)) {
        stop("Could not remove the stale SVA association figure: ",
        file,
        call. = FALSE
        )
    }
    }
    emitLogMinfiEwasWater(c(
    "SVA plot type:             association",
    paste0(
        "SVA association figure:   omitted; no estimable ",
        "technical-factor associations."
    ),
    "============================================================"
    ), verbose = verbose, log_path = logPath)
    TRUE
}

runSelectedSvaPlotDnaEpico <- function(
    analysisData, plot, state, display, file, res
) {
    if (identical(plot, "matrix")) {
    return(svaEnmixRunMatrixPlot(
        analysisData = analysisData, state = state$matrix,
        display = display, file = file, res = res
    ))
    }
    draw_fun <- switch(plot,
    sentrix_id = function() svaEnmixPlotSentrixID(analysisData),
    sentrix_position = function() svaEnmixPlotSentrixPosition(analysisData),
    association = function() {
        drawPlotObjectMinfiEwasWater(
        state$association
        )
    }
    )
    runPlotMinfiEwasWater(
    draw_fun = draw_fun, display = display, file = file,
    width = state$matrix$width, height = state$matrix$height, res = res
    )
}

#' Plot surrogate variables for svaEnmix
#'
#' @param analysisData Object returned by `analyzeSvaEnmix()`.
#' @param plot Character. Plot type: `'sentrix_id'`, `'sentrix_position'`, or
#'   `'matrix'`, or `'association'`.
#' @param display Logical. If `TRUE`, draw the plot on the active graphics
#'   device.
#' @param file Character or `NULL`. TIFF file path used for saved output.
#' @param width Integer. Plot width in pixels when `file` is supplied.
#' @param height Integer. Plot height in pixels when `file` is supplied.
#' @param res Integer. TIFF resolution in DPI when `file` is supplied.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return Invisibly returns saved TIFF path(s), otherwise `NULL`. Association
#'   output is omitted when no technical-factor association is estimable.
#'
#' @description
#' Draw one of the standard surrogate-variable plots used by `svaEnmix()`.
#' An empty association result is recorded in the log without creating a
#' placeholder figure.
#'
#' @examples
#' ex <- dnaEPICO:::exampleSvaAnalysisStateDnaEpico()
#' plotSvaEnmix(
#'   analysisData = ex$analysisData,
#'   plot = "sentrix_id",
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#'
#' @export
plotSvaEnmix <- function(
    analysisData,
    plot = c("sentrix_id", "sentrix_position", "matrix", "association"),
    display = FALSE, file = NULL, width = 2000L, height = 1000L,
    res = 150L, verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_plotSvaEnmix.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs, log_dir = log_dir, log_file = log_file
    )
    plot <- match.arg(plot)
    state <- prepareSvaPlotDnaEpico(
    analysisData, plot, file, display, width, height, res,
    verbose, log_path
    )
    if (omitEmptySvaAssociationDnaEpico(
    plot, state$association, file, verbose, log_path
    )) {
    return(invisible(NULL))
    }
    output <- runSelectedSvaPlotDnaEpico(
    analysisData, plot, state, display, file, res
    )
    emitLogMinfiEwasWater(c(
    paste("SVA plot type:            ", plot),
    paste("Display plot:             ", display),
    svaEnmixSavedPlotLogLines(output),
    "============================================================"
    ), verbose = verbose, log_path = log_path)
    invisible(output)
}

#' @keywords internal
#' @noRd
readPhenotypeFileSvaEnmix <- function(file, sepType = NULL) {
    separator <- resolveSeparatorMinfiEwasWater(sepType)
    if (is.null(separator)) {
    return(utils::read.csv(file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    ))
    }
    utils::read.table(file,
    header = TRUE, sep = separator, stringsAsFactors = FALSE,
    check.names = FALSE, quote = "\"", comment.char = ""
    )
}

validateReplacementFramesSvaEnmix <- function(original, updated, SampleID) {
    if (!is.data.frame(original) || !is.data.frame(updated)) {
    stop("The original and updated phenotype objects must be data frames.",
        call. = FALSE
    )
    }
    if (!(SampleID %in% colnames(original)) ||
    !(SampleID %in% colnames(updated))) {
    stop("SampleID column is missing from the phenotype replacement: ",
        SampleID,
        call. = FALSE
    )
    }
    if (nrow(original) != nrow(updated)) {
    stop("The updated phenotype has ", nrow(updated), " rows; expected ",
        nrow(original), ". The original file was not changed.",
        call. = FALSE
    )
    }
    original_ids <- validateSampleIdentifiersDnaEpico(
    original[[SampleID]], paste0(
        "Original phenotype column '",
        SampleID, "'"
    )
    )
    updated_ids <- validateSampleIdentifiersDnaEpico(
    updated[[SampleID]], paste0(
        "Updated phenotype column '",
        SampleID, "'"
    )
    )
    if (!identical(original_ids, updated_ids)) {
    stop("Sample identifiers or their order changed; the original file ",
        "was not changed.",
        call. = FALSE
    )
    }
    invisible(NULL)
}

validateOriginalColumnsSvaEnmix <- function(original, updated) {
    if (anyDuplicated(colnames(updated))) {
    stop("The updated phenotype contains duplicated column names.",
        call. = FALSE
    )
    }
    missing <- setdiff(colnames(original), colnames(updated))
    if (length(missing)) {
    missing_text <- paste(missing, collapse = ", ")
    stop(sprintf(
        "The updated phenotype is missing original column(s): %s",
        missing_text
    ), call. = FALSE)
    }
    for (column in colnames(original)) {
    if (!isTRUE(all.equal(
        original[[column]], updated[[column]],
        check.attributes = FALSE
    ))) {
        stop("Original phenotype column changed: ", column,
        call. = FALSE
        )
    }
    }
    invisible(NULL)
}

validateReplacementPcColumnsSvaEnmix <- function(original, updated, pcColumns) {
    if (!length(pcColumns) || anyNA(pcColumns) || any(!nzchar(pcColumns))) {
    stop("The SVA output does not contain named PC columns.",
        call. = FALSE
    )
    }
    if (anyDuplicated(pcColumns) ||
    !all(grepl("^PC[1-9][0-9]*$", pcColumns))) {
    stop("SVA columns must have unique names PC1, PC2, and so on.",
        call. = FALSE
    )
    }
    missing <- setdiff(pcColumns, colnames(updated))
    if (length(missing)) {
    missing_text <- paste(missing, collapse = ", ")
    stop(sprintf(
        "The updated phenotype is missing PC column(s): %s",
        missing_text
    ), call. = FALSE)
    }
    if (any(pcColumns %in% colnames(original))) {
    stop("PC column names already exist in the original phenotype file.",
        call. = FALSE
    )
    }
    if (!all(vapply(updated[pcColumns], is.numeric, logical(1)))) {
    stop("All appended PC columns must be numeric.", call. = FALSE)
    }
    if (!all(is.finite(as.matrix(updated[pcColumns])))) {
    stop("Appended PC columns must contain only finite values.",
        call. = FALSE
    )
    }
    invisible(NULL)
}

#' @keywords internal
#' @noRd
validatePhenotypeReplacementSvaEnmix <- function(
    original, updated, SampleID, pcColumns
) {
    validateReplacementFramesSvaEnmix(original, updated, SampleID)
    validateOriginalColumnsSvaEnmix(original, updated)
    validateReplacementPcColumnsSvaEnmix(original, updated, pcColumns)
    invisible(TRUE)
}

#' @keywords internal
#' @noRd
writePhenotypeFileSvaEnmix <- function(data, file, sepType = NULL) {
    separator <- resolveSeparatorMinfiEwasWater(sepType)
    if (is.null(separator)) {
    utils::write.csv(data, file = file, row.names = FALSE)
    } else {
    utils::write.table(data,
        file = file, sep = separator,
        row.names = FALSE, col.names = TRUE, quote = TRUE,
        na = "NA"
    )
    }
    invisible(file)
}

#' @keywords internal
#' @noRd
stagePhenotypeReplacementSvaEnmix <- function(
    mergedPheno, original, tempFile, SampleID, pcColumns, sepType
) {
    writePhenotypeFileSvaEnmix(mergedPheno, tempFile, sepType = sepType)
    written <- readPhenotypeFileSvaEnmix(tempFile, sepType = sepType)
    validatePhenotypeReplacementSvaEnmix(
    original = original, updated = written,
    SampleID = SampleID, pcColumns = pcColumns
    )
    invisible(NULL)
}

activatePhenotypeReplacementSvaEnmix <- function(
    phenoFile, tempFile, backupFile
) {
    if (!file.rename(phenoFile, backupFile)) {
    stop("Could not create a rollback copy of the phenotype file.",
        call. = FALSE
    )
    }
    if (file.rename(tempFile, phenoFile)) {
    return(invisible(NULL))
    }
    restored <- file.rename(backupFile, phenoFile)
    if (!isTRUE(restored)) {
    failure_message <- paste(
        "Phenotype replacement failed and the original could not be",
        "restored from"
    )
    stop(sprintf(
        "%s: %s",
        failure_message,
        backupFile
    ), call. = FALSE)
    }
    stop("Phenotype replacement failed; the original file was restored.",
    call. = FALSE
    )
}

validateCommittedPhenotypeSvaEnmix <- function(
    phenoFile, backupFile, original, SampleID, pcColumns, sepType
) {
    error <- tryCatch(
    {
        final <- readPhenotypeFileSvaEnmix(phenoFile, sepType = sepType)
        validatePhenotypeReplacementSvaEnmix(
        original = original, updated = final,
        SampleID = SampleID, pcColumns = pcColumns
        )
        NULL
    },
    error = identity
    )
    if (is.null(error)) {
    return(invisible(NULL))
    }
    unlink(phenoFile)
    restored <- file.rename(backupFile, phenoFile)
    if (!isTRUE(restored)) {
    failure_message <- paste(
        "Final phenotype validation failed and the original could",
        "not be restored from"
    )
    stop(sprintf(
        "%s: %s",
        failure_message,
        backupFile
    ), call. = FALSE)
    }
    restored_message <- paste(
        "Final phenotype validation failed; the original file was",
        "restored"
    )
    stop(sprintf(
    "%s: %s",
    restored_message, conditionMessage(error)
    ), call. = FALSE)
}

#' @keywords internal
#' @noRd
replacePhenotypeFileSvaEnmix <- function(mergedPheno,
    phenoFile, SampleID, pcColumns, sepType = NULL) {
    if (!file.exists(phenoFile)) {
        stop("Phenotype file does not exist: ",
            phenoFile, call. = FALSE)
    }
    original <- readPhenotypeFileSvaEnmix(phenoFile,
        sepType = sepType)
    validatePhenotypeReplacementSvaEnmix(original = original,
        updated = mergedPheno, SampleID = SampleID,
        pcColumns = pcColumns)
    output_dir <- dirname(phenoFile)
    temp_file <- tempfile(pattern = paste0(".",
        basename(phenoFile), ".sva-"), tmpdir = output_dir,
        fileext = ".tmp")
    backup_file <- tempfile(pattern = paste0(".",
        basename(phenoFile), ".backup-"),
        tmpdir = output_dir, fileext = ".tmp")
    backup_active <- FALSE
    replacement_complete <- FALSE
    on.exit({
        if (file.exists(temp_file)) unlink(temp_file)
        if (backup_active && file.exists(backup_file)) {
            if (file.exists(phenoFile)) unlink(phenoFile)
            file.rename(backup_file, phenoFile)
        }
        if (replacement_complete && file.exists(backup_file)) {
            unlink(backup_file)
        }
    }, add = TRUE)
    stagePhenotypeReplacementSvaEnmix(mergedPheno,
        original, temp_file, SampleID, pcColumns,
        sepType)
    activatePhenotypeReplacementSvaEnmix(phenoFile,
        temp_file, backup_file)
    backup_active <- TRUE
    validateCommittedPhenotypeSvaEnmix(phenoFile,
        backup_file, original, SampleID,
        pcColumns, sepType)
    if (file.exists(backup_file) && unlink(backup_file) !=
        0L) {
        stop("The phenotype file was replaced, but its temporary backup could ",
            "not be removed.", call. = FALSE)
    }
    backup_active <- FALSE
    replacement_complete <- TRUE
    invisible(phenoFile)
}

#' Write svaEnmix outputs to disk
#'
#' @param svaData Object returned by `estimateSvaEnmixControls()`.
#' @param mergedPheno Phenotype data frame returned by `mergeSvaTargetsEnmix()`.
#' @param analysisData Optional object returned by `analyzeSvaEnmix()`.
#' @param phenoFile Character or `NULL`. Existing phenotype file replaced with
#'   `mergedPheno` after validation.
#' @param SampleID Character. Phenotype column containing sample identifiers.
#' @param sepType Character or `NULL`. Field separator used by `phenoFile`.
#' @param dataBaseDir Character. Base directory used for saved data outputs.
#' @param rBaseDir Character. Base directory used for saved `.RData` outputs.
#' @param scriptLabel Character. Label used to create the output subdirectory.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_svaEnmix_paths'` containing the paths
#'   written to disk.
#'
#' @description
#' Write the SVA matrix, phenotype, and model-summary outputs.
#'
#' @examples
#' ex <- dnaEPICO:::exampleSvaAnalysisStateDnaEpico()
#' temp_dir <- tempdir()
#' pheno_file <- file.path(temp_dir, "phenoLC.csv")
#' utils::write.csv(ex$targets, pheno_file, row.names = FALSE)
#' output_paths <- writeSvaEnmixOutputs(
#'   svaData = list(sva = ex$sva),
#'   mergedPheno = ex$mergedPheno,
#'   analysisData = ex$analysisData,
#'   phenoFile = pheno_file,
#'   SampleID = "Sample_Name",
#'   dataBaseDir = file.path(temp_dir, "data"),
#'   rBaseDir = file.path(temp_dir, "rData"),
#'   scriptLabel = "svaEnmixExample",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(output_paths)
#'
#' @export
writeSvaEnmixOutputs <- function(svaData, mergedPheno, analysisData = NULL,
    phenoFile = NULL, SampleID = "Sample_Name", sepType = NULL, dataBaseDir =
        "data",
    rBaseDir = "rData", scriptLabel = "svaEnmix", verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_writeSvaEnmixOutputs.txt") {
    log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir,
        log_file = log_file)
    data_dir <- file.path(dataBaseDir, scriptLabel)
    r_dir <- file.path(rBaseDir, scriptLabel)
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(r_dir, recursive = TRUE, showWarnings = FALSE)
    sva_rdata_path <- file.path(r_dir, "svaMatrix.RData")
    sva_csv_path <- file.path(data_dir, "svaMatrix.csv")
    pheno_output_path <- phenoFile
    saveNamedObjectMinfiEwasWater(svaData$sva, "sva", sva_rdata_path)
    utils::write.csv(svaData$sva, sva_csv_path, row.names = TRUE)
    if (!is.null(analysisData)) { for (i in seq_len(analysisData$K)) {
            utils::capture.output(summary(analysisData$fullModels[[i]]),
                file = file.path(data_dir, paste0("summary_full_sva",
                    i, ".txt")))
            if (length(analysisData$droptermSteps[[i]]) > 0L) {
                for (step in analysisData$droptermSteps[[i]]) {
                    utils::capture.output(step$dropterm, file = file.path(
            data_dir,
                    paste0("dropterm_step_sva", i, ".txt")), append = TRUE)
                    utils::capture.output(step$summary, file = file.path(
            data_dir,
                    paste0("dropterm_model_sva", i, ".txt")), append = TRUE)
                } }
            utils::capture.output(analysisData$anovaFull[[i]], file =
            file.path(data_dir, paste0("anova_full_sva", i, ".txt")))
            utils::capture.output(analysisData$anovaReduced[[i]], file =
            file.path(data_dir,
                paste0("anova_reduced_sva", i, ".txt")))
        } }; if (!is.null(pheno_output_path)) {
        replacePhenotypeFileSvaEnmix(mergedPheno = mergedPheno, phenoFile =
            pheno_output_path,
            SampleID = SampleID, pcColumns = colnames(svaData$sva), sepType =
            sepType) }
    emitLogMinfiEwasWater(c(paste("SVA Matrix RData saved to: ",
        sva_rdata_path),
        paste("SVA Matrix CSV saved to:   ", sva_csv_path), paste(
            "Updated phenotype + SVA:   ",
            if (is.null(pheno_output_path)) {
                "not saved" } else { pheno_output_path
            }), "============================================================"),
        verbose = verbose, log_path = log_path)
    structure(list(svaRData = sva_rdata_path, svaCSV = sva_csv_path,
        phenoWithSva = pheno_output_path, dataDir = data_dir, rDir = r_dir),
        class = "dnaEPICO_svaEnmix_paths") }
