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
#'     RGSet = ex$RGSet,
#'     ctrlSvaPercVar = 0.5,
#'     ctrlSvaFlag = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' sva_data$K
#'
#' @export
estimateSvaEnmixControls <- function(
    RGSet, ctrlSvaPercVar = 0.9,
    ctrlSvaFlag = 1, verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_estimateSvaEnmixControls.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    ctrlSvaPercVar <- validateProbabilityDnaEpico(
        ctrlSvaPercVar,
        "ctrlSvaPercVar"
    )
    if (length(ctrlSvaFlag) != 1L || !is.numeric(ctrlSvaFlag) ||
        is.na(ctrlSvaFlag) || !is.finite(ctrlSvaFlag) || !(ctrlSvaFlag %in%
        c(1, 2))) {
        stop("ctrlSvaFlag must be 1 (percvar) or 2 (one requested component).",
            call. = FALSE
        )
    }

    sva <- ENmix::ctrlsva(
        rgSet = RGSet, percvar = ctrlSvaPercVar,
        flag = ctrlSvaFlag
    )
    sva <- as.matrix(sva)
    nan_sva_converted <- sum(is.nan(sva))
    if (nan_sva_converted > 0L) {
        sva[is.nan(sva)] <- NA_real_
    }
    sample_names <- validateSampleIdentifiersDnaEpico(
        svaEnmixGetRGSetSampleNames(RGSet),
        "RGSet sample identifiers"
    )
    if (!is.numeric(sva) || nrow(sva) != length(sample_names) ||
        ncol(sva) == 0L || anyNA(sva) || any(!is.finite(sva)) ||
        is.null(rownames(sva))) {
        if (nan_sva_converted > 0L) {
            stop("ENmix::ctrlsva() must return a finite numeric matrix with one named row per RGSet sample; ",
                nan_sva_converted,
                    " upstream NaN value(s) were converted to NA.",
                call. = FALSE
            )
        }
        stop("ENmix::ctrlsva() must return a finite numeric matrix with one named row per RGSet sample.",
            call. = FALSE
        )
    }
    sample_match <- matchSampleIdentifiersDnaEpico(
        query = sample_names,
        reference = rownames(sva), queryLabel = "RGSet sample identifiers",
        referenceLabel = "ENmix SVA row names", requireSameSet = TRUE
    )
    sva <- sva[sample_match, , drop = FALSE]
    colnames(sva) <- paste0("PC", seq_len(ncol(sva)))

    emitLogMinfiEwasWater(
        c(
            paste(
                "ctrlSva percvar:          ",
                ctrlSvaPercVar
            ), paste(
                "ctrlSva flag:             ",
                ctrlSvaFlag
            ), paste("ctrlSva NaN to NA:        ", nan_sva_converted),
            paste("Number of surrogate variables:", ncol(sva)),
                "Surrogate variables matrix (first few rows):",
            previewLinesMinfiEwasWater(sva),
                "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(
        list(
            sva = sva, K = ncol(sva), ctrlSvaPercVar = ctrlSvaPercVar,
            ctrlSvaFlag = ctrlSvaFlag, nanConverted = nan_sva_converted
        ),
        class = "dnaEPICO_svaEnmix_sva"
    )
}

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
#'     targets = ex$targets,
#'     sva = ex$sva,
#'     SampleID = "Sample_Name",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' colnames(merged_pheno)[seq_len(4)]
#'
#' @export
mergeSvaTargetsEnmix <- function(
    targets, sva, SampleID = "Sample_Name",
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_mergeSvaTargetsEnmix.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    if (!(SampleID %in% colnames(targets))) {
        stop("SampleID column not found in targets: ", SampleID,
            call. = FALSE
        )
    }

    if (!is.matrix(sva) || !is.numeric(sva) || ncol(sva) == 0L) {
        stop("sva must be a numeric matrix with at least one column.",
            call. = FALSE
        )
    }
    if (is.null(rownames(sva))) {
        stop("sva must have row names that match the phenotype SampleID column.",
            call. = FALSE
        )
    }
    if (anyNA(sva) || any(!is.finite(sva))) {
        stop("sva contains missing or non-finite values.", call. = FALSE)
    }

    sva_data <- data.frame(
        sample_id = rownames(sva), sva, row.names = NULL,
        check.names = FALSE, stringsAsFactors = FALSE
    )
    colnames(sva_data)[1] <- SampleID
    surrogate_columns <- setdiff(colnames(sva_data), SampleID)
    duplicate_columns <- intersect(colnames(targets), surrogate_columns)
    if (length(duplicate_columns) > 0L) {
        stop("The phenotype table already contains output columns that would be duplicated: ",
            paste(duplicate_columns, collapse = ", "),
                ". Start from the pre-SVA phenotype file or remove the existing columns explicitly.",
            call. = FALSE
        )
    }

    match_idx <- matchSampleIdentifiersDnaEpico(
        query = targets[[SampleID]],
        reference = sva_data[[SampleID]], queryLabel = paste0(
            "Phenotype column '",
            SampleID, "'"
        ), referenceLabel = "SVA row names",
        requireSameSet = TRUE
    )
    merged_pheno <- cbind(targets, sva_data[match_idx, surrogate_columns,
        drop = FALSE
    ])

    emitLogMinfiEwasWater(
        c(paste(
            "Merged phenotype rows:     ",
            nrow(merged_pheno)
        ), paste(
            "Merged phenotype columns:  ",
            ncol(merged_pheno)
        ), "============================================================"),
        verbose = verbose, log_path = log_path
    )

    merged_pheno
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
#'     sva = ex$sva,
#'     RGSet = ex$RGSet,
#'     SentrixIDColumn = "Sentrix_ID",
#'     SentrixPositionColumn = "Sentrix_Position",
#'     verbose = FALSE,
#'     logs = FALSE
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
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    sample_names <- svaEnmixGetRGSetSampleNames(RGSet)

    if (length(sample_names) == 0L) {
        stop("Could not determine sample names from the loaded RGSet.",
            call. = FALSE
        )
    }

    if (!is.matrix(sva) || !is.numeric(sva) || ncol(sva) == 0L) {
        stop("sva must be a numeric matrix with at least one surrogate variable.",
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

    match_idx <- matchSampleIdentifiersDnaEpico(
        query = sample_names,
        reference = rownames(sva), queryLabel = "RGSet sample identifiers",
        referenceLabel = "SVA row names", requireSameSet = TRUE
    )

    sva <- sva[match_idx, , drop = FALSE]
    col_data <- svaEnmixGetRGSetColData(RGSet)
    missing_sentrix <- setdiff(
        c(SentrixIDColumn, SentrixPositionColumn),
        colnames(col_data)
    )
    if (length(missing_sentrix) > 0L) {
        stop("RGSet column data is missing: ", paste(missing_sentrix,
            collapse = ", "
        ), call. = FALSE)
    }
    sentrix_id <- as.factor(col_data[[SentrixIDColumn]])
    sentrix_position <- as.factor(col_data[[SentrixPositionColumn]])
    if (anyNA(sentrix_id) || anyNA(sentrix_position)) {
        stop("Sentrix chip and position columns cannot contain missing values.",
            call. = FALSE
        )
    }
    K <- ncol(sva)

    technical_data <- data.frame(SentrixID = sentrix_id,
        SentrixPosition = sentrix_position)
    technical_terms <- colnames(technical_data)[vapply(
        technical_data,
        function(values) nlevels(droplevels(values)) > 1L, logical(1)
    )]
    technical_formula <- if (length(technical_terms) == 0L) {
        stats::as.formula("sva_value ~ 1")
    } else {
        stats::reformulate(termlabels = technical_terms, response = "sva_value")
    }

    full_models <- lapply(seq_len(K), function(i) {
        model_data <- technical_data
        model_data$sva_value <- sva[, i]
        stats::lm(technical_formula, data = model_data)
    })

    reduced_models <- vector("list", K)
    dropterm_steps <- vector("list", K)
    warning_state <- new.env(parent = emptyenv())
    warning_state$messages <- character(0)
    capture_analysis_warnings <- function(expression) {
        withCallingHandlers(expression, warning = function(condition) {
            warning_state$messages <- unique(c(
                warning_state$messages,
                conditionMessage(condition)
            ))
            invokeRestart("muffleWarning")
        })
    }
    for (i in seq_len(K)) {
        model_tmp <- full_models[[i]]
        model_steps <- list()

        repeat {
            drop_table <- capture_analysis_warnings(stats::drop1(model_tmp,
                test = "F"
            ))
            drop_p <- drop_table$`Pr(F)`
            valid_idx <- which(!is.na(drop_p))

            if (length(valid_idx) == 0L) {
                break
            }

            max_idx <- valid_idx[which.max(drop_p[valid_idx])]
            max_p <- drop_p[[max_idx]]
            term_to_drop <- rownames(drop_table)[[max_idx]]

            if (!is.finite(max_p) || max_p <= 0.05 || identical(
                term_to_drop,
                "<none>"
            )) {
                break
            }

            model_tmp <- stats::update(model_tmp, paste(
                ". ~ . -",
                term_to_drop
            ))
            model_steps[[length(model_steps) + 1L]] <- list(
                dropterm = drop_table,
                summary = summary(model_tmp)
            )
        }

        reduced_models[[i]] <- model_tmp
        dropterm_steps[[i]] <- model_steps
    }

    anova_full <- lapply(full_models,
        function(model) capture_analysis_warnings(stats::anova(model)))
    anova_reduced <- lapply(reduced_models,
        function(model) capture_analysis_warnings(stats::anova(model)))
    analysis_warnings <- warning_state$messages

    emitLogMinfiEwasWater(
        c(
            paste(
                "Number of surrogate variables (K):",
                K
            ), paste("SentrixID class:            ", class(sentrix_id)[1]),
            paste("SentrixID unique levels:    ", length(unique(sentrix_id))),
            previewLinesMinfiEwasWater(table(sentrix_id)), paste(
                "SentrixPosition class:      ",
                class(sentrix_position)[1]
            ), paste(
                "SentrixPosition unique levels:",
                length(unique(sentrix_position))
            ), previewLinesMinfiEwasWater(table(sentrix_position)),
            paste("Technical terms modelled:     ",
                if (length(technical_terms)) {
                paste(technical_terms,
                    collapse = ", "
                )
            } else {
                "intercept only"
            }), "First row of SVA matrix:",
            previewLinesMinfiEwasWater(sva[1, , drop = FALSE]), paste(
                "Sample names in SVA matrix:",
                paste(rownames(sva)[seq_len(min(5L, nrow(sva)))],
                    collapse = ", "
                )
            ), paste(
                "Sample names in colData(RGSet):",
                paste(sample_names[seq_len(min(5L, length(sample_names)))],
                    collapse = ", "
                )
            ), paste(
                "Model-selection warnings:    ",
                length(analysis_warnings)
            ), if (length(analysis_warnings) >
                0L) {
                analysis_warnings
            } else {
                character(0)
            }, "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(
        list(
            sva = sva, K = K, sentrixID = sentrix_id,
            sentrixPosition = sentrix_position,
                technicalTerms = technical_terms,
            fullModels = full_models, reducedModels = reduced_models,
            droptermSteps = dropterm_steps, anovaFull = anova_full,
            anovaReduced = anova_reduced, warnings = analysis_warnings
        ),
        class = "dnaEPICO_svaEnmix_analysis"
    )
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
        return(c(width = as.integer(ceiling(width)),
            height = as.integer(ceiling(height))))
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
    color_map <- grDevices::rainbow(length(levels(analysisData$sentrixID)))
    sva <- analysisData$sva
    graphics::plot(sva[, 1], sva[, 2],
        col = color_map[analysisData$sentrixID],
        pch = 16, xlab = "Surrogate Variable 1 (PC1)",
            ylab = "Surrogate Variable 2 (PC2)",
        main = "Surrogate Variables Colored by Chip (SentrixID)"
    )
    graphics::legend("topright",
        legend = levels(analysisData$sentrixID),
        col = color_map, pch = 16, title = "SentrixID", cex = 0.6
    )
}

#' @keywords internal
#' @noRd
svaEnmixPlotSentrixPosition <- function(analysisData) {
    color_map <-
        grDevices::rainbow(length(levels(analysisData$sentrixPosition)))
    sva <- analysisData$sva
    graphics::plot(sva[, 1], sva[, 2],
        col = color_map[analysisData$sentrixPosition],
        pch = 16, xlab = "Surrogate Variable 1 (PC1)",
            ylab = "Surrogate Variable 2 (PC2)",
        main = "Surrogate Variables Colored by Sentrix Position"
    )
    graphics::legend("topright",
        legend = levels(analysisData$sentrixPosition),
        col = color_map, pch = 16, title = "SentrixPosition",
        cex = 0.6
    )
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
    color_map <- grDevices::rainbow(length(state$sentrixIDLevels))
    pch_map <- seq_along(sentrix_position_levels)
    scales <- svaEnmixMatrixPlotScales(page_panel_count)

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(
        mfrow = c(length(row_indices), length(col_indices)),
        family = "Times", las = 1,
            mar = svaEnmixMatrixPanelMargin(page_panel_count),
        oma = c(3.2, 3.2, 2, 0.5), mgp = c(1.1, 0.35, 0), tcl = -0.18,
        cex.axis = scales$axis, cex.lab = scales$label
    )

    for (row_pos in seq_along(row_indices)) {
        for (col_pos in seq_along(col_indices)) {
            row_index <- row_indices[[row_pos]]
            col_index <- col_indices[[col_pos]]
            graphics::plot(sva[, col_index], sva[, row_index],
                col = color_map[analysisData$sentrixID],
                    pch = pch_map[analysisData$sentrixPosition],
                cex = scales$point, xlab = if (row_pos == length(row_indices)) {
                    paste("SV", col_index)
                } else {
                    ""
                }, ylab = if (col_pos == 1L) {
                    paste("SV", row_index)
                } else {
                    ""
                }, main = "", xaxt = if (row_pos == length(row_indices)) {
                    "s"
                } else {
                    "n"
                }, yaxt = if (col_pos == 1L) {
                    "s"
                } else {
                    "n"
                }
            )
            if (row_pos == 1L && col_pos == 1L) {
                svaEnmixPlotMatrixLegends(
                    state, sentrix_position_levels,
                    color_map, pch_map, scales$legend
                )
            }
        }
    }

    graphics::mtext(svaEnmixPlotMatrixPageTitle(page, state$pageCount),
        side = 3, outer = TRUE, line = 0.5, cex = 0.9
    )
}

#' @keywords internal
#' @noRd
svaEnmixSaveMatrixPage <- function(
    file, draw_fun, width, height,
    res
) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::tiff(
        filename = file, width = width, height = height,
        res = res, type = "cairo"
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

#' Plot surrogate variables for svaEnmix
#'
#' @param analysisData Object returned by `analyzeSvaEnmix()`.
#' @param plot Character. Plot type: `'sentrix_id'`, `'sentrix_position'`, or
#'   `'matrix'`.
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
#' @return Invisibly returns saved TIFF path(s), otherwise `NULL`.
#'
#' @description
#' Draw one of the standard surrogate-variable plots used by `svaEnmix()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleSvaAnalysisStateDnaEpico()
#' plotSvaEnmix(
#'     analysisData = ex$analysisData,
#'     plot = "sentrix_id",
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#'
#' @export
plotSvaEnmix <- function(
    analysisData, plot = c(
        "sentrix_id",
        "sentrix_position", "matrix"
    ), display = FALSE, file = NULL,
    width = 2000L, height = 1000L, res = 150L, verbose = FALSE,
    logs = FALSE, log_dir = NULL, log_file = "log_plotSvaEnmix.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )
    plot <- match.arg(plot)
    sva <- analysisData$sva
    matrix_state <- list(width = width, height = height)
    will_draw <- !is.null(file) || isTRUE(display)

    if (plot %in% c("sentrix_id", "sentrix_position") && ncol(sva) <
        2L) {
        stop("At least two surrogate variables are required for this plot.",
            call. = FALSE
        )
    }

    if (identical(plot, "matrix")) {
        matrix_state <- svaEnmixMatrixPlotState(
            analysisData,
            width, height, res
        )
        if (will_draw) {
            svaEnmixLogMatrixPlotState(
                matrix_state, res, verbose,
                log_path
            )
        }
    }

    draw_fun <- switch(plot,
        sentrix_id = function() svaEnmixPlotSentrixID(analysisData),
        sentrix_position = function() svaEnmixPlotSentrixPosition(analysisData),
        matrix = NULL
    )

    if (identical(plot, "matrix")) {
        output <- svaEnmixRunMatrixPlot(
            analysisData = analysisData,
            state = matrix_state, display = display, file = file,
            res = res
        )
    } else {
        output <- runPlotMinfiEwasWater(
            draw_fun = draw_fun,
            display = display, file = file, width = matrix_state$width,
            height = matrix_state$height, res = res
        )
    }

    emitLogMinfiEwasWater(
        c(
            paste(
                "SVA plot type:            ",
                plot
            ), paste("Display plot:             ", display),
            svaEnmixSavedPlotLogLines(output),
                "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

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

#' @keywords internal
#' @noRd
validatePhenotypeReplacementSvaEnmix <- function(
    original, updated,
    SampleID, pcColumns
) {
    if (!is.data.frame(original) || !is.data.frame(updated)) {
        stop("The original and updated phenotype objects must be data frames.",
            call. = FALSE
        )
    }
    if (!(SampleID %in% colnames(original)) || !(SampleID %in%
        colnames(updated))) {
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
        original[[SampleID]],
        paste0("Original phenotype column '", SampleID, "'")
    )
    updated_ids <- validateSampleIdentifiersDnaEpico(
        updated[[SampleID]],
        paste0("Updated phenotype column '", SampleID, "'")
    )
    if (!identical(original_ids, updated_ids)) {
        stop("Sample identifiers or their order changed; the original file was not changed.",
            call. = FALSE
        )
    }
    if (anyDuplicated(colnames(updated))) {
        stop("The updated phenotype contains duplicated column names.",
            call. = FALSE
        )
    }
    missing_original <- setdiff(colnames(original), colnames(updated))
    if (length(missing_original) > 0L) {
        stop("The updated phenotype is missing original column(s): ",
            paste(missing_original, collapse = ", "),
            call. = FALSE
        )
    }
    for (column in colnames(original)) {
        if (!isTRUE(all.equal(original[[column]], updated[[column]],
            check.attributes = FALSE
        ))) {
            stop("Original phenotype column changed: ", column,
                call. = FALSE
            )
        }
    }
    if (length(pcColumns) == 0L || anyNA(pcColumns) ||
        any(!nzchar(pcColumns))) {
        stop("The SVA output does not contain named PC columns.",
            call. = FALSE
        )
    }
    if (anyDuplicated(pcColumns) || !all(grepl(
        "^PC[1-9][0-9]*$",
        pcColumns
    ))) {
        stop("SVA columns must have unique names PC1, PC2, and so on.",
            call. = FALSE
        )
    }
    missing_pc <- setdiff(pcColumns, colnames(updated))
    if (length(missing_pc) > 0L) {
        stop("The updated phenotype is missing PC column(s): ",
            paste(missing_pc, collapse = ", "),
            call. = FALSE
        )
    }
    if (any(pcColumns %in% colnames(original))) {
        stop("PC column names already exist in the original phenotype file.",
            call. = FALSE
        )
    }
    pc_numeric <- vapply(updated[pcColumns], is.numeric, logical(1))
    if (!all(pc_numeric)) {
        stop("All appended PC columns must be numeric.", call. = FALSE)
    }
    if (!all(is.finite(as.matrix(updated[pcColumns])))) {
        stop("Appended PC columns must contain only finite values.",
            call. = FALSE
        )
    }
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
replacePhenotypeFileSvaEnmix <- function(
    mergedPheno, phenoFile,
    SampleID, pcColumns, sepType = NULL
) {
    if (!file.exists(phenoFile)) {
        stop("Phenotype file does not exist: ", phenoFile, call. = FALSE)
    }
    original <- readPhenotypeFileSvaEnmix(phenoFile, sepType = sepType)
    validatePhenotypeReplacementSvaEnmix(
        original = original,
        updated = mergedPheno, SampleID = SampleID, pcColumns = pcColumns
    )

    output_dir <- dirname(phenoFile)
    temp_file <- tempfile(pattern = paste0(
        ".", basename(phenoFile),
        ".sva-"
    ), tmpdir = output_dir, fileext = ".tmp")
    backup_file <- tempfile(pattern = paste0(
        ".", basename(phenoFile),
        ".backup-"
    ), tmpdir = output_dir, fileext = ".tmp")
    backup_active <- FALSE
    replacement_complete <- FALSE
    on.exit(
        {
            if (file.exists(temp_file)) unlink(temp_file)
            if (isTRUE(backup_active) && file.exists(backup_file)) {
                if (file.exists(phenoFile)) unlink(phenoFile)
                file.rename(backup_file, phenoFile)
            }
            if (isTRUE(replacement_complete) &&
                file.exists(backup_file)) unlink(backup_file)
        },
        add = TRUE
    )

    writePhenotypeFileSvaEnmix(mergedPheno, temp_file, sepType = sepType)
    written <- readPhenotypeFileSvaEnmix(temp_file, sepType = sepType)
    validatePhenotypeReplacementSvaEnmix(
        original = original,
        updated = written, SampleID = SampleID, pcColumns = pcColumns
    )

    if (!file.rename(phenoFile, backup_file)) {
        stop("Could not create a rollback copy of the phenotype file.",
            call. = FALSE
        )
    }
    backup_active <- TRUE
    if (!file.rename(temp_file, phenoFile)) {
        restored <- file.rename(backup_file, phenoFile)
        backup_active <- FALSE
        if (!isTRUE(restored)) {
            stop("Phenotype replacement failed and the original could not be restored from: ",
                backup_file,
                call. = FALSE
            )
        }
        stop("Phenotype replacement failed; the original file was restored.",
            call. = FALSE
        )
    }
    final_validation <- tryCatch(
        {
            final <- readPhenotypeFileSvaEnmix(phenoFile, sepType = sepType)
            validatePhenotypeReplacementSvaEnmix(
                original = original,
                updated = final, SampleID = SampleID, pcColumns = pcColumns
            )
            NULL
        },
        error = identity
    )
    if (inherits(final_validation, "error")) {
        unlink(phenoFile)
        restored <- file.rename(backup_file, phenoFile)
        backup_active <- FALSE
        if (!isTRUE(restored)) {
            stop("Final phenotype validation failed and the original could not be restored from: ",
                backup_file,
                call. = FALSE
            )
        }
        stop("Final phenotype validation failed; the original file was restored: ",
            conditionMessage(final_validation),
            call. = FALSE
        )
    }
    if (file.exists(backup_file) && unlink(backup_file) != 0L) {
        stop("The phenotype file was replaced, but its temporary backup could not be removed.",
            call. = FALSE
        )
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
#'     svaData = list(sva = ex$sva),
#'     mergedPheno = ex$mergedPheno,
#'     analysisData = ex$analysisData,
#'     phenoFile = pheno_file,
#'     SampleID = "Sample_Name",
#'     dataBaseDir = file.path(temp_dir, "data"),
#'     rBaseDir = file.path(temp_dir, "rData"),
#'     scriptLabel = "svaEnmixExample",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(output_paths)
#'
#' @export
writeSvaEnmixOutputs <- function(
    svaData, mergedPheno, analysisData = NULL,
    phenoFile = NULL, SampleID = "Sample_Name", sepType = NULL,
    dataBaseDir = "data", rBaseDir = "rData", scriptLabel = "svaEnmix",
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_writeSvaEnmixOutputs.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    data_dir <- file.path(dataBaseDir, scriptLabel)
    r_dir <- file.path(rBaseDir, scriptLabel)
    dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(r_dir, recursive = TRUE, showWarnings = FALSE)

    sva_rdata_path <- file.path(r_dir, "svaMatrix.RData")
    sva_csv_path <- file.path(data_dir, "svaMatrix.csv")
    pheno_output_path <- phenoFile

    saveNamedObjectMinfiEwasWater(svaData$sva, "sva", sva_rdata_path)
    utils::write.csv(svaData$sva, sva_csv_path, row.names = TRUE)

    if (!is.null(analysisData)) {
        for (i in seq_len(analysisData$K)) {
            utils::capture.output(summary(analysisData$fullModels[[i]]),
                file = file.path(data_dir, paste0(
                    "summary_full_sva",
                    i, ".txt"
                ))
            )

            if (length(analysisData$droptermSteps[[i]]) > 0L) {
                for (step in analysisData$droptermSteps[[i]]) {
                    utils::capture.output(step$dropterm,
                        file = file.path(
                            data_dir,
                            paste0("dropterm_step_sva", i, ".txt")
                        ),
                        append = TRUE
                    )
                    utils::capture.output(step$summary,
                        file = file.path(
                            data_dir,
                            paste0("dropterm_model_sva", i, ".txt")
                        ),
                        append = TRUE
                    )
                }
            }

            utils::capture.output(analysisData$anovaFull[[i]],
                file = file.path(data_dir, paste0(
                    "anova_full_sva",
                    i, ".txt"
                ))
            )
            utils::capture.output(analysisData$anovaReduced[[i]],
                file = file.path(data_dir, paste0(
                    "anova_reduced_sva",
                    i, ".txt"
                ))
            )
        }
    }

    if (!is.null(pheno_output_path)) {
        replacePhenotypeFileSvaEnmix(
            mergedPheno = mergedPheno,
            phenoFile = pheno_output_path, SampleID = SampleID,
            pcColumns = colnames(svaData$sva), sepType = sepType
        )
    }

    emitLogMinfiEwasWater(
        c(
            paste(
                "SVA Matrix RData saved to: ",
                sva_rdata_path
            ), paste(
                "SVA Matrix CSV saved to:   ",
                sva_csv_path
            ), paste("Updated phenotype + SVA:   ",
                if (is.null(pheno_output_path)) "not saved" else pheno_output_path),
            "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(
        svaRData = sva_rdata_path, svaCSV = sva_csv_path,
        phenoWithSva = pheno_output_path, dataDir = data_dir,
        rDir = r_dir
    ), class = "dnaEPICO_svaEnmix_paths")
}
