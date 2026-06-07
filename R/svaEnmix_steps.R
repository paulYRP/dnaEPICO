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
#' @return A list with class `"dnaEPICO_svaEnmix_sva"` containing the surrogate
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
estimateSvaEnmixControls <- function(
    RGSet,
    ctrlSvaPercVar = 0.90,
    ctrlSvaFlag = 1,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_estimateSvaEnmixControls.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  sva <- ENmix::ctrlsva(
    rgSet = RGSet,
    percvar = ctrlSvaPercVar,
    flag = ctrlSvaFlag
  )
  sva <- as.matrix(sva)

  emitLogMinfiEwasWater(
    c(
      paste("ctrlSva percvar:          ", ctrlSvaPercVar),
      paste("ctrlSva flag:             ", ctrlSvaFlag),
      paste("Number of surrogate variables:", ncol(sva)),
      "Surrogate variables matrix (first few rows):",
      previewLinesMinfiEwasWater(sva),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      sva = sva,
      K = ncol(sva),
      ctrlSvaPercVar = ctrlSvaPercVar,
      ctrlSvaFlag = ctrlSvaFlag
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
#'   targets = ex$targets,
#'   sva = ex$sva,
#'   SampleID = "Sample_Name",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' colnames(merged_pheno)[seq_len(4)]
#'
#' @export
mergeSvaTargetsEnmix <- function(
    targets,
    sva,
    SampleID = "Sample_Name",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_mergeSvaTargetsEnmix.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  if (!(SampleID %in% colnames(targets))) {
    stop("SampleID column not found in targets: ", SampleID, call. = FALSE)
  }

  if (is.null(rownames(sva))) {
    stop("sva must have row names that match the phenotype SampleID column.", call. = FALSE)
  }

  sva_data <- data.frame(
    sample_id = rownames(sva),
    sva,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(sva_data)[1] <- SampleID

  match_idx <- match(as.character(targets[[SampleID]]), as.character(sva_data[[SampleID]]))
  merged_pheno <- cbind(
    targets,
    sva_data[match_idx, setdiff(colnames(sva_data), SampleID), drop = FALSE]
  )

  emitLogMinfiEwasWater(
    c(
      paste("Merged phenotype rows:     ", nrow(merged_pheno)),
      paste("Merged phenotype columns:  ", ncol(merged_pheno)),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
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
#' @return A list with class `"dnaEPICO_svaEnmix_analysis"` containing the
#'   aligned Sentrix factors, full and reduced linear models, and ANOVA tables.
#'
#' @description
#' Fit linear models for each surrogate variable against Sentrix chip and
#' Sentrix position, perform backward elimination with `MASS::dropterm()`, and
#' return the in-memory analysis objects.
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
    sva,
    RGSet,
    SentrixIDColumn = "Sentrix_ID",
    SentrixPositionColumn = "Sentrix_Position",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_analyzeSvaEnmix.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  sample_names <- svaEnmixGetRGSetSampleNames(RGSet)

  if (length(sample_names) == 0L) {
    stop("Could not determine sample names from the loaded RGSet.", call. = FALSE)
  }

  match_idx <- match(sample_names, rownames(sva))

  if (anyNA(match_idx)) {
    stop("The row names of sva do not align with the RGSet sample names.", call. = FALSE)
  }

  sva <- sva[match_idx, , drop = FALSE]
  col_data <- svaEnmixGetRGSetColData(RGSet)
  sentrix_id <- as.factor(col_data[[SentrixIDColumn]])
  sentrix_position <- as.factor(col_data[[SentrixPositionColumn]])
  K <- ncol(sva)

  full_models <- lapply(
    seq_len(K),
    function(i) {
      stats::lm(
        sva[, i] ~ SentrixID + SentrixPosition,
        data = data.frame(
          SentrixID = sentrix_id,
          SentrixPosition = sentrix_position
        )
      )
    }
  )

  reduced_models <- vector("list", K)
  dropterm_steps <- vector("list", K)

  for (i in seq_len(K)) {
    model_tmp <- full_models[[i]]
    model_steps <- list()

    repeat {
      drop_table <- MASS::dropterm(model_tmp, test = "F")
      drop_p <- drop_table$`Pr(F)`
      valid_idx <- which(!is.na(drop_p))

      if (length(valid_idx) == 0L) {
        break
      }

      max_idx <- valid_idx[which.max(drop_p[valid_idx])]
      max_p <- drop_p[[max_idx]]
      term_to_drop <- rownames(drop_table)[[max_idx]]

      if (!is.finite(max_p) || max_p <= 0.05 || identical(term_to_drop, "<none>")) {
        break
      }

      model_tmp <- stats::update(model_tmp, paste(". ~ . -", term_to_drop))
      model_steps[[length(model_steps) + 1L]] <- list(
        dropterm = drop_table,
        summary = summary(model_tmp)
      )
    }

    reduced_models[[i]] <- model_tmp
    dropterm_steps[[i]] <- model_steps
  }

  anova_full <- lapply(full_models, stats::anova)
  anova_reduced <- lapply(reduced_models, stats::anova)

  emitLogMinfiEwasWater(
    c(
      paste("Number of surrogate variables (K):", K),
      paste("SentrixID class:            ", class(sentrix_id)[1]),
      paste("SentrixID unique levels:    ", length(unique(sentrix_id))),
      previewLinesMinfiEwasWater(table(sentrix_id)),
      paste("SentrixPosition class:      ", class(sentrix_position)[1]),
      paste("SentrixPosition unique levels:", length(unique(sentrix_position))),
      previewLinesMinfiEwasWater(table(sentrix_position)),
      "First row of SVA matrix:",
      previewLinesMinfiEwasWater(sva[1, , drop = FALSE]),
      paste(
        "Sample names in SVA matrix:",
        paste(rownames(sva)[seq_len(min(5L, nrow(sva)))], collapse = ", ")
      ),
      paste(
        "Sample names in colData(RGSet):",
        paste(sample_names[seq_len(min(5L, length(sample_names)))], collapse = ", ")
      ),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      sva = sva,
      K = K,
      sentrixID = sentrix_id,
      sentrixPosition = sentrix_position,
      fullModels = full_models,
      reducedModels = reduced_models,
      droptermSteps = dropterm_steps,
      anovaFull = anova_full,
      anovaReduced = anova_reduced
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
  page_limit <- as.integer(svaEnmixPositiveNumber(page_limit, 6))
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
        page = page,
        rowIndices = row_block,
        colIndices = col_block
      )
      page <- page + 1L
    }
  }

  pages
}

#' @keywords internal
#' @noRd
svaEnmixMatrixPlotDimensions <- function(width, height, res, panel_count) {
  fallback_width <- 2000
  fallback_height <- 1000
  width <- svaEnmixPositiveNumber(width, fallback_width)
  height <- svaEnmixPositiveNumber(height, fallback_height)
  res <- svaEnmixPositiveNumber(res, NA_real_)
  panel_count <- svaEnmixPositiveNumber(panel_count, NA_real_)

  if (
    !is.finite(res) || res <= 0 ||
      !is.finite(panel_count) || panel_count <= 0L
  ) {
    return(c(
      width = as.integer(ceiling(width)),
      height = as.integer(ceiling(height))
    ))
  }

  minimum_pixels <- ceiling(
    round(((panel_count * 1.6) + 1.0) * res, digits = 6L)
  )
  c(
    width = as.integer(max(width, minimum_pixels)),
    height = as.integer(max(height, minimum_pixels))
  )
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
svaEnmixMatrixPlotState <- function(analysisData, width, height, res) {
  sentrix_id_levels <- levels(analysisData$sentrixID)
  sentrix_id_legend_limit <- 20L
  page_panel_limit <- 6L
  pages <- svaEnmixMatrixPages(analysisData$K, page_panel_limit)
  page_panel_count <- min(analysisData$K, page_panel_limit)
  dimensions <- svaEnmixMatrixPlotDimensions(
    width = width,
    height = height,
    res = res,
    panel_count = page_panel_count
  )

  list(
    width = dimensions[["width"]],
    height = dimensions[["height"]],
    pages = pages,
    pageCount = length(pages),
    pagePanelLimit = page_panel_limit,
    sentrixIDLevels = sentrix_id_levels,
    sentrixIDLegendLimit = sentrix_id_legend_limit,
    showSentrixIDLegend = length(sentrix_id_levels) <= sentrix_id_legend_limit
  )
}

#' @keywords internal
#' @noRd
svaEnmixLogMatrixPlotState <- function(state, res, verbose, log_path) {
  log_lines <- c(
    paste(
      "SVA matrix plot dimensions:",
      state$width,
      "x",
      state$height,
      "@",
      res
    )
  )

  if (state$pageCount > 1L) {
    log_lines <- c(
      log_lines,
      paste(
        "SVA matrix plot pages:",
        state$pageCount,
        "pages with up to",
        state$pagePanelLimit,
        "x",
        state$pagePanelLimit,
        "panels per page"
      )
    )
  }

  if (!isTRUE(state$showSentrixIDLegend)) {
    log_lines <- c(
      log_lines,
      paste(
        "SentrixID legend suppressed for matrix plot:",
        length(state$sentrixIDLevels),
        "levels exceed limit of",
        state$sentrixIDLegendLimit
      )
    )
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
  graphics::plot(
    sva[, 1],
    sva[, 2],
    col = color_map[analysisData$sentrixID],
    pch = 16,
    xlab = "Surrogate Variable 1 (PC1)",
    ylab = "Surrogate Variable 2 (PC2)",
    main = "Surrogate Variables Colored by Chip (SentrixID)"
  )
  graphics::legend(
    "topright",
    legend = levels(analysisData$sentrixID),
    col = color_map,
    pch = 16,
    title = "SentrixID",
    cex = 0.6
  )
}

#' @keywords internal
#' @noRd
svaEnmixPlotSentrixPosition <- function(analysisData) {
  color_map <- grDevices::rainbow(length(levels(analysisData$sentrixPosition)))
  sva <- analysisData$sva
  graphics::plot(
    sva[, 1],
    sva[, 2],
    col = color_map[analysisData$sentrixPosition],
    pch = 16,
    xlab = "Surrogate Variable 1 (PC1)",
    ylab = "Surrogate Variable 2 (PC2)",
    main = "Surrogate Variables Colored by Sentrix Position"
  )
  graphics::legend(
    "topright",
    legend = levels(analysisData$sentrixPosition),
    col = color_map,
    pch = 16,
    title = "SentrixPosition",
    cex = 0.6
  )
}

#' @keywords internal
#' @noRd
svaEnmixMatrixPlotScales <- function(K) {
  list(
    axis = svaEnmixMatrixScale(
      K,
      multiplier = 3.4,
      minimum = 0.35,
      maximum = 0.7
    ),
    label = svaEnmixMatrixScale(
      K,
      multiplier = 3.8,
      minimum = 0.45,
      maximum = 0.8
    ),
    point = svaEnmixMatrixScale(
      K,
      multiplier = 4.0,
      minimum = 0.45,
      maximum = 0.8
    ),
    legend = svaEnmixMatrixScale(
      K,
      multiplier = 3.2,
      minimum = 0.35,
      maximum = 0.6
    )
  )
}

#' @keywords internal
#' @noRd
svaEnmixPlotMatrixLegends <- function(
    state,
    sentrix_position_levels,
    color_map,
    pch_map,
    legend_cex
) {
  if (isTRUE(state$showSentrixIDLegend)) {
    graphics::legend(
      "topright",
      legend = state$sentrixIDLevels,
      col = color_map,
      pch = 15,
      title = "SentrixID",
      cex = legend_cex,
      bty = "n"
    )
  }

  graphics::legend(
    if (isTRUE(state$showSentrixIDLegend)) "bottomright" else "topright",
    legend = sentrix_position_levels,
    pch = pch_map,
    title = "SentrixPosition",
    cex = legend_cex,
    bty = "n"
  )
}

#' @keywords internal
#' @noRd
svaEnmixPlotMatrixPageTitle <- function(page, page_count) {
  title <- "Effects of Sentrix ID (color) and Sentrix Position (shape)"
  if (page_count <= 1L) {
    return(title)
  }

  paste(
    title,
    sprintf(
      "(page %d of %d: SV%d-SV%d by SV%d-SV%d)",
      page$page,
      page_count,
      min(page$rowIndices),
      max(page$rowIndices),
      min(page$colIndices),
      max(page$colIndices)
    )
  )
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
    family = "Times",
    las = 1,
    mar = svaEnmixMatrixPanelMargin(page_panel_count),
    oma = c(3.2, 3.2, 2.0, 0.5),
    mgp = c(1.1, 0.35, 0),
    tcl = -0.18,
    cex.axis = scales$axis,
    cex.lab = scales$label
  )

  for (row_pos in seq_along(row_indices)) {
    for (col_pos in seq_along(col_indices)) {
      row_index <- row_indices[[row_pos]]
      col_index <- col_indices[[col_pos]]
      graphics::plot(
        sva[, col_index],
        sva[, row_index],
        col = color_map[analysisData$sentrixID],
        pch = pch_map[analysisData$sentrixPosition],
        cex = scales$point,
        xlab = if (row_pos == length(row_indices)) {
          paste("SV", col_index)
        } else {
          ""
        },
        ylab = if (col_pos == 1L) paste("SV", row_index) else "",
        main = "",
        xaxt = if (row_pos == length(row_indices)) "s" else "n",
        yaxt = if (col_pos == 1L) "s" else "n"
      )
      if (row_pos == 1L && col_pos == 1L) {
        svaEnmixPlotMatrixLegends(
          state, sentrix_position_levels, color_map, pch_map, scales$legend
        )
      }
    }
  }

  graphics::mtext(
    svaEnmixPlotMatrixPageTitle(page, state$pageCount),
    side = 3, outer = TRUE, line = 0.5, cex = 0.9
  )
}

#' @keywords internal
#' @noRd
svaEnmixSaveMatrixPage <- function(file, draw_fun, width, height, res) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  grDevices::tiff(
    filename = file,
    width = width,
    height = height,
    res = res,
    type = "cairo"
  )
  tryCatch(
    draw_fun(),
    finally = grDevices::dev.off()
  )

  invisible(file)
}

#' @keywords internal
#' @noRd
svaEnmixRunMatrixPlot <- function(analysisData, state, display, file, res) {
  if (is.null(file) && !isTRUE(display)) {
    return(invisible(NULL))
  }

  output_files <- NULL
  if (!is.null(file)) {
    output_files <- vapply(
      state$pages,
      function(page) svaEnmixMatrixPageFile(file, page$page, state$pageCount),
      character(1)
    )
    output_files <- unname(output_files)
  }

  if (!is.null(file)) {
    for (page in state$pages) {
      page_file <- svaEnmixMatrixPageFile(file, page$page, state$pageCount)
      svaEnmixSaveMatrixPage(
        file = page_file,
        draw_fun = function() svaEnmixPlotMatrixPage(analysisData, state, page),
        width = state$width,
        height = state$height,
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
#' @param plot Character. Plot type: `"sentrix_id"`, `"sentrix_position"`, or
#'   `"matrix"`.
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
    plot = c("sentrix_id", "sentrix_position", "matrix"),
    display = FALSE,
    file = NULL,
    width = 2000L,
    height = 1000L,
    res = 150L,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_plotSvaEnmix.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )
  plot <- match.arg(plot)
  sva <- analysisData$sva
  matrix_state <- list(width = width, height = height)
  will_draw <- !is.null(file) || isTRUE(display)

  if (plot %in% c("sentrix_id", "sentrix_position") && ncol(sva) < 2L) {
    stop(
      "At least two surrogate variables are required for this plot.",
      call. = FALSE
    )
  }

  if (identical(plot, "matrix")) {
    matrix_state <- svaEnmixMatrixPlotState(analysisData, width, height, res)
    if (will_draw) {
      svaEnmixLogMatrixPlotState(matrix_state, res, verbose, log_path)
    }
  }

  draw_fun <- switch(
    plot,
    sentrix_id = function() svaEnmixPlotSentrixID(analysisData),
    sentrix_position = function() svaEnmixPlotSentrixPosition(analysisData),
    matrix = NULL
  )

  if (identical(plot, "matrix")) {
    output <- svaEnmixRunMatrixPlot(
      analysisData = analysisData,
      state = matrix_state,
      display = display,
      file = file,
      res = res
    )
  } else {
    output <- runPlotMinfiEwasWater(
      draw_fun = draw_fun,
      display = display,
      file = file,
      width = matrix_state$width,
      height = matrix_state$height,
      res = res
    )
  }

  emitLogMinfiEwasWater(
    c(
      paste("SVA plot type:            ", plot),
      paste("Display plot:             ", display),
      svaEnmixSavedPlotLogLines(output),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  invisible(output)
}

#' Write svaEnmix outputs to disk
#'
#' @param svaData Object returned by `estimateSvaEnmixControls()`.
#' @param mergedPheno Phenotype data frame returned by `mergeSvaTargetsEnmix()`.
#' @param analysisData Optional object returned by `analyzeSvaEnmix()`.
#' @param phenoFile Character or `NULL`. When supplied, `mergedPheno` is written
#'   back to this path for legacy compatibility.
#' @param dataBaseDir Character. Base directory used for saved data outputs.
#' @param rBaseDir Character. Base directory used for saved `.RData` outputs.
#' @param scriptLabel Character. Label used to create the output subdirectory.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_svaEnmix_paths"` containing the paths
#'   written to disk.
#'
#' @description
#' Write the legacy CSV, `.RData`, and text-summary outputs used by the original
#' `svaEnmix()` workflow.
#'
#' @examples
#' ex <- dnaEPICO:::exampleSvaAnalysisStateDnaEpico()
#' temp_dir <- tempdir()
#' output_paths <- writeSvaEnmixOutputs(
#'   svaData = list(sva = ex$sva),
#'   mergedPheno = ex$mergedPheno,
#'   analysisData = ex$analysisData,
#'   phenoFile = file.path(temp_dir, "phenoLC.csv"),
#'   dataBaseDir = file.path(temp_dir, "data"),
#'   rBaseDir = file.path(temp_dir, "rData"),
#'   scriptLabel = "svaEnmixExample",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(output_paths)
#'
#' @export
writeSvaEnmixOutputs <- function(
    svaData,
    mergedPheno,
    analysisData = NULL,
    phenoFile = NULL,
    dataBaseDir = "data",
    rBaseDir = "rData",
    scriptLabel = "svaEnmix",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_writeSvaEnmixOutputs.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
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

  if (!is.null(pheno_output_path)) {
    utils::write.csv(mergedPheno, file = pheno_output_path, row.names = FALSE)
  }

  if (!is.null(analysisData)) {
    for (i in seq_len(analysisData$K)) {
      utils::capture.output(
        summary(analysisData$fullModels[[i]]),
        file = file.path(data_dir, paste0("summary_full_sva", i, ".txt"))
      )

      if (length(analysisData$droptermSteps[[i]]) > 0L) {
        for (step in analysisData$droptermSteps[[i]]) {
          utils::capture.output(
            step$dropterm,
            file = file.path(data_dir, paste0("dropterm_step_sva", i, ".txt")),
            append = TRUE
          )
          utils::capture.output(
            step$summary,
            file = file.path(data_dir, paste0("dropterm_model_sva", i, ".txt")),
            append = TRUE
          )
        }
      }

      utils::capture.output(
        analysisData$anovaFull[[i]],
        file = file.path(data_dir, paste0("anova_full_sva", i, ".txt"))
      )
      utils::capture.output(
        analysisData$anovaReduced[[i]],
        file = file.path(data_dir, paste0("anova_reduced_sva", i, ".txt"))
      )
    }
  }

  emitLogMinfiEwasWater(
    c(
      paste("SVA Matrix RData saved to: ", sva_rdata_path),
      paste("SVA Matrix CSV saved to:   ", sva_csv_path),
      paste(
        "Saved phenoLC + SVA:       ",
        if (is.null(pheno_output_path)) "not saved" else pheno_output_path
      ),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      svaRData = sva_rdata_path,
      svaCSV = sva_csv_path,
      phenoWithSva = pheno_output_path,
      dataDir = data_dir,
      rDir = r_dir
    ),
    class = "dnaEPICO_svaEnmix_paths"
  )
}
