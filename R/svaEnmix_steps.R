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
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
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
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
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

  sample_names <- Biobase::sampleNames(RGSet)
  match_idx <- match(sample_names, rownames(sva))

  if (anyNA(match_idx)) {
    stop("The row names of sva do not align with the RGSet sample names.", call. = FALSE)
  }

  sva <- sva[match_idx, , drop = FALSE]
  col_data <- SummarizedExperiment::colData(RGSet)
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
#' @return Invisibly returns `file` when a TIFF is written, otherwise `NULL`.
#'
#' @description
#' Draw one of the standard surrogate-variable plots used by `svaEnmix()`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
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

  if (plot %in% c("sentrix_id", "sentrix_position") && ncol(sva) < 2L) {
    stop("At least two surrogate variables are required for this plot.", call. = FALSE)
  }

  draw_fun <- switch(
    plot,
    sentrix_id = function() {
      color_map <- grDevices::rainbow(length(levels(analysisData$sentrixID)))
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
    },
    sentrix_position = function() {
      color_map <- grDevices::rainbow(length(levels(analysisData$sentrixPosition)))
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
    },
    matrix = function() {
      K <- analysisData$K
      color_map <- grDevices::rainbow(length(levels(analysisData$sentrixID)))
      pch_map <- seq_along(levels(analysisData$sentrixPosition))

      old_par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old_par), add = TRUE)
      graphics::par(mfrow = c(K, K), family = "Times", las = 1)

      for (i in seq_len(K)) {
        for (j in seq_len(K)) {
          graphics::plot(
            sva[, j],
            sva[, i],
            col = color_map[analysisData$sentrixID],
            pch = pch_map[analysisData$sentrixPosition],
            xlab = paste("SV", j),
            ylab = paste("SV", i),
            main = "Effects of Sentrix ID (color) & Sentrix Position (shape)"
          )

          if (i == 1L && j == 1L) {
            graphics::legend(
              "topright",
              legend = levels(analysisData$sentrixID),
              col = color_map,
              pch = 15,
              title = "SentrixID",
              cex = 0.6
            )
            graphics::legend(
              "bottomright",
              legend = levels(analysisData$sentrixPosition),
              pch = pch_map,
              title = "SentrixPosition",
              cex = 0.6
            )
          }
        }
      }
    }
  )

  output <- runPlotMinfiEwasWater(
    draw_fun = draw_fun,
    display = display,
    file = file,
    width = width,
    height = height,
    res = res
  )

  emitLogMinfiEwasWater(
    c(
      paste("SVA plot type:            ", plot),
      paste("Display plot:             ", display),
      paste(
        "Saved plot path:          ",
        if (is.null(output)) "not saved" else output
      ),
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
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
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
