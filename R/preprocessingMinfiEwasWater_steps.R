#' Read IDAT files into an annotated RGChannelSet
#'
#' Read methylation-array IDAT files with `minfi::read.metharray.exp()`, set
#' sample names from the phenotype table, apply the requested annotation, and
#' return the resulting `RGChannelSet`.
#'
#' @param idatFolder Character. Directory containing the IDAT files.
#' @param targets Data frame returned by `readPhenotypeTargets()` or an
#'   equivalent phenotype table.
#' @param SampleID Character. Name of the phenotype column containing sample
#'   identifiers used to label the `RGChannelSet`.
#' @param arrayType Character. Array name passed to
#'   `Biobase::annotation(RGSet)`, for example
#'   `"IlluminaHumanMethylationEPICv2"`.
#' @param annotationVersion Character. Annotation build passed to
#'   `Biobase::annotation(RGSet)`, for example `"20a1.hg38"` for EPIC v2 hg38
#'   annotations or `"ilmn12.hg19"` for 450K hg19 annotations.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return An annotated `RGChannelSet`.
#'
#' @examples
#' if (requireNamespace("minfiData", quietly = TRUE) &&
#'     requireNamespace("IlluminaHumanMethylation450kmanifest", quietly = TRUE) &&
#'     requireNamespace("IlluminaHumanMethylation450kanno.ilmn12.hg19", quietly = TRUE)) {
#'   ex <- dnaEPICO:::exampleMinfiIdatInputsDnaEpico(n = 4)
#'   rgset <- readRGSetMinfiEwasWater(
#'     idatFolder = ex$idatFolder,
#'     targets = ex$targets,
#'     SampleID = "Sample_Name",
#'     arrayType = ex$arrayType,
#'     annotationVersion = ex$annotationVersion
#'   )
#'   class(rgset)
#' }
#'
#' @export
readRGSetMinfiEwasWater <- function(
    idatFolder,
    targets,
    SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
    annotationVersion = "20a1.hg38",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_readRGSetMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  if (!dir.exists(idatFolder)) {
    stop("idatFolder does not exist: ", idatFolder, call. = FALSE)
  }

  if (!(SampleID %in% colnames(targets))) {
    stop("SampleID column not found in targets: ", SampleID, call. = FALSE)
  }

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      paste("IDAT folder:              ", idatFolder),
      paste("Array type:               ", arrayType),
      paste("Annotation version:       ", annotationVersion)
    ),
    verbose = verbose,
    log_path = log_path
  )

  RGSet <- minfi::read.metharray.exp(
    base = idatFolder,
    targets = targets,
    extended = FALSE,
    recursive = FALSE,
    verbose = FALSE
  )

  colnames(RGSet) <- targets[[SampleID]]
  Biobase::annotation(RGSet) <- c(
    array = arrayType,
    annotation = annotationVersion
  )

  manifest_lines <- utils::capture.output(methods::show(minfi::getManifest(RGSet)))

  emitLogMinfiEwasWater(
    c(
      paste(
        "RGSet loaded with",
        ncol(RGSet),
        "samples."
      ),
      paste(
        "Applied annotation:       ",
        paste(Biobase::annotation(RGSet), collapse = ", ")
      ),
      "Manifest used:",
      manifest_lines,
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  RGSet
}

#' Build raw minfi preprocessing objects
#'
#' Create a raw `MethylSet`, `RatioSet`, and genome-mapped object from an
#' `RGChannelSet`, and return them together in a single structured object.
#'
#' @param RGSet An `RGChannelSet`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_minfiEwasWater_raw"` containing `MSet`,
#'   `RatioSet`, and `GSet`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' raw_data <- buildRawMinfiEwasWater(
#'   RGSet = ex$RGSet,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(raw_data)
#'
#' @export
buildRawMinfiEwasWater <- function(
    RGSet,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_buildRawMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  emitLogMinfiEwasWater(
    "Running preprocessRaw(), ratioConvert(), and mapToGenome().",
    verbose = verbose,
    log_path = log_path
  )

  MSet <- minfi::preprocessRaw(RGSet)
  RatioSet <- minfi::ratioConvert(MSet, what = "both", keepCN = TRUE)
  GSet <- minfi::mapToGenome(RatioSet)

  meth_cols <- seq_len(min(ncol(MSet), 3L))
  gset_cols <- seq_len(min(ncol(GSet), 5L))

  emitLogMinfiEwasWater(
    c(
      paste("MSet dimensions:          ", paste(dim(MSet), collapse = " x ")),
      paste(
        "RatioSet dimensions:      ",
        paste(dim(RatioSet), collapse = " x ")
      ),
      paste("GSet dimensions:          ", paste(dim(GSet), collapse = " x ")),
      "Preview of methylated intensities:",
      previewLinesMinfiEwasWater(minfi::getMeth(MSet)[, meth_cols, drop = FALSE]),
      "Preview of unmethylated intensities:",
      previewLinesMinfiEwasWater(minfi::getUnmeth(MSet)[, meth_cols, drop = FALSE]),
      "Preview of beta values:",
      previewLinesMinfiEwasWater(minfi::getBeta(GSet)[, gset_cols, drop = FALSE]),
      "Preview of M-values:",
      previewLinesMinfiEwasWater(minfi::getM(GSet)[, gset_cols, drop = FALSE]),
      "Preview of copy-number values:",
      previewLinesMinfiEwasWater(minfi::getCN(GSet)[, gset_cols, drop = FALSE]),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      MSet = MSet,
      RatioSet = RatioSet,
      GSet = GSet
    ),
    class = "dnaEPICO_minfiEwasWater_raw"
  )
}

#' Assess sample quality before sample filtering
#'
#' Compute minfi QC metrics and detection P values, identify failed samples
#' using the requested threshold, and return the assessment as a single object.
#'
#' @param rawData Object returned by `buildRawMinfiEwasWater()`.
#' @param RGSet An `RGChannelSet` aligned with `rawData`.
#' @param qcCutoff Numeric. Cutoff passed to `minfi::plotQC()` when the QC plot
#'   is drawn.
#' @param detPtype Character. Detection P-value mode passed to
#'   `minfi::detectionP()`. Common values in minfi workflows are `"m+u"` and
#'   `"negative"`. The default used here is `"m+u"`.
#' @param detPThreshold Numeric. Samples with mean detection P value above this
#'   threshold are marked as failed.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_minfiEwasWater_assessment"` containing
#'   the QC object, detection P matrix, mean detection P values, and failed
#'   sample identifiers.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' raw_data <- buildRawMinfiEwasWater(ex$RGSet, verbose = FALSE, logs = FALSE)
#' assessment <- assessSamplesMinfiEwasWater(
#'   rawData = raw_data,
#'   RGSet = ex$RGSet,
#'   detPThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(assessment)
#'
#' @export
assessSamplesMinfiEwasWater <- function(
    rawData,
    RGSet,
    qcCutoff = 10.5,
    detPtype = "m+u",
    detPThreshold = 0.05,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_assessSamplesMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  qc <- minfi::getQC(rawData$MSet)
  detP <- minfi::detectionP(RGSet, type = detPtype)
  meanDetP <- colMeans(detP)
  failedSamples <- names(meanDetP[meanDetP > detPThreshold])

  preview_cols <- seq_len(min(ncol(detP), 5L))

  emitLogMinfiEwasWater(
    c(
      paste("QC cutoff (median):       ", qcCutoff),
      paste("Detection P type:         ", detPtype),
      paste("Detection P threshold:    ", detPThreshold),
      "Preview of detection P values:",
      previewLinesMinfiEwasWater(detP[, preview_cols, drop = FALSE]),
      paste("Number of failed samples: ", length(failedSamples)),
      if (length(failedSamples) > 0L) {
        paste("Failed sample IDs:        ", paste(failedSamples, collapse = ", "))
      } else {
        "Failed sample IDs:         none"
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      qc = qc,
      qcCutoff = qcCutoff,
      detP = detP,
      detPtype = detPtype,
      detPThreshold = detPThreshold,
      meanDetP = meanDetP,
      failedSamples = failedSamples
    ),
    class = "dnaEPICO_minfiEwasWater_assessment"
  )
}

#' Plot quality-assessment outputs for preprocessingMinfiEwasWater
#'
#' Draw either the minfi QC plot or the detection P-value plot from an
#' assessment object returned by `assessSamplesMinfiEwasWater()`.
#'
#' @param assessment Object returned by `assessSamplesMinfiEwasWater()`.
#' @param plot Character. Plot type: `"qc"` or `"detection"`.
#' @param display Logical. If `TRUE`, draw the plot on the active graphics
#'   device.
#' @param file Character or `NULL`. TIFF file written when supplied.
#' @param width Integer. TIFF width in pixels when `file` is supplied.
#' @param height Integer. TIFF height in pixels when `file` is supplied.
#' @param res Integer. TIFF resolution in DPI when `file` is supplied.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return Invisibly returns the saved TIFF path when `file` is supplied,
#'   otherwise `NULL`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' raw_data <- buildRawMinfiEwasWater(ex$RGSet, verbose = FALSE, logs = FALSE)
#' assessment <- assessSamplesMinfiEwasWater(
#'   rawData = raw_data,
#'   RGSet = ex$RGSet,
#'   detPThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' plotAssessmentMinfiEwasWater(
#'   assessment = assessment,
#'   plot = "qc",
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#'
#' @export
plotAssessmentMinfiEwasWater <- function(
    assessment,
    plot = c("qc", "detection"),
    display = FALSE,
    file = NULL,
    width = 2000L,
    height = 1000L,
    res = 150L,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_plotAssessmentMinfiEwasWater.txt"
) {
  plot <- match.arg(plot)

  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  draw_fun <- switch(
    plot,
    qc = function() {
      minfi::plotQC(assessment$qc, badSampleCutoff = assessment$qcCutoff)
    },
    detection = function() {
      graphics::barplot(
        assessment$meanDetP,
        las = 3,
        cex.names = 0.8,
        ylab = "Mean detection p-values"
      )
      graphics::abline(
        h = assessment$detPThreshold,
        col = "red",
        lwd = 2,
        lty = 2
      )
    }
  )

  runPlotMinfiEwasWater(
    draw_fun = draw_fun,
    display = display,
    file = file,
    width = width,
    height = height,
    res = res
  )

  emitLogMinfiEwasWater(
    c(
      paste("Assessment plot created:   ", plot),
      if (is.null(file)) {
        "Assessment plot file:      none"
      } else {
        paste("Assessment plot file:      ", file)
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  invisible(file)
}

#' Filter failed samples from an RGSet and phenotype table
#'
#' Remove failed samples identified during sample assessment and return the
#' filtered `RGChannelSet` together with the aligned phenotype table.
#'
#' @param RGSet An `RGChannelSet`.
#' @param targets Data frame containing phenotype information.
#' @param failedSamples Character vector of sample identifiers to remove.
#' @param SampleID Character. Name of the sample identifier column in `targets`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_minfiEwasWater_samples"` containing the
#'   filtered `RGSet`, aligned phenotype table, and failed sample identifiers.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' filtered_samples <- filterSamplesMinfiEwasWater(
#'   RGSet = ex$RGSet,
#'   targets = ex$targets,
#'   failedSamples = ex$targets$Sample_Name[1],
#'   SampleID = "Sample_Name",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' nrow(filtered_samples$targets)
#'
#' @export
filterSamplesMinfiEwasWater <- function(
    RGSet,
    targets,
    failedSamples = character(0),
    SampleID = "Sample_Name",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_filterSamplesMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  keep_samples <- !(colnames(RGSet) %in% failedSamples)
  RGSet_filtered <- RGSet[, keep_samples]

  sample_names <- colnames(RGSet_filtered)
  matched <- match(sample_names, targets[[SampleID]])
  targets_filtered <- targets[matched, , drop = FALSE]
  rownames(targets_filtered) <- NULL

  emitLogMinfiEwasWater(
    c(
      paste("Samples before filtering: ", ncol(RGSet)),
      paste("Samples after filtering:  ", ncol(RGSet_filtered)),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      RGSet = RGSet_filtered,
      targets = targets_filtered,
      failedSamples = failedSamples
    ),
    class = "dnaEPICO_minfiEwasWater_samples"
  )
}

#' Predict biological sex from a filtered raw-data object
#'
#' Predict sample sex from a genome-mapped methylation object, align the
#' predictions with phenotype data, and return a structured object that can be
#' plotted or merged into downstream phenotype tables.
#'
#' @param rawData Object returned by `buildRawMinfiEwasWater()`.
#' @param targets Filtered phenotype data frame aligned with `rawData`.
#' @param SampleID Character. Name of the sample identifier column in `targets`.
#' @param sexColumn Character. Name of the phenotype column containing reported
#'   sex.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_minfiEwasWater_sex"` containing the sex
#'   prediction result, aligned phenotype data, plotting data, and mismatch
#'   table.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiWorkflowStateDnaEpico()
#' sex_data <- predictSexMinfiEwasWater(
#'   rawData = ex$rawFiltered,
#'   targets = ex$sampleData$targets,
#'   SampleID = "Sample_Name",
#'   sexColumn = "Sex",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(sex_data)
#'
#' @export
predictSexMinfiEwasWater <- function(
    rawData,
    targets,
    SampleID = "Sample_Name",
    sexColumn = "Sex",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_predictSexMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  if (!(SampleID %in% colnames(targets))) {
    stop("SampleID column not found in targets: ", SampleID, call. = FALSE)
  }

  if (!(sexColumn %in% colnames(targets))) {
    stop("sexColumn not found in targets: ", sexColumn, call. = FALSE)
  }

  pSex <- minfi::getSex(rawData$GSet)
  sex_plot_data <- as.data.frame(pSex, stringsAsFactors = FALSE)
  sex_plot_data$SampleID <- rownames(sex_plot_data)

  matched <- match(sex_plot_data$SampleID, targets[[SampleID]])
  targets_aligned <- targets[matched, , drop = FALSE]
  rownames(targets_aligned) <- NULL

  recoded_sex <- targets_aligned[[sexColumn]]
  if (is.character(recoded_sex) || is.factor(recoded_sex)) {
    recoded_sex <- as.character(recoded_sex)
    female_values <- c("F", "Female", "f", "female", "FEMALE")
    male_values <- c("M", "Male", "m", "male", "MALE")
    recoded <- rep(NA_integer_, length(recoded_sex))
    recoded[recoded_sex %in% female_values] <- 0L
    recoded[recoded_sex %in% male_values] <- 1L
    recoded_sex <- recoded
  } else {
    recoded_sex <- as.integer(recoded_sex)
  }

  targets_aligned[[sexColumn]] <- recoded_sex
  targets_aligned$PredSex <- ifelse(pSex$predictedSex == "F", 0L, 1L)

  sex_plot_data[[sexColumn]] <- targets_aligned[[sexColumn]]
  mismatches <- targets_aligned[
    !is.na(targets_aligned[[sexColumn]]) &
      !is.na(targets_aligned$PredSex) &
      targets_aligned[[sexColumn]] != targets_aligned$PredSex,
    ,
    drop = FALSE
  ]

  sex_values <- targets[[sexColumn]]
  nSexNA <- sum(is.na(sex_values))
  if (is.character(sex_values) || is.factor(sex_values)) {
    known_sex <- c(
      "F", "Female", "f", "female", "FEMALE",
      "M", "Male", "m", "male", "MALE"
    )
    unknown_sex <- setdiff(unique(as.character(sex_values)), known_sex)
  } else {
    unknown_sex <- character(0)
  }

  emitLogMinfiEwasWater(
    c(
      "Sex column integrity check:",
      paste("  Sex column used:        ", sexColumn),
      paste("  NA values:              ", nSexNA),
      paste(
        "  Unknown labels:         ",
        if (length(unknown_sex) == 0L) {
          "none"
        } else {
          paste(unknown_sex, collapse = ", ")
        }
      ),
      paste("Mismatches found:         ", nrow(mismatches)),
      if (nrow(mismatches) > 0L) {
        previewLinesMinfiEwasWater(
          mismatches[, seq_len(min(3L, ncol(mismatches))), drop = FALSE]
        )
      } else {
        character(0)
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      pSex = pSex,
      targets = targets_aligned,
      sexPlotData = sex_plot_data,
      mismatches = mismatches,
      SampleID = SampleID,
      sexColumn = sexColumn
    ),
    class = "dnaEPICO_minfiEwasWater_sex"
  )
}

#' Plot predicted or clinical sex from `predictSexMinfiEwasWater()`
#'
#' @param sexData Object returned by `predictSexMinfiEwasWater()`.
#' @param type Character. Plot type: `"predicted"` for methylation-predicted
#'   sex or `"clinical"` for reported sex.
#' @param display Logical. If `TRUE`, draw the plot on the active graphics
#'   device.
#' @param file Character or `NULL`. TIFF file written when supplied.
#' @param width Integer. TIFF width in pixels when `file` is supplied.
#' @param height Integer. TIFF height in pixels when `file` is supplied.
#' @param res Integer. TIFF resolution in DPI when `file` is supplied.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return Invisibly returns the saved TIFF path when `file` is supplied,
#'   otherwise `NULL`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiWorkflowStateDnaEpico()
#' plotSexMinfiEwasWater(
#'   sexData = ex$sexData,
#'   type = "predicted",
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#'
#' @export
plotSexMinfiEwasWater <- function(
    sexData,
    type = c("predicted", "clinical"),
    display = FALSE,
    file = NULL,
    width = 2000L,
    height = 1000L,
    res = 70L,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_plotSexMinfiEwasWater.txt"
) {
  type <- match.arg(type)

  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  draw_fun <- switch(
    type,
    predicted = function() {
      graphics::plot(
        x = sexData$pSex$xMed,
        y = sexData$pSex$yMed,
        type = "n",
        xlab = "X chr, median total intensity (log2)",
        ylab = "Y chr, median total intensity (log2)"
      )
      graphics::text(
        x = sexData$pSex$xMed,
        y = sexData$pSex$yMed,
        labels = sexData$targets[[sexData$SampleID]],
        col = ifelse(
          sexData$pSex$predictedSex == "M",
          "deepskyblue",
          "deeppink3"
        )
      )
      graphics::legend(
        "bottomleft",
        c("M", "F"),
        col = c("deepskyblue", "deeppink3"),
        pch = 16
      )
    },
    clinical = function() {
      graphics::plot(
        x = sexData$sexPlotData$xMed,
        y = sexData$sexPlotData$yMed,
        type = "n",
        xlab = "X chr, median total intensity (log2)",
        ylab = "Y chr, median total intensity (log2)"
      )
      graphics::text(
        x = sexData$sexPlotData$xMed,
        y = sexData$sexPlotData$yMed,
        labels = sexData$sexPlotData$SampleID,
        col = ifelse(
          sexData$sexPlotData[[sexData$sexColumn]] == 1L,
          "deepskyblue",
          "deeppink3"
        )
      )
      graphics::legend(
        "bottomleft",
        c("M", "F"),
        col = c("deepskyblue", "deeppink3"),
        pch = 16
      )
    }
  )

  runPlotMinfiEwasWater(
    draw_fun = draw_fun,
    display = display,
    file = file,
    width = width,
    height = height,
    res = res
  )

  emitLogMinfiEwasWater(
    c(
      paste("Sex plot created:         ", type),
      if (is.null(file)) {
        "Sex plot file:            none"
      } else {
        paste("Sex plot file:            ", file)
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  invisible(file)
}

#' Normalize filtered samples with minfi and wateRmelon methods
#'
#' Apply one or more supported normalization methods to a filtered `RGSet` and
#' return all normalized objects together in a single result object.
#'
#' @param sampleData Object returned by `filterSamplesMinfiEwasWater()`.
#' @param sexColumn Character. Name of the phenotype column used as the optional
#'   sex covariate for normalization methods that support it.
#' @param normMethods Character vector or semicolon-separated string of
#'   normalization methods. Supported values are `"adjustedfunnorm"`,
#'   `"funnorm"`, `"illumina"`, `"quantile"`, and `"swan"`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_minfiEwasWater_norm"` containing the
#'   requested normalized objects and the first method as `primary`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' sample_data <- filterSamplesMinfiEwasWater(
#'   RGSet = ex$RGSet,
#'   targets = ex$targets,
#'   failedSamples = character(0),
#'   SampleID = "Sample_Name",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' norm_data <- normalizeMinfiEwasWater(
#'   sampleData = sample_data,
#'   sexColumn = "Sex",
#'   normMethods = "quantile",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(norm_data$normalized)
#'
#' @export
normalizeMinfiEwasWater <- function(
    sampleData,
    sexColumn = "Sex",
    normMethods = "adjustedfunnorm",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_normalizeMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  method_list <- splitOptionMinfiEwasWater(normMethods, sep = ";")
  normalized <- vector("list", length(method_list))
  names(normalized) <- method_list

  col_data <- SummarizedExperiment::colData(sampleData$RGSet)
  sex_vec <- NULL
  if (!is.null(sexColumn) && sexColumn %in% colnames(col_data)) {
    sex_vec <- col_data[[sexColumn]]
    sex_vec <- as.character(sex_vec)
    sex_vec[sex_vec %in% c("0", "F", "Female", "female", "FEMALE")] <- "F"
    sex_vec[sex_vec %in% c("1", "M", "Male", "male", "MALE")] <- "M"
    if (!all(sex_vec %in% c("F", "M"))) {
      sex_vec <- NULL
    }
  }

  emitLogMinfiEwasWater(
    c(
      paste(
        "Normalization methods:    ",
        paste(method_list, collapse = ", ")
      ),
      paste("Sex column:               ", sexColumn)
    ),
    verbose = verbose,
    log_path = log_path
  )

  for (i in seq_along(method_list)) {
    method <- method_list[[i]]
    emitLogMinfiEwasWater(
      paste("Applying normalization:   ", method),
      verbose = verbose,
      log_path = log_path
    )

    normalized[[i]] <- switch(
      method,
      adjustedfunnorm = wateRmelon::adjustedFunnorm(
        sampleData$RGSet,
        sex = sex_vec
      ),
      funnorm = minfi::preprocessFunnorm(sampleData$RGSet, sex = sex_vec),
      illumina = minfi::preprocessIllumina(sampleData$RGSet),
      quantile = minfi::preprocessQuantile(sampleData$RGSet, sex = sex_vec),
      swan = minfi::preprocessSWAN(sampleData$RGSet),
      stop("Unknown normalization method: ", method, call. = FALSE)
    )
  }

  emitLogMinfiEwasWater(
    "=======================================================================",
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      primary = normalized[[1L]],
      normalized = normalized,
      methods = method_list,
      sexColumn = sexColumn
    ),
    class = "dnaEPICO_minfiEwasWater_norm"
  )
}

#' Plot raw and normalized methylation distributions
#'
#' Draw the density comparison plot used to inspect raw versus normalized data.
#'
#' @param RGSet An `RGChannelSet` aligned with `targets`.
#' @param normData Object returned by `normalizeMinfiEwasWater()`.
#' @param targets Filtered phenotype data aligned with `RGSet`.
#' @param sexColumn Character. Name of the phenotype column used to colour the
#'   density curves.
#' @param display Logical. If `TRUE`, draw the plot on the active graphics
#'   device.
#' @param file Character or `NULL`. TIFF file written when supplied.
#' @param width Integer. TIFF width in pixels when `file` is supplied.
#' @param height Integer. TIFF height in pixels when `file` is supplied.
#' @param res Integer. TIFF resolution in DPI when `file` is supplied.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return Invisibly returns the saved TIFF path when `file` is supplied,
#'   otherwise `NULL`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiWorkflowStateDnaEpico()
#' plotNormalizationMinfiEwasWater(
#'   RGSet = ex$sampleData$RGSet,
#'   normData = ex$normData,
#'   targets = ex$sampleData$targets,
#'   sexColumn = "Sex",
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#'
#' @export
plotNormalizationMinfiEwasWater <- function(
    RGSet,
    normData,
    targets,
    sexColumn = "Sex",
    display = FALSE,
    file = NULL,
    width = 2000L,
    height = 1000L,
    res = 150L,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_plotNormalizationMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  draw_fun <- function() {
    graphics::par(mfrow = c(1, 2))
    minfi::densityPlot(
      RGSet,
      sampGroups = targets[[sexColumn]],
      main = "Raw",
      legend = FALSE
    )
    graphics::legend(
      "top",
      legend = levels(factor(targets[[sexColumn]])),
      text.col = RColorBrewer::brewer.pal(8, "Dark2")
    )

    minfi::densityPlot(
      minfi::getBeta(normData$primary),
      sampGroups = targets[[sexColumn]],
      main = "Normalized",
      legend = FALSE
    )
    graphics::legend(
      "top",
      legend = levels(factor(targets[[sexColumn]])),
      text.col = RColorBrewer::brewer.pal(8, "Dark2")
    )
  }

  runPlotMinfiEwasWater(
    draw_fun = draw_fun,
    display = display,
    file = file,
    width = width,
    height = height,
    res = res
  )

  emitLogMinfiEwasWater(
    c(
      if (is.null(file)) {
        "Raw-vs-normalized plot:   not written to file"
      } else {
        paste("Raw-vs-normalized plot:   ", file)
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  invisible(file)
}

#' Plot raw beta-value density from a raw preprocessing object
#'
#' Draw the pre-normalization beta density plot from a raw minfi object and a
#' grouping variable in the phenotype table.
#'
#' @param rawData Object returned by `buildRawMinfiEwasWater()`.
#' @param targets Filtered phenotype data aligned with `rawData`.
#' @param plotGroupVar Character. Phenotype column used to group samples in the
#'   density plot.
#' @param display Logical. If `TRUE`, draw the plot on the active graphics
#'   device.
#' @param file Character or `NULL`. TIFF file written when supplied.
#' @param width Integer. TIFF width in pixels when `file` is supplied.
#' @param height Integer. TIFF height in pixels when `file` is supplied.
#' @param res Integer. TIFF resolution in DPI when `file` is supplied.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return Invisibly returns the saved TIFF path when `file` is supplied,
#'   otherwise `NULL`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiWorkflowStateDnaEpico()
#' plotRawDensityMinfiEwasWater(
#'   rawData = ex$rawFiltered,
#'   targets = ex$sampleData$targets,
#'   plotGroupVar = "Sex",
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#'
#' @export
plotRawDensityMinfiEwasWater <- function(
    rawData,
    targets,
    plotGroupVar = "Sex",
    display = FALSE,
    file = NULL,
    width = 2000L,
    height = 1000L,
    res = 150L,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_plotRawDensityMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  draw_fun <- function() {
    minfi::densityPlot(
      minfi::getBeta(rawData$MSet),
      sampGroups = targets[[plotGroupVar]],
      pal = RColorBrewer::brewer.pal(8, "Dark2"),
      main = paste("Density Plot of Beta Values by", plotGroupVar),
      legend = TRUE
    )
  }

  runPlotMinfiEwasWater(
    draw_fun = draw_fun,
    display = display,
    file = file,
    width = width,
    height = height,
    res = res
  )

  emitLogMinfiEwasWater(
    c(
      if (is.null(file)) {
        "Raw beta density plot:    not written to file"
      } else {
        paste("Raw beta density plot:    ", file)
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  invisible(file)
}

#' Filter probes from a normalized methylation object
#'
#' Apply detection P-value, chromosome, SNP, and cross-reactive probe filters
#' to the primary normalized object and return the filtered result.
#'
#' @param normData Object returned by `normalizeMinfiEwasWater()`.
#' @param RGSet Filtered `RGChannelSet` aligned with `normData`.
#' @param pvalThreshold Numeric. Probes must have detection P values below this
#'   threshold in all samples to be retained.
#' @param chrToRemove Character vector or comma-separated string of chromosome
#'   names to remove, for example `"chrX,chrY"`.
#' @param snpsToRemove Character vector or comma-separated string of SNP probe
#'   types to remove, for example `"SBE,CpG"`.
#' @param mafThreshold Numeric. Minor allele frequency threshold passed to
#'   `minfi::dropLociWithSnps()`.
#' @param crossReactivePath Character. Path to a CSV file containing a `ProbeID`
#'   column of cross-reactive probes to remove.
#' @param detPtype Character. Detection P-value mode passed to
#'   `minfi::detectionP()` for the probe filter. Common values in minfi
#'   workflows are `"m+u"` and `"negative"`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_minfiEwasWater_filter"` containing the
#'   filtered object and counts for each filtering stage.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiWorkflowStateDnaEpico()
#' filtered_data <- filterProbesMinfiEwasWater(
#'   normData = ex$normData,
#'   RGSet = ex$sampleData$RGSet,
#'   pvalThreshold = 1,
#'   chrToRemove = "chrY",
#'   snpsToRemove = "SBE",
#'   mafThreshold = 1,
#'   crossReactivePath = ex$crossReactivePath,
#'   detPtype = "m+u",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' filtered_data$counts[["crossReactive"]]
#'
#' @export
filterProbesMinfiEwasWater <- function(
    normData,
    RGSet,
    pvalThreshold = 0.01,
    chrToRemove = "chrX,chrY",
    snpsToRemove = "SBE,CpG",
    mafThreshold = 0.1,
    crossReactivePath,
    detPtype = "m+u",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_filterProbesMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  if (!file.exists(crossReactivePath)) {
    stop(
      "crossReactivePath does not exist: ",
      crossReactivePath,
      call. = FALSE
    )
  }

  chr_list <- splitOptionMinfiEwasWater(chrToRemove, sep = ",")
  snp_list <- splitOptionMinfiEwasWater(snpsToRemove, sep = ",")

  detP <- minfi::detectionP(RGSet, type = detPtype)
  detP <- detP[
    match(Biobase::featureNames(normData$primary), rownames(detP)),
    ,
    drop = FALSE
  ]

  keep_detp <- rowSums(detP < pvalThreshold) == ncol(normData$primary)
  filtered_detp <- normData$primary[keep_detp, ]

  ann <- minfi::getAnnotation(RGSet)
  remove_probes <- ann$Name[ann$chr %in% chr_list]
  keep_chr <- !(Biobase::featureNames(filtered_detp) %in% remove_probes)
  filtered_chr <- filtered_detp[keep_chr, ]

  filtered_snp <- minfi::dropLociWithSnps(
    filtered_chr,
    snps = snp_list,
    maf = mafThreshold
  )

  cross_reactive <- utils::read.csv(
    crossReactivePath,
    stringsAsFactors = FALSE
  )
  keep_cross <- !(Biobase::featureNames(filtered_snp) %in% cross_reactive$ProbeID)
  filtered_final <- filtered_snp[keep_cross, ]

  counts <- c(
    start = nrow(normData$primary),
    detP = nrow(filtered_detp),
    chromosome = nrow(filtered_chr),
    snp = nrow(filtered_snp),
    crossReactive = nrow(filtered_final)
  )

  emitLogMinfiEwasWater(
    c(
      paste("Probe filter threshold:    ", pvalThreshold),
      paste("Chromosomes removed:       ", paste(chr_list, collapse = ", ")),
      paste("SNP filters removed:       ", paste(snp_list, collapse = ", ")),
      paste("MAF threshold:             ", mafThreshold),
      paste("Cross-reactive file:       ", crossReactivePath),
      paste("Probes after detP filter:  ", counts[["detP"]]),
      paste("Probes after chr filter:   ", counts[["chromosome"]]),
      paste("Probes after SNP filter:   ", counts[["snp"]]),
      paste("Probes after cross filter: ", counts[["crossReactive"]]),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      filtered = filtered_final,
      detPFiltered = filtered_detp,
      chrFiltered = filtered_chr,
      snpFiltered = filtered_snp,
      detP = detP,
      counts = counts,
      pvalThreshold = pvalThreshold,
      chrToRemove = chr_list,
      snpsToRemove = snp_list,
      mafThreshold = mafThreshold,
      crossReactivePath = crossReactivePath
    ),
    class = "dnaEPICO_minfiEwasWater_filter"
  )
}

#' Extract beta, M, and copy-number matrices from a filtered object
#'
#' @param filteredData Object returned by `filterProbesMinfiEwasWater()`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_minfiEwasWater_metrics"` containing
#'   `beta`, `m`, and `cn`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiWorkflowStateDnaEpico()
#' metrics_data <- extractMetricsMinfiEwasWater(
#'   filteredData = ex$filteredData,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(metrics_data)
#'
#' @export
extractMetricsMinfiEwasWater <- function(
    filteredData,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_extractMetricsMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  beta <- minfi::getBeta(filteredData$filtered)
  m <- minfi::getM(filteredData$filtered)
  cn <- minfi::getCN(filteredData$filtered)

  preview_cols <- seq_len(min(ncol(beta), 5L))

  emitLogMinfiEwasWater(
    c(
      "Extracting final DNAm matrices (M, Beta, CN).",
      "Preview of beta values:",
      previewLinesMinfiEwasWater(beta[, preview_cols, drop = FALSE]),
      "Preview of M-values:",
      previewLinesMinfiEwasWater(m[, preview_cols, drop = FALSE]),
      "Preview of copy-number values:",
      previewLinesMinfiEwasWater(cn[, preview_cols, drop = FALSE]),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      beta = beta,
      m = m,
      cn = cn
    ),
    class = "dnaEPICO_minfiEwasWater_metrics"
  )
}

#' Plot multidimensional scaling or density summaries from final metrics
#'
#' @param metricsData Object returned by `extractMetricsMinfiEwasWater()`.
#' @param targets Filtered phenotype data aligned with `metricsData`.
#' @param plot Character. Plot type: `"mds"` or `"density"`.
#' @param plotGroupVar Character. Phenotype column used for the main grouping.
#' @param sexColumn Character. Phenotype column used for the sex grouping in the
#'   MDS plot.
#' @param display Logical. If `TRUE`, draw the plot on the active graphics
#'   device.
#' @param file Character or `NULL`. TIFF file written when supplied.
#' @param width Integer. TIFF width in pixels when `file` is supplied.
#' @param height Integer. TIFF height in pixels when `file` is supplied.
#' @param res Integer. TIFF resolution in DPI when `file` is supplied.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return Invisibly returns the saved TIFF path when `file` is supplied,
#'   otherwise `NULL`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiWorkflowStateDnaEpico()
#' plotMetricsMinfiEwasWater(
#'   metricsData = ex$metricsData,
#'   targets = ex$sampleData$targets,
#'   plot = "density",
#'   plotGroupVar = "Sex",
#'   sexColumn = "Sex",
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#'
#' @export
plotMetricsMinfiEwasWater <- function(
    metricsData,
    targets,
    plot = c("mds", "density"),
    plotGroupVar = "Sex",
    sexColumn = "Sex",
    display = FALSE,
    file = NULL,
    width = 2000L,
    height = 1000L,
    res = 150L,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_plotMetricsMinfiEwasWater.txt"
) {
  plot <- match.arg(plot)

  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  draw_fun <- switch(
    plot,
    mds = function() {
      group_factor <- factor(targets[[plotGroupVar]])
      group_sex <- factor(targets[[sexColumn]])
      pal <- RColorBrewer::brewer.pal(8, "Dark2")

      graphics::par(mfrow = c(1, 2))
      limma::plotMDS(
        metricsData$m,
        main = plotGroupVar,
        top = 1000,
        gene.selection = "common",
        col = pal[group_factor],
        dim = c(1, 2)
      )
      graphics::legend(
        "right",
        legend = levels(group_factor),
        text.col = pal,
        cex = 0.7,
        bg = "white"
      )

      limma::plotMDS(
        metricsData$m,
        main = sexColumn,
        top = 1000,
        gene.selection = "common",
        col = pal[group_sex],
        dim = c(2, 3)
      )
      graphics::legend(
        "topright",
        legend = levels(group_sex),
        text.col = pal,
        cex = 0.7,
        bg = "white"
      )
    },
    density = function() {
      group_factor <- factor(targets[[plotGroupVar]])
      graphics::par(mfrow = c(1, 2))
      minfi::densityPlot(
        metricsData$beta,
        sampGroups = group_factor,
        main = "Beta values",
        legend = FALSE,
        xlab = "Beta values"
      )
      graphics::legend(
        "top",
        legend = levels(group_factor),
        text.col = RColorBrewer::brewer.pal(8, "Dark2")
      )

      minfi::densityPlot(
        metricsData$m,
        sampGroups = group_factor,
        main = "M-values",
        legend = FALSE,
        xlab = "M values"
      )
      graphics::legend(
        "topleft",
        legend = levels(group_factor),
        text.col = RColorBrewer::brewer.pal(8, "Dark2")
      )
    }
  )

  runPlotMinfiEwasWater(
    draw_fun = draw_fun,
    display = display,
    file = file,
    width = width,
    height = height,
    res = res
  )

  emitLogMinfiEwasWater(
    c(
      paste("Metrics plot created:     ", plot),
      if (is.null(file)) {
        "Metrics plot file:        none"
      } else {
        paste("Metrics plot file:        ", file)
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  invisible(file)
}

#' Plot ENmix control images from an RGSet
#'
#' Call `ENmix::plotCtrl()` for a supplied `RGSet`. This function only writes
#' files when `output_dir` is provided because `ENmix::plotCtrl()` produces JPG
#' files on disk rather than returning a plot object.
#'
#' @param RGSet An `RGChannelSet`.
#' @param output_dir Character or `NULL`. Directory where ENmix control JPG
#'   files should be written. If `NULL`, the function returns without writing
#'   files.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return Invisibly returns `output_dir`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' output_dir <- file.path(tempdir(), "enmix-control-plots")
#' plotCtrlMinfiEwasWater(
#'   RGSet = ex$RGSet,
#'   output_dir = output_dir,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' dir.exists(output_dir)
#'
#' @export
plotCtrlMinfiEwasWater <- function(
    RGSet,
    output_dir = NULL,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_plotCtrlMinfiEwasWater.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = log_dir,
    log_file = log_file
  )

  if (is.null(output_dir)) {
    emitLogMinfiEwasWater(
      c(
        "ENmix control plots skipped because output_dir is NULL.",
        "======================================================================="
      ),
      verbose = verbose,
      log_path = log_path
    )

    return(invisible(NULL))
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(output_dir)

  old_options <- options(bitmapType = "cairo")
  on.exit(options(old_options), add = TRUE)

  ENmix::plotCtrl(RGSet)

  emitLogMinfiEwasWater(
    c(
      paste("Generated ENmix control JPGs in:", output_dir),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  invisible(output_dir)
}
