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
#'   `'IlluminaHumanMethylationEPICv2'`.
#' @param annotationVersion Character. Annotation build passed to
#'   `Biobase::annotation(RGSet)`, for example `'20a1.hg38'` for EPIC v2 hg38
#'   annotations or `'ilmn12.hg19'` for 450K hg19 annotations.
#' @param force Logical. Passed to `minfi::read.metharray.exp()`. Use `TRUE`
#'   only after confirming that the selected IDAT files should be read together.
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
#'     ex <- dnaEPICO:::exampleMinfiIdatInputsDnaEpico(n = 4)
#'     rgset <- readRGSetMinfiEwasWater(
#'         idatFolder = ex$idatFolder,
#'         targets = ex$targets,
#'         SampleID = "Sample_Name",
#'         arrayType = ex$arrayType,
#'         annotationVersion = ex$annotationVersion
#'     )
#'     class(rgset)
#' }
#'
#' @export
readRGSetMinfiEwasWater <- function(
    idatFolder, targets, SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
        annotationVersion = "20a1.hg38",
    force = FALSE, verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_readRGSetMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    if (!dir.exists(idatFolder)) {
        stop("idatFolder does not exist: ", idatFolder, call. = FALSE)
    }

    if (!(SampleID %in% colnames(targets))) {
        stop("SampleID column not found in targets: ", SampleID,
            call. = FALSE
        )
    }
    force <- validateLogicalScalarDnaEpico(force, "force")

    emitLogMinfiEwasWater(
        c(
            "============================================================",
            paste("IDAT folder:              ", idatFolder), paste(
                "Array type:               ",
                arrayType
            ), paste("Annotation version:       ", annotationVersion),
            paste("Force IDAT read:          ", force)
        ),
        verbose = verbose,
        log_path = log_path
    )

    RGSet <- minfi::read.metharray.exp(
        base = idatFolder, targets = targets,
        extended = FALSE, recursive = FALSE, verbose = FALSE,
        force = force
    )

    sample_ids <- validateSampleIdentifiersDnaEpico(
        targets[[SampleID]],
        paste0("Phenotype column '", SampleID, "'")
    )
    if (ncol(RGSet) != length(sample_ids)) {
        stop("The IDAT reader returned ", ncol(RGSet), " samples for ",
            length(sample_ids), " phenotype rows.",
            call. = FALSE
        )
    }

    colnames(RGSet) <- sample_ids
    Biobase::annotation(RGSet) <- c(array = arrayType,
        annotation = annotationVersion)

    manifest_lines <-
        utils::capture.output(methods::show(minfi::getManifest(RGSet)))

    emitLogMinfiEwasWater(
        c(
            paste(
                "RGSet loaded with", ncol(RGSet),
                "samples."
            ), paste("Applied annotation:       ",
                paste(Biobase::annotation(RGSet),
                collapse = ", "
            )), "Manifest used:", manifest_lines,
            "============================================================"
        ),
        verbose = verbose, log_path = log_path
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
#' @return A list with class `'dnaEPICO_minfiEwasWater_raw'` containing `MSet`,
#'   `RatioSet`, and `GSet`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' raw_data <- buildRawMinfiEwasWater(
#'     RGSet = ex$RGSet,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(raw_data)
#'
#' @export
buildRawMinfiEwasWater <- function(
    RGSet, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_buildRawMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    emitLogMinfiEwasWater("Running preprocessRaw(), ratioConvert(), and mapToGenome().",
        verbose = verbose, log_path = log_path
    )

    MSet <- minfi::preprocessRaw(RGSet)
    RatioSet <- minfi::ratioConvert(MSet, what = "both", keepCN = TRUE)
    GSet <- minfi::mapToGenome(RatioSet)

    meth_cols <- seq_len(min(ncol(MSet), 3L))
    gset_cols <- seq_len(min(ncol(GSet), 5L))

    emitLogMinfiEwasWater(
        c(
            paste(
                "MSet dimensions:          ",
                paste(dim(MSet), collapse = " x ")
            ), paste(
                "RatioSet dimensions:      ",
                paste(dim(RatioSet), collapse = " x ")
            ), paste(
                "GSet dimensions:          ",
                paste(dim(GSet), collapse = " x ")
            ), "Preview of methylated intensities:",
            previewLinesMinfiEwasWater(minfi::getMeth(MSet)[, meth_cols,
                drop = FALSE
            ]), "Preview of unmethylated intensities:",
            previewLinesMinfiEwasWater(minfi::getUnmeth(MSet)[, meth_cols,
                drop = FALSE
            ]), "Preview of beta values:",
                previewLinesMinfiEwasWater(minfi::getBeta(GSet)[,
                gset_cols,
                drop = FALSE
            ]), "Preview of M-values:",
            previewLinesMinfiEwasWater(minfi::getM(GSet)[, gset_cols,
                drop = FALSE
            ]), "Preview of copy-number values:",
            previewLinesMinfiEwasWater(minfi::getCN(GSet)[, gset_cols,
                drop = FALSE
            ]), "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(MSet = MSet, RatioSet = RatioSet, GSet = GSet),
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
#'   `minfi::detectionP()`. Common values in minfi workflows are `'m+u'` and
#'   `'negative'`. The default used here is `'m+u'`.
#' @param detPThreshold Numeric. Samples with mean detection P value above this
#'   threshold are marked as failed.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_minfiEwasWater_assessment'` containing
#'   the QC object, detection P matrix, mean detection P values, and failed
#'   sample identifiers.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' raw_data <- buildRawMinfiEwasWater(ex$RGSet, verbose = FALSE, logs = FALSE)
#' assessment <- assessSamplesMinfiEwasWater(
#'     rawData = raw_data,
#'     RGSet = ex$RGSet,
#'     detPThreshold = 1,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(assessment)
#'
#' @export
assessSamplesMinfiEwasWater <- function(
    rawData, RGSet, qcCutoff = 10.5,
    detPtype = "m+u", detPThreshold = 0.05, verbose = FALSE,
    logs = FALSE, log_dir = NULL,
        log_file = "log_assessSamplesMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    qc <- minfi::getQC(rawData$MSet)
    detP <- minfi::detectionP(RGSet, type = detPtype)
    detPThreshold <- validateProbabilityDnaEpico(
        detPThreshold,
        "detPThreshold"
    )
    nanDetPConverted <- sum(is.nan(detP))
    if (nanDetPConverted > 0L) {
        detP[is.nan(detP)] <- NA_real_
    }

    missingDetP <- colSums(!is.finite(detP))
    observedDetP <- colSums(is.finite(detP))
    meanDetP <- vapply(seq_len(ncol(detP)), function(index) {
        meanFiniteOrNADnaEpico(detP[
            ,
            index
        ])
    }, numeric(1))
    names(meanDetP) <- colnames(detP)
    failedSamples <- names(meanDetP[!is.finite(meanDetP) | meanDetP >
        detPThreshold])

    preview_cols <- seq_len(min(ncol(detP), 5L))

    emitLogMinfiEwasWater(
        c(
            paste(
                "QC cutoff (median):       ",
                qcCutoff
            ), paste("Detection P type:         ", detPtype),
            paste("Detection P threshold:    ", detPThreshold), paste(
                "Detection P NaN to NA:     ",
                nanDetPConverted
            ), paste(
                "Missing detection values: ",
                sum(missingDetP)
            ), "Preview of detection P values:",
            previewLinesMinfiEwasWater(detP[, preview_cols, drop = FALSE]),
            paste("Number of failed samples: ", length(failedSamples)),
            if (length(failedSamples) > 0L) {
                paste("Failed sample IDs:        ", paste(failedSamples,
                    collapse = ", "
                ))
            } else {
                "Failed sample IDs:         none"
            }, "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(
        list(
            qc = qc, qcCutoff = qcCutoff, detP = detP,
            detPtype = detPtype, detPThreshold = detPThreshold,
                nanDetPConverted = nanDetPConverted,
            observedDetP = observedDetP, missingDetP = missingDetP,
            meanDetP = meanDetP, failedSamples = failedSamples
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
#' @param plot Character. Plot type: `'qc'` or `'detection'`.
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
#' @examples
#' assessment <- list(
#'     meanDetP = c(S1 = 0.01, S2 = 0.02, S3 = 0.04),
#'     detPThreshold = 0.05
#' )
#' plotAssessmentMinfiEwasWater(
#'     assessment = assessment,
#'     plot = "detection",
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#'
#' @export
plotAssessmentMinfiEwasWater <- function(
    assessment, plot = c(
        "qc",
        "detection"
    ), display = FALSE, file = NULL, width = 2000L,
    height = 1000L, res = 150L, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_plotAssessmentMinfiEwasWater.txt"
) {
    plot <- match.arg(plot)

    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    draw_fun <- switch(plot,
        qc = function() {
            minfi::plotQC(assessment$qc, badSampleCutoff = assessment$qcCutoff)
        },
        detection = function() {
            graphics::barplot(assessment$meanDetP,
                las = 3, cex.names = 0.8,
                ylab = "Mean detection p-values"
            )
            graphics::abline(
                h = assessment$detPThreshold, col = "red",
                lwd = 2, lty = 2
            )
        }
    )

    runPlotMinfiEwasWater(
        draw_fun = draw_fun, display = display,
        file = file, width = width, height = height, res = res
    )

    emitLogMinfiEwasWater(
        c(paste(
            "Assessment plot created:   ",
            plot
        ), if (is.null(file)) {
            "Assessment plot file:      none"
        } else {
            paste("Assessment plot file:      ", file)
        }, "============================================================"),
        verbose = verbose, log_path = log_path
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
#' @return A list with class `'dnaEPICO_minfiEwasWater_samples'` containing the
#'   filtered `RGSet`, aligned phenotype table, and failed sample identifiers.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' filtered_samples <- filterSamplesMinfiEwasWater(
#'     RGSet = ex$RGSet,
#'     targets = ex$targets,
#'     failedSamples = ex$targets$Sample_Name[1],
#'     SampleID = "Sample_Name",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' nrow(filtered_samples$targets)
#'
#' @export
filterSamplesMinfiEwasWater <- function(
    RGSet, targets, failedSamples = character(0),
    SampleID = "Sample_Name", verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_filterSamplesMinfiEwasWater.txt"
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
    sample_names <- validateSampleIdentifiersDnaEpico(
        colnames(RGSet),
        "RGSet sample identifiers"
    )
    target_sample_ids <- validateSampleIdentifiersDnaEpico(
        targets[[SampleID]],
        paste0("phenotype column '", SampleID, "'")
    )
    matchSampleIdentifiersDnaEpico(
        query = sample_names, reference = target_sample_ids,
        queryLabel = "RGSet sample identifiers", referenceLabel = paste0(
            "phenotype column '",
            SampleID, "'"
        ), requireSameSet = TRUE
    )
    if (length(failedSamples) > 0L) {
        failedSamples <- validateSampleIdentifiersDnaEpico(
            failedSamples,
            "failedSamples"
        )
        unknown_failed_samples <- setdiff(failedSamples, sample_names)
        if (length(unknown_failed_samples) > 0L) {
            stop("failedSamples are not present in RGSet: ",
                paste(unknown_failed_samples, collapse = ", "),
                call. = FALSE
            )
        }
    }

    keep_samples <- !(sample_names %in% failedSamples)
    if (!any(keep_samples)) {
        stop("Sample filtering would remove every sample.", call. = FALSE)
    }
    RGSet_filtered <- RGSet[, keep_samples]

    sample_names <- colnames(RGSet_filtered)
    matched <- matchSampleIdentifiersDnaEpico(
        query = sample_names,
        reference = targets[[SampleID]],
            queryLabel = "RGSet sample identifiers",
        referenceLabel = paste0(
            "phenotype column '", SampleID,
            "'"
        )
    )
    targets_filtered <- targets[matched, , drop = FALSE]
    rownames(targets_filtered) <- NULL

    emitLogMinfiEwasWater(
        c(
            paste(
                "Samples before filtering: ",
                ncol(RGSet)
            ), paste("Samples after filtering:  ", ncol(RGSet_filtered)),
            "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(
        RGSet = RGSet_filtered, targets = targets_filtered,
        failedSamples = failedSamples
    ), class = "dnaEPICO_minfiEwasWater_samples")
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
#' @return A list with class `'dnaEPICO_minfiEwasWater_sex'` containing the sex
#'   prediction result, aligned phenotype data, plotting data, and mismatch
#'   table.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiWorkflowStateDnaEpico()
#' sex_data <- predictSexMinfiEwasWater(
#'     rawData = ex$rawFiltered,
#'     targets = ex$sampleData$targets,
#'     SampleID = "Sample_Name",
#'     sexColumn = "Sex",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(sex_data)
#'
#' @export
predictSexMinfiEwasWater <- function(
    rawData, targets, SampleID = "Sample_Name",
    sexColumn = "Sex", verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_predictSexMinfiEwasWater.txt"
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

    if (!(sexColumn %in% colnames(targets))) {
        stop("sexColumn not found in targets: ", sexColumn, call. = FALSE)
    }

    pSex <- as.data.frame(minfi::getSex(rawData$GSet), stringsAsFactors = FALSE)
    numeric_sex_columns <- names(pSex)[vapply(
        pSex, is.numeric,
        logical(1)
    )]
    nanSexConverted <- sum(vapply(
        pSex[numeric_sex_columns],
        function(values) sum(is.nan(values)), integer(1)
    ))
    for (column in numeric_sex_columns) {
        pSex[[column]][is.nan(pSex[[column]])] <- NA_real_
    }
    sex_plot_data <- as.data.frame(pSex, stringsAsFactors = FALSE)
    sex_plot_data$SampleID <- rownames(sex_plot_data)

    matched <- matchSampleIdentifiersDnaEpico(
        query = sex_plot_data$SampleID,
        reference = targets[[SampleID]],
            queryLabel = "Predicted-sex sample identifiers",
        referenceLabel = paste0(
            "phenotype column '", SampleID,
            "'"
        ), requireSameSet = TRUE
    )
    targets_aligned <- targets[matched, , drop = FALSE]
    rownames(targets_aligned) <- NULL

    reported_sex <- canonicalizeSexDnaEpico(targets_aligned[[sexColumn]])
    predicted_sex <- canonicalizeSexDnaEpico(pSex$predictedSex)
    targets_aligned$PredSex <- predicted_sex$code
    targets_aligned$SexMismatch <- !is.na(reported_sex$code) &
        !is.na(predicted_sex$code) & reported_sex$code != predicted_sex$code

    sex_plot_data[[sexColumn]] <- reported_sex$code
    sex_plot_data$PredSex <- predicted_sex$code
    sex_plot_data$SexMismatch <- targets_aligned$SexMismatch
    mismatches <- targets_aligned[targets_aligned$SexMismatch, ,
        drop = FALSE
    ]

    sex_values <- targets[[sexColumn]]
    nSexNA <- sum(is.na(sex_values))
    unknown_sex <- reported_sex$unknown

    emitLogMinfiEwasWater(
        c(
            "Sex column integrity check:", paste(
                "  Sex column used:        ",
                sexColumn
            ), paste("  NA values:              ", nSexNA),
            paste("  NaN values converted:   ", nanSexConverted),
            paste("  Unknown labels:         ", if (length(unknown_sex) ==
                0L) {
                "none"
            } else {
                paste(unknown_sex, collapse = ", ")
            }), paste("Mismatches found:         ", nrow(mismatches)),
            if (nrow(mismatches) > 0L) {
                previewLinesMinfiEwasWater(mismatches[, seq_len(min(
                    3L,
                    ncol(mismatches)
                )), drop = FALSE])
            } else {
                character(0)
            }, "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(
        list(
            pSex = pSex, targets = targets_aligned, sexPlotData = sex_plot_data,
            mismatches = mismatches,
                reportedSexOriginal = reported_sex$original,
            reportedSexCode = reported_sex$code,
                reportedSexLabel = reported_sex$label,
            predictedSexCode = predicted_sex$code,
                predictedSexLabel = predicted_sex$label,
            unknownReportedSex = reported_sex$unknown,
                unknownPredictedSex = predicted_sex$unknown,
            nanConverted = nanSexConverted, removedSampleIDs = character(0),
            removeSexMismatch = FALSE, SampleID = SampleID,
                sexColumn = sexColumn
        ),
        class = "dnaEPICO_minfiEwasWater_sex"
    )
}

#' Plot predicted or clinical sex from `predictSexMinfiEwasWater()`
#'
#' @param sexData Object returned by `predictSexMinfiEwasWater()`.
#' @param type Character. Plot type: `'predicted'` for methylation-predicted
#'   sex or `'clinical'` for reported sex.
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
#' @examples
#' ex <- dnaEPICO:::exampleSexPlotStateDnaEpico()
#' plotSexMinfiEwasWater(
#'     sexData = ex,
#'     type = "predicted",
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#'
#' @export
plotSexMinfiEwasWater <- function(
    sexData, type = c(
        "predicted",
        "clinical"
    ), display = FALSE, file = NULL, width = 2000L,
    height = 1000L, res = 70L, verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_plotSexMinfiEwasWater.txt"
) {
    type <- match.arg(type)

    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    draw_fun <- switch(type,
        predicted = function() {
            graphics::plot(
                x = sexData$pSex$xMed, y = sexData$pSex$yMed,
                type = "n", xlab = "X chr, median total intensity (log2)",
                ylab = "Y chr, median total intensity (log2)"
            )
            graphics::text(
                x = sexData$pSex$xMed, y = sexData$pSex$yMed,
                labels = sexData$targets[[sexData$SampleID]],
                    col = ifelse(sexData$pSex$predictedSex ==
                    "M", "deepskyblue", "deeppink3")
            )
            graphics::legend("bottomleft", c("M", "F"), col = c(
                "deepskyblue",
                "deeppink3"
            ), pch = 16)
        },
        clinical = function() {
            graphics::plot(
                x = sexData$sexPlotData$xMed, y = sexData$sexPlotData$yMed,
                type = "n", xlab = "X chr, median total intensity (log2)",
                ylab = "Y chr, median total intensity (log2)"
            )
            graphics::text(
                x = sexData$sexPlotData$xMed, y = sexData$sexPlotData$yMed,
                labels = sexData$sexPlotData$SampleID,
                    col = ifelse(sexData$sexPlotData[[sexData$sexColumn]] ==
                    1L, "deepskyblue", "deeppink3")
            )
            graphics::legend("bottomleft", c("M", "F"), col = c(
                "deepskyblue",
                "deeppink3"
            ), pch = 16)
        }
    )

    runPlotMinfiEwasWater(
        draw_fun = draw_fun, display = display,
        file = file, width = width, height = height, res = res
    )

    emitLogMinfiEwasWater(
        c(paste(
            "Sex plot created:         ",
            type
        ), if (is.null(file)) {
            "Sex plot file:            none"
        } else {
            paste("Sex plot file:            ", file)
        }, "============================================================"),
        verbose = verbose, log_path = log_path
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
#'   normalization methods. Supported values are `'adjustedfunnorm'`,
#'   `'funnorm'`, `'illumina'`, `'quantile'`, and `'swan'`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_minfiEwasWater_norm'` containing the
#'   requested normalized objects and the first method as `primary`.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
#' sample_data <- filterSamplesMinfiEwasWater(
#'     RGSet = ex$RGSet,
#'     targets = ex$targets,
#'     failedSamples = character(0),
#'     SampleID = "Sample_Name",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' norm_data <- normalizeMinfiEwasWater(
#'     sampleData = sample_data,
#'     sexColumn = "Sex",
#'     normMethods = "quantile",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(norm_data$normalized)
#'
#' @export
normalizeMinfiEwasWater <- function(
    sampleData, sexColumn = "Sex",
    normMethods = "adjustedfunnorm", verbose = FALSE, logs = FALSE,
    log_dir = NULL, log_file = "log_normalizeMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    method_list <- tolower(splitOptionMinfiEwasWater(normMethods,
        sep = ";"
    ))
    supported_methods <- c(
        "adjustedfunnorm", "funnorm", "illumina",
        "quantile", "swan"
    )
    if (length(method_list) == 0L) {
        stop("At least one normalization method must be supplied.",
            call. = FALSE
        )
    }
    unknown_methods <- setdiff(method_list, supported_methods)
    if (length(unknown_methods) > 0L) {
        stop("Unknown normalization method(s): ", paste(unknown_methods,
            collapse = ", "
        ), call. = FALSE)
    }
    normalized <- vector("list", length(method_list))
    names(normalized) <- method_list

    col_data <- SummarizedExperiment::colData(sampleData$RGSet)
    sex_vec <- resolveNormalizationSexDnaEpico(col_data, sexColumn)

    emitLogMinfiEwasWater(c(paste(
        "Normalization methods:    ",
        paste(method_list, collapse = ", ")
    ), paste(
        "Sex column:               ",
        sexColumn
    )), verbose = verbose, log_path = log_path)

    for (i in seq_along(method_list)) {
        method <- method_list[[i]]
        emitLogMinfiEwasWater(paste(
            "Applying normalization:   ",
            method
        ), verbose = verbose, log_path = log_path)

        normalized[[i]] <- switch(method,
            adjustedfunnorm = wateRmelon::adjustedFunnorm(sampleData$RGSet,
                sex = sex_vec
            ),
            funnorm = minfi::preprocessFunnorm(sampleData$RGSet,
                sex = sex_vec
            ),
            illumina = minfi::preprocessIllumina(sampleData$RGSet),
            quantile = minfi::preprocessQuantile(sampleData$RGSet,
                sex = sex_vec
            ),
            swan = minfi::preprocessSWAN(sampleData$RGSet),
            stop("Unknown normalization method: ", method, call. = FALSE)
        )
    }

    emitLogMinfiEwasWater("============================================================",
        verbose = verbose, log_path = log_path
    )

    structure(list(
        primary = normalized[[1L]], normalized = normalized,
        methods = method_list, sexColumn = sexColumn
    ), class = "dnaEPICO_minfiEwasWater_norm")
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
#' @examples
#' ex <- dnaEPICO:::exampleMinfiMetricsStateDnaEpico()
#' plotNormalizationMinfiEwasWater(
#'     RGSet = ex$beta,
#'     normData = ex$normData,
#'     targets = ex$targets,
#'     sexColumn = "Sex",
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#'
#' @export
plotNormalizationMinfiEwasWater <- function(
    RGSet, normData,
    targets, sexColumn = "Sex", display = FALSE, file = NULL,
    width = 2000L, height = 1000L, res = 150L, verbose = FALSE,
    logs = FALSE, log_dir = NULL,
        log_file = "log_plotNormalizationMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    draw_fun <- function() {
        graphics::par(mfrow = c(1, 2))
        minfi::densityPlot(RGSet,
            sampGroups = targets[[sexColumn]],
            main = "Raw", legend = FALSE
        )
        graphics::legend("top",
            legend = levels(factor(targets[[sexColumn]])),
            text.col = RColorBrewer::brewer.pal(8, "Dark2")
        )

        minfi::densityPlot(minfi::getBeta(normData$primary),
            sampGroups = targets[[sexColumn]], main = "Normalized",
            legend = FALSE
        )
        graphics::legend("top",
            legend = levels(factor(targets[[sexColumn]])),
            text.col = RColorBrewer::brewer.pal(8, "Dark2")
        )
    }

    runPlotMinfiEwasWater(
        draw_fun = draw_fun, display = display,
        file = file, width = width, height = height, res = res
    )

    emitLogMinfiEwasWater(
        c(if (is.null(file)) {
            "Raw-vs-normalized plot:   not written to file"
        } else {
            paste("Raw-vs-normalized plot:   ", file)
        }, "============================================================"),
        verbose = verbose, log_path = log_path
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
#' @examples
#' ex <- dnaEPICO:::exampleMinfiMetricsStateDnaEpico()
#' plotRawDensityMinfiEwasWater(
#'     rawData = ex$rawData,
#'     targets = ex$targets,
#'     plotGroupVar = "Sex",
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#'
#' @export
plotRawDensityMinfiEwasWater <- function(
    rawData, targets, plotGroupVar = "Sex",
    display = FALSE, file = NULL, width = 2000L, height = 1000L,
    res = 150L, verbose = FALSE, logs = FALSE, log_dir = NULL,
    log_file = "log_plotRawDensityMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    draw_fun <- function() {
        minfi::densityPlot(minfi::getBeta(rawData$MSet),
            sampGroups = targets[[plotGroupVar]],
            pal = RColorBrewer::brewer.pal(8, "Dark2"), main = paste(
                "Density Plot of Beta Values by",
                plotGroupVar
            ), legend = TRUE
        )
    }

    runPlotMinfiEwasWater(
        draw_fun = draw_fun, display = display,
        file = file, width = width, height = height, res = res
    )

    emitLogMinfiEwasWater(
        c(if (is.null(file)) {
            "Raw beta density plot:    not written to file"
        } else {
            paste("Raw beta density plot:    ", file)
        }, "============================================================"),
        verbose = verbose, log_path = log_path
    )

    invisible(file)
}

#' Read probe-exclusion identifiers from one CSV file
#'
#' @param probeExclusionPath Character. Path to a CSV file containing probe
#'   identifiers to exclude.
#' @param probeExclusionIdColumn Character or `NULL`. Column containing probe
#'   identifiers. When `NULL` or `''`, common names are auto-detected before
#'   falling back to an unlabeled or probe-like first column.
#' @param featureNames Character vector or `NULL`. Optional array feature names
#'   used to count overlap.
#'
#' @return A list with class `'dnaEPICO_probeExclusion_file_ids'`.
#'
#' @keywords internal
#' @noRd
readProbeExclusionFileIdsMinfiEwasWater <- function(
    probeExclusionPath,
    probeExclusionIdColumn = NULL, featureNames = NULL
) {
    if (!file.exists(probeExclusionPath)) {
        stop("probeExclusionPath does not exist: ", probeExclusionPath,
            call. = FALSE
        )
    }

    probe_exclusion <- utils::read.csv(probeExclusionPath,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (ncol(probe_exclusion) == 0L) {
        stop("probeExclusionPath contains no columns: ", probeExclusionPath,
            call. = FALSE
        )
    }

    id_column <- probeExclusionIdColumn
    if (length(id_column) == 0L || is.na(id_column[[1L]])) {
        id_column <- NULL
    } else {
        id_column <- trimws(as.character(id_column[[1L]]))
        if (!nzchar(id_column) || identical(
            toupper(id_column),
            "NULL"
        )) {
            id_column <- NULL
        }
    }

    available_columns <- colnames(probe_exclusion)
    standard_columns <- c("ProbeID", "TargetID", "IlmnID", "Name")

    if (!is.null(id_column)) {
        if (!(id_column %in% available_columns)) {
            stop("probeExclusionIdColumn not found in probe-exclusion file: ",
                id_column, ". Available columns: ",
                    paste(ifelse(nzchar(available_columns),
                    available_columns, "<blank>"
                ), collapse = ", "),
                call. = FALSE
            )
        }
        selected_column <- id_column
        selected_index <- match(id_column, available_columns)
        source <- "configured"
    } else {
        standard_match <- standard_columns[standard_columns %in%
            available_columns]
        if (length(standard_match) > 0L) {
            selected_column <- standard_match[[1L]]
            selected_index <- match(selected_column, available_columns)
            source <- "standard"
        } else {
            first_column <- available_columns[[1L]]
            first_values <- trimws(as.character(probe_exclusion[[1L]]))
            first_values <- first_values[!is.na(first_values) &
                nzchar(first_values)]
            first_column_is_unlabeled <- !nzchar(first_column)
            first_column_is_probe_like <- length(first_values) >
                0L && all(grepl("^(cg|ch\\.)", first_values))

            if (isTRUE(first_column_is_unlabeled) ||
                isTRUE(first_column_is_probe_like)) {
                selected_column <- first_column
                selected_index <- 1L
                source <- if (isTRUE(first_column_is_unlabeled)) {
                    "unlabeled first column"
                } else {
                    "probe-like first column"
                }
            } else {
                stop("Could not detect a probe-exclusion ID column in ",
                    probeExclusionPath,
                        ". Set probeExclusionIdColumn to one of: ",
                    paste(ifelse(nzchar(available_columns), available_columns,
                        "<blank>"
                    ), collapse = ", "),
                    call. = FALSE
                )
            }
        }
    }

    ids <- trimws(as.character(probe_exclusion[[selected_index]]))
    ids <- unique(ids[!is.na(ids) & nzchar(ids)])

    if (length(ids) == 0L) {
        stop("Probe-exclusion ID column contains no usable IDs: ",
            ifelse(nzchar(selected_column), selected_column,
                "<blank>"
            ),
            call. = FALSE
        )
    }

    overlap <- NA_integer_
    if (!is.null(featureNames)) {
        overlap <- sum(ids %in% featureNames)
    }

    structure(list(
        ids = ids, path = probeExclusionPath, column = selected_column,
        source = source, nRows = nrow(probe_exclusion), nIds = length(ids),
        overlap = overlap
    ), class = "dnaEPICO_probeExclusion_file_ids")
}

normalizeProbeExclusionPathsMinfiEwasWater <- function(probeExclusionPath) {
    if (length(probeExclusionPath) == 0L || all(is.na(probeExclusionPath))) {
        return(character(0))
    }

    paths <- splitOptionMinfiEwasWater(probeExclusionPath, sep = ";")
    paths[toupper(paths) != "NULL"]
}

normalizeProbeExclusionColumnsMinfiEwasWater <- function(
    probeExclusionIdColumn,
    nFiles
) {
    if (length(probeExclusionIdColumn) == 0L ||
        all(is.na(probeExclusionIdColumn))) {
        return(rep(list(NULL), nFiles))
    }

    columns <- splitOptionMinfiEwasWater(probeExclusionIdColumn,
        sep = ";"
    )
    if (length(columns) == 0L || identical(
        toupper(columns[[1L]]),
        "NULL"
    )) {
        return(rep(list(NULL), nFiles))
    }

    if (length(columns) == 1L) {
        return(rep(as.list(columns), nFiles))
    }

    if (length(columns) != nFiles) {
        stop("probeExclusionIdColumn must be NULL, a single column name, or one ",
            "semicolon-separated column name per probe-exclusion file.",
            call. = FALSE
        )
    }

    as.list(columns)
}

#' Read probe-exclusion identifiers from one or more CSV files
#'
#' @param probeExclusionPath Character vector or semicolon-separated string of
#'   CSV files containing probe identifiers to exclude.
#' @param probeExclusionIdColumn Character, `NULL`, or semicolon-separated
#'   string of column names.
#' @param featureNames Character vector or `NULL`. Optional array feature names
#'   used to count overlap.
#'
#' @return A list with class `'dnaEPICO_probeExclusion_ids'`.
#'
#' @keywords internal
#' @noRd
readProbeExclusionIdsMinfiEwasWater <- function(
    probeExclusionPath,
    probeExclusionIdColumn = NULL, featureNames = NULL
) {
    paths <- normalizeProbeExclusionPathsMinfiEwasWater(probeExclusionPath)
    columns <- normalizeProbeExclusionColumnsMinfiEwasWater(
        probeExclusionIdColumn = probeExclusionIdColumn,
        nFiles = length(paths)
    )

    if (length(paths) == 0L) {
        return(structure(
            list(
                ids = character(0), files = data.frame(
                    path = character(0),
                    column = character(0), source = character(0),
                        nRows = integer(0),
                    nIds = integer(0), overlap = integer(0),
                        stringsAsFactors = FALSE
                ),
                nIds = 0L,
                    overlap = if (is.null(featureNames)) NA_integer_ else 0L
            ),
            class = "dnaEPICO_probeExclusion_ids"
        ))
    }

    file_results <- Map(f = function(path, column) {
        readProbeExclusionFileIdsMinfiEwasWater(
            probeExclusionPath = path,
            probeExclusionIdColumn = column, featureNames = featureNames
        )
    }, path = paths, column = columns)

    ids <- unique(unlist(lapply(file_results, `[[`, "ids"), use.names = FALSE))
    overlap <- NA_integer_
    if (!is.null(featureNames)) {
        overlap <- sum(ids %in% featureNames)
    }

    files <- do.call(rbind, lapply(file_results, function(result) {
        data.frame(
            path = result$path, column = ifelse(nzchar(result$column),
                result$column, "<blank>"
            ), source = result$source,
            nRows = result$nRows, nIds = result$nIds, overlap = result$overlap,
            stringsAsFactors = FALSE
        )
    }))

    structure(list(
        ids = ids, files = files, nIds = length(ids),
        overlap = overlap
    ), class = "dnaEPICO_probeExclusion_ids")
}

normalizeEpicV2ManifestFlagsMinfiEwasWater <- function(epicV2ManifestFlags) {
    default_flags <- c(
        CH_WGBS_evidence = TRUE, CH_BLAT = TRUE,
        MissingPos = TRUE, MismatchPos = FALSE
    )

    if (length(epicV2ManifestFlags) == 0L || all(is.na(epicV2ManifestFlags))) {
        return(default_flags)
    }

    if (is.logical(epicV2ManifestFlags)) {
        if (is.null(names(epicV2ManifestFlags))) {
            if (length(epicV2ManifestFlags) != length(default_flags)) {
                stop("Unnamed epicV2ManifestFlags must have ",
                    length(default_flags), " values.",
                    call. = FALSE
                )
            }
            names(epicV2ManifestFlags) <- names(default_flags)
        }

        unknown_flags <- setdiff(
            names(epicV2ManifestFlags),
            names(default_flags)
        )
        if (length(unknown_flags) > 0L) {
            stop("Unknown EPICv2 manifest flag(s): ", paste(unknown_flags,
                collapse = ", "
            ), call. = FALSE)
        }

        if (anyNA(epicV2ManifestFlags)) {
            stop("epicV2ManifestFlags must contain only TRUE or FALSE values.",
                call. = FALSE
            )
        }

        default_flags[names(epicV2ManifestFlags)] <- epicV2ManifestFlags
        return(default_flags)
    }

    flag_values <- splitOptionMinfiEwasWater(epicV2ManifestFlags,
        sep = ";"
    )
    if (length(flag_values) == 0L) {
        return(default_flags)
    }

    flag_names <- sub("=.*$", "", flag_values)
    flag_enabled <- rep(TRUE, length(flag_values))
    has_equals <- grepl("=", flag_values, fixed = TRUE)
    flag_enabled[has_equals] <- as.logical(toupper(sub(
        "^.*=",
        "", flag_values[has_equals]
    )))
    if (anyNA(flag_enabled)) {
        stop("epicV2ManifestFlags must contain only TRUE or FALSE values.",
            call. = FALSE
        )
    }

    names(flag_enabled) <- flag_names
    normalizeEpicV2ManifestFlagsMinfiEwasWater(flag_enabled)
}

extractEpicV2ManifestExclusionIdsMinfiEwasWater <- function(
    manifest,
    epicV2ManifestFlags = c(
        CH_WGBS_evidence = TRUE, CH_BLAT = TRUE,
        MissingPos = TRUE, MismatchPos = FALSE
    ), featureNames = NULL
) {
    flags <- normalizeEpicV2ManifestFlagsMinfiEwasWater(epicV2ManifestFlags)
    enabled_flags <- names(flags)[!is.na(flags) & flags]

    if (length(enabled_flags) == 0L) {
        return(structure(
            list(
                ids = character(0), flags = flags,
                flagCounts = stats::setNames(integer(0), character(0)),
                nIds = 0L,
                    overlap = if (is.null(featureNames)) NA_integer_ else 0L
            ),
            class = "dnaEPICO_epicV2Manifest_ids"
        ))
    }

    manifest <- as.data.frame(manifest, stringsAsFactors = FALSE)
    missing_flags <- setdiff(enabled_flags, colnames(manifest))
    if (length(missing_flags) > 0L) {
        stop("EPICv2 manifest is missing requested flag column(s): ",
            paste(missing_flags, collapse = ", "),
            call. = FALSE
        )
    }

    ids <- rownames(manifest)
    default_row_names <- identical(ids, as.character(seq_len(nrow(manifest))))
    if (is.null(ids) || all(!nzchar(ids)) || isTRUE(default_row_names)) {
        if (!("IlmnID" %in% colnames(manifest))) {
            stop("EPICv2 manifest must have IlmnID row names or an IlmnID column.",
                call. = FALSE
            )
        }
        ids <- as.character(manifest$IlmnID)
    }

    selected <- rep(FALSE, nrow(manifest))
    flag_counts <- stats::setNames(
        integer(length(enabled_flags)),
        enabled_flags
    )
    for (flag in enabled_flags) {
        flag_values <- toupper(trimws(as.character(manifest[[flag]])))
        flag_selected <- !is.na(flag_values) & flag_values %in%
            c("Y", "TRUE", "1")
        flag_counts[[flag]] <- sum(flag_selected)
        selected <- selected | flag_selected
    }

    selected_ids <- unique(ids[selected & !is.na(ids) & nzchar(ids)])
    overlap <- NA_integer_
    if (!is.null(featureNames)) {
        overlap <- sum(selected_ids %in% featureNames)
    }

    structure(list(
        ids = selected_ids, flags = flags, flagCounts = flag_counts,
        nIds = length(selected_ids), overlap = overlap
    ), class = "dnaEPICO_epicV2Manifest_ids")
}

readEpicV2ManifestExclusionIdsMinfiEwasWater <- function(
    useEpicV2Manifest = FALSE,
    epicV2ManifestFlags = c(
        CH_WGBS_evidence = TRUE, CH_BLAT = TRUE,
        MissingPos = TRUE, MismatchPos = FALSE
    ), featureNames = NULL
) {
    if (!isTRUE(useEpicV2Manifest)) {
        return(extractEpicV2ManifestExclusionIdsMinfiEwasWater(
            manifest = data.frame(row.names = character(0)),
            epicV2ManifestFlags = stats::setNames(rep(
                FALSE,
                4L
            ), c(
                "CH_WGBS_evidence", "CH_BLAT", "MissingPos",
                "MismatchPos"
            )), featureNames = featureNames
        ))
    }

    if (!requireNamespace("AnnotationHub", quietly = TRUE)) {
        stop("AnnotationHub is required when useEpicV2Manifest = TRUE. ",
            "Install AnnotationHub or set useEpicV2Manifest = FALSE.",
            call. = FALSE
        )
    }

    hub <- AnnotationHub::AnnotationHub()
    manifest <- hub[["AH116484"]]

    extractEpicV2ManifestExclusionIdsMinfiEwasWater(
        manifest = manifest,
        epicV2ManifestFlags = epicV2ManifestFlags, featureNames = featureNames
    )
}

#' Filter probes from a normalized methylation object
#'
#' Apply detection P-value, chromosome, SNP, and probe-exclusion filters to the
#' primary normalized object and return the filtered result.
#'
#' @param normData Object returned by `normalizeMinfiEwasWater()`.
#' @param RGSet Filtered `RGChannelSet` aligned with `normData`.
#' @param pvalThreshold Numeric. Probes must have detection P values below this
#'   threshold in all samples to be retained.
#' @param chrToRemove Character vector or comma-separated string of chromosome
#'   names to remove, for example `'chrX,chrY'`.
#' @param snpsToRemove Character vector or comma-separated string of SNP probe
#'   types to remove, for example `'SBE,CpG'`.
#' @param mafThreshold Numeric. Minor allele frequency threshold passed to
#'   `minfi::dropLociWithSnps()`.
#' @param probeExclusionPath Character vector or semicolon-separated string of
#'   CSV files containing probe IDs to remove.
#' @param probeExclusionIdColumn Character or `NULL`. Column containing probe
#'   IDs. When `NULL` or `''`, each file is auto-detected using `ProbeID`,
#'   `TargetID`, `IlmnID`, or `Name`, then falling back to an unlabeled or
#'   probe-like first column.
#' @param useEpicV2Manifest Logical. If `TRUE`, also remove EPICv2 probes
#'   flagged in the Peters et al. expanded manifest from AnnotationHub resource
#'   `AH116484`.
#' @param epicV2ManifestFlags Named logical vector controlling which EPICv2
#'   manifest flags are removed. Defaults remove `CH_WGBS_evidence`, `CH_BLAT`,
#'   and `MissingPos`, but not `MismatchPos`.
#' @param crossReactivePath Deprecated alias for `probeExclusionPath`.
#' @param crossReactiveIdColumn Deprecated alias for `probeExclusionIdColumn`.
#' @param detPtype Character. Detection P-value mode passed to
#'   `minfi::detectionP()` for the probe filter. Common values in minfi
#'   workflows are `'m+u'` and `'negative'`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `'dnaEPICO_minfiEwasWater_filter'` containing the
#'   filtered object and counts for each filtering stage.
#'
#' @examplesIf requireNamespace("minfiData", quietly = TRUE)
#' ex <- dnaEPICO:::exampleMinfiWorkflowStateDnaEpico()
#' filtered_data <- filterProbesMinfiEwasWater(
#'     normData = ex$normData,
#'     RGSet = ex$sampleData$RGSet,
#'     pvalThreshold = 1,
#'     chrToRemove = "chrY",
#'     snpsToRemove = "SBE",
#'     mafThreshold = 1,
#'     probeExclusionPath = ex$probeExclusionPath,
#'     detPtype = "m+u",
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' filtered_data$counts[["probeExclusion"]]
#'
#' @export
filterProbesMinfiEwasWater <- function(
    normData, RGSet, pvalThreshold = 0.01,
    chrToRemove = "chrX,chrY", snpsToRemove = "SBE,CpG", mafThreshold = 0.1,
    probeExclusionPath, probeExclusionIdColumn = NULL,
        useEpicV2Manifest = FALSE,
    epicV2ManifestFlags = c(
        CH_WGBS_evidence = TRUE, CH_BLAT = TRUE,
        MissingPos = TRUE, MismatchPos = FALSE
    ), detPtype = "m+u",
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_filterProbesMinfiEwasWater.txt",
    crossReactivePath = NULL, crossReactiveIdColumn = NULL
) {
    if (!is.null(crossReactivePath)) {
        probeExclusionPath <- crossReactivePath
    }
    if (!is.null(crossReactiveIdColumn)) {
        probeExclusionIdColumn <- crossReactiveIdColumn
    }

    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    chr_list <- splitOptionMinfiEwasWater(chrToRemove, sep = ",")
    snp_list <- splitOptionMinfiEwasWater(snpsToRemove, sep = ",")

    pvalThreshold <- validateProbabilityDnaEpico(
        pvalThreshold,
        "pvalThreshold"
    )
    mafThreshold <- validateProbabilityDnaEpico(
        mafThreshold,
        "mafThreshold"
    )
    useEpicV2Manifest <- validateLogicalScalarDnaEpico(
        useEpicV2Manifest,
        "useEpicV2Manifest"
    )

    detP <- minfi::detectionP(RGSet, type = detPtype)
    nan_detp_converted <- sum(is.nan(detP))
    if (nan_detp_converted > 0L) {
        detP[is.nan(detP)] <- NA_real_
    }

    normalized_samples <- colnames(normData$primary)
    rgset_samples <- colnames(RGSet)
    matchSampleIdentifiersDnaEpico(
        query = normalized_samples,
        reference = rgset_samples, queryLabel = "Normalized sample identifiers",
        referenceLabel = "RGSet sample identifiers", requireSameSet = TRUE
    )

    feature_match <- match(
        Biobase::featureNames(normData$primary),
        rownames(detP)
    )
    if (anyNA(feature_match)) {
        stop("Detection P values are missing for ", sum(is.na(feature_match)),
            " normalized probes.",
            call. = FALSE
        )
    }
    detP <- detP[feature_match, normalized_samples, drop = FALSE]

    keep_detp <- rowSums(!is.na(detP) & detP < pvalThreshold) ==
        ncol(detP)
    filtered_detp <- normData$primary[keep_detp, ]

    ann <- minfi::getAnnotation(RGSet)
    remove_probes <- ann$Name[ann$chr %in% chr_list]
    keep_chr <- !(Biobase::featureNames(filtered_detp) %in% remove_probes)
    filtered_chr <- filtered_detp[keep_chr, ]

    filtered_snp <- minfi::dropLociWithSnps(filtered_chr,
        snps = snp_list,
        maf = mafThreshold
    )

    probe_file_ids <- readProbeExclusionIdsMinfiEwasWater(
        probeExclusionPath = probeExclusionPath,
        probeExclusionIdColumn = probeExclusionIdColumn,
            featureNames = Biobase::featureNames(filtered_snp)
    )
    epicv2_manifest_ids <- readEpicV2ManifestExclusionIdsMinfiEwasWater(
        useEpicV2Manifest = useEpicV2Manifest,
        epicV2ManifestFlags = epicV2ManifestFlags,
            featureNames = Biobase::featureNames(filtered_snp)
    )

    probe_exclusion_ids <- unique(c(probe_file_ids$ids,
        epicv2_manifest_ids$ids))
    probe_exclusion_overlap <- NA_integer_
    if (length(probe_exclusion_ids) > 0L) {
        probe_exclusion_overlap <- sum(probe_exclusion_ids %in%
            Biobase::featureNames(filtered_snp))
    }

    if (length(probe_exclusion_ids) > 0L &&
        isTRUE(!is.na(probe_exclusion_overlap)) &&
        identical(probe_exclusion_overlap, 0L)) {
        warning("No probe-exclusion IDs overlap the filtered array feature names. ",
            "Check that probeExclusionPath, probeExclusionIdColumn, and manifest ",
            "settings match the array platform.",
            call. = FALSE
        )
    }

    keep_probe_exclusion <- !(Biobase::featureNames(filtered_snp) %in%
        probe_exclusion_ids)
    filtered_final <- filtered_snp[keep_probe_exclusion, ]

    counts <- c(
        start = nrow(normData$primary), detP = nrow(filtered_detp),
        chromosome = nrow(filtered_chr), snp = nrow(filtered_snp),
        probeExclusion = nrow(filtered_final)
    )
    if (nrow(filtered_final) == 0L) {
        stop("Probe filtering removed every probe. Review the detection P-value, chromosome, SNP, and probe-exclusion settings.",
            call. = FALSE
        )
    }

    probe_file_lines <- if (nrow(probe_file_ids$files) > 0L) {
        apply(probe_file_ids$files, 1L, function(row) {
            paste0(
                "  - ", row[["path"]], " [", row[["column"]],
                "; ", row[["source"]], "; IDs=", row[["nIds"]],
                "; overlap=", row[["overlap"]], "]"
            )
        })
    } else {
        "  - none"
    }
    epicv2_flag_lines <- if (isTRUE(useEpicV2Manifest)) {
        paste0(
            "  - ", names(epicv2_manifest_ids$flagCounts),
            ": ", epicv2_manifest_ids$flagCounts
        )
    } else {
        "  - disabled"
    }

    emitLogMinfiEwasWater(
        c(
            paste(
                "Probe filter threshold:    ",
                pvalThreshold
            ), paste(
                "Detection P NaN to NA:     ",
                nan_detp_converted
            ), paste(
                "Chromosomes removed:       ",
                paste(chr_list, collapse = ", ")
            ), paste(
                "SNP filters removed:       ",
                paste(snp_list, collapse = ", ")
            ), paste(
                "MAF threshold:             ",
                mafThreshold
            ), "Probe-exclusion files:", probe_file_lines,
            paste("Probe-exclusion file IDs: ", probe_file_ids$nIds),
            paste("EPICv2 manifest enabled:  ", useEpicV2Manifest),
            "EPICv2 manifest flag counts:", epicv2_flag_lines, paste(
                "EPICv2 manifest IDs:      ",
                epicv2_manifest_ids$nIds
            ), paste(
                "Probe-exclusion IDs total:",
                length(probe_exclusion_ids)
            ), paste(
                "Probe-exclusion overlap:  ",
                probe_exclusion_overlap
            ), paste(
                "Probes after detP filter:  ",
                counts[["detP"]]
            ), paste(
                "Probes after chr filter:   ",
                counts[["chromosome"]]
            ), paste(
                "Probes after SNP filter:   ",
                counts[["snp"]]
            ), paste(
                "Probes after exclusion:    ",
                counts[["probeExclusion"]]
            ), "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(
        filtered = filtered_final, detPFiltered = filtered_detp,
        chrFiltered = filtered_chr, snpFiltered = filtered_snp,
        detP = detP, counts = counts, pvalThreshold = pvalThreshold,
        nanDetPConverted = nan_detp_converted, chrToRemove = chr_list,
        snpsToRemove = snp_list, mafThreshold = mafThreshold,
        probeExclusionPath = probeExclusionPath,
            probeExclusionIdColumn = probeExclusionIdColumn,
        probeExclusionIds = probe_exclusion_ids,
            probeExclusionFileIds = probe_file_ids,
        epicV2ManifestIds = epicv2_manifest_ids,
            useEpicV2Manifest = useEpicV2Manifest,
        epicV2ManifestFlags = normalizeEpicV2ManifestFlagsMinfiEwasWater(epicV2ManifestFlags),
        crossReactivePath = probeExclusionPath,
            crossReactiveIdColumn = probeExclusionIdColumn,
        crossReactiveIds = probe_file_ids
    ), class = "dnaEPICO_minfiEwasWater_filter")
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
#' @return A list with class `'dnaEPICO_minfiEwasWater_metrics'` containing
#'   `beta`, `m`, and `cn`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMinfiMetricsStateDnaEpico()
#' metrics_data <- extractMetricsMinfiEwasWater(
#'     filteredData = ex$filteredData,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' names(metrics_data)
#'
#' @export
extractMetricsMinfiEwasWater <- function(
    filteredData, verbose = FALSE,
    logs = FALSE, log_dir = NULL,
        log_file = "log_extractMetricsMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    beta <- minfi::getBeta(filteredData$filtered)
    m <- minfi::getM(filteredData$filtered)
    cn <- minfi::getCN(filteredData$filtered)

    if (is.null(rownames(beta)) || is.null(colnames(beta))) {
        stop("The extracted beta matrix must have probe and sample identifiers.",
            call. = FALSE
        )
    }
    validateMethylationProbeIdentifiersDnaEpico(
        rownames(beta),
        "Extracted methylation probe identifiers"
    )
    validateSampleIdentifiersDnaEpico(colnames(beta),
        "Extracted methylation sample identifiers")
    metrics_match <- vapply(list(m = m, cn = cn), function(metric) {
        identical(dim(metric), dim(beta)) && identical(
            rownames(metric),
            rownames(beta)
        ) && identical(colnames(metric), colnames(beta))
    }, logical(1))
    if (!all(metrics_match)) {
        stop("Extracted Beta, M, and CN matrices must have identical probe and sample order.",
            call. = FALSE
        )
    }

    validations <- list(
        beta = inspectMethylationMatrixDnaEpico(
            beta,
            "beta"
        ), m = inspectMethylationMatrixDnaEpico(m, "m"),
        cn = inspectMethylationMatrixDnaEpico(cn, "cn")
    )
    beta <- validations$beta$values
    m <- validations$m$values
    cn <- validations$cn$values

    converted_nan <- sum(vapply(
        validations, function(validation) validation$boundaries$NaN.Converted,
        numeric(1)
    ))
    invalid_cpgs <- sum(vapply(
        validations, function(validation) validation$boundaries$Invalid.CpGs,
        numeric(1)
    ))
    if (converted_nan > 0L) {
        warning(converted_nan,
            " NaN values returned by methylation accessors were converted to NA.",
            call. = FALSE
        )
    }
    if (invalid_cpgs > 0L) {
        warning(invalid_cpgs,
            " invalid scale-specific CpG records were retained for reporting.",
            call. = FALSE
        )
    }

    preview_cols <- seq_len(min(ncol(beta), 5L))

    emitLogMinfiEwasWater(
        c(
            "Extracting final DNAm matrices (M, Beta, CN).",
            unlist(lapply(validations, function(validation) {
                formatMethylationBoundariesLogDnaEpico(validation$boundaries)
            }), use.names = FALSE), "Preview of beta values:",
                previewLinesMinfiEwasWater(beta[,
                preview_cols,
                drop = FALSE
            ]), "Preview of M-values:",
            previewLinesMinfiEwasWater(m[, preview_cols, drop = FALSE]),
            "Preview of copy-number values:", previewLinesMinfiEwasWater(cn[,
                preview_cols,
                drop = FALSE
            ]), "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(beta = beta, m = m, cn = cn, validation = validations),
        class = "dnaEPICO_minfiEwasWater_metrics"
    )
}

#' Plot multidimensional scaling or density summaries from final metrics
#'
#' @param metricsData Object returned by `extractMetricsMinfiEwasWater()`.
#' @param targets Filtered phenotype data aligned with `metricsData`.
#' @param plot Character. Plot type: `'mds'` or `'density'`.
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
#' @examples
#' ex <- dnaEPICO:::exampleMinfiMetricsStateDnaEpico()
#' plotMetricsMinfiEwasWater(
#'     metricsData = ex$metricsData,
#'     targets = ex$targets,
#'     plot = "density",
#'     plotGroupVar = "Sex",
#'     sexColumn = "Sex",
#'     display = FALSE,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#'
#' @export
plotMetricsMinfiEwasWater <- function(
    metricsData, targets, plot = c(
        "mds",
        "density"
    ), plotGroupVar = "Sex", sexColumn = "Sex", display = FALSE,
    file = NULL, width = 2000L, height = 1000L, res = 150L, verbose = FALSE,
    logs = FALSE, log_dir = NULL, log_file = "log_plotMetricsMinfiEwasWater.txt"
) {
    plot <- match.arg(plot)

    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    draw_fun <- switch(plot,
        mds = function() {
            group_factor <- factor(targets[[plotGroupVar]])
            group_sex <- factor(targets[[sexColumn]])
            pal <- RColorBrewer::brewer.pal(8, "Dark2")

            graphics::par(mfrow = c(1, 2))
            limma::plotMDS(metricsData$m,
                main = plotGroupVar, top = 1000,
                gene.selection = "common", col = pal[group_factor],
                dim = c(1, 2)
            )
            graphics::legend("right",
                legend = levels(group_factor),
                text.col = pal, cex = 0.7, bg = "white"
            )

            limma::plotMDS(metricsData$m,
                main = sexColumn, top = 1000,
                gene.selection = "common", col = pal[group_sex],
                dim = c(2, 3)
            )
            graphics::legend("topright",
                legend = levels(group_sex),
                text.col = pal, cex = 0.7, bg = "white"
            )
        },
        density = function() {
            group_factor <- factor(targets[[plotGroupVar]])
            graphics::par(mfrow = c(1, 2))
            minfi::densityPlot(metricsData$beta,
                sampGroups = group_factor,
                main = "Beta values", legend = FALSE, xlab = "Beta values"
            )
            graphics::legend("top",
                legend = levels(group_factor),
                text.col = RColorBrewer::brewer.pal(8, "Dark2")
            )

            minfi::densityPlot(metricsData$m,
                sampGroups = group_factor,
                main = "M-values", legend = FALSE, xlab = "M values"
            )
            graphics::legend("topleft",
                legend = levels(group_factor),
                text.col = RColorBrewer::brewer.pal(8, "Dark2")
            )
        }
    )

    runPlotMinfiEwasWater(
        draw_fun = draw_fun, display = display,
        file = file, width = width, height = height, res = res
    )

    emitLogMinfiEwasWater(
        c(paste(
            "Metrics plot created:     ",
            plot
        ), if (is.null(file)) {
            "Metrics plot file:        none"
        } else {
            paste("Metrics plot file:        ", file)
        }, "============================================================"),
        verbose = verbose, log_path = log_path
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
#'     RGSet = ex$RGSet,
#'     output_dir = output_dir,
#'     verbose = FALSE,
#'     logs = FALSE
#' )
#' dir.exists(output_dir)
#'
#' @export
plotCtrlMinfiEwasWater <- function(
    RGSet, output_dir = NULL,
    verbose = FALSE, logs = FALSE, log_dir = NULL,
        log_file = "log_plotCtrlMinfiEwasWater.txt"
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = log_dir,
        log_file = log_file
    )

    if (is.null(output_dir)) {
        emitLogMinfiEwasWater(
            c(
                "ENmix control plots skipped because output_dir is NULL.",
                "============================================================"
            ),
            verbose = verbose, log_path = log_path
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
        c(paste(
            "Generated ENmix control JPGs in:",
            output_dir
        ), "============================================================"),
        verbose = verbose, log_path = log_path
    )

    invisible(output_dir)
}
