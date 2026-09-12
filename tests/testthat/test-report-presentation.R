test_that("figure exports are real PNG, JPEG and TIFF files with matching dimensions", {
    skip_if_not_installed("magick")
    tmp <- withr::local_tempdir()
    source <- file.path(tmp, "source.png")
    magick::image_write(magick::image_blank(120, 80, "transparent"), source)
    dir.create(file.path(tmp, "exports"))
    result <- dnaEPICO:::.copyReportFigureFileDnaEpico(source, "figure",
        file.path(tmp, "exports"), TRUE)
    expect_true(all(file.exists(result$formats)))
    info <- lapply(result$formats, function(path) magick::image_info(magick::image_read(path)))
    expect_identical(vapply(info, function(x) x$format, character(1)),
        c(png = "PNG", jpeg = "JPEG", tiff = "TIFF"))
    expect_true(all(vapply(info, function(x) x$width == 120 && x$height == 80, logical(1))))
    expect_false(info$jpeg$matte)
    for (format in c("png", "jpeg")) {
        path <- result$formats[[format]]
        payload <- readLines(paste0(path, ".download.js"))
        encoded <- sub('^document.currentScript.dataset.downloadData = "(.*)";$',
            '\\1', payload[2])
        expect_identical(base64enc::base64decode(encoded),
            readBin(path, "raw", n = file.info(path)$size))
        expect_match(payload[1], paste0('"image/', format, '"'), fixed = TRUE)
    }
    expect_false(file.exists(paste0(result$formats[["tiff"]], ".download.js")))
})

test_that("original PNG and JPEG downloads work without image conversion", {
    tmp <- withr::local_tempdir()
    for (extension in c("png", "jpg", "jpeg")) {
        source <- file.path(tmp, paste0("source.", extension))
        bytes <- as.raw(0:255)
        writeBin(bytes, source)
        result <- dnaEPICO:::.copyReportFigureFileDnaEpico(source, "export", tmp, FALSE)
        payload <- readLines(paste0(result$download, ".download.js"))
        encoded <- sub('^document.currentScript.dataset.downloadData = "(.*)";$',
            '\\1', payload[2])
        expect_identical(base64enc::base64decode(encoded), bytes)
    }
})

test_that("identifier distributions are omitted but observation summaries remain", {
    items <- data.frame(original_name = c("distribution_subject_key_categorical.tiff",
        "participantObservationCount_subject_key.tiff", "distribution_age_continuous.tiff"))
    kept <- dnaEPICO:::.excludeParticipantFiguresDnaEpico(items, "subject_key")
    expect_equal(nrow(kept), 2L)
    expect_match(kept$original_name[1], "participantObservationCount")
    prepared <- list(data = data.frame(subject_key = factor(c("A", "A", "B")), age = 1:3),
        personVar = "subject_key", phenotypes = "age", factorVars = "subject_key",
        covariates = "subject_key")
    plots <- dnaEPICO::plotMethylationGLMDistributions(prepared)
    expect_false("subject_key" %in% names(plots$factors))
    expect_true("age" %in% names(plots$phenotypes))
})

test_that("participant heatmaps retain counts without overlapping cell labels", {
    p <- dnaEPICO:::categoricalModelAssociationPlotDnaEpico(
        rep(c("A", "B"), each = 80), rep(paste0("P", 1:80), each = 2),
        "Profession", "Participant", identifier = TRUE)
    expect_equal(sum(p$data$Freq), 160)
    expect_false(any(vapply(p$layers, function(x) inherits(x$geom, "GeomText"), logical(1))))
    expect_equal(nlevels(p$data$panel), 2)
    short <- dnaEPICO:::createModelAssociationPlotDnaEpico(data.frame(
        Profession = factor(c("A", NA, "B")), subject_key = c(1, 1, 2)),
        "Profession", "subject_key", identifierColumns = "subject_key")
    expect_equal(sum(short$data$Freq), 3)
    expect_true("Missing" %in% short$data$phenotype)
    expect_false(any(vapply(short$layers, function(x) inherits(x$geom, "GeomText"), logical(1))))
})

test_that("array-sized residual diagnostics retain every CpG without a count legend", {
    data <- data.frame(meanMethylation = seq(0, 1, length.out = 10001), diagnosticY = 0.03)
    p <- dnaEPICO:::createDiagnosticMeanPlotDnaEpico(data, list(label = "Average beta"), "Residual SD")
    built <- ggplot2::ggplot_build(p)
    expect_equal(nrow(built$data[[1]]), nrow(data))
    expect_equal(p$theme$legend.position, "none")
    expect_true(all(is.finite(built$data[[2]]$y)))
})

test_that("detection bars use sample means in input order", {
    p <- dnaEPICO:::sampleDetectionBarPlotDnaEpico(list(
        meanDetP = c(S2 = 0.03, S1 = 0.001, S3 = 0.01), detPThreshold = 0.05))
    data <- ggplot2::ggplot_build(p)$data[[1]]
    expect_equal(data$y, c(0.03, 0.001, 0.01))
    expect_null(p$labels$caption)
})

test_that("all numeric and named timepoints have individual highlights", {
    env <- new.env(parent = asNamespace("dnaEPICO"))
    env$html_escape <- dnaEPICO:::.htmlEscapeDnamReport
    dnaEPICO:::.installDrHelpers(env)
    for (values in list(c("1", "2", "3"), c("Baseline", "Follow-up"))) {
        text <- env$collapse_values(values, highlight = TRUE)
        html <- env$report_inline_markup(paste0("Timepoint values ", text, "."))
        for (value in values) expect_match(html,
            paste0('class="dnaepico-data-value">', value, '</span>'), fixed = TRUE)
    }
})
