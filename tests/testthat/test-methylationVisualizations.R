test_that("model distribution and diagnostic plots omit embedded descriptions", {
    histogram <- dnaEPICO:::createDistributionPlotMethylationGLM(
        c(1, 2, 2, 3, NA), "age", type = "hist"
    )
    bars <- dnaEPICO:::createDistributionPlotMethylationGLM(
        c("Case", "Control", "Case", NA), "status", type = "bar"
    )

    expect_s3_class(histogram, "ggplot")
    expect_s3_class(bars, "ggplot")
    expect_null(histogram$labels$title)
    expect_null(bars$labels$title)
    expect_null(histogram$labels$caption)
    expect_null(bars$labels$caption)
    expect_true(any(vapply(histogram$layers, function(layer) {
        inherits(layer$geom, "GeomPoint")
    }, logical(1))))
    expect_true(any(vapply(histogram$layers, function(layer) {
        inherits(layer$geom, "GeomBoxplot")
    }, logical(1))))
    expect_true(any(vapply(bars$layers, function(layer) {
        inherits(layer$geom, "GeomText")
    }, logical(1))))

    cpgs <- paste0("cg", sprintf("%08d", seq_len(20)))
    diagnostics <- dnaEPICO:::buildMethylationTermDiagnosticsDnaEpico(
        summaryData = data.frame(
            CpG = cpgs,
            Term = rep("age", length(cpgs)),
            P.value = seq(0.001, 0.2, length.out = length(cpgs)),
            FDR = seq(0.01, 0.4, length.out = length(cpgs)),
            Estimate = seq(-0.2, 0.2, length.out = length(cpgs)),
            Std.Error = rep(0.03, length(cpgs)),
            ResidualSD = seq(0.1, 0.3, length.out = length(cpgs)),
            check.names = FALSE
        ),
        phenotype = "age", term = "age", termColumn = "Term",
        pValueColumn = "P.value", yColumn = "ResidualSD",
        yLabel = "Residual standard deviation",
        diagnosticMean = list(
            values = stats::setNames(seq(0.2, 0.8, length.out = 20), cpgs),
            label = "Mean beta value"
        ),
        fdrThreshold = 0.05, estimateColumn = "Estimate",
        standardErrorColumn = "Std.Error"
    )

    expect_null(diagnostics$plots$qqplot$labels$title)
    expect_null(diagnostics$plots$qqplot$labels$caption)
    expect_true(is.na(diagnostics$lambda))
    expect_true(any(vapply(diagnostics$plots$qqplot$layers, function(layer) {
        inherits(layer$geom, "GeomRibbon")
    }, logical(1))))
    expect_s3_class(diagnostics$plots$volcano, "ggplot")
    expect_s3_class(diagnostics$plots$effectForest, "ggplot")
    expect_identical(
        diagnostics$plots$effectForest$labels$y,
        "CpG (top 20 ranked by p-value)"
    )
    expect_null(diagnostics$plots$volcano$labels$title)
    label_layers <- Filter(function(layer) {
        inherits(layer$geom, "GeomTextRepel")
    }, diagnostics$plots$residualSignificance$layers)
    expect_length(label_layers, 1L)
    expect_equal(label_layers[[1L]]$geom_params$max.overlaps, Inf)
    expect_true(all(vapply(diagnostics$plots, function(plot) {
        is.null(plot) || is.null(plot$labels$caption)
    }, logical(1))))
})

test_that("four-set overlap has a ranked intersection companion", {
    sets <- list(
        Profession1 = c("cg1", "cg2", "cg5"),
        Profession2 = c("cg2", "cg3", "cg5"),
        Age = c("cg1", "cg4", "cg5"),
        Sex = c("cg2", "cg4", "cg5")
    )
    genes <- stats::setNames(paste0("GENE", 1:5), paste0("cg", 1:5))

    plot <- dnaEPICO:::createIntersectionPlotDnaEpico(sets, genes)

    expect_s3_class(plot, "ggplot")
    expect_null(plot$labels$title)
    expect_null(plot$labels$caption)
    expect_identical(plot$labels$x, "CpGs")
    expect_identical(plot$labels$y, "Intersection membership")
    expect_true(any(vapply(plot$layers, function(layer) {
        inherits(layer$geom, "GeomCol")
    }, logical(1))))
})

test_that("LME visualisations include paired changes for two timepoints", {
    prepared <- list(
        data = data.frame(
            person = rep(paste0("P", 1:4), each = 2),
            Timepoint = rep(c("T1", "T2"), 4),
            score = c(1, 2, 2, 3, 4, 3, 5, 6),
            sex = rep(c("F", "M"), each = 4),
            stringsAsFactors = FALSE
        ), phenotypes = "score", covariates = "sex", factorVars = "sex",
        interactionTerm = NULL, timeVar = "Timepoint", personVar = "person"
    )

    result <- dnaEPICO:::plotMethylationLMEDistributions(
        prepared, outputDir = NULL, display = FALSE
    )

    expect_s3_class(result, "dnaEPICO_methylationLME_distribution_plots")
    expect_true("change_score" %in% names(result$plots))
    expect_null(result$plots$change_score$labels$title)
    expect_null(result$plots$change_score$labels$caption)
})

test_that("Manhattan plots use annotated raw p-value columns and versioned names", {
    tmp <- withr::local_tempdir()
    annotated <- data.frame(
        IlmnID = paste0("cg", sprintf("%08d", seq_len(44))),
        chr = rep(paste0("chr", 1:22), each = 2),
        pos = rep(c(1000, 2000), 22),
        `statusCaseP.Value` = c(0, 1e-9, rep(0.1, 42)),
        `statusCaseAdjusted.P.Value` = rep(0.2, 44),
        check.names = FALSE
    )

    expect_identical(
        dnaEPICO:::identifyAnnotatedPvalueColumnsDnaEpico(annotated),
        "statusCaseP.Value"
    )
    result <- dnaEPICO:::plotAnnotatedManhattanDnaEpico(
        annotatedResults = annotated, analysis = "GLM", outputDir = tmp,
        plotWidth = 600, plotHeight = 400, plotDPI = 72, display = FALSE
    )

    expect_s3_class(result, "dnaEPICO_manhattan_plots")
    expect_equal(dnaEPICO:::genomeWideThresholdDnaEpico(), 9e-8)
    expect_true(all(result$manifest$genomeWideThreshold == 9e-8))
    expect_setequal(
        result$manifest$filename,
        c(
            "manhattan_statusCaseP.Value_v1.tiff",
            "manhattan_statusCaseP.Value_v2.tiff"
        )
    )
    expect_true(all(file.exists(file.path(tmp, result$manifest$filename))))
    expect_false(any(grepl("model", result$manifest$filename,
        ignore.case = TRUE
    )))
    expect_true(all(vapply(result$plots[["statusCaseP.Value"]], function(plot) {
        is.null(plot$labels$title)
    }, logical(1))))
    expect_true(all(vapply(result$plots[["statusCaseP.Value"]], function(plot) {
        is.null(plot$labels$caption)
    }, logical(1))))
    expect_true(all(vapply(result$plots[["statusCaseP.Value"]], function(plot) {
        identical(plot$theme$legend.position, "none")
    }, logical(1))))
    expect_s3_class(
        result$plots[["statusCaseP.Value"]][["v1"]]$coordinates,
        "CoordPolar"
    )
    expect_true(any(vapply(
        result$plots[["statusCaseP.Value"]][["v1"]]$layers,
        function(layer) inherits(layer$geom, "GeomRect"), logical(1)
    )))
})

test_that("model-level Venn expands factor terms and keeps label mappings", {
    testthat::skip_if_not_installed("ggVennDiagram")

    tmp <- withr::local_tempdir()
    ids <- paste0("cg", sprintf("%08d", seq_len(6)))
    annotated <- data.frame(
        IlmnID = ids,
        Profession1P.Value = c(1e-9, 0.01, 0.2, 1e-6, 0.8, 0.6),
        Profession2P.Value = c(0.02, 1e-10, 0.3, 0.7, 1e-6, 0.9),
        AgeP.Value = c(0.03, 0.04, 1e-11, 0.8, 0.9, 1e-6),
        Profession_Omnibus_P.Value = c(1e-9, 0.02, 0.4, 0.8, 0.9, 0.7),
        Age_Omnibus_P.Value = c(0.03, 1e-10, 0.5, 0.7, 0.8, 0.6),
        UCSC_RefGene_Name = paste0("U", seq_along(ids)),
        GencodeV50_RefGene_Name = paste0("G", seq_along(ids)),
        GencodeV50_NonAnnotated_RefGene_Name = rep("", length(ids)),
        check.names = FALSE
    )
    summaries <- list(
        diagnosticSummaries = list(
            Profession = data.frame(
                Coefficient = c("Profession1", "Profession2")
            ),
            Age = data.frame(Coefficient = "Age")
        ),
        phenotypes = c("Profession", "Age"),
        settings = list(
            interactionTerm = NULL, omnibusTest = TRUE,
            factorLevels = list(Profession = c("0", "1", "2"))
        )
    )

    result <- dnaEPICO:::generateModelVennDDnaEpico(
        annotatedResults = annotated, modelSummaries = summaries,
        analysis = "GLM", vennDPhenotypes = "Profession,Age",
        vennDLabels = "EMD,FIRE,Age",
        vennDOmnibusPhenotypes = "Profession,Age",
        vennDOmnibusLabels = "Profession,Age", outputDir = tmp,
        plotWidth = 500, plotHeight = 450, plotDPI = 72
    )

    expect_s3_class(result, "dnaEPICO_vennD_plots")
    expect_equal(nrow(result$files), 8L)
    expect_true(all(file.exists(result$files$path)))
    expect_setequal(unique(result$files$threshold), c("nominal", "genomeWide"))
    expect_setequal(unique(result$files$annotation), c("UCSC", "GENCODEv50"))
    expect_false(any(grepl("suggestive|combined", result$files$filename)))
    expect_true(any(grepl(
        "vennD_GLM_EMD_FIRE_Age_nominal_UCSC.tiff",
        result$files$filename, fixed = TRUE
    )))
    expect_setequal(names(result$sheets), c(
        "vennDNOM", "vennDSUGG", "vennDGW",
        "vennDOmnibusNOM", "vennDOmnibusSUGG", "vennDOmnibusGW"
    ))
    expect_identical(
        result$mappings$coefficient$pValueColumn,
        c("Profession1P.Value", "Profession2P.Value", "AgeP.Value")
    )
    expect_identical(
        result$mappings$coefficient$displayLabel, c("EMD", "FIRE", "Age")
    )
    expect_true(any(result$metadataRows$Value == "0"))
    expect_false(file.exists(file.path(tmp, "figure_manifest.tsv")))
})

test_that("Venn gene counts keep UCSC and GENCODE sources separate", {
    testthat::skip_if_not_installed("ggVennDiagram")
    annotation <- data.frame(
        UCSC_RefGene_Name = c("", "UCSC2"),
        GencodeV50_RefGene_Name = c("GENCODE1", ""),
        GencodeV50_NonAnnotated_RefGene_Name = c("NEAREST1", "NEAREST2"),
        stringsAsFactors = FALSE, check.names = FALSE
    )

    expect_identical(
        dnaEPICO:::assignVennGenesDnaEpico(
            annotation, "UCSC_RefGene_Name"
        ),
        c(NA_character_, "UCSC2")
    )
    expect_identical(
        dnaEPICO:::assignVennGenesDnaEpico(
            annotation, c(
                "GencodeV50_RefGene_Name",
                "GencodeV50_NonAnnotated_RefGene_Name"
            )
        ),
        c("GENCODE1", "NEAREST2")
    )

    labels <- dnaEPICO:::vennRegionLabelsDnaEpico(
        sets = list(A = c("cg1", "cg2"), B = "cg2"),
        geneMap = c(cg1 = "GENE1", cg2 = "GENE2")
    )
    expect_false(any(grepl("\\(\\s|\\s\\)", labels$displayLabel)))
})

test_that("model-level Venn validates labels, omnibus settings, and annotation", {
    data <- data.frame(
        IlmnID = c("cg1", "cg2"), AP.Value = c(0.01, 0.2),
        BP.Value = c(0.02, 0.3), UCSC_RefGene_Name = c("A", "B"),
        GencodeV50_RefGene_Name = c("A", "B"),
        GencodeV50_NonAnnotated_RefGene_Name = c("", ""),
        check.names = FALSE
    )
    summaries <- list(
        diagnosticSummaries = list(
            A = data.frame(Coefficient = "A"),
            B = data.frame(Coefficient = "B")
        ), phenotypes = c("A", "B"),
        settings = list(
            interactionTerm = NULL, omnibusTest = FALSE,
            factorLevels = list()
        )
    )
    expect_error(
        dnaEPICO:::generateModelVennDDnaEpico(
            data, summaries, analysis = "GLM", vennDPhenotypes = "A,B",
            vennDLabels = "same,SAME"
        ), "unique without regard"
    )
    expect_error(
        dnaEPICO:::generateModelVennDDnaEpico(
            data, summaries, analysis = "GLM",
            vennDOmnibusPhenotypes = "A,B"
        ), "omnibusTest = TRUE"
    )
    expect_error(
        dnaEPICO:::generateModelVennDDnaEpico(
            data[, !grepl("Gencode", names(data)), drop = FALSE], summaries,
            analysis = "GLM", vennDPhenotypes = "A,B"
        ), "release-matched pair of GENCODE"
    )
})

test_that("NULL model-level Venn selection removes only stale analysis figures", {
    tmp <- withr::local_tempdir()
    glm_file <- file.path(tmp, "vennD_GLM_old_nominal_UCSC.tiff")
    lme_file <- file.path(tmp, "vennD_LME_old_nominal_UCSC.tiff")
    writeLines("old", glm_file)
    writeLines("keep", lme_file)

    result <- dnaEPICO:::generateModelVennDDnaEpico(
        annotatedResults = data.frame(), modelSummaries = list(),
        analysis = "GLM", outputDir = tmp
    )
    expect_s3_class(result, "dnaEPICO_vennD_plots")
    expect_false(file.exists(glm_file))
    expect_true(file.exists(lme_file))
})

test_that("ordinary model-level Venn resolves nlme coefficient columns", {
    data <- data.frame(
        IlmnID = c("cg1", "cg2", "cg3"),
        Profession_1_P.Value = c(0.01, 0.2, 0.3),
        Profession_2_P.Value = c(0.02, 0.3, 0.4),
        Age_Age_P.Value = c(0.03, 0.4, 0.5),
        UCSC_RefGene_Name = c("A", "B", "C"),
        GencodeV50_RefGene_Name = c("A", "B", "C"),
        GencodeV50_NonAnnotated_RefGene_Name = c("", "", ""),
        check.names = FALSE
    )
    summaries <- list(
        diagnosticSummaries = list(
            Profession = data.frame(Interaction.Term = c(
                "Profession.1", "Profession.2"
            )), Age = data.frame(Interaction.Term = "Age")
        ), phenotypes = c("Profession", "Age"),
        settings = list(
            interactionTerm = NULL, omnibusTest = FALSE,
            lmeEngine = "nlme", factorLevels = list()
        )
    )
    mapping <- dnaEPICO:::resolveCoefficientVennColumnsDnaEpico(
        data, summaries, c("Profession", "Age"), analysis = "LME"
    )
    expect_identical(mapping$pValueColumn, c(
        "Profession_1_P.Value", "Profession_2_P.Value", "Age_Age_P.Value"
    ))
})

test_that("annotated workbook places Venn sheets before the dictionary", {
    path <- file.path(withr::local_tempdir(), "annotatedGLM.xlsx")
    dnaEPICO:::writeAnnotatedWorkbookMethylationGLM(
        annotated_df = data.frame(IlmnID = "cg1"), file = path,
        resultSheet = "annotatedGLM",
        metadata = data.frame(Key = "analysis", Value = "methylationGLM"),
        extraSheets = list(
            vennDNOM = data.frame(A_ONLY = "cg1"),
            vennDSUGG = data.frame(A_ONLY = character()),
            vennDGW = data.frame(A_ONLY = character())
        ),
        dictionary = data.frame(
            Column = "IlmnID", Description = "Probe", Formula = ""
        )
    )
    expect_identical(openxlsx::getSheetNames(path), c(
        "annotatedGLM", "metadata", "vennDNOM", "vennDSUGG", "vennDGW",
        "dictionary"
    ))
})
