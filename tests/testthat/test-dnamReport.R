create_dnam_report_example <- function(path) {
    dirs <- list(
        qcDir = file.path(path, "figures", "preprocessingMinfiEwasWater", "enmix"),
        preprocessingDir = file.path(path, "figures", "preprocessingMinfiEwasWater", "qc"),
        postprocessingDir = file.path(path, "figures", "preprocessingMinfiEwasWater", "metrics"),
        svaDir = file.path(path, "figures", "svaEnmix"),
        glmDir = file.path(path, "figures", "methylationGLM"),
        lmeDir = file.path(path, "figures", "methylationLME"),
        glmTableDir = file.path(path, "data", "qpasst1", "methylationGLM"),
        lmeTableDir = file.path(path, "data", "qpasst1", "methylationLME"),
        phenoTabDir = file.path(path, "data", "qpasst1", "preprocessingMinfiEwasWater"),
        figDir = file.path(path, "reports", "figures"),
        outputDir = file.path(path, "reports"),
        glmReportDir = file.path(path, "reports", "assets", "results", "glm_results"),
        lmeReportDir = file.path(path, "reports", "assets", "results", "lme_results"),
        logDir = file.path(path, "logs")
    )

    invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

    glm_subdir <- file.path(dirs$glmDir, "status")
    lme_subdir <- file.path(dirs$lmeDir, "score")
    dir.create(glm_subdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(lme_subdir, recursive = TRUE, showWarnings = FALSE)

    file.create(file.path(dirs$qcDir, "qc_1.jpg"))
    file.create(file.path(dirs$preprocessingDir, "pre_1.tiff"))
    file.create(file.path(dirs$postprocessingDir, "post_1.tiff"))
    file.create(file.path(dirs$svaDir, "sva_1.tiff"))
    file.create(file.path(
        glm_subdir, "qqplot_status_coefficientPvalue.tiff"
    ))
    file.create(file.path(
        lme_subdir, "qqplot_score_coefficientPvalue.tiff"
    ))

    dirs$glmTablePath <- file.path(dirs$glmTableDir, "annotatedGLM.xlsx")
    dirs$lmeTablePath <- file.path(dirs$lmeTableDir, "annotatedLME.xlsx")
    dirs$phenoTab <- file.path(dirs$phenoTabDir, "phenoLC.csv")
    utils::write.csv(
        data.frame(
            person = c("P1", "P2"),
            Sentrix_ID = c("Chip1", "Chip1"),
            Timepoint = c(1, 2),
            check.names = FALSE
        ),
        dirs$phenoTab,
        row.names = FALSE
    )
    glm_dictionary <- data.frame(
        Column = c("P.Value", "IlmnID", "Name", "chr", "pos"),
        Description = c(
            "Pvalue from GLM model",
            "CpG probe identifier",
            "CpG probe identifier",
            "Genomic annotation or supporting result column",
            "Genomic annotation or supporting result column"
        ),
        Formula = c(
            "GLM: Beta values ~ `status` + `sex`",
            "",
            "",
            "",
            ""
        ),
        check.names = FALSE
    )
    glm_data <- data.frame(
            IlmnID = "cg00000029",
            Name = "cg00000029",
            P.Value = 0.001,
            chr = "chr16",
            pos = 53434200,
            UCSC_RefGene_Group = "TSS1500",
            UCSC_RefGene_Name = "RBL2",
            Relation_to_Island = "OpenSea",
            GencodeV41_Group = "TSS",
            check.names = FALSE
    )
    openxlsx::write.xlsx(
        list(annotatedGLM = glm_data, dictionary = glm_dictionary),
        dirs$glmTablePath,
        overwrite = TRUE
    )
    dnaEPICO:::writeReportTableSidecarDnaEpico(
        tableData = glm_data,
        workbookFile = dirs$glmTablePath,
        sidecarDir = dirs$glmReportDir,
        sheet = "annotatedGLM",
        idColumn = "IlmnID",
        dictionary = glm_dictionary
    )
    lme_dictionary <- data.frame(
        Column = c("phenotype:Timepoint P.Value", "IlmnID", "Name", "chr", "pos"),
        Description = c(
            "Pvalue from LME model",
            "CpG probe identifier",
            "CpG probe identifier",
            "Genomic annotation or supporting result column",
            "Genomic annotation or supporting result column"
        ),
        Formula = c(
            "LME: Beta values ~ `score` + `Timepoint` + `sex` + (1 | `person` )",
            "",
            "",
            "",
            ""
        ),
        check.names = FALSE
    )
    lme_data <- data.frame(
            IlmnID = "cg00000108",
            Name = "cg00000108",
            "phenotype:Timepoint P.Value" = 0.002,
            chr = "chr3",
            pos = 37417715,
            UCSC_RefGene_Group = "Body",
            UCSC_RefGene_Name = "C3orf35",
            Relation_to_Island = "Island",
            GencodeV41_Group = "gene body",
            check.names = FALSE
    )
    openxlsx::write.xlsx(
        list(annotatedLME = lme_data, dictionary = lme_dictionary),
        dirs$lmeTablePath,
        overwrite = TRUE
    )
    dnaEPICO:::writeReportTableSidecarDnaEpico(
        tableData = lme_data,
        workbookFile = dirs$lmeTablePath,
        sidecarDir = dirs$lmeReportDir,
        sheet = "annotatedLME",
        idColumn = "IlmnID",
        dictionary = lme_dictionary
    )

    dirs
}

test_that("report stage runner preserves direct early returns", {
    expected <- list(ok = FALSE, error = "unavailable")
    stage <- function() expected

    observed <- dnaEPICO:::.runDrStages(
        arguments = list(),
        stages = list(stage),
        parent = environment()
    )

    expect_identical(observed, expected)
})

write_multiple_formula_workbooks <- function(example_dirs, backend = "lme4") {
    phenotypes <- c("T1D", "HbA1c", "C-peptide")
    glm_formulas <- c(
        "GLM: M-values ~ `T1D` + `Age` + `Sex` + `PC1` + `PC2` + `PC3`",
        "GLM: M-values ~ `HbA1c` + `Age` + `Sex` + `PC1` + `PC2` + `PC3`",
        "GLM: M-values ~ `C-peptide` + `Age` + `Sex` + `PC1` + `PC2` + `PC3`"
    )
    lme_formulas <- c(
        "LME: M-values ~ `T1D` + `Age` + `Sex` + `PC1` + `PC2` + `PC3` + (1 | `SID` )",
        "LME: M-values ~ `HbA1c` + `Age` + `Sex` + `PC1` + `PC2` + `PC3` + (1 | `SID` )",
        "LME: M-values ~ `C-peptide` + `Age` + `Sex` + `PC1` + `PC2` + `PC3` + (1 | `SID` )"
    )
    pvalue_columns <- paste0(phenotypes, "_P.Value")
    annotation_columns <- c("IlmnID", "Name", "chr", "pos")

    glm_data <- data.frame(IlmnID = "cg00000029", check.names = FALSE)
    glm_data[pvalue_columns] <- list(0.001, 0.02, 0.03)
    glm_data$Name <- "cg00000029"
    glm_data$chr <- "chr16"
    glm_data$pos <- 53434200
    glm_dictionary <- data.frame(
        Column = c(pvalue_columns, annotation_columns),
        Description = c(
            rep("Pvalue from GLM model", length(pvalue_columns)),
            "CpG probe identifier",
            "CpG probe identifier",
            "Genomic annotation or supporting result column",
            "Genomic annotation or supporting result column"
        ),
        Formula = c(glm_formulas, rep("", length(annotation_columns))),
        check.names = FALSE
    )
    openxlsx::write.xlsx(
        list(annotatedGLM = glm_data, dictionary = glm_dictionary),
        example_dirs$glmTablePath,
        overwrite = TRUE
    )

    lme_data <- data.frame(IlmnID = "cg00000108", check.names = FALSE)
    lme_data[pvalue_columns] <- list(0.002, 0.01, 0.04)
    lme_data$Name <- "cg00000108"
    lme_data$chr <- "chr3"
    lme_data$pos <- 37417715
    lme_dictionary <- data.frame(
        Column = c(pvalue_columns, annotation_columns),
        Description = c(
            rep("Pvalue from LME model", length(pvalue_columns)),
            "CpG probe identifier",
            "CpG probe identifier",
            "Genomic annotation or supporting result column",
            "Genomic annotation or supporting result column"
        ),
        Formula = c(lme_formulas, rep("", length(annotation_columns))),
        check.names = FALSE
    )
    lme_metadata <- data.frame(
        Key = c(
            "backend",
            "fitting_function",
            "libraries",
            "correlation_structure",
            "correlation_variable",
            "interaction_term",
            "phenotypes"
        ),
        Value = c(
            backend,
            if (identical(backend, "nlme")) "nlme::lme" else "lmerTest::lmer",
            if (identical(backend, "nlme")) "nlme" else "lme4,lmerTest",
            if (identical(backend, "nlme")) "AR1" else "none",
            if (identical(backend, "nlme")) "Timepoint" else "None",
            "None",
            paste(phenotypes, collapse = ",")
        ),
        stringsAsFactors = FALSE
    )
    openxlsx::write.xlsx(
        list(annotatedLME = lme_data, dictionary = lme_dictionary, metadata = lme_metadata),
        example_dirs$lmeTablePath,
        overwrite = TRUE
    )
    dnaEPICO:::writeReportTableSidecarDnaEpico(
        tableData = glm_data,
        workbookFile = example_dirs$glmTablePath,
        sidecarDir = example_dirs$glmReportDir,
        sheet = "annotatedGLM",
        idColumn = "IlmnID",
        dictionary = glm_dictionary
    )
    dnaEPICO:::writeReportTableSidecarDnaEpico(
        tableData = lme_data,
        workbookFile = example_dirs$lmeTablePath,
        sidecarDir = example_dirs$lmeReportDir,
        sheet = "annotatedLME",
        idColumn = "IlmnID",
        dictionary = lme_dictionary,
        workbookMetadata = lme_metadata
    )

    invisible(example_dirs)
}

test_that("report table sidecars are read in bounded chunks without losing rows", {
    tmp <- withr::local_tempdir()
    table_path <- file.path(tmp, "annotatedLME.tsv.gz")
    table_data <- data.frame(
        IlmnID = sprintf("cg%08d", seq_len(10025)),
        P.Value = seq_len(10025) / 10026,
        Profession_Model.Message = NA_character_,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    data.table::fwrite(
        table_data,
        file = table_path,
        sep = "\t",
        quote = TRUE,
        compress = "gzip"
    )
    identifiers <- character(0)
    stream_result <- dnaEPICO:::streamReportTableDnaEpico(
        tableFile = table_path,
        chunkSize = 5000,
        expectedRows = nrow(table_data),
        expectedColumns = ncol(table_data),
        chunkHandler = function(chunk, chunk_number) {
            identifiers <<- c(identifiers, chunk$IlmnID)
        }
    )

    expect_equal(stream_result$nRows, 10025)
    expect_equal(stream_result$nChunks, 3)
    expect_equal(stream_result$maximumChunkRows, 5000)
    expect_identical(identifiers, table_data$IlmnID)
})

test_that("prepareDnamReportInputs returns a structured inventory quietly", {
    tmp <- withr::local_tempdir()
    example_dirs <- create_dnam_report_example(tmp)

    expect_message(
        result <- prepareDnamReportInputs(
            outputDir = example_dirs$outputDir,
            qcDir = example_dirs$qcDir,
            preprocessingDir = example_dirs$preprocessingDir,
            postprocessingDir = example_dirs$postprocessingDir,
            svaDir = example_dirs$svaDir,
            glmDir = example_dirs$glmDir,
            lmeDir = example_dirs$lmeDir,
            figDir = example_dirs$figDir,
            verbose = FALSE,
            logs = FALSE
        ),
        NA
    )

    expect_s3_class(result, "dnaEPICO_dnamReport_prepared")
    expect_equal(result$figureInventory$qc$count, 1)
    expect_equal(result$figureInventory$glm$count, 1)
    expect_equal(result$figureInventory$lme$count, 1)
    expect_length(result$missingFigureDirectories, 0)
    expect_match(result$outputFile, "docs/index\\.html$")
})

test_that("dnamReport renders the dashboard and writes logs on request", {
    tmp <- withr::local_tempdir()
    example_dirs <- create_dnam_report_example(tmp)
    inferred_model <- basename(example_dirs$outputDir)
    report_glm_dir <- file.path(
        tmp, "figures", inferred_model, "methylationGLM"
    )
    report_lme_dir <- file.path(
        tmp, "figures", inferred_model, "methylationLME"
    )
    dir.create(report_glm_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(report_lme_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(
        file.path(
            example_dirs$glmDir, "status",
            "qqplot_status_coefficientPvalue.tiff"
        ),
        report_glm_dir
    )
    file.copy(
        file.path(
            example_dirs$lmeDir, "score",
            "qqplot_score_coefficientPvalue.tiff"
        ),
        report_lme_dir
    )
    file.create(file.path(
        report_lme_dir,
        paste0(
            "vennD_LME_Time_Time_x_Diagnosis_Treatment_x_Time_",
            "genomeWide_GENCODEv50.tiff"
        )
    ))
    cpg_detection_path <- file.path(tmp, "cpg-detection.csv")
    sample_detection_path <- file.path(tmp, "sample-detection.csv")
    utils::write.csv(
        data.frame(
            metric = c(
                "Total CpGs assessed",
                "CpGs detected in at least one sample",
                "CpGs detected in all samples",
                "CpGs never detected",
                "CpGs with no valid detection P values"
            ),
            nCpGs = c(10L, 10L, 9L, 0L, 0L)
        ),
        cpg_detection_path,
        row.names = FALSE
    )
    utils::write.csv(
        data.frame(
            UID = c("Sample_1", "Sample_2"),
            nAssessed = c(10L, 10L),
            nDetected = c(10L, 9L),
            pDetected = c(100, 90)
        ),
        sample_detection_path,
        row.names = FALSE
    )

    expect_message(
        result <- dnamReport(
            outputDir = example_dirs$outputDir,
            phenoTab = example_dirs$phenoTab,
            enmixTab = example_dirs$qcDir,
            qcTab = example_dirs$preprocessingDir,
            svaTab = example_dirs$svaDir,
            metricTab = example_dirs$postprocessingDir,
            glmTab = example_dirs$glmTablePath,
            lmeTab = example_dirs$lmeTablePath,
            cpgDetectionPath = cpg_detection_path,
            sampleDetectionPath = sample_detection_path,
            verbose = TRUE,
            logs = TRUE,
            logTab = example_dirs$logDir
        ),
        "Report path:"
    )

    expect_s3_class(result, "dnaEPICO_dnamReport")
    expect_true(result$status %in% c("rendered", "skipped", "failed"))
    expect_s3_class(result$preparedReport, "dnaEPICO_dnamReport_prepared")
    expect_s3_class(result$renderResult, "dnaEPICO_dnamReport_render")
    expect_true(file.exists(result$logFile))
    report_log <- paste(readLines(result$logFile, warn = FALSE), collapse = "\n")
    expect_match(report_log, "GLM report table source: streamed_sidecar", fixed = TRUE)
    expect_match(report_log, "LME report table source: streamed_sidecar", fixed = TRUE)
    if (identical(result$status, "rendered")) {
        expect_true(file.exists(result$outputFile))
    }
    expect_true(any(grepl("glm\\.qmd$", result$sourceFiles)))
    expect_true(any(grepl("lme\\.qmd$", result$sourceFiles)))
    expect_true(any(grepl("glm-visualisations\\.qmd$", result$sourceFiles)))
    expect_true(any(grepl("lme-visualisations\\.qmd$", result$sourceFiles)))
    expect_equal(
        unname(result$resultTableSources[c("GLM", "LME")]),
        c("streamed_sidecar", "streamed_sidecar")
    )

    glm_qmd <- readLines(file.path(result$projectDir, "glm.qmd"), warn = FALSE)
    lme_qmd <- readLines(file.path(result$projectDir, "lme.qmd"), warn = FALSE)
    data_qmd <- readLines(file.path(result$projectDir, "index.qmd"), warn = FALSE)
    expect_true(any(grepl('data-result-key="phenotype_data"', data_qmd, fixed = TRUE)))
    expect_true(any(grepl('data-role="filter-column"', data_qmd, fixed = TRUE)))
    expect_equal(sum(grepl("dnaepico-viewer-pagination", data_qmd, fixed = TRUE)), 2)
    expect_false(any(grepl("DT::datatable", data_qmd, fixed = TRUE)))
    expect_false(any(grepl("htmlwidget", data_qmd, fixed = TRUE)))
    expect_true(any(grepl(
        "Table 1. Generalised Linear Model Results and CpG Annotation by Phenotype",
        glm_qmd,
        fixed = TRUE
    )))
    expect_false(any(grepl("paged viewer", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("fitted with `glm2`", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("The recorded model formula is", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("GLM: Beta values ~ status + sex", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("data-result-key=\"glm_results\"", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl('data-role="workbook-sheet"', glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("Sheet: dictionary", glm_qmd, fixed = TRUE)))
    expect_false(any(grepl("DT::datatable", glm_qmd, fixed = TRUE)))
    expect_false(any(grepl("openxlsx::read.xlsx", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl(
        "Table 1. Linear Mixed-Effects Model Results and CpG Annotation by Phenotype",
        lme_qmd,
        fixed = TRUE
    )))
    expect_false(any(grepl("paged viewer", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("`lmerTest::lmer()`", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("No interaction was fitted", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("The recorded model formula is", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("data-result-key=\"lme_results\"", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl('data-role="workbook-sheet"', lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("Sheet: dictionary", lme_qmd, fixed = TRUE)))
    expect_false(any(grepl("DT::datatable", lme_qmd, fixed = TRUE)))
    expect_false(any(grepl("openxlsx::read.xlsx", lme_qmd, fixed = TRUE)))

    result_assets <- file.path(result$projectDir, "assets", "results")
    expect_true(file.exists(file.path(result_assets, "phenotype_data", "manifest.js")))
    expect_true(file.exists(file.path(result_assets, "phenotype_data", "chunk-0001.js")))
    expect_true(file.exists(file.path(result_assets, "phenotype_data", basename(example_dirs$phenoTab))))
    expect_true(file.exists(file.path(result_assets, "glm_results", "manifest.js")))
    expect_true(file.exists(file.path(result_assets, "glm_results", "chunk-0001.js")))
    expect_false(file.exists(file.path(result_assets, "glm_results", "annotatedGLM.xlsx")))
    expect_true(file.exists(file.path(result_assets, "glm_results", "annotatedGLM.tsv.gz")))
    expect_true(file.exists(file.path(result_assets, "lme_results", "manifest.js")))
    expect_true(file.exists(file.path(result_assets, "lme_results", "chunk-0001.js")))
    expect_false(file.exists(file.path(result_assets, "lme_results", "annotatedLME.xlsx")))
    expect_true(file.exists(file.path(result_assets, "lme_results", "annotatedLME.tsv.gz")))
    expect_true(file.exists(file.path(
        result_assets, "glm_results_dictionary", "manifest.js"
    )))
    expect_true(file.exists(file.path(
        result_assets, "lme_results_dictionary", "manifest.js"
    )))
    expect_true(file.exists(file.path(
        result$projectDir, "assets", "workbooks", "glm", "annotatedGLM.xlsx"
    )))
    expect_true(file.exists(file.path(
        result$projectDir, "assets", "workbooks", "lme", "annotatedLME.xlsx"
    )))
    expect_true(file.exists(file.path(result$projectDir, "assets", "cpg-viewer.js")))
    expect_true(file.exists(file.path(
        result$projectDir, "assets", "figure-viewer.js"
    )))
    expect_true(file.exists(file.path(
        result$projectDir, "assets", "report-interactions.js"
    )))
    viewer_script <- readLines(
        file.path(result$projectDir, "assets", "cpg-viewer.js"),
        warn = FALSE
    )
    expect_true(any(grepl("maxCachedChunks", viewer_script, fixed = TRUE)))
    expect_true(any(grepl("dnaepico:activate-viewer", viewer_script,
        fixed = TRUE
    )))
    expect_true(any(grepl("matchedRow.scrollIntoView", viewer_script,
        fixed = TRUE
    )))
    expect_true(any(grepl("dnaepico-viewer-match", viewer_script,
        fixed = TRUE
    )))
    expect_true(any(grepl("function commitPageNumber", viewer_script,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        '.dnaepico-result-content[data-result-key]', viewer_script,
        fixed = TRUE
    )))
    figure_viewer_script <- readLines(
        file.path(result$projectDir, "assets", "figure-viewer.js"),
        warn = FALSE
    )
    report_interaction_script <- readLines(
        file.path(result$projectDir, "assets", "report-interactions.js"),
        warn = FALSE
    )
    expect_true(any(grepl(
        'window.location.protocol !== "file:"',
        report_interaction_script,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        'document.documentElement.classList.add("dnaepico-file-protocol")',
        report_interaction_script,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "search.remove()", report_interaction_script, fixed = TRUE
    )))
    expect_true(any(grepl(
        'canvas.addEventListener("wheel"',
        figure_viewer_script, fixed = TRUE
    )))
    expect_true(any(grepl(
        'if (!event.ctrlKey && !event.metaKey)',
        figure_viewer_script, fixed = TRUE
    )))
    expect_true(any(grepl(
        'canvas.addEventListener("keydown"',
        figure_viewer_script, fixed = TRUE
    )))
    expect_true(any(grepl(
        "dnaepico.figureZoom", figure_viewer_script, fixed = TRUE
    )))
    expect_true(any(grepl(
        "dnaepico:figure-zoom-change", figure_viewer_script, fixed = TRUE
    )))
    expect_true(any(grepl(
        "window.sessionStorage.setItem", figure_viewer_script,
        fixed = TRUE
    )))
    expect_false(any(grepl(
        "setNormalCanvasHeight", figure_viewer_script, fixed = TRUE
    )))
    expect_true(any(grepl(
        'var expansionObserver = new MutationObserver(scheduleGeometry)',
        figure_viewer_script, fixed = TRUE
    )))
    expect_true(any(grepl(
        "ResizeObserver", figure_viewer_script,
        fixed = TRUE
    )))
    interaction_script <- readLines(file.path(
        result$projectDir, "assets", "report-interactions.js"
    ), warn = FALSE)
    expect_true(any(grepl(
        'setAttribute("role", "button")',
        interaction_script, fixed = TRUE
    )))
    expect_true(any(grepl(
        'event.key !== "Enter"', interaction_script,
        fixed = TRUE
    )))
    expect_false(any(grepl("Download all results (XLSX)", glm_qmd, fixed = TRUE)))
    expect_false(any(grepl("Download all results (XLSX)", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("Download complete workbook (XLSX)", glm_qmd,
        fixed = TRUE
    )))
    expect_true(any(grepl("Download complete workbook (XLSX)", lme_qmd,
        fixed = TRUE
    )))

    quarto_yml <- readLines(file.path(result$projectDir, "_quarto.yml"),
        warn = FALSE
    )
    expect_true(any(grepl("    right:", quarto_yml, fixed = TRUE)))
    expect_false(any(grepl("    left:", quarto_yml, fixed = TRUE)))
    expect_true(any(quarto_yml == "  search:"))
    expect_true(any(quarto_yml == "    location: navbar"))
    expect_true(any(quarto_yml == "    type: overlay"))
    expect_true(any(quarto_yml == "    expandable: false"))
    expect_false(any(quarto_yml == "    search: true"))
    expected_tabs <- c(
        "Data", "Quality Control", "Batch effect", "Metrics",
        "GLM Analysis", "GLM Visualisations", "LME Analysis",
        "LME Visualisations", "Report", "Logs"
    )
    for (tab in expected_tabs) {
        expect_true(any(grepl(paste0('text: "', tab, '"'), quarto_yml,
            fixed = TRUE
        )))
    }

    quality_control_qmd <- readLines(
        file.path(result$projectDir, "quality-control.qmd"), warn = FALSE
    )
    metrics_qmd <- readLines(file.path(result$projectDir, "metrics.qmd"), warn = FALSE)
    logs_qmd <- readLines(file.path(result$projectDir, "logs.qmd"), warn = FALSE)
    expect_false(file.exists(file.path(result$projectDir, "enmix-qc.qmd")))
    expect_true(any(grepl(
        "ENmix produces control plots", quality_control_qmd, fixed = TRUE
    )))
    expect_true(any(grepl("## ENmix QC", quality_control_qmd, fixed = TRUE)))
    expect_true(any(grepl(
        'data-role="figure-select"', quality_control_qmd, fixed = TRUE
    )))
    expect_false(any(grepl(
        'data-role="figure-zoom"', quality_control_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        'data-role="figure-stage"', quality_control_qmd, fixed = TRUE
    )))
    expect_equal(sum(startsWith(
        trimws(quality_control_qmd),
        "::: {.card .dnaepico-summary-table-card"
    )), 2L)
    expect_equal(sum(grepl(
        "dnaepico-summary-table-wrap",
        quality_control_qmd, fixed = TRUE
    )), 2L)
    expect_true(any(grepl(
        'data-role="figure-description"', metrics_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        "`Methylation Analysis` was not found.",
        logs_qmd,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "`Data Preparation` was not found.",
        logs_qmd,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "`Batch Effect` was not found.",
        logs_qmd,
        fixed = TRUE
    )))
    expect_true(any(grepl('::: {.card title="GLM Analysis"}', logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("methylationGLM.txt", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl('::: {.card title="LME Analysis"}', logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("methylationLME.txt", logs_qmd, fixed = TRUE)))

    report_qmd <- readLines(file.path(result$projectDir, "report.qmd"), warn = FALSE)
    expect_true(any(grepl("format:", report_qmd, fixed = TRUE)))
    expect_true(any(grepl("  dashboard:", report_qmd, fixed = TRUE)))
    expect_true(any(grepl("    scrolling: true", report_qmd, fixed = TRUE)))
    expect_true(any(grepl(
        "    orientation: columns", report_qmd, fixed = TRUE
    )))
    expect_true(any(report_qmd == "## Column"))
    expect_false(any(grepl("## {.sidebar", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("dnaepico-page-heading", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("dnaEPICO Report", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("dnaepico-context-note", report_qmd, fixed = TRUE)))
    generated_pages <- list.files(
        result$projectDir, pattern = "[.]qmd$", full.names = TRUE
    )
    page_content <- unlist(lapply(generated_pages, readLines, warn = FALSE),
        use.names = FALSE
    )
    for (page in generated_pages) {
        page_lines <- readLines(page, warn = FALSE)
        expect_true(any(grepl(
            "    orientation: columns", page_lines, fixed = TRUE
        )), info = basename(page))
        expect_true(any(grepl(
            '<script src="assets/report-interactions.js"></script>',
            page_lines,
            fixed = TRUE
        )), info = basename(page))
    }
    expect_false(any(grepl("## {.sidebar", page_content, fixed = TRUE)))
    expect_false(any(grepl("dnaepico-notes-panel", page_content,
        fixed = TRUE
    )))
    expect_false(any(grepl(
        "dnaEPICO analysis report", page_content, fixed = TRUE
    )))
    glm_visualisations_qmd <- readLines(
        file.path(result$projectDir, "glm-visualisations.qmd"),
        warn = FALSE
    )
    lme_visualisations_qmd <- readLines(
        file.path(result$projectDir, "lme-visualisations.qmd"),
        warn = FALSE
    )
    lme_qmd <- readLines(file.path(result$projectDir, "lme.qmd"),
        warn = FALSE
    )
    expect_true(any(grepl("## EWAS QQ Plots", glm_visualisations_qmd,
        fixed = TRUE
    )))
    expect_true(any(grepl("## EWAS QQ Plots", lme_visualisations_qmd,
        fixed = TRUE
    )))
    expect_false(any(grepl("## Model Variables", glm_visualisations_qmd,
        fixed = TRUE
    )))
    expect_false(any(grepl("## Model Diagnostics", glm_visualisations_qmd,
        fixed = TRUE
    )))
    data_qmd <- readLines(file.path(result$projectDir, "index.qmd"),
        warn = FALSE
    )
    glm_qmd <- readLines(file.path(result$projectDir, "glm.qmd"),
        warn = FALSE
    )
    expect_true(any(grepl(
        'class="dnaepico-viewer-controls"', data_qmd, fixed = TRUE
    )))
    expect_false(any(grepl(
        'class="dnaepico-viewer-toolbar"', page_content, fixed = TRUE
    )))
    expect_false(any(grepl(
        'class="dnaepico-viewer-filter"', page_content, fixed = TRUE
    )))
    expect_lt(
        which(grepl("dnaepico-viewer-controls", data_qmd,
            fixed = TRUE
        ))[[1L]],
        which(grepl("dnaepico-content-description", data_qmd,
            fixed = TRUE
        ))[[1L]]
    )
    expect_lt(
        which(grepl("dnaepico-content-description", data_qmd,
            fixed = TRUE
        ))[[1L]],
        which(grepl("dnaepico-viewer-table-wrap", data_qmd,
            fixed = TRUE
        ))[[1L]]
    )
    expect_lt(
        which(grepl("dnaepico-content-notes", glm_qmd,
            fixed = TRUE
        ))[[1L]],
        which(grepl("dnaepico-workbook-selector", glm_qmd,
            fixed = TRUE
        ))[[1L]]
    )
    expect_lt(
        which(grepl("dnaepico-content-notes", lme_qmd,
            fixed = TRUE
        ))[[1L]],
        which(grepl("dnaepico-workbook-selector", lme_qmd,
            fixed = TRUE
        ))[[1L]]
    )
    expect_true(any(grepl(
        '<script src="assets/figure-viewer.js"></script>',
        glm_visualisations_qmd,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        '<script src="assets/report-interactions.js"></script>',
        glm_visualisations_qmd,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        '"previewPath":', glm_visualisations_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        '"downloadPath":', glm_visualisations_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        '"downloadName":', glm_visualisations_qmd, fixed = TRUE
    )))
    expect_false(any(grepl(
        "figure-fullscreen", page_content, fixed = TRUE
    )))
    expect_false(any(grepl("Full screen", page_content, fixed = TRUE)))
    expect_false(any(grepl(
        "requestFullscreen", page_content, fixed = TRUE
    )))
    expect_false(any(grepl(
        "exitFullscreen", page_content, fixed = TRUE
    )))
    expect_true(any(grepl(
        ".card .dnaepico-selected-figure",
        glm_visualisations_qmd, fixed = TRUE
    )))
    expect_false(any(grepl(
        'class="dnaepico-figure-card-toolbar"',
        glm_visualisations_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        'class="dnaepico-figure-actions"',
        glm_visualisations_qmd, fixed = TRUE
    )))
    figure_action_positions <- which(grepl(
        "dnaepico-figure-actions", glm_visualisations_qmd,
        fixed = TRUE
    ))
    figure_card_positions <- which(grepl(
        ".card .dnaepico-selected-figure", glm_visualisations_qmd,
        fixed = TRUE
    ))
    expect_lt(
        max(figure_action_positions),
        max(figure_card_positions)
    )
    expect_true(any(grepl(
        "data-browser-id=", glm_visualisations_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        "expandable=\"true\" fill=\"false\"",
        glm_visualisations_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        ".card .dnaepico-result-card",
        glm_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        "expandable=\"false\" fill=\"false\"",
        glm_qmd, fixed = TRUE
    )))
    control_position <- which(grepl(
        "dnaepico-viewer-controls", glm_qmd, fixed = TRUE
    ))[[1L]]
    result_card_positions <- which(grepl(
        "dnaepico-result-card", glm_qmd, fixed = TRUE
    ))
    expect_true(any(result_card_positions > control_position))
    batch_qmd <- readLines(
        file.path(result$projectDir, "batch-effect.qmd"), warn = FALSE
    )
    for (figure_page in list(metrics_qmd, batch_qmd)) {
        expect_true(any(grepl(
            ".card .dnaepico-selected-figure",
            figure_page, fixed = TRUE
        )))
        expect_true(any(grepl(
            'class="dnaepico-figure-actions"',
            figure_page, fixed = TRUE
        )))
        expect_false(any(grepl("Full screen", figure_page, fixed = TRUE)))
    }
    expect_false(any(grepl(
        "function updateFigureGeometry()", page_content,
        fixed = TRUE
    )))
    site_css <- readLines(
        file.path(result$projectDir, "assets", "qpasst.css"), warn = FALSE
    )
    expect_true(any(grepl("--dna-nav-width", site_css, fixed = TRUE)))
    expect_true(any(grepl(
        "html.dnaepico-file-protocol #quarto-search",
        site_css,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "body.qpasst-report-page .dnaepico-data-value",
        site_css,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "body.qpasst-report-page .dnaepico-report-section",
        site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        "body.qpasst-report-page .dnaepico-report-section li::marker",
        site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        "background: #dff4fd !important", site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        "color: #10202a !important", site_css, fixed = TRUE
    )))
    expect_true(any(grepl("width: 8.5rem", site_css, fixed = TRUE)))
    expect_true(any(grepl(
        "grid-template-rows: none !important", site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        "body.bslib-has-full-screen #quarto-header", site_css,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "visibility: hidden !important", site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        "[data-full-screen='true']", site_css, fixed = TRUE
    )))
    expect_false(any(grepl("data-zoom-mode", site_css, fixed = TRUE)))
    expect_true(any(grepl(
        ".dnaepico-figure-stage img", site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        ".dnaepico-selected-figure[data-full-screen='false']",
        site_css,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "height: clamp(24rem, 62vh, 52rem) !important",
        site_css,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "body.quarto-dashboard div.dnaepico-selected-figure",
        site_css,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "object-fit: contain !important", site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        ".dnaepico-result-card[data-full-screen='true']",
        site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        ".dnaepico-summary-table-card[data-full-screen='true']",
        site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        "body.bslib-has-full-screen .bslib-card > .card-body",
        site_css, fixed = TRUE
    )))
    expect_false(any(grepl(":fullscreen", site_css, fixed = TRUE)))
    expect_false(any(grepl("100vh", site_css, fixed = TRUE)))
    expect_true(any(grepl(
        ".dnaepico-viewer-match > td", site_css, fixed = TRUE
    )))
    expect_false(any(grepl("scheduleZoom", page_content, fixed = TRUE)))
    expect_false(any(grepl("ResizeObserver", page_content, fixed = TRUE)))
    expect_true(any(grepl(
        "#quarto-header .quarto-navbar-tools",
        site_css, fixed = TRUE
    )))
    expect_true(any(grepl(
        "href: https://github.com/paulYRP/dnaEPICO",
        quarto_yml, fixed = TRUE
    )))
    expect_true(any(grepl(
        "href: https://bioconductor.org/packages/dnaEPICO/",
        quarto_yml, fixed = TRUE
    )))
    expect_equal(sum(grepl("target: _blank", quarto_yml, fixed = TRUE)), 2L)
    expect_true(file.exists(file.path(
        result$projectDir, "assets", "dnaEPICORM.svg"
    )))
    expect_true(any(grepl("The table has 2 rows and 3 columns.", report_qmd, fixed = TRUE)))
    expect_true(any(grepl(
        'class="dnaepico-data-value"', report_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        '>glm2</span>', report_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        "Using person as the participant identifier",
        report_qmd,
        fixed = TRUE
    )))
    expect_false(any(grepl(
        "Using Sentrix_ID as the participant identifier",
        report_qmd,
        fixed = TRUE
    )))
    expect_false(any(grepl("The Data tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("The ENmix QC tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("The GLM Analysis tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("The LME Analysis tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("The Logs tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("Generated Outputs", report_qmd, fixed = TRUE)))
    expect_true(any(grepl("<h3>Logs</h3>", report_qmd, fixed = TRUE)))
    expected_report_sections <- c(
        "Summary", "Data and study design", "Methylation preprocessing",
        "Quality control", "Methylation metrics", "Batch-effect assessment",
        "Generalised linear model analysis",
        "Linear mixed-effects analysis", "CpG and gene-overlap figures",
        "Logs", "Observations"
    )
    expect_true(any(grepl("<h2>Analysis Report</h2>", report_qmd,
        fixed = TRUE
    )))
    for (section in expected_report_sections) {
        expect_true(any(grepl(
            paste0("<h3>", section, "</h3>"), report_qmd,
            fixed = TRUE
        )), info = section)
    }
    expect_true(any(grepl("Study data:", report_qmd, fixed = TRUE)))
    expect_true(any(grepl(
        "Statistical analyses:", report_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        "No log file was available for", report_qmd, fixed = TRUE
    )))
    expect_false(any(grepl("At a glance", report_qmd, fixed = TRUE)))
    expect_false(any(grepl(
        "Workflow and provenance", report_qmd, fixed = TRUE
    )))
    expect_false(any(grepl(
        "Points requiring attention", report_qmd, fixed = TRUE
    )))
    expect_false(any(grepl(
        paste0(
            "Technical details are retained, with short explanations for ",
            "readers without a methylation-analysis background."
        ),
        report_qmd, fixed = TRUE
    )))
    expect_false(any(grepl("Workflow Status", report_qmd, fixed = TRUE)))

    venn_title <- paste0(
        "Genome-wide CpG and GENCODE v50 Gene Overlap Across ",
        "Time-related Model Effects"
    )
    venn_description <- paste0(
        "This Venn diagram compares genome-wide significant CpGs and their ",
        "GENCODE v50 gene annotations for the Time main effect, the ",
        "interaction between Time and Diagnosis, and the interaction between ",
        "Treatment and Time."
    )
    expect_true(any(grepl(
        venn_title, lme_visualisations_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(
        venn_description, lme_visualisations_qmd, fixed = TRUE
    )))
    expect_true(any(grepl(venn_title, report_qmd, fixed = TRUE)))
    expect_false(any(grepl(
        "CpG and Gene Overlap: vennD", lme_visualisations_qmd, fixed = TRUE
    )))

    if (identical(result$status, "rendered")) {
        expect_true(file.exists(file.path(result$projectDir, "docs", "search.json")))
        expect_lt(file.info(file.path(result$projectDir, "docs", "index.html"))$size, 100000)
        expect_lt(file.info(file.path(result$projectDir, "docs", "lme.html"))$size, 500000)
        expect_lt(file.info(file.path(result$projectDir, "docs", "glm.html"))$size, 500000)
        data_html <- readLines(file.path(result$projectDir, "docs", "index.html"), warn = FALSE)
        glm_html <- readLines(file.path(result$projectDir, "docs", "glm.html"), warn = FALSE)
        lme_html <- readLines(file.path(result$projectDir, "docs", "lme.html"), warn = FALSE)
        metrics_html <- readLines(
            file.path(result$projectDir, "docs", "metrics.html"), warn = FALSE
        )
        quality_control_html <- readLines(
            file.path(result$projectDir, "docs", "quality-control.html"),
            warn = FALSE
        )
        batch_html <- readLines(
            file.path(result$projectDir, "docs", "batch-effect.html"),
            warn = FALSE
        )
        expect_true(any(grepl('data-role="body"', data_html, fixed = TRUE)))
        expect_false(any(grepl("html-widget", data_html, fixed = TRUE)))
        expect_false(any(grepl("jquery.dataTables", data_html, fixed = TRUE)))
        expect_true(any(grepl('data-role="body"', glm_html, fixed = TRUE)))
        expect_true(any(grepl('data-role="body"', lme_html, fixed = TRUE)))
        expect_false(any(grepl("&lt;table class=", glm_html, fixed = TRUE)))
        expect_false(any(grepl("&lt;table class=", lme_html, fixed = TRUE)))
        expect_false(any(grepl("&lt;label&gt;Rows per page", glm_html, fixed = TRUE)))
        expect_false(any(grepl("&lt;label&gt;Rows per page", lme_html, fixed = TRUE)))
        expect_true(any(grepl(
            'Table 1. Generalised Linear Model Results',
            glm_html,
            fixed = TRUE
        )))
        expect_true(any(grepl(
            'Table 1. Linear Mixed-Effects Model Results',
            lme_html,
            fixed = TRUE
        )))
        expect_true(any(grepl(
            'Data Preview',
            data_html,
            fixed = TRUE
        )))
        expect_true(any(grepl("card-header", glm_html, fixed = TRUE)))
        expect_true(any(grepl("card-header", lme_html, fixed = TRUE)))
        expect_true(any(grepl("card-header", data_html, fixed = TRUE)))
        for (model_html in list(glm_html, lme_html)) {
            expect_true(any(grepl(
                "dnaepico-result-card", model_html, fixed = TRUE
            )))
            expect_true(any(grepl(
                "bslib-full-screen-enter", model_html, fixed = TRUE
            )))
            expect_false(any(grepl(
                "<p>::: {.card .dnaepico-result-card", model_html,
                fixed = TRUE
            )))
        }
        for (figure_html in list(metrics_html, batch_html)) {
            expect_true(any(grepl("card-header", figure_html, fixed = TRUE)))
            expect_true(any(grepl(
                "bslib-full-screen-enter", figure_html, fixed = TRUE
            )))
        }
        expect_true(any(grepl(
            "dnaepico-summary-table-card",
            quality_control_html, fixed = TRUE
        )))
        expect_true(any(grepl(
            "dnaepico-summary-table-wrap",
            quality_control_html, fixed = TRUE
        )))
        expect_true(any(grepl(
            "bslib-full-screen-enter",
            quality_control_html, fixed = TRUE
        )))
    }

    printed <- utils::capture.output(print(result))
    expect_equal(
        printed,
        c(
            "Class type: dnaEPICO_dnamReport",
            paste("Log output path:", normalizePath(result$logFile, winslash = "/", mustWork = FALSE)),
            paste("Report output path:", normalizePath(result$outputFile, winslash = "/", mustWork = FALSE))
        )
    )
    expect_false(any(grepl("\\$preparedReport|\\$figureInventory|\\$renderResult", printed)))
})

test_that("dnamReport limits model pages to the requested workflow scope", {
    tmp <- withr::local_tempdir()
    example_dirs <- create_dnam_report_example(tmp)
    project_dir <- example_dirs$outputDir
    docs_dir <- file.path(project_dir, "docs")
    dir.create(file.path(docs_dir, "assets", "results", "glm_results"),
        recursive = TRUE
    )
    dir.create(file.path(docs_dir, "assets", "results", "lme_results"),
        recursive = TRUE
    )
    writeLines("stale", file.path(project_dir, "glm.qmd"))
    writeLines("stale", file.path(project_dir, "lme.qmd"))
    writeLines("stale", file.path(docs_dir, "glm.html"))
    writeLines("stale", file.path(docs_dir, "lme.html"))
    writeLines("stale", file.path(
        docs_dir, "assets", "results", "glm_results", "stale.txt"
    ))
    writeLines("stale", file.path(
        docs_dir, "assets", "results", "lme_results", "stale.txt"
    ))
    glm_sidecar <- file.path(
        project_dir, "assets", "results", "glm_results", "annotatedGLM.tsv.gz"
    )
    lme_sidecar <- file.path(
        project_dir, "assets", "results", "lme_results", "annotatedLME.tsv.gz"
    )
    expect_true(file.exists(glm_sidecar))
    expect_true(file.exists(lme_sidecar))

    preprocessing_report <- dnamReport(
        outputDir = project_dir,
        phenoTab = example_dirs$phenoTab,
        enmixTab = example_dirs$qcDir,
        qcTab = example_dirs$preprocessingDir,
        svaTab = example_dirs$svaDir,
        metricTab = example_dirs$postprocessingDir,
        glmTab = example_dirs$glmTablePath,
        lmeTab = example_dirs$lmeTablePath,
        modelSections = character(0),
        verbose = FALSE, logs = FALSE, logTab = example_dirs$logDir
    )

    quarto_yml <- readLines(file.path(project_dir, "_quarto.yml"), warn = FALSE)
    logs_qmd <- readLines(file.path(project_dir, "logs.qmd"), warn = FALSE)
    report_qmd <- readLines(file.path(project_dir, "report.qmd"), warn = FALSE)
    expect_identical(preprocessing_report$modelSections, character(0))
    expect_identical(
        unname(preprocessing_report$resultTableSources),
        c("not_requested", "not_requested")
    )
    expect_false(file.exists(file.path(project_dir, "glm.qmd")))
    expect_false(file.exists(file.path(project_dir, "lme.qmd")))
    expect_false(file.exists(file.path(project_dir, "glm-visualisations.qmd")))
    expect_false(file.exists(file.path(project_dir, "lme-visualisations.qmd")))
    expect_false(file.exists(file.path(docs_dir, "glm.html")))
    expect_false(file.exists(file.path(docs_dir, "lme.html")))
    expect_false(dir.exists(file.path(
        docs_dir, "assets", "results", "glm_results"
    )))
    expect_false(dir.exists(file.path(
        docs_dir, "assets", "results", "lme_results"
    )))
    expect_false(any(grepl("href: glm.qmd", quarto_yml, fixed = TRUE)))
    expect_false(any(grepl("href: lme.qmd", quarto_yml, fixed = TRUE)))
    expect_false(any(grepl("MethylationGLM|GLM Analysis", logs_qmd)))
    expect_false(any(grepl("MethylationLME|LME Analysis", logs_qmd)))
    expect_false(any(grepl("<h3>GLM</h3>", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("<h3>LME</h3>", report_qmd, fixed = TRUE)))
    expect_true(file.exists(glm_sidecar))
    expect_true(file.exists(lme_sidecar))

    glm_report <- dnamReport(
        outputDir = project_dir,
        phenoTab = example_dirs$phenoTab,
        enmixTab = example_dirs$qcDir,
        qcTab = example_dirs$preprocessingDir,
        svaTab = example_dirs$svaDir,
        metricTab = example_dirs$postprocessingDir,
        glmTab = example_dirs$glmTablePath,
        lmeTab = example_dirs$lmeTablePath,
        modelSections = "glm",
        verbose = FALSE, logs = FALSE, logTab = example_dirs$logDir
    )
    expect_identical(glm_report$modelSections, "glm")
    expect_true(file.exists(file.path(project_dir, "glm.qmd")))
    expect_true(file.exists(file.path(project_dir, "glm-visualisations.qmd")))
    expect_false(file.exists(file.path(project_dir, "lme.qmd")))
    expect_false(file.exists(file.path(project_dir, "lme-visualisations.qmd")))

    lme_report <- dnamReport(
        outputDir = project_dir,
        phenoTab = example_dirs$phenoTab,
        enmixTab = example_dirs$qcDir,
        qcTab = example_dirs$preprocessingDir,
        svaTab = example_dirs$svaDir,
        metricTab = example_dirs$postprocessingDir,
        glmTab = example_dirs$glmTablePath,
        lmeTab = example_dirs$lmeTablePath,
        modelSections = "lme",
        verbose = FALSE, logs = FALSE, logTab = example_dirs$logDir
    )
    expect_identical(lme_report$modelSections, "lme")
    expect_false(file.exists(file.path(project_dir, "glm.qmd")))
    expect_false(file.exists(file.path(project_dir, "glm-visualisations.qmd")))
    expect_true(file.exists(file.path(project_dir, "lme.qmd")))
    expect_true(file.exists(file.path(project_dir, "lme-visualisations.qmd")))
})

test_that("dnamReport validates requested model sections", {
    expect_identical(
        dnaEPICO:::normalizeModelSectionsDnamReport(c("LME", "glm", "glm")),
        c("glm", "lme")
    )
    expect_identical(
        dnaEPICO:::normalizeModelSectionsDnamReport(character(0)),
        character(0)
    )
    expect_error(
        dnaEPICO:::normalizeModelSectionsDnamReport("unsupported"),
        "Unsupported modelSections"
    )
    expect_error(
        dnaEPICO:::normalizeModelSectionsDnamReport(NA_character_),
        "cannot contain missing"
    )
})

test_that("dnamReport labels and describes nlme reports from workbook metadata", {
    tmp <- withr::local_tempdir()
    example_dirs <- create_dnam_report_example(tmp)
    workbook <- openxlsx::loadWorkbook(example_dirs$lmeTablePath)
    openxlsx::addWorksheet(workbook, "metadata")
    openxlsx::writeData(
        workbook,
        "metadata",
        data.frame(
            Key = c(
                "backend",
                "fitting_function",
                "libraries",
                "correlation_structure",
                "correlation_variable",
                "interaction_term"
            ),
            Value = c("nlme", "nlme::lme", "nlme", "AR1", "Timepoint", "None"),
            stringsAsFactors = FALSE
        )
    )
    openxlsx::saveWorkbook(workbook, example_dirs$lmeTablePath, overwrite = TRUE)
    dnaEPICO:::writeReportTableSidecarDnaEpico(
        tableData = openxlsx::read.xlsx(
            example_dirs$lmeTablePath,
            sheet = "annotatedLME",
            check.names = FALSE
        ),
        workbookFile = example_dirs$lmeTablePath,
        sidecarDir = file.path(tmp, "quarto", "assets", "results", "lme_results"),
        sheet = "annotatedLME",
        idColumn = "IlmnID",
        dictionary = openxlsx::read.xlsx(
            example_dirs$lmeTablePath,
            sheet = "dictionary",
            check.names = FALSE
        ),
        workbookMetadata = openxlsx::read.xlsx(
            example_dirs$lmeTablePath,
            sheet = "metadata",
            check.names = FALSE
        )
    )

    result <- dnamReport(
        outputDir = file.path(tmp, "quarto"),
        phenoTab = example_dirs$phenoTab,
        enmixTab = example_dirs$qcDir,
        qcTab = example_dirs$preprocessingDir,
        svaTab = example_dirs$svaDir,
        metricTab = example_dirs$postprocessingDir,
        glmTab = example_dirs$glmTablePath,
        lmeTab = example_dirs$lmeTablePath,
        verbose = FALSE,
        logs = FALSE,
        logTab = example_dirs$logDir
    )

    quarto_yml <- readLines(file.path(result$projectDir, "_quarto.yml"), warn = FALSE)
    lme_qmd <- readLines(file.path(result$projectDir, "lme.qmd"), warn = FALSE)
    logs_qmd <- readLines(file.path(result$projectDir, "logs.qmd"), warn = FALSE)
    expect_true(any(grepl('text: "LME Analysis"', quarto_yml, fixed = TRUE)))
    expect_true(any(grepl('title: "nlme Analysis"', lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("analysed using `nlme::lme()`", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("`AR1` residual correlation structure", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl('class="dnaepico-model-formula"', lme_qmd, fixed = TRUE)))
    expect_false(any(grepl("~ `", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl('::: {.card title="nlme Analysis"}', logs_qmd, fixed = TRUE)))
})

test_that("dnamReport presents multiple GLM, LME, and nlme formulas by phenotype", {
    tmp <- withr::local_tempdir()
    example_dirs <- create_dnam_report_example(tmp)
    write_multiple_formula_workbooks(example_dirs, backend = "lme4")

    lme_result <- dnamReport(
        outputDir = file.path(tmp, "quarto-lme"),
        phenoTab = example_dirs$phenoTab,
        enmixTab = example_dirs$qcDir,
        qcTab = example_dirs$preprocessingDir,
        svaTab = example_dirs$svaDir,
        metricTab = example_dirs$postprocessingDir,
        glmTab = example_dirs$glmTablePath,
        lmeTab = example_dirs$lmeTablePath,
        verbose = FALSE,
        logs = FALSE,
        logTab = example_dirs$logDir
    )
    glm_qmd <- paste(readLines(file.path(lme_result$projectDir, "glm.qmd"), warn = FALSE), collapse = "\n")
    lme_qmd <- paste(readLines(file.path(lme_result$projectDir, "lme.qmd"), warn = FALSE), collapse = "\n")
    report_qmd <- paste(readLines(file.path(lme_result$projectDir, "report.qmd"), warn = FALSE), collapse = "\n")

    for (page in list(glm_qmd, lme_qmd)) {
        expect_match(page, "3 phenotype-specific models were fitted\\.")
        expect_match(page, "View recorded model formulas \\(3\\)")
        expect_match(page, "<strong>T1D</strong>", fixed = TRUE)
        expect_match(page, "<strong>HbA1c</strong>", fixed = TRUE)
        expect_match(page, "<strong>C-peptide</strong>", fixed = TRUE)
        expect_false(grepl("The recorded model formula is", page, fixed = TRUE))
    }
    expect_match(glm_qmd, "GLM: M-values ~ T1D + Age + Sex + PC1 + PC2 + PC3", fixed = TRUE)
    expect_match(lme_qmd, "LME: M-values ~ T1D + Age + Sex + PC1 + PC2 + PC3 + (1 | SID)", fixed = TRUE)
    expect_false(grepl("View recorded model formulas", report_qmd, fixed = TRUE))
    expect_false(grepl("dnaepico-model-formula-item", report_qmd, fixed = TRUE))

    write_multiple_formula_workbooks(example_dirs, backend = "nlme")
    nlme_result <- dnamReport(
        outputDir = file.path(tmp, "quarto-nlme"),
        phenoTab = example_dirs$phenoTab,
        enmixTab = example_dirs$qcDir,
        qcTab = example_dirs$preprocessingDir,
        svaTab = example_dirs$svaDir,
        metricTab = example_dirs$postprocessingDir,
        glmTab = example_dirs$glmTablePath,
        lmeTab = example_dirs$lmeTablePath,
        verbose = FALSE,
        logs = FALSE,
        logTab = example_dirs$logDir
    )
    nlme_qmd <- paste(readLines(file.path(nlme_result$projectDir, "lme.qmd"), warn = FALSE), collapse = "\n")
    expect_match(nlme_qmd, 'title: "nlme Analysis"', fixed = TRUE)
    expect_match(nlme_qmd, "3 phenotype-specific models were fitted.", fixed = TRUE)
    expect_match(nlme_qmd, "nlme: M-values ~ T1D + Age + Sex + PC1 + PC2 + PC3 + (1 | SID)", fixed = TRUE)
})

test_that("dnamReport aggregates recognised workflow details into the Report", {
    tmp <- withr::local_tempdir()
    example_dirs <- create_dnam_report_example(tmp)
    log_content <- list(
        log_preprocessingMinfiEwasWater.txt = c(
            "Array type: IlluminaHumanMethylation450k",
            "Annotation version: ilmn12.hg19",
            "Normalization methods: adjustedfunnorm",
            "Detection p-value threshold: 0.05",
            "Reference: saliva",
            "Phenotype file loaded with 6 samples and 13 columns.",
            "Samples after filtering: 6"
        ),
        log_preprocessingPheno.txt = c(
            "Using all 6 samples.",
            "Identifier column: person",
            "Timepoints: 1, 2",
            "Beta dimensions: 451304 x 6"
        ),
        log_svaEnmix.txt = c(
            "Using all 6 samples.",
            "Number of surrogate variables (K): 3",
            "Technical terms modelled: SentrixID, SentrixPosition"
        ),
        log_methylationGLM.txt = c(
            "Phenotypes: status",
            "Covariates: age, sex",
            "Interaction term: None",
            "CpG columns retained: 451304",
            "P-value adjustment method: fdr"
        ),
        log_methylationLME.txt = c(
            "Participant key: person",
            paste0(
                "Outputs: model design, longitudinal distributions, ",
                "Manhattan and Venn figures"
            )
        )
    )
    for (filename in names(log_content)) {
        writeLines(
            log_content[[filename]],
            file.path(example_dirs$logDir, filename),
            useBytes = TRUE
        )
    }

    result <- dnamReport(
        outputDir = file.path(tmp, "report-output"),
        phenoTab = example_dirs$phenoTab,
        enmixTab = example_dirs$qcDir,
        qcTab = example_dirs$preprocessingDir,
        svaTab = example_dirs$svaDir,
        metricTab = example_dirs$postprocessingDir,
        glmTab = example_dirs$glmTablePath,
        lmeTab = example_dirs$lmeTablePath,
        verbose = FALSE,
        logs = FALSE,
        logTab = example_dirs$logDir
    )
    report_qmd <- paste(readLines(
        file.path(result$projectDir, "report.qmd"), warn = FALSE
    ), collapse = "\n")

    expected_details <- c(
        "IlluminaHumanMethylation450k",
        "adjustedfunnorm normalisation",
        "beta-matrix dimensions",
        "surrogate variables",
        "phenotypes status",
        "participant key person",
        "The Data tab contains",
        "workflow logs most commonly record"
    )
    for (detail in expected_details) {
        expect_match(report_qmd, detail, fixed = TRUE)
    }
    for (stage in c(
        "Methylation Analysis", "Data Preparation", "Batch Effect",
        "GLM Analysis", "LME Analysis"
    )) {
        expect_match(report_qmd, stage, fixed = TRUE)
    }
    expect_false(grepl(
        "The Logs tab contains one section per analysis stage",
        report_qmd, fixed = TRUE
    ))
})

test_that("dnamReport validates the detection P-value threshold before writing", {
    output_dir <- file.path(withr::local_tempdir(), "quarto")

    expect_error(
        dnamReport(
            outputDir = output_dir,
            detPThreshold = 2,
            verbose = FALSE,
            logs = FALSE
        ),
        "detPThreshold must be a single number between 0 and 1",
        fixed = TRUE
    )
    expect_false(dir.exists(output_dir))
})
