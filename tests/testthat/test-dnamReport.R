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
    file.create(file.path(glm_subdir, "qqplot_status.tiff"))
    file.create(file.path(lme_subdir, "qqplot_score.tiff"))

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
        Fit.Status = "fitted",
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
        "Table 1. Generalised Linear Model Results and Genomic Annotation of CpG Sites by Phenotype(s)",
        glm_qmd,
        fixed = TRUE
    )))
    expect_false(any(grepl("paged viewer", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("fitted with `glm2`", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("The recorded model formula is", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("GLM: Beta values ~ status + sex", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("data-result-key=\"glm_results\"", glm_qmd, fixed = TRUE)))
    expect_false(any(grepl("DT::datatable", glm_qmd, fixed = TRUE)))
    expect_false(any(grepl("openxlsx::read.xlsx", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl(
        "Table 1. Linear Mixed-Effects Model Results and Genomic Annotation of CpG Sites by Phenotype(s)",
        lme_qmd,
        fixed = TRUE
    )))
    expect_false(any(grepl("paged viewer", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("`lmerTest::lmer()`", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("No phenotype-by-time interaction was fitted", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("The recorded model formula is", lme_qmd, fixed = TRUE)))
    expect_true(any(grepl("data-result-key=\"lme_results\"", lme_qmd, fixed = TRUE)))
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
    expect_true(file.exists(file.path(result$projectDir, "assets", "cpg-viewer.js")))
    viewer_script <- readLines(
        file.path(result$projectDir, "assets", "cpg-viewer.js"),
        warn = FALSE
    )
    expect_true(any(grepl("maxCachedChunks", viewer_script, fixed = TRUE)))
    expect_false(any(grepl("Download all results (XLSX)", glm_qmd, fixed = TRUE)))
    expect_false(any(grepl("Download all results (XLSX)", lme_qmd, fixed = TRUE)))

    enmix_qmd <- readLines(file.path(result$projectDir, "enmix-qc.qmd"), warn = FALSE)
    metrics_qmd <- readLines(file.path(result$projectDir, "metrics.qmd"), warn = FALSE)
    logs_qmd <- readLines(file.path(result$projectDir, "logs.qmd"), warn = FALSE)
    expect_true(any(grepl("`ENmix` produces control plots", enmix_qmd, fixed = TRUE)))
    expect_true(any(grepl("`Figure 1`", metrics_qmd, fixed = TRUE)))
    expect_true(any(grepl("Displays the methylation preprocessing log", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("Displays the phenotype preparation log", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("Displays the hidden-effect and surrogate-variable analysis log", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl('::: {.card title="GLM Analysis"}', logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("methylationGLM.txt", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl('::: {.card title="LME Analysis"}', logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("methylationLME.txt", logs_qmd, fixed = TRUE)))

    report_qmd <- readLines(file.path(result$projectDir, "report.qmd"), warn = FALSE)
    expect_true(any(grepl("The table has 2 rows and 3 columns.", report_qmd, fixed = TRUE)))
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

    if (identical(result$status, "rendered")) {
        expect_lt(file.info(file.path(result$projectDir, "docs", "index.html"))$size, 100000)
        expect_lt(file.info(file.path(result$projectDir, "docs", "lme.html"))$size, 500000)
        expect_lt(file.info(file.path(result$projectDir, "docs", "glm.html"))$size, 500000)
        data_html <- readLines(file.path(result$projectDir, "docs", "index.html"), warn = FALSE)
        glm_html <- readLines(file.path(result$projectDir, "docs", "glm.html"), warn = FALSE)
        lme_html <- readLines(file.path(result$projectDir, "docs", "lme.html"), warn = FALSE)
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
    expect_true(any(grepl('text: "nlme Analysis"', quarto_yml, fixed = TRUE)))
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
