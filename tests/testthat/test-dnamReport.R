create_dnam_report_example <- function(path) {
    dirs <- list(
        qcDir = file.path(path, "figures", "preprocessingMinfiEwasWater", "enmix"),
        preprocessingDir = file.path(path, "figures", "preprocessingMinfiEwasWater", "qc"),
        postprocessingDir = file.path(path, "figures", "preprocessingMinfiEwasWater", "metrics"),
        svaDir = file.path(path, "figures", "svaEnmix"),
        glmDir = file.path(path, "figures", "methylationGLM_T1"),
        glmmDir = file.path(path, "figures", "methylationGLMM_T1T2"),
        glmTableDir = file.path(path, "data", "qpasst1", "methylationGLM_T1"),
        lmerTableDir = file.path(path, "data", "qpasst1", "methylationGLMM_T1T2"),
        phenoTabDir = file.path(path, "data", "qpasst1", "preprocessingMinfiEwasWater"),
        figDir = file.path(path, "reports", "figures"),
        outputDir = file.path(path, "reports"),
        logDir = file.path(path, "logs")
    )

    invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

    glm_subdir <- file.path(dirs$glmDir, "status")
    glmm_subdir <- file.path(dirs$glmmDir, "score")
    dir.create(glm_subdir, recursive = TRUE, showWarnings = FALSE)
    dir.create(glmm_subdir, recursive = TRUE, showWarnings = FALSE)

    file.create(file.path(dirs$qcDir, "qc_1.jpg"))
    file.create(file.path(dirs$preprocessingDir, "pre_1.tiff"))
    file.create(file.path(dirs$postprocessingDir, "post_1.tiff"))
    file.create(file.path(dirs$svaDir, "sva_1.tiff"))
    file.create(file.path(glm_subdir, "qqplot_status.tiff"))
    file.create(file.path(glmm_subdir, "qqplot_score.tiff"))

    dirs$glmTablePath <- file.path(dirs$glmTableDir, "annotatedGLM.xlsx")
    dirs$lmerTablePath <- file.path(dirs$lmerTableDir, "annotatedLME.xlsx")
    dirs$phenoTab <- file.path(dirs$phenoTabDir, "phenoLC.csv")
    utils::write.csv(
        data.frame(
            UID = c("S1", "S2"),
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
    openxlsx::write.xlsx(
        list(
            annotatedGLM = data.frame(
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
            ),
            dictionary = glm_dictionary
        ),
        dirs$glmTablePath,
        overwrite = TRUE
    )
    lmer_dictionary <- data.frame(
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
    openxlsx::write.xlsx(
        list(
            annotatedLME = data.frame(
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
            ),
            dictionary = lmer_dictionary
        ),
        dirs$lmerTablePath,
        overwrite = TRUE
    )

    dirs
}

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
            glmmDir = example_dirs$glmmDir,
            figDir = example_dirs$figDir,
            verbose = FALSE,
            logs = FALSE
        ),
        NA
    )

    expect_s3_class(result, "dnaEPICO_dnamReport_prepared")
    expect_equal(result$figureInventory$qc$count, 1)
    expect_equal(result$figureInventory$glm$count, 1)
    expect_equal(result$figureInventory$glmm$count, 1)
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
            lmerTab = example_dirs$lmerTablePath,
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
    if (identical(result$status, "rendered")) {
        expect_true(file.exists(result$outputFile))
    }
    expect_true(any(grepl("glm\\.qmd$", result$sourceFiles)))
    expect_true(any(grepl("lmer\\.qmd$", result$sourceFiles)))

    glm_qmd <- readLines(file.path(result$projectDir, "glm.qmd"), warn = FALSE)
    lmer_qmd <- readLines(file.path(result$projectDir, "lmer.qmd"), warn = FALSE)
    expect_true(any(grepl(
        "Table 1. Generalised Linear Model Results and Genomic Annotation of CpG Sites by Phenotype(s)",
        glm_qmd,
        fixed = TRUE
    )))
    expect_true(any(grepl("`Table 1` presents", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("`glm2` package", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("annotatedGLM.xlsx", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl("annotatedGLM", glm_qmd, fixed = TRUE)))
    expect_true(any(grepl(
        "Table 1. Linear Mixed-Effects Model Results and Genomic Annotation of CpG Sites by Phenotype(s) and Timepoint",
        lmer_qmd,
        fixed = TRUE
    )))
    expect_true(any(grepl("`Table 1` presents", lmer_qmd, fixed = TRUE)))
    expect_true(any(grepl("`lmer` package", lmer_qmd, fixed = TRUE)))
    expect_true(any(grepl("annotatedLME", lmer_qmd, fixed = TRUE)))
    expect_true(any(grepl("annotatedLME.xlsx", lmer_qmd, fixed = TRUE)))

    enmix_qmd <- readLines(file.path(result$projectDir, "enmix-qc.qmd"), warn = FALSE)
    metrics_qmd <- readLines(file.path(result$projectDir, "metrics.qmd"), warn = FALSE)
    logs_qmd <- readLines(file.path(result$projectDir, "logs.qmd"), warn = FALSE)
    expect_true(any(grepl("`ENmix` produces control plots", enmix_qmd, fixed = TRUE)))
    expect_true(any(grepl("`Figure 1`", metrics_qmd, fixed = TRUE)))
    expect_true(any(grepl("Displays the methylation preprocessing log", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("Displays the phenotype preparation log", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("Displays the hidden-effect and surrogate-variable analysis log", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("### GLM Analysis", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("methylationGLM_T1.txt", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("### LMER Analysis", logs_qmd, fixed = TRUE)))
    expect_true(any(grepl("methylationGLMM_T1T2.txt", logs_qmd, fixed = TRUE)))

    report_qmd <- readLines(file.path(result$projectDir, "report.qmd"), warn = FALSE)
    expect_true(any(grepl("The table has 2 rows and 2 columns.", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("The Data tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("The ENmix QC tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("The GLM Analysis tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("The LMER Analysis tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("The Logs tab summarises", report_qmd, fixed = TRUE)))
    expect_false(any(grepl("Generated Outputs", report_qmd, fixed = TRUE)))

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
