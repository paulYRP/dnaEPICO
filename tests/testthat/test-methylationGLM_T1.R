create_methylation_glm_t1_example <- function(path) {
    phenoBT1 <- data.frame(
        Sample_Name = c("S1", "S2", "S3", "S4"),
        status = factor(c("Case", "Case", "Control", "Control")),
        sex = factor(c("F", "M", "F", "M")),
        cg00000029 = c(0.20, 0.25, 0.22, 0.27),
        cg00000108 = c(0.60, 0.55, 0.52, 0.58),
        check.names = FALSE
    )

    input_path <- file.path(path, "phenoBetaT1.RData")
    save(phenoBT1, file = input_path)

    list(
        inputPheno = input_path,
        phenoBT1 = phenoBT1
    )
}

test_that("methylationGLM_T1 returns in-memory results quietly by default", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_glm_t1_example(tmp)

    expect_message(
        result <- methylationGLM_T1(
            inputPheno = example_data$inputPheno,
            outputLogs = file.path(tmp, "logs"),
            outputRData = file.path(tmp, "rData", "methylationGLM_T1", "models"),
            outputPlots = file.path(tmp, "figures", "methylationGLM_T1"),
            phenotypes = "status",
            covariates = "sex",
            factorVars = "status,sex",
            cpgLimit = 2,
            nCores = 1,
            summaryPval = 1,
            saveSignificantCpGs = TRUE,
            significantCpGPval = 1,
            annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            annotationCols = "Name,chr,pos",
            saveTxtSummaries = TRUE,
            display = FALSE,
            verbose = FALSE,
            logs = FALSE,
            saveOutputs = FALSE
        ),
        NA
    )

    expect_s3_class(result, "dnaEPICO_methylationGLM_T1")
    expect_null(result$savedFiles)
    expect_true("status" %in% names(result$modelFits$fits))
    expect_true("status" %in% names(result$modelSummaries$summaries))
    expect_s3_class(result$distributionPlots$phenotypes$status, "ggplot")
    expect_s3_class(result$diagnosticPlots$plots$status$qqplot, "ggplot")
    expect_true("IlmnID" %in% colnames(result$annotation$data))
    expect_false(dir.exists(file.path(tmp, "figures", "methylationGLM_T1")))
})

test_that("methylationGLM_T1 can write outputs and logs on request", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_glm_t1_example(tmp)

    expect_message(
        result <- methylationGLM_T1(
            inputPheno = example_data$inputPheno,
            outputLogs = file.path(tmp, "logs"),
            outputRData = file.path(tmp, "rData", "methylationGLM_T1", "models"),
            outputPlots = file.path(tmp, "figures", "methylationGLM_T1"),
            phenotypes = "status",
            covariates = "sex",
            factorVars = "status,sex",
            cpgLimit = 2,
            nCores = 1,
            summaryPval = 1,
            saveSignificantCpGs = TRUE,
            significantCpGDir = file.path(tmp, "results", "cpgs", "methylationGLM_T1"),
            significantCpGPval = 1,
            saveTxtSummaries = TRUE,
            summaryTxtDir = file.path(tmp, "results", "summary", "methylationGLM_T1"),
            annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            annotationCols = "Name,chr,pos",
            annotatedGLMOut = file.path(tmp, "data", "methylationGLM_T1"),
            display = FALSE,
            verbose = TRUE,
            logs = TRUE,
            saveOutputs = TRUE
        ),
        "Starting DNAm GLM Analysis"
    )

    expect_s3_class(result$savedFiles, "dnaEPICO_methylationGLM_T1_paths")
    expect_true(file.exists(file.path(tmp, "logs", "log_methylationGLM_T1.txt")))
    expect_true(file.exists(result$savedFiles$modelFiles[["status"]]))
    expect_true(file.exists(result$savedFiles$summaryFiles[["status"]]))
    expect_true(file.exists(result$savedFiles$summaryTxtFiles[["status"]]))
    expect_true(file.exists(result$savedFiles$annotatedGLM))
    expect_true(file.exists(file.path(tmp, "figures", "methylationGLM_T1", "bar_status.tiff")))
    expect_true(file.exists(file.path(tmp, "figures", "methylationGLM_T1", "qqplot_status.tiff")))
    expect_true(length(result$savedFiles$significantCpGFiles$status) >= 1)
    expect_true(all(file.exists(result$savedFiles$significantCpGFiles$status)))
})

test_that("annotateMethylationGLM_T1Summaries accepts annotation package names", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()

    annotation_data <- annotateMethylationGLM_T1Summaries(
        modelSummaries = ex$modelSummaries,
        annotationObject = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        annotationCols = "Name,chr,pos",
        verbose = FALSE,
        logs = FALSE
    )

    expect_s3_class(annotation_data, "dnaEPICO_methylationGLM_T1_annotation")
    expect_true("IlmnID" %in% colnames(annotation_data$data))
})
