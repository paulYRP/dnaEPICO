create_methylation_glmm_t1t2_example <- function(path, include_person = TRUE) {
    phenoBT1T2 <- data.frame(
        SID = c("P1A", "P1B", "P2A", "P2B", "P3A", "P3B", "P4A", "P4B"),
        Timepoint = factor(c("1", "2", "1", "2", "1", "2", "1", "2")),
        score = c(10, 12, 9, 11, 13, 14, 8, 9),
        sex = factor(c("F", "F", "M", "M", "F", "F", "M", "M")),
        cg00000029 = c(0.25, 0.27, 0.20, 0.22, 0.30, 0.31, 0.18, 0.20),
        cg00000108 = c(0.50, 0.53, 0.55, 0.57, 0.48, 0.49, 0.60, 0.61),
        check.names = FALSE
    )

    if (isTRUE(include_person)) {
        phenoBT1T2$person <- c(1, 1, 2, 2, 3, 3, 4, 4)
    }

    input_path <- file.path(path, "phenoBetaT1T2.RData")
    save(phenoBT1T2, file = input_path)

    list(
        inputPheno = input_path,
        phenoBT1T2 = phenoBT1T2
    )
}

test_that("methylationGLMM_T1T2 returns in-memory results quietly by default", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    testthat::skip_if_not_installed("lmerTest")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_glmm_t1t2_example(tmp, include_person = TRUE)

    expect_message(
        result <- methylationGLMM_T1T2(
            inputPheno = example_data$inputPheno,
            outputLogs = file.path(tmp, "logs"),
            outputRData = file.path(tmp, "rData", "methylationGLMM_T1T2", "models"),
            outputPlots = file.path(tmp, "figures", "methylationGLMM_T1T2"),
            personVar = "person",
            timeVar = "Timepoint",
            phenotypes = "score",
            covariates = "sex",
            factorVars = "sex,Timepoint",
            cpgLimit = 2,
            nCores = 1,
            summaryPval = 1,
            saveSignificantInteractions = TRUE,
            significantInteractionPval = 1,
            annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            annotationCols = "Name,chr,pos",
            display = FALSE,
            verbose = FALSE,
            logs = FALSE,
            saveOutputs = FALSE
        ),
        NA
    )

    expect_s3_class(result, "dnaEPICO_methylationGLMM_T1T2")
    expect_null(result$savedFiles)
    expect_true("score" %in% names(result$modelFits$fits))
    expect_true("score" %in% names(result$modelSummaries$summaries))
    expect_s3_class(result$diagnosticPlots$plots$score$qqplot, "ggplot")
    expect_true("IlmnID" %in% colnames(result$annotation$data))
    expect_false(isTRUE(result$preparedData$personCreated))
    expect_false(dir.exists(file.path(tmp, "figures", "methylationGLMM_T1T2")))
})

test_that("methylationGLMM_T1T2 can write outputs and derive person IDs on request", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    testthat::skip_if_not_installed("lmerTest")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_glmm_t1t2_example(tmp, include_person = FALSE)

    expect_message(
        result <- methylationGLMM_T1T2(
            inputPheno = example_data$inputPheno,
            outputLogs = file.path(tmp, "logs"),
            outputRData = file.path(tmp, "rData", "methylationGLMM_T1T2", "models"),
            outputPlots = file.path(tmp, "figures", "methylationGLMM_T1T2"),
            personVar = "person",
            timeVar = "Timepoint",
            phenotypes = "score",
            covariates = "sex",
            factorVars = "sex,Timepoint",
            cpgLimit = 2,
            nCores = 1,
            summaryPval = 1,
            saveSignificantInteractions = TRUE,
            significantInteractionDir = file.path(tmp, "results", "cpgs", "methylationGLMM_T1T2"),
            significantInteractionPval = 1,
            saveTxtSummaries = TRUE,
            summaryTxtDir = file.path(tmp, "results", "summary", "methylationGLMM_T1T2"),
            annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            annotationCols = "Name,chr,pos",
            annotatedLMEOut = file.path(tmp, "data", "methylationGLMM_T1T2"),
            display = FALSE,
            verbose = TRUE,
            logs = TRUE,
            saveOutputs = TRUE
        ),
        "Starting DNAm LME Analysis"
    )

    expect_s3_class(result$savedFiles, "dnaEPICO_methylationGLMM_T1T2_paths")
    expect_true(file.exists(file.path(tmp, "logs", "log_methylationGLMM_T1T2.txt")))
    expect_true(file.exists(result$savedFiles$modelFiles[["score"]]))
    expect_true(file.exists(result$savedFiles$summaryFiles[["score"]]))
    expect_true(file.exists(result$savedFiles$summaryTxtFiles[["score"]]))
    expect_true(file.exists(result$savedFiles$annotatedLME))
    expect_true(file.exists(file.path(tmp, "figures", "methylationGLMM_T1T2", "qqplot_score.tiff")))
    expect_true(isTRUE(result$preparedData$personCreated))
    expect_true(length(result$savedFiles$significantInteractionFiles$score) >= 1)
    expect_true(all(file.exists(result$savedFiles$significantInteractionFiles$score)))
})
