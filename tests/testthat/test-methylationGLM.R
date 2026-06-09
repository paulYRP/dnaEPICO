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

test_that("methylationGLM returns in-memory results quietly by default", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_glm_t1_example(tmp)

    expect_message(
        result <- methylationGLM(
            inputPheno = example_data$inputPheno,
            outputLogs = file.path(tmp, "logs"),
            outputRData = file.path(tmp, "rData", "methylationGLM", "models"),
            outputPlots = file.path(tmp, "figures", "methylationGLM"),
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

    expect_s3_class(result, "dnaEPICO_methylationGLM")
    expect_null(result$savedFiles)
    expect_true("status" %in% names(result$modelFits$fits))
    expect_true("status" %in% names(result$modelSummaries$summaries))
    expect_equal(result$modelFits$settings$parallelBackend, "serial")
    expect_equal(nrow(result$modelFits$summaryCache$status), 2)
    expect_s3_class(result$distributionPlots$phenotypes$status, "ggplot")
    expect_s3_class(result$diagnosticPlots$plots$status$qqplot, "ggplot")
    expect_true("IlmnID" %in% colnames(result$annotation$data))
    expect_equal(result$runSettings$analysisLabel, "methylationGLM")
    expect_equal(result$runSettings$internalResponseColumn, "beta")
    expect_false(dir.exists(file.path(tmp, "figures", "methylationGLM")))
})

test_that("methylationGLM updates the internal response column for M-value runs", {
    tmp <- withr::local_tempdir()
    example_data <- create_methylation_glm_t1_example(tmp)

    prepared <- prepareMethylationGLMData(
        inputPheno = example_data$inputPheno,
        phenotypes = "status",
        covariates = "sex",
        factorVars = "status,sex",
        cpgLimit = 1,
        methylationScale = "m",
        verbose = FALSE,
        logs = FALSE
    )
    fits <- fitMethylationGLMModels(
        preparedData = prepared,
        nCores = 1,
        verbose = FALSE,
        logs = FALSE
    )

    expect_equal(prepared$internalResponseColumn, "m")
    expect_equal(prepared$responseLabel, "M-values")
    expect_equal(fits$settings$internalResponseColumn, "m")
    expect_equal(unname(fits$formulas["status"]), "m ~ `status` + `sex`")
})

test_that("methylationGLM can write outputs and logs on request", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_glm_t1_example(tmp)

    expect_message(
        result <- methylationGLM(
            inputPheno = example_data$inputPheno,
            outputLogs = file.path(tmp, "logs"),
            outputRData = file.path(tmp, "rData", "methylationGLM", "models"),
            outputPlots = file.path(tmp, "figures", "methylationGLM"),
            phenotypes = "status",
            covariates = "sex",
            factorVars = "status,sex",
            cpgLimit = 2,
            nCores = 1,
            summaryPval = 1,
            saveSignificantCpGs = TRUE,
            significantCpGDir = file.path(tmp, "results", "cpgs", "methylationGLM"),
            significantCpGPval = 1,
            saveTxtSummaries = TRUE,
            summaryTxtDir = file.path(tmp, "results", "summary", "methylationGLM"),
            annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            annotationCols = "Name,chr,pos",
            annotatedGLMOut = file.path(tmp, "data", "methylationGLM"),
            display = FALSE,
            verbose = TRUE,
            logs = TRUE,
            saveOutputs = TRUE
        ),
        "Starting DNAm GLM Analysis"
    )

    expect_s3_class(result$savedFiles, "dnaEPICO_methylationGLM_paths")
    log_file <- file.path(tmp, "logs", "log_methylationGLM.txt")
    expect_true(file.exists(log_file))
    log_text <- paste(readLines(log_file, warn = FALSE), collapse = "\n")
    expect_match(log_text, "Fit-time summary rows cached")
    expect_match(log_text, "Summary source:\\s+fit-time cache")
    expect_equal(result$runSettings$methylationObjectPrefix, "phenoBeta")
    expect_equal(result$modelFits$settings$internalResponseColumn, "beta")
    expect_true(file.exists(result$savedFiles$modelFiles[["status"]]))
    expect_true(file.exists(result$savedFiles$summaryFiles[["status"]]))
    expect_true(file.exists(result$savedFiles$summaryTxtFiles[["status"]]))
    expect_true(file.exists(result$savedFiles$annotatedGLM))
    expect_match(result$savedFiles$annotatedGLM, "annotatedGLM\\.xlsx$")
    expect_equal(
        openxlsx::getSheetNames(result$savedFiles$annotatedGLM),
        c("annotatedGLM", "dictionary")
    )
    dictionary <- openxlsx::read.xlsx(
        result$savedFiles$annotatedGLM,
        sheet = "dictionary",
        check.names = FALSE
    )
    expect_equal(colnames(dictionary), c("Column", "Description", "Formula"))
    expect_true(any(dictionary$Description == "Pvalue from GLM model"))
    expect_true(any(dictionary$Formula == "GLM: Beta values ~ `status` + `sex`"))
    expect_true(file.exists(file.path(tmp, "figures", "methylationGLM", "bar_status.tiff")))
    expect_true(file.exists(file.path(tmp, "figures", "methylationGLM", "qqplot_status.tiff")))
    expect_true(length(result$savedFiles$significantCpGFiles$status) >= 1)
    expect_true(all(file.exists(result$savedFiles$significantCpGFiles$status)))
})

test_that("annotated workbook dictionary reports M-values when the fitted response uses M scale", {
    model_results <- list(
        settings = list(
            methylationScale = "m",
            internalResponseColumn = "m"
        ),
        fits = list(
            status = list(
                cg00000029 = list(
                    fitted = c(-0.4, 0.2, 1.3),
                    residuals = c(0, 0, 0)
                )
            )
        ),
        formulas = c(status = "m ~ `status` + `sex`")
    )

    dictionary <- dnaEPICO:::buildAnnotatedWorkbookDictionaryMethylationGLM(
        columns = c("CpG", "statusCaseP.Value"),
        modelDescription = "Pvalue from GLM model",
        formulaText = model_results$formulas,
        modelLabel = "GLM",
        responseLabel = dnaEPICO:::inferMethylationValueLabelMethylationGLM(model_results)
    )

    expect_equal(
        dictionary$Formula[dictionary$Column == "statusCaseP.Value"],
        "GLM: M-values ~ `status` + `sex`"
    )
})

test_that("annotated workbook dictionary reports copy number when requested", {
    model_results <- list(
        settings = list(
            methylationScale = "cn",
            internalResponseColumn = "cn"
        ),
        formulas = c(status = "cn ~ `status` + `sex`")
    )

    dictionary <- dnaEPICO:::buildAnnotatedWorkbookDictionaryMethylationGLM(
        columns = c("CpG", "statusCaseP.Value"),
        modelDescription = "Pvalue from GLM model",
        formulaText = model_results$formulas,
        modelLabel = "GLM",
        responseLabel = dnaEPICO:::inferMethylationValueLabelMethylationGLM(model_results)
    )

    expect_equal(
        dictionary$Formula[dictionary$Column == "statusCaseP.Value"],
        "GLM: Copy number values ~ `status` + `sex`"
    )
})

test_that("annotateMethylationGLMSummaries accepts annotation package names", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")

    ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()

    annotation_data <- annotateMethylationGLMSummaries(
        modelSummaries = ex$modelSummaries,
        annotationObject = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        annotationCols = "Name,chr,pos",
        verbose = FALSE,
        logs = FALSE
    )

    expect_s3_class(annotation_data, "dnaEPICO_methylationGLM_annotation")
    expect_true("IlmnID" %in% colnames(annotation_data$data))
})
