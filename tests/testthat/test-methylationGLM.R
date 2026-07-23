create_methylation_glm_t1_example <- function(path) {
    phenoBT1 <- data.frame(
        Sample_Name = paste0("S", seq_len(8)),
        status = factor(rep(c("Case", "Case", "Control", "Control"), 2)),
        sex = factor(rep(c("F", "M", "F", "M"), 2)),
        cg00000029 = c(0.20, 0.25, 0.22, 0.27, 0.23, 0.29, 0.26, 0.30),
        cg00000108 = c(0.60, 0.55, 0.52, 0.58, 0.62, 0.57, 0.50, 0.56),
        check.names = FALSE
    )

    input_path <- file.path(path, "phenoBT1.RData")
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
    report_assets <- file.path(
        tmp, "reports", "model1", "assets", "results", "glm_results"
    )

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
            reportAssetsDir = report_assets,
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
    expect_equal(result$runSettings$methylationObjectPrefix, "phenoB")
    expect_equal(result$modelFits$settings$internalResponseColumn, "beta")
    expect_length(result$savedFiles$modelFiles, 0L)
    expect_true(file.exists(result$savedFiles$summaryFiles[["status"]]))
    compact_summary <- readRDS(result$savedFiles$summaryFiles[["status"]])
    expect_s3_class(
        compact_summary,
        "dnaEPICO_methylation_phenotype_summary"
    )
    expect_true(compact_summary$complete)
    expect_identical(compact_summary$cpgOrder, result$preparedData$cpgColumns)
    expect_false(any(vapply(
        compact_summary,
        inherits,
        logical(1),
        what = "glm"
    )))
    expect_true(file.exists(result$savedFiles$summaryTxtFiles[["status"]]))
    expect_true(file.exists(result$savedFiles$annotatedGLM))
    expect_true(file.exists(result$savedFiles$annotatedGLMText))
    expect_true(file.exists(result$savedFiles$annotatedGLMReportMetadata))
    expect_true(file.exists(result$savedFiles$annotatedGLMDictionary))
    expect_true(file.exists(result$savedFiles$annotatedGLMMetadata))
    expect_match(result$savedFiles$annotatedGLM, "annotatedGLM\\.xlsx$")
    expect_equal(
        normalizePath(
            dirname(result$savedFiles$annotatedGLMText),
            winslash = "/", mustWork = FALSE
        ),
        normalizePath(report_assets, winslash = "/", mustWork = FALSE)
    )
    expect_false(dir.exists(file.path(
        dirname(result$savedFiles$annotatedGLM), "report-assets"
    )))
    expect_equal(
        sort(basename(unlist(result$savedFiles[c(
            "annotatedGLMText",
            "annotatedGLMReportMetadata",
            "annotatedGLMDictionary",
            "annotatedGLMMetadata"
        )]))),
        sort(c(
            "annotatedGLM.tsv.gz",
            "annotatedGLM.report.tsv",
            "annotatedGLM.dictionary.tsv",
            "annotatedGLM.metadata.tsv"
        ))
    )
    expect_equal(
        openxlsx::getSheetNames(result$savedFiles$annotatedGLM),
        c("annotatedGLM", "metadata", "dictionary")
    )
    dictionary <- openxlsx::read.xlsx(
        result$savedFiles$annotatedGLM,
        sheet = "dictionary",
        check.names = FALSE
    )
    expect_equal(colnames(dictionary), c("Column", "Description", "Formula"))
    expect_true(any(dictionary$Description == "Pvalue from GLM model"))
    expect_true(any(dictionary$Formula == "GLM: Beta values ~ `status` + `sex`"))
    metadata <- openxlsx::read.xlsx(
        result$savedFiles$annotatedGLM,
        sheet = "metadata",
        check.names = FALSE
    )
    expect_equal(metadata$Value[metadata$Key == "backend"], "glm2")
    expect_equal(metadata$Value[metadata$Key == "fitting_function"], "glm2::glm2")
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

test_that("multiple GLM formulas use one consistent response label", {
    dictionary <- dnaEPICO:::buildAnnotatedWorkbookDictionaryMethylationGLM(
        columns = c("CpG", "unmappedP.Value"),
        modelDescription = "Pvalue from GLM model",
        formulaText = c(
            status = "beta ~ `status`",
            age = "beta ~ `age`"
        ),
        modelLabel = "GLM",
        responseLabel = "Beta values"
    )

    expect_equal(
        dictionary$Formula[dictionary$Column == "unmappedP.Value"],
        "GLM: Beta values ~ `status`; Beta values ~ `age`"
    )

    mapped_dictionary <- dnaEPICO:::buildAnnotatedWorkbookDictionaryMethylationGLM(
        columns = c("CpG", "statusCaseP.Value", "ageP.Value"),
        modelDescription = "Pvalue from GLM model",
        formulaText = c(
            status = "beta ~ `status`",
            age = "beta ~ `age`"
        ),
        modelLabel = "GLM",
        responseLabel = "Beta values"
    )

    expect_equal(
        mapped_dictionary$Formula[
            mapped_dictionary$Column == "statusCaseP.Value"
        ],
        "GLM: Beta values ~ `status`"
    )
    expect_equal(
        mapped_dictionary$Formula[
            mapped_dictionary$Column == "ageP.Value"
        ],
        "GLM: Beta values ~ `age`"
    )
})

test_that("GLM native model conditions use one message column", {
    tmp <- withr::local_tempdir()
    example_data <- create_methylation_glm_t1_example(tmp)
    prepared <- prepareMethylationGLMData(
        inputPheno = example_data$inputPheno,
        phenotypes = "status",
        covariates = "sex",
        factorVars = "status,sex",
        cpgLimit = 1,
        logs = FALSE
    )
    fits <- fitMethylationGLMModels(prepared, nCores = 1, logs = FALSE)
    summaries <- summarizeMethylationGLMModels(
        fits, prepared, logs = FALSE
    )
    summary_table <- summaries$summaries$status
    expect_true("Model.Message" %in% names(summary_table))
    expect_false(any(c(
        "Fit.Status", "Converged", "Convergence.Message",
        "Fit.Warning", "Inference.Included", "Exclusion.Reason"
    ) %in% names(summary_table)))

    annotation <- annotateMethylationGLMSummaries(
        modelSummaries = summaries,
        annotationObject = data.frame(
            CpG = "cg00000029",
            chr = "chr1",
            stringsAsFactors = FALSE
        ),
        annotationCols = "chr",
        logs = FALSE
    )
    expect_true("status_Model.Message" %in% names(annotation$data))
    expect_equal(names(annotation$data)[[1L]], "IlmnID")
    expect_lt(
        match("chr", names(annotation$data)),
        match("status_Model.Message", names(annotation$data))
    )
})

test_that("GLM failures are retained internally and omitted without a p-value", {
    failed_fit <- dnaEPICO:::newMethylationFitErrorDnaEpico(
        reason = "The GLM did not produce a result.",
        errorClass = "dnaEPICO_methylationGLM_fit_error"
    )
    fits <- list(status = list(cg_failed = failed_fit))
    model_messages <- dnaEPICO:::collectModelMessagesMethylationGLM(fits)
    failures <- dnaEPICO:::collectFitFailuresMethylationModels(
        fits,
        "dnaEPICO_methylationGLM_fit_error"
    )
    expect_false(model_messages$P.Value.Available)
    expect_match(model_messages$Model.Message, "^ERROR:")

    annotation <- annotateMethylationGLMSummaries(
        modelSummaries = list(
            diagnosticSummaries = list(status = data.frame()),
            summaries = list(status = data.frame()),
            phenotypes = "status",
            fitFailures = failures,
            modelMessages = model_messages
        ),
        annotationObject = data.frame(
            CpG = "cg_failed",
            chr = "chr1",
            stringsAsFactors = FALSE
        ),
        annotationCols = "chr",
        logs = FALSE
    )
    expect_equal(nrow(annotation$data), 0L)
})
test_that("scaleVars standardizes only selected GLM fixed effects", {
    tmp <- withr::local_tempdir()
    example_data <- create_methylation_glm_t1_example(tmp)
    phenoBT1 <- example_data$phenoBT1
    phenoBT1$age <- seq(20, 55, length.out = nrow(phenoBT1))
    save(phenoBT1, file = example_data$inputPheno)

    prepared <- prepareMethylationGLMData(
        inputPheno = example_data$inputPheno,
        phenotypes = "status",
        covariates = "sex,age",
        factorVars = "status,sex",
        scaleVars = "age",
        cpgLimit = 1,
        logs = FALSE
    )

    expect_equal(prepared$data$age, phenoBT1$age)
    expect_equal(mean(prepared$modelData$age), 0, tolerance = 1e-12)
    expect_equal(stats::sd(prepared$modelData$age), 1, tolerance = 1e-12)
    expect_equal(prepared$scaleVars, "age")
    expect_equal(prepared$scalingMetadata$Variable, "age")
    expect_error(
        prepareMethylationGLMData(
            inputPheno = example_data$inputPheno,
            phenotypes = "status",
            covariates = "sex,age",
            factorVars = "status,sex",
            scaleVars = "sex",
            cpgLimit = 1,
            logs = FALSE
        ),
        "both factorVars and scaleVars"
    )
    expect_error(
        prepareMethylationGLMData(
            inputPheno = example_data$inputPheno,
            phenotypes = "status",
            covariates = "sex",
            factorVars = "status,sex",
            scaleVars = "age",
            cpgLimit = 1,
            logs = FALSE
        ),
        "fixed-effect variables"
    )
})

test_that("parallel plans use engine crossover and resource caps", {
    expect_equal(dnaEPICO:::parseFiniteNumericMethylationModels("64"), 64)
    expect_equal(
        dnaEPICO:::parseFiniteNumericMethylationModels(" 1.5e3 "),
        1500
    )
    expect_true(is.na(
        dnaEPICO:::parseFiniteNumericMethylationModels("not-a-number")
    ))

    withr::local_envvar(c(
        DNAEPICO_PARALLEL_BACKEND = "auto",
        DNAEPICO_MAX_WORKERS = "2"
    ))
    expect_equal(dnaEPICO:::parallelCrossoverMethylationModels("glm2"), 25000L)
    expect_equal(dnaEPICO:::parallelCrossoverMethylationModels("lme4"), 1500L)
    expect_equal(dnaEPICO:::parallelCrossoverMethylationModels("nlme"), 5000L)

    small <- dnaEPICO:::resolveParallelPlanMethylationModels(
        engine = "glm2",
        nCores = 8,
        nCpGs = 100
    )
    expect_equal(small$backend, "serial")
    expect_equal(small$workerCount, 1L)
    expect_match(small$reason, "below")

    withr::local_envvar(DNAEPICO_PARALLEL_BACKEND = "psock")
    forced <- dnaEPICO:::resolveParallelPlanMethylationModels(
        engine = "glm2",
        nCores = 8,
        nCpGs = 3
    )
    expect_lte(forced$workerCount, 2L)
    expect_lte(forced$workerCount, 3L)
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
