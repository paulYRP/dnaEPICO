create_methylation_lme_example <- function(path, include_person = TRUE) {
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

create_methylation_lme_ar1_example <- function(path) {
    set.seed(1)
    pheno <- data.frame(
        SID = paste0("S", rep(seq_len(8), each = 3), "_", rep(seq_len(3), times = 8)),
        person = rep(seq_len(8), each = 3),
        Timepoint = factor(rep(c("T1", "T2", "T3"), times = 8), levels = c("T1", "T2", "T3")),
        VisitOrder = rep(seq_len(3), times = 8),
        status = factor(rep(rep(c("Control", "Case"), each = 4), each = 3)),
        sex = factor(rep(c("F", "M"), length.out = 24)),
        cg00000029 =
            0.2 +
                rep(rep(c(0, 0.03), each = 4), each = 3) +
                rep(seq(0, 0.02, length.out = 3), times = 8) +
                stats::rnorm(24, 0, 0.005),
        cg00000108 =
            0.5 +
                rep(rep(c(0, -0.02), each = 4), each = 3) +
                rep(seq(0, 0.01, length.out = 3), times = 8) +
                stats::rnorm(24, 0, 0.005),
        check.names = FALSE
    )

    input_path <- file.path(path, "phenoMT1T2T3.RData")
    save(pheno, file = input_path)

    list(
        inputPheno = input_path,
        annotation = data.frame(
            CpG = c("cg00000029", "cg00000108"),
            chr = c("chr1", "chr1"),
            pos = c(1, 2),
            stringsAsFactors = FALSE
        )
    )
}

test_that("methylationLME returns in-memory results quietly by default", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    testthat::skip_if_not_installed("lmerTest")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_lme_example(tmp, include_person = TRUE)

    expect_message(
        result <- methylationLME(
            inputPheno = example_data$inputPheno,
            outputLogs = file.path(tmp, "logs"),
            outputRData = file.path(tmp, "rData", "methylationLME", "models"),
            outputPlots = file.path(tmp, "figures", "methylationLME"),
            personVar = "person",
            timeVar = "Timepoint",
            phenotypes = "score",
            covariates = "sex",
            factorVars = "sex",
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

    expect_s3_class(result, "dnaEPICO_methylationLME")
    expect_null(result$savedFiles)
    expect_true("score" %in% names(result$modelFits$fits))
    expect_true("score" %in% names(result$modelSummaries$summaries))
    expect_equal(result$modelFits$settings$parallelBackend, "serial")
    expect_equal(nrow(result$modelFits$summaryCache$score), 2)
    expect_s3_class(result$diagnosticPlots$plots$score$qqplot, "ggplot")
    expect_true("IlmnID" %in% colnames(result$annotation$data))
    expect_false(isTRUE(result$preparedData$personCreated))
    expect_equal(result$runSettings$analysisLabel, "methylationLME")
    expect_equal(result$runSettings$internalResponseColumn, "beta")
    expect_false(dir.exists(file.path(tmp, "figures", "methylationLME")))
})

test_that("methylationLME updates the internal response column for copy-number runs", {
    testthat::skip_if_not_installed("lmerTest")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_lme_example(tmp, include_person = TRUE)

    prepared <- prepareMethylationLMEData(
        inputPheno = example_data$inputPheno,
        personVar = "person",
        timeVar = "Timepoint",
        phenotypes = "score",
        covariates = "sex",
        factorVars = "sex",
        cpgLimit = 1,
        methylationScale = "cn",
        verbose = FALSE,
        logs = FALSE
    )
    fits <- fitMethylationLMEModels(
        preparedData = prepared,
        lmeLibs = "lme4,lmerTest",
        nCores = 1,
        verbose = FALSE,
        logs = FALSE
    )

    expect_equal(prepared$internalResponseColumn, "cn")
    expect_equal(prepared$responseLabel, "Copy number values")
    expect_equal(fits$settings$internalResponseColumn, "cn")
    expect_equal(unname(fits$formulas["score"]), "cn ~ `score` + `sex` + (1 | `person` )")
})

test_that("methylationLME can write outputs and derive person IDs on request", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    testthat::skip_if_not_installed("lmerTest")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_lme_example(tmp, include_person = FALSE)

    expect_message(
        result <- methylationLME(
            inputPheno = example_data$inputPheno,
            outputLogs = file.path(tmp, "logs"),
            outputRData = file.path(tmp, "rData", "methylationLME", "models"),
            outputPlots = file.path(tmp, "figures", "methylationLME"),
            personVar = "person",
            timeVar = "Timepoint",
            phenotypes = "score",
            covariates = "sex",
            factorVars = "sex",
            cpgLimit = 2,
            nCores = 1,
            summaryPval = 1,
            saveSignificantInteractions = TRUE,
            significantInteractionDir = file.path(tmp, "results", "cpgs", "methylationLME"),
            significantInteractionPval = 1,
            saveTxtSummaries = TRUE,
            summaryTxtDir = file.path(tmp, "results", "summary", "methylationLME"),
            annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            annotationCols = "Name,chr,pos",
            annotatedLMEOut = file.path(tmp, "data", "methylationLME"),
            display = FALSE,
            verbose = TRUE,
            logs = TRUE,
            saveOutputs = TRUE
        ),
        "Starting DNAm LME Analysis"
    )

    expect_s3_class(result$savedFiles, "dnaEPICO_methylationLME_paths")
    log_file <- file.path(tmp, "logs", "log_methylationLME.txt")
    expect_true(file.exists(log_file))
    log_text <- paste(readLines(log_file, warn = FALSE), collapse = "\n")
    expect_match(log_text, "Fit-time summary rows cached")
    expect_match(log_text, "Summary source:\\s+fit-time cache")
    expect_equal(result$runSettings$methylationObjectPrefix, "phenoBeta")
    expect_equal(result$modelFits$settings$internalResponseColumn, "beta")
    expect_true(file.exists(result$savedFiles$modelFiles[["score"]]))
    expect_true(file.exists(result$savedFiles$summaryFiles[["score"]]))
    expect_true(file.exists(result$savedFiles$summaryTxtFiles[["score"]]))
    expect_true(file.exists(result$savedFiles$annotatedLME))
    expect_match(result$savedFiles$annotatedLME, "annotatedLME\\.xlsx$")
    expect_equal(
        openxlsx::getSheetNames(result$savedFiles$annotatedLME),
        c("annotatedLME", "dictionary", "metadata")
    )
    dictionary <- openxlsx::read.xlsx(
        result$savedFiles$annotatedLME,
        sheet = "dictionary",
        check.names = FALSE
    )
    expect_equal(colnames(dictionary), c("Column", "Description", "Formula"))
    expect_true(any(dictionary$Description == "Pvalue from LME model"))
    expect_true(any(dictionary$Formula == "LME: Beta values ~ `score` + `sex` + (1 | `person` )"))
    metadata <- openxlsx::read.xlsx(
        result$savedFiles$annotatedLME,
        sheet = "metadata",
        check.names = FALSE
    )
    expect_true(all(c("Key", "Value") %in% colnames(metadata)))
    expect_equal(metadata$Value[metadata$Key == "backend"], "lme4")
    expect_equal(metadata$Value[metadata$Key == "fitting_function"], "lmerTest::lmer")
    expect_true(file.exists(file.path(tmp, "figures", "methylationLME", "qqplot_score.tiff")))
    expect_true(isTRUE(result$preparedData$personCreated))
    expect_true(length(result$savedFiles$significantInteractionFiles$score) >= 1)
    expect_true(all(file.exists(result$savedFiles$significantInteractionFiles$score)))
})

test_that("methylationLME supports nlme AR1 models through lmeLibs", {
    testthat::skip_if_not_installed("nlme")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_lme_ar1_example(tmp)

    expect_message(
        result <- methylationLME(
            inputPheno = example_data$inputPheno,
            outputLogs = file.path(tmp, "logs"),
            outputRData = file.path(tmp, "rData", "methylationLME", "models"),
            outputPlots = file.path(tmp, "figures", "methylationLME"),
            personVar = "person",
            timeVar = "Timepoint",
            phenotypes = "status",
            covariates = "sex",
            factorVars = "status,sex",
            lmeLibs = "nlme",
            correlationStructure = "AR1",
            correlationVar = "VisitOrder",
            methylationScale = "m",
            cpgLimit = 2,
            nCores = 1,
            summaryPval = 1,
            saveSignificantInteractions = TRUE,
            significantInteractionPval = 1,
            annotationPackage = example_data$annotation,
            annotationCols = "chr,pos",
            display = FALSE,
            verbose = FALSE,
            logs = TRUE,
            saveOutputs = FALSE
        ),
        NA
    )

    expect_s3_class(result, "dnaEPICO_methylationLME")
    expect_equal(result$modelFits$settings$lmeLibs, "nlme")
    expect_equal(result$modelFits$settings$correlationStructure, "AR1")
    expect_equal(result$modelFits$settings$correlationVar, "VisitOrder")
    expect_equal(result$runSettings$correlationStructure, "AR1")
    expect_equal(result$runSettings$correlationVar, "VisitOrder")
    expect_equal(result$runSettings$methylationScale, "m")
    expect_equal(nrow(result$modelFits$summaryCache$status), 2)
    expect_true(all(c("CpG", "Interaction.Term", "Estimate", "Std.Error", "t.value", "P.value") %in%
        colnames(result$modelSummaries$summaries$status)))

    log_text <- paste(readLines(file.path(tmp, "logs", "log_methylationLME.txt"), warn = FALSE), collapse = "\n")
    expect_match(log_text, "LME libraries:\\s+nlme")
    expect_match(log_text, "Correlation structure:\\s+AR1")
    expect_match(log_text, "Correlation variable:\\s+VisitOrder")
    expect_false(grepl("Resolved LME engine", log_text, fixed = TRUE))
})

test_that("methylationLME requires a correlation variable for AR1 models", {
    testthat::skip_if_not_installed("nlme")

    tmp <- withr::local_tempdir()
    example_data <- create_methylation_lme_ar1_example(tmp)

    expect_error(
        methylationLME(
            inputPheno = example_data$inputPheno,
            personVar = "person",
            timeVar = "Timepoint",
            phenotypes = "status",
            covariates = "sex",
            factorVars = "status,sex",
            lmeLibs = "nlme",
            correlationStructure = "AR1",
            methylationScale = "m",
            annotationPackage = example_data$annotation,
            display = FALSE,
            verbose = FALSE,
            logs = FALSE,
            saveOutputs = FALSE
        ),
        "correlationVar must be supplied"
    )

})

test_that("methylationLME formulas support direct interaction variables", {
    formula_text <- dnaEPICO:::buildFormulaMethylationLME(
        phenotype = "caseStatus",
        personVar = "person",
        covariates = "Sex",
        interactionTerm = "Age",
        includeRandomTerm = TRUE
    )

    expect_equal(
        formula_text,
        "beta ~ `caseStatus` * `Age` + `Sex` + (1 | `person` )"
    )
})

test_that("annotateMethylationLMESummaries accepts annotation package names", {
    testthat::skip_if_not_installed("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    testthat::skip_if_not_installed("lmerTest")

    ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()

    annotation_data <- annotateMethylationLMESummaries(
        modelSummaries = ex$modelSummaries,
        annotationObject = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
        annotationCols = "Name,chr,pos",
        verbose = FALSE,
        logs = FALSE
    )

    expect_s3_class(annotation_data, "dnaEPICO_methylationLME_annotation")
    expect_true("IlmnID" %in% colnames(annotation_data$data))
})
