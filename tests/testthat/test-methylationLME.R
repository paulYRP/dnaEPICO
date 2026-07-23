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

    input_path <- file.path(path, "phenoBT1T2.RData")
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

test_that("lme4 native singular messages are retained without changing p-values", {
    model_data <- data.frame(
        person = factor(rep(seq_len(8), each = 3)),
        visit = rep(c(-1, 0, 1), times = 8)
    )
    cpg_values <-
        0.4 +
        (0.03 * model_data$visit) +
        rep(c(0.001, -0.002, 0.001), times = 8)

    fit <- dnaEPICO:::fitCpGModelMethylationLME(
        cpg = "cg_singular",
        cpgValues = cpg_values,
        modelData = model_data,
        formulaText = "beta ~ visit + (1 | person)",
        personVar = "person",
        lmeEngine = "lme4",
        responseVar = "beta"
    )
    expect_false(inherits(fit, "dnaEPICO_methylationLME_fit_error"))
    expect_match(fit$modelMessage, "singular", ignore.case = TRUE)

    summary_table <- dnaEPICO:::summarizeCpGFitMethylationLME(
        cpg = "cg_singular",
        modelObj = fit,
        phenotype = "visit"
    )
    expect_true(is.finite(summary_table$P.value))
    expect_match(summary_table$Model.Message, "singular", ignore.case = TRUE)
    expect_false(any(c(
        "Fit.Status", "Singular.Fit", "Converged",
        "Convergence.Message", "Fit.Warning",
        "Inference.Included", "Exclusion.Reason"
    ) %in% names(summary_table)))

    fit$pValueAvailable <- TRUE
    fits <- list(visit = list(cg_singular = fit))
    model_messages <- dnaEPICO:::collectModelMessagesMethylationLME(fits)
    annotation <- annotateMethylationLMESummaries(
        modelSummaries = list(
            diagnosticSummaries = list(visit = summary_table),
            summaries = list(visit = summary_table),
            phenotypes = "visit",
            fitFailures = data.frame(),
            modelMessages = model_messages,
            omnibusTests = list(),
            settings = list(interactionTerm = NULL)
        ),
        annotationObject = data.frame(
            CpG = "cg_singular",
            chr = "chr1",
            pos = 1,
            stringsAsFactors = FALSE
        ),
        annotationCols = c("chr", "pos"),
        verbose = FALSE,
        logs = FALSE
    )
    expect_equal(annotation$data$IlmnID, "cg_singular")
    expect_match(
        annotation$data$visit_Model.Message,
        "singular",
        ignore.case = TRUE
    )
    expect_lt(
        match("pos", names(annotation$data)),
        match("visit_Model.Message", names(annotation$data))
    )
})

test_that("mixed-model messages do not override returned coefficient p-values", {
    coefficient_table <- matrix(
        c(0.4, 0.01, 40, 1e-08, 0.03, 0.01, 3, 0.01),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(
            c("(Intercept)", "visit"),
            c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
        )
    )
    fit <- list(
        coef = coefficient_table,
        coefficientTerms = c("(Intercept)" = "(Intercept)", visit = "visit"),
        modelMessage = "WARNING: failed to converge with max|grad| = 0.01"
    )
    summary_table <- dnaEPICO:::summarizeCpGFitMethylationLME(
        cpg = "cg_message",
        modelObj = fit,
        phenotype = "visit"
    )
    expect_equal(summary_table$P.value, 0.01)
    expect_match(summary_table$Model.Message, "failed to converge")
})

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
            summaryTxtDir = file.path(tmp, "results", "summary"),
            significantInteractionDir = file.path(tmp, "results", "significant"),
            annotatedLMEOut = file.path(tmp, "data", "methylationLME"),
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
    report_assets <- file.path(
        tmp, "reports", "model1", "assets", "results", "lme_results"
    )

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
            reportAssetsDir = report_assets,
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
    expect_equal(result$runSettings$methylationObjectPrefix, "phenoB")
    expect_equal(result$modelFits$settings$internalResponseColumn, "beta")
    expect_length(result$savedFiles$modelFiles, 0L)
    expect_true(file.exists(result$savedFiles$summaryFiles[["score"]]))
    compact_summary <- readRDS(result$savedFiles$summaryFiles[["score"]])
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
        what = "merMod"
    )))
    expect_true(file.exists(result$savedFiles$summaryTxtFiles[["score"]]))
    expect_true(file.exists(result$savedFiles$annotatedLME))
    expect_true(file.exists(result$savedFiles$annotatedLMEText))
    expect_true(file.exists(result$savedFiles$annotatedLMEReportMetadata))
    expect_true(file.exists(result$savedFiles$annotatedLMEDictionary))
    expect_true(file.exists(result$savedFiles$annotatedLMEMetadata))
    expect_match(result$savedFiles$annotatedLME, "annotatedLME\\.xlsx$")
    expect_equal(
        normalizePath(
            dirname(result$savedFiles$annotatedLMEText),
            winslash = "/", mustWork = FALSE
        ),
        normalizePath(report_assets, winslash = "/", mustWork = FALSE)
    )
    expect_false(dir.exists(file.path(
        dirname(result$savedFiles$annotatedLME), "report-assets"
    )))
    expect_equal(
        sort(basename(unlist(result$savedFiles[c(
            "annotatedLMEText",
            "annotatedLMEReportMetadata",
            "annotatedLMEDictionary",
            "annotatedLMEMetadata"
        )]))),
        sort(c(
            "annotatedLME.tsv.gz",
            "annotatedLME.report.tsv",
            "annotatedLME.dictionary.tsv",
            "annotatedLME.metadata.tsv"
        ))
    )
    expect_equal(
        openxlsx::getSheetNames(result$savedFiles$annotatedLME),
        c("annotatedLME", "metadata", "dictionary")
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
    report_assets <- file.path(
        tmp, "reports", "nlme", "assets", "results", "lme_results"
    )

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
            summaryTxtDir = file.path(tmp, "results", "summary", "nlme"),
            significantInteractionDir = file.path(tmp, "results", "significant", "nlme"),
            annotatedLMEOut = file.path(tmp, "data", "nlme"),
            reportAssetsDir = report_assets,
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
            saveOutputs = TRUE
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
    expect_true(all(c(
        "CpG", "Interaction.Term", "Estimate", "Std.Error", "t.value", "P.value",
        "Model.Message"
    ) %in%
        colnames(result$modelSummaries$summaries$status)))
    expect_true(all(c(
        "Phenotype", "CpG", "Model.Message", "P.Value.Available"
    ) %in% colnames(result$modelFits$modelMessages)))
    expect_true(all(result$modelFits$modelMessages$P.Value.Available))
    expect_equal(
        openxlsx::getSheetNames(result$savedFiles$annotatedLME),
        c("annotatedLME", "metadata", "dictionary")
    )
    metadata <- openxlsx::read.xlsx(
        result$savedFiles$annotatedLME,
        sheet = "metadata",
        check.names = FALSE
    )
    expect_equal(metadata$Value[metadata$Key == "backend"], "nlme")
    expect_equal(metadata$Value[metadata$Key == "fitting_function"], "nlme::lme")
    annotated <- openxlsx::read.xlsx(
        result$savedFiles$annotatedLME,
        sheet = "annotatedLME",
        check.names = FALSE
    )
    expect_true("status_Model.Message" %in% colnames(annotated))
    expect_true(file.exists(file.path(report_assets, "annotatedLME.tsv.gz")))
    expect_false(dir.exists(file.path(
        dirname(result$savedFiles$annotatedLME), "report-assets"
    )))

    log_text <- paste(readLines(file.path(tmp, "logs", "log_methylationLME.txt"), warn = FALSE), collapse = "\n")
    expect_match(log_text, "LME libraries:\\s+nlme")
    expect_match(log_text, "Correlation structure:\\s+AR1")
    expect_match(log_text, "Correlation variable:\\s+VisitOrder")
    expect_false(grepl("Resolved LME engine", log_text, fixed = TRUE))
})

test_that("scaleVars standardizes selected longitudinal fixed effects", {
    tmp <- withr::local_tempdir()
    example_data <- create_methylation_lme_example(tmp)
    prepared <- prepareMethylationLMEData(
        inputPheno = example_data$inputPheno,
        personVar = "person",
        timeVar = "Timepoint",
        phenotypes = "score",
        covariates = "sex",
        factorVars = "sex",
        scaleVars = "score",
        cpgLimit = 1,
        logs = FALSE
    )

    expect_equal(prepared$data$score, example_data$phenoBT1T2$score)
    expect_equal(mean(prepared$modelData$score), 0, tolerance = 1e-12)
    expect_equal(stats::sd(prepared$modelData$score), 1, tolerance = 1e-12)
    expect_equal(prepared$scalingMetadata$Variable, "score")

    prepared_with_correlation_scaled <- prepared
    prepared_with_correlation_scaled$scaleVars <- c(
        prepared_with_correlation_scaled$scaleVars,
        "Timepoint"
    )
    expect_error(
        fitMethylationLMEModels(
            preparedData = prepared_with_correlation_scaled,
            nCores = 1,
            lmeLibs = "nlme",
            correlationStructure = "AR1",
            correlationVar = "Timepoint",
            logs = FALSE
        ),
        "correlationVar cannot also"
    )
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
