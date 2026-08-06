create_omnibus_glm_example <- function(path) {
    set.seed(20260805)
    n <- 72L
    phenoBT1 <- data.frame(
        SID = paste0("S", seq_len(n)),
        group = factor(rep(c("A", "B", "C"), each = n / 3L)),
        modifier = factor(rep(c("M1", "M2", "M3"), times = n / 3L)),
        trend = seq(-1, 1, length.out = n),
        age = stats::rnorm(n, 40, 8),
        stringsAsFactors = TRUE
    )
    design <- stats::model.matrix(~ group + age, data = phenoBT1)
    coefficients <- c(
        "(Intercept)" = 0.45, groupB = 0.05, groupC = 0.10, age = 0.001
    )
    phenoBT1$cg00000029 <- as.numeric(
        design %*% coefficients + stats::rnorm(n, 0, 0.025)
    )
    phenoBT1$cg00000108 <- as.numeric(
        0.55 + 0.04 * phenoBT1$trend + stats::rnorm(n, 0, 0.025)
    )

    input_path <- file.path(path, "phenoBT1.RData")
    save(phenoBT1, file = input_path)
    list(
        data = phenoBT1, input = input_path,
        annotation = data.frame(
            CpG = c("cg00000029", "cg00000108"),
            Name = c("cg00000029", "cg00000108"),
            chr = c("chr16", "chr16"),
            pos = c(1001L, 1002L),
            stringsAsFactors = FALSE
        )
    )
}

test_that("linearHypothesis handles one-degree-of-freedom numeric terms", {
    example <- create_omnibus_glm_example(withr::local_tempdir())
    model_data <- example$data[, c("trend", "age")]
    formula_text <- "beta ~ trend + age"
    target <- dnaEPICO:::resolveOmnibusTargetTermMethylationGLM(
        formula_text, model_data, "trend"
    )
    fit <- dnaEPICO:::fitCpGModelMethylationGLM(
        "cg00000108", example$data$cg00000108, model_data,
        formula_text, omnibusTest = TRUE, omnibusTerm = target
    )

    expect_false(inherits(fit, "dnaEPICO_methylationGLM_fit_error"))
    expect_identical(fit$omnibus$status, "tested")
    expect_equal(fit$omnibus$numeratorDf, 1)
    expect_equal(
        fit$omnibus$fValue,
        unname(fit$coef["trend", "t value"]^2),
        tolerance = 1e-10
    )
    expect_equal(
        fit$omnibus$pValue,
        unname(fit$coef["trend", "Pr(>|t|)"]),
        tolerance = 1e-10
    )

    keep <- example$data$group != "C"
    factor_data <- data.frame(
        group = droplevels(example$data$group[keep]),
        age = example$data$age[keep]
    )
    factor_fit <- dnaEPICO:::fitCpGModelMethylationGLM(
        "cg00000029", example$data$cg00000029[keep], factor_data,
        "beta ~ group + age", omnibusTest = TRUE, omnibusTerm = "group"
    )
    expect_equal(factor_fit$omnibus$numeratorDf, 1)
    expect_equal(
        factor_fit$omnibus$pValue,
        unname(factor_fit$coef["groupB", "Pr(>|t|)"]),
        tolerance = 1e-10
    )
})

test_that("multi-level GLM omnibus tests agree with a partial F test", {
    example <- create_omnibus_glm_example(withr::local_tempdir())
    model_data <- example$data[, c("group", "age")]
    formula_text <- "beta ~ group + age"
    target <- dnaEPICO:::resolveOmnibusTargetTermMethylationGLM(
        formula_text, model_data, "group"
    )
    fit <- dnaEPICO:::fitCpGModelMethylationGLM(
        "cg00000029", example$data$cg00000029, model_data,
        formula_text, omnibusTest = TRUE, omnibusTerm = target
    )

    direct_data <- model_data
    direct_data$beta <- example$data$cg00000029
    reduced <- glm2::glm2(beta ~ age,
        data = direct_data,
        family = stats::gaussian()
    )
    full <- glm2::glm2(beta ~ group + age,
        data = direct_data,
        family = stats::gaussian()
    )
    partial <- stats::anova(reduced, full, test = "F")
    tested_row <- nrow(partial)

    expect_identical(target, "group")
    expect_identical(fit$omnibus$status, "tested")
    expect_equal(fit$omnibus$numeratorDf, 2)
    expect_equal(
        fit$omnibus$fValue,
        partial[["F"]][[tested_row]],
        tolerance = 1e-10
    )
    expect_equal(
        fit$omnibus$pValue,
        partial[["Pr(>F)"]][[tested_row]],
        tolerance = 1e-10
    )

    releveled_data <- model_data
    releveled_data$group <- stats::relevel(releveled_data$group, ref = "C")
    releveled_fit <- dnaEPICO:::fitCpGModelMethylationGLM(
        "cg00000029", example$data$cg00000029, releveled_data,
        formula_text, omnibusTest = TRUE, omnibusTerm = target
    )
    expect_equal(
        releveled_fit$omnibus$pValue,
        fit$omnibus$pValue,
        tolerance = 1e-10
    )
})

test_that("GLM omnibus interaction tests use the complete interaction term", {
    example <- create_omnibus_glm_example(withr::local_tempdir())
    model_data <- example$data[, c("trend", "modifier", "age")]
    formula_text <- "beta ~ trend * modifier + age"
    target <- dnaEPICO:::resolveOmnibusTargetTermMethylationGLM(
        formula_text, model_data, "trend", "modifier"
    )
    fit <- dnaEPICO:::fitCpGModelMethylationGLM(
        "cg00000108", example$data$cg00000108, model_data,
        formula_text, omnibusTest = TRUE, omnibusTerm = target
    )

    expect_identical(target, "trend:modifier")
    expect_identical(fit$omnibus$status, "tested")
    expect_equal(fit$omnibus$numeratorDf, 2)
    selected_terms <- fit$coefficientTerms[
        fit$coefficientTerms == target
    ]
    expect_length(selected_terms, 2L)
})

test_that("an unavailable GLM omnibus test retains coefficient inference", {
    example <- create_omnibus_glm_example(withr::local_tempdir())
    model_data <- example$data[, c("trend", "age")]
    fit <- dnaEPICO:::fitCpGModelMethylationGLM(
        "cg00000108", example$data$cg00000108, model_data,
        "beta ~ trend + age", omnibusTest = TRUE,
        omnibusTerm = "term_not_in_model"
    )

    expect_false(inherits(fit, "dnaEPICO_methylationGLM_fit_error"))
    expect_identical(fit$omnibus$status, "not_estimable")
    expect_match(fit$omnibus$reason, "no estimable")
    expect_true(is.finite(fit$coef["trend", "Pr(>|t|)"]))
    expect_match(fit$modelMessage, "no estimable")
})

test_that("GLM omnibus configuration controls persistence and adjustment", {
    tmp <- withr::local_tempdir()
    example <- create_omnibus_glm_example(tmp)
    prepared <- prepareMethylationGLMData(
        inputPheno = example$input, phenotypes = "group",
        covariates = "age", factorVars = "group", cpgLimit = 2,
        logs = FALSE
    )
    summary_dir <- file.path(tmp, "summaries")

    without_omnibus <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 1, omnibusTest = FALSE,
        summaryDir = summary_dir, resumeFromSummary = TRUE, logs = FALSE
    )
    expect_equal(nrow(without_omnibus$omnibusTests$group), 0L)

    with_omnibus <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 1, omnibusTest = TRUE,
        summaryDir = summary_dir, resumeFromSummary = TRUE, logs = FALSE
    )
    expect_identical(with_omnibus$fittedPhenotypes, "group")
    expect_length(with_omnibus$resumedPhenotypes, 0L)
    expect_identical(unname(with_omnibus$omnibusTargets), "group")
    expect_equal(nrow(with_omnibus$omnibusTests$group), 2L)

    resumed <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 1, omnibusTest = TRUE,
        summaryDir = summary_dir, resumeFromSummary = TRUE, logs = FALSE
    )
    expect_identical(resumed$resumedPhenotypes, "group")
    expect_length(resumed$fittedPhenotypes, 0L)

    summaries <- summarizeMethylationGLMModels(
        modelResults = resumed, preparedData = prepared,
        padjmethod = "fdr", nCores = 1, logs = FALSE
    )
    omnibus <- summaries$omnibusTests$group
    valid <- is.finite(omnibus$Omnibus.P.Value)
    expect_equal(
        omnibus$Omnibus.Adjusted.P.Value[valid],
        stats::p.adjust(omnibus$Omnibus.P.Value[valid], method = "fdr")
    )
})

test_that("GLM significant CpGs follow the configured test", {
    tmp <- withr::local_tempdir()
    example <- create_omnibus_glm_example(tmp)
    prepared <- prepareMethylationGLMData(
        inputPheno = example$input, phenotypes = "group",
        covariates = "age", factorVars = "group", cpgLimit = 2,
        logs = FALSE
    )
    fits <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 1, omnibusTest = TRUE,
        logs = FALSE
    )
    cpgs <- prepared$cpgColumns

    fits$summaryCache$group[["Pr(>|t|)"]] <- 1
    fits$summaryCache$group[["Pr(>|t|)"]][
        fits$summaryCache$group$CpG == cpgs[[1L]]
    ] <- 0.001
    fits$omnibusTests$group$Omnibus.Status <- "tested"
    fits$omnibusTests$group$Omnibus.P.Value <- c(0.20, 0.01)

    omnibus_hits <- collectSignificantCpGsMethylationGLM(
        modelResults = fits, pvalThreshold = 0.05, logs = FALSE
    )
    expect_identical(names(omnibus_hits$group), cpgs[[2L]])

    fits$omnibusTests$group$Omnibus.Status[[2L]] <- "not_estimable"
    unavailable_hits <- collectSignificantCpGsMethylationGLM(
        modelResults = fits, pvalThreshold = 0.05, logs = FALSE
    )
    expect_length(unavailable_hits$group, 0L)

    fits$settings$omnibusTest <- FALSE
    coefficient_hits <- collectSignificantCpGsMethylationGLM(
        modelResults = fits, pvalThreshold = 0.05, logs = FALSE
    )
    expect_identical(names(coefficient_hits$group), cpgs[[1L]])
})

test_that("GLM omnibus collection checks the fitted interaction term", {
    tmp <- withr::local_tempdir()
    example <- create_omnibus_glm_example(tmp)
    prepared <- prepareMethylationGLMData(
        inputPheno = example$input, phenotypes = "group",
        covariates = "age", factorVars = "group", cpgLimit = 1,
        logs = FALSE
    )
    fits <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 1, omnibusTest = TRUE,
        logs = FALSE
    )

    expect_error(
        collectSignificantCpGsMethylationGLM(
            modelResults = fits, pvalThreshold = 0.05,
            interactionTerm = "modifier", logs = FALSE
        ),
        "interactionTerm does not match"
    )
})

test_that("GLM omnibus results propagate to annotated outputs", {
    tmp <- withr::local_tempdir()
    example <- create_omnibus_glm_example(tmp)
    report_assets <- file.path(
        tmp, "reports", "model1", "assets", "results", "glm_results"
    )
    result <- methylationGLM(
        inputPheno = example$input,
        outputLogs = file.path(tmp, "logs"),
        outputRData = file.path(tmp, "rData", "models"),
        outputPlots = file.path(tmp, "figures"),
        phenotypes = "group", covariates = "age", factorVars = "group",
        omnibusTest = TRUE, cpgLimit = 2, nCores = 1, summaryPval = NA,
        saveSignificantCpGs = FALSE,
        summaryTxtDir = file.path(tmp, "results", "summary", "methylationGLM"),
        annotationPackage = example$annotation,
        annotationCols = "Name,chr,pos",
        annotatedGLMOut = file.path(tmp, "data", "methylationGLM"),
        reportAssetsDir = report_assets,
        display = FALSE, verbose = FALSE, logs = TRUE,
        saveOutputs = TRUE
    )

    expected <- paste0(
        "group_Omnibus_",
        c(
            "F.Value", "Num.DF", "Den.DF", "P.Value",
            "Adjusted.P.Value", "Method"
        )
    )
    expect_true(all(expected %in% colnames(result$annotation$data)))
    workbook <- result$savedFiles$annotatedGLM
    workbook_data <- openxlsx::read.xlsx(workbook, sheet = "annotatedGLM")
    expect_true(all(expected %in% colnames(workbook_data)))
    metadata <- openxlsx::read.xlsx(workbook, sheet = "metadata")
    metadata_values <- stats::setNames(metadata$Value, metadata$Key)
    expect_identical(metadata_values[["omnibus_test"]], "TRUE")
    expect_identical(
        metadata_values[["omnibus_method"]],
        "car::linearHypothesis Wald F"
    )
    expect_identical(metadata_values[["p_adjust_method"]], "fdr")
    expect_true(nzchar(metadata_values[["car_version"]]))
    expect_true(file.exists(file.path(report_assets, "annotatedGLM.tsv.gz")))
})

test_that("GLM omnibus testing runs on PSOCK workers", {
    tmp <- withr::local_tempdir()
    example <- create_omnibus_glm_example(tmp)
    prepared <- prepareMethylationGLMData(
        inputPheno = example$input,
        phenotypes = c("group", "trend"), covariates = "age",
        factorVars = "group", cpgLimit = 2, logs = FALSE
    )
    withr::local_envvar(c(
        DNAEPICO_PARALLEL_BACKEND = "psock",
        DNAEPICO_MAX_WORKERS = "2"
    ))

    fits <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 2, omnibusTest = TRUE,
        logs = FALSE
    )
    expect_identical(fits$settings$parallelBackend, "psock")
    expect_true(fits$settings$clusterReusedAcrossPhenotypes)
    expect_setequal(names(fits$omnibusTests), c("group", "trend"))
    expect_true(all(vapply(fits$omnibusTests, nrow, integer(1)) == 2L))
    expect_true(all(vapply(fits$omnibusTests, function(table) {
        all(table$Omnibus.Status == "tested")
    }, logical(1))))
})
