contains_native_model <- function(x) {
    if (inherits(x, c("glm", "merMod", "lme"))) {
        return(TRUE)
    }
    if (!is.list(x)) {
        return(FALSE)
    }
    any(vapply(x, contains_native_model, logical(1)))
}

test_that("GLM phenotype summaries are complete and resumable", {
    tmp <- withr::local_tempdir()
    phenoBT1 <- data.frame(
        status = factor(rep(c("Control", "Case"), each = 4)),
        sex = factor(rep(c("F", "M"), 4)),
        cg00000029 = c(0.20, 0.21, 0.19, 0.22, 0.40, 0.41, 0.39, 0.42),
        cg00000108 = c(0.60, 0.61, 0.59, 0.62, 0.50, 0.51, 0.49, 0.52),
        check.names = FALSE
    )
    input <- file.path(tmp, "phenoBT1.RData")
    save(phenoBT1, file = input)
    prepared <- prepareMethylationGLMData(
        inputPheno = input, phenotypes = "status", covariates = "sex",
        factorVars = "status,sex", logs = FALSE
    )
    summary_dir <- file.path(tmp, "summaries")

    first <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 1,
        summaryDir = summary_dir, resumeFromSummary = TRUE,
        logs = FALSE
    )
    summary_file <- first$summaryFiles$status
    expect_true(file.exists(summary_file))
    artifact <- readRDS(summary_file)
    expect_s3_class(
        artifact,
        "dnaEPICO_methylation_phenotype_summary"
    )
    expect_identical(artifact$cpgOrder, prepared$cpgColumns)
    expect_identical(as.character(artifact$modelMessages$CpG),
        prepared$cpgColumns)
    expect_false(contains_native_model(artifact))
    expect_false(any(grepl(
        "chunk|checkpoint",
        list.files(summary_dir, all.files = TRUE),
        ignore.case = TRUE
    )))

    poisoned <- prepared
    poisoned$data[prepared$cpgColumns] <- "not numeric"
    resumed <- fitMethylationGLMModels(
        preparedData = poisoned, nCores = 2,
        summaryDir = summary_dir, resumeFromSummary = TRUE,
        logs = FALSE
    )
    expect_identical(
        resumed$summaryCache$status,
        first$summaryCache$status
    )
    expect_identical(
        resumed$coefficientResults$status$pValue,
        first$coefficientResults$status$pValue
    )
})

test_that("GLM signatures ignore one-level factors outside the formula", {
    tmp <- withr::local_tempdir()
    phenoBT1 <- data.frame(
        status = factor(rep(c("Control", "Case"), each = 4)),
        Timepoint = factor(rep("1", 8)),
        cg00000029 = c(0.20, 0.21, 0.19, 0.22, 0.40, 0.41, 0.39, 0.42),
        check.names = FALSE
    )
    input <- file.path(tmp, "phenoBT1.RData")
    save(phenoBT1, file = input)
    prepared <- prepareMethylationGLMData(
        inputPheno = input, phenotypes = "status",
        covariates = character(0), factorVars = "status,Timepoint",
        logs = FALSE
    )
    summary_dir <- file.path(tmp, "summaries")

    expect_error(
        result <- fitMethylationGLMModels(
            preparedData = prepared, nCores = 1,
            summaryDir = summary_dir, logs = FALSE
        ),
        NA
    )
    artifact <- readRDS(result$summaryFiles$status)

    expect_identical(names(artifact$signature$factors), "status")
    expect_false("Timepoint" %in% names(artifact$signature$factors))
    expect_false(inherits(
        result$fits$status[["cg00000029"]],
        "dnaEPICO_methylationGLM_fit_error"
    ))
})

test_that("an unreadable phenotype summary refits the complete phenotype", {
    tmp <- withr::local_tempdir()
    phenoBT1 <- data.frame(
        status = rep(c(0, 1), each = 4),
        cg00000029 = c(0.20, 0.21, 0.19, 0.22, 0.40, 0.41, 0.39, 0.42),
        cg00000108 = c(0.60, 0.61, 0.59, 0.62, 0.50, 0.51, 0.49, 0.52),
        check.names = FALSE
    )
    input <- file.path(tmp, "phenoBT1.RData")
    save(phenoBT1, file = input)
    prepared <- prepareMethylationGLMData(
        inputPheno = input, phenotypes = "status",
        covariates = character(0), factorVars = character(0),
        logs = FALSE
    )
    summary_dir <- file.path(tmp, "summaries")
    first <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 1,
        summaryDir = summary_dir, logs = FALSE
    )
    writeLines("incomplete", first$summaryFiles$status)

    refitted <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 1,
        summaryDir = summary_dir, resumeFromSummary = TRUE,
        logs = FALSE
    )
    artifact <- readRDS(refitted$summaryFiles$status)
    expect_true(artifact$complete)
    expect_identical(artifact$cpgOrder, prepared$cpgColumns)
    expect_equal(nrow(artifact$modelMessages), length(prepared$cpgColumns))
})

test_that("LME compact summaries retain coefficients and omnibus results", {
    testthat::skip_if_not_installed("lmerTest")
    tmp <- withr::local_tempdir()
    set.seed(14)
    phenoBT1T2 <- data.frame(
        person = factor(rep(seq_len(12), each = 2)),
        Timepoint = factor(rep(c("T1", "T2"), 12)),
        group = factor(rep(rep(c("A", "B"), each = 6), each = 2)),
        cg00000029 = rep(seq(0.2, 0.4, length.out = 12), each = 2) +
            rep(c(0, 0.02), 12) + stats::rnorm(24, 0, 0.002),
        cg00000108 = rep(seq(0.5, 0.7, length.out = 12), each = 2) +
            rep(c(0, -0.01), 12) + stats::rnorm(24, 0, 0.002),
        check.names = FALSE
    )
    input <- file.path(tmp, "phenoBT1T2.RData")
    save(phenoBT1T2, file = input)
    prepared <- prepareMethylationLMEData(
        inputPheno = input, personVar = "person", timeVar = "Timepoint",
        phenotypes = "Timepoint", covariates = character(0),
        factorVars = "Timepoint,group", interactionTerm = "group",
        logs = FALSE
    )
    summary_dir <- file.path(tmp, "summaries")
    result <- fitMethylationLMEModels(
        preparedData = prepared, nCores = 1,
        omnibusTest = TRUE, omnibusDdf = "Satterthwaite",
        summaryDir = summary_dir, resumeFromSummary = TRUE,
        logs = FALSE
    )
    artifact <- readRDS(result$summaryFiles$Timepoint)
    expect_false(contains_native_model(artifact))
    expect_identical(artifact$cpgOrder, prepared$cpgColumns)
    expect_equal(nrow(artifact$omnibusTests), length(prepared$cpgColumns))
    expect_true(all(c(
        "estimate", "stdError", "df", "statistic", "pValue"
    ) %in% names(artifact$coefficientResults)))
})

test_that("completed phenotypes are reused independently", {
    tmp <- withr::local_tempdir()
    phenoBT1 <- data.frame(
        status = factor(rep(c("Control", "Case"), each = 5)),
        exposure = seq_len(10),
        cg00000029 = seq(0.20, 0.47, length.out = 10),
        cg00000108 = seq(0.60, 0.42, length.out = 10),
        check.names = FALSE
    )
    input <- file.path(tmp, "phenoBT1.RData")
    save(phenoBT1, file = input)
    prepared <- prepareMethylationGLMData(
        inputPheno = input, phenotypes = c("status", "exposure"),
        covariates = character(0), factorVars = "status",
        logs = FALSE
    )
    summary_dir <- file.path(tmp, "summaries")
    first <- fitMethylationGLMModels(
        preparedData = prepared, nCores = 1,
        summaryDir = summary_dir, logs = FALSE
    )
    status_before <- readRDS(first$summaryFiles$status)
    unlink(first$summaryFiles$exposure)

    poisoned <- prepared
    poisoned$data[prepared$cpgColumns] <- "not numeric"
    resumed <- suppressWarnings(fitMethylationGLMModels(
        preparedData = poisoned, nCores = 1,
        summaryDir = summary_dir, resumeFromSummary = TRUE,
        logs = FALSE
    ))

    expect_identical(
        resumed$summaryCache$status,
        status_before$targetSummary
    )
    expect_equal(resumed$failureCounts[["status"]],
        status_before$failureCount)
    expect_equal(resumed$failureCounts[["exposure"]],
        length(prepared$cpgColumns))
    expect_identical(resumed$resumedPhenotypes, "status")
    expect_identical(resumed$fittedPhenotypes, "exposure")
    expect_true(file.exists(resumed$summaryFiles$exposure))
})


test_that("serial and PSOCK compact results agree for GLM and nlme", {
    testthat::skip_if_not_installed("nlme")
    testthat::skip_if(dnaEPICO:::availableWorkersMethylationModels() < 2L)
    tmp <- withr::local_tempdir()
    set.seed(21)
    person <- factor(rep(seq_len(8), each = 2))
    score <- rep(seq(-1, 1, length.out = 8), each = 2) +
        rep(c(0, 0.2), 8)
    phenoBT1T2 <- data.frame(
        person = person,
        Timepoint = factor(rep(c("T1", "T2"), 8)),
        score = score,
        cg00000029 = 0.3 + 0.04 * score + stats::rnorm(16, 0, 0.01),
        cg00000108 = 0.6 - 0.03 * score + stats::rnorm(16, 0, 0.01),
        check.names = FALSE
    )
    input <- file.path(tmp, "phenoBT1T2.RData")
    save(phenoBT1T2, file = input)

    glm_prepared <- prepareMethylationGLMData(
        inputPheno = input, phenotypes = "score",
        covariates = character(0), factorVars = character(0),
        logs = FALSE
    )
    lme_prepared <- prepareMethylationLMEData(
        inputPheno = input, personVar = "person", timeVar = "Timepoint",
        phenotypes = "score", covariates = character(0),
        factorVars = character(0), logs = FALSE
    )

    withr::local_envvar(DNAEPICO_PARALLEL_BACKEND = "serial")
    glm_serial <- fitMethylationGLMModels(
        glm_prepared, nCores = 2, logs = FALSE
    )
    nlme_serial <- fitMethylationLMEModels(
        lme_prepared, nCores = 2, lmeLibs = "nlme", logs = FALSE
    )

    withr::local_envvar(c(
        DNAEPICO_PARALLEL_BACKEND = "psock",
        DNAEPICO_MAX_WORKERS = "2"
    ))
    glm_psock <- fitMethylationGLMModels(
        glm_prepared, nCores = 2, logs = FALSE
    )
    nlme_psock <- fitMethylationLMEModels(
        lme_prepared, nCores = 2, lmeLibs = "nlme", logs = FALSE
    )

    expect_equal(
        glm_psock$coefficientResults$score,
        glm_serial$coefficientResults$score,
        tolerance = 1e-12
    )
    expect_equal(
        nlme_psock$coefficientResults$score,
        nlme_serial$coefficientResults$score,
        tolerance = 1e-12
    )
})
