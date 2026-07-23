create_omnibus_lme_example <- function(path, groups = 2L) {
    set.seed(20260722)
    n_people <- 24L
    visits <- c("T1", "T2", "T3")
    profession_levels <- LETTERS[seq_len(groups)]
    profession_by_person <- rep(
        profession_levels,
        length.out = n_people
    )
    phenoBT1T2T3 <- data.frame(
        SID = paste0(
            "P", rep(seq_len(n_people), each = length(visits)),
            "_", rep(visits, times = n_people)
        ),
        person = factor(rep(seq_len(n_people), each = length(visits))),
        Timepoint = factor(
            rep(visits, times = n_people),
            levels = visits
        ),
        Profession = factor(
            rep(profession_by_person, each = length(visits)),
            levels = profession_levels
        ),
        stringsAsFactors = TRUE
    )
    design <- stats::model.matrix(
        ~ Timepoint * Profession,
        data = phenoBT1T2T3
    )
    coefficients <- rep(0, ncol(design))
    names(coefficients) <- colnames(design)
    coefficients[["(Intercept)"]] <- 0.45
    interaction_columns <- grepl(":", names(coefficients), fixed = TRUE)
    coefficients[interaction_columns] <- seq(
        0.04, 0.08,
        length.out = sum(interaction_columns)
    )
    participant_effect <- rep(
        stats::rnorm(n_people, 0, 0.02),
        each = length(visits)
    )
    phenoBT1T2T3$cg00000029 <- as.numeric(
        design %*% coefficients + participant_effect +
            stats::rnorm(nrow(phenoBT1T2T3), 0, 0.01)
    )
    phenoBT1T2T3$cg00000108 <- as.numeric(
        0.55 + participant_effect +
            stats::rnorm(nrow(phenoBT1T2T3), 0, 0.01)
    )

    input_path <- file.path(path, "phenoBT1T2T3.RData")
    save(phenoBT1T2T3, file = input_path)
    list(
        data = phenoBT1T2T3,
        input = input_path,
        annotation = data.frame(
            CpG = c("cg00000029", "cg00000108"),
            Name = c("cg00000029", "cg00000108"),
            chr = c("chr16", "chr16"),
            pos = c(1001L, 1002L),
            stringsAsFactors = FALSE
        )
    )
}

test_that("omnibus interaction tests agree with contestMD", {
    example <- create_omnibus_lme_example(withr::local_tempdir())
    model_data <- example$data[, c(
        "person", "Timepoint", "Profession"
    )]
    formula_text <-
        "beta ~ Timepoint * Profession + (1 | person)"
    target <- dnaEPICO:::resolveOmnibusTargetTermMethylationLME(
        formulaText = formula_text, data = model_data,
        phenotype = "Timepoint", interactionTerm = "Profession"
    )
    fit <- dnaEPICO:::fitCpGModelMethylationLME(
        cpg = "cg00000029", cpgValues = example$data$cg00000029,
        modelData = model_data, formulaText = formula_text,
        personVar = "person", omnibusTest = TRUE,
        omnibusDdf = "Satterthwaite", omnibusTerm = target
    )
    expect_false(inherits(fit, "dnaEPICO_methylationLME_fit_error"))
    expect_identical(target, "Timepoint:Profession")
    expect_identical(fit$omnibus$status, "tested")
    expect_equal(fit$omnibus$numeratorDf, 2)

    direct_data <- model_data
    direct_data$beta <- example$data$cg00000029
    direct_fit <- lmerTest::lmer(
        beta ~ Timepoint * Profession + (1 | person),
        data = direct_data, REML = TRUE
    )
    interaction_names <- names(fit$coefficientTerms)[
        fit$coefficientTerms == target
    ]
    fixed_names <- names(lme4::fixef(direct_fit))
    contrast <- matrix(
        0, nrow = length(interaction_names), ncol = length(fixed_names)
    )
    selected <- match(interaction_names, fixed_names)
    contrast[cbind(seq_along(selected), selected)] <- 1
    direct <- lmerTest::contestMD(
        direct_fit, contrast,
        ddf = "Satterthwaite"
    )
    expect_equal(fit$omnibus$fValue, direct[["F value"]][[1L]])
    expect_equal(fit$omnibus$denominatorDf, direct$DenDF[[1L]])
    expect_equal(fit$omnibus$pValue, direct[["Pr(>F)"]][[1L]])
})

test_that("omnibus tests are invariant to factor reference levels", {
    example <- create_omnibus_lme_example(withr::local_tempdir())
    fit_test <- function(data) {
        formula_text <-
            "beta ~ Timepoint * Profession + (1 | person)"
        target <- dnaEPICO:::resolveOmnibusTargetTermMethylationLME(
            formula_text, data, "Timepoint", "Profession"
        )
        dnaEPICO:::fitCpGModelMethylationLME(
            "cg00000029", example$data$cg00000029, data,
            formula_text, "person", omnibusTest = TRUE,
            omnibusTerm = target
        )$omnibus
    }
    first <- example$data[, c("person", "Timepoint", "Profession")]
    second <- first
    second$Timepoint <- stats::relevel(second$Timepoint, ref = "T3")
    second$Profession <- stats::relevel(second$Profession, ref = "B")

    expect_equal(fit_test(first)$pValue, fit_test(second)$pValue,
        tolerance = 1e-10
    )
})

test_that("contestMD handles one-degree-of-freedom omnibus tests", {
    example <- create_omnibus_lme_example(withr::local_tempdir())
    model_data <- example$data[, c("person", "Timepoint")]
    model_data$trend <- as.numeric(model_data$Timepoint) - 1
    formula_text <- "beta ~ trend + (1 | person)"
    target <- dnaEPICO:::resolveOmnibusTargetTermMethylationLME(
        formula_text, model_data, "trend"
    )
    fit <- dnaEPICO:::fitCpGModelMethylationLME(
        "cg00000029", example$data$cg00000029, model_data,
        formula_text, "person", omnibusTest = TRUE,
        omnibusTerm = target
    )

    expect_identical(fit$omnibus$status, "tested")
    expect_equal(fit$omnibus$numeratorDf, 1)
    coefficient_p <- fit$coef["trend", "Pr(>|t|)"]
    expect_equal(fit$omnibus$pValue, coefficient_p, tolerance = 1e-10)
})

test_that("omnibus configuration is validated before fitting", {
    expect_error(
        dnaEPICO:::normalizeOmnibusDdfMethylationLME("invalid"),
        "Satterthwaite.*Kenward-Roger"
    )
    expect_error(
        fitMethylationLMEModels(
            preparedData = list(), lmeLibs = "nlme",
            omnibusTest = TRUE
        ),
        "only.*lmerTest/lme4.*not nlme"
    )
})

test_that("non-estimable omnibus tests retain the CpG and reason", {
    example <- create_omnibus_lme_example(withr::local_tempdir())
    model_data <- example$data[, c(
        "person", "Timepoint", "Profession"
    )]
    fit <- dnaEPICO:::fitCpGModelMethylationLME(
        "cg00000029", example$data$cg00000029, model_data,
        "beta ~ Timepoint * Profession + (1 | person)",
        "person", omnibusTest = TRUE,
        omnibusTerm = "term_not_in_model"
    )
    expect_false(inherits(fit, "dnaEPICO_methylationLME_fit_error"))
    expect_identical(fit$omnibus$status, "not_estimable")
    expect_match(fit$omnibus$reason, "no estimable")

    summary <- dnaEPICO:::summarizeCpGFitMethylationLME(
        "cg00000029", fit, "Timepoint", "Profession"
    )
    expect_true(all(summary$CpG == "cg00000029"))
    expect_true(any(is.finite(summary$P.value)))
    expect_match(summary$Model.Message[[1L]], "no estimable")
})

test_that("Kenward-Roger omnibus testing is available through pbkrtest", {
    skip_if_not_installed("pbkrtest")
    example <- create_omnibus_lme_example(withr::local_tempdir())
    model_data <- example$data[, c(
        "person", "Timepoint", "Profession"
    )]
    formula_text <-
        "beta ~ Timepoint * Profession + (1 | person)"
    target <- dnaEPICO:::resolveOmnibusTargetTermMethylationLME(
        formula_text, model_data, "Timepoint", "Profession"
    )
    fit <- dnaEPICO:::fitCpGModelMethylationLME(
        "cg00000029", example$data$cg00000029, model_data,
        formula_text, "person", omnibusTest = TRUE,
        omnibusDdf = "Kenward-Roger", omnibusTerm = target
    )
    expect_identical(fit$omnibus$status, "tested")
    expect_identical(fit$omnibus$method, "Kenward-Roger")
    expect_equal(fit$omnibus$numeratorDf, 2)
})

test_that("significant interaction collection uses omnibus p-values", {
    make_fit <- function(coefficient_p, omnibus_p) {
        coefficient_table <- matrix(
            c(
                0.4, 0.01, 40, 1e-12,
                0.02, 0.01, 2, coefficient_p
            ),
            nrow = 2, byrow = TRUE,
            dimnames = list(
                c("(Intercept)", "TimepointT2:ProfessionB"),
                c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
            )
        )
        list(
            coef = coefficient_table,
            coefficientTerms = c(
                "(Intercept)" = "(Intercept)",
                "TimepointT2:ProfessionB" = "Timepoint:Profession"
            ),
            omnibus = list(
                status = "tested", pValue = omnibus_p,
                term = "Timepoint:Profession"
            )
        )
    }
    model_results <- list(
        fits = list(Timepoint = list(
            cg_component_only = make_fit(0.001, 0.5),
            cg_omnibus = make_fit(0.5, 0.001)
        )),
        settings = list(
            omnibusTest = TRUE,
            interactionTerm = "Profession"
        )
    )
    significant <- collectSignificantInteractionsMethylationLME(
        modelResults = model_results, pvalThreshold = 0.05,
        interactionTerm = "Profession", logs = FALSE
    )
    expect_setequal(names(significant$Timepoint), "cg_omnibus")
})

test_that("multiple phenotypes receive separate omnibus tests", {
    tmp <- withr::local_tempdir()
    example <- create_omnibus_lme_example(tmp)
    phenoBT1T2T3 <- example$data
    set.seed(17)
    person_exposure <- stats::rnorm(length(levels(phenoBT1T2T3$person)))
    phenoBT1T2T3$Exposure <- rep(person_exposure, each = 3)
    input <- file.path(tmp, "phenoBT1T2T3_multi.RData")
    save(phenoBT1T2T3, file = input)
    prepared <- prepareMethylationLMEData(
        inputPheno = input, personVar = "person", timeVar = "Timepoint",
        phenotypes = c("Timepoint", "Exposure"), covariates = NULL,
        factorVars = "Timepoint,Profession",
        interactionTerm = "Profession", cpgLimit = 2,
        logs = FALSE
    )
    withr::local_envvar(c(
        DNAEPICO_PARALLEL_BACKEND = "psock",
        DNAEPICO_MAX_WORKERS = "2"
    ))
    fits <- suppressWarnings(fitMethylationLMEModels(
        preparedData = prepared, nCores = 2,
        omnibusTest = TRUE, omnibusDdf = "Satterthwaite",
        logs = FALSE
    ))
    expect_setequal(names(fits$omnibusTests), c("Timepoint", "Exposure"))
    expect_identical(
        unname(fits$omnibusTargets),
        c("Timepoint:Profession", "Exposure:Profession")
    )
    expect_equal(nrow(fits$omnibusTests$Timepoint), 2L)
    expect_equal(nrow(fits$omnibusTests$Exposure), 2L)
    expect_identical(fits$settings$parallelBackend, "psock")
    expect_true(fits$settings$clusterReusedAcrossPhenotypes)

    summaries <- summarizeMethylationLMEModels(
        modelResults = fits, preparedData = prepared,
        padjmethod = "fdr", nCores = 1, logs = FALSE
    )
    annotated <- annotateMethylationLMESummaries(
        modelSummaries = summaries,
        annotationObject = example$annotation,
        annotationCols = "Name,chr,pos", logs = FALSE
    )
    columns <- colnames(annotated$data)
    timepoint_coefficients <- grep(
        "^Timepoint_.*_P[.]Value$", columns
    )
    timepoint_coefficients <- timepoint_coefficients[!grepl(
        "_Omnibus_", columns[timepoint_coefficients], fixed = TRUE
    )]
    exposure_coefficients <- grep(
        "^Exposure_.*_P[.]Value$", columns
    )
    exposure_coefficients <- exposure_coefficients[!grepl(
        "_Omnibus_", columns[exposure_coefficients], fixed = TRUE
    )]
    timepoint_raw <- match(
        "Timepoint_Profession_Omnibus_P.Value", columns
    )
    exposure_raw <- match(
        "Exposure_Profession_Omnibus_P.Value", columns
    )
    omnibus_details <- grep("_Omnibus_", columns, fixed = TRUE)
    omnibus_details <- setdiff(
        omnibus_details, c(timepoint_raw, exposure_raw)
    )
    expect_lt(
        max(c(timepoint_coefficients, timepoint_raw)),
        min(c(exposure_coefficients, exposure_raw))
    )
    expect_lt(
        max(c(exposure_coefficients, exposure_raw)),
        match("Name", columns)
    )
    expect_lt(match("pos", columns), min(omnibus_details))
    expect_lt(max(omnibus_details), match("Timepoint_Model.Message", columns))
})

test_that("omnibus results propagate to workbook and report sidecars", {
    example <- create_omnibus_lme_example(withr::local_tempdir())
    tmp <- withr::local_tempdir()
    file.copy(example$input, file.path(tmp, basename(example$input)))
    input <- file.path(tmp, basename(example$input))
    report_assets <- file.path(
        tmp, "reports", "model1", "assets", "results", "lme_results"
    )

    result <- suppressWarnings(methylationLME(
        inputPheno = input,
        outputLogs = file.path(tmp, "logs"),
        outputRData = file.path(tmp, "rData", "models"),
        outputPlots = file.path(tmp, "figures"),
        personVar = "person", timeVar = "Timepoint",
        phenotypes = "Timepoint", covariates = NULL,
        factorVars = "Timepoint,Profession",
        interactionTerm = "Profession",
        omnibusTest = TRUE, omnibusDdf = "Satterthwaite",
        cpgLimit = 2, nCores = 1, summaryPval = NA,
        significantInteractionDir = file.path(
            tmp, "preliminaryResults", "cpgs", "methylationLME"
        ),
        significantInteractionPval = 1,
        summaryTxtDir = file.path(
            tmp, "preliminaryResults", "summary", "methylationLME"
        ),
        annotationPackage = example$annotation,
        annotationCols = "Name,chr,pos",
        annotatedLMEOut = file.path(tmp, "data", "methylationLME"),
        reportAssetsDir = report_assets,
        display = FALSE, verbose = FALSE, logs = TRUE,
        saveOutputs = TRUE
    ))

    expect_true(result$modelFits$settings$omnibusTest)
    omnibus <- result$modelSummaries$omnibusTests$Timepoint
    expect_equal(nrow(omnibus), 2L)
    valid <- is.finite(omnibus$Omnibus.P.Value)
    expect_equal(
        omnibus$Omnibus.Adjusted.P.Value[valid],
        stats::p.adjust(omnibus$Omnibus.P.Value[valid], method = "fdr")
    )

    expected_omnibus <- paste0(
        "Timepoint_Profession_Omnibus_",
        c(
            "F.Value", "Num.DF", "Den.DF", "P.Value",
            "Adjusted.P.Value", "Method"
        )
    )
    workbook <- result$savedFiles$annotatedLME
    result_sheet <- openxlsx::read.xlsx(
        workbook, sheet = "annotatedLME"
    )
    columns <- colnames(result_sheet)
    expect_true(all(expected_omnibus %in% columns))
    coefficient_columns <- grep(
        "^Timepoint_.*:.*_P\\.Value$", columns,
        value = TRUE
    )
    expect_true(length(coefficient_columns) >= 2L)
    raw_omnibus <- "Timepoint_Profession_Omnibus_P.Value"
    detail_omnibus <- setdiff(expected_omnibus, raw_omnibus)
    expect_lt(
        max(match(coefficient_columns, columns)),
        match(raw_omnibus, columns)
    )
    expect_lt(
        match(raw_omnibus, columns),
        match("Name", columns)
    )
    expect_lt(
        match("pos", columns),
        min(match(detail_omnibus, columns))
    )
    expect_lt(
        max(match(detail_omnibus, columns)),
        match("Timepoint_Model.Message", columns)
    )
    expect_equal(
        openxlsx::getSheetNames(workbook),
        c("annotatedLME", "metadata", "dictionary")
    )
    dictionary <- openxlsx::read.xlsx(workbook, sheet = "dictionary")
    expect_true(any(grepl("joint fixed-effect omnibus", dictionary$Description)))
    expect_true(all(
        dictionary$Formula[dictionary$Column %in% expected_omnibus] != ""
    ))
    metadata <- openxlsx::read.xlsx(workbook, sheet = "metadata")
    metadata_values <- stats::setNames(metadata$Value, metadata$Key)
    expect_identical(metadata_values[["omnibus_test"]], "TRUE")
    expect_identical(metadata_values[["omnibus_joint"]], "TRUE")
    expect_identical(
        metadata_values[["random_effect_structure"]],
        "(1 | person)"
    )
    expect_identical(
        metadata_values[["omnibus_ddf"]],
        "Satterthwaite"
    )
    expect_identical(metadata_values[["p_adjust_method"]], "fdr")
    expect_true(file.exists(file.path(
        report_assets,
        "annotatedLME.tsv.gz"
    )))
    expect_false(dir.exists(file.path(
        dirname(workbook),
        "report-assets"
    )))
})
