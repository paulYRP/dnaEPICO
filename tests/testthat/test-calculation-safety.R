test_that("sample identifiers are reordered explicitly and duplicates fail", {
    expect_identical(
        dnaEPICO:::matchSampleIdentifiersDnaEpico(
            query = c("S2", "S1"),
            reference = c("S1", "S2"),
            requireSameSet = TRUE
        ),
        c(2L, 1L)
    )

    expect_error(
        dnaEPICO:::matchSampleIdentifiersDnaEpico(
            query = c("S1", "S1"),
            reference = c("S1", "S2")
        ),
        "duplicate sample identifiers"
    )
    expect_error(
        dnaEPICO:::matchSampleIdentifiersDnaEpico(
            query = c("S1", "S3"),
            reference = c("S1", "S2")
        ),
        "S3"
    )
})

test_that("sex labels are canonicalized without changing reported values", {
    original <- c("Female", "M", "0", "1", "unknown", NA)
    result <- dnaEPICO:::canonicalizeSexDnaEpico(original)

    expect_identical(result$original, original)
    expect_identical(result$code, c(0L, 1L, 0L, 1L, NA_integer_, NA_integer_))
    expect_identical(result$label, c("F", "M", "F", "M", NA, NA))
    expect_identical(result$unknown, "unknown")
    expect_identical(formals(preprocessingMinfiEwasWater)$removeSexMismatch, FALSE)

    expect_identical(
        dnaEPICO:::resolveNormalizationSexDnaEpico(
            data.frame(
                Sex = c("F", "unknown", NA),
                PredSex = c(0L, 1L, 0L)
            ),
            "Sex"
        ),
        c("F", "M", "F")
    )
})

test_that("sex mismatch removal keeps RGSet and phenotype rows aligned", {
    rgset <- matrix(seq_len(12), nrow = 3)
    colnames(rgset) <- c("S3", "S1", "S2", "S4")
    targets <- data.frame(
        Sample_Name = c("S1", "S2", "S3", "S4"),
        Sex = c("F", "M", "F", "M"),
        stringsAsFactors = FALSE
    )

    result <- filterSamplesMinfiEwasWater(
        RGSet = rgset,
        targets = targets,
        failedSamples = c("S2", "S4"),
        verbose = FALSE,
        logs = FALSE
    )

    expect_identical(colnames(result$RGSet), c("S3", "S1"))
    expect_identical(result$targets$Sample_Name, c("S3", "S1"))
})

test_that("coefficient extraction follows formula terms exactly", {
    model_data <- data.frame(
        status = factor(rep(c("Control", "Case"), each = 4)),
        sex = factor(rep(c("F", "M"), times = 4))
    )
    formula_text <- "beta ~ status * sex"
    term_map <- dnaEPICO:::buildCoefficientTermMapMethylationModels(
        formulaText = formula_text,
        data = model_data
    )

    interaction_rows <- dnaEPICO:::findCoefficientRowsMethylationGLM(
        coefNames = names(term_map),
        variable = "status",
        interactionTerm = "sex",
        coefficientTerms = term_map
    )
    main_rows <- dnaEPICO:::findCoefficientRowsMethylationGLM(
        coefNames = names(term_map),
        variable = "status",
        coefficientTerms = term_map
    )

    expect_true(length(interaction_rows) > 0L)
    expect_true(all(term_map[interaction_rows] == "status:sex"))
    expect_true(length(main_rows) > 0L)
    expect_true(all(term_map[main_rows] == "status"))
    expect_false(any(interaction_rows %in% main_rows))
})

test_that("rank-deficient fixed-effect designs fail before CpG fitting", {
    model_data <- data.frame(x = 1:8, duplicate_x = 1:8)
    expect_error(
        dnaEPICO:::validateFixedEffectDesignMethylationModels(
            formulaText = "beta ~ x + duplicate_x",
            data = model_data
        ),
        "rank deficient"
    )
    expect_error(
        dnaEPICO:::validateFixedEffectDesignMethylationModels(
            formulaText = "beta ~ x",
            data = data.frame(x = c(1, Inf, 3))
        ),
        "non-finite numeric values"
    )
})

test_that("a phenotype cannot be used as its own interaction term", {
    expect_error(
        dnaEPICO:::buildFormulaMethylationGLM(
            phenotype = "status",
            interactionTerm = "status"
        ),
        "must differ"
    )
    expect_error(
        dnaEPICO:::buildFormulaMethylationLME(
            phenotype = "score",
            personVar = "person",
            interactionTerm = "score"
        ),
        "must differ"
    )
})

test_that("phenotype-only GLM and LME formulas contain no empty term", {
    expect_identical(
        dnaEPICO:::buildFormulaMethylationGLM(
            phenotype = "status",
            covariates = character(0)
        ),
        "beta ~ `status`"
    )
    expect_identical(
        dnaEPICO:::buildFormulaMethylationLME(
            phenotype = "score",
            personVar = "person",
            covariates = character(0)
        ),
        "beta ~ `score` + (1 | `person` )"
    )
})

test_that("longitudinal subject derivation accepts only explicit A/B visits", {
    valid <- data.frame(SID = c("P1A", "P1B", "P2A", "P2B"))
    derived <- dnaEPICO:::ensurePersonColumnMethylationLME(valid)
    expect_identical(derived$data$person, c("P1", "P1", "P2", "P2"))

    invalid <- data.frame(SID = c("P1_visit1", "P1_visit2", "P2_visit1"))
    expect_error(
        dnaEPICO:::ensurePersonColumnMethylationLME(invalid),
        "Supply an explicit subject identifier"
    )
})

test_that("probability and CpG-limit options reject invalid values", {
    expect_equal(
        dnaEPICO:::validateProbabilityDnaEpico(0.05, "threshold"),
        0.05
    )
    expect_error(
        dnaEPICO:::validateProbabilityDnaEpico(1.2, "threshold"),
        "between 0 and 1"
    )
    expect_error(
        dnaEPICO:::validateCpgLimitMethylationModels(2.5),
        "positive whole number"
    )
    expect_error(
        dnaEPICO:::validateCpgPrefixDnaEpico(""),
        "cpgPrefix"
    )
    expect_error(
        dnaEPICO:::normalizeOptionalNumericMethylationGLM(c(1, 2)),
        "one value"
    )
    expect_error(
        dnaEPICO:::normalizeOptionalNumericMethylationGLM("not-a-number"),
        "must be numeric"
    )
})

test_that("phenotype-to-PRS mappings are unambiguous", {
    expect_identical(
        dnaEPICO:::parsePrsMapMethylationGLM("T1D:PRS_T1D,HbA1c:PRS_HbA1c"),
        c(T1D = "PRS_T1D", HbA1c = "PRS_HbA1c")
    )
    expect_error(
        dnaEPICO:::parsePrsMapMethylationGLM("T1D:PRS1,T1D:PRS2"),
        "duplicate phenotype mappings"
    )
    expect_error(
        dnaEPICO:::parsePrsMapMethylationGLM("T1D:"),
        "format 'Phenotype:PRS'"
    )
    expect_error(
        dnaEPICO:::parsePrsMapMethylationGLM("T1D:PRS:extra"),
        "format 'Phenotype:PRS'"
    )
})

test_that("constant CpG responses are reported as failed GLM fits", {
    result <- dnaEPICO:::fitCpGModelMethylationGLM(
        cpg = "cg00000001",
        cpgValues = rep(0.5, 8),
        modelData = data.frame(status = rep(c(0, 1), each = 4)),
        formulaText = "beta ~ status"
    )

    expect_s3_class(result, "dnaEPICO_methylationGLM_fit_error")
    expect_match(result$error, "no observed variation")
})

test_that("constant CpG responses are reported before mixed-model fitting", {
    result <- dnaEPICO:::fitCpGModelMethylationLME(
        cpg = "cg00000001",
        cpgValues = rep(0.5, 8),
        modelData = data.frame(
            person = factor(rep(seq_len(4), each = 2)),
            score = rep(c(0, 1), 4)
        ),
        formulaText = "beta ~ score + (1 | person)",
        personVar = "person",
        lmeEngine = "lme4"
    )

    expect_s3_class(result, "dnaEPICO_methylationLME_fit_error")
    expect_match(result$error, "no observed variation")
})

test_that("external annotation identifiers cannot be implicit or duplicated", {
    expect_error(
        dnaEPICO:::coerceAnnotationDataMethylationGLM(data.frame(chr = "chr1")),
        "must include a CpG or IlmnID column"
    )
    expect_error(
        dnaEPICO:::coerceAnnotationDataMethylationGLM(data.frame(
            CpG = c("cg00000001", "cg00000001"),
            chr = c("chr1", "chr1")
        )),
        "duplicate CpG identifiers"
    )
})

test_that("finite summaries return NA instead of creating NaN", {
    expect_equal(dnaEPICO:::meanFiniteOrNADnaEpico(c(NA_real_, NaN, Inf)), NA_real_)
    expect_equal(dnaEPICO:::meanFiniteOrNADnaEpico(c(NA_real_, 0.2, 0.4)), 0.3)

    timepoint_summary <- summarizeTimepointsMethylationLME(
        data = data.frame(
            Timepoint = c("T1", "T1", "T2", "T2"),
            score = c(NA_real_, NaN, 1, Inf)
        ),
        timeVar = "Timepoint",
        phenotypes = "score",
        logs = FALSE
    )
    expect_true(is.na(timepoint_summary$score_mean[timepoint_summary$Timepoint == "T1"]))
    expect_false(any(is.nan(timepoint_summary$score_mean)))
    expect_identical(timepoint_summary$score_n, c(0L, 1L))
})

test_that("diagnostic means stay on the declared methylation scale", {
    prepared <- list(
        data = data.frame(
            cg_boundary = c(0, 0),
            cg_valid = c(0.2, 0.4),
            cg_invalid = c(Inf, NA_real_)
        ),
        cpgColumns = c("cg_boundary", "cg_valid", "cg_invalid"),
        methylationScale = "beta",
        responseLabel = "Beta values",
        invalidCpGs = data.frame(CpG = "cg_invalid")
    )
    result <- dnaEPICO:::diagnosticMeanMethylationModels(prepared)

    expect_identical(result$label, "Average Beta")
    expect_equal(unname(result$values[c("cg_boundary", "cg_valid")]), c(0, 0.3))
    expect_true(is.na(result$values[["cg_invalid"]]))
    expect_false(any(is.nan(result$values)))
})

test_that("saved plot devices close before an optional display draw", {
    tmp <- withr::local_tempdir()
    plot_file <- file.path(tmp, "device-test.tiff")
    draws <- 0L
    grDevices::pdf(file.path(tmp, "display-device.pdf"))
    on.exit(grDevices::dev.off(), add = TRUE)
    device_before <- grDevices::dev.cur()

    dnaEPICO:::runPlotMinfiEwasWater(
        draw_fun = function() {
            draws <<- draws + 1L
            graphics::plot.new()
            graphics::text(0.5, 0.5, "device test")
        },
        display = TRUE,
        file = plot_file,
        width = 400,
        height = 400,
        res = 72
    )

    expect_equal(draws, 2L)
    expect_true(file.exists(plot_file))
    expect_identical(grDevices::dev.cur(), device_before)
})

test_that("GLM filtering does not truncate diagnostics or annotated CpGs", {
    summary_cache <- data.frame(
        CpG = rep(c("cg1", "cg2"), each = 2L),
        Coefficient = rep(c("statusB", "statusC"), times = 2L),
        Estimate = c(0.2, -0.1, 0.05, 0.08),
        `Std. Error` = rep(0.05, 4L),
        `t value` = c(4, -2, 1, 1.6),
        `Pr(>|t|)` = c(0.001, 0.04, 0.5, 0.2),
        ResidualSD = rep(0.1, 4L),
        check.names = FALSE
    )
    prepared <- list(
        data = data.frame(cg1 = c(0.2, 0.3), cg2 = c(0.4, 0.5)),
        cpgColumns = c("cg1", "cg2"),
        methylationScale = "beta",
        responseLabel = "Beta values",
        interactionTerm = NULL,
        invalidCpGs = data.frame(CpG = character(0))
    )
    summaries <- summarizeMethylationGLMModels(
        modelResults = list(
            fits = list(status = list()),
            summaryCache = list(status = summary_cache),
            fitFailures = data.frame()
        ),
        preparedData = prepared,
        summaryPval = 0.05,
        logs = FALSE
    )

    expect_equal(nrow(summaries$summaries$status), 2L)
    expect_equal(nrow(summaries$diagnosticSummaries$status), 4L)

    diagnostics <- plotMethylationGLMDiagnostics(
        modelSummaries = summaries,
        preparedData = prepared,
        display = FALSE,
        logs = FALSE
    )
    expect_setequal(names(diagnostics$plots$status), c("statusB", "statusC"))
    expect_setequal(names(diagnostics$inflationFactors$status), c("statusB", "statusC"))

    annotation <- annotateMethylationGLMSummaries(
        modelSummaries = summaries,
        annotationObject = data.frame(CpG = c("cg1", "cg2"), chr = c("chr1", "chr2")),
        annotationCols = "chr",
        logs = FALSE
    )
    expect_setequal(annotation$data$IlmnID, c("cg1", "cg2"))
})

test_that("LME filtering does not truncate diagnostics or annotated CpGs", {
    summary_cache <- data.frame(
        CpG = rep(c("cg1", "cg2"), each = 2L),
        Interaction.Term = rep(c("score", "score:visit"), times = 2L),
        Estimate = c(0.2, -0.1, 0.05, 0.08),
        Std.Error = rep(0.05, 4L),
        t.value = c(4, -2, 1, 1.6),
        P.value = c(0.001, 0.04, 0.5, 0.2),
        stringsAsFactors = FALSE
    )
    prepared <- list(
        data = data.frame(cg1 = c(0.2, 0.3), cg2 = c(0.4, 0.5)),
        cpgColumns = c("cg1", "cg2"),
        methylationScale = "beta",
        responseLabel = "Beta values",
        interactionTerm = "visit",
        invalidCpGs = data.frame(CpG = character(0))
    )
    summaries <- summarizeMethylationLMEModels(
        modelResults = list(
            fits = list(score = list()),
            summaryCache = list(score = summary_cache),
            fitFailures = data.frame()
        ),
        preparedData = prepared,
        summaryPval = 0.05,
        logs = FALSE
    )

    expect_equal(nrow(summaries$summaries$score), 2L)
    expect_equal(nrow(summaries$diagnosticSummaries$score), 4L)

    diagnostics <- plotMethylationLMEDiagnostics(
        modelSummaries = summaries,
        preparedData = prepared,
        display = FALSE,
        logs = FALSE
    )
    expect_setequal(names(diagnostics$plots$score), c("score", "score:visit"))
    expect_setequal(names(diagnostics$inflationFactors$score), c("score", "score:visit"))

    annotation <- annotateMethylationLMESummaries(
        modelSummaries = summaries,
        annotationObject = data.frame(CpG = c("cg1", "cg2"), chr = c("chr1", "chr2")),
        annotationCols = "chr",
        logs = FALSE
    )
    expect_setequal(annotation$data$IlmnID, c("cg1", "cg2"))
})
