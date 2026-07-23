test_that("methylation scales are normalized and ranges are read-only", {
    expect_identical(
        vapply(
            c("beta", "Beta", "BETA", "m", "M", "cn", "CN", "cN"),
            dnaEPICO:::normalizeMethylationScaleDnaEpico,
            character(1)
        ),
        c(
            beta = "beta", Beta = "beta", BETA = "beta",
            m = "m", M = "m", cn = "cn", CN = "cn", cN = "cn"
        )
    )
    expect_error(
        dnaEPICO:::normalizeMethylationScaleDnaEpico("unsupported"),
        "Beta, M, CN"
    )

    beta <- rbind(
        cg_boundary = c(0, 1),
        cg_nan = c(NaN, 0.5),
        cg_outside = c(-0.01, 1.01),
        cg_inf = c(Inf, 0.3)
    )
    original <- beta
    range <- dnaEPICO:::summarizeMethylationRangeDnaEpico(beta, "Beta")
    expect_identical(beta, original)
    expect_identical(
        names(range),
        c("Scale", "Observed.Minimum", "Observed.Maximum")
    )
    expect_identical(range$Scale, "beta")
    expect_equal(range$Observed.Minimum, -0.01)
    expect_equal(range$Observed.Maximum, 1.01)
    expect_false("Defined.Domain" %in% names(range))
})

test_that("methylation probe identifiers must be present and unique", {
    expect_identical(
        dnaEPICO:::validateMethylationProbeIdentifiersDnaEpico(
            c("cg00000001", "cg00000002")
        ),
        c("cg00000001", "cg00000002")
    )
    expect_error(
        dnaEPICO:::validateMethylationProbeIdentifiersDnaEpico(
            c("cg00000001", "cg00000001")
        ),
        "duplicate CpG identifiers"
    )
    expect_error(
        dnaEPICO:::validateMethylationProbeIdentifiersDnaEpico(c("cg00000001", "")),
        "missing or blank CpG identifiers"
    )
})

test_that("GLM passes numeric CpG values to glm2 without range rejection", {
    tmp <- withr::local_tempdir()
    phenoBT1 <- data.frame(
        Sample_Name = paste0("S", seq_len(12)),
        phenotype = rep(c(0, 1), 6),
        cg_valid = seq(0.1, 0.65, length.out = 12),
        cg_nan = c(NaN, seq(0.2, 0.7, length.out = 11)),
        cg_outside = c(1.1, seq(0.2, 0.7, length.out = 11)),
        check.names = FALSE
    )
    input_file <- file.path(tmp, "phenoBT1.RData")
    save(phenoBT1, file = input_file)

    prepared <- prepareMethylationGLMData(
        inputPheno = input_file,
        phenotypes = "phenotype",
        covariates = character(0),
        factorVars = character(0),
        methylationScale = "beta",
        verbose = FALSE,
        logs = FALSE
    )
    expect_true(is.nan(prepared$data$cg_nan[[1L]]))
    expect_equal(prepared$data$cg_outside[[1L]], 1.1)

    fits <- fitMethylationGLMModels(
        preparedData = prepared,
        nCores = 1,
        verbose = FALSE,
        logs = FALSE
    )
    summaries <- summarizeMethylationGLMModels(
        modelResults = fits,
        preparedData = prepared,
        verbose = FALSE,
        logs = FALSE
    )
    annotation <- annotateMethylationGLMSummaries(
        modelSummaries = summaries,
        annotationObject = data.frame(
            CpG = c("cg_valid", "cg_nan", "cg_outside"),
            chr = "chr1",
            stringsAsFactors = FALSE
        ),
        annotationCols = "chr",
        verbose = FALSE,
        logs = FALSE
    )

    expect_equal(nrow(fits$fitFailures), 0L)
    expect_setequal(
        summaries$summaries$phenotype$CpG,
        c("cg_valid", "cg_nan", "cg_outside")
    )
    expect_setequal(
        annotation$data$IlmnID,
        c("cg_valid", "cg_nan", "cg_outside")
    )
})

test_that("GLM retains p-values returned by glm2 for a constant response", {
    tmp <- withr::local_tempdir()
    phenoBT1 <- data.frame(
        Sample_Name = paste0("S", seq_len(8)),
        phenotype = rep(c(0, 1), 4),
        cg_constant = rep(1.2, 8),
        check.names = FALSE
    )
    input_file <- file.path(tmp, "phenoBT1.RData")
    save(phenoBT1, file = input_file)
    prepared <- prepareMethylationGLMData(
        inputPheno = input_file,
        phenotypes = "phenotype",
        covariates = character(0),
        factorVars = character(0),
        verbose = FALSE,
        logs = FALSE
    )

    fits <- suppressWarnings(
        fitMethylationGLMModels(prepared, nCores = 1, logs = FALSE)
    )
    expect_true(fits$modelMessages$P.Value.Available)
    expect_equal(fits$failureCounts[["phenotype"]], 0L)
})

test_that("lme4 and nlme do not reject numeric CpGs by methylation range", {
    testthat::skip_if_not_installed("lmerTest")
    testthat::skip_if_not_installed("nlme")

    tmp <- withr::local_tempdir()
    set.seed(11)
    phenoBT1T2 <- data.frame(
        SID = paste0("P", rep(seq_len(10), each = 2), c("A", "B")),
        person = rep(seq_len(10), each = 2),
        Timepoint = factor(rep(c("T1", "T2"), 10)),
        phenotype = rep(c(0, 1), each = 10),
        cg_valid = 0.3 + rep(c(0, 0.04), 10) + stats::rnorm(20, 0, 0.01),
        cg_outside = c(-0.1, 0.3 + rep(c(0, 0.04), 9), 0.36),
        check.names = FALSE
    )
    input_file <- file.path(tmp, "phenoBT1T2.RData")
    save(phenoBT1T2, file = input_file)
    prepared <- prepareMethylationLMEData(
        inputPheno = input_file,
        personVar = "person",
        timeVar = "Timepoint",
        phenotypes = "phenotype",
        covariates = character(0),
        factorVars = "phenotype",
        verbose = FALSE,
        logs = FALSE
    )

    lme4_fits <- fitMethylationLMEModels(
        preparedData = prepared,
        lmeLibs = "lme4,lmerTest",
        nCores = 1,
        logs = FALSE
    )
    nlme_fits <- fitMethylationLMEModels(
        preparedData = prepared,
        lmeLibs = "nlme",
        nCores = 1,
        logs = FALSE
    )

    expect_true("cg_outside" %in%
        lme4_fits$coefficientResults$phenotype$cpgOrder)
    expect_true("cg_outside" %in%
        nlme_fits$coefficientResults$phenotype$cpgOrder)
    expect_true(any(is.finite(
        lme4_fits$coefficientResults$phenotype$estimate["cg_outside", ]
    )))
    expect_true(any(is.finite(
        nlme_fits$coefficientResults$phenotype$estimate["cg_outside", ]
    )))
})

test_that("model preparation selects scale-specific data from multi-object external files", {
    tmp <- withr::local_tempdir()
    phenoMT1 <- data.frame(
        Sample_Name = paste0("S", seq_len(8)),
        phenotype = rep(c(0, 1), 4),
        cg_external = seq(-4, 3),
        check.names = FALSE
    )
    distractor <- data.frame(not_the_model_input = 1)
    glm_file <- file.path(tmp, "phenoMT1.RData")
    save(phenoMT1, distractor, file = glm_file)

    glm_data <- prepareMethylationGLMData(
        inputPheno = glm_file,
        phenotypes = "phenotype",
        covariates = character(0),
        factorVars = character(0),
        methylationScale = "m",
        verbose = FALSE,
        logs = FALSE
    )
    expect_equal(glm_data$data$cg_external, phenoMT1$cg_external)
    expect_identical(glm_data$methylationScale, "m")

    phenoCNT1T2 <- data.frame(
        SID = paste0("P", rep(seq_len(4), each = 2), rep(c("A", "B"), 4)),
        person = rep(seq_len(4), each = 2),
        Timepoint = rep(c(1, 2), 4),
        phenotype = rep(c(0, 1), 4),
        cg_external = seq(-20, 15, length.out = 8),
        check.names = FALSE
    )
    lme_file <- file.path(tmp, "phenoCNT1T2.RData")
    save(phenoCNT1T2, distractor, file = lme_file)

    lme_data <- prepareMethylationLMEData(
        inputPheno = lme_file,
        personVar = "person",
        timeVar = "Timepoint",
        phenotypes = "phenotype",
        covariates = character(0),
        factorVars = character(0),
        methylationScale = "cn",
        verbose = FALSE,
        logs = FALSE
    )
    expect_identical(lme_data$data$cg_external, phenoCNT1T2$cg_external)
    expect_identical(lme_data$methylationScale, "cn")
})
