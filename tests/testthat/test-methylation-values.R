test_that("methylation-scale validation follows Beta, M, and CN definitions", {
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
        cg_low = c(-0.01, 0.2),
        cg_inf = c(Inf, 0.3),
        cg_missing = c(NA_real_, NA_real_)
    )
    beta_result <- dnaEPICO:::inspectMethylationMatrixDnaEpico(beta, "beta")

    expect_true(is.na(beta_result$values["cg_nan", 1]))
    expect_equal(beta_result$boundaries$Observed.Minimum, -0.01)
    expect_equal(beta_result$boundaries$Observed.Maximum, 1)
    expect_equal(beta_result$boundaries$At.Lower.Boundary, 1L)
    expect_equal(beta_result$boundaries$At.Upper.Boundary, 1L)
    expect_equal(beta_result$boundaries$NaN.Converted, 1L)
    expect_equal(beta_result$boundaries$Invalid.CpGs, 3L)
    expect_setequal(
        beta_result$invalidCpGs$CpG,
        c("cg_low", "cg_inf", "cg_missing")
    )
    expect_identical(
        beta_result$issues$Status[beta_result$issues$CpG == "cg_nan"],
        "NaN converted to NA"
    )

    m_values <- rbind(cg_finite = c(-100, 100), cg_inf = c(-Inf, 0))
    m_result <- dnaEPICO:::inspectMethylationMatrixDnaEpico(m_values, "m")
    expect_equal(m_result$boundaries$Out.Of.Range, 0L)
    expect_identical(m_result$invalidCpGs$CpG, "cg_inf")

    cn_values <- rbind(cg_finite = c(-20, 30), cg_inf = c(Inf, 0))
    cn_result <- dnaEPICO:::inspectMethylationMatrixDnaEpico(cn_values, "cn")
    expect_equal(cn_result$boundaries$Out.Of.Range, 0L)
    expect_identical(cn_result$invalidCpGs$CpG, "cg_inf")
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

test_that("GLM skips and reports invalid CpGs without stopping valid fits", {
    tmp <- withr::local_tempdir()
    phenoBT1 <- data.frame(
        Sample_Name = paste0("S", seq_len(12)),
        phenotype = rep(c(0, 1), 6),
        cg_valid = seq(0.1, 0.65, length.out = 12),
        cg_nan = c(NaN, seq(0.2, 0.7, length.out = 11)),
        cg_invalid = c(1.1, seq(0.2, 0.7, length.out = 11)),
        check.names = FALSE
    )
    input_file <- file.path(tmp, "phenoBT1.RData")
    save(phenoBT1, file = input_file)

    prepared <- suppressWarnings(prepareMethylationGLMData(
        inputPheno = input_file,
        phenotypes = "phenotype",
        covariates = character(0),
        factorVars = character(0),
        methylationScale = "beta",
        verbose = FALSE,
        logs = FALSE
    ))
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
            CpG = c("cg_valid", "cg_nan", "cg_invalid"),
            chr = "chr1",
            stringsAsFactors = FALSE
        ),
        annotationCols = "chr",
        verbose = FALSE,
        logs = FALSE
    )

    expect_identical(fits$fitFailures$CpG, "cg_invalid")
    expect_identical(fits$fitFailures$Status, "invalid")
    expect_match(fits$fitFailures$Error, "outside \\[0, 1\\]")
    expect_setequal(summaries$summaries$phenotype$CpG, c("cg_valid", "cg_nan"))
    expect_setequal(annotation$data$IlmnID, c("cg_valid", "cg_nan", "cg_invalid"))
    expect_identical(
        annotation$data$phenotype_Fit.Status[
            annotation$data$IlmnID == "cg_invalid"
        ],
        "invalid"
    )
    expect_match(
        annotation$data$phenotype_Exclusion.Reason[
            annotation$data$IlmnID == "cg_invalid"
        ],
        "outside \\[0, 1\\]"
    )
})

test_that("GLM retains an inventory when every CpG is invalid", {
    tmp <- withr::local_tempdir()
    phenoBT1 <- data.frame(
        Sample_Name = paste0("S", seq_len(8)),
        phenotype = rep(c(0, 1), 4),
        cg_invalid = rep(1.2, 8),
        check.names = FALSE
    )
    input_file <- file.path(tmp, "phenoBT1.RData")
    save(phenoBT1, file = input_file)
    prepared <- suppressWarnings(prepareMethylationGLMData(
        inputPheno = input_file,
        phenotypes = "phenotype",
        covariates = character(0),
        factorVars = character(0),
        verbose = FALSE,
        logs = FALSE
    ))

    expect_warning(
        fits <- fitMethylationGLMModels(prepared, nCores = 1, logs = FALSE),
        "failure inventory was retained"
    )
    expect_equal(nrow(fits$fitFailures), 1L)
    expect_identical(fits$fitFailures$Status, "invalid")
})

test_that("lme4 and nlme skip the same invalid longitudinal CpG", {
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
        cg_invalid = c(-0.1, rep(0.4, 19)),
        check.names = FALSE
    )
    input_file <- file.path(tmp, "phenoBT1T2.RData")
    save(phenoBT1T2, file = input_file)
    prepared <- suppressWarnings(prepareMethylationLMEData(
        inputPheno = input_file,
        personVar = "person",
        timeVar = "Timepoint",
        phenotypes = "phenotype",
        covariates = character(0),
        factorVars = "phenotype",
        verbose = FALSE,
        logs = FALSE
    ))

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

    expect_identical(lme4_fits$fitFailures$CpG, "cg_invalid")
    expect_identical(nlme_fits$fitFailures$CpG, "cg_invalid")
    expect_identical(lme4_fits$fitFailures$Status, "invalid")
    expect_identical(nlme_fits$fitFailures$Status, "invalid")
    expect_false(inherits(lme4_fits$fits$phenotype$cg_valid, "dnaEPICO_methylationLME_fit_error"))
    expect_false(inherits(nlme_fits$fits$phenotype$cg_valid, "dnaEPICO_methylationLME_fit_error"))
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
