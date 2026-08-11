create_preprocessing_pheno_example <- function(path) {
    pheno <- data.frame(
        Sample_Name = c("S1", "S2", "S3"),
        Timepoint = c("1", "1", "2"),
        Sex = c(0, 1, 0),
        stringsAsFactors = FALSE
    )

    beta <- matrix(
        c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60),
        nrow = 2,
        dimnames = list(c("cg1", "cg2"), pheno$Sample_Name)
    )
    m <- beta * 10
    cn <- beta * 100

    pheno_file <- file.path(path, "pheno.csv")
    beta_path <- file.path(path, "beta.RData")
    m_path <- file.path(path, "m.RData")
    cn_path <- file.path(path, "cn.RData")

    utils::write.csv(pheno, pheno_file, row.names = FALSE)
    save(beta, file = beta_path)
    save(m, file = m_path)
    save(cn, file = cn_path)

    list(
        pheno = pheno,
        beta = beta,
        m = m,
        cn = cn,
        phenoFile = pheno_file,
        betaPath = beta_path,
        mPath = m_path,
        cnPath = cn_path
    )
}

test_that("preprocessingPheno returns in-memory objects quietly by default", {
    tmp <- withr::local_tempdir()
    example_data <- create_preprocessing_pheno_example(tmp)

    expect_message(
        result <- preprocessingPheno(
            phenoFile = example_data$phenoFile,
            betaPath = example_data$betaPath,
            mPath = example_data$mPath,
            cnPath = example_data$cnPath,
            SampleID = "Sample_Name",
            timeVar = "Timepoint",
            timepoints = "1,2",
            combineTimepoints = "1,2",
            outputPheno = file.path(tmp, "data", "preprocessingPheno"),
            outputRData = file.path(tmp, "rData", "preprocessingPheno", "metrics"),
            outputRDataMerge = file.path(tmp, "rData", "preprocessingPheno", "mergeData"),
            sexColumn = "Sex",
            outputLogs = file.path(tmp, "logs"),
            outputDir = file.path(tmp, "clockFoundation"),
            saveOutputs = FALSE
        ),
        NA
    )

    expect_s3_class(result, "dnaEPICO_preprocessingPheno")
    expect_null(result$savedFiles)
    expect_identical(result$timepointData$timepoints, c("1", "2"))
    expect_equal(nrow(result$combinedData$pheno), 3)
    expect_equal(nrow(result$combinedData$phenoB), 3)
    expect_true(all(c("ProbeID", "S1", "S2", "S3") %in% colnames(result$clockFoundation$betaCSV)))
    expect_true("id" %in% colnames(result$clockFoundation$phenoCF))
    expect_false(dir.exists(file.path(tmp, "data", "preprocessingPheno")))
})

test_that("preprocessingPheno can write legacy outputs and logs on request", {
    tmp <- withr::local_tempdir()
    example_data <- create_preprocessing_pheno_example(tmp)

    expect_message(
        result <- preprocessingPheno(
            phenoFile = example_data$phenoFile,
            betaPath = example_data$betaPath,
            mPath = example_data$mPath,
            cnPath = example_data$cnPath,
            SampleID = "Sample_Name",
            timeVar = "Timepoint",
            timepoints = "1,2",
            combineTimepoints = "1,2",
            outputPheno = file.path(tmp, "data", "preprocessingPheno"),
            outputRData = file.path(tmp, "rData", "preprocessingPheno", "metrics"),
            outputRDataMerge = file.path(tmp, "rData", "preprocessingPheno", "mergeData"),
            sexColumn = "Sex",
            outputLogs = file.path(tmp, "logs"),
            outputDir = file.path(tmp, "clockFoundation"),
            verbose = TRUE,
            logs = TRUE,
            saveOutputs = TRUE
        ),
        "Starting preprocessingPheno"
    )

    expect_s3_class(result$savedFiles, "dnaEPICO_preprocessingPheno_paths")
    expect_true(file.exists(file.path(tmp, "logs", "log_preprocessingPheno.txt")))
    expect_true(file.exists(file.path(tmp, "data", "preprocessingPheno", "phenoT1.csv")))
    expect_true(file.exists(file.path(tmp, "rData", "preprocessingPheno", "metrics", "betaT1.RData")))
    pheno_b_t1 <- file.path(
        tmp,
        "rData",
        "preprocessingPheno",
        "mergeData",
        "phenoBT1.RData"
    )
    expect_true(file.exists(pheno_b_t1))
    expect_equal(
        result$savedFiles$combinedPhenoB,
        file.path(
            tmp,
            "rData",
            "preprocessingPheno",
            "mergeData",
            "phenoBT1T2.RData"
        )
    )
    expect_true(file.exists(file.path(tmp, "clockFoundation", "beta.csv")))
    expect_true(file.exists(file.path(tmp, "clockFoundation", "phenoCF.csv")))
})

test_that("preprocessingPheno skips combined outputs when requested", {
    tmp <- withr::local_tempdir()
    example_data <- create_preprocessing_pheno_example(tmp)

    result <- preprocessingPheno(
        phenoFile = example_data$phenoFile,
        betaPath = example_data$betaPath,
        mPath = example_data$mPath,
        cnPath = example_data$cnPath,
        SampleID = "Sample_Name",
        timeVar = "Timepoint",
        timepoints = "1,2",
        combineTimepoints = NULL,
        outputPheno = file.path(tmp, "data", "preprocessingPheno"),
        outputRData = file.path(tmp, "rData", "preprocessingPheno", "metrics"),
        outputRDataMerge = file.path(tmp, "rData", "preprocessingPheno", "mergeData"),
        sexColumn = "Sex",
        outputLogs = file.path(tmp, "logs"),
        outputDir = file.path(tmp, "clockFoundation"),
        saveOutputs = TRUE
    )

    expect_true("combinedData" %in% names(result))
    expect_null(result$combinedData)
    expect_identical(result$timepointData$timepoints, c("1", "2"))
    expect_true(file.exists(file.path(
        tmp, "data", "preprocessingPheno", "phenoT1.csv"
    )))
    expect_true(file.exists(file.path(
        tmp, "data", "preprocessingPheno", "phenoT2.csv"
    )))
    expect_true(file.exists(file.path(
        tmp, "rData", "preprocessingPheno", "mergeData", "phenoBT1.RData"
    )))
    expect_true(file.exists(file.path(
        tmp, "rData", "preprocessingPheno", "mergeData", "phenoBT2.RData"
    )))
    expect_false(file.exists(file.path(
        tmp, "data", "preprocessingPheno", "phenoT1T2.csv"
    )))
    expect_false(file.exists(file.path(
        tmp, "rData", "preprocessingPheno", "mergeData", "phenoBT1T2.RData"
    )))
    expect_true(all(c(
        "combinedPheno", "combinedPhenoMethylation", "combinedPhenoB"
    ) %in% names(result$savedFiles)))
    expect_null(result$savedFiles$combinedPheno)
    expect_null(result$savedFiles$combinedPhenoMethylation)
    expect_null(result$savedFiles$combinedPhenoB)
    expect_true(file.exists(file.path(tmp, "clockFoundation", "beta.csv")))
    expect_true(file.exists(file.path(tmp, "clockFoundation", "phenoCF.csv")))
})

test_that("preprocessingPheno selects M and CN modeling scales without changing Clock Foundation beta", {
    tmp <- withr::local_tempdir()
    example_data <- create_preprocessing_pheno_example(tmp)

    result_m <- preprocessingPheno(
        phenoFile = example_data$phenoFile,
        betaPath = example_data$betaPath,
        mPath = example_data$mPath,
        cnPath = example_data$cnPath,
        SampleID = "Sample_Name",
        timeVar = "Timepoint",
        timepoints = "1,2",
        combineTimepoints = "1,2",
        methylationScale = "m",
        outputPheno = file.path(tmp, "data", "preprocessingPhenoM"),
        outputRData = file.path(tmp, "rData", "preprocessingPhenoM", "metrics"),
        outputRDataMerge = file.path(tmp, "rData", "preprocessingPhenoM", "mergeData"),
        sexColumn = "Sex",
        outputLogs = file.path(tmp, "logs"),
        outputDir = file.path(tmp, "clockFoundationM"),
        saveOutputs = FALSE
    )

    expect_equal(result_m$methylationScale, "m")
    expect_true("phenoM" %in% names(result_m$combinedData))
    expect_equal(result_m$combinedData$phenoM$cg1, c(1, 3, 5))
    expect_equal(
        result_m$clockFoundation$betaCSV$S1,
        c(0.10, 0.20)
    )

    result_cn <- preprocessingPheno(
        phenoFile = example_data$phenoFile,
        betaPath = example_data$betaPath,
        mPath = example_data$mPath,
        cnPath = example_data$cnPath,
        SampleID = "Sample_Name",
        timeVar = "Timepoint",
        timepoints = "1,2",
        combineTimepoints = "1,2",
        methylationScale = "cn",
        outputPheno = file.path(tmp, "data", "preprocessingPhenoCN"),
        outputRData = file.path(tmp, "rData", "preprocessingPhenoCN", "metrics"),
        outputRDataMerge = file.path(tmp, "rData", "preprocessingPhenoCN", "mergeData"),
        sexColumn = "Sex",
        outputLogs = file.path(tmp, "logs"),
        outputDir = file.path(tmp, "clockFoundationCN"),
        saveOutputs = TRUE
    )

    pheno_cn_t1 <- file.path(tmp, "rData", "preprocessingPhenoCN", "mergeData", "phenoCNT1.RData")
    pheno_cn_t1t2 <- file.path(tmp, "rData", "preprocessingPhenoCN", "mergeData", "phenoCNT1T2.RData")
    expect_true(file.exists(pheno_cn_t1))
    expect_true(file.exists(pheno_cn_t1t2))
    merge_files <- basename(list.files(file.path(tmp, "rData", "preprocessingPhenoCN", "mergeData")))
    expect_false("phenoCnT1.RData" %in% merge_files)
    expect_equal(result_cn$savedFiles$combinedPhenoCN, pheno_cn_t1t2)

    loaded <- new.env(parent = emptyenv())
    load(pheno_cn_t1, envir = loaded)
    expect_true("phenoCNT1" %in% ls(loaded))
    expect_equal(loaded$phenoCNT1$cg1, c(10, 30))
})

test_that("Clock Foundation input preserves numeric beta values", {
    beta <- rbind(
        cg_boundary = c(0, 1),
        cg_nan = c(NaN, 0.4),
        cg_outside = c(-0.01, 0.5),
        cg_missing = c(NA_real_, NA_real_)
    )
    colnames(beta) <- c("S1", "S2")
    pheno <- data.frame(
        Sample_Name = c("S2", "S1"),
        Sex = c("M", "F"),
        stringsAsFactors = FALSE
    )

    result <- buildClockFoundationInputsPreprocessingPheno(
        beta = beta,
        pheno = pheno,
        verbose = FALSE,
        logs = FALSE
    )

    expect_identical(
        result$betaCSV$ProbeID,
        c("cg_boundary", "cg_nan", "cg_outside", "cg_missing")
    )
    expect_true(is.nan(result$betaCSV$S1[result$betaCSV$ProbeID == "cg_nan"]))
    expect_identical(result$phenoCF$id, c("S1", "S2"))
    expect_equal(result$methylationRange$Observed.Minimum, -0.01)
    expect_equal(result$methylationRange$Observed.Maximum, 1)
})

test_that("Clock Foundation preserves reported sex values without fallback", {
    sample_ids <- paste0("S", seq_len(6L))
    beta <- matrix(
        seq(0.1, 0.6, length.out = 12L),
        nrow = 2L,
        dimnames = list(c("cg1", "cg2"), sample_ids)
    )
    pheno <- data.frame(
        Sample_Name = rev(sample_ids),
        Gender = c(NA, "", "Unknown", "not recorded", "F", "1"),
        PredSex = c(0L, 1L, 0L, 1L, 1L, 0L),
        stringsAsFactors = FALSE
    )
    expected <- pheno$Gender[match(sample_ids, pheno$Sample_Name)]
    expected_prediction <- pheno$PredSex[match(sample_ids, pheno$Sample_Name)]

    result <- buildClockFoundationInputsPreprocessingPheno(
        beta = beta,
        pheno = pheno,
        SampleID = "Sample_Name",
        sexColumn = "Gender",
        verbose = FALSE,
        logs = FALSE
    )

    expect_identical(result$phenoCF$id, sample_ids)
    expect_identical(result$phenoCF$Gender, expected)
    expect_identical(result$phenoCF$PredSex, expected_prediction)
})

test_that("Clock Foundation input rejects ambiguous CpG identifiers", {
    beta <- matrix(
        c(0.1, 0.2, 0.3, 0.4),
        nrow = 2,
        dimnames = list(c("cg1", "cg1"), c("S1", "S2"))
    )
    pheno <- data.frame(
        Sample_Name = c("S1", "S2"),
        Sex = c("F", "M"),
        stringsAsFactors = FALSE
    )

    expect_error(
        buildClockFoundationInputsPreprocessingPheno(
            beta = beta,
            pheno = pheno,
            logs = FALSE
        ),
        "duplicate CpG identifiers"
    )
})

test_that("timepoint parsing ignores missing rows and removes duplicate requests", {
    beta <- matrix(
        c(0.1, 0.2, 0.3),
        nrow = 1,
        dimnames = list("cg1", c("S1", "S2", "S3"))
    )
    metrics <- list(beta = beta, m = beta, cn = beta)
    pheno <- data.frame(
        Sample_Name = c("S1", "S2", "S3"),
        Timepoint = c("1", NA_character_, "2"),
        stringsAsFactors = FALSE
    )

    result <- splitTimepointsPreprocessingPheno(
        pheno = pheno,
        metricsData = metrics,
        timepoints = "1,1,2",
        logs = FALSE
    )

    expect_identical(result$timepoints, c("1", "2"))
    expect_identical(result$data[["1"]]$pheno$Sample_Name, "S1")
    expect_identical(result$data[["2"]]$pheno$Sample_Name, "S3")
})
