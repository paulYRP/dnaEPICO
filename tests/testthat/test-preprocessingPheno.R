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
    expect_equal(nrow(result$combinedData$phenoBeta), 3)
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
    expect_true(file.exists(file.path(tmp, "rData", "preprocessingPheno", "mergeData", "phenoBetaT1.RData")))
    expect_true(file.exists(file.path(tmp, "clockFoundation", "beta.csv")))
    expect_true(file.exists(file.path(tmp, "clockFoundation", "phenoCF.csv")))
})
