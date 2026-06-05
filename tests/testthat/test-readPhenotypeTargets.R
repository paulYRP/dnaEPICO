test_that("readPhenotypeTargets is quiet by default and subsets rows", {
    pheno <- data.frame(
        Sample_Name = c("S1", "S2", "S3"),
        Sex = c("F", "M", "F"),
        stringsAsFactors = FALSE
    )

    pheno_file <- tempfile(fileext = ".csv")
    utils::write.csv(pheno, pheno_file, row.names = FALSE)

    targets <- NULL

    expect_silent({
        targets <- dnaEPICO::readPhenotypeTargets(
            phenoFile = pheno_file,
            SampleID = "Sample_Name",
            nSamples = 2
        )
    })

    expect_s3_class(targets, "data.frame")
    expect_equal(nrow(targets), 2L)
    expect_identical(targets$Sample_Name, c("S1", "S2"))
})

test_that("readPhenotypeTargets treats absent separator values as comma defaults", {
    pheno <- data.frame(
        Sample_Name = c("S1", "S2"),
        Sex = c("F", "M"),
        stringsAsFactors = FALSE
    )

    pheno_file <- tempfile(fileext = ".csv")
    utils::write.csv(pheno, pheno_file, row.names = FALSE)

    for (sep_type in list(NULL, "", "NULL")) {
        targets <- dnaEPICO::readPhenotypeTargets(
            phenoFile = pheno_file,
            sepType = sep_type,
            SampleID = "Sample_Name"
        )

        expect_identical(targets$Sample_Name, pheno$Sample_Name)
        expect_identical(targets$Sex, pheno$Sex)
    }
})

test_that("readPhenotypeTargets preserves explicit tab separators", {
    pheno <- data.frame(
        Sample_Name = c("S1", "S2"),
        Sex = c("F", "M"),
        stringsAsFactors = FALSE
    )

    pheno_file <- tempfile(fileext = ".tsv")
    utils::write.table(
        pheno,
        pheno_file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
    )

    for (sep_type in list("\t", "\\t")) {
        targets <- dnaEPICO::readPhenotypeTargets(
            phenoFile = pheno_file,
            sepType = sep_type,
            SampleID = "Sample_Name"
        )

        expect_identical(targets$Sample_Name, pheno$Sample_Name)
        expect_identical(targets$Sex, pheno$Sex)
    }
})

test_that("readPhenotypeTargets can emit verbose messages and write logs", {
    pheno <- data.frame(
        Sample_Name = c("S1", "S2"),
        Sex = c("F", "M"),
        stringsAsFactors = FALSE
    )

    pheno_file <- tempfile(fileext = ".csv")
    log_dir <- tempfile("readTargets-log-")
    utils::write.csv(pheno, pheno_file, row.names = FALSE)

    expect_message(
        dnaEPICO::readPhenotypeTargets(
            phenoFile = pheno_file,
            SampleID = "Sample_Name",
            verbose = TRUE,
            logs = TRUE,
            log_dir = log_dir
        ),
        "Phenotype file loaded with 2 samples and 2 columns."
    )

    log_file <- file.path(log_dir, "log_readPhenotypeTargets.txt")
    expect_true(file.exists(log_file))

    log_lines <- readLines(log_file, warn = FALSE)
    expect_true(any(grepl("Preview of targets:", log_lines, fixed = TRUE)))
})

test_that("readPhenotypeTargets validates the SampleID column", {
    pheno <- data.frame(
        sample_name = c("S1", "S2"),
        Sex = c("F", "M"),
        stringsAsFactors = FALSE
    )

    pheno_file <- tempfile(fileext = ".csv")
    utils::write.csv(pheno, pheno_file, row.names = FALSE)

    expect_error(
        dnaEPICO::readPhenotypeTargets(
            phenoFile = pheno_file,
            SampleID = "Sample_Name"
        ),
        "SampleID column not found in phenotype data"
    )
})
