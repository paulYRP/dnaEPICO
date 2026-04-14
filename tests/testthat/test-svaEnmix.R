create_toy_sva_analysis <- function() {
    sva <- matrix(
        c(-1.0, 1.0, 0.0, 0.5, -0.5, 0.0),
        nrow = 3,
        dimnames = list(c("S1", "S2", "S3"), c("PC1", "PC2"))
    )
    sentrix_id <- factor(c("Chip1", "Chip1", "Chip2"))
    sentrix_position <- factor(c("R02C02", "R04C01", "R04C02"))

    full_models <- lapply(
        seq_len(ncol(sva)),
        function(i) {
            stats::lm(
                sva[, i] ~ sentrix_id + sentrix_position,
                data = data.frame(
                    sentrix_id = sentrix_id,
                    sentrix_position = sentrix_position
                )
            )
        }
    )

    structure(
        list(
            sva = sva,
            K = ncol(sva),
            sentrixID = sentrix_id,
            sentrixPosition = sentrix_position,
            fullModels = full_models,
            reducedModels = full_models,
            droptermSteps = list(list(), list()),
            anovaFull = lapply(full_models, stats::anova),
            anovaReduced = lapply(full_models, stats::anova)
        ),
        class = "dnaEPICO_svaEnmix_analysis"
    )
}

test_that("mergeSvaTargetsEnmix preserves target order", {
    targets <- data.frame(
        Sample_Name = c("S2", "S1", "S3"),
        Group = c("A", "A", "B"),
        stringsAsFactors = FALSE
    )
    sva <- create_toy_sva_analysis()$sva

    merged <- mergeSvaTargetsEnmix(
        targets = targets,
        sva = sva,
        SampleID = "Sample_Name"
    )

    expect_identical(merged$Sample_Name, targets$Sample_Name)
    expect_true(all(c("PC1", "PC2") %in% colnames(merged)))
})

test_that("plotSvaEnmix saves a sentrix-id plot", {
    tmp <- withr::local_tempdir()
    analysis_data <- create_toy_sva_analysis()
    plot_file <- file.path(tmp, "sva_SentrixID.tiff")

    output <- plotSvaEnmix(
        analysisData = analysis_data,
        plot = "sentrix_id",
        file = plot_file
    )

    expect_identical(output, plot_file)
    expect_true(file.exists(plot_file))
})

test_that("writeSvaEnmixOutputs writes matrix and summary files", {
    tmp <- withr::local_tempdir()
    analysis_data <- create_toy_sva_analysis()
    sva_data <- structure(
        list(
            sva = analysis_data$sva,
            K = analysis_data$K,
            ctrlSvaPercVar = 0.9,
            ctrlSvaFlag = 1
        ),
        class = "dnaEPICO_svaEnmix_sva"
    )
    merged_pheno <- data.frame(
        Sample_Name = c("S1", "S2", "S3"),
        PC1 = analysis_data$sva[, 1],
        PC2 = analysis_data$sva[, 2],
        stringsAsFactors = FALSE
    )
    pheno_file <- file.path(tmp, "pheno.csv")

    paths <- writeSvaEnmixOutputs(
        svaData = sva_data,
        mergedPheno = merged_pheno,
        analysisData = analysis_data,
        phenoFile = pheno_file,
        dataBaseDir = file.path(tmp, "data"),
        rBaseDir = file.path(tmp, "rData"),
        scriptLabel = "svaEnmix"
    )

    expect_s3_class(paths, "dnaEPICO_svaEnmix_paths")
    expect_true(file.exists(paths$svaRData))
    expect_true(file.exists(paths$svaCSV))
    expect_true(file.exists(paths$phenoWithSva))
    expect_true(file.exists(file.path(paths$dataDir, "summary_full_sva1.txt")))
    expect_true(file.exists(file.path(paths$dataDir, "anova_full_sva1.txt")))
})

test_that("svaEnmix logs fatal errors before stopping", {
    tmp <- withr::local_tempdir()
    pheno_file <- file.path(tmp, "pheno.csv")
    log_dir <- file.path(tmp, "logs")
    pheno <- data.frame(
        UID = "S1",
        Sentrix_ID = "Chip1",
        Sentrix_Position = "R01C01",
        stringsAsFactors = FALSE
    )
    utils::write.csv(pheno, pheno_file, row.names = FALSE)

    expect_error(
        svaEnmix(
            phenoFile = pheno_file,
            rgsetData = file.path(tmp, "missing_RGSet.RData"),
            SampleID = "UID",
            SentrixIDColumn = "Sentrix_ID",
            SentrixPositionColumn = "Sentrix_Position",
            outputLogs = log_dir,
            logs = TRUE,
            verbose = FALSE,
            saveOutputs = FALSE
        ),
        "Input file does not exist"
    )

    log_file <- file.path(log_dir, "log_svaEnmix.txt")
    expect_true(file.exists(log_file))
    log_lines <- readLines(log_file, warn = FALSE)
    expect_true(any(grepl("ERROR in svaEnmix", log_lines, fixed = TRUE)))
    expect_true(any(grepl("Input file does not exist", log_lines, fixed = TRUE)))
    expect_true(any(grepl("Failing call:", log_lines, fixed = TRUE)))
    expect_true(any(grepl("Call stack:", log_lines, fixed = TRUE)))
    expect_true(any(grepl("svaEnmix", log_lines, fixed = TRUE)))
})

test_that("svaEnmix accepts saved RGSet inputs", {
    testthat::skip_if_not_installed("minfiData")

    tmp <- withr::local_tempdir()
    ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
    pheno_file <- file.path(tmp, "pheno.csv")
    rgset_path <- file.path(tmp, "RGSet.RData")
    RGSet <- ex$RGSet
    utils::write.csv(ex$targets, pheno_file, row.names = FALSE)
    save(RGSet, file = rgset_path)

    result <- svaEnmix(
        phenoFile = pheno_file,
        rgsetData = rgset_path,
        SampleID = "Sample_Name",
        arrayType = "IlluminaHumanMethylation450k",
        annotationVersion = "ilmn12.hg19",
        SentrixIDColumn = "Sentrix_ID",
        SentrixPositionColumn = "Sentrix_Position",
        outputLogs = file.path(tmp, "logs"),
        figureBaseDir = file.path(tmp, "figures"),
        dataBaseDir = file.path(tmp, "data"),
        rBaseDir = file.path(tmp, "rData"),
        verbose = FALSE,
        logs = FALSE,
        saveOutputs = FALSE
    )

    expect_s3_class(result, "dnaEPICO_svaEnmix")
    expect_equal(nrow(result$targets), ncol(result$RGSet))
})
