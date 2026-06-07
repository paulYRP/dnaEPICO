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

create_large_sva_analysis <- function(
    sample_count = 48L,
    K = 6L,
    sentrix_id_count = 44L
) {
    sample_ids <- paste0("S", seq_len(sample_count))
    sentrix_id_levels <- paste0("Chip", seq_len(sentrix_id_count))
    sentrix_position_levels <- paste0("R", sprintf("%02d", seq_len(8L)), "C01")
    sva <- matrix(
        seq_len(sample_count * K) / 100,
        nrow = sample_count,
        dimnames = list(sample_ids, paste0("PC", seq_len(K)))
    )

    structure(
        list(
            sva = sva,
            K = K,
            sentrixID = factor(
                sentrix_id_levels[
                    ((seq_len(sample_count) - 1L) %% sentrix_id_count) + 1L
                ],
                levels = sentrix_id_levels
            ),
            sentrixPosition = factor(
                sentrix_position_levels[
                    ((seq_len(sample_count) - 1L) %%
                        length(sentrix_position_levels)) + 1L
                ],
                levels = sentrix_position_levels
            )
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

test_that("plotSvaEnmix adapts matrix plots and SentrixID legends", {
    tmp <- withr::local_tempdir()
    log_dir <- file.path(tmp, "logs")
    analysis_data <- create_large_sva_analysis()
    plot_file <- file.path(tmp, "sva_SentrixIDPosition.tiff")

    output <- plotSvaEnmix(
        analysisData = analysis_data,
        plot = "matrix",
        file = plot_file,
        width = 600,
        height = 400,
        res = 150,
        logs = TRUE,
        log_dir = log_dir
    )

    expect_identical(output, plot_file)
    expect_true(file.exists(plot_file))

    log_lines <- readLines(file.path(log_dir, "log_plotSvaEnmix.txt"))
    expect_true(any(grepl(
        "SVA matrix plot dimensions:",
        log_lines,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "SentrixID legend suppressed for matrix plot",
        log_lines,
        fixed = TRUE
    )))
})

test_that("plotSvaEnmix paginates large matrix plots", {
    tmp <- withr::local_tempdir()
    log_dir <- file.path(tmp, "logs")
    analysis_data <- create_large_sva_analysis(sample_count = 56L, K = 7L)
    plot_file <- file.path(tmp, "sva_SentrixIDPosition.tiff")

    output <- plotSvaEnmix(
        analysisData = analysis_data,
        plot = "matrix",
        file = plot_file,
        width = 600,
        height = 400,
        res = 150,
        logs = TRUE,
        log_dir = log_dir
    )

    expect_length(output, 4L)
    expect_true(all(file.exists(output)))
    expect_false(file.exists(plot_file))
    expect_true(all(grepl("_page0[1-4]\\.tiff$", output)))

    log_lines <- readLines(file.path(log_dir, "log_plotSvaEnmix.txt"))
    expect_true(any(grepl(
        "SVA matrix plot pages: 4 pages",
        log_lines,
        fixed = TRUE
    )))
    expect_true(any(grepl("Saved plot paths:", log_lines, fixed = TRUE)))
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

test_that("svaEnmix accepts wrapper files that contain an RGSet element", {
    testthat::skip_if_not_installed("minfiData")

    tmp <- withr::local_tempdir()
    ex <- dnaEPICO:::exampleMinfiBaseDataDnaEpico()
    pheno_file <- file.path(tmp, "pheno.csv")
    rgset_path <- file.path(tmp, "RGSet_wrapper.RData")
    rgset_state <- list(
        RGSet = ex$RGSet,
        note = "wrapper object"
    )
    utils::write.csv(ex$targets, pheno_file, row.names = FALSE)
    save(rgset_state, file = rgset_path)

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

test_that("svaEnmix reports invalid saved RGSet objects clearly", {
    tmp <- withr::local_tempdir()
    pheno_file <- file.path(tmp, "pheno.csv")
    bad_path <- file.path(tmp, "bad_RGSet.RData")
    bad_object <- list(not_rgset = 1)
    pheno <- data.frame(
        UID = "S1",
        Sentrix_ID = "Chip1",
        Sentrix_Position = "R01C01",
        stringsAsFactors = FALSE
    )
    utils::write.csv(pheno, pheno_file, row.names = FALSE)
    save(bad_object, file = bad_path)

    expect_error(
        svaEnmix(
            phenoFile = pheno_file,
            rgsetData = bad_path,
            SampleID = "UID",
            SentrixIDColumn = "Sentrix_ID",
            SentrixPositionColumn = "Sentrix_Position",
            outputLogs = file.path(tmp, "logs"),
            logs = TRUE,
            verbose = FALSE,
            saveOutputs = FALSE
        ),
        "does not expose a usable sample count"
    )
})
