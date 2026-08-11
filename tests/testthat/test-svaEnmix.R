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
            anovaFull = suppressWarnings(lapply(full_models, stats::anova)),
            anovaReduced = suppressWarnings(lapply(full_models, stats::anova))
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

test_that("mergeSvaTargetsEnmix prevents repeated surrogate-variable columns", {
    targets <- data.frame(
        Sample_Name = c("S1", "S2", "S3"),
        PC1 = c(0, 0, 0),
        stringsAsFactors = FALSE
    )
    expect_error(
        mergeSvaTargetsEnmix(
            targets = targets,
            sva = create_toy_sva_analysis()$sva,
            SampleID = "Sample_Name"
        ),
        "would be duplicated"
    )
})

test_that("ENmix control-SVA options follow the documented flag contract", {
    expect_error(
        estimateSvaEnmixControls(NULL, ctrlSvaPercVar = -0.01),
        "between 0 and 1"
    )
    expect_error(
        estimateSvaEnmixControls(NULL, ctrlSvaFlag = 0),
        "must be 1"
    )
    expect_error(
        estimateSvaEnmixControls(NULL, ctrlSvaFlag = 3),
        "must be 1"
    )
})

test_that("Sentrix analysis omits technical factors without variation", {
    testthat::skip_if_not_installed("SummarizedExperiment")
    testthat::skip_if_not_installed("S4Vectors")

    sample_ids <- paste0("S", seq_len(6))
    sva <- matrix(
        seq_len(12) / 10,
        nrow = 6,
        dimnames = list(sample_ids, c("PC1", "PC2"))
    )
    rgset <- SummarizedExperiment::SummarizedExperiment(
        assays = list(signal = matrix(1, nrow = 2, ncol = 6)),
        colData = S4Vectors::DataFrame(
            Sentrix_ID = rep("Chip1", 6),
            Sentrix_Position = rep("R01C01", 6),
            row.names = sample_ids
        )
    )

    result <- analyzeSvaEnmix(
        sva = sva,
        RGSet = rgset,
        logs = FALSE
    )

    expect_length(result$technicalTerms, 0L)
    expect_true(all(vapply(result$fullModels, function(model) {
        identical(attr(stats::terms(model), "term.labels"), character(0))
    }, logical(1))))
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

test_that("plotSvaEnmix saves a technical-factor association summary", {
    tmp <- withr::local_tempdir()
    plot_file <- file.path(tmp, "sva_technicalFactorAssociations.tiff")
    analysis_data <- create_large_sva_analysis(
        sample_count = 48L, K = 3L, sentrix_id_count = 4L
    )
    technical_data <- data.frame(
        sentrix_id = analysis_data$sentrixID,
        sentrix_position = analysis_data$sentrixPosition
    )
    models <- lapply(seq_len(ncol(analysis_data$sva)), function(index) {
        stats::lm(
            analysis_data$sva[, index] ~ sentrix_id + sentrix_position,
            data = technical_data
        )
    })
    analysis_data$anovaFull <- suppressWarnings(lapply(models, stats::anova))

    output <- plotSvaEnmix(
        analysisData = analysis_data,
        plot = "association", file = plot_file,
        width = 700, height = 500, res = 72
    )

    expect_identical(output, plot_file)
    expect_true(file.exists(plot_file))
    association_plot <- dnaEPICO:::svaEnmixPlotAssociations(
        analysis_data
    )
    expect_null(association_plot$labels$title)
    expect_null(association_plot$labels$caption)
})

test_that("plotSvaEnmix omits empty technical-factor association figures", {
    tmp <- withr::local_tempdir()
    plot_file <- file.path(tmp, "stale_association.tiff")
    writeLines("stale", plot_file)
    analysis_data <- create_large_sva_analysis()
    analysis_data$anovaFull <- list()

    output <- plotSvaEnmix(
        analysisData = analysis_data, plot = "association",
        file = plot_file, width = 700, height = 500, res = 72
    )

    expect_null(output)
    expect_false(file.exists(plot_file))
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
        CD8T = c(0.20, 0.25, 0.30),
        PC1 = analysis_data$sva[, 1],
        PC2 = analysis_data$sva[, 2],
        stringsAsFactors = FALSE
    )
    pheno_file <- file.path(tmp, "pheno.csv")
    utils::write.csv(
        merged_pheno[c("Sample_Name", "CD8T")],
        pheno_file,
        row.names = FALSE
    )

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
    expect_identical(
        names(utils::read.csv(pheno_file, check.names = FALSE)),
        c("Sample_Name", "CD8T", "PC1", "PC2")
    )
})

test_that("writeSvaEnmixOutputs transactionally updates one phenotype file", {
    tmp <- withr::local_tempdir()
    analysis_data <- create_toy_sva_analysis()
    pheno_file <- file.path(tmp, "phenoLC.csv")
    input_data <- data.frame(
        Sample_Name = rownames(analysis_data$sva),
        CD8T = c(0.20, 0.25, 0.30),
        status = c("case", "control", "case"),
        stringsAsFactors = FALSE
    )
    utils::write.csv(input_data, pheno_file, row.names = FALSE)
    merged <- cbind(input_data, as.data.frame(analysis_data$sva))

    paths <- writeSvaEnmixOutputs(
        svaData = list(sva = analysis_data$sva),
        mergedPheno = merged,
        phenoFile = pheno_file,
        SampleID = "Sample_Name",
        dataBaseDir = file.path(tmp, "data"),
        rBaseDir = file.path(tmp, "rData"),
        scriptLabel = "svaEnmix"
    )

    updated <- utils::read.csv(pheno_file, check.names = FALSE)
    expect_identical(paths$phenoWithSva, pheno_file)
    expect_identical(updated$Sample_Name, input_data$Sample_Name)
    expect_equal(updated[c("CD8T", "status")], input_data[c("CD8T", "status")])
    expect_equal(
        unname(as.matrix(updated[c("PC1", "PC2")])),
        unname(as.matrix(merged[c("PC1", "PC2")]))
    )
    expect_identical(
        list.files(tmp, pattern = "^\\.phenoLC\\.csv\\.(sva|backup)-", all.files = TRUE),
        character(0)
    )
    expect_identical(
        list.files(tmp, pattern = "^phenoLC", all.files = TRUE),
        "phenoLC.csv"
    )
})

test_that("failed phenotype validation leaves the original file unchanged", {
    tmp <- withr::local_tempdir()
    pheno_file <- file.path(tmp, "phenoLC.csv")
    original <- data.frame(
        Sample_Name = c("S1", "S2", "S3"),
        CD8T = c(0.20, 0.25, 0.30),
        stringsAsFactors = FALSE
    )
    utils::write.csv(original, pheno_file, row.names = FALSE)
    checksum <- tools::md5sum(pheno_file)
    invalid <- cbind(original, PC1 = c(0.1, Inf, 0.3))

    expect_error(
        getFromNamespace("replacePhenotypeFileSvaEnmix", "dnaEPICO")(
            mergedPheno = invalid,
            phenoFile = pheno_file,
            SampleID = "Sample_Name",
            pcColumns = "PC1"
        ),
        "only finite values"
    )

    expect_identical(unname(tools::md5sum(pheno_file)), unname(checksum))
    expect_equal(utils::read.csv(pheno_file), original)
    expect_identical(
        list.files(tmp, pattern = "^\\.phenoLC\\.csv\\.(sva|backup)-", all.files = TRUE),
        character(0)
    )
})

test_that("svaEnmix exposes one phenotype-file interface", {
    expect_false("outputPhenoFile" %in% names(formals(svaEnmix)))
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
    expect_equal(colnames(result$svaData$sva), paste0("PC", seq_len(result$svaData$K)))
    expect_false(anyDuplicated(colnames(result$mergedPheno)) > 0L)
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
