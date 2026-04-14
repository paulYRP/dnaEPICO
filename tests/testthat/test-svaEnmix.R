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
