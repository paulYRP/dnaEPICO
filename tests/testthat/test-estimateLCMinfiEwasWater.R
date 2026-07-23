build_estimate_lc_example_data <- function() {
    ref_file <- system.file("extdata", "saliva.txt", package = "dnaEPICO")
    beta <- as.matrix(utils::read.table(ref_file))[1:20, , drop = FALSE]
    colnames(beta) <- c("sample1", "sample2")

    targets <- data.frame(
        Sample_Name = colnames(beta),
        Timepoint = c("T1", "T2"),
        stringsAsFactors = FALSE
    )

    list(beta = beta, targets = targets)
}

test_that("estimateLC returns the expected cell-type columns", {
    example_data <- build_estimate_lc_example_data()

    result <- estimateLC(
        meth = example_data$beta,
        ref = "saliva",
        constrained = FALSE
    )

    expect_s3_class(result, "data.table")
    expect_identical(
        colnames(result),
        c("Leukocytes", "Epithelial.cells")
    )
    expect_equal(as.numeric(result[1, ]), c(1, 0), tolerance = 1e-8)
    expect_equal(as.numeric(result[2, ]), c(0, 1), tolerance = 1e-8)
})

test_that("estimateLC averages duplicate EPICv2 loci instead of choosing one", {
    example_data <- build_estimate_lc_example_data()
    epic <- example_data$beta
    rownames(epic) <- paste0(rownames(epic), "_TC11")
    duplicate <- epic[1, , drop = FALSE]
    rownames(duplicate) <- sub("_TC11$", "_BC11", rownames(duplicate))
    duplicate[,] <- c(0.2, 0.8)
    epic_with_duplicate <- rbind(epic, duplicate)

    collapsed <- epic
    collapsed[1, ] <- colMeans(epic_with_duplicate[c(1, nrow(epic_with_duplicate)), , drop = FALSE])

    expect_equal(
        as.matrix(estimateLC(epic_with_duplicate, ref = "saliva", constrained = FALSE)),
        as.matrix(estimateLC(collapsed, ref = "saliva", constrained = FALSE)),
        tolerance = 1e-10
    )
})

test_that("estimateLC rejects duplicate exact probe identifiers", {
    example_data <- build_estimate_lc_example_data()
    duplicated_beta <- rbind(example_data$beta, example_data$beta[1, , drop = FALSE])

    expect_error(
        estimateLC(duplicated_beta, ref = "saliva", constrained = FALSE),
        "duplicate CpG identifiers"
    )
})

test_that("estimateLCMinfiEwasWater merges and orders phenoLC columns", {
    example_data <- build_estimate_lc_example_data()

    result <- estimateLCMinfiEwasWater(
        beta = example_data$beta,
        targets = example_data$targets,
        lcRef = "saliva",
        phenoOrder = "Sample_Name;Timepoint"
    )

    expect_s3_class(result, "dnaEPICO_minfiEwasWater_lc")
    expect_identical(
        colnames(result$phenoLC)[1:4],
        c("Sample_Name", "Timepoint", "Leukocytes", "Epithelial.cells")
    )
    expect_equal(result$phenoLC$Leukocytes, c(1, 0), tolerance = 1e-8)
    expect_equal(result$phenoLC$Epithelial.cells, c(0, 1), tolerance = 1e-8)
})

test_that("estimateLCMinfiEwasWater validates reference and beta inputs", {
    example_data <- build_estimate_lc_example_data()
    targets <- data.frame(
        Sample_Name = colnames(example_data$beta),
        stringsAsFactors = FALSE
    )

    expect_error(
        estimateLCMinfiEwasWater(
            beta = example_data$beta,
            targets = targets,
            lcRef = c("saliva", "salivaEPIC")
        ),
        "one non-empty"
    )
    expect_error(
        estimateLCMinfiEwasWater(
            beta = example_data$beta,
            targets = targets,
            lcRef = "saliva",
            constrained = "FALSE"
        ),
        "TRUE or FALSE"
    )

    result <- estimateLCMinfiEwasWater(
        beta = example_data$beta,
        targets = targets,
        lcRef = "saliva"
    )
    expect_identical(
        names(result$methylationRange),
        c("Scale", "Observed.Minimum", "Observed.Maximum")
    )
    expect_equal(
        result$methylationRange$Observed.Minimum,
        min(example_data$beta)
    )
    expect_equal(
        result$methylationRange$Observed.Maximum,
        max(example_data$beta)
    )
    expect_false(any(c(
        "methylationBoundaries", "methylationIssues", "invalidCpGs"
    ) %in% names(result)))
})
