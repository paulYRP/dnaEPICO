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
