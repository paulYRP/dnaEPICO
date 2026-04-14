create_dnam_report_example <- function(path) {
  dirs <- list(
    qcDir = file.path(path, "figures", "preprocessingMinfiEwasWater", "enMix"),
    preprocessingDir = file.path(path, "figures", "preprocessingMinfiEwasWater", "qc"),
    postprocessingDir = file.path(path, "figures", "preprocessingMinfiEwasWater", "metrics"),
    svaDir = file.path(path, "figures", "svaEnmix"),
    glmDir = file.path(path, "figures", "methylationGLM_T1"),
    glmmDir = file.path(path, "figures", "methylationGLMM_T1T2"),
    figDir = file.path(path, "reports", "figures"),
    outputDir = file.path(path, "reports"),
    logDir = file.path(path, "logs")
  )

  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  glm_subdir <- file.path(dirs$glmDir, "status")
  glmm_subdir <- file.path(dirs$glmmDir, "score")
  dir.create(glm_subdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(glmm_subdir, recursive = TRUE, showWarnings = FALSE)

  file.create(file.path(dirs$qcDir, "qc_1.jpg"))
  file.create(file.path(dirs$preprocessingDir, "pre_1.tiff"))
  file.create(file.path(dirs$postprocessingDir, "post_1.tiff"))
  file.create(file.path(dirs$svaDir, "sva_1.tiff"))
  file.create(file.path(glm_subdir, "qqplot_status.tiff"))
  file.create(file.path(glmm_subdir, "qqplot_score.tiff"))

  dirs
}

test_that("prepareDnamReportInputs returns a structured inventory quietly", {
  tmp <- withr::local_tempdir()
  example_dirs <- create_dnam_report_example(tmp)

  expect_message(
    result <- prepareDnamReportInputs(
      output = "DNAm_Report.pdf",
      outputDir = example_dirs$outputDir,
      qcDir = example_dirs$qcDir,
      preprocessingDir = example_dirs$preprocessingDir,
      postprocessingDir = example_dirs$postprocessingDir,
      svaDir = example_dirs$svaDir,
      glmDir = example_dirs$glmDir,
      glmmDir = example_dirs$glmmDir,
      figDir = example_dirs$figDir,
      reportTitle = "DNA methylation",
      author = "School of Biomedical Sciences",
      date = format(Sys.Date(), "%B %d, %Y"),
      verbose = FALSE,
      logs = FALSE
    ),
    NA
  )

  expect_s3_class(result, "dnaEPICO_dnamReport_prepared")
  expect_match(basename(result$templatePath), "dnamReport\\.Rmd$")
  expect_equal(result$figureInventory$qc$count, 1)
  expect_equal(result$figureInventory$glm$count, 1)
  expect_equal(result$figureInventory$glmm$count, 1)
  expect_length(result$missingFigureDirectories, 0)
  expect_match(result$outputFile, "DNAm_Report\\.pdf$")
})

test_that("dnamReport can prepare without rendering and write logs on request", {
  tmp <- withr::local_tempdir()
  example_dirs <- create_dnam_report_example(tmp)

  expect_message(
    result <- dnamReport(
      output = "DNAm_Report.pdf",
      outputDir = example_dirs$outputDir,
      qcDir = example_dirs$qcDir,
      preprocessingDir = example_dirs$preprocessingDir,
      postprocessingDir = example_dirs$postprocessingDir,
      svaDir = example_dirs$svaDir,
      glmDir = example_dirs$glmDir,
      glmmDir = example_dirs$glmmDir,
      figDir = example_dirs$figDir,
      reportTitle = "DNA methylation",
      author = "School of Biomedical Sciences",
      date = format(Sys.Date(), "%B %d, %Y"),
      render = FALSE,
      verbose = TRUE,
      logs = TRUE,
      logDir = example_dirs$logDir
    ),
    "Starting DNA Methylation Report Step"
  )

  expect_s3_class(result, "dnaEPICO_dnamReport")
  expect_identical(result$status, "prepared")
  expect_s3_class(result$preparedReport, "dnaEPICO_dnamReport_prepared")
  expect_true(file.exists(result$logFile))
  expect_false(file.exists(result$outputFile))
  expect_null(result$renderResult)
})
