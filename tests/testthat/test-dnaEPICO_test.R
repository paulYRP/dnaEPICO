test_that("preprocessingMinfiEwasWater runs on GSE285813 example data", {

  # Locate extdata
  baseDir <- system.file("extdata", "GSE285813", package = "dnaEPICO")
  expect_true(dir.exists(baseDir))

  phenoFile <- file.path(baseDir, "data/preprocessingMinfiEwasWater/pheno.csv")
  idatFolder <- file.path(baseDir, "data/preprocessingMinfiEwasWater/idats")
  crossReactivePath <-
    file.path(baseDir,
              "data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv")
  lcPhenoDir <- file.path(baseDir, "data/preprocessingMinfiEwasWater")

  expect_true(file.exists(phenoFile))
  expect_true(dir.exists(idatFolder))
  expect_true(file.exists(crossReactivePath))
  expect_true(dir.exists(lcPhenoDir))

  tmpDir <- tempdir()
  oldWd <- getwd()
  setwd(tmpDir)
  on.exit(setwd(oldWd), add = TRUE)

  # Testing preprocessingMinfiEwasWater function
  # nSamples = 5, lower than this will not work due to visualisationn
  expect_error(
    preprocessingMinfiEwasWater(
      phenoFile = phenoFile,
      idatFolder = idatFolder,
      crossReactivePath = crossReactivePath,
      lcPhenoDir = lcPhenoDir,
      nSamples = 5,
      SampleID = "Geo_Accession",
      phenoOrder = paste(
        "Geo_Accession",
        "Timepoint",
        "Sex",
        "PredSex",
        "Basename",
        "Sentrix_ID",
        "Sentrix_Position",
        sep = ";"
      )
    ),
    NA
  )
})
