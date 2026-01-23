test_that("preprocessingMinfiEwasWater runs using staged minfiData IDATs", {

  # Skip test if example data are not available
  skip_if_not_installed("minfiData")

  # Load required package
  library(minfiData)

  # ------------------------------------------------------------------
  # 1. Locate minfiData example files
  # ------------------------------------------------------------------

  baseDirMinfi <- system.file("extdata", package = "minfiData")

  idatFiles <- list.files(
    baseDirMinfi,
    pattern = "\\.idat$",
    recursive = TRUE,
    full.names = TRUE
  )

  expect_true(length(idatFiles) > 0)

  # ------------------------------------------------------------------
  # 2. Locate dnaEPICO packaged extdata (cross-reactive probes)
  # ------------------------------------------------------------------

  baseDirDnaEPICO <- system.file("extdata", package = "dnaEPICO")
  expect_true(dir.exists(baseDirDnaEPICO))

  crossReactivePath <- file.path(
    baseDirDnaEPICO,
    "12864_2024_10027_MOESM8_ESM.csv"
  )

  expect_true(file.exists(crossReactivePath))

  # ------------------------------------------------------------------
  # 3. Create temporary working directory
  # ------------------------------------------------------------------

  tmpDir <- withr::local_tempdir()

  tmpIdatDir <- file.path(tmpDir, "idats")
  dir.create(tmpIdatDir)

  file.copy(idatFiles, tmpIdatDir, overwrite = TRUE)

  # ------------------------------------------------------------------
  # 4. Create phenotype file (pheno_minfiData.csv)
  # ------------------------------------------------------------------

  targets <- read.csv(
    file.path(baseDirMinfi, "SampleSheet.csv"),
    stringsAsFactors = FALSE,
    skip = 7
  )

  targets$Basename <- paste0(
    targets$Sentrix_ID,
    "_",
    targets$Sentrix_Position
  )

  phenoFile <- file.path(tmpDir, "pheno_minfiData.csv")

  write.csv(
    targets,
    phenoFile,
    row.names = FALSE
  )

  expect_true(file.exists(phenoFile))

  # ------------------------------------------------------------------
  # 5. Set working directory safely
  # ------------------------------------------------------------------

  oldWd <- getwd()
  setwd(tmpDir)
  on.exit(setwd(oldWd), add = TRUE)

  # ------------------------------------------------------------------
  # 6. Run preprocessing pipeline
  # ------------------------------------------------------------------

  expect_error(
    preprocessingMinfiEwasWater(
      phenoFile = phenoFile,
      idatFolder = tmpIdatDir,
      nSamples = 5,
      SampleID = "Sample_Name",
      phenoOrder = "Sample_Name;sex;Basename;Sentrix_ID;Sentrix_Position",
      arrayType = "IlluminaHumanMethylation450k",
      annotationVersion = "ilmn12.hg19",
      sexColumn = "sex",
      plotGroupVar = "sex",
      outputLogs = file.path(tmpDir, "logs"),
      baseDataFolder = file.path(tmpDir, "rData"),
      crossReactivePath = crossReactivePath,
      lcPhenoDir = tmpDir
    ),
    NA
  )

  # ------------------------------------------------------------------
  # 7. Minimal output sanity checks
  # ------------------------------------------------------------------

  expect_true(dir.exists(file.path(tmpDir, "logs")))
  expect_true(dir.exists(file.path(tmpDir, "rData")))
})
