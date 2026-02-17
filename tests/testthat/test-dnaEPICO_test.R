test_that("preprocessingMinfiEwasWater runs using staged minfiData IDATs", {
  
  library(minfi)
  library(minfiData)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

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

  targets$Timepoint <- 1

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
    dnaEPICO::preprocessingMinfiEwasWater(
      phenoFile = phenoFile,
      idatFolder = tmpIdatDir,
      nSamples = 6,
      SampleID = "Sample_Name",
      phenoOrder = "Sample_Name;sex;Basename;Sentrix_ID;Sentrix_Position",
      arrayType = "IlluminaHumanMethylation450k",
      annotationVersion = "ilmn12.hg19",
      sexColumn = "sex",
      plotGroupVar = "sex",
      lcRef = "saliva",
      outputLogs = file.path(tmpDir, "logs"),
      baseDataFolder = file.path(tmpDir, "rData"),
      crossReactivePath = crossReactivePath,
      lcPhenoDir = tmpDir
    ),
    NA
  )

  # ------------------------------------------------------------------
  # 7. Run sva
  # ------------------------------------------------------------------

  rgsetPath <- file.path(
    tmpDir,
    "rData",
    "preprocessingMinfiEwasWater",
    "objects",
    "RGSet.RData"
  )

  expect_true(file.exists(rgsetPath))

  expect_error(
      dnaEPICO::svaEnmix(
        phenoFile = phenoFile,
        rgsetData = rgsetPath,
        sepType = "",
        outputLogs = file.path(tmpDir, "logs"),
        nSamples = NA,
        SampleID = "Sample_Name",
        arrayType = "IlluminaHumanMethylation450k",
        annotationVersion = "ilmn12.hg19",
        SentrixIDColumn = "Sentrix_ID",
        SentrixPositionColumn = "Sentrix_Position",
        ctrlSvaPercVar = 0.90,
        ctrlSvaFlag = 1,
        scriptLabel = "svaEnmix",
        tiffWidth = 2000,
        tiffHeight = 1000,
        tiffRes = 150,
        dataBaseDir = file.path(tmpDir, "data") ,
        figureBaseDir = file.path(tmpDir, "figures")
        
      ),
      NA
    )

  # ------------------------------------------------------------------
  # 8. Run preprocessing phenotype files
  # ------------------------------------------------------------------

    betaPath <- file.path(
      tmpDir,
      "rData",
      "preprocessingMinfiEwasWater",
      "metrics",
      "beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData"
    )

    mPath <- file.path(
      tmpDir,
      "rData",
      "preprocessingMinfiEwasWater",
      "metrics",
      "m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData"
    )

    cnPath <- file.path(
      tmpDir,
      "rData",
      "preprocessingMinfiEwasWater",
      "metrics",
      "cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData"
    )

    expect_true(file.exists(betaPath))
    expect_true(file.exists(mPath))
    expect_true(file.exists(cnPath))

    expect_error(
        dnaEPICO::preprocessingPheno(
          phenoFile = phenoFile,
          sepType = "",
          betaPath = betaPath,
          mPath = mPath,
          cnPath = cnPath,
          SampleID = "Sample_Name",
          timeVar = "Timepoint",
          timepoints = "1",
          combineTimepoints = "1",
          sexColumn = "sex",
          outputPheno = file.path(tmpDir, "data", "preprocessingPheno"),
          outputRData = file.path(tmpDir, "rData", "preprocessingPheno", "metrics"),
          outputRDataMerge = file.path(tmpDir, "rData", "preprocessingPheno", "mergeData"),
          outputLogs = file.path(tmpDir, "logs"),
          outputDir = file.path(tmpDir, "data", "preprocessingPheno")
        ),
        NA
      )

  # ------------------------------------------------------------------
  # 8. Generate Report
  # ------------------------------------------------------------------
    dnamReport(
      output = "DNAm_Report.pdf",
      outputDir = file.path(tmpDir, "reports"),

      qcDir = file.path(
        tmpDir, "figures", "preprocessingMinfiEwasWater", "enMix"
      ),
      preprocessingDir = file.path(
        tmpDir, "figures", "preprocessingMinfiEwasWater", "qc"
      ),
      postprocessingDir = file.path(
        tmpDir, "figures", "preprocessingMinfiEwasWater", "metrics"
      ),
      svaDir = file.path(
        tmpDir, "figures", "svaEnmix"
      ),
      glmDir = file.path(
        tmpDir, "figures", "methylationGLM_T1"
      ),
      glmmDir = file.path(
        tmpDir, "figures", "methylationGLMM_T1T2"
      ),

      reportTitle = "DNA methylation",
      author = "School of Biomedical Sciences",
      date = format(Sys.Date(), "%B %d, %Y")
    )
})
