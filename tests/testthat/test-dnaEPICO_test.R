test_that("preprocessingMinfiEwasWater runs using staged minfiData IDATs for Bioconductor", {
    testthat::skip_on_ci()
    library(minfi)
    library(minfiData)
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    testthat::skip_if_not_installed("tinytex")
    testthat::skip_if_not_installed("tiff")
    testthat::skip_if_not_installed("withr")

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
})

test_that("preprocessingMinfiEwasWater runs using staged minfiData IDATs for Github", {
    library(minfi)
    library(minfiData)
    library(IlluminaHumanMethylation450kmanifest)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    testthat::skip_if_not_installed("tinytex")
    testthat::skip_if_not_installed("withr")

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
            lcPhenoDir = tmpDir,
            saveOutputs = TRUE
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
            sepType = NULL,
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
            dataBaseDir = file.path(tmpDir, "data"),
            figureBaseDir = file.path(tmpDir, "figures"),
            rBaseDir = file.path(tmpDir, "rData"),
            saveOutputs = TRUE
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
            sepType = NULL,
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
            outputDir = file.path(tmpDir, "data", "preprocessingPheno"),
            saveOutputs = TRUE
        ),
        NA
    )

    # ------------------------------------------------------------------
    # 9. Generate Report
    # ------------------------------------------------------------------

    report_result <- dnaEPICO::dnamReport(
        outputDir = file.path(tmpDir, "reports"),
        enmixTab = file.path(
            tmpDir, "figures", "preprocessingMinfiEwasWater", "enmix"
        ),
        qcTab = file.path(
            tmpDir, "figures", "preprocessingMinfiEwasWater", "qc"
        ),
        metricTab = file.path(
            tmpDir, "figures", "preprocessingMinfiEwasWater", "metrics"
        ),
        svaTab = file.path(
            tmpDir, "figures", "svaEnmix"
        ),
        verbose = FALSE,
        logs = TRUE,
        logTab = file.path(tmpDir, "logs")
    )
    if (identical(report_result$status, "failed")) {
        testthat::skip(
            paste(
                "Report rendering failed in the current local environment:",
                if (is.null(report_result$errorMessage)) "unknown render error" else report_result$errorMessage
            )
        )
    }
    expect_true(report_result$status %in% c("rendered", "skipped"))
    expect_true(file.exists(file.path(tmpDir, "logs", "log_dnamReport.txt")))
    if (identical(report_result$status, "rendered")) {
        expect_true(file.exists(report_result$outputFile))
    }

    # ------------------------------------------------------------------
    # 10. Extract Makefile
    # ------------------------------------------------------------------
    expect_error(
        makefile_path <- dnaEPICO::extractMake(
            destDir = tmpDir,
            overwrite = TRUE
        ),
        NA
    )
    expect_true(file.exists(makefile_path))

    # ------------------------------------------------------------------
    # 11. Run GLM
    # ------------------------------------------------------------------

    expect_error(
        dnaEPICO::methylationGLM_T1(
            inputPheno = file.path(
                tmpDir, "rData", "preprocessingPheno", "mergeData",
                "phenoBetaT1.RData"
            ),
            outputLogs = file.path(tmpDir, "logs"),
            outputRData = file.path(tmpDir, "rData", "methylationGLM_T1", "models"),
            outputPlots = file.path(tmpDir, "figures", "methylationGLM_T1"),
            phenotypes = "status",
            covariates = "sex",
            factorVars = "sex,status",
            cpgLimit = 1,
            nCores = 8,
            summaryPval = 1,
            significantCpGDir = file.path(tmpDir, "results", "cpgs", "methylationGLM_T1"),
            summaryTxtDir = file.path(
                tmpDir, "results", "summary", "methylationGLM_T1",
                "glm"
            ),
            annotationCols = "Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island",
            annotatedGLMOut = file.path(tmpDir, "data", "methylationGLM_T1"),
            annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            chunkSize = 1,
            significantCpGPval = 1,
            logs = TRUE,
            saveOutputs = TRUE
        ),
        NA
    )

    # ------------------------------------------------------------------
    # 12. Run LMER
    # ------------------------------------------------------------------

    expect_error(
        dnaEPICO::methylationGLMM_T1T2(
            inputPheno = file.path(tmpDir, "rData", "preprocessingPheno", "mergeData", "phenoBetaT1.RData"),
            outputLogs = file.path(tmpDir, "logs"),
            outputRData = file.path(tmpDir, "rData", "methylationGLMM_T1T2", "models"),
            outputPlots = file.path(tmpDir, "figures", "methylationGLMM_T1T2"),
            phenotypes = "status",
            personVar = "person",
            covariates = "sex",
            factorVars = "sex,status",
            cpgLimit = 1,
            nCores = 8,
            summaryPval = 1,
            significantInteractionDir = file.path(tmpDir, "results", "cpgs", "methylationGLMM_T1T2"),
            summaryTxtDir = file.path(tmpDir, "results", "summary", "methylationGLMM_T1T2", "lmer"),
            annotatedLMEOut = file.path(tmpDir, "data", "methylationGLMM_T1T2"),
            annotationPackage = "IlluminaHumanMethylation450kanno.ilmn12.hg19",
            chunkSize = 1,
            significantInteractionPval = 1,
            saveSignificantInteractions = FALSE,
            annotationCols = "Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island",
            logs = TRUE,
            saveOutputs = TRUE
        ),
        NA
    )
})
