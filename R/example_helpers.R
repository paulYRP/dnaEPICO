#' Internal example helpers for documentation
#'
#' These helpers create small, reproducible objects used by runnable examples in
#' the package documentation.
#'
#' @keywords internal
#' @noRd

requireExamplePackageDnaEpico <- function(package_name) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
        stop(
            "Package '",
            package_name,
            "' is required for this example.",
            call. = FALSE
        )
    }
}

exampleCacheDnaEpico <- new.env(parent = emptyenv())

getCachedExampleDnaEpico <- function(key, builder) {
    if (exists(key, envir = exampleCacheDnaEpico, inherits = FALSE)) {
        return(get(key, envir = exampleCacheDnaEpico, inherits = FALSE))
    }

    value <- builder()
    assign(key, value, envir = exampleCacheDnaEpico)
    value
}

exampleTempDirDnaEpico <- function(prefix) {
    temp_dir <- tempfile(pattern = prefix)
    dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
    temp_dir
}

createDummyImageDnaEpico <- function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(filename = path, width = 200, height = 200)
    graphics::plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    graphics::text(1, 1, labels = basename(path))
    grDevices::dev.off()
    invisible(path)
}

exampleMinfiBaseDataDnaEpico <- function() {
    getCachedExampleDnaEpico(
        key = "minfi_base",
        builder = function() {
            requireExamplePackageDnaEpico("minfiData")

            example_env <- new.env(parent = emptyenv())
            utils::data("RGsetEx", package = "minfiData", envir = example_env)

            rgset <- get("RGsetEx", envir = example_env)
            targets <- as.data.frame(SummarizedExperiment::colData(rgset))

            if (!("Sample_Name" %in% colnames(targets))) {
                targets$Sample_Name <- Biobase::sampleNames(rgset)
            }

            if (!("Sex" %in% colnames(targets))) {
                if ("sex" %in% colnames(targets)) {
                    sex_values <- as.character(targets$sex)
                    targets$Sex <- ifelse(sex_values %in% c("0", "F", "Female"), "F", "M")
                } else {
                    targets$Sex <- rep(c("F", "M"), length.out = nrow(targets))
                }
            }

            if (!("Sentrix_ID" %in% colnames(targets))) {
                if ("Slide" %in% colnames(targets)) {
                    targets$Sentrix_ID <- as.character(targets$Slide)
                } else {
                    targets$Sentrix_ID <- paste0("Slide", seq_len(nrow(targets)))
                }
            }

            if (!("Sentrix_Position" %in% colnames(targets))) {
                if ("Array" %in% colnames(targets)) {
                    targets$Sentrix_Position <- as.character(targets$Array)
                } else {
                    targets$Sentrix_Position <- paste0("Position", seq_len(nrow(targets)))
                }
            }

            if (!("Timepoint" %in% colnames(targets))) {
                targets$Timepoint <- rep(c("1", "2"), length.out = nrow(targets))
            }

            rgset_col_data <- SummarizedExperiment::colData(rgset)
            rgset_col_data$Sex <- targets$Sex
            rgset_col_data$Sentrix_ID <- targets$Sentrix_ID
            rgset_col_data$Sentrix_Position <- targets$Sentrix_Position
            rgset_col_data$Timepoint <- targets$Timepoint
            SummarizedExperiment::colData(rgset) <- rgset_col_data

            list(
                RGSet = rgset,
                targets = targets,
                crossReactivePath = system.file(
                    "extdata",
                    "12864_2024_10027_MOESM8_ESM.csv",
                    package = "dnaEPICO"
                )
            )
        }
    )
}

exampleMinfiIdatInputsDnaEpico <- function(n = 6L) {
    getCachedExampleDnaEpico(
        key = paste0("minfi_idat_", as.integer(n)),
        builder = function() {
            requireExamplePackageDnaEpico("minfiData")
            requireExamplePackageDnaEpico("IlluminaHumanMethylation450kmanifest")
            requireExamplePackageDnaEpico("IlluminaHumanMethylation450kanno.ilmn12.hg19")

            base_dir <- system.file("extdata", package = "minfiData")
            targets <- utils::read.csv(
                file.path(base_dir, "SampleSheet.csv"),
                stringsAsFactors = FALSE,
                skip = 7
            )
            targets <- utils::head(targets, n = as.integer(n))
            targets$Basename <- paste0(targets$Sentrix_ID, "_", targets$Sentrix_Position)

            if (!("Sex" %in% colnames(targets)) && "sex" %in% colnames(targets)) {
                targets$Sex <- targets$sex
            }

            temp_dir <- exampleTempDirDnaEpico("dnaEPICO-idat-example-")
            idat_dir <- file.path(temp_dir, "idats")
            dir.create(idat_dir, recursive = TRUE, showWarnings = FALSE)

            available_idats <- list.files(
                base_dir,
                pattern = "[.]idat$",
                recursive = TRUE,
                full.names = TRUE
            )

            for (basename_value in targets$Basename) {
                matched_files <- grep(
                    paste0(basename_value, ".*[.]idat$"),
                    available_idats,
                    value = TRUE
                )
                file.copy(matched_files, idat_dir, overwrite = TRUE)
            }

            pheno_file <- file.path(temp_dir, "pheno.csv")
            utils::write.csv(targets, pheno_file, row.names = FALSE)

            list(
                tempDir = temp_dir,
                idatFolder = idat_dir,
                phenoFile = pheno_file,
                targets = targets,
                arrayType = "IlluminaHumanMethylation450k",
                annotationVersion = "ilmn12.hg19",
                crossReactivePath = system.file(
                    "extdata",
                    "12864_2024_10027_MOESM8_ESM.csv",
                    package = "dnaEPICO"
                )
            )
        }
    )
}

exampleMinfiWorkflowStateDnaEpico <- function() {
    getCachedExampleDnaEpico(
        key = "minfi_workflow",
        builder = function() {
            example_data <- exampleMinfiBaseDataDnaEpico()

            rgset <- example_data$RGSet
            targets <- example_data$targets

            raw_data <- buildRawMinfiEwasWater(
                RGSet = rgset,
                verbose = FALSE,
                logs = FALSE
            )
            assessment <- assessSamplesMinfiEwasWater(
                rawData = raw_data,
                RGSet = rgset,
                detPThreshold = 1,
                verbose = FALSE,
                logs = FALSE
            )
            sample_data <- filterSamplesMinfiEwasWater(
                RGSet = rgset,
                targets = targets,
                failedSamples = character(0),
                SampleID = "Sample_Name",
                verbose = FALSE,
                logs = FALSE
            )
            raw_filtered <- buildRawMinfiEwasWater(
                RGSet = sample_data$RGSet,
                verbose = FALSE,
                logs = FALSE
            )
            sex_data <- predictSexMinfiEwasWater(
                rawData = raw_filtered,
                targets = sample_data$targets,
                SampleID = "Sample_Name",
                sexColumn = "Sex",
                verbose = FALSE,
                logs = FALSE
            )

            rgset_col_data <- SummarizedExperiment::colData(sample_data$RGSet)
            rgset_col_data$Sex <- sample_data$targets$Sex
            rgset_col_data$PredSex <- sample_data$targets$PredSex
            SummarizedExperiment::colData(sample_data$RGSet) <- rgset_col_data

            norm_data <- normalizeMinfiEwasWater(
                sampleData = sample_data,
                sexColumn = "Sex",
                normMethods = "quantile",
                verbose = FALSE,
                logs = FALSE
            )
            filtered_data <- filterProbesMinfiEwasWater(
                normData = norm_data,
                RGSet = sample_data$RGSet,
                pvalThreshold = 1,
                chrToRemove = "chrY",
                snpsToRemove = "SBE",
                mafThreshold = 1,
                crossReactivePath = example_data$crossReactivePath,
                detPtype = "m+u",
                verbose = FALSE,
                logs = FALSE
            )
            metrics_data <- extractMetricsMinfiEwasWater(
                filteredData = filtered_data,
                verbose = FALSE,
                logs = FALSE
            )

            list(
                RGSet = rgset,
                targets = targets,
                rawData = raw_data,
                assessment = assessment,
                sampleData = sample_data,
                rawFiltered = raw_filtered,
                sexData = sex_data,
                normData = norm_data,
                filteredData = filtered_data,
                metricsData = metrics_data,
                crossReactivePath = example_data$crossReactivePath
            )
        }
    )
}

examplePreprocessingPhenoStateDnaEpico <- function() {
    getCachedExampleDnaEpico(
        key = "preprocessing_pheno",
        builder = function() {
            temp_dir <- exampleTempDirDnaEpico("dnaEPICO-preprocessingPheno-example-")

            pheno <- data.frame(
                Sample_Name = c("S1", "S2", "S3", "S4"),
                Timepoint = c("1", "1", "2", "2"),
                Sex = c("Female", "Male", "Female", "Male"),
                Age = c(20, 22, 21, 23),
                stringsAsFactors = FALSE
            )

            beta <- matrix(
                c(0.20, 0.25, 0.22, 0.27, 0.60, 0.55, 0.52, 0.58, 0.10, 0.15, 0.12, 0.16),
                nrow = 3,
                byrow = TRUE,
                dimnames = list(c("cg00000029", "cg00000108", "cg00000109"), pheno$Sample_Name)
            )
            m <- log2(beta / (1 - beta))
            cn <- matrix(
                1,
                nrow = nrow(beta),
                ncol = ncol(beta),
                dimnames = dimnames(beta)
            )

            pheno_path <- file.path(temp_dir, "phenoLC.csv")
            beta_path <- file.path(temp_dir, "beta.RData")
            m_path <- file.path(temp_dir, "m.RData")
            cn_path <- file.path(temp_dir, "cn.RData")

            utils::write.csv(pheno, pheno_path, row.names = FALSE)
            saveNamedObjectMinfiEwasWater(beta, "beta", beta_path)
            saveNamedObjectMinfiEwasWater(m, "m", m_path)
            saveNamedObjectMinfiEwasWater(cn, "cn", cn_path)

            metrics_data <- loadMetricsPreprocessingPheno(
                betaPath = beta_path,
                mPath = m_path,
                cnPath = cn_path,
                verbose = FALSE,
                logs = FALSE
            )
            timepoint_data <- splitTimepointsPreprocessingPheno(
                pheno = pheno,
                metricsData = metrics_data,
                SampleID = "Sample_Name",
                timeVar = "Timepoint",
                timepoints = "1,2",
                verbose = FALSE,
                logs = FALSE
            )
            combined_data <- combineTimepointsPreprocessingPheno(
                timepointData = timepoint_data,
                combineTimepoints = "1,2",
                verbose = FALSE,
                logs = FALSE
            )
            clock_foundation <- buildClockFoundationInputsPreprocessingPheno(
                beta = timepoint_data$data[["1"]]$beta,
                pheno = timepoint_data$data[["1"]]$pheno,
                SampleID = "Sample_Name",
                sexColumn = "Sex",
                verbose = FALSE,
                logs = FALSE
            )

            preprocessing_data <- list(
                pheno = pheno,
                metrics = metrics_data,
                timepointData = timepoint_data,
                combinedData = combined_data,
                clockFoundation = clock_foundation
            )

            list(
                tempDir = temp_dir,
                pheno = pheno,
                phenoPath = pheno_path,
                betaPath = beta_path,
                mPath = m_path,
                cnPath = cn_path,
                metricsData = metrics_data,
                timepointData = timepoint_data,
                combinedData = combined_data,
                clockFoundation = clock_foundation,
                preprocessingData = preprocessing_data
            )
        }
    )
}

exampleSvaAnalysisStateDnaEpico <- function() {
    getCachedExampleDnaEpico(
        key = "sva_analysis",
        builder = function() {
            example_data <- exampleMinfiBaseDataDnaEpico()
            sample_names <- Biobase::sampleNames(example_data$RGSet)

            sva <- cbind(
                sva1 = seq_along(sample_names) / 10,
                sva2 = rev(seq_along(sample_names)) / 10
            )
            rownames(sva) <- sample_names

            merged_pheno <- mergeSvaTargetsEnmix(
                targets = example_data$targets,
                sva = sva,
                SampleID = "Sample_Name",
                verbose = FALSE,
                logs = FALSE
            )
            analysis_data <- analyzeSvaEnmix(
                sva = sva,
                RGSet = example_data$RGSet,
                SentrixIDColumn = "Sentrix_ID",
                SentrixPositionColumn = "Sentrix_Position",
                verbose = FALSE,
                logs = FALSE
            )

            list(
                RGSet = example_data$RGSet,
                targets = example_data$targets,
                sva = sva,
                mergedPheno = merged_pheno,
                analysisData = analysis_data
            )
        }
    )
}

exampleMethylationGLMStateDnaEpico <- function() {
    getCachedExampleDnaEpico(
        key = "methylation_glm",
        builder = function() {
            temp_dir <- exampleTempDirDnaEpico("dnaEPICO-glm-example-")
            input_path <- file.path(temp_dir, "phenoBT1.RData")

            phenoBT1 <- data.frame(
                Sample_Name = c("S1", "S2", "S3", "S4"),
                status = factor(c("Case", "Case", "Control", "Control")),
                sex = factor(c("F", "M", "F", "M")),
                age = c(20, 22, 23, 24),
                cg00000029 = c(0.20, 0.25, 0.22, 0.27),
                cg00000108 = c(0.60, 0.55, 0.52, 0.58),
                check.names = FALSE
            )
            save(phenoBT1, file = input_path)

            prepared_data <- prepareMethylationGLM_T1Data(
                inputPheno = input_path,
                phenotypes = "status",
                covariates = "sex,age",
                factorVars = "status,sex",
                cpgLimit = 2,
                verbose = FALSE,
                logs = FALSE
            )
            model_results <- fitMethylationGLM_T1Models(
                preparedData = prepared_data,
                nCores = 1,
                verbose = FALSE,
                logs = FALSE
            )
            model_summaries <- summarizeMethylationGLM_T1Models(
                modelResults = model_results,
                preparedData = prepared_data,
                summaryResidualSD = TRUE,
                summaryPval = NA,
                nCores = 1,
                verbose = FALSE,
                logs = FALSE
            )

            list(
                tempDir = temp_dir,
                inputPath = input_path,
                preparedData = prepared_data,
                modelResults = model_results,
                modelSummaries = model_summaries
            )
        }
    )
}

exampleMethylationGLMMStateDnaEpico <- function() {
    getCachedExampleDnaEpico(
        key = "methylation_glmm",
        builder = function() {
            temp_dir <- exampleTempDirDnaEpico("dnaEPICO-glmm-example-")
            input_path <- file.path(temp_dir, "phenoBT1T2.RData")

            phenoBT1T2 <- data.frame(
                Sample_Name = paste0("S", seq_len(8)),
                person = rep(seq_len(4), each = 2),
                Timepoint = rep(c("T1", "T2"), 4),
                score = c(10, 11, 8, 9, 7, 8, 12, 13),
                sex = factor(rep(c("F", "M"), each = 4)),
                cg00000029 = c(0.20, 0.22, 0.25, 0.24, 0.18, 0.19, 0.28, 0.29),
                cg00000108 = c(0.60, 0.59, 0.55, 0.56, 0.52, 0.53, 0.58, 0.60),
                check.names = FALSE
            )
            save(phenoBT1T2, file = input_path)

            prepared_data <- prepareMethylationGLMM_T1T2Data(
                inputPheno = input_path,
                personVar = "person",
                timeVar = "Timepoint",
                phenotypes = "score",
                covariates = "sex",
                factorVars = "sex,Timepoint",
                cpgLimit = 2,
                verbose = FALSE,
                logs = FALSE
            )
            model_results <- fitMethylationGLMM_T1T2Models(
                preparedData = prepared_data,
                nCores = 1,
                verbose = FALSE,
                logs = FALSE
            )
            model_summaries <- summarizeMethylationGLMM_T1T2Models(
                modelResults = model_results,
                preparedData = prepared_data,
                summaryPval = NA,
                nCores = 1,
                verbose = FALSE,
                logs = FALSE
            )

            list(
                tempDir = temp_dir,
                inputPath = input_path,
                preparedData = prepared_data,
                modelResults = model_results,
                modelSummaries = model_summaries
            )
        }
    )
}

exampleDnamReportStateDnaEpico <- function() {
    getCachedExampleDnaEpico(
        key = "dnam_report",
        builder = function() {
            temp_dir <- exampleTempDirDnaEpico("dnaEPICO-report-example-")
            qc_dir <- file.path(temp_dir, "qc")
            preprocessing_dir <- file.path(temp_dir, "preprocessing")
            postprocessing_dir <- file.path(temp_dir, "postprocessing")
            sva_dir <- file.path(temp_dir, "sva")
            glm_dir <- file.path(temp_dir, "glm")
            glmm_dir <- file.path(temp_dir, "glmm")
            fig_dir <- file.path(temp_dir, "report-figures")

            createDummyImageDnaEpico(file.path(qc_dir, "qc.png"))
            createDummyImageDnaEpico(file.path(preprocessing_dir, "preprocessing.png"))
            createDummyImageDnaEpico(file.path(postprocessing_dir, "postprocessing.png"))
            createDummyImageDnaEpico(file.path(sva_dir, "sva.png"))
            createDummyImageDnaEpico(file.path(glm_dir, "glm.png"))
            createDummyImageDnaEpico(file.path(glmm_dir, "glmm.png"))

            prepared_report <- prepareDnamReportInputs(
                output = "DNAm_Report.pdf",
                outputDir = temp_dir,
                qcDir = qc_dir,
                preprocessingDir = preprocessing_dir,
                postprocessingDir = postprocessing_dir,
                svaDir = sva_dir,
                glmDir = glm_dir,
                glmmDir = glmm_dir,
                figDir = fig_dir,
                verbose = FALSE,
                logs = FALSE
            )

            list(
                tempDir = temp_dir,
                preparedReport = prepared_report
            )
        }
    )
}
