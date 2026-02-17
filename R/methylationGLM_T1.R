#' Run methylationGLM_T1.R
#' @import parallel
#' @import glm2
#' @import ggplot2
#' @import ggrepel
#' @import minfi
#'
#' @param inputPheno Character. Path to merged phenotype and beta-value RData file.
#' @param outputLogs Character. Directory for log files.
#' @param outputRData Character. Directory for model RData outputs.
#' @param outputPlots Character. Directory for output plots.
#' @param phenotypes Character. Comma-separated phenotype variables to model.
#' @param covariates Character. Comma-separated covariate variables.
#' @param factorVars Character. Comma-separated factor variables.
#' @param cpgPrefix Character. CpG identifier prefix.
#' @param cpgLimit Integer or NA. Maximum number of CpGs to analyse.
#' @param nCores Integer. Number of CPU cores to use.
#' @param plotWidth Integer. Plot width in pixels.
#' @param plotHeight Integer. Plot height in pixels.
#' @param plotDPI Integer. Plot resolution in DPI.
#' @param interactionTerm Character or NULL. Optional interaction term.
#' @param libPath Character or NULL. Optional library path.
#' @param glmLibs Character. GLM libraries to use.
#' @param prsMap Character or NULL. Optional PRS mapping file.
#' @param summaryPval Numeric or NA. P-value threshold for summaries.
#' @param summaryResidualSD Logical. Whether to summarise residual SD.
#' @param saveSignificantCpGs Logical. Save significant CpGs.
#' @param significantCpGDir Character. Directory for significant CpGs.
#' @param significantCpGPval Numeric. P-value threshold for significant CpGs.
#' @param saveTxtSummaries Logical. Save text summaries.
#' @param chunkSize Integer. Number of CpGs processed per chunk.
#' @param summaryTxtDir Character. Directory for summary text files.
#' @param fdrThreshold Numeric. FDR threshold.
#' @param padjmethod Character. P-value adjustment method.
#' @param annotationPackage Character. Annotation package name.
#' @param annotationCols Character. Annotation columns to extract.
#' @param annotatedGLMOut Character. Output directory for annotated results.
#'
#' @return
#' Invisibly returns \code{NULL}. This function is called for its side effects,
#' running CpG-level GLM analyses and writing model results, summaries,
#' plots, and annotated tables to disk.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' methylationGLM_T1(
#'   inputPheno = "rData/preprocessingPheno/mergeData/phenoBetaT1.RData",
#'   outputLogs = "logs",
#'   outputRData = "rData/methylationGLM_T1/models",
#'   outputPlots = "figures/methylationGLM_T1",
#'   phenotypes = "DASS_Depression,DASS_Anxiety,DASS_Stress,
#'                 PCL5_TotalScore,MHCSF_TotalScore,BRS_TotalScore",
#'   covariates = "Sex,Age,Ethnicity,TraumaDefinition,
#'                 Leukocytes,Epithelial.cells",
#'   factorVars = "Sex,Ethnicity,TraumaDefinition",
#'   cpgPrefix = "cg",
#'   cpgLimit = NA,
#'   nCores = 4,
#'   plotWidth = 2000,
#'   plotHeight = 1000,
#'   plotDPI = 150,
#'   interactionTerm = NULL,
#'   libPath = NULL,
#'   glmLibs = "glm2",
#'   prsMap = NULL,
#'   summaryPval = NA,
#'   summaryResidualSD = TRUE,
#'   saveSignificantCpGs = FALSE,
#'   significantCpGDir = "preliminaryResults/cpgs/methylationGLM_T1",
#'   significantCpGPval = 0.05,
#'   saveTxtSummaries = TRUE,
#'   chunkSize = 10000,
#'   summaryTxtDir = "preliminaryResults/summary/methylationGLM_T1/glm",
#'   fdrThreshold = 0.05,
#'   padjmethod = "fdr",
#'   annotationPackage =
#'     "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
#'   annotationCols =
#'     "Name,chr,pos,UCSC_RefGene_Group,
#'      UCSC_RefGene_Name,Relation_to_Island,
#'      GencodeV41_Group",
#'   annotatedGLMOut = "data/methylationGLM_T1"
#' )
#' }
#'
#' @export
methylationGLM_T1 <- function(
    inputPheno = "rData/preprocessingPheno/mergeData/phenoBetaT1.RData",
    outputLogs = "logs",
    outputRData = "rData/methylationGLM_T1/models",
    outputPlots = "figures/methylationGLM_T1",
    phenotypes = "DASS_Depression,DASS_Anxiety,DASS_Stress,PCL5_TotalScore,
                  MHCSF_TotalScore,BRS_TotalScore",
    covariates = "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells",
    factorVars = "Sex,Ethnicity,TraumaDefinition",
    cpgPrefix = "cg",
    cpgLimit = NA,
    nCores = 32,
    plotWidth = 2000,
    plotHeight = 1000,
    plotDPI = 150,
    interactionTerm = NULL,
    libPath = NULL,
    glmLibs = "glm2",
    prsMap = NULL,
    summaryPval = NA,
    summaryResidualSD = TRUE,
    saveSignificantCpGs = FALSE,
    significantCpGDir = "preliminaryResults/cpgs/methylationGLM_T1",
    significantCpGPval = 0.05,
    saveTxtSummaries = TRUE,
    chunkSize = NULL,
    summaryTxtDir = "preliminaryResults/summary/methylationGLM_T1/glm",
    fdrThreshold = 0.05,
    padjmethod = "fdr",
    annotationPackage = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    annotationCols = "Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,
                      Relation_to_Island,GencodeV41_Group",
    annotatedGLMOut = "data/methylationGLM_T1"
) {

# Set default values for optional arguments
if (is.character(cpgLimit) && tolower(cpgLimit) == "na") {
        cpgLimit <- NA
}

if (is.character(summaryPval) && tolower(summaryPval) == "na") {
        summaryPval <- NA
}

if (!is.null(prsMap)) {
        prsMapList <- setNames(
                sapply(strsplit(unlist(strsplit(prsMap, ",")), ":"), `[`, 2),
                sapply(strsplit(unlist(strsplit(prsMap, ",")), ":"), `[`, 1)
        )
} else {
        prsMapList <- list()
}


dir.create(annotatedGLMOut, recursive = TRUE, showWarnings = FALSE)
dir.create(summaryTxtDir, recursive = TRUE, showWarnings = FALSE)
dir.create(significantCpGDir, recursive = TRUE, showWarnings = FALSE)
dir.create(outputRData, recursive = TRUE, showWarnings = FALSE)
dir.create(outputPlots, recursive = TRUE, showWarnings = FALSE)

glmLibList   <- strsplit(glmLibs, ",")[[1]]
covariateList <- strsplit(covariates, ",")[[1]]
phenotypeList <- strsplit(phenotypes, ",")[[1]]
factorVarsList <- strsplit(factorVars, ",")[[1]]
# ==============================================================================

# ----------- Setup Logging -----------
dir.create(outputLogs, recursive = TRUE, showWarnings = FALSE)
logFilePath <- file.path(outputLogs, "log_methylationGLM_T1.txt")
logCon <- file(logFilePath, open = "wt")
sink(logCon, split = TRUE)
sink(logCon, type = "message")
# ==============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting DNAm GLM Analysis (Timepoint 1) ====\n")
cat("Start time: ", format(Sys.time()), "\n")
cat("Input phenotype + beta file: ", inputPheno, "\n")
cat("Output RData folder: ", outputRData, "\n")
cat("Output logs folder: ", outputLogs, "\n")
cat("Output plots folder: ", outputPlots, "\n")
cat("Phenotypes: ", phenotypes, "\n")
cat("Covariates: ", covariates, "\n")
cat("Factor variables: ", factorVars, "\n")
cat("CpG column prefix: ", cpgPrefix, "\n")
cat("CpG limit: ", ifelse(is.na(cpgLimit), "All", cpgLimit), "\n")
cat("Number of cores: ", nCores, "\n")
cat("Interaction term: ", interactionTerm, "\n")
cat("Library path: ", libPath, "\n")
cat("GLM libraries: ", glmLibs, "\n")
cat("PRS mapping: ", ifelse(is.null(prsMap), "None", paste(prsMap, collapse = ", ")), "\n")
cat("Include Residual SD in summary: ", summaryResidualSD, "\n")
cat("Summary p-value filter: ", ifelse(is.na(summaryPval), "None", summaryPval), "\n")
cat("Save summary tables: ", saveTxtSummaries, "\n")
cat("Chunk size for parallel processing: ", ifelse(is.null(chunkSize), "Auto", chunkSize), "\n")
cat("Summary output folder: ", summaryTxtDir, "\n")
cat("FDR threshold: ", fdrThreshold, "\n")
cat("P-value adjustment method: ", padjmethod, "\n")
cat("Save significant CpGs: ", saveSignificantCpGs, "\n")
cat("Significant CpG output folder: ", significantCpGDir, "\n")
cat("Significance p-value cutoff: ", significantCpGPval, "\n")
cat("Annotation package: ", annotationPackage, "\n")
cat("Annotation columns: ", annotationCols, "\n")
cat("Annotated output CSV: ", annotatedGLMOut, "\n")
cat("=======================================================================\n")

# ----------- Load Data -----------
inpt <- load(inputPheno)
phenoBT1 <- get(inpt)

cat("Loaded phenotype + beta data from:", inputPheno, "\n")
phenotypes <- strsplit(phenotypes, ",")[[1]]
covariates <- strsplit(covariates, ",")[[1]]
factorVars <- strsplit(factorVars, ",")[[1]]

# ----------- Analysis Info -----------
cat("phenoBT1 dimensions:", dim(phenoBT1), "\n")
cgCount <- sum(grepl("^cg", colnames(phenoBT1)))
chCount <- sum(grepl("^ch", colnames(phenoBT1)))
cat("Number of CpG columns (cg):", cgCount, "\n")
cat("Number of CpG columns (ch):", chCount, "\n")

cat("Checking structure of the dataset...\n")
if (!is.null(interactionTerm) && interactionTerm != "") {
  if (interactionTerm %in% names(phenoBT1)) {
    print(table(phenoBT1[[interactionTerm]], useNA = "ifany"))
  } else {
    cat("Warning: interactionTerm not found in phenoBT1 columns.\n")
  }
} else {
  cat("No interaction term specified; skipping table summary.\n")
}
cat("=======================================================================\n")

# ----------- Convert Factors -----------
for (var in factorVars) {
        if (var %in% colnames(phenoBT1)) {
                phenoBT1[[var]] <- as.factor(phenoBT1[[var]])
        }
}

# ----------- Missing Summary & Distribution -----------
cat("Missing summary:\n")
print(sapply(phenoBT1[, c(phenotypes, covariates)], function(x) sum(is.na(x))))
cat("=======================================================================\n")

cat("Summary statistics:\n")
print(summary(phenoBT1[, c(phenotypes, covariates)]))

for (var in phenotypes) {

        if (is.numeric(phenoBT1[[var]])) {
          p <- ggplot(phenoBT1, aes_string(x = var)) +
            geom_histogram(bins = 30, fill = "steelblue", color = "white") +
            labs(title = paste("Distribution of", var),
                 x = var, y = "Frequency") +
            theme_minimal()

        } else {
          p <- ggplot(phenoBT1, aes_string(x = var)) +
            geom_bar(fill = "steelblue") +
            labs(title = paste("Distribution of", var), x = var, y = "Count") +
            theme_minimal()
        }

        tiff(filename = file.path(outputPlots, paste0("hist_", var, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")

        print(p)
        dev.off()
}

catVars <- intersect(factorVars, colnames(phenoBT1))
for (var in catVars) {
        p <- ggplot(phenoBT1, aes_string(x = var)) +
                geom_bar(fill = "darkorange") +
                labs(title = paste("Distribution of",
                                   var), x = var, y = "Count") +
                theme_minimal()
        tiff(filename = file.path(outputPlots, paste0("bar_", var, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")
        print(p)
        dev.off()
}

numVars <- setdiff(covariates, factorVars)
numVars <- intersect(numVars, colnames(phenoBT1))

for (var in numVars) {
        p <- ggplot(phenoBT1, aes_string(x = var)) +
                geom_histogram(bins = 30, fill = "darkgreen", color = "white") +
                labs(title = paste("Distribution of", var),
                     x = var, y = "Frequency") +
                theme_minimal()

        tiff(filename = file.path(outputPlots, paste0("hist_", var, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")
        print(p)
        dev.off()
}

cat("Plots saved to:", outputPlots, "\n")
cat("=======================================================================\n")

# ----------- GLM Function -----------
glm <- function(
                phenoScore,
                merge,
                covariates = covariateList,
                factorVars = factorVarsList,
                cpgPrefix = cpgPrefix,
                cpgLimit = cpgLimit,
                nCore = nCores,
                libPath = libPath,
                glmLibList = glmLibList,
                interactionTerm = NULL

) {
        if (!phenoScore %in% colnames(merge)) {
                stop(paste("Phenotype", phenoScore, "not found in the dataset"))
        }

        if (!is.null(interactionTerm) && interactionTerm != "" && interactionTerm %in% colnames(merge)) {
          interactionPart <- paste(phenoScore, "*", interactionTerm)
          fixedTerms <- setdiff(c(covariates), c(phenoScore, interactionTerm))
          fullFormula <- paste("~", paste(c(interactionPart, fixedTerms), collapse = " + "))
        } else {
          fullFormula <- paste("~", paste(c(phenoScore, covariates), collapse = " + "))
        }

        cpgCol <- grep(paste0("^", cpgPrefix), colnames(merge), value = TRUE)
        if (!is.na(cpgLimit)) {
                cpgCol <- head(cpgCol, as.numeric(cpgLimit))
        }
        cl <- makeCluster(nCore)
        clusterExport(
                cl,
                varlist = c("phenoScore", "interactionTerm",
                            "covariates",
                            "factorVars",
                            "libPath",
                            "glmLibList", "fullFormula", "chunkSize"),
                envir = environment()
        )

        # clusterEvalQ(cl, {
        #         if (!is.null(libPath)) {
        #                 .libPaths(libPath)
        #         }
        #         sapply(glmLibList, function(pkg) {
        #                 if (!require(pkg, character.only = TRUE)) {
        #                         stop(paste("Failed to load package:", pkg))
        #                 }
        #         })
        # })

        clusterEvalQ(cl, {
                if (!is.null(libPath)) {
                        .libPaths(libPath)
                }

                for (pkg in glmLibList) {
                        if (!requireNamespace(pkg, quietly = TRUE)) {
                        stop(paste("Failed to load package:", pkg))
                        }
                }

                NULL
        })


        fit <- function(cpg) {
                tryCatch({
                        vars <- unique(c(phenoScore,
                                                 covariates,
                                                 cpg,
                                                 if (!is.null(interactionTerm) && interactionTerm != "") interactionTerm))

                        model <- merge[, vars, drop = FALSE]
                        colnames(model)[ncol(model)] <- "beta"

                        for (var in factorVars) {
                                if (var %in% colnames(model)) {
                                        model[[var]] <- as.factor(model[[var]])
                                }
                        }

                        fit <- glm2(
                                formula = as.formula(paste("beta", fullFormula)),
                                data = model,
                                family = gaussian(),
                                na.action = na.exclude
                        )

                        return(list(
                                coef = summary(fit)$coefficients,
                                residuals = residuals(fit),
                                fitted = fitted(fit)
                        ))
                }, error = function(e) NULL)
        }

        fitList <- parLapply(cl, cpgCol, fit)
        names(fitList) <- cpgCol
        stopCluster(cl)

        return(fitList)
}

# ----------- Run GLMs -----------
cat("Running GLMs for all phenotypes...\n")
for (pheno in phenotypeList) {
        outputFile <- file.path(outputRData, paste0(pheno, "GLM.RData"))

        if (file.exists(outputFile)) {
                cat("Loading existing GLM for:", pheno, "\n")
                load(outputFile)
                assign(paste0(pheno, "GLM"), fitResult, envir = .GlobalEnv)
                next
        }

        cat("Running GLM for:", pheno, "\n")

        # Include PRS if defined for this phenotype
        prsVar <- if (pheno %in% names(prsMapList)) prsMapList[[pheno]] else NULL
        allCovariates <- if (!is.null(prsVar)) c(covariateList, prsVar) else covariateList

        # Interaction
        modelFormula <- if (!is.null(interactionTerm) && interactionTerm != "") {
          paste("~", paste(c(paste0(pheno, "*", interactionTerm), allCovariates), collapse = " + "))
        } else {
          paste("~", paste(c(pheno, allCovariates), collapse = " + "))
        }
        cat("Formula:", modelFormula, "\n")

        fitResult <- glm(
                phenoScore = pheno,
                merge = phenoBT1,
                covariates = allCovariates,
                factorVars = factorVarsList,
                cpgPrefix = cpgPrefix,
                cpgLimit = cpgLimit,
                nCore = nCores,
                libPath = libPath,
                glmLibList = glmLibList,
                interactionTerm = interactionTerm
        )

        save(fitResult, file = outputFile)

        assign(paste0(pheno, "GLM"), fitResult, envir = .GlobalEnv)

}
cat("Finished running GLMs for all phenotypes.\n")
cat("=======================================================================\n")

# ----------- Extract GLM Summary in Parallel -----------
cpgsGLM <- function(
                fitList,
                variable,
                interactionTerm = NULL,
                includeResidualSD = summaryResidualSD,
                pValue = summaryPval,
                nCore = nCores,
                libPath = libPath,
                glmLibList = glmLibList,
                chunkSize = NULL

) {

  cat("Starting extraction for variable:", variable, "\n")
  cpgNames <- names(fitList)
  if (is.null(chunkSize)) {
    chunkSize <- max(10, floor(length(cpgNames) / (nCore * 4)))
  }
  cat("Total CpGs:", length(cpgNames), "| Using chunkSize:", chunkSize, "\n")

  splitIntoChunks <- function(x, size) {
    split(x, ceiling(seq_along(x) / size))
  }

  cpgNames <- names(fitList)
  cpgChunks <- splitIntoChunks(cpgNames, chunkSize)

  cl <- makeCluster(nCore)
  clusterExport(
    cl,
    varlist = c("fitList",
                "variable",
                "interactionTerm",
                "includeResidualSD",
                "pValue",
                "libPath",
                "glmLibList"),
    envir = environment()
  )

  clusterEvalQ(cl, {
    if (!is.null(libPath)) {
        .libPaths(libPath)
    }

    for (pkg in glmLibList) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("Failed to load package:", pkg))
        }
    }

    NULL
  })
 

  results <- parLapplyLB(cl, cpgChunks, function(chunk) {
    outList <- vector("list", length(chunk))
    idx <- 1
    for (cpg in chunk) {
      modelObj <- fitList[[cpg]]
      if (is.null(modelObj) || is.null(modelObj$coef)) next

      coefTable <- modelObj$coef
      if (!is.null(interactionTerm) && interactionTerm != "") {
        pattern <- paste0("^",
                          variable, ".*:", interactionTerm)
        matchedRows <- grep(pattern, rownames(coefTable), value = TRUE)

        if (length(matchedRows) == 0) {
          matchedRows <- grep(paste0("^", variable),
                              rownames(coefTable), value = TRUE)
        }
      } else {
        matchedRows <- grep(paste0("^", variable),
                            rownames(coefTable), value = TRUE)
      }

      if (length(matchedRows) == 0) next

      tmp <- coefTable[matchedRows, , drop = FALSE]
      tmp <- as.data.frame(tmp)
      tmp$CpG <- cpg
      tmp$Coefficient <- rownames(tmp)

      if (includeResidualSD && !is.null(modelObj$residuals)) {
        tmp$ResidualSD <- sd(modelObj$residuals, na.rm = TRUE)
      }

      if (!is.null(tmp)) {
        outList[[idx]] <- tmp
        idx <- idx + 1
      }
    }
    do.call(rbind, outList)
  })

  stopCluster(cl)

  summary <- do.call(rbind, results)

  if (is.null(summary) || nrow(summary) == 0) {
    warning("No CpG-level results extracted for variable:", variable)
    return(NULL)
  }

  cols <- c("CpG", "Coefficient", "Estimate", "Std. Error", "t value",
            "Pr(>|t|)")
  if (includeResidualSD) cols <- c(cols, "ResidualSD")
  summary <- summary[, cols, drop = FALSE]

  if (!is.na(pValue)) {
    summary <- subset(summary, `Pr(>|t|)` < pValue)
  }

  rownames(summary) <- NULL

  cat("Completed GLM summary extraction for:", variable, "\n")
  return(summary)
}

cat("=======================================================================\n")
phenotypeList <- strsplit(phenotypes, ",")[[1]]

# ----------- Run and Save GLM Summaries -----------
for (pheno in phenotypeList) {
        summaryFile <- file.path(outputRData, paste0(pheno,
                                                         "SummaryGLM.RData"))

        if (file.exists(summaryFile)) {
                cat("Loading existing summary for:", pheno, "\n")
                load(summaryFile)
                assign(paste0(pheno, "SummaryGLM"), fitResult,
                       envir = .GlobalEnv)
                next
        }

        cat("Running GLM summary extraction for:", pheno, "\n")

        fitObjectName <- paste0(pheno, "GLM")
        fitObject <- get(fitObjectName)

        fitResult <- cpgsGLM(
                fitList = fitObject,
                variable = pheno,
                includeResidualSD = summaryResidualSD,
                pValue = summaryPval,
                nCore = nCores,
                libPath = libPath,
                glmLibList = glmLibList,
                chunkSize = chunkSize
        )

        save(fitResult, file = summaryFile)

        assign(paste0(pheno, "SummaryGLM"), fitResult, envir = .GlobalEnv)

        cat("Saved:", paste0(pheno, "SummaryGLM.RData"), "\n")
}
cat("=======================================================================\n")

# ----------- Save Significant CpGs -----------
saveSignificantCpGs <- function(
                resultList,
                resultName,
                baseDir = significantCpGDir,
                pvalThreshold = significantCpGPval,
                interactionTerm = NULL
) {
  resultDir <- file.path(baseDir, resultName)
  if (!dir.exists(resultDir)) dir.create(resultDir, recursive = TRUE)

  for (i in seq_along(resultList)) {
    coefTable <- resultList[[i]]$coef
    cpgName <- names(resultList)[i]

    if (!is.null(interactionTerm) && interactionTerm != "") {
      pattern <- paste0("^", resultName, ".*:", interactionTerm)
      matchedRows <- grep(pattern, rownames(coefTable), value = TRUE)

      if (length(matchedRows) == 0) {
        matchedRows <- grep(paste0("^", resultName),
                            rownames(coefTable), value = TRUE)
        if (length(matchedRows) > 0) {
          cat("Interaction term", interactionTerm, "dropped for CpG", cpgName,
              "- extracting main effect for", resultName, "\n")
        }
      }
    } else {
      matchedRows <- grep(paste0("^", resultName),
                          rownames(coefTable), value = TRUE)
    }

    if (length(matchedRows) > 0) {
      matchedPvals <- coefTable[matchedRows, "Pr(>|t|)"]

      if (any(matchedPvals < pvalThreshold, na.rm = TRUE)) {
        cpgDir <- file.path(resultDir, cpgName)
        if (!dir.exists(cpgDir)) dir.create(cpgDir)
        outputFile <- file.path(cpgDir, paste0(cpgName, ".txt"))
        write.table(coefTable, file = outputFile, sep = "\t", quote = FALSE)
      }
    }
  }
}

# ---------- Save Significant CpGs to Directory -----------
if (isTRUE(saveSignificantCpGs)) {
        cat("Saving significant CpGs to:", significantCpGDir, "\n")

        for (pheno in strsplit(phenotypes, ",")[[1]]) {
                modelObjName <- paste0(pheno, "GLM")
                if (exists(modelObjName)) {
                        resultObj <- get(modelObjName)
                        saveSignificantCpGs(
                                resultList = resultObj,
                                resultName = pheno
                        )
                }
        }

        cat("Finished saving significant CpG model outputs.\n")
}
cat("=======================================================================\n")

# ----------- Save GLM Summary Tables -----------
saveSummaryToTxt <- function(
                summaryDF,
                outputFile,
                sortByPval = TRUE
) {
        if (sortByPval) {
                if ("P.value" %in% colnames(summaryDF)) {
                        summaryDF <- summaryDF[order(summaryDF$P.value), ]
                } else if ("Pr(>|t|)" %in% colnames(summaryDF)) {
                        summaryDF <- summaryDF[order(summaryDF[["Pr(>|t|)"]]), ]
                }
        }

        dir.create(dirname(outputFile), recursive = TRUE, showWarnings = FALSE)

        write.table(
                summaryDF,
                file = outputFile,
                sep = "\t",
                row.names = FALSE,
                quote = FALSE
        )

        cat("Saved summary table to:", outputFile, "\n")
}

# ----------- Save All Summaries from Phenotype List -----------

if (saveTxtSummaries) {
        cat("Saving summaries to TXT files...\n")

        phenoList <- strsplit(phenotypes, ",")[[1]]

        for (pheno in phenoList) {

                outputFile <- file.path(
                        summaryTxtDir,
                        paste0(pheno, "SummaryGLM.txt")
                )

                summaryObjName <- paste0(pheno, "SummaryGLM")

                if (exists(summaryObjName)) {
                        summaryObj <- get(summaryObjName)

                        saveSummaryToTxt(
                                summaryDF = summaryObj,
                                outputFile = outputFile
                        )
                }
        }

        cat("Finished writing all GLM summaries to TXT files.\n")
}
cat("=======================================================================\n")

# ----------- Diagnostic Plot Function -----------
diagnosticPlots <- function(
                summary,
                betaMatrix,
                variable,
                cpgPrefix = cpgPrefix,
                padjmethod = padjmethod,
                fdrThreshold = fdrThreshold,
                outputDir = outputPlots,
                plotWidth = plotWidth,
                plotHeight = plotHeight,
                plotDPI = plotDPI
) {
        cat("Generating diagnostic plots for:", variable, "\n")

        dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

        summary$FDR <- p.adjust(summary$`Pr(>|t|)`, method = padjmethod)

        chisq <- qchisq(1 - summary$`Pr(>|t|)`, df = 1)
        lambda <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1), 3)
        message("Genomic inflation factor: ", lambda)

        # -- Q-Q Plot of P-values --
        pvals <- summary$`Pr(>|t|)`
        pvals <- pvals[!is.na(pvals)]
        tiff(filename = file.path(outputDir,
                                  paste0("qqplot_", variable, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")
        qqplot(
                -log10(ppoints(length(pvals))),
                -log10(sort(pvals)),
                main = paste("Q-Q Plot of p-values for", variable,
                             "\nGenomic Inflation Factor= ", lambda),
                xlab = "Expected -log10(p)",
                ylab = "Observed -log10(p)",
                pch = 16,
                col = "black"
        )
        abline(0, 1, col = "red")
        dev.off()

        # -- Residual SD vs log2(mean Beta) --
        cpgCols <- grep(paste0("^", cpgPrefix),
                        colnames(betaMatrix), value = TRUE)
        betaMat <- betaMatrix[, cpgCols, drop = FALSE]
        meanBeta <- colMeans(betaMat, na.rm = TRUE)
        meanBetaLog <- log2(meanBeta)
        summary$log2meanBeta <- meanBetaLog[summary$CpG]
        summary <- summary[!is.na(summary$log2meanBeta), ]

        tiff(filename = file.path(outputDir,
                                  paste0("residualSD_", variable, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")
        plot(summary$log2meanBeta,
             summary$ResidualSD,
             pch = 20,
             col = rgb(0, 0, 0, 0.4),
             xlab = "log2(Average Beta)",
             ylab = "Residual Standard Deviation",
             main = "plotSA-style: Residual SD vs Average Beta")
        lines(lowess(summary$log2meanBeta, summary$ResidualSD),
              col = "red", lwd = 2)
        dev.off()

        # -- Residual SD vs Significance --
        p <- ggplot(summary,
                    aes(x = -log10(`Pr(>|t|)`),
                        y = ResidualSD,
                        color = FDR < fdrThreshold)) +
                geom_point(alpha = 0.6) +
                geom_text_repel(data = subset(summary,
                                              FDR < fdrThreshold),
                                aes(label = CpG)) +
                scale_color_manual(values = c("FALSE" = "grey50",
                                              "TRUE" = "firebrick")) +
                labs(
                        title = paste("Residual SD vs Significance for",
                                      variable),
                        x = "-log10(p-value)",
                        y = "Residual SD",
                        color = paste("FDR <", fdrThreshold)
                ) +
                theme_minimal()

        # Save as TIFF using tiff() + print() + dev.off()
        tiff(filename = file.path(outputDir,
                                     paste0("residualSignificance_",
                                            variable, ".tiff")),
                width = plotWidth,
                height = plotHeight,
                res = plotDPI, type = "cairo"
        )
        print(p)
        dev.off()


}

# ----------- Generate Diagnostic Plots -----------
cat("Generating diagnostic plots for all phenotypes...\n")
for (pheno in strsplit(phenotypes, ",")[[1]]) {
        summaryName <- paste0(pheno, "SummaryGLM")
        if (exists(summaryName)) {
                diagnosticPlots(
                        summary = get(summaryName),
                        betaMatrix = phenoBT1,
                        variable = pheno,
                        fdrThreshold = fdrThreshold,
                        cpgPrefix = cpgPrefix,
                        padjmethod = padjmethod,
                        outputDir = outputPlots,
                        plotWidth = plotWidth,
                        plotHeight = plotHeight,
                        plotDPI = plotDPI
                )
        }
}
cat("Finished generating all diagnostic plots.\n")
cat("Plots saved to:", outputPlots, "\n")
cat("=======================================================================\n")

# ----------- Load Annotation Data -----------
cat("Loading annotation object:", annotationPackage, "\n")
annotationObject <- getAnnotation(get(annotationPackage))

cat("Annotation loaded with", nrow(annotationObject), "probes\n")
cat("Annotation columns available:\n")
print(colnames(annotationObject))
cat("=======================================================================\n")

# ----------- Apply Annotation -----------
annotateGLM <- function(
                summaryList,
                annotationObject,
                annotationCols = strsplit(annotationCols, ",")[[1]]
) {
        cat("Starting annotation of GLM summaries...\n")
        print(annotationCols)
        annotationCols <- trimws(annotationCols)
        annotationCols <- gsub("\n", "", annotationCols)
        cat("\nCorrected annotation columns:\n")
        print(annotationCols)

        modelNames <- names(summaryList)

        cat("Merging GLM summaries...\n")
        cleanedSummaries <- lapply(modelNames, function(modelName) {
          df <- summaryList[[modelName]]

          coefNames <- unique(df$Coefficient)
          dfList <- lapply(coefNames, function(coefName) {
            subdf <- df[df$Coefficient == coefName, ]
            subdf[[paste0(coefName, "P.Value")]] <- subdf$`Pr(>|t|)`
            subdf <- subdf[, c("CpG", paste0(coefName, "P.Value"))]
            return(subdf)
          })

        # Merge multiple coefficient data.frames by CpG
        mergedCoef <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE),
                             dfList)
        return(mergedCoef)

        })

        mergedSummary <- Reduce(function(x, y) merge(x,
                                                     y,
                                                     by = "CpG",
                                                     all = TRUE),
                                cleanedSummaries)

        annDF <- as.data.frame(annotationObject)
        annDF$CpG <- rownames(annDF)

        cat("Merging annotation with GLM summary...\n")
        annotatedResults <- merge(mergedSummary,
                                  annDF, by = "CpG", all.x = TRUE) 

        if (is.null(annotatedResults) || nrow(annotatedResults) == 0) {
                stop("Annotation merge produced empty result")
                }
                                 
        finalCols <- c("CpG",
                       unlist(lapply(cleanedSummaries, function(df) colnames(df)[-1])),
                       annotationCols)

        cat("\nColumns in annotatedResults:\n")
        print(colnames(annotatedResults))

        cat("\nRequested finalCols:\n")
        print(finalCols)

        missingCols <- setdiff(finalCols, colnames(annotatedResults))
        cat("\nMissing columns:\n")
        print(missingCols)

        annotatedResults <- annotatedResults[, finalCols]
        colnames(annotatedResults)[1] <- "IlmnID"

        cat("Annotation completed. Annotated CpGs:",
            nrow(annotatedResults), "\n")
        return(annotatedResults)
}

# ----------- Execute Annotation and Save Results -----------

cat("Running annotation of GLM summary results...\n")

# Split phenotype names and fetch each corresponding summary object
phenotypeNames <- strsplit(phenotypes, ",")[[1]]
summaryList <- setNames(
        lapply(phenotypeNames, function(pheno) {
                get(paste0(pheno, "SummaryGLM"))
        }),
        phenotypeNames
)

# Convert annotationCols from comma-separated string to vector
annotationColsVec <- strsplit(annotationCols, ",")[[1]]

# Run the annotation function
annotatedGLM <- annotateGLM(
        summaryList = summaryList,
        annotationObject = getAnnotation(get(annotationPackage)),
        annotationCols = annotationColsVec
)

# Save annotated results
cat("Saving annotated GLM summary to:", annotatedGLMOut, "\n")
write.csv(
        annotatedGLM,
        file = file.path(annotatedGLMOut, "annotatedGLM.csv"),
        row.names = FALSE)
cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")
sink()
close(logCon)
}