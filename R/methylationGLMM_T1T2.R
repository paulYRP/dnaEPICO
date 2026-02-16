#' Run methylationGLMM_T1T2.R
#' @import parallel
#' @import ggplot2
#' @import ggrepel
#' @import minfi
#' @importFrom data.table rbindlist
#' @importFrom dplyr select group_by summarise across all_of
#' @import tidyr
#' @import readr
#' @import stringr
#' @importFrom lme4 ranef fixef
#' @import lmerTest
#'
#' @param inputPheno Character. Path to merged longitudinal phenotype and beta-value RData file.
#' @param outputLogs Character. Directory for log files.
#' @param outputRData Character. Directory for model RData outputs.
#' @param outputPlots Character. Directory for output plots.
#' @param personVar Character. Subject identifier variable.
#' @param timeVar Character. Time variable name.
#' @param phenotypes Character. Comma-separated phenotype variables to model.
#' @param covariates Character. Comma-separated covariate variables.
#' @param factorVars Character. Comma-separated factor variables.
#' @param lmeLibs Character. Mixed-effects modelling libraries.
#' @param prsMap Character or NULL. Optional PRS mapping file.
#' @param libPath Character or NULL. Optional library path.
#' @param cpgPrefix Character. CpG identifier prefix.
#' @param cpgLimit Integer or NA. Maximum number of CpGs to analyse.
#' @param nCores Integer. Number of CPU cores to use.
#' @param summaryPval Numeric or NA. P-value threshold for summaries.
#' @param plotWidth Integer. Plot width in pixels.
#' @param plotHeight Integer. Plot height in pixels.
#' @param plotDPI Integer. Plot resolution in DPI.
#' @param interactionTerm Character or NULL. Optional interaction term.
#' @param saveSignificantInteractions Logical. Save significant interaction results.
#' @param significantInteractionDir Character. Directory for significant interactions.
#' @param significantInteractionPval Numeric. P-value threshold for significant interactions.
#' @param saveTxtSummaries Logical. Save text summaries.
#' @param chunkSize Integer. Number of CpGs processed per chunk.
#' @param summaryTxtDir Character. Directory for summary text files.
#' @param fdrThreshold Numeric. FDR threshold.
#' @param padjmethod Character. P-value adjustment method.
#' @param annotationPackage Character. Annotation package name.
#' @param annotationCols Character. Annotation columns to extract.
#' @param annotatedLMEOut Character. Output directory for annotated results.
#'
#' @return
#' Invisibly returns \code{NULL}. This function is called for its side effects,
#' running CpG-level longitudinal mixed-effects models and writing model
#' results, interaction summaries, plots, and annotated tables to disk.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' methylationGLMM_T1T2(
#'   inputPheno = "rData/preprocessingPheno/mergeData/phenoBetaT1T2.RData",
#'   outputLogs = "logs",
#'   outputRData = "rData/methylationGLMM_T1T2/models",
#'   outputPlots = "figures/methylationGLMM_T1T2",
#'   personVar = "person",
#'   timeVar = "Timepoint",
#'   phenotypes = "DASS_Depression,DASS_Anxiety,DASS_Stress,
#'                 PCL5_TotalScore,MHCSF_TotalScore,BRS_TotalScore",
#'   covariates = "Sex,Age,Ethnicity,TraumaDefinition,
#'                 Leukocytes,Epithelial.cells",
#'   factorVars = "Sex,Ethnicity,TraumaDefinition,Timepoint",
#'   lmeLibs = "lme4,lmerTest",
#'   prsMap = NULL,
#'   libPath = NULL,
#'   cpgPrefix = "cg",
#'   cpgLimit = NA,
#'   nCores = 4,
#'   summaryPval = NA,
#'   plotWidth = 2000,
#'   plotHeight = 1000,
#'   plotDPI = 150,
#'   interactionTerm = "TreatvControl_Time1_vs_Time2",
#'   saveSignificantInteractions = TRUE,
#'   significantInteractionDir =
#'     "preliminaryResults/cpgs/methylationGLMM_T1T2",
#'   significantInteractionPval = 0.05,
#'   saveTxtSummaries = TRUE,
#'   chunkSize = 10000,
#'   summaryTxtDir =
#'     "preliminaryResults/summary/methylationGLMM_T1T2/lmer",
#'   fdrThreshold = 0.05,
#'   padjmethod = "fdr",
#'   annotationPackage =
#'     "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
#'   annotationCols =
#'     "Name,chr,pos,UCSC_RefGene_Group,
#'      UCSC_RefGene_Name,Relation_to_Island,
#'      GencodeV41_Group",
#'   annotatedLMEOut = "data/methylationGLMM_T1T2"
#' )
#' }
#'
#' @export
methylationGLMM_T1T2 <- function(
    inputPheno = "rData/preprocessingPheno/mergeData/phenoBetaT1T2.RData",
    outputLogs = "logs/",
    outputRData = "rData/methylationGLMM_T1T2/models",
    outputPlots = "figures/methylationGLMM_T1T2",
    personVar = "person",
    timeVar = "Timepoint",
    phenotypes = "DASS_Depression,DASS_Anxiety,DASS_Stress,PCL5_TotalScore,
                  MHCSF_TotalScore,BRS_TotalScore",
    covariates = "Sex,Age,Ethnicity,TraumaDefinition,Leukocytes,Epithelial.cells",
    factorVars = "Sex,Ethnicity,TraumaDefinition,Timepoint",
    lmeLibs = "lme4,lmerTest",
    prsMap = NULL,
    libPath = NULL,
    cpgPrefix = "cg",
    cpgLimit = NA,
    nCores = 32,
    summaryPval = NA,
    plotWidth = 2000,
    plotHeight = 1000,
    plotDPI = 150,
    interactionTerm = NULL,
    saveSignificantInteractions = TRUE,
    significantInteractionDir = "preliminaryResults/cpgs/methylationGLMM_T1T2",
    significantInteractionPval = 0.05,
    saveTxtSummaries = TRUE,
    chunkSize = 10000,
    summaryTxtDir = "preliminaryResults/summary/methylationGLMM_T1T2/lmer",
    fdrThreshold = 0.05,
    padjmethod = "fdr",
    annotationPackage = "IlluminaHumanMethylationEPICv2anno.20a1.hg38",
    annotationCols = "Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,
                      Relation_to_Island,GencodeV41_Group",
    annotatedLMEOut = "data/methylationGLMM_T1T2"
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

lmeLibList   <- strsplit(lmeLibs, ",")[[1]]
covariateList <- strsplit(covariates, ",")[[1]]
phenotypeList <- strsplit(phenotypes, ",")[[1]]
factorVarsList <- strsplit(factorVars, ",")[[1]]

dir.create(summaryTxtDir, recursive = TRUE, showWarnings = FALSE)
dir.create(significantInteractionDir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotatedLMEOut, recursive = TRUE, showWarnings = FALSE)
dir.create(outputRData, recursive = TRUE, showWarnings = FALSE)
dir.create(outputPlots, recursive = TRUE, showWarnings = FALSE)
# ==============================================================================

# ----------- Setup Logging -----------
dir.create(outputLogs, recursive = TRUE, showWarnings = FALSE)
logFilePath <- file.path(outputLogs, "log_methylationGLMM_T1T2.txt")
logCon <- file(logFilePath, open = "wt")
sink(logCon, split = TRUE)
sink(logCon, type = "message")
# ==============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting DNAm LME Analysis (Timepoint 1 vs 3) ====\n")
cat("Start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input phenotype + beta file: ", inputPheno, "\n")
cat("Person ID variable: ", personVar, "\n")
cat("Timepoint variable: ", timeVar, "\n")
cat("Phenotypes: ", phenotypes, "\n")
cat("Covariates: ", covariates, "\n")
cat("PRS mapping: ", if (is.null(prsMap)) "None" else paste(unlist(prsMap), collapse = ", "), "\n")
cat("Factor variables: ", factorVars, "\n")
cat("LME libraries: ", lmeLibs, "\n")
cat("Library path: ", ifelse(is.null(libPath), "Default (R system)", libPath), "\n")
cat("CpG prefix pattern: ", cpgPrefix, "\n")
cat("CpG limit: ", ifelse(is.na(cpgLimit), "All CpGs", cpgLimit), "\n")
cat("Number of cores: ", nCores, "\n")
cat("Interaction term: ", interactionTerm, "\n")
cat("Save significant interactions: ", saveSignificantInteractions, "\n")
cat("Significant interaction p-value cutoff: ", significantInteractionPval, "\n")
cat("Directory to save significant interactions: ", significantInteractionDir, "\n")
cat("Save summary TXT files: ", saveTxtSummaries, "\n")
cat("Summary p-value threshold: ", summaryPval, "\n")
cat("Padjmethod: ", padjmethod, "\n")
cat("FDR threshold for significance: ", fdrThreshold, "\n")
cat("Multiple testing correction method: ", padjmethod, "\n")
cat("Save significant interaction CpGs: ", saveSignificantInteractions, "\n")
cat("Significant interaction p-value cutoff: ", significantInteractionPval, "\n")
cat("Directory to save significant interaction CpGs: ", significantInteractionDir, "\n")
cat("Summary TXT output directory: ", summaryTxtDir, "\n")
cat("Annotation package: ", annotationPackage, "\n")
cat("Annotation columns: ", annotationCols, "\n")
cat("FDR threshold for significance: ", fdrThreshold, "\n")
cat("Output logs directory: ", outputLogs, "\n")
cat("Directory to save annotated LME results: ", annotatedLMEOut, "\n")
cat("=======================================================================\n\n")

# ----------- Load Data -----------
cat("Loading input phenotype + beta file...\n")
inpt <- load(inputPheno)
phenoBT1T2 <- get(inpt)

# ----------- Check and Create person Column if Missing -----------
if (!(personVar %in% colnames(phenoBT1T2))) {
  cat(paste0("Column '",
             personVar, "' not found. Creating it from 'SID'...\n"))

  phenoBT1T2[[personVar]] <- as.numeric(factor(gsub("[AB]$", "",
                                                        phenoBT1T2$SID)))

  cat("Example mapping of SID to person ID:\n")
  print(head(phenoBT1T2[order(phenoBT1T2[[personVar]],
                              phenoBT1T2$SID),
                        c("SID", personVar)], 20
  ))

  cat("Count of records per person ID:\n")
  print(table(phenoBT1T2[[personVar]]))
} else {
  cat(paste0("Column '",
             personVar, "' already exists. Skipping creation.\n"))
}

cat("Checking structure of merged longitudinal dataset...\n")
print(table(table(phenoBT1T2[[personVar]])))
print(table(phenoBT1T2[[timeVar]]))
print(unique(phenoBT1T2[[timeVar]]))
print(dim(phenoBT1T2))
cat("=======================================================================\n")

# ----------- Summary Stats by Timepoint -----------
cat("Summary statistics for phenotype scores by timepoint:\n")
phenoBT1T2 %>%
        dplyr::select(all_of(timeVar), all_of(phenotypeList)) %>%
        group_by(.data[[timeVar]]) %>%
        summarise(across(everything(), list(
                mean = ~mean(., na.rm = TRUE),
                sd   = ~sd(., na.rm = TRUE),
                n    = ~sum(!is.na(.))
        ))) %>%
        print(width = Inf)
cat("=======================================================================\n")

# ----------- Run LME Function -----------
lme <- function(
                phenoScore,
                merge,
                personVar = personVar,
                timeVar = timeVar,
                covariates = covariateList,
                factorVars = factorVarsList,
                interactionTerm = interactionTerm,
                cpgPrefix = cpgPrefix,
                cpgLimit = cpgLimit,
                nCore = nCores,
                libPath = libPath,
                lmeLibs = lmeLibList
) {
        if (!phenoScore %in% colnames(merge)) {
                stop(paste("Phenotype", phenoScore, "not found"))
        }

        covariateNames <- c(timeVar, phenoScore, covariates)
        cpgCol <- grep(paste0("^", cpgPrefix), colnames(merge), value = TRUE)
        if (!is.na(cpgLimit)) {
                cpgCol <- head(cpgCol, as.numeric(cpgLimit))
        }

        cl <- makeCluster(nCore)
        clusterExport(cl, varlist = c("merge", "phenoScore", "personVar",
                                      "timeVar", "covariateNames", "factorVars",
                                      "interactionTerm", "libPath", "lmeLibs"),
                      envir = environment())

        # clusterEvalQ(cl, {
        #         if (!is.null(libPath)) {
        #                 .libPaths(libPath)
        #         }
        #         sapply(lmeLibs, function(pkg) {
        #                 if (!require(pkg, character.only = TRUE)) {
        #                         stop(paste("Failed to load package:", pkg))
        #                 }
        #         })
        # })

        clusterEvalQ(cl, {
                if (!is.null(libPath)) {
                        .libPaths(libPath)
                }

                for (pkg in lmeLibs) {
                        if (!requireNamespace(pkg, quietly = TRUE)) {
                        stop(paste("Package not installed:", pkg))
                        }
                }

                NULL
        })


        fit <- function(cpg) {
                tryCatch({
                        modelData <- merge[, c(personVar, covariateNames, cpg)]
                        colnames(modelData)[ncol(modelData)] <- "beta"
                        for (var in factorVars) {
                                if (var %in% colnames(modelData)) {
                                        modelData[[var]] <- as.factor(modelData[[var]])
                                }
                        }

                        if (!is.null(interactionTerm) && interactionTerm != "") {
                          intFormula <- paste(phenoScore, "*", interactionTerm)
                          fixedTerms <- setdiff(covariateNames, c(phenoScore,
                                                                  interactionTerm))
                          fixedPart <- paste(c(intFormula, fixedTerms), collapse = " + ")
                        } else {
                          fixedTerms <- setdiff(covariateNames, phenoScore)
                          fixedPart <- paste(c(phenoScore, fixedTerms), collapse = " + ")
                        }

                        form <- as.formula(paste("beta ~",
                                                 fixedPart,
                                                 "+ (1|",
                                                 personVar, ")"))

                        model <- lmer(form,
                                      data = modelData,
                                      na.action = na.exclude, REML = TRUE)

                        list(
                                coef = summary(model)$coefficients,
                                residuals = residuals(model),
                                fitted = fitted(model),
                                ranef = ranef(model),
                                fixef = fixef(model)
                        )
                }, error = function(e) NULL)
        }

        fitList <- parLapply(cl, cpgCol, fit)
        names(fitList) <- cpgCol
        stopCluster(cl)
        return(fitList)
}

# ----------- Run LMEs for All Phenotypes -----------
cat("Running LMEs for all phenotypes...\n")
for (pheno in phenotypeList) {
        outputFile <- file.path(outputRData, paste0(pheno, "LME.RData"))

        if (file.exists(outputFile)) {
                cat("Loading existing LME for:", pheno, "\n")
                load(outputFile)
                assign(paste0(pheno, "LME"), fitResult, envir = .GlobalEnv)
                next
        }

        cat("Running LME for:", pheno, "\n")

        prsVar <- if (pheno %in% names(prsMapList)) prsMapList[[pheno]] else NULL
        allCovariates <- if (!is.null(prsVar)) c(covariateList,
                                                 prsVar) else covariateList

        if (!is.null(interactionTerm) && interactionTerm != "") {
          fixedTerms <- setdiff(allCovariates, c(pheno, interactionTerm))
          modelFormula <- paste("~",
                                paste(c(paste0(pheno, "*", interactionTerm),
                                        fixedTerms),
                                      collapse = " + "),
                                "+ (1|", personVar, ")")
        } else {
          fixedTerms <- setdiff(allCovariates, pheno)
          modelFormula <- paste("~",
                                paste(c(pheno, fixedTerms),
                                      collapse = " + "),
                                "+ (1|", personVar, ")")
        }



        cat("Formula:", modelFormula, "\n")


        fitResult <- lme(
                phenoScore = pheno,
                merge = phenoBT1T2,
                personVar = personVar,
                timeVar = timeVar,
                covariates = allCovariates,
                factorVars = factorVarsList,
                interactionTerm = interactionTerm,
                cpgPrefix = cpgPrefix,
                cpgLimit = cpgLimit,
                nCore = nCores,
                libPath = libPath,
                lmeLibs = lmeLibList
        )

        save(fitResult, file = file.path(outputRData,
                                         paste0(pheno, "LME.RData")))

        assign(paste0(pheno, "LME"), fitResult, envir = .GlobalEnv)

}
cat("Finished running LMEs for all phenotypes.\n")
cat("=======================================================================\n")

# ----------- Extract LME Interaction Summary -----------
cpgsLME <- function(
                fitList,
                phenotype,
                interactionTerm = interactionTerm,
                pValue = summaryPval,
                nCore = nCores,
                libPath = libPath,
                lmeLibs = lmeLibList,
                chunkSize = chunkSize

) {
        cat("Extracting LME interaction effects for:", phenotype, "\n")

        if (is.null(interactionTerm) || interactionTerm == "") {
          cat("No interaction term provided, extracting main effects only.\n")
        } else {
          cat("Interaction term detected:", interactionTerm, "\n")
        }


        cpgNames <- names(fitList)
        if (is.null(chunkSize)) {
          chunkSize <- max(1000, floor(length(cpgNames) / (nCore * 4)))
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
                varlist = c("fitList", "pValue", "interactionTerm",
                            "phenotype", "libPath", "lmeLibs"),
                envir = environment()
        )

        # clusterEvalQ(cl, {
        #         if (!is.null(libPath)) {
        #                 .libPaths(libPath)
        #         }
        #         sapply(lmeLibs, function(pkg) {
        #                 if (!require(pkg, character.only = TRUE)) {
        #                         stop(paste("Failed to load package:", pkg))
        #                 }
        #         })
        # })

        clusterEvalQ(cl, {
                if (!is.null(libPath)) {
                        .libPaths(libPath)
                }

                for (pkg in lmeLibs) {
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
            modelOutput <- fitList[[cpg]]
            if (is.null(modelOutput) || is.null(modelOutput$coef)) next

            coefTable <- modelOutput$coef

            if (!is.null(interactionTerm) && interactionTerm != "") {
              pattern <- paste0("^", phenotype, ".*:", interactionTerm)
              matchedTerms <- grep(pattern, rownames(coefTable), value = TRUE)
            } else {
              pattern <- paste0("^", phenotype)
              matchedTerms <- grep(pattern, rownames(coefTable), value = TRUE)
            }

            if (length(matchedTerms) == 0) next

            tmp <- tryCatch({
              do.call(rbind, lapply(matchedTerms, function(term) {
                coefRow <- coefTable[term, ]
                data.frame(
                  CpG = cpg,
                  Interaction.Term = term,
                  Estimate = coefRow["Estimate"],
                  Std.Error = coefRow["Std. Error"],
                  t.value = coefRow["t value"],
                  P.value = coefRow["Pr(>|t|)"],
                  row.names = NULL
                )
              }))})
            if (!is.null(tmp)) {
              outList[[idx]] <- tmp
              idx <- idx + 1
              }
          }
          if (idx > 1) {
            return(data.table::rbindlist(outList[1:(idx-1)],
                                         use.names = TRUE, fill = TRUE))
          } else {
            return(NULL)
          }
        })

        stopCluster(cl)

        summary <- do.call(rbind, results)

        if (!is.na(pValue)) {
                summary <- subset(summary, `P.value` < pValue)
        }


        cat("Finished extracting interaction for:", phenotype, "\n")
        return(summary)
}

# ----------- Run and Save LME Summaries -----------
cat("Running LME summary extraction for all phenotypes...\n")

for (pheno in phenotypeList) {
        summaryFile <- file.path(outputRData, paste0(pheno,
                                                         "SummaryLME.RData"))

        if (file.exists(summaryFile)) {
                cat("Loading existing summary for:", pheno, "\n")
                load(summaryFile)
                assign(paste0(pheno, "SummaryLME"), fitResult,
                       envir = .GlobalEnv)
                next
        }

        cat("Processing interaction effects for:", pheno, "\n")

        fitObjectName <- paste0(pheno, "LME")
        if (!exists(fitObjectName)) {
                cat("Object not found:", fitObjectName, "- skipping\n")
                next
        }

        fitObject <- get(fitObjectName)

        fitResult <- cpgsLME(
                fitList = fitObject,
                phenotype = pheno,
                nCore = nCores,
                pValue = summaryPval,
                libPath = libPath,
                lmeLibs = lmeLibList,
                interactionTerm = interactionTerm
        )

        save(fitResult, file = summaryFile)

        assign(paste0(pheno, "SummaryLME"), fitResult, envir = .GlobalEnv)

        cat("Saved:", paste0(pheno, "SummaryLME.RData"), "\n")
}

cat("Finished extracting all LME summaries.\n")
cat("=======================================================================\n")

# ----------- Save Significant LME Interaction Results -----------
saveSignificantInteractions <- function(
                resultList,
                resultName = deparse(substitute(resultList)),
                baseDir = significantInteractionDir,
                pvalThreshold = significantInteractionPval,
                interactionTerm = interactionTerm
) {
        resultDir <- file.path(baseDir, resultName)
        if (!dir.exists(resultDir)) dir.create(resultDir, recursive = TRUE)

        if (is.null(interactionTerm) || interactionTerm == "") {
          cat("No interaction term detected, extracting main effects for",
              resultName, "\n")
          pattern <- paste0("^", resultName)
        } else {
          cat("Interaction term detected:", interactionTerm,
              "- extracting interaction effects for", resultName, "\n")
          pattern <- paste0("^", resultName, ".*:", interactionTerm)
        }

        for (i in seq_along(resultList)) {
          coefTable <- resultList[[i]]$coef
          cpgName <- names(resultList)[i]
          matchedRows <- grep(pattern, rownames(coefTable), value = FALSE)

          if (length(matchedRows) > 0) {
            termPvals <- coefTable[matchedRows, "Pr(>|t|)"]
            if (any(termPvals  < pvalThreshold, na.rm = TRUE)) {
              cpgDir <- file.path(resultDir, cpgName)
              if (!dir.exists(cpgDir)) dir.create(cpgDir)

              outputFile <- file.path(cpgDir, paste0(cpgName, ".txt"))
              write.table(coefTable, file = outputFile,
                          sep = "\t", quote = FALSE)
            }
          }
        }
}

cat("Saving significant LME interaction results...\n")

cat("Saving significant interaction CpGs to:", significantInteractionDir, "\n")
for (pheno in phenotypeList) {
        modelObjName <- paste0(pheno, "LME")
        if (exists(modelObjName)) {
                resultObj <- get(modelObjName)
                saveSignificantInteractions(resultList = resultObj, resultName = pheno)
        }
}
cat("Finished saving significant CpG interaction outputs.\n")
cat("=======================================================================\n")

# ----------- Save Summary Tables as TXT -----------
saveSummaryToTxt <- function(
                summaryDF,
                outputFile
) {
        if ("P.value" %in% colnames(summaryDF)) {
                summaryDF <- summaryDF[order(summaryDF$P.value), ]
        } else if ("Pr(>|t|)" %in% colnames(summaryDF)) {
                summaryDF <- summaryDF[order(summaryDF[["Pr(>|t|)"]]), ]
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

# ----------- Run and Save TXT Summaries -----------

if (saveTxtSummaries) {
        cat("Saving summaries to TXT files...\n")

        phenoList <- strsplit(phenotypes, ",")[[1]]

        for (pheno in phenoList) {

                outputFile <- file.path(
                        summaryTxtDir,
                        paste0(pheno, "SummaryLME.txt")
                )

                summaryObjName <- paste0(pheno, "SummaryLME")

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

# ----------- Diagnostic Plots for LME Results -----------
diagnosticPlotsLME <- function(summary,
                               phenoBeta,
                               variable,
                               fdrThreshold = fdrThreshold,
                               padjmethod = padjmethod,
                               cpgPrefix = cpgPrefix,
                               outputDir = outputPlots,
                               plotWidth = plotWidth,
                               plotHeight = plotHeight,
                               plotDPI = plotDPI
                               ) {

        ## Multiple Testing Correction
        summary$FDR <- p.adjust(summary$P.value, method = padjmethod)

        ## Genomic Inflation Factor (λ)
        chisq <- qchisq(1 - summary$P.value, df = 1)
        lambda <- round(median(chisq, na.rm = TRUE) / qchisq(0.5, df = 1), 3)
        message("Genomic inflation factor: ", lambda)

        ## Q-Q Plot
        pvals <- summary$P.value
        pvals <- pvals[!is.na(pvals)]

        tiff(filename = file.path(outputDir,
                                  paste0("qqplot_", variable, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")

        qqplot(-log10(ppoints(length(pvals))), -log10(sort(pvals)),
               main = paste("Q-Q Plot of p-values for", variable,
                            "\nGenomic Inflation Factor= ", lambda),
               xlab = "Expected -log10(p)",
               ylab = "Observed -log10(p)",
               pch = 16, col = "black")
        abline(0, 1, col = "red")
        dev.off()

        ## Mean Beta Values
        cpgCols <- grep(paste0("^", cpgPrefix),
                        colnames(phenoBeta), value = TRUE)
        betas <- phenoBeta[, cpgCols]
        meanBeta <- colMeans(betas, na.rm = TRUE)
        meanBetaLog <- log2(meanBeta)

        # Merge log2meanBeta with summary
        summary$log2meanBeta <- meanBetaLog[summary$CpG]
        summary <- summary[!is.na(summary$log2meanBeta), ]

        ## Plot: Residual Proxy (SD estimate) vs Mean Methylation
        tiff(filename = file.path(outputDir,
                                  paste0("residualSD_", variable, ".tiff")),
             width = plotWidth,
             height = plotHeight,
             res = plotDPI, type = "cairo")

        plot(summary$log2meanBeta,
             summary$Std.Error,
             pch = 20,
             col = rgb(0, 0, 0, 0.4),
             xlab = "log2(Average Beta)",
             ylab = "Standard Error",
             main = "SD vs Average Beta (proxy from Std. Error)")
        lines(lowess(summary$log2meanBeta, summary$Std.Error),
              col = "red", lwd = 2)
        dev.off()

        ## Plot: Significance vs Variability (colored by FDR)
        p <- ggplot(summary, aes(x = -log10(P.value),
                            y = Std.Error, color = FDR < fdrThreshold)) +
                geom_point(alpha = 0.6) +
                geom_text_repel(data = subset(summary,
                                              FDR < fdrThreshold),
                                aes(label = CpG)) +
                scale_color_manual(values = c("FALSE" = "grey50",
                                              "TRUE" = "firebrick")) +
                labs(title = paste("Standard Error vs Significance for", variable),
                     x = "-log10(p-value)",
                     y = "Standard Error",
                     color = paste("FDR <", fdrThreshold)) +
                theme_minimal()

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
# ==============================================================================

# ----------- Generate Diagnostic Plots for All Phenotypes -----------
cat("Generating diagnostic plots for all phenotypes...\n")

# Parse phenotype list
phenotypeNames <- strsplit(phenotypes, ",")[[1]]

for (pheno in phenotypeNames) {

        summaryName <- paste0(pheno, "SummaryLME")

        if (exists(summaryName)) {
                cat("Generating diagnostic plots for:", pheno, "\n")

                diagnosticPlotsLME(
                        summary = get(summaryName),
                        phenoBeta = phenoBT1T2,
                        variable = pheno,
                        cpgPrefix = cpgPrefix,
                        fdrThreshold = fdrThreshold,
                        padjmethod = padjmethod,
                        outputDir = outputPlots,
                        plotWidth = plotWidth,
                        plotHeight = plotHeight,
                        plotDPI = plotDPI
                )
        } else {
                cat("Warning: Summary object not found for phenotype:", pheno, "\n")
        }
}

cat("Finished generating all diagnostic plots.\n")
cat("Plots saved to:", outputPlots, "\n")
cat("=======================================================================\n")

# ----------- Load Annotation Object -----------
cat("Loading annotation object:", annotationPackage, "\n")
annotationObject <- getAnnotation(get(annotationPackage))

cat("Annotation loaded with", nrow(annotationObject), "probes\n")
cat("Available columns:\n")
print(colnames(annotationObject))
cat("=======================================================================\n")

# ----------- Function: Annotate LME Summaries -----------
annotateLME <- function(
                summaryList,
                annotationObject,
                annotationCols = strsplit(annotationCols, ",")[[1]]
) {
        modelNames <- names(summaryList)

        cat("Merging LME summaries...\n")

        cleanedSummaries <- list()

        for (modelName in modelNames) {
                df <- summaryList[[modelName]]

                if (!all(c("CpG", "Interaction.Term", "P.value") %in% colnames(df))) {
                        warning(paste("Skipping", modelName, ": required columns not found"))
                        next
                }

                dfSplit <- split(df, df$Interaction.Term)

                modelTables <- lapply(names(dfSplit), function(term) {
                        subDf <- dfSplit[[term]]

                        interactionSuffix <- gsub(paste0("^", modelName, "\\."), "", term)
                        pCol <- paste0(modelName, "_", interactionSuffix, "_P.Value")

                        subDf[[pCol]] <- subDf$P.value
			          subDf <- as.data.frame(subDf)
                        subDf <- subDf[, c("CpG", pCol)]
                        return(subDf)
                })

                mergedModel <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE),
                                      modelTables)
                cleanedSummaries[[modelName]] <- mergedModel
        }

        mergedSummary <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE),
                                cleanedSummaries)

        annDF <- as.data.frame(annotationObject)

	annDF$CpG <- rownames(annDF)

        if (!is.null(annotationCols)) {
                annDF <- annDF[, unique(c("CpG", annotationCols)), drop = FALSE]
        }

        annotatedResults <- merge(mergedSummary, annDF, by = "CpG", all.x = TRUE)
        colnames(annotatedResults)[1] <- "IlmnID"

        cat("Annotation completed. Annotated CpGs:", nrow(annotatedResults), "\n")
        return(annotatedResults)
}

# ----------- Execute Annotation and Save CSV -----------
cat("Running annotation of LME summary results...\n")

# Build summaryList dynamically
phenotypeNames <- strsplit(phenotypes, ",")[[1]]
summaryList <- setNames(
        lapply(phenotypeNames, function(pheno) {
                get(paste0(pheno, "SummaryLME"))
        }),
        phenotypeNames
)

cat("\n Phenotypes to annotate:\n")
print(names(summaryList))

for (nm in names(summaryList)) {
  cat("First few rows of", nm, "summary:\n")
  print(head(summaryList[[nm]]))
  cat("Column names in", nm, "summary:\n")
  print(colnames(summaryList[[nm]]))
}

# Annotate
annotatedLME <- annotateLME(
        summaryList = summaryList,
        annotationObject = annotationObject,
        annotationCols = strsplit(annotationCols, ",")[[1]]
)

# Save as CSV inside the output directory
annotatedLMEPath <- file.path(annotatedLMEOut, "annotatedLME.csv")

cat("Saving annotated LME results to:", annotatedLMEPath, "\n")

write.csv(
        annotatedLME,
        file = annotatedLMEPath,
        row.names = FALSE
)

cat("Finished writing annotated LME results.\n")
cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")
sink()
close(logCon)
}