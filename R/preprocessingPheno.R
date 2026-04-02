#' Run preprocessingPheno.R
#' @importFrom magrittr %>%
#'
#' @param phenoFile Character. Path to phenotype file.
#' @param sepType Character. Field separator for phenotype file.
#' @param betaPath Character. Path to beta-values RData file.
#' @param mPath Character. Path to M-values RData file.
#' @param cnPath Character. Path to copy number RData file.
#' @param SampleID Character. Sample identifier column name.
#' @param timeVar Character. Time variable column name.
#' @param timepoints Character. Timepoints to retain.
#' @param combineTimepoints Character. Timepoints to combine.
#' @param outputPheno Character. Output directory for processed phenotype files.
#' @param outputRData Character. Output directory for metric RData files.
#' @param outputRDataMerge Character. Output directory for merged RData files.
#' @param sexColumn Character. Sex column name.
#' @param outputLogs Character. Directory for log files.
#' @param outputDir Character. Base output directory.
#'
#' @return
#' Invisibly returns \code{NULL}. This function is called for its side effects,
#' preparing analysis-ready phenotype-methylation datasets by subsetting,
#' merging timepoints, aligning samples with beta, M, and CN matrices, and
#' writing processed phenotype tables and RData objects to disk.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' preprocessingPheno(
#'   phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
#'   sepType = "",
#'   betaPath = "rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
#'   mPath = "rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
#'   cnPath = "rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
#'   SampleID = "Sample_Name",
#'   timeVar = "Timepoint",
#'   timepoints = "1,2",
#'   combineTimepoints = "1,2",
#'   outputPheno = "data/preprocessingPheno",
#'   outputRData = "rData/preprocessingPheno/metrics",
#'   outputRDataMerge = "rData/preprocessingPheno/mergeData",
#'   sexColumn = "Sex",
#'   outputLogs = "logs",
#'   outputDir = "data/preprocessingPheno"
#' )
#' }
#'
#' @export
preprocessingPheno <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
    sepType = "",
    betaPath =
      "rData/preprocessingMinfiEwasWater/metrics/beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    mPath =
      "rData/preprocessingMinfiEwasWater/metrics/m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    cnPath =
      "rData/preprocessingMinfiEwasWater/metrics/cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData",
    SampleID = "Sample_Name",
    timeVar = "Timepoint",
    timepoints = "1,2",
    combineTimepoints = "1,2",
    outputPheno = "data/preprocessingPheno",
    outputRData = "rData/preprocessingPheno/metrics",
    outputRDataMerge = "rData/preprocessingPheno/mergeData",
    sexColumn = "Sex",
    outputLogs = "logs",
    outputDir = "data/preprocessingPheno"
) {

dir.create(outputRData, recursive = TRUE, showWarnings = FALSE)
dir.create(outputRDataMerge, recursive = TRUE, showWarnings = FALSE)
dir.create(outputLogs, recursive = TRUE, showWarnings = FALSE)
dir.create(outputPheno, recursive = TRUE, showWarnings = FALSE)
dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

#===============================================================================

# ----------- Logging Setup -----------
logFilePath <- file.path(outputLogs, "log_preprocessingPheno.txt")
logCon <- file(logFilePath, open = "wt")

sink(logCon, split = TRUE)
sink(logCon, type = "message")
#===============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting Phenotype Preprocessing ====\n")
cat("Start Time:               ", format(Sys.time()), "\n")
cat("Log file path:            ", logFilePath, "\n\n")
cat("Phenotype file:           ", phenoFile, "\n")
cat("Beta path:                ", betaPath, "\n")
cat("M-values path:            ", mPath, "\n")
cat("CN path:                  ", cnPath, "\n\n")

cat("Identifier column:        ", SampleID, "\n")
cat("Timepoint column:        ", timeVar, "\n")
cat("Timepoints (if present):  ", timepoints, "\n")
cat("Combine timepoints:       ", combineTimepoints, "\n\n")
cat("Sex column:               ", sexColumn, "\n")

cat("Output phenotype dir:     ", outputPheno, "\n")
cat("RData metrics dir:        ", outputRData, "\n")
cat("RData merge dir:          ", outputRDataMerge, "\n\n")
cat("=======================================================================\n")

# ----------- Load Data -----------
load(betaPath)
load(mPath)
load(cnPath)

cat("Beta dimensions: ", dim(beta), "\n")
cat("M dimensions: ", dim(m), "\n")
cat("CN dimensions: ", dim(cn), "\n")

# ----------- Read Phenotype File -----------
if (sepType == "\\t") {
  sepChar <- "\t"
} else if (sepType == "") {
  sepChar <- NULL
} else {
  sepChar <- sepType
}

# Now read the phenotype file
if (!is.null(sepChar)) {
  pheno <- utils::read.csv(phenoFile, sep = sepChar)
} else {
  pheno <- utils::read.csv(phenoFile)
}

cat("Phenotype file loaded with",
    nrow(pheno), "samples and", ncol(pheno), "columns.\n")
cat("Preview of phenoLC:\n")
print(utils::head(pheno[, 1:5]))
cat("=======================================================================\n")

# ----------- Subsetting Timepoints & Data Splitting -----------
timepoints <- as.numeric(strsplit(timepoints, ",")[[1]])
cat("Subsetting to timepoints:", paste(timepoints, collapse = ", "), "\n")

# Print available timepoints in the phenotype
cat("Available values in", timeVar, "column:\n")
print(table(pheno[[timeVar]], useNA = "ifany"))

for (tp in timepoints) {
  # Subset phenotype by Timepoint
  phenoSub <- subset(pheno, pheno[[timeVar]] == tp)
  assign(paste0("phenoT", tp), phenoSub)

  # Subset matrices using SID (SampleID) from phenoSub
  sids <- as.character(phenoSub[[SampleID]])

  assign(paste0("betaT", tp), beta[, sids])
  assign(paste0("mT", tp),    m[,    sids])
  assign(paste0("cnT", tp),   cn[,   sids])
}

# Save each subset
for (tp in timepoints) {
        utils::write.csv(get(paste0("phenoT", tp)), file = file.path(outputPheno,
                                                              paste0("phenoT",
                                                                     tp, ".csv")),
                  row.names = FALSE)
        save(list = paste0("betaT", tp), file = file.path(outputRData,
                                                          paste0("betaT",
                                                                 tp, ".RData")))
        save(list = paste0("mT", tp), file = file.path(outputRData,
                                                       paste0("mT",
                                                              tp, ".RData")))
}

# ----------- Merge Combined Timepoints for Longitudinal Analysis -----------
cat("Combining timepoints:", combineTimepoints, "\n")
combineTPs <- as.numeric(strsplit(combineTimepoints, ",")[[1]])

combinedPhenoList <- lapply(combineTPs, function(tp) get(paste0("phenoT", tp)))
phenoCombined <- do.call(rbind, combinedPhenoList)

combineSuffix <- paste0("T", paste(combineTPs, collapse = "T"))
utils::write.csv(phenoCombined,
          file = file.path(outputPheno,
                                          paste0("pheno",
                                                 combineSuffix, ".csv")),
          row.names = FALSE)

cat("Saved combined phenotype file for T1T2 at:", outputPheno, "\n")
cat("=======================================================================\n")

# ----------- Merge Beta Matrix with Phenotype ----------
mergeBeta <- function(phenoFrame, betaMatrix, id = SampleID) {
        rownames(phenoFrame) <- phenoFrame[[id]]
        matched <- intersect(rownames(phenoFrame), colnames(betaMatrix))
        phenoFrame <- phenoFrame[matched, ]
        betaMatrix <- betaMatrix[, matched]

        betaTranp <- as.data.frame(t(betaMatrix))
        mergedData <- cbind(phenoFrame, betaTranp)
        return(mergedData)
}

# Perform merge for each timepoint

mergedList <- list()
for (tp in timepoints) {
        cat("Processing merge for timepoint:", tp, "\n")

        phenoObj <- paste0("phenoT", tp)
        betaObj <- paste0("betaT", tp)

        if (!exists(phenoObj) || !exists(betaObj)) {
                cat("Warning: One or both objects not found for T", tp, "\n", sep = "")
                next
        }

        phenoTemp <- get(phenoObj)
        betaTemp <- get(betaObj)

        cat("  - pheno rows:", nrow(phenoTemp), "\n")
        cat("  - beta cols:", ncol(betaTemp), "\n")

        mergedTemp <- tryCatch({
                mergeBeta(phenoTemp, betaTemp)
        }, error = function(e) {
                cat("[ERROR] mergeBeta failed for timepoint", tp, ":\n", conditionMessage(e), "\n")
                return(NULL)
        })

        if (!is.null(mergedTemp)) {
                mergedList[[as.character(tp)]] <- mergedTemp
                save(mergedTemp, file = file.path(outputRDataMerge,
                                                  paste0("phenoBetaT", tp, ".RData")))
                cat("Saved merged object for T", tp, "\n", sep = "")
        } else {
                cat("Skipping save for T", tp, " due to error\n")
        }
}


# ----------- Combine merged phenotype + beta matrix -----------
combined <- do.call(rbind, mergedList[as.character(combineTPs)])
save(combined,
     file = file.path(outputRDataMerge,
                      paste0("phenoBeta", combineSuffix, ".RData")))

cat("Combined data saved for timepoints:",
    paste(combineTPs, collapse = ", "), "\n")

cat("Merged data saved to:", outputRDataMerge, "\n")
cat("=======================================================================\n")
# ----------- Preprocessing Betas for Horvath Calculator -----------

betaCSV <- as.data.frame(beta)
betaCSV <- tibble::rownames_to_column(betaCSV, var = "ProbeID")

# Inspect changes
dim(betaCSV)
print(utils::head(betaCSV)[1:5, 1:5])

betaCSVPath <- file.path(outputDir, "beta.csv")
utils::write.csv(betaCSV, file = betaCSVPath, row.names = FALSE)
cat("Beta CSV file for ClockFundation saved to:", betaCSVPath, "\n")

zipFile <- file.path(outputDir, "beta.zip")

utils::zip(zipfile = zipFile, files = betaCSVPath, flags = "-j")

cat("Beta ZIP file for ClockFundation saved to:", zipFile, "\n")

# ----------- Preprocessing CSV for Horvath Calculator -----------

# Rename the column "SampleName" to "id"
pheno <- pheno %>%
  dplyr::rename(id = SampleID)

# Recode Sex only if values are not already "Male"/"Female"
uniqueSex <- unique(pheno[[sexColumn]])

if (!all(uniqueSex %in% c("Male", "Female"))) {
  cat("Re-encoding Sex: 0 = Female, 1 = Male\n")
  pheno[[sexColumn]] <- ifelse(pheno[[sexColumn]] == 0, "Female", "Male")
} else {
  cat("Sex column already contains 'Male' and 'Female'. Skipping recoding.\n")
}

phenoCSVPath <- file.path(outputDir, "phenoCF.csv")

utils::write.csv(pheno, file = phenoCSVPath, row.names = FALSE)
cat("Sample file for ClockFundation saved to:", phenoCSVPath, "\n")

cat("=======================================================================\n")

# ----------- Close Logging -----------
cat("\nSession Info:\n")
print(utils::sessionInfo())
# =============================================================================
sink(type = "message")
sink()
close(logCon)
}
