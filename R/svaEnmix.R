#' Run svaEnmix.R
#' @import minfi
#' @import Gviz
#' @import ENmix
#' @import ggplot2
#' @import ggpubr
#' @importFrom MASS dropterm
#'
#' @param phenoFile Character. Path to phenotype file with cell composition data.
#' @param rgsetData Character. Path to RGSet RData file.
#' @param sepType Character. Field separator for phenotype file.
#' @param outputLogs Character. Directory for log files.
#' @param nSamples Integer or NA. Number of samples to subset for testing.
#' @param SampleID Character. Sample identifier column name.
#' @param arrayType Character. Illumina array type.
#' @param annotationVersion Character. Annotation version.
#' @param SentrixIDColumn Character. Sentrix ID column name.
#' @param SentrixPositionColumn Character. Sentrix position column name.
#' @param ctrlSvaPercVar Numeric. Proportion of variance explained by control probes.
#' @param ctrlSvaFlag Integer. Flag indicating use of control probes.
#' @param scriptLabel Character. Label used in output naming.
#' @param tiffWidth Integer. Width of TIFF plots in pixels.
#' @param tiffHeight Integer. Height of TIFF plots in pixels.
#' @param tiffRes Integer. Resolution (DPI) for TIFF plots.
#' @param figureBaseDir Character. Base directory for Figures outputs.
#' @param dataBaseDir Character. Base directory for Data outputs.
#'
#' @return
#' Invisibly returns \code{NULL}. This function is called for its side effects,
#' executing the external \code{svaEnmix.R} script and writing results, figures,
#' and logs to disk.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' svaEnmix(
#'   phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
#'   rgsetData = "rData/preprocessingMinfiEwasWater/objects/RGSet.RData",
#'   sepType = "",
#'   outputLogs = "logs",
#'   nSamples = 5,
#'   SampleID = "Sample_Name",
#'   arrayType = "IlluminaHumanMethylationEPICv2",
#'   annotationVersion = "20a1.hg38",
#'   SentrixIDColumn = "Sentrix_ID",
#'   SentrixPositionColumn = "Sentrix_Position",
#'   ctrlSvaPercVar = 0.90,
#'   ctrlSvaFlag = 1,
#'   scriptLabel = "svaEnmix",
#'   tiffWidth = 2000,
#'   tiffHeight = 1000,
#'   tiffRes = 150,
#'   figureBaseDir = "figures",
#'   dataBaseDir = "data"
#' )
#' }
#'
#' @export
svaEnmix <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/phenoLC.csv",
    rgsetData = "rData/preprocessingMinfiEwasWater/objects/RGSet.RData",
    sepType = "",
    outputLogs = "logs",
    nSamples = NA,
    SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
    annotationVersion = "20a1.hg38",
    SentrixIDColumn = "Sentrix_ID",
    SentrixPositionColumn = "Sentrix_Position",
    ctrlSvaPercVar = 0.90,
    ctrlSvaFlag = 1,
    scriptLabel = "svaEnmix",
    tiffWidth = 2000,
    tiffHeight = 1000,
    tiffRes = 150,
    figureBaseDir = "figures",
    dataBaseDir = "data",
    rBaseDir = "rData"
) {

# ----------- Logging Setup -----------
dir.create(outputLogs, recursive = TRUE, showWarnings = FALSE)

logFilePath <- file.path(outputLogs,"log_svaEnmix.txt")
logCon <- file(logFilePath, open = "wt")

sink(logCon, split = TRUE)
sink(logCon, type = "message")
# ==============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting SVA Estimation with Enmix ====\n")
cat("Start time: ", format(Sys.time()), "\n\n")
cat("Log file path: ", logFilePath, "\n\n")
cat("Pheno file: ", phenoFile, "\n")
cat("Separator type: ", ifelse(is.null(sepType), "default, No separator", sepType), "\n")
cat("Log directory: ", outputLogs, "\n")
cat("Sample limit: ", ifelse(is.na(nSamples), "All", nSamples), "\n")
cat("SampleID column: ", SampleID, "\n")
cat("Sentrix ID column: ", SentrixIDColumn, "\n")
cat("Sentrix Position column: ", SentrixPositionColumn, "\n")
cat("Script label: ", scriptLabel, "\n")
cat("ctrlSva percvar: ", ctrlSvaPercVar, "\n")
cat("ctrlSva flag: ", ctrlSvaFlag, "\n")
cat("TIFF dimensions (WxH): ", tiffWidth, "x", tiffHeight, " at", tiffRes, "dpi\n")
# =============================================================================

# ----------- Directory Setup for Figures and Data -----------
dir.create(file.path(figureBaseDir, scriptLabel), showWarnings = FALSE,
           recursive = TRUE)
dir.create(file.path(dataBaseDir, scriptLabel), showWarnings = FALSE,
           recursive = TRUE)
dir.create(file.path(rBaseDir, scriptLabel), showWarnings = FALSE,
           recursive = TRUE)
cat("=======================================================================\n")

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
  targets <- read.csv(phenoFile, sep = sepChar)
} else {
  targets <- read.csv(phenoFile)
}

if (!is.na(nSamples) && nSamples < nrow(targets)) {
  targets <- targets[1:nSamples, ]
  cat("Subsetting to", nSamples, "samples for testing.\n")
} else {
  cat("Using all", nrow(targets), "samples.\n")
}

cat("Phenotype file loaded with",
    nrow(targets), "samples and", ncol(targets), "columns.\n")
cat("Preview of targets:\n")
print(head(targets[, 1:6]))
cat("=======================================================================\n")

# ----------- Load IDAT Files into RGSet -----------
load(rgsetData)

# Assign custom sample names
sampleNames(RGSet) <- targets[[SampleID]]
cat("RGSet loaded with", length(sampleNames(RGSet)), "samples.\n")
cat("=======================================================================\n")

# ----------- Estimate Surrogate Variables from Control Probes -----------
sva <- ctrlsva(
  rgSet = RGSet,
  percvar = ctrlSvaPercVar,
  flag = ctrlSvaFlag
)
cat("Surrogate variables matrix (first few rows):\n")
print(head(sva))

# ---- Save SVA matrix ----
svaSentrixRDataPath <- file.path(rBaseDir,
                                 scriptLabel, "svaMatrix.RData")
save(sva, file = svaSentrixRDataPath)
cat("SVA Matrix RData saved to: ", svaSentrixRDataPath, "\n")

svaSentrixRDataCSV <- file.path(dataBaseDir,
                                scriptLabel, "svaMatrix.csv")
write.csv(sva, svaSentrixRDataCSV, row.names = TRUE)
cat("SVA Matrix CSV saved to: ", svaSentrixRDataCSV, "\n")

# ---- Prepare SVA for merge ----
svaD <- data.frame(SID = rownames(sva), sva, row.names = NULL)
names(svaD)[1] <- SampleID

# ---- Merge with phenotype ----
pheno <- merge(targets, svaD, by = SampleID, all.x = TRUE)

write.csv(pheno,
          file = phenoFile,
          row.names = FALSE)
cat("Saved phenoLC + SVA:", phenoFile, "\n")

# ----------- Plot SVA Colored by SentrixID -----------
sentrixID <- as.factor(pData(RGSet)[[SentrixIDColumn]])

# Create TIFF output
svaSentrixPath <- file.path(figureBaseDir,
                        scriptLabel, "sva_SentrixID.tiff")
tiff(filename = svaSentrixPath,
     width = tiffWidth,
     height = tiffHeight,
     res = tiffRes, type = "cairo")

plot(sva[, 1], sva[, 2],
     col = rainbow(length(levels(sentrixID)))[sentrixID],
     pch = 16,
     xlab = "Surrogate Variable 1 (PC1)",
     ylab = "Surrogate Variable 2 (PC2)",
     main = "Surrogate Variables Colored by Chip (SentrixID)")
legend("topright", legend = levels(sentrixID),
       col = rainbow(length(levels(sentrixID))),
       pch = 16, title = "SentrixID", cex = 0.6)

dev.off()

cat("SVA Sentrix plot saved to: ", svaSentrixPath, "\n")

# ----------- Plot SVA Colored by SentrixPosition -----------
sentrixPos <- as.factor(pData(RGSet)[[SentrixPositionColumn]])

svaPositionpath <- file.path(figureBaseDir, scriptLabel, "sva_SentrixPosition.tiff")

tiff(filename = svaPositionpath,
     width = tiffWidth,
     height = tiffHeight,
     res = tiffRes, type = "cairo")

plot(sva[, 1], sva[, 2],
     col = rainbow(length(levels(sentrixPos)))[sentrixPos],
     pch = 16,
     xlab = "Surrogate Variable 1 (PC1)",
     ylab = "Surrogate Variable 2 (PC2)",
     main = "Surrogate Variables Colored by Sentrix Position")
legend("topright", legend = levels(sentrixPos),
       col = rainbow(length(levels(sentrixPos))),
       pch = 16, title = "SentrixPosition", cex = 0.6)

dev.off()

cat("SVA Position plot saved to: ", svaPositionpath, "\n")
cat("=======================================================================\n")

# ----------- Linear Models for Surrogate Variables (ANOVA) -----------
K <- ncol(sva)
cat("Number of surrogate variables (K):", K, "\n")

# Print class and levels of SentrixID
cat("SentrixID class:", class(sentrixID), "\n")
cat("SentrixID unique levels:", length(unique(sentrixID)), "\n")
print(table(sentrixID))

# Print class and levels of SentrixPosition
cat("SentrixPosition class:", class(sentrixPos), "\n")
cat("SentrixPosition unique levels:", length(unique(sentrixPos)), "\n")
print(table(sentrixPos))

# Print example row of SVA matrix
cat("First row of SVA matrix:\n")
print(sva[1, ])

# Confirm if sample names align
cat("Sample names in SVA matrix:",
    paste(rownames(sva)[1:5], collapse = ", "), "\n")
cat("Sample names in pData(RGSet):",
    paste(rownames(pData(RGSet))[1:5], collapse = ", "), "\n")

# Fit linear models for each surrogate variable
lmsvaFull <- lapply(1:K, function(i)
  lm(sva[, i] ~ SentrixID + SentrixPosition,
     data.frame("SentrixID" = sentrixID,
                "SentrixPosition" = sentrixPos))
)

lmsvaRed <- vector("list", K)

# ----------- Save summaries of full models -----------
capture.output(summary(lmsvaFull[[1]]),
               file = file.path(dataBaseDir, scriptLabel, "summary_full_sva1.txt"))

if (K >= 2) {
  capture.output(summary(lmsvaFull[[2]]),
                 file = file.path(dataBaseDir, scriptLabel, "summary_full_sva2.txt"))
}

# Perform backward elimination and write ANOVA output
for(i in 1:K){
  lmtmp = lmsvaFull[[i]]
  while(1){
    dttmp = dropterm(lmtmp, test = "F")
    if(max(dttmp$`Pr(F)`, na.rm = TRUE) > (0.05))
      ttmp = rownames(dttmp)[which.max(dttmp$`Pr(F)`)]
    else break
    lmtmp = update(lmtmp, paste(".~. - ", ttmp) )
    capture.output(dttmp,
                   file = file.path(dataBaseDir,
                                    scriptLabel, paste0("dropterm_step_sva", i, ".txt")),
                   append = TRUE)
    capture.output(summary(lmtmp),
                   file = file.path(dataBaseDir,
                                    scriptLabel, paste0("dropterm_model_sva", i, ".txt")),
                   append = TRUE)
  }

  lmsvaRed[[i]] = lmtmp
}

# ----------- Save ANOVA summaries for full and reduced models -----------
for (i in 1:K) {
  capture.output(anova(lmsvaFull[[i]]),
                 file = file.path(dataBaseDir,
                                  scriptLabel, paste0("anova_full_sva", i, ".txt")))

  capture.output(anova(lmsvaRed[[i]]),
                 file = file.path(dataBaseDir,
                                  scriptLabel, paste0("anova_reduced_sva", i, ".txt")))
}
cat("=======================================================================\n")

# ----------- Plot Matrix of Surrogate Variables Colored by SentrixID and Shape by Position -----------

# Prepare output TIFF file

svaSentrixPositionPath <- file.path(figureBaseDir,
                                    scriptLabel, "sva_SentrixIDPosition.tiff")

tiff(filename = svaSentrixPositionPath,
     width = tiffWidth,
     height = tiffHeight,
     res = tiffRes,
     type = "cairo")

# Prepare plotting layout
par(mfrow = c(K, K), family = "Times", las = 1)

# Extract and map IDs/positions
colorMap <- rainbow(length(levels(sentrixID)))
pchMap <- 1:length(levels(sentrixPos))

# Plot matrix
for (i in 1:K) {
  for (j in 1:K) {
    plot(sva[, j], sva[, i],
         col = colorMap[sentrixID],
         pch = pchMap[sentrixPos],
         xlab = paste("SV", j),
         ylab = paste("SV", i),
         main = "Effects of Sentrix ID (color) & Sentrix Position (shape)")

    # Legend only in top-left panel
    if (i == 1 && j == 1) {
      legend("topright",
             legend = levels(sentrixID),
             col = colorMap,
             pch = 15,
             title = "SentrixID",
             cex = 0.6)
      legend("bottomright",
             legend = levels(sentrixPos),
             pch = pchMap,
             title = "SentrixPosition",
             cex = 0.6)
    }
  }
}

# Close plotting device
dev.off()

cat("SVA Sentrix/Position plot saved to: ", svaSentrixPositionPath, "\n")
cat("=======================================================================\n")

cat("Session info:\n")
print(sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")
sink()
close(logCon)
}