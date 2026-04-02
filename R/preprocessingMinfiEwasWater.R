#' Run preprocessingMinfiEwasWater.R
#'
#' @param phenoFile Character. Path to phenotype CSV file.
#' @param idatFolder Character. Path to IDAT files directory.
#' @param outputLogs Character. Directory for log files.
#' @param nSamples Integer or NA. Number of samples to subset for testing.
#' @param SampleID Character. Sample identifier column name.
#' @param arrayType Character. Illumina array type.
#' @param annotationVersion Character. Annotation version.
#' @param scriptLabel Character. Label used in output naming.
#' @param baseDataFolder Character. Base directory for RData outputs.
#' @param figureBaseDir Character. Base directory for Figures outputs.
#' @param sepType Character. Field separator for phenotype file.
#' @param tiffWidth Integer. Width of TIFF plots in pixels.
#' @param tiffHeight Integer. Height of TIFF plots in pixels.
#' @param tiffRes Integer. Resolution (DPI) for TIFF plots.
#' @param qcCutoff Numeric. Quality-control cutoff threshold.
#' @param detPtype Character. Detection p-value calculation type.
#' @param detPThreshold Numeric. Detection p-value threshold.
#' @param normMethods Character. Normalization method(s).
#' @param sexColumn Character. Sex column name in phenotype data.
#' @param pvalThreshold Numeric. Probe-level p-value threshold.
#' @param chrToRemove Character. Chromosomes to remove.
#' @param snpsToRemove Character. SNP probe types to remove.
#' @param mafThreshold Numeric. Minor allele frequency threshold.
#' @param crossReactivePath Character. Path to cross-reactive probe file.
#' @param plotGroupVar Character. Variable used for grouping plots.
#' @param lcRef Character. Reference panel for cell composition.
#' @param phenoOrder Character. Semicolon-separated phenotype column order.
#' @param lcPhenoDir Character. Output directory for cell composition phenotype.
#'
#' @return
#' Invisibly returns \code{NULL}. This function is called for its side effects,
#' performing Illumina EPICv2 preprocessing, quality control, normalisation,
#' probe filtering, cell composition estimation, and writing plots, logs,
#' and RData objects to disk.
#'
#' @examples
#' tmp <- tempdir()
#' stopifnot(dir.exists(tmp))
#'
#' \donttest{
#' preprocessingMinfiEwasWater(
#'   phenoFile = "data/preprocessingMinfiEwasWater/pheno.csv",
#'   idatFolder = "data/preprocessingMinfiEwasWater/idats",
#'   outputLogs = "logs",
#'   nSamples = NA,
#'   SampleID = "Sample_Name",
#'   arrayType = "IlluminaHumanMethylationEPICv2",
#'   annotationVersion = "20a1.hg38",
#'   scriptLabel = "preprocessingMinfiEwasWater",
#'   baseDataFolder = "rData",
#'   figureBaseDir = "figures",
#'   sepType = "",
#'   tiffWidth = 2000,
#'   tiffHeight = 1000,
#'   tiffRes = 150,
#'   qcCutoff = 10.5,
#'   detPtype = "m+u",
#'   detPThreshold = 0.05,
#'   normMethods = "adjustedfunnorm",
#'   sexColumn = "Sex",
#'   pvalThreshold = 0.01,
#'   chrToRemove = "chrX,chrY",
#'   snpsToRemove = "SBE,CpG",
#'   mafThreshold = 0.1,
#'   crossReactivePath = "data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv",
#'   plotGroupVar = "Sex",
#'   lcRef = "salivaEPIC",
#'   phenoOrder = "Sample_Name;Timepoint;Sex;PredSex;Basename;Sentrix_ID;Sentrix_Position",
#'   lcPhenoDir = "data/preprocessingMinfiEwasWater"
#' )
#' }
#'
#' @export
preprocessingMinfiEwasWater <- function(
    phenoFile = "data/preprocessingMinfiEwasWater/pheno.csv",
    idatFolder = "data/preprocessingMinfiEwasWater/idats",
    outputLogs = "logs",
    nSamples = NA,
    SampleID = "Sample_Name",
    arrayType = "IlluminaHumanMethylationEPICv2",
    annotationVersion = "20a1.hg38",
    scriptLabel = "preprocessingMinfiEwasWater",
    baseDataFolder = "rData",
    figureBaseDir = "figures",
    sepType = "",
    tiffWidth = 2000,
    tiffHeight = 1000,
    tiffRes = 150,
    qcCutoff = 10.5,
    detPtype = "m+u",
    detPThreshold = 0.05,
    normMethods = "adjustedfunnorm",
    sexColumn = "Sex",
    pvalThreshold = 0.01,
    chrToRemove = "chrX,chrY",
    snpsToRemove = "SBE,CpG",
    mafThreshold = 0.1,
    crossReactivePath =
      "data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv",
    plotGroupVar = "Sex",
    lcRef = "salivaEPIC",
    phenoOrder = "Sample_Name;Timepoint;Sex;PredSex;Basename;Sentrix_ID;Sentrix_Position",
    lcPhenoDir = "data/preprocessingMinfiEwasWater"
) {

# Split comma/semicolon lists
chrToRemoveList <- strsplit(chrToRemove, ",")[[1]]
snpList         <- strsplit(snpsToRemove, ",")[[1]]
normMethodList  <- strsplit(normMethods, ";")[[1]]
# ==============================================================================

# ----------- Logging Setup -----------
logFilePath <- file.path(outputLogs, "log_preprocessingMinfiEwasWater.txt")
dir.create(outputLogs, recursive = TRUE, showWarnings = FALSE)

logCon <- file(logFilePath, open = "wt")

sink(logCon, split = TRUE)
sink(logCon, type = "message")
# ==============================================================================

# ----------- Logging Start Info -----------
cat("==== Starting", scriptLabel, "====\n")
cat("Start Time:               ", format(Sys.time()), "\n")
cat("Log file path:            ", logFilePath, "\n\n")
cat("Phenotype file:           ", phenoFile, "\n")
cat("Separator type: ", ifelse(is.null(sepType), "default (',')", sepType), "\n")
cat("IDAT folder:              ", idatFolder, "\n")
cat("nSamples limit:           ", ifelse(is.na(nSamples), "all", nSamples), "\n")
cat("SampleID column:          ", SampleID, "\n")
cat("Array type:               ", arrayType, "\n")
cat("Annotation version:       ", annotationVersion, "\n")
cat("Base RData folder:        ", baseDataFolder, "\n")
cat("Base Figure folder:        ", figureBaseDir, "\n")
cat("TIFF size (w x h @ dpi):  ", tiffWidth, " x ", tiffHeight, " @ ", tiffRes, "\n")
cat("QC cutoff (median):       ", qcCutoff, "\n")
cat("Detection P-value type:   ", detPtype, "\n\n")
cat("Detection p-value threshold:", detPThreshold, "\n")
cat("Normalization methods:    ", paste(normMethodList, collapse = ", "), "\n")
cat("Sex column:               ", sexColumn, "\n")
cat("Plot grouping variable:   ", plotGroupVar, "\n\n")
cat("Probe filtering:\n")
cat("  P-value threshold:      ", pvalThreshold, "\n")
cat("  Chromosomes to remove:  ", chrToRemove, "\n")
cat("  SNP positions filter:   ", snpsToRemove, "\n")
cat("  MAF threshold:          ", mafThreshold, "\n")
cat("  Cross-reactive file:    ", ifelse(is.null(crossReactivePath),
                                         "data/preprocessingMinfiEwasWater/12864_2024_10027_MOESM8_ESM.csv",
                                         crossReactivePath), "\n\n")
cat("Cell composition (estimateLC):\n")
cat("  Reference:              ", lcRef, "\n")
cat("  Leading pheno order:    ", phenoOrder, "\n")
# =============================================================================

# ----------- Prepare Subfolders for metrics -----------
objectDir  <- file.path(baseDataFolder, scriptLabel, "objects")
normDir    <- file.path(baseDataFolder, scriptLabel, "normObjects")
metricsDir <- file.path(baseDataFolder, scriptLabel, "metrics")
filterDir  <- file.path(baseDataFolder, scriptLabel, "filterObjects")

dir.create(objectDir, recursive = TRUE, showWarnings = FALSE)
dir.create(normDir, recursive = TRUE, showWarnings = FALSE)
dir.create(metricsDir, recursive = TRUE, showWarnings = FALSE)
dir.create(filterDir, recursive = TRUE, showWarnings = FALSE)

# ----------- Prepare Subfolders for rData/qc -----------

qcDir <- file.path(baseDataFolder, scriptLabel, "qc")

if (!dir.exists(qcDir)) {
  dir.create(qcDir, recursive = TRUE, showWarnings = FALSE)
}

# ----------- Prepare Subfolders for figures/metrics -----------

metricsFigDir <- file.path(figureBaseDir, scriptLabel, "metrics")

if (!dir.exists(metricsFigDir)) {
  dir.create(metricsFigDir, recursive = TRUE, showWarnings = FALSE)
}
# ----------- Prepare Subfolders for figures/qc -----------

qcFigDir <- file.path(figureBaseDir, scriptLabel, "qc")

if (!dir.exists(qcFigDir)) {
  dir.create(qcFigDir, recursive = TRUE, showWarnings = FALSE)
}
# ----------- Prepare Subfolders for figures/enmix -----------
# Target folder
enmixDir <- file.path(figureBaseDir, scriptLabel, "enMix")

# Create directory if missing
if (!dir.exists(enmixDir)) {
  dir.create(enmixDir, recursive = TRUE, showWarnings = FALSE)
}
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
  targets <- utils::read.csv(phenoFile, sep = sepChar)
} else {
  targets <- utils::read.csv(phenoFile)
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
print(utils::head(targets[, 1:5]))
cat("=======================================================================\n")

# ----------- Load IDAT Files into RGSet -----------
RGSet <- minfi::read.metharray.exp(
        base = idatFolder,
        targets = targets,
        extended = FALSE,
        recursive = FALSE,
        verbose = FALSE
)

# Assign custom sample names
Biobase::sampleNames(RGSet) <- targets[[SampleID]]
cat("RGSet loaded with", length(Biobase::sampleNames(RGSet)), "samples.\n")
cat("=======================================================================\n")

owd <- getwd(); on.exit(setwd(owd), add = TRUE)
setwd(enmixDir)

op <- options(bitmapType = "cairo")
on.exit(options(op), add = TRUE)

# Generate ENmix control plots (JPGs will be created in enmixDir)
ENmix::plotCtrl(RGSet)

setwd(owd)

cat("Generated ENmix control JPGs in:", enmixDir, "\n")

cat("=======================================================================\n")

# ----------- Apply Annotation -----------
Biobase::annotation(RGSet) <- c(
        array = arrayType,
        annotation = annotationVersion
)
cat("Applied annotation: ", paste(Biobase::annotation(RGSet), collapse = ", "), "\n")
cat("Manifest used:\n")
methods::show(minfi::getManifest(RGSet))
cat("=======================================================================\n")

# ----------- Save RGSet -----------
RGSetPath <- file.path(objectDir, "RGSet.RData")
save(RGSet, file = RGSetPath)
cat("RGSet saved to: ", RGSetPath, "\n")
cat("=======================================================================\n")

# ----------- Preprocess Raw (create MSet) -----------
cat("Running preprocessRaw() to generate MSet...\n")
MSet <- minfi::preprocessRaw(RGSet)
cat("MSet created with", ncol(MSet), "samples and", nrow(MSet), "probes.\n")
cat("=======================================================================\n")

# Save MSet object
MSetPath <- file.path(objectDir, "MSet.RData")
save(MSet, file = MSetPath)
cat("MSet saved to:", MSetPath, "\n")
cat("=======================================================================\n")

# Display methylated and unmethylated intensity
cat("Preview of methylated intensities:\n")
print(utils::head(minfi::getMeth(MSet)[, 1:3]))
cat("=======================================================================\n")
cat("Preview of unmethylated intensities:\n")
print(utils::head(minfi::getUnmeth(MSet)[, 1:3]))
cat("=======================================================================\n")

# ----------- Ratio Conversion and Genome Mapping -----------
cat("Converting MSet to RatioSet and GSet...\n")
RatioSet <- minfi::ratioConvert(MSet, what = "both", keepCN = TRUE)
cat("RatioSet created.\n")
print(RatioSet)
cat("=======================================================================\n")
GSet <- minfi::mapToGenome(RatioSet)
cat("GSet created.\n")
print(GSet)
cat("=======================================================================\n")

# Save RatioSet and GSet
RatioSetPath <- file.path(objectDir, "RatioSet.RData")
GSetPath <- file.path(objectDir, "GSet.RData")
save(RatioSet, file = RatioSetPath)
save(GSet, file = GSetPath)
cat("=======================================================================\n")

# ----------- Extract Methylation Metrics -----------
cat("Extracting methylation raw level metrics from GSet, these metrics are not saved...
    \n")

beta <- minfi::getBeta(GSet)
cat("Preview of beta values:\n")
print(utils::head(beta[, 1:5]))
cat("=======================================================================\n")

m <- minfi::getM(GSet)
cat("Preview of M-values:\n")
print(utils::head(m[, 1:5]))
cat("=======================================================================\n")

cn <- minfi::getCN(GSet)
cat("Preview of copy number values:\n")
print(utils::head(cn[, 1:5]))

cat("=======================================================================\n")

# ----------- Quality Control Plot (from MSet) -----------
cat("Running QC plotting from MSet object...\n")
qc <- minfi::getQC(MSet)

qcPath <- file.path(figureBaseDir, scriptLabel, "qc", "quality_control(MSet).tiff")
grDevices::tiff(filename = qcPath,
     width = tiffWidth,
     height = tiffHeight,
     res = tiffRes, type = "cairo")

minfi::plotQC(qc, badSampleCutoff = qcCutoff)
grDevices::dev.off()

cat("QC plot saved to: ", qcPath, "\n")
cat("=======================================================================\n")

# ----------- Calculate Detection P-values -----------
cat("Calculating detection p-values...\n")

detP <- minfi::detectionP(RGSet, type = detPtype)
cat("Detection p-values calculated using type: ", detPtype, "\n")

cat("Preview of detection p-values:\n")
print(utils::head(detP[, 1:5]))

detPpath <- file.path(qcDir, "detP_RGSet.RData")
save(detP, file = detPpath)
cat("Detection RData p-values saved to: ", detPpath, "\n")

detPlotPath <- file.path(figureBaseDir, scriptLabel, "qc", "detection_pvalues(RGSet).tiff")

grDevices::tiff(filename = detPlotPath,
     width = tiffWidth,
     height = tiffHeight,
     res = tiffRes, type = "cairo")
graphics::barplot(colMeans(detP),
        las=3,
        cex.names=0.8,
        ylab="Mean detection p-values")
graphics::abline(h=0.05,col="red", lwd = 2, lty = 2)
grDevices::dev.off()

cat("Detection plot p-values saved to: ", detPlotPath, "\n")
cat("=======================================================================\n")

# ----------- Remove samples based on detection P-values -----------
cat("Calculate the mean detection p-values across all samples...\n")
meanDetP <- colMeans(detP)

# === Identify Failed Samples ===
failedSamples <- names(meanDetP[meanDetP > detPThreshold])
nFailed <- length(failedSamples)
nBefore <- ncol(RGSet)

cat("Number of failed samples:", nFailed, "\n")
if (nFailed > 0) {
  cat("Failed sample IDs:\n")
  cat(paste(failedSamples, collapse = ", "), "\n")
}

# === Remove Failed Samples from RGSet ===
RGSet <- RGSet[, !(colnames(RGSet) %in% failedSamples)]
nAfter <- ncol(RGSet)

cat("Samples before filtering:", nBefore, "\n")
cat("Samples after filtering:", nAfter, "\n")

# ----------- Save RGSet -----------
RGSetPath <- file.path(objectDir, "RGSet.RData")
save(RGSet, file = RGSetPath)
cat("RGSet saved after removing the failed samples to: ", RGSetPath, "\n")
cat("=======================================================================\n")

# ----------- Density Plot of Beta Values from MSet -----------
cat("Generating density plot of Beta values...\n")

phenoData <- Biobase::pData(MSet)

# Ensure output directory exists
denBetaPath <- file.path(figureBaseDir,
                         scriptLabel, "qc", "densityBeta(MSet).tiff")

grDevices::tiff(filename = denBetaPath,
     width = tiffWidth,
     height = tiffHeight,
     res = tiffRes, type = "cairo")

minfi::densityPlot(MSet,
            sampGroups = phenoData[[plotGroupVar]],
            pal = RColorBrewer::brewer.pal(8, "Dark2"),
            main = paste("Density Plot of Beta Values by", plotGroupVar),
            add = TRUE,
            legend = TRUE)

grDevices::dev.off()

cat("Density plot saved to: ", denBetaPath, "\n")
cat("=======================================================================\n")

cat("Predicting sex based on Beta values...\n")
pSex <- minfi::getSex(GSet)
utils::head(pSex)

# -------------- Plot Sex predictions --------------
pSexPath <- file.path(figureBaseDir,
                         scriptLabel, "qc", "sexPrediction(GSet).tiff")

grDevices::tiff(filename = pSexPath,
     width = tiffWidth,
     height = tiffHeight,
     res = 70, type = "cairo")

graphics::plot(x = pSex$xMed,
     y = pSex$yMed,
     type = "n",
     xlab = "X chr, median total intensity (log2)",
     ylab = "Y chr, median total intensity (log2)")
graphics::text(x = pSex$xMed, y = pSex$yMed, labels = targets[[SampleID]],
     col = ifelse(pSex$predictedSex == "M", "deepskyblue", "deeppink3"))
graphics::legend("bottomleft", c("M", "F"), col = c("deepskyblue", "deeppink3"), pch = 16)
grDevices::dev.off()

cat("Predicted Sex plot saved to: ", pSexPath, "\n")
cat("=======================================================================\n")

# Create clinical sex plot
cat("Clinical sex values...\n")
pSexD <- as.data.frame(pSex)
pSexD <- merge(pSexD, targets, by.x="row.names", by.y = SampleID)
utils::head(pSexD[, 1:4])

# Extract sex column dynamically
sexVec <- targets[[sexColumn]]

# Identify NA values
nSexNA <- sum(is.na(sexVec))

# Identify unknown / unexpected values (character or factor only)
if (is.character(sexVec) || is.factor(sexVec)) {
  knownSex <- c("F", "Female", "f", "female", "FEMALE",
                "M", "Male", "m", "male", "MALE")
  unknownSex <- setdiff(unique(sexVec), knownSex)
} else {
  unknownSex <- character(0)
}

# Log sex integrity summary
cat("Sex column integrity check:\n")
cat("  Sex column used: ", sexColumn, "\n")
cat("  NA values:       ", nSexNA, "\n")
cat("  Unknown labels:  ",
    if (length(unknownSex) == 0) "None" else paste(unknownSex, collapse = ", "),
    "\n")

# Recode sex (F = 0, M = 1) controlled by opt

if (is.character(targets[[sexColumn]]) || is.factor(targets[[sexColumn]])) {
  targets[[sexColumn]] <-
    ifelse(targets[[sexColumn]] %in%
             c("F", "Female", "f", "female", "FEMALE"),
           0, 1)
}

# Synchronise recoded sex back into plotting data
pSexD[[sexColumn]] <- targets[[sexColumn]][
    match(pSexD$Row.names, targets[[SampleID]])
  ]

# -------------- Plot clinical sex --------------
pSexClPath <- file.path(figureBaseDir,
                      scriptLabel, "qc", "sexClinical(GSet).tiff")

grDevices::tiff(filename = pSexClPath,
     width = tiffWidth,
     height = tiffHeight,
     res = 70, type = "cairo")

graphics::plot(x = pSexD$xMed, y = pSexD$yMed, type = "n", xlab = "X chr, median total intensity (log2)", ylab = "Y chr, median total intensity (log2)")
graphics::text(x = pSexD$xMed, y = pSexD$yMed, labels = pSexD$Row.names,
     col = ifelse(pSexD[[sexColumn]] == "1", "deepskyblue", "deeppink3"))
graphics::legend("bottomleft", c("M", "F"), col = c("deepskyblue", "deeppink3"), pch = 16)
grDevices::dev.off()

cat("Clinical Sex plot saved to: ", pSexClPath, "\n")
cat("=======================================================================\n")

# Bind the predicted sex to the targets file and identify any mismatches
targets$PredSex <- pSex$predictedSex
# Convert F = 0 and M = 1 in the column predSex
targets$PredSex <- ifelse(targets$PredSex == "F", 0, 1)

# === Remove Failed Samples from targets ===
targets <- targets[!(targets[[SampleID]] %in% failedSamples), ]

# Add PredSex to pData
Biobase::pData(RGSet)$PredSex <- targets$PredSex

cat("Mistmaches found")
print(targets[targets[[sexColumn]] != targets$PredSex, 1:3])
cat("=======================================================================\n")

cat("Running normalization methods using Minfi and WateRmelon: ",
    paste(normMethodList, collapse = ", "), "\n")

sexVec <- NULL
if (!is.null(sexColumn) && sexColumn %in% colnames(Biobase::pData(RGSet))) {
  sexVec <- Biobase::pData(RGSet)[, sexColumn]
} else {
  cat("Note: sexColumn not found in Biobase::pData(RGSet).
      Fallback to NULL; funnorm/adjustedfunnorm will run without sex covariate.\n")
}

normPaths <- c(); firstMethod <- TRUE
for (method in normMethodList) {
        cat("Applying normalization:", method, "\n")       

        normObj <- switch(
                method,
                "adjustedfunnorm" = wateRmelon::adjustedFunnorm(RGSet, sex = sexVec),
                "funnorm"         = minfi::preprocessFunnorm(RGSet, sex = sexVec),
                "illumina"        = minfi::preprocessIllumina(RGSet),
                "quantile"        = minfi::preprocessQuantile(RGSet, sex = sexVec),
                "swan"            = minfi::preprocessSWAN(RGSet),
                stop(paste("Unknown normalization method:", method))
        )
        if (method %in% c("funnorm","adjustedfunnorm") && is.null(sexVec)) {
          cat("Requested method uses sex, but sex not provided;
              proceeded with sex = NULL.\n")
        }

        if (firstMethod) {
                MSetF <- normObj
                firstMethod <- FALSE
        }

        normPath <- file.path(normDir, paste0("norm_", method, "_RGSet.RData"))
        save(normObj, file = normPath)
        normPaths <- c(normPaths, normPath)
        cat("Saved normalized object: ", normPath, "\n")
}

# -------------- Plot Row vs Normalise data --------------
rawNormlPath <- file.path(figureBaseDir,
                        scriptLabel, "qc", "sexComparison_RawNorm(MSetF).tiff")

grDevices::tiff(filename = rawNormlPath,
     width = tiffWidth,
     height = tiffHeight,
     res = tiffRes, type = "cairo")

graphics::par(mfrow=c(1,2))
minfi::densityPlot(RGSet,
            sampGroups=targets[[sexColumn]],
            main="Raw",
            legend=FALSE)
graphics::legend("top",
       legend = levels(factor(targets[[sexColumn]])),
       text.col=RColorBrewer::brewer.pal(8,"Dark2"))

minfi::densityPlot(minfi::getBeta(MSetF),
            sampGroups=targets[[sexColumn]],
            main="Normalized",
            legend=FALSE)
graphics::legend("top",
       legend = levels(factor(targets[[sexColumn]])),
       text.col=RColorBrewer::brewer.pal(8,"Dark2"))
grDevices::dev.off()

cat("Plot Raw vs Normalisation data saved to: ", rawNormlPath, "\n")
cat("=======================================================================\n")

# ----------- Probe Filtering Based on Detection P-values -----------
cat("Filtering probes with detection p-values: ",
    pvalThreshold, "...\n")

# Recompute detection p-values
detP <- minfi::detectionP(RGSet)

# Align detection p-values with normalized probes
detP <- detP[match(Biobase::featureNames(MSetF), rownames(detP)), ]

# Identify probes retained across all samples
keep <- rowSums(detP < pvalThreshold) == ncol(MSetF)
cat("Probes retained: ", sum(keep), "/", length(keep), "\n")

MSetF_Flt <- MSetF[keep, ]
MSetFfltPath <- file.path(filterDir, "removProbes_MSetF_Flt.RData")
save(MSetF_Flt, file = MSetFfltPath)
cat("Filtered object saved to: ", MSetFfltPath, "\n")
cat("=======================================================================\n")

# ----------- Filter Probes on Sex Chromosomes -----------
cat("Removing probes on chromosomes: ", paste(chrToRemoveList,
                                              collapse = ", "), "\n")
# Identify probes to remove
ann <- minfi::getAnnotation(RGSet)
removeProbes <- ann$Name[ann$chr %in% chrToRemoveList]
keepChr <- !(Biobase::featureNames(MSetF_Flt) %in% removeProbes)

MSetF_Flt_Rxy <- MSetF_Flt[keepChr, ]

cat("Remaining probes after removing selected chromosomes:\n")
print(table(minfi::getAnnotation(MSetF_Flt_Rxy)$chr))

rxyPath <- file.path(filterDir, "removChrXY_MSetF_Flt_Rxy.RData")
save(MSetF_Flt_Rxy, file = rxyPath)
cat("Sex chromosome-filtered object saved to: ", rxyPath, "\n")
cat("=======================================================================\n")

# ----------- Remove Probes with SNPs -----------
cat("Removing probes with SNPs at: ", paste(snpList, collapse = ", "),
    " with MAF >=", mafThreshold, "\n")

# Apply SNP probe filtering
MSetF_Flt_Rxy_Ds <- minfi::dropLociWithSnps(
        MSetF_Flt_Rxy,
        snps = snpList,
        maf = mafThreshold
)
cat("Remaining probes after SNP filtering: ", nrow(MSetF_Flt_Rxy_Ds), "\n")

snpPath <- file.path(filterDir, paste0("removSNPs_MAF", mafThreshold,
                                       "_MSetF_Flt_Rxy_Ds.RData"))
save(MSetF_Flt_Rxy_Ds, file = snpPath)
cat("SNP-filtered object saved to: ", snpPath, "\n")
cat("=======================================================================\n")

# ----------- Remove Cross-Reactive Probes -----------
cat("Loading cross-reactive probe list from:\n", crossReactivePath, "\n")

xReactiveProbes <- utils::read.csv(crossReactivePath, stringsAsFactors = FALSE)

# Filter out cross-reactive probes
keepCr <- !(Biobase::featureNames(MSetF_Flt_Rxy_Ds) %in% xReactiveProbes$ProbeID)
cat("Probes retained after cross-reactive filter: ", sum(keepCr), "\n")

MSetF_Flt_Rxy_Ds_Rc <- MSetF_Flt_Rxy_Ds[keepCr, ]
rcPath <- file.path(filterDir, "removCrossReactive_MSetF_Flt_Rxy_Ds_Rc.RData")
save(MSetF_Flt_Rxy_Ds_Rc, file = rcPath)
cat("Cross-reactive-filtered object saved to: ", rcPath, "\n")
cat("=======================================================================\n")

# ----------- Final DNAm Matrices from Filtered Data -----------
cat("Extracting final DNAm matrices (M, Beta, CN)...\n")

# M-values
m <- minfi::getM(MSetF_Flt_Rxy_Ds_Rc)
mOut <- file.path(metricsDir, "m_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
save(m, file = mOut)
cat("M-values saved to: ", mOut, "\n")
print(utils::head(m[, 1:5]))

# Beta-values
beta <- minfi::getBeta(MSetF_Flt_Rxy_Ds_Rc)
betaOut <- file.path(metricsDir, "beta_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
save(beta, file = betaOut)
cat("Beta-values saved to: ", betaOut, "\n")
print(utils::head(beta[, 1:5]))

# CN-values
cn <- minfi::getCN(MSetF_Flt_Rxy_Ds_Rc)
cnOut <- file.path(metricsDir, "cn_NomFilt_MSetF_Flt_Rxy_Ds_Rc.RData")
save(cn, file = cnOut)
cat("CN-values saved to: ", cnOut, "\n")
print(utils::head(cn[, 1:5]))
cat("=======================================================================\n")

# ----- Examine higher dimensions to look at other sources of variation -----

groupFactor <- factor(targets[[plotGroupVar]])
groupSex <- factor(targets[[sexColumn]])

mdsPath <- file.path(figureBaseDir,
                          scriptLabel,
                          "metrics",
                          "examineMDS_PostFilteringCrossRect(MSetF_Flt_Rxy_Ds_Rc).tiff")

grDevices::tiff(filename = mdsPath,
     width = tiffWidth,
     height = tiffHeight,
     res = tiffRes, type = "cairo")

pal <- RColorBrewer::brewer.pal(8,"Dark2")
graphics::par(mfrow=c(1,2))
limma::plotMDS(minfi::getM(MSetF_Flt_Rxy_Ds_Rc),
        main="Timepoint",
        top=1000, gene.selection="common",
        col=pal[groupFactor], dim=c(1,2))
graphics::legend("right", legend=levels(groupFactor),
       text.col = RColorBrewer::brewer.pal(8,"Dark2"),
       cex=0.7, bg="white")
limma::plotMDS(minfi::getM(MSetF_Flt_Rxy_Ds_Rc),
        main="Sex",
        top=1000, gene.selection="common",
        col=pal[groupSex], dim=c(2,3))
graphics::legend("topright", legend=levels(groupSex),
       text.col = RColorBrewer::brewer.pal(8,"Dark2"),
       cex=0.7, bg="white")

grDevices::dev.off()

cat("Plot examineMDS_PostFilteringCrossRect saved to: ", mdsPath, "\n")
cat("=======================================================================\n")

# ----------- Plot Density of Final Beta & M Values by Group Variable -----------
cat("Plotting final density plots for grouping variable: ",
    plotGroupVar, "\n")

betaMPlotPath <- file.path(figureBaseDir,
                     scriptLabel,
                     "metrics",
                     "densityBeta&M(MSetF_Flt_Rxy_Ds_Rc).tiff")

# Create TIFF output
grDevices::tiff(filename = betaMPlotPath,
     width = tiffWidth,
     height = tiffHeight,
     res = tiffRes, type = "cairo")
graphics::par(mfrow = c(1, 2))

# Beta plot
minfi::densityPlot(beta,
            sampGroups = groupFactor,
            main = "Beta values",
            legend = FALSE,
            xlab = "Beta values")
graphics::legend("top", legend = levels(groupFactor), text.col = RColorBrewer::brewer.pal(8,"Dark2"))

# M-value plot
minfi::densityPlot(m,
            sampGroups = groupFactor,
            main = "M-values",
            legend = FALSE,
            xlab = "M values")
graphics::legend("topleft", legend = levels(groupFactor), text.col = RColorBrewer::brewer.pal(8,"Dark2"))

grDevices::dev.off()
cat("Density plots saved to: ", betaMPlotPath, "\n")

cat("=======================================================================\n")

# ------ Cell Type Estimation (Reference-driven, auto-install) -------

cat("Cell composition reference selected:", lcRef, "\n")

ewasRefs <- c(
  "saliva","salivaEPIC"
)

useLC <- any(sapply(ewasRefs, grepl, x = lcRef, fixed = TRUE))

if (useLC) {

  cat("Using internal Houseman implementation (estimateLC)\n")

  lc <- estimateLC(
    meth = beta,
    ref = lcRef,
    constrained = FALSE
  )

} else {

  cat("Using ENmix Houseman-based cell composition\n")

  lc <- ENmix::estimateCellProp(
    userdata = beta,
    refdata = lcRef,
    nonnegative = TRUE,
    normalize = FALSE,
    nProbes = 50,
    refplot = FALSE
  )
}

phenoLC <- phenoLC <- cbind(targets, lc)[, !duplicated(colnames(cbind(targets, lc)))]

leadCols <- strsplit(phenoOrder, ";", fixed = TRUE)[[1]]
leadCols <- leadCols[leadCols %in% colnames(phenoLC)]
phenoLC <- dplyr::select(phenoLC, dplyr::all_of(leadCols), dplyr::everything())

if (!dir.exists(lcPhenoDir)) dir.create(lcPhenoDir, recursive = TRUE)
lcPhenoOut <- file.path(lcPhenoDir, "phenoLC.csv")
utils::write.csv(phenoLC,
          file = lcPhenoOut,
          row.names = FALSE)
cat("Saved phenoLC:", lcPhenoOut, "\n")
cat("=======================================================================\n")

cat("Session info:\n")
print(utils::sessionInfo())
# ==============================================================================

# ----------- Close Logging -----------
sink(type = "message")
sink()
close(logCon)

}
