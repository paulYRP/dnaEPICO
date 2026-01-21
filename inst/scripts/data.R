# data.R
#
# Purpose:
#   Document how the example DNA methylation data shipped with dnaEPICO
#   (inst/extdata/GSE285813) was obtained and prepared.
#
# This script is for documentation and reproducibility purposes only.
# It is NOT executed during package installation.

## -----------------------------------------------------------
## Data source
## -----------------------------------------------------------
S
# GEO accession:
#   GSE285813
#
# Description:
#   DNA methylation profiling of whole blood samples from septic and
#   control patients using Illumina EPIC arrays.
#
# Source URL:
#   https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285813

## -----------------------------------------------------------
## Raw data retrieval (manual steps)
## -----------------------------------------------------------

# Raw IDAT files and phenotype metadata were downloaded manually from GEO
# using the supplementary files associated with GSE285813.

## -----------------------------------------------------------
## Subsetting for dnaEPICO
## -----------------------------------------------------------

# Bioconductor testing:
#
# - A small subset of samples was selected (5 samples)
# - Only the corresponding IDAT files were retained
# - A phenotype CSV file was preprocessed containing:
#     * Geo_Accession
#     * Timepoint
#     * Sex
#     * Basename
#     * Sentrix_ID
#     * Sentrix_Position
#
# This is sufficient to:
#   - Test preprocessingMinfiEwasWater()

## -----------------------------------------------------------
## Additional resources
## -----------------------------------------------------------

# Cross-reactive probe list:
#   12864_2024_10027_MOESM8_ESM.csv
#
# This file was obtained from the supplementary material of the
# corresponding publication and is included unchanged for filtering
# demonstration purposes. Reduced the size below 100 MB.
#
# Source URL:
#   https://link.springer.com/article/10.1186/s12864-024-10027-5

## -----------------------------------------------------------
## Final location in the package
## -----------------------------------------------------------

# The prepared files were placed in:
#
# inst/extdata/GSE285813/data/preprocessingMinfiEwasWater/
#
# These files are used exclusively for:
#   - Examples
#   - Unit tests
#   - Documentation
#
# They are NOT modified at runtime.
