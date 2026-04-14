# Data Description

# Purpose:
#   Document how external and internal data files used by the dnaEPICO
# package were obtained, prepared, and stored for reproducibility.
#
# This script is for documentation and reproducibility purposes only.
# It is NOT executed during package installation.

# -----------------------------------------------------------
# Data sources
# -----------------------------------------------------------

# 1. Cross-reactive probe list:
#   12864_2024_10027_MOESM8_ESM.csv
#
# This file was obtained from the supplementary material of the
# corresponding publication and is included for probe filtering
# demonstration purposes. The file size was reduced (<100 MB) to
# comply with package distribution requirements.
#
# Source URL:
#   https://link.springer.com/article/10.1186/s12864-024-10027-5
#
#
# 2. Cell composition reference coefficients (Houseman method):
#
#   saliva.txt
# salivaEPIC.txt
#
# These reference coefficient matrices are used by the function
# `estimateLC()` to estimate cell type proportions from DNA
# methylation beta values using the Houseman algorithm.
#
# The reference data originate from the ewastools pipeline:
#
#   Murat K, et al. (2020)
# "Ewastools: Infinium Human Methylation BeadChip pipeline for
#    population epigenetics integrated into Galaxy."
# GigaScience, 9(5):giaa049.
# https://academic.oup.com/gigascience/article/9/5/giaa049/5836679
#
# The saliva reference panel corresponds to:
#
#   Middleton LYM, et al. (2020)
# "Saliva cell type DNA methylation reference panel for
#    epidemiology studies in children."
# https://pubmed.ncbi.nlm.nih.gov/33588693/

# -----------------------------------------------------------
# Make/ directory
# -----------------------------------------------------------

# The `make/` directory contains workflow orchestration files used to
# run the complete dnaEPICO analysis pipeline in a reproducible and
# automated manner. These files are not executed during package
# installation and are not required for normal package usage.
#
# Files in the `make/` directory include:
#
#   - Makefile rules that define parameters between preprocessing,
# modelling, and reporting steps
# - Makefile model pipeline configuration variables and models shared
# across scripts
#
# These files allow users and developers to reproduce the full analysis
# performed by dnaEPICO outside of the package runtime environment,
# ensuring transparency and methodological reproducibility.

# -----------------------------------------------------------
# dnamReport.Rmd template
# -----------------------------------------------------------

# `dnamReport.Rmd` is the packaged report template used by `dnamReport()`

# -----------------------------------------------------------
# Final location in the package
# -----------------------------------------------------------

# The prepared files were placed in:
#
#   inst/extdata/
#
#   - 12864_2024_10027_MOESM8_ESM.csv
# - saliva.txt
# - salivaEPIC.txt
# - dnamReport.Rmd
#
# inst/extdata/make/
#
#   - Workflow reproducible pipeline execution files
