# ===============================================
# USER CONFIGURATION
# ===============================================
MAKEFLAGS += --output-sync

# ===============================================
# MODELS SELECTION PARALLEL RUN
# ===============================================
MODEL ?= model1 # Single model
MODELS = modelA modelB modelC # Multiple models for parallel run

# ===============================================
# PER-MODEL OVERRIDES
# ===============================================
ifeq ($(MODEL), modelA)
  PHENO_FILE = $(DATA_DIR)/$(STEP1)/phenoA.csv


else ifeq ($(MODEL), modelB)
  PHENO_FILE = $(DATA_DIR)/$(STEP1)/phenoB.csv

else ifeq ($(MODEL), modelC)
  PHENO_FILE = $(DATA_DIR)/$(STEP1)/phenoC.csv

else # default (model1)
  PHENO_FILE = $(DATA_DIR)/$(STEP1)/pheno.csv

endif

# ===============================================
# DIRECTORIES
# ===============================================
LOGS_DIR = logs
DATA_DIR = data
RDATA_DIR = rData
RESULTS_DIR = results
FIGURES_DIR = figures
STEP1 = preprocessingMinfiEwasWater
STEP2 = svaEnmix
STEP3 = preprocessingPheno
STEP4 = methylationGLM_T1
STEP5 = methylationGLMM_T1T2
STEP6 = reports
METRICS_DIR = metrics
IDAT_DIR = idats
OBJ_DIR = objects
MERGE_DIR = mergeData
MODEL_DIR = models
CPG_DIR = cpgs
SUMMARY_DIR = summary
ENMIX_DIR = enmix
SVA_DIR = sva
QC_DIR = qc
R_DIR = ~/R/x86_64-pc-linux-gnu-library/4.4 # NULL is default R library path

# ===============================================
# GLOBAL PARAMETERS
# ===============================================
# STEP 1-5 PARAMETERS
PHENO_FILE = $(DATA_DIR)/$(STEP1)/phenoEMD.csv
SEP_TYPE = ""
SAMPLE_ID = Sample_Name
N_SAMPLES = 30
ARRAY_TYPE = IlluminaHumanMethylationEPICv2
ANNOTATION_VERSION = 20a1.hg38
TIFF_WIDTH = 2000
TIFF_HEIGHT = 1000
TIFF_RES = 150
SEX_COLUMN = Sex
SENTRIX_ID_COLUMN = Sentrix_ID
SENTRIX_POSITION_COLUMN = Sentrix_Position
BASENAME_COLUMN = Basename
TIME_VAR = Timepoint
PHENO_ORDER = $(SAMPLE_ID);$(TIME_VAR);$(SEX_COLUMN);PredSex;$(BASENAME_COLUMN);$(SENTRIX_ID_COLUMN);$(SENTRIX_POSITION_COLUMN)

# STEP 4-5 PARAMETERS
PHENOTYPES = TreatmentGroup
COVARIATES = Sex
FACTOR_VARS = Sex,TreatmentGroup
CPG_PREFIX = cg
CPG_LIMIT = NA
PRS_MAP = NULL
SUMMARY_PVAL = NA
N_CORES = 64
SAVE_TXT_SUMMARIES = TRUE
CHUNK_SIZE = 10000
FDR_THRESHOLD = 0.05
PADJ_METHOD = fdr
ANNOTATION_PACKAGE = IlluminaHumanMethylationEPICv2anno.20a1.hg38
ANNOTATION_COLS = Name,chr,pos,UCSC_RefGene_Group,UCSC_RefGene_Name,Relation_to_Island,GencodeV41_Group

# ===============================================
# STEP 1 PARAMETERS
# ===============================================
QC_CUTOFF = 10.5
DET_PTYPE = m+u
DET_PTHRESHOLD = 0.05
FUNNORM_SEED = 123
PVAL_THRESHOLD = 0.01
CHR_TO_REMOVE = chrX,chrY
SNPS_TO_REMOVE = SBE,CpG
CROSS_REACTIVE_FILE = $(DATA_DIR)/$(STEP1)/12864_2024_10027_MOESM8_ESM.csv
MAF_THRESHOLD = 0.1
PLOT_GROUP_VAR = TreatvControl_Time1_vs_Time2
LC_REF = salivaEPIC

# ===============================================
# STEP 2 PARAMETERS
# ===============================================
CTRL_SVA_PERC_VAR = 0.90
CTRL_SVA_FLAG = 1

# ===============================================
# STEP 3 PARAMETERS
# ===============================================
TIMEPOINTS = 1,2
COMBINE_TIMEPOINTS = 1,2

# ===============================================
# STEP 4 PARAMETERS
# ===============================================
GLM_LIBS = glm2
INTERACTION_GLM = NULL
SUMMARY_RESIDUAL_SD = TRUE
SAVE_SIGNIFICANT_CPGS = TRUE
SIGNIFICANT_CPG_PVAL = 0.1

# ==============================================
# STEP 5 PARAMETERS
# ==============================================
PERSON_VAR = person
LME_LIBS = lme4,lmerTest
INTERACTION_LME = TreatvControl_Time1_vs_Time2
SAVE_SIGNIFICANT_INTERACTIONS = TRUE
SIGNIFICANT_INTERACTION_PVAL = 0.1

# ===============================================
# ARGUMENT NORMALISATION
# ===============================================

ifeq ($(PRS_MAP),NULL)
  PRS_MAP_ARG = NULL
else
  PRS_MAP_ARG = '$(PRS_MAP)'
endif

ifeq ($(R_DIR),NULL)
  R_DIR_ARG = NULL
else
  R_DIR_ARG = '$(R_DIR)'
endif

ifeq ($(INTERACTION_GLM),NULL)
  INTERACTION_GLM_ARG = NULL
else
  INTERACTION_GLM_ARG = '$(INTERACTION_GLM)'
endif

ifeq ($(INTERACTION_LME),NULL)
  INTERACTION_LME_ARG = NULL
else
  INTERACTION_LME_ARG = '$(INTERACTION_LME)'
endif

# ===============================================
# INCLUDE PIPELINE FROM dnapipeR
# ===============================================
DNAPIPE_MK := $(shell Rscript -e "cat(system.file('extdata','make','Makefile.pipeline',package='dnaEPICO'))")
include $(DNAPIPE_MK)

