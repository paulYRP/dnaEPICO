# ================================================================
# dnaEPICO - Global variable declarations for NSE compatibility
# ================================================================
#
# Purpose:
# Declare variables used via Non-Standard Evaluation (NSE) so that
# R CMD check does not produce "no visible binding for global variable"
#
# This file does not modify runtime behaviour. It only informs the
# static code analyser about variables used inside:
#   - data.table operations
#   - modelling outputs
#   - annotation tables
#
# ================================================================

utils::globalVariables(c(
    # ---------- Illumina annotation / manifest ----------
    "IlmnID", "Name", "AddressA_ID", "AddressB_ID", "Next_Base",
    "Color_Channel", "Probe_Type", "Infinium_Design_Type",
    "CHR", "MAPINFO", "Rep_Num",
    "CHR_37", "MAPINFO_37",
    "CHR_hg38", "Start_hg38", "Strand_hg38",

    # ---------- Chromosome / probe mapping ----------
    "chr37", "chr38", "probe_id",
    "probe_design", "channel", "probe_type", "next_base",

    # ---------- DNAm metrics ----------
    "m", "cn",

    # ---------- Cell composition / estimateLC ----------
    "CpG", "FDR", "P.value", "Std.Error",

    # ---------- GLM / LME outputs ----------
    "Pr(>|t|)", "ResidualSD",
    "value", "expected", "observed", "log2meanBeta",

    # ---------- Pipeline variables ----------
    "idatFolder", "rBaseDir"
))
