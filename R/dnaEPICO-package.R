#' dnaEPICO: DNA methylation preprocessing and modeling workflows
#'
#' `dnaEPICO` provides workflows for preprocessing and analyzing Illumina DNA
#' methylation arrays. It supports quality control, normalization, cell-type
#' estimation, surrogate-variable analysis, phenotype preparation, CpG-wise
#' generalised linear models, longitudinal mixed-effects models using
#' `lmerTest`/`lme4` or `nlme`, and web reporting.
#'
#' The package supports two usage styles:
#' \itemize{
#'   \item interactive use, where functions return structured in-memory result
#'   objects for inspection and composition; and
#'   \item file-based pipeline use, where the same functions can write logs,
#'   plots, tables, and serialized objects when `saveOutputs = TRUE`.
#' }
#'
#' The main high-level entry points are:
#' \itemize{
#'   \item [preprocessingMinfiEwasWater()]
#'   \item [svaEnmix()]
#'   \item [preprocessingPheno()]
#'   \item [methylationGLM()]
#'   \item [methylationLME()]
#'   \item [dnamReport()]
#' }
#'
#' @name dnaEPICO
#' @docType package
#' @keywords package
"_PACKAGE"

utils::globalVariables(c(
    "SV1", "SV2", "X", "Y", "angle", "cell", "chromosome", "count",
    "cumulativePosition", "density", "diagnosticEstimate", "displayCpG",
    "displayLabel", "end", "failed", "group", "intersectionLabel", "label",
    "labelOffset", "lower", "meanDetectionP", "midpoint", "minusLog10P",
    "observations", "percentage", "proportion", "radius", "sex",
    "significance", "stage", "start", "surrogateVariable",
    "technicalFactor", "threshold", "time", "upper", "variable1",
    "variable2", "x", "xMed", "y", "yMed"
))
