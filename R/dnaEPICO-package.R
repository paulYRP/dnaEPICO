#' dnaEPICO: DNA methylation preprocessing and modeling workflows
#'
#' The `dnaEPICO` package provides a structured workflow for preprocessing and
#' analyzing Illumina DNA methylation array data, including quality control,
#' normalization, cell-type estimation, surrogate-variable analysis, phenotype
#' preparation, CpG-wise generalised linear models, longitudinal mixed-effects
#' models using `lmerTest`/`lme4` or `nlme` with optional residual correlation
#' structures, and interactive web reporting.
#'
#' The package supports two complementary usage styles:
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
