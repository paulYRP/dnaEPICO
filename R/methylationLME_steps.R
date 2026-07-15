#' Ensure a person identifier column exists for longitudinal LME analyses
#'
#' @param data Data frame containing the longitudinal phenotype-plus-beta data.
#' @param personVar Character. Name of the subject identifier column.
#' @param sidVar Character. Name of the fallback sample identifier column used to
#'   derive `personVar` when it is missing.
#'
#' @return A list containing the updated data frame, whether `personVar` was
#'   created, an optional SID-to-person preview, and counts per person.
#'
#' @description
#' Internal helper that preserves the package's legacy behavior of deriving a
#' subject identifier from `SID` when the requested `personVar` column is not
#' present.
#'
#' @keywords internal
#' @noRd
ensurePersonColumnMethylationLME <- function(
    data,
    personVar = "person",
    sidVar = "SID"
) {
  person_created <- FALSE
  mapping_preview <- NULL

  if (!(personVar %in% colnames(data))) {
    if (!(sidVar %in% colnames(data))) {
      stop(
        "Column '",
        personVar,
        "' was not found and cannot be created because '",
        sidVar,
        "' is missing.",
        call. = FALSE
      )
    }

    data[[personVar]] <- as.numeric(
      factor(gsub("[AB]$", "", as.character(data[[sidVar]])))
    )
    person_created <- TRUE
    mapping_preview <- utils::head(
      data[order(data[[personVar]], data[[sidVar]]), c(sidVar, personVar), drop = FALSE],
      20L
    )
  }

  person_counts <- table(data[[personVar]], useNA = "ifany")

  list(
    data = data,
    personCreated = person_created,
    mappingPreview = mapping_preview,
    personCounts = person_counts
  )
}

#' Build a mixed-effects formula for methylationLME helpers
#'
#' @param phenotype Character. Phenotype variable of interest.
#' @param personVar Character. Subject identifier variable used as a random
#'   intercept.
#' @param covariates Character vector of additional fixed-effect covariates.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return Character scalar containing a model formula.
#'
#' @description
#' Internal helper that builds the per-CpG linear mixed-effects formula used by
#' `methylationLME()` and its modular helpers.
#'
#' @keywords internal
#' @noRd
buildFormulaMethylationLME <- function(
    phenotype,
    personVar,
    covariates = character(0),
    interactionTerm = NULL,
    includeRandomTerm = TRUE,
    responseVar = "beta"
) {
  quoted_phenotype <- quoteNamesMethylationGLM(phenotype)
  quoted_person <- quoteNamesMethylationGLM(personVar)
  fixed_terms <- unique(c(covariates))

  if (!is.null(interactionTerm) && nzchar(interactionTerm)) {
    quoted_interaction <- quoteNamesMethylationGLM(interactionTerm)
    interaction_part <- paste(quoted_phenotype, quoted_interaction, sep = " * ")
    fixed_terms <- setdiff(fixed_terms, c(interactionTerm, phenotype))
    fixed_formula_terms <- c(interaction_part, quoteNamesMethylationGLM(fixed_terms))
  } else {
    fixed_terms <- setdiff(fixed_terms, phenotype)
    fixed_formula_terms <- c(quoted_phenotype, quoteNamesMethylationGLM(fixed_terms))
  }

  fixed_formula_terms <- unique(fixed_formula_terms[nzchar(fixed_formula_terms)])
  if (length(fixed_formula_terms) == 0L) {
    stop("At least one fixed-effect term is required.", call. = FALSE)
  }

  formula_text <- paste(
    responseVar,
    "~",
    paste(fixed_formula_terms, collapse = " + ")
  )

  if (isTRUE(includeRandomTerm)) {
    formula_text <- paste(
      formula_text,
      "+ (1 |",
      quoted_person,
      ")"
    )
  }

  formula_text
}

normalizeCorrelationStructureMethylationLME <- function(
    correlationStructure = "none"
) {
  if (length(correlationStructure) == 0L || all(is.na(correlationStructure))) {
    return("none")
  }

  structure_value <- tolower(trimws(as.character(correlationStructure[[1L]])))
  if (structure_value %in% c("", "none", "null", "na")) {
    return("none")
  }
  if (identical(structure_value, "ar1")) {
    return("AR1")
  }
  if (identical(structure_value, "car1")) {
    return("CAR1")
  }

  stop(
    "correlationStructure must be one of: none, AR1, CAR1.",
    call. = FALSE
  )
}

normalizeCorrelationVariableMethylationLME <- function(
    correlationVar = NULL
) {
  if (is.null(correlationVar) || length(correlationVar) == 0L || all(is.na(correlationVar))) {
    return(NULL)
  }

  parsed <- splitOptionMinfiEwasWater(correlationVar, sep = ",")
  if (length(parsed) == 0L) {
    return(NULL)
  }
  if (length(parsed) > 1L) {
    stop("correlationVar must contain a single variable name.", call. = FALSE)
  }

  parsed <- trimws(as.character(parsed[[1L]]))
  if (!nzchar(parsed) || tolower(parsed) %in% c("null", "na")) {
    return(NULL)
  }

  parsed
}

resolveLmeLibrariesMethylationLME <- function(
    lmeLibs = "lme4,lmerTest"
) {
  requested <- splitOptionMinfiEwasWater(lmeLibs, sep = ",")
  if (length(requested) == 0L) {
    requested <- c("lme4", "lmerTest")
  }

  requested_lower <- tolower(requested)
  has_nlme <- "nlme" %in% requested_lower
  has_lme4 <- any(requested_lower %in% c("lme4", "lmertest"))

  if (isTRUE(has_nlme) && isTRUE(has_lme4)) {
    stop(
      "lmeLibs must choose either 'lme4,lmerTest' or 'nlme', not both.",
      call. = FALSE
    )
  }

  if (isTRUE(has_nlme)) {
    return(list(
      engine = "nlme",
      requestedPackages = requested,
      requiredPackages = "nlme"
    ))
  }

  if (isTRUE(has_lme4)) {
    return(list(
      engine = "lme4",
      requestedPackages = requested,
      requiredPackages = unique(c(requested, "lme4", "lmerTest"))
    ))
  }

  stop(
    "lmeLibs must contain either 'lme4'/'lmerTest' or 'nlme'.",
    call. = FALSE
  )
}

coerceCorrelationTimeMethylationLME <- function(x) {
  if (is.numeric(x)) {
    return(as.numeric(x))
  }

  character_x <- trimws(as.character(x))
  non_missing <- !is.na(x)
  numeric_pattern <- "^[+-]?(?:\\d+(?:\\.\\d*)?|\\.\\d+)(?:[eE][+-]?\\d+)?$"
  is_numeric_text <- grepl(numeric_pattern, character_x)

  if (all(is_numeric_text[non_missing])) {
    numeric_x <- rep(NA_real_, length(x))
    numeric_x[non_missing] <- as.numeric(character_x[non_missing])
    return(numeric_x)
  }

  as.numeric(factor(x, levels = unique(x)))
}

addCorrelationTimeVariableMethylationLME <- function(
    modelData,
    correlationVar
) {
  correlation_time_var <- "dnaEPICO_lme_correlation_time"
  while (correlation_time_var %in% colnames(modelData)) {
    correlation_time_var <- paste0(correlation_time_var, "_")
  }

  modelData[[correlation_time_var]] <- coerceCorrelationTimeMethylationLME(
    modelData[[correlationVar]]
  )

  list(
    data = modelData,
    correlationTimeVar = correlation_time_var
  )
}

buildCorrelationMethylationLME <- function(
    correlationStructure = "none",
    correlationTimeVar = NULL,
    personVar
) {
  if (identical(correlationStructure, "none")) {
    return(NULL)
  }

  if (is.null(correlationTimeVar) || !nzchar(correlationTimeVar)) {
    stop(
      "An internal correlation time variable is required for AR1/CAR1 models.",
      call. = FALSE
    )
  }

  correlation_formula <- stats::as.formula(paste(
    "~",
    quoteNamesMethylationGLM(correlationTimeVar),
    "|",
    quoteNamesMethylationGLM(personVar)
  ))

  switch(
    correlationStructure,
    AR1 = nlme::corAR1(form = correlation_formula),
    CAR1 = nlme::corCAR1(form = correlation_formula),
    stop("Unsupported correlation structure.", call. = FALSE)
  )
}

coerceCoefficientTableMethylationLME <- function(
    coefTable,
    engine = "lme4"
) {
  coef_table <- as.data.frame(coefTable, check.names = FALSE)

  if (identical(engine, "nlme")) {
    rename_map <- c(
      Value = "Estimate",
      Std.Error = "Std. Error",
      DF = "df",
      `t-value` = "t value",
      `p-value` = "Pr(>|t|)"
    )
    for (old_name in names(rename_map)) {
      if (old_name %in% colnames(coef_table)) {
        colnames(coef_table)[colnames(coef_table) == old_name] <- rename_map[[old_name]]
      }
    }
  }

  as.matrix(coef_table)
}

#' Find coefficient rows for a longitudinal phenotype effect or interaction
#'
#' @param coefNames Character vector of coefficient names.
#' @param phenotype Character. Phenotype variable of interest.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return Character vector of matching coefficient names.
#'
#' @description
#' Internal helper that identifies the longitudinal fixed-effect terms that
#' should be extracted from each fitted mixed-effects model.
#'
#' @keywords internal
#' @noRd
findCoefficientRowsMethylationLME <- function(
    coefNames,
    phenotype,
    interactionTerm = NULL
) {
  normalized_names <- gsub("`", "", coefNames, fixed = TRUE)
  phenotype_pattern <- escapeRegexMethylationGLM(phenotype)

  if (!is.null(interactionTerm) && nzchar(interactionTerm)) {
    interaction_pattern <- escapeRegexMethylationGLM(interactionTerm)
    matches <- grepl(
      paste0("^", phenotype_pattern, ".*:", interaction_pattern),
      normalized_names
    )
  } else {
    matches <- grepl(paste0("^", phenotype_pattern), normalized_names)
  }

  coefNames[matches]
}

#' Fit a single CpG-level mixed-effects model for methylationLME helpers
#'
#' @param cpg Character. CpG column name.
#' @param data Data frame containing phenotype and beta columns.
#' @param modelVars Character vector of variables to retain from `data`.
#' @param personVar Character. Subject identifier variable.
#' @param formulaText Character scalar containing the model formula.
#' @param factorVars Character vector of variables to coerce to factors.
#'
#' @return A list containing coefficient, fitted-value, residual, random-effect,
#'   and fixed-effect components, or an error object when fitting fails.
#'
#' @description
#' Internal helper that fits one per-CpG linear mixed-effects model with a
#' subject-level random intercept.
#'
#' @keywords internal
#' @noRd
fitCpGModelMethylationLME <- function(
    cpg,
    data,
    modelVars,
    personVar,
    formulaText,
    factorVars = character(0),
    lmeEngine = "lme4",
    correlationStructure = "none",
    correlationTimeVar = NULL,
    responseVar = "beta"
) {
  tryCatch(
    {
      model_data <- data[, modelVars, drop = FALSE]
      model_data[[responseVar]] <- as.numeric(data[[cpg]])
      model_data[[personVar]] <- as.factor(model_data[[personVar]])

      for (var in intersect(factorVars, colnames(model_data))) {
        model_data[[var]] <- as.factor(model_data[[var]])
      }

      if (identical(lmeEngine, "nlme")) {
        fit <- nlme::lme(
          fixed = stats::as.formula(formulaText),
          random = stats::as.formula(paste("~ 1 |", quoteNamesMethylationGLM(personVar))),
          correlation = buildCorrelationMethylationLME(
            correlationStructure = correlationStructure,
            correlationTimeVar = correlationTimeVar,
            personVar = personVar
          ),
          data = model_data,
          na.action = stats::na.exclude,
          method = "REML",
          control = nlme::lmeControl(returnObject = TRUE)
        )
        coef_table <- coerceCoefficientTableMethylationLME(
          summary(fit)$tTable,
          engine = lmeEngine
        )
        ranef_values <- nlme::ranef(fit)
        fixef_values <- nlme::fixef(fit)
      } else {
        fit <- lmerTest::lmer(
          formula = stats::as.formula(formulaText),
          data = model_data,
          na.action = stats::na.exclude,
          REML = TRUE
        )
        coef_table <- coerceCoefficientTableMethylationLME(
          summary(fit)$coefficients,
          engine = lmeEngine
        )
        ranef_values <- lme4::ranef(fit)
        fixef_values <- lme4::fixef(fit)
      }

      list(
        coef = coef_table,
        residuals = stats::residuals(fit),
        fitted = stats::fitted(fit),
        ranef = ranef_values,
        fixef = fixef_values
      )
    },
    error = function(error) {
      structure(
        list(error = conditionMessage(error)),
        class = "dnaEPICO_methylationLME_fit_error"
      )
    }
  )
}

#' Summarize a single CpG-level mixed-effects model fit
#'
#' @param cpg Character. CpG identifier.
#' @param modelObj List returned by `fitCpGModelMethylationLME()`.
#' @param phenotype Character. Phenotype variable of interest.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return A data frame for the requested CpG, or `NULL`.
#'
#' @description
#' Internal helper that extracts phenotype-specific fixed-effect rows from a
#' single CpG-level mixed-effects model fit.
#'
#' @keywords internal
#' @noRd
summarizeCpGFitMethylationLME <- function(
    cpg,
    modelObj,
    phenotype,
    interactionTerm = NULL
) {
  if (is.null(modelObj) || inherits(modelObj, "dnaEPICO_methylationLME_fit_error")) {
    return(NULL)
  }

  coef_table <- modelObj$coef
  if (is.null(coef_table)) {
    return(NULL)
  }

  matched_terms <- findCoefficientRowsMethylationLME(
    coefNames = rownames(coef_table),
    phenotype = phenotype,
    interactionTerm = interactionTerm
  )
  if (length(matched_terms) == 0L) {
    return(NULL)
  }

  do.call(
    rbind,
    lapply(
      matched_terms,
      function(term) {
        coef_row <- coef_table[term, ]
        data.frame(
          CpG = cpg,
          Interaction.Term = term,
          Estimate = unname(coef_row["Estimate"]),
          Std.Error = unname(coef_row["Std. Error"]),
          t.value = unname(coef_row["t value"]),
          P.value = unname(coef_row["Pr(>|t|)"]),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }
    )
  )
}

filterSummaryByPvalueMethylationLME <- function(
    summaryDf,
    pValueFilter
) {
  summary_df <- summaryDf
  if (is.null(summary_df) || nrow(summary_df) == 0L) {
    return(data.frame())
  }

  if (nrow(summary_df) > 0L && !is.na(pValueFilter)) {
    summary_df <- summary_df[summary_df$P.value < pValueFilter, , drop = FALSE]
  }
  rownames(summary_df) <- NULL

  summary_df
}

fitCpGModelMethylationLMEPrepared <- function(
    cpg,
    cpgValues,
    modelData,
    formulaText,
    personVar,
    lmeEngine = "lme4",
    correlationStructure = "none",
    correlationTimeVar = NULL,
    responseVar = "beta"
) {
  tryCatch(
    {
      model_data <- modelData
      model_data[[responseVar]] <- as.numeric(cpgValues)

      if (identical(lmeEngine, "nlme")) {
        fit <- nlme::lme(
          fixed = stats::as.formula(formulaText),
          random = stats::as.formula(paste("~ 1 |", quoteNamesMethylationGLM(personVar))),
          correlation = buildCorrelationMethylationLME(
            correlationStructure = correlationStructure,
            correlationTimeVar = correlationTimeVar,
            personVar = personVar
          ),
          data = model_data,
          na.action = stats::na.exclude,
          method = "REML",
          control = nlme::lmeControl(returnObject = TRUE)
        )
        coef_table <- coerceCoefficientTableMethylationLME(
          summary(fit)$tTable,
          engine = lmeEngine
        )
        ranef_values <- nlme::ranef(fit)
        fixef_values <- nlme::fixef(fit)
      } else {
        fit <- lmerTest::lmer(
          formula = stats::as.formula(formulaText),
          data = model_data,
          na.action = stats::na.exclude,
          REML = TRUE
        )
        coef_table <- coerceCoefficientTableMethylationLME(
          summary(fit)$coefficients,
          engine = lmeEngine
        )
        ranef_values <- lme4::ranef(fit)
        fixef_values <- lme4::fixef(fit)
      }

      list(
        coef = coef_table,
        residuals = stats::residuals(fit),
        fitted = stats::fitted(fit),
        ranef = ranef_values,
        fixef = fixef_values
      )
    },
    error = function(error) {
      structure(
        list(error = conditionMessage(error)),
        class = "dnaEPICO_methylationLME_fit_error"
      )
    }
  )
}

fitMethylationLMEBatch <- function(
    cpgBatch,
    data,
    modelData,
    formulaText,
    personVar,
    lmeEngine = "lme4",
    correlationStructure = "none",
    correlationTimeVar = NULL,
    phenotype,
    interactionTerm = NULL,
    responseVar = "beta"
) {
  fits <- vector("list", length(cpgBatch))
  names(fits) <- cpgBatch
  summaries <- vector("list", length(cpgBatch))
  names(summaries) <- cpgBatch

  for (cpg in cpgBatch) {
    model_obj <- fitCpGModelMethylationLMEPrepared(
      cpg = cpg,
      cpgValues = data[[cpg]],
      modelData = modelData,
      formulaText = formulaText,
      personVar = personVar,
      lmeEngine = lmeEngine,
      correlationStructure = correlationStructure,
      correlationTimeVar = correlationTimeVar,
      responseVar = responseVar
    )
    fits[[cpg]] <- model_obj
    summaries[[cpg]] <- summarizeCpGFitMethylationLME(
      cpg = cpg,
      modelObj = model_obj,
      phenotype = phenotype,
      interactionTerm = interactionTerm
    )
  }

  summaries <- Filter(Negate(is.null), summaries)
  summary_df <- if (length(summaries) == 0L) {
    data.frame()
  } else {
    out <- do.call(rbind, summaries)
    rownames(out) <- NULL
    out
  }

  list(fits = fits, summaries = summary_df)
}

#' Summarize phenotype values by timepoint for longitudinal methylation analyses
#'
#' @param data Data frame containing the longitudinal phenotype-plus-beta data.
#' @param timeVar Character. Name of the time variable.
#' @param phenotypes Character vector or comma-separated string of phenotype
#'   variables to summarize.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A data frame with one row per timepoint and summary columns for each
#'   requested phenotype.
#'
#' @description
#' Summarize the requested phenotype variables by timepoint. Numeric phenotypes
#' are reported with mean, standard deviation, and non-missing counts; non-numeric
#' phenotypes are reported with non-missing counts and the observed levels.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' timepoint_summary <- summarizeTimepointsMethylationLME(
#'   data = ex$preparedData$data,
#'   timeVar = "Timepoint",
#'   phenotypes = "score",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' nrow(timepoint_summary)
#'
#' @export
summarizeTimepointsMethylationLME <- function(
    data,
    timeVar = "Timepoint",
    phenotypes,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  phenotype_list <- splitOptionMinfiEwasWater(phenotypes, sep = ",")

  if (!(timeVar %in% colnames(data))) {
    stop("timeVar column not found: ", timeVar, call. = FALSE)
  }

  missing_phenotypes <- setdiff(phenotype_list, colnames(data))
  if (length(missing_phenotypes) > 0L) {
    stop(
      "Phenotype columns not found for timepoint summary: ",
      paste(missing_phenotypes, collapse = ", "),
      call. = FALSE
    )
  }

  split_data <- split(data, as.character(data[[timeVar]]))
  summaries <- lapply(
    split_data,
    function(df) {
      out <- list()
      out[[timeVar]] <- as.character(df[[timeVar]][[1L]])

      for (phenotype in phenotype_list) {
        values <- df[[phenotype]]
        if (is.numeric(values)) {
          out[[paste0(phenotype, "_mean")]] <- mean(values, na.rm = TRUE)
          out[[paste0(phenotype, "_sd")]] <- stats::sd(values, na.rm = TRUE)
          out[[paste0(phenotype, "_n")]] <- sum(!is.na(values))
        } else {
          observed_levels <- unique(as.character(values[!is.na(values)]))
          out[[paste0(phenotype, "_n")]] <- sum(!is.na(values))
          out[[paste0(phenotype, "_levels")]] <- paste(observed_levels, collapse = ",")
        }
      }

      as.data.frame(out, stringsAsFactors = FALSE)
    }
  )
  summary_df <- do.call(rbind, summaries)
  rownames(summary_df) <- NULL

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      "Summary statistics for phenotype scores by timepoint:",
      previewLinesMinfiEwasWater(summary_df),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  summary_df
}

#' Prepare longitudinal phenotype-plus-methylation data for mixed-effects analyses
#'
#' @param inputPheno Character. Path to the merged longitudinal phenotype-plus-
#'   methylation object created by `preprocessingPheno()`.
#' @param personVar Character. Name of the subject identifier column.
#' @param timeVar Character. Name of the time variable.
#' @param phenotypes Character vector or comma-separated string of phenotype
#'   variables to model.
#' @param covariates Character vector or comma-separated string of covariate
#'   variables to adjust for.
#' @param factorVars Character vector or comma-separated string of variables that
#'   should be converted to factors before modeling.
#' @param prsMap Character vector or comma-separated string of phenotype-to-PRS
#'   mappings in the form `"Phenotype:PRS"`.
#' @param cpgPrefix Character. Prefix used to identify methylation columns.
#' @param cpgLimit Integer or `NA`. Maximum number of CpGs to retain. `NA`
#'   keeps all matching CpGs.
#' @param methylationScale Character. Methylation metric represented by the CpG
#'   columns. One of `"beta"`, `"m"`, or `"cn"`.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationLME_data"` containing
#'   the prepared analysis data, parsed variable selections, CpG columns,
#'   timepoint summaries, and subject-ID diagnostics.
#'
#' @description
#' Load the merged longitudinal phenotype-plus-methylation object, ensure that a subject
#' identifier column is available, validate the requested modeling variables,
#' convert selected variables to factors, and return a single in-memory object
#' for downstream mixed-effects modeling helpers.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' prepared_data <- prepareMethylationLMEData(
#'   inputPheno = ex$inputPath,
#'   personVar = "person",
#'   timeVar = "Timepoint",
#'   phenotypes = "score",
#'   covariates = "sex",
#'   factorVars = "sex",
#'   cpgLimit = 2,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(prepared_data)
#'
#' @export
prepareMethylationLMEData <- function(
    inputPheno,
    personVar = "person",
    timeVar = "Timepoint",
    phenotypes,
    covariates,
    factorVars,
    prsMap = NULL,
    cpgPrefix = "cg",
    cpgLimit = NA,
    methylationScale = "beta",
    interactionTerm = NULL,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  phenotype_list <- splitOptionMinfiEwasWater(phenotypes, sep = ",")
  covariate_list <- splitOptionMinfiEwasWater(covariates, sep = ",")
  factor_list <- splitOptionMinfiEwasWater(factorVars, sep = ",")
  prs_map <- parsePrsMapMethylationGLM(prsMap)
  cpg_limit <- normalizeOptionalNumericMethylationGLM(cpgLimit)
  methylation_scale <- normalizeMethylationScaleDnaEpico(methylationScale)
  methylation_label <- methylationScaleResponseLabelDnaEpico(methylation_scale)
  methylation_prefix <- methylationScaleObjectPrefixDnaEpico(methylation_scale)
  response_column <- methylationScaleResponseColumnDnaEpico(methylation_scale)
  analysis_data <- loadSavedObjectPreprocessingPheno(inputPheno, preferred_name = "phenoBT1T2")

  if (!is.data.frame(analysis_data)) {
    analysis_data <- as.data.frame(analysis_data, stringsAsFactors = FALSE)
  }

  person_data <- ensurePersonColumnMethylationLME(
    data = analysis_data,
    personVar = personVar
  )
  analysis_data <- person_data$data

  if (!(timeVar %in% colnames(analysis_data))) {
    stop("timeVar column not found in inputPheno: ", timeVar, call. = FALSE)
  }

  missing_phenotypes <- setdiff(phenotype_list, colnames(analysis_data))
  if (length(missing_phenotypes) > 0L) {
    stop(
      "Phenotype columns not found in inputPheno: ",
      paste(missing_phenotypes, collapse = ", "),
      call. = FALSE
    )
  }

  missing_covariates <- setdiff(covariate_list, colnames(analysis_data))
  if (length(missing_covariates) > 0L) {
    stop(
      "Covariate columns not found in inputPheno: ",
      paste(missing_covariates, collapse = ", "),
      call. = FALSE
    )
  }

  mapped_prs <- unname(prs_map[names(prs_map) %in% phenotype_list])
  missing_prs <- setdiff(mapped_prs, colnames(analysis_data))
  if (length(missing_prs) > 0L) {
    stop("PRS columns not found in inputPheno: ", paste(unique(missing_prs), collapse = ", "), call. = FALSE)
  }

  resolved_interaction <- interactionTerm
  if (!is.null(resolved_interaction) && nzchar(resolved_interaction)) {
    if (!(resolved_interaction %in% colnames(analysis_data))) {
      emitLogMinfiEwasWater(
        c(
          "=======================================================================",
          paste("Requested interactionTerm was not found and will be ignored:", resolved_interaction),
          "======================================================================="
        ),
        verbose = verbose,
        log_path = log_path
      )
      resolved_interaction <- NULL
    }
  }

  for (var in intersect(c(personVar, factor_list), colnames(analysis_data))) {
    if (identical(var, personVar)) {
      next
    }
    analysis_data[[var]] <- as.factor(analysis_data[[var]])
  }

  cpg_columns <- grep(
    paste0("^", escapeRegexMethylationGLM(cpgPrefix)),
    colnames(analysis_data),
    value = TRUE
  )
  if (!is.na(cpg_limit)) {
    cpg_columns <- utils::head(cpg_columns, as.integer(cpg_limit))
  }
  if (length(cpg_columns) == 0L) {
    stop("No CpG columns were found with prefix '", cpgPrefix, "'.", call. = FALSE)
  }

  timepoint_summary <- summarizeTimepointsMethylationLME(
    data = analysis_data,
    timeVar = timeVar,
    phenotypes = phenotype_list,
    verbose = FALSE,
    logs = FALSE
  )

  log_lines <- c(
    "=======================================================================",
    paste("Loaded longitudinal phenotype + methylation data from:", inputPheno),
    paste("Merged modeling object:            ", methylation_prefix, "*"),
    paste("Data dimensions:                  ", paste(dim(analysis_data), collapse = " x ")),
    paste("Person variable:                  ", personVar),
    paste("Time variable:                    ", timeVar),
    paste("Phenotypes:                       ", paste(phenotype_list, collapse = ", ")),
    paste("Covariates:                       ", paste(covariate_list, collapse = ", ")),
    paste("Factor variables:                 ", paste(factor_list, collapse = ", ")),
    paste("CpG columns retained:             ", length(cpg_columns)),
    if (isTRUE(person_data$personCreated)) {
      paste("Created person variable from SID: ", personVar)
    } else {
      paste("Person variable already present:  ", personVar)
    },
    "Count of records per person ID:",
    previewLinesMinfiEwasWater(person_data$personCounts),
    paste("Values observed in", timeVar, ":"),
    previewLinesMinfiEwasWater(table(analysis_data[[timeVar]], useNA = "ifany"))
  )

  if (!is.null(person_data$mappingPreview)) {
    log_lines <- c(
      log_lines,
      "Example mapping of SID to person ID:",
      previewLinesMinfiEwasWater(person_data$mappingPreview)
    )
  }

  log_lines <- c(
    log_lines,
    "Summary statistics for phenotype scores by timepoint:",
    previewLinesMinfiEwasWater(timepoint_summary),
    "======================================================================="
  )
  emitLogMinfiEwasWater(log_lines, verbose = verbose, log_path = log_path)

  structure(
    list(
      data = analysis_data,
      personVar = personVar,
      timeVar = timeVar,
      phenotypes = phenotype_list,
      covariates = covariate_list,
      factorVars = factor_list,
      prsMap = prs_map,
      cpgColumns = cpg_columns,
      cpgPrefix = cpgPrefix,
      cpgLimit = cpg_limit,
      methylationScale = methylation_scale,
      responseLabel = methylation_label,
      methylationObjectPrefix = methylation_prefix,
      internalResponseColumn = response_column,
      interactionTerm = resolved_interaction,
      requestedInteractionTerm = interactionTerm,
      personCreated = person_data$personCreated,
      personCounts = person_data$personCounts,
      personMappingPreview = person_data$mappingPreview,
      timepointSummary = timepoint_summary
    ),
    class = "dnaEPICO_methylationLME_data"
  )
}
#' Fit CpG-wise mixed-effects models for longitudinal methylation analyses
#'
#' @param preparedData Object returned by `prepareMethylationLMEData()`.
#' @param nCores Integer. Number of worker processes to use.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded
#'   to worker processes.
#' @param lmeLibs Character vector or comma-separated string of package names to
#'   check on worker processes. The default is `"lme4,lmerTest"`.
#' @param correlationStructure Character. Residual correlation structure used
#'   when `lmeLibs` selects `"nlme"`. One of `"none"`, `"AR1"`, or `"CAR1"`.
#' @param correlationVar Character or `NULL`. Variable used to order repeated
#'   observations within `personVar` for `AR1` or `CAR1` residual correlation
#'   structures. Must be supplied explicitly for `AR1` or `CAR1`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationLME_models"`
#'   containing fitted model lists, model formulas, and counts of failed CpG
#'   fits.
#'
#' @description
#' Fit one linear mixed-effects model per CpG for each phenotype requested in the
#' object returned by `prepareMethylationLMEData()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' model_results <- fitMethylationLMEModels(
#'   preparedData = ex$preparedData,
#'   nCores = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(model_results$fits)
#'
#' @export
fitMethylationLMEModels <- function(
    preparedData,
    nCores = 1L,
    libPath = NULL,
    lmeLibs = "lme4,lmerTest",
    correlationStructure = "none",
    correlationVar = NULL,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)

  if (is.null(libPath)) {
    libPath <- .libPaths()
  }

  lme_config <- resolveLmeLibrariesMethylationLME(lmeLibs)
  lme_lib_list <- lme_config$requestedPackages
  required_lme_lib_list <- lme_config$requiredPackages
  lme_engine <- lme_config$engine
  correlation_structure <- normalizeCorrelationStructureMethylationLME(correlationStructure)
  correlation_var <- normalizeCorrelationVariableMethylationLME(
    correlationVar = correlationVar
  )
  if (!identical(lme_engine, "nlme") && !identical(correlation_structure, "none")) {
    stop(
      "correlationStructure can only be AR1 or CAR1 when lmeLibs selects 'nlme'.",
      call. = FALSE
    )
  }
  if (!identical(correlation_structure, "none")) {
    if (is.null(correlation_var)) {
      stop(
        "correlationVar must be supplied when correlationStructure is AR1 or CAR1.",
        call. = FALSE
      )
    }
  }

  analysis_data <- preparedData$data
  cpg_columns <- preparedData$cpgColumns
  n_cores <- max(1L, as.integer(nCores))
  fits <- list()
  summary_cache <- list()
  formulas <- stats::setNames(character(length(preparedData$phenotypes)), preparedData$phenotypes)
  failure_counts <- stats::setNames(integer(length(preparedData$phenotypes)), preparedData$phenotypes)
  backend <- resolveParallelBackendMethylationModels(n_cores)

  for (phenotype in preparedData$phenotypes) {
    prs_var <- character(0)
    if (phenotype %in% names(preparedData$prsMap)) {
      prs_var <- unname(preparedData$prsMap[[phenotype]])
    }
    covariates <- unique(c(preparedData$covariates, prs_var))
    model_vars <- unique(c(
      preparedData$personVar,
      phenotype,
      covariates,
      preparedData$interactionTerm,
      if (!identical(correlation_structure, "none")) {
        correlation_var
      } else {
        character(0)
      }
    ))
    model_vars <- model_vars[!is.na(model_vars) & nzchar(model_vars)]

    missing_vars <- setdiff(model_vars, colnames(analysis_data))
    if (length(missing_vars) > 0L) {
      stop(
        "Model variables not found for phenotype ",
        phenotype,
        ": ",
        paste(missing_vars, collapse = ", "),
        call. = FALSE
      )
    }

    display_formula_text <- buildFormulaMethylationLME(
      phenotype = phenotype,
      personVar = preparedData$personVar,
      covariates = covariates,
      interactionTerm = preparedData$interactionTerm,
      includeRandomTerm = TRUE,
      responseVar = preparedData$internalResponseColumn
    )
    formula_text <- buildFormulaMethylationLME(
      phenotype = phenotype,
      personVar = preparedData$personVar,
      covariates = covariates,
      interactionTerm = preparedData$interactionTerm,
      includeRandomTerm = !identical(lme_engine, "nlme"),
      responseVar = preparedData$internalResponseColumn
    )

    base_model_data <- analysis_data[, model_vars, drop = FALSE]
    factor_vars <- preparedData$factorVars
    person_var <- preparedData$personVar
    base_model_data[[person_var]] <- as.factor(base_model_data[[person_var]])
    for (var in intersect(factor_vars, colnames(base_model_data))) {
      base_model_data[[var]] <- as.factor(base_model_data[[var]])
    }
    correlation_time_var <- NULL
    if (!identical(correlation_structure, "none")) {
      correlation_data <- addCorrelationTimeVariableMethylationLME(
        modelData = base_model_data,
        correlationVar = correlation_var
      )
      base_model_data <- correlation_data$data
      correlation_time_var <- correlation_data$correlationTimeVar
    }
    cpg_batches <- chunkCpGColumnsMethylationModels(
      cpgColumns = cpg_columns,
      nCores = n_cores,
      batchesPerCore = 8L
    )
    batch_worker <- fitMethylationLMEBatch
    resolved_interaction <- preparedData$interactionTerm
    response_var <- preparedData$internalResponseColumn

    if (!identical(backend, "serial") && length(cpg_batches) > 1L) {
      cluster_size <- min(n_cores, length(cpg_batches))

      if (identical(backend, "fork")) {
        batch_results <- parallel::mclapply(
          cpg_batches,
          function(batch) {
            validateWorkerPackagesMethylationModels(
              libPath = libPath,
              packages = required_lme_lib_list
            )
            batch_worker(
              cpgBatch = batch,
              data = analysis_data,
              modelData = base_model_data,
              formulaText = formula_text,
              personVar = person_var,
              lmeEngine = lme_engine,
              correlationStructure = correlation_structure,
              correlationTimeVar = correlation_time_var,
              phenotype = phenotype,
              interactionTerm = resolved_interaction,
              responseVar = response_var
            )
          },
          mc.cores = cluster_size,
          mc.preschedule = FALSE
        )
      } else {
        cl <- makePsockClusterMethylationModels(cluster_size)
        batch_results <- tryCatch(
          {
            parallel::clusterExport(
              cl,
              varlist = c(
                "analysis_data",
                "base_model_data",
                "formula_text",
                "phenotype",
                "resolved_interaction",
                "person_var",
                "lme_engine",
                "correlation_structure",
                "correlation_time_var",
                "response_var",
                "batch_worker",
                "libPath",
                "required_lme_lib_list",
                "validateWorkerPackagesMethylationModels",
                "fitCpGModelMethylationLMEPrepared",
                "buildCorrelationMethylationLME",
                "coerceCoefficientTableMethylationLME",
                "summarizeCpGFitMethylationLME",
                "findCoefficientRowsMethylationLME",
                "quoteNamesMethylationGLM",
                "escapeRegexMethylationGLM"
              ),
              envir = environment()
            )

            parallel::clusterEvalQ(
              cl,
              validateWorkerPackagesMethylationModels(
                libPath = libPath,
                packages = required_lme_lib_list
              )
            )

            parallel::parLapplyLB(
              cl,
              cpg_batches,
              function(batch) {
                batch_worker(
                  cpgBatch = batch,
                  data = analysis_data,
                  modelData = base_model_data,
                  formulaText = formula_text,
                  personVar = person_var,
                  lmeEngine = lme_engine,
                  correlationStructure = correlation_structure,
                  correlationTimeVar = correlation_time_var,
                  phenotype = phenotype,
                  interactionTerm = resolved_interaction,
                  responseVar = response_var
                )
              }
            )
          },
          finally = {
            parallel::stopCluster(cl)
          }
        )
      }
    } else {
      validateWorkerPackagesMethylationModels(
        libPath = libPath,
        packages = required_lme_lib_list
      )
      batch_results <- lapply(
        cpg_batches,
        function(batch) {
          batch_worker(
            cpgBatch = batch,
            data = analysis_data,
            modelData = base_model_data,
            formulaText = formula_text,
            personVar = person_var,
            lmeEngine = lme_engine,
            correlationStructure = correlation_structure,
            correlationTimeVar = correlation_time_var,
            phenotype = phenotype,
            interactionTerm = resolved_interaction,
            responseVar = response_var
          )
        }
      )
    }

    combined_results <- combineFitBatchResultsMethylationModels(
      batchResults = batch_results,
      cpgColumns = cpg_columns
    )
    fit_list <- combined_results$fits
    phenotype_summary_cache <- filterSummaryByPvalueMethylationLME(
      summaryDf = combined_results$summaries,
      pValueFilter = NA_real_
    )
    failures <- vapply(
      fit_list,
      function(x) inherits(x, "dnaEPICO_methylationLME_fit_error"),
      logical(1)
    )

    fits[[phenotype]] <- fit_list
    summary_cache[[phenotype]] <- phenotype_summary_cache
    formulas[[phenotype]] <- display_formula_text
    failure_counts[[phenotype]] <- sum(failures)

    emitLogMinfiEwasWater(
      c(
        "=======================================================================",
        paste("Fitted phenotype:            ", phenotype),
        paste("Formula:                     ", display_formula_text),
        paste("Correlation structure:       ", correlation_structure),
        paste(
          "Correlation variable:        ",
          if (identical(correlation_structure, "none")) "None" else correlation_var
        ),
        paste("CpGs attempted:              ", length(cpg_columns)),
        paste("Failed CpG fits:             ", failure_counts[[phenotype]]),
        paste("Parallel backend:            ", backend),
        paste("Fit batches:                 ", length(cpg_batches)),
        paste("Fit batch size:              ", if (length(cpg_batches) == 0L) 0L else max(vapply(cpg_batches, length, integer(1)))),
        paste("Fit-time summary rows cached:", nrow(phenotype_summary_cache)),
        "======================================================================="
      ),
      verbose = verbose,
      log_path = log_path
    )
  }

  structure(
    list(
      fits = fits,
      summaryCache = summary_cache,
      formulas = formulas,
      phenotypes = names(fits),
      failureCounts = failure_counts,
      settings = list(
        nCores = n_cores,
        parallelBackend = backend,
        fitBatchCount = length(chunkCpGColumnsMethylationModels(cpg_columns, n_cores, 8L)),
        libPath = libPath,
        lmeLibs = lme_lib_list,
        lmeEngine = lme_engine,
        correlationStructure = correlation_structure,
        correlationVar = if (identical(correlation_structure, "none")) {
          NULL
        } else {
          correlation_var
        },
        methylationScale = preparedData$methylationScale,
        methylationObjectPrefix = preparedData$methylationObjectPrefix,
        responseLabel = preparedData$responseLabel,
        internalResponseColumn = preparedData$internalResponseColumn,
        interactionTerm = preparedData$interactionTerm
      ),
      responseLabel = preparedData$responseLabel
    ),
    class = "dnaEPICO_methylationLME_models"
  )
}

#' Summarize CpG-wise mixed-effects model fits for longitudinal analyses
#'
#' @param modelResults Object returned by `fitMethylationLMEModels()`.
#' @param preparedData Object returned by `prepareMethylationLMEData()`.
#' @param summaryPval Numeric or `NA`. Optional p-value filter applied to the
#'   returned summary tables. `NA` keeps all rows.
#' @param nCores Integer. Number of worker processes to use while extracting
#'   summary rows.
#' @param chunkSize Integer or `NULL`. Number of CpGs processed per parallel
#'   chunk. `NULL` chooses a value automatically.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationLME_summaries"`
#'   containing one CpG-level summary data frame per phenotype.
#'
#' @description
#' Extract phenotype-specific fixed-effect tables from the fitted mixed-effects
#' model object returned by `fitMethylationLMEModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' summary_results <- summarizeMethylationLMEModels(
#'   modelResults = ex$modelResults,
#'   preparedData = ex$preparedData,
#'   summaryPval = NA,
#'   nCores = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(summary_results$summaries)
#'
#' @export
summarizeMethylationLMEModels <- function(
    modelResults,
    preparedData,
    summaryPval = NA,
    nCores = 1L,
    chunkSize = NULL,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  p_value_filter <- normalizeOptionalNumericMethylationGLM(summaryPval)
  chunk_size <- normalizeChunkSizeMethylationGLM(chunkSize)
  n_cores <- max(1L, as.integer(nCores))
  summaries <- list()

  for (phenotype in names(modelResults$fits)) {
    if (
      !is.null(modelResults$summaryCache) &&
        !is.null(modelResults$summaryCache[[phenotype]])
    ) {
      summary_df <- filterSummaryByPvalueMethylationLME(
        summaryDf = modelResults$summaryCache[[phenotype]],
        pValueFilter = p_value_filter
      )
      summaries[[phenotype]] <- summary_df

      emitLogMinfiEwasWater(
        c(
          "=======================================================================",
          paste("Summarized phenotype:        ", phenotype),
          paste("LME summary rows returned:   ", nrow(summary_df)),
          "Summary source:              fit-time cache",
          if (is.na(p_value_filter)) {
            "P-value filter:              none"
          } else {
            paste("P-value filter:              ", p_value_filter)
          },
          "======================================================================="
        ),
        verbose = verbose,
        log_path = log_path
      )
      next
    }

    fit_list <- modelResults$fits[[phenotype]]
    cpg_names <- names(fit_list)

    local_chunk_size <- chunk_size
    if (is.null(local_chunk_size)) {
      local_chunk_size <- max(10L, floor(length(cpg_names) / max(n_cores * 4L, 1L)))
    }
    local_chunk_size <- max(1L, as.integer(local_chunk_size))
    cpg_chunks <- split(cpg_names, ceiling(seq_along(cpg_names) / local_chunk_size))

    summary_worker <- summarizeCpGFitMethylationLME
    resolved_interaction <- preparedData$interactionTerm

    if (n_cores > 1L && length(cpg_chunks) > 1L) {
      cluster_size <- min(n_cores, length(cpg_chunks))
      cl <- parallel::makeCluster(cluster_size)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      parallel::clusterExport(
        cl,
        varlist = c(
          "fit_list",
          "phenotype",
          "resolved_interaction",
          "summary_worker"
        ),
        envir = environment()
      )

      result_chunks <- parallel::parLapplyLB(
        cl,
        cpg_chunks,
        function(chunk) {
          rows <- lapply(
            chunk,
            function(cpg) {
              summary_worker(
                cpg = cpg,
                modelObj = fit_list[[cpg]],
                phenotype = phenotype,
                interactionTerm = resolved_interaction
              )
            }
          )
          rows <- Filter(Negate(is.null), rows)
          if (length(rows) == 0L) {
            return(NULL)
          }
          do.call(rbind, rows)
        }
      )
      parallel::stopCluster(cl)
      on.exit(NULL, add = FALSE)
    } else {
      result_chunks <- lapply(
        cpg_chunks,
        function(chunk) {
          rows <- lapply(
            chunk,
            function(cpg) {
              summary_worker(
                cpg = cpg,
                modelObj = fit_list[[cpg]],
                phenotype = phenotype,
                interactionTerm = resolved_interaction
              )
            }
          )
          rows <- Filter(Negate(is.null), rows)
          if (length(rows) == 0L) {
            return(NULL)
          }
          do.call(rbind, rows)
        }
      )
    }

    result_chunks <- Filter(Negate(is.null), result_chunks)
    if (length(result_chunks) == 0L) {
      summary_df <- data.frame()
    } else {
      summary_df <- do.call(rbind, result_chunks)
      rownames(summary_df) <- NULL
    }

    if (nrow(summary_df) > 0L && !is.na(p_value_filter)) {
      summary_df <- summary_df[summary_df$P.value < p_value_filter, , drop = FALSE]
      rownames(summary_df) <- NULL
    }

    summaries[[phenotype]] <- summary_df

    emitLogMinfiEwasWater(
      c(
        "=======================================================================",
        paste("Summarized phenotype:        ", phenotype),
        paste("LME summary rows returned:   ", nrow(summary_df)),
        paste("Summary chunk size:          ", local_chunk_size),
        if (is.na(p_value_filter)) {
          "P-value filter:              none"
        } else {
          paste("P-value filter:              ", p_value_filter)
        },
        "======================================================================="
      ),
      verbose = verbose,
      log_path = log_path
    )
  }

  structure(
    list(
      summaries = summaries,
      phenotypes = names(summaries),
      settings = list(
        summaryPval = p_value_filter,
        chunkSize = chunk_size,
        interactionTerm = preparedData$interactionTerm
      )
    ),
    class = "dnaEPICO_methylationLME_summaries"
  )
}

#' Collect significant longitudinal model terms from fitted mixed-effects models
#'
#' @param modelResults Object returned by `fitMethylationLMEModels()`.
#' @param pvalThreshold Numeric. Threshold applied to the extracted phenotype or
#'   interaction p-values.
#' @param interactionTerm Character or `NULL`. Optional interaction term. When
#'   `NULL`, phenotype main effects are used.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationLME_significant"`
#'   containing the retained coefficient tables for each phenotype.
#'
#' @description
#' Collect raw coefficient tables for CpGs whose phenotype main effect or
#' requested interaction p-value passes the chosen threshold.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' significant_hits <- collectSignificantInteractionsMethylationLME(
#'   modelResults = ex$modelResults,
#'   pvalThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(significant_hits)
#'
#' @export
collectSignificantInteractionsMethylationLME <- function(
    modelResults,
    pvalThreshold = 0.05,
    interactionTerm = NULL,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  threshold <- as.numeric(pvalThreshold[[1L]])
  retained <- list()

  for (phenotype in names(modelResults$fits)) {
    fit_list <- modelResults$fits[[phenotype]]
    phenotype_hits <- list()
    if (
      !is.null(modelResults$summaryCache) &&
        !is.null(modelResults$summaryCache[[phenotype]]) &&
        optionalTermMatchesMethylationModels(
          requested = interactionTerm,
          cached = modelResults$settings$interactionTerm
        )
    ) {
      cached_summary <- modelResults$summaryCache[[phenotype]]
      if (nrow(cached_summary) > 0L && !is.na(threshold)) {
        hit_cpgs <- unique(cached_summary$CpG[cached_summary$P.value < threshold])
        hit_cpgs <- hit_cpgs[!is.na(hit_cpgs) & hit_cpgs %in% names(fit_list)]
        for (cpg in hit_cpgs) {
          model_obj <- fit_list[[cpg]]
          if (is.null(model_obj) || inherits(model_obj, "dnaEPICO_methylationLME_fit_error")) {
            next
          }
          if (!is.null(model_obj$coef)) {
            phenotype_hits[[cpg]] <- as.data.frame(model_obj$coef)
          }
        }
      }

      retained[[phenotype]] <- phenotype_hits
      next
    }

    for (cpg in names(fit_list)) {
      model_obj <- fit_list[[cpg]]
      if (is.null(model_obj) || inherits(model_obj, "dnaEPICO_methylationLME_fit_error")) {
        next
      }

      coef_table <- model_obj$coef
      if (is.null(coef_table)) {
        next
      }

      matched_rows <- findCoefficientRowsMethylationLME(
        coefNames = rownames(coef_table),
        phenotype = phenotype,
        interactionTerm = interactionTerm
      )
      if (length(matched_rows) == 0L) {
        next
      }

      matched_pvals <- coef_table[matched_rows, "Pr(>|t|)", drop = TRUE]
      if (any(matched_pvals < threshold, na.rm = TRUE)) {
        phenotype_hits[[cpg]] <- as.data.frame(coef_table)
      }
    }

    retained[[phenotype]] <- phenotype_hits
  }

  hit_counts <- vapply(retained, length, integer(1))
  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      paste("Significant longitudinal terms retained at p <", threshold, ":"),
      paste(names(hit_counts), hit_counts, sep = ": ", collapse = "; "),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(retained, class = "dnaEPICO_methylationLME_significant")
}
#' Plot longitudinal mixed-effects model diagnostics
#'
#' @param modelSummaries Object returned by `summarizeMethylationLMEModels()`.
#' @param preparedData Object returned by `prepareMethylationLMEData()`.
#' @param fdrThreshold Numeric. False-discovery-rate threshold used to highlight
#'   CpGs in the diagnostic plots.
#' @param padjmethod Character. P-value adjustment method passed to
#'   `stats::p.adjust()`.
#' @param outputDir Character or `NULL`. Directory used for TIFF files. When
#'   `NULL`, plots are returned in memory only.
#' @param plotWidth Integer. TIFF width in pixels when plots are written to disk.
#' @param plotHeight Integer. TIFF height in pixels when plots are written to
#'   disk.
#' @param plotDPI Integer. TIFF resolution in DPI when plots are written to
#'   disk.
#' @param display Logical. If `TRUE`, draw plots on the active graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationLME_diagnostic_plots"`
#'   containing the generated `ggplot2` objects, genomic inflation factors, and
#'   any saved TIFF file paths.
#'
#' @description
#' Create Q-Q and standard-error diagnostic plots from the mixed-effects summary
#' tables returned by `summarizeMethylationLMEModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' diagnostic_plots <- plotMethylationLMEDiagnostics(
#'   modelSummaries = ex$modelSummaries,
#'   preparedData = ex$preparedData,
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(diagnostic_plots$plots)
#'
#' @export
plotMethylationLMEDiagnostics <- function(
    modelSummaries,
    preparedData,
    fdrThreshold = 0.05,
    padjmethod = "fdr",
    outputDir = NULL,
    plotWidth = 2000L,
    plotHeight = 1000L,
    plotDPI = 150L,
    display = FALSE,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  summary_list <- modelSummaries$summaries
  beta_matrix <- preparedData$data[, preparedData$cpgColumns, drop = FALSE]
  mean_beta <- colMeans(beta_matrix, na.rm = TRUE)
  plot_list <- list()
  inflation_factors <- list()
  saved_files <- list()

  for (phenotype in names(summary_list)) {
    summary_df <- summary_list[[phenotype]]
    if (is.null(summary_df) || nrow(summary_df) == 0L) {
      next
    }

    summary_df$FDR <- stats::p.adjust(summary_df$P.value, method = padjmethod)
    pvals <- summary_df$P.value
    pvals <- pvals[!is.na(pvals)]
    if (length(pvals) == 0L) {
      next
    }

    chisq <- stats::qchisq(1 - pvals, df = 1)
    lambda <- round(stats::median(chisq, na.rm = TRUE) / stats::qchisq(0.5, df = 1), 3)
    inflation_factors[[phenotype]] <- lambda

    qq_data <- data.frame(
      expected = -log10(stats::ppoints(length(pvals))),
      observed = -log10(sort(pvals))
    )
    qq_plot <- ggplot2::ggplot(qq_data, ggplot2::aes(x = expected, y = observed)) +
      ggplot2::geom_point(color = "black") +
      ggplot2::geom_abline(intercept = 0, slope = 1, color = "red") +
      ggplot2::labs(
        title = paste(
          "Q-Q Plot of p-values for",
          phenotype,
          "\nGenomic Inflation Factor =",
          lambda
        ),
        x = "Expected -log10(p)",
        y = "Observed -log10(p)"
      ) +
      ggplot2::theme_minimal()

    summary_df$log2meanBeta <- log2(pmax(mean_beta[summary_df$CpG], .Machine$double.xmin))
    residual_plot <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = log2meanBeta, y = Std.Error)
    ) +
      ggplot2::geom_point(alpha = 0.6, color = "black") +
      ggplot2::labs(
        title = "Standard Error vs Average Beta",
        x = "log2(Average Beta)",
        y = "Standard Error"
      ) +
      ggplot2::theme_minimal()

    if (nrow(summary_df) >= 3L) {
      residual_plot <- residual_plot +
        ggplot2::geom_smooth(
          method = "loess",
          formula = y ~ x,
          se = FALSE,
          color = "red"
        )
    }

    residual_significance_plot <- ggplot2::ggplot(
      summary_df,
      ggplot2::aes(x = -log10(P.value), y = Std.Error, color = FDR < fdrThreshold)
    ) +
      ggplot2::geom_point(alpha = 0.6) +
      ggrepel::geom_text_repel(
        data = subset(summary_df, FDR < fdrThreshold),
        ggplot2::aes(label = CpG)
      ) +
      ggplot2::scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick")) +
      ggplot2::labs(
        title = paste("Standard Error vs Significance for", phenotype),
        x = "-log10(p-value)",
        y = "Standard Error",
        color = paste("FDR <", fdrThreshold)
      ) +
      ggplot2::theme_minimal()

    qq_file <- NULL
    residual_file <- NULL
    significance_file <- NULL
    if (!is.null(outputDir)) {
      qq_file <- file.path(outputDir, paste0("qqplot_", phenotype, ".tiff"))
      residual_file <- file.path(outputDir, paste0("residualSD_", phenotype, ".tiff"))
      significance_file <- file.path(outputDir, paste0("residualSignificance_", phenotype, ".tiff"))
    }

    runPlotMinfiEwasWater(
      draw_fun = function() drawPlotObjectMinfiEwasWater(qq_plot),
      display = display,
      file = qq_file,
      width = plotWidth,
      height = plotHeight,
      res = plotDPI
    )
    runPlotMinfiEwasWater(
      draw_fun = function() drawPlotObjectMinfiEwasWater(residual_plot),
      display = display,
      file = residual_file,
      width = plotWidth,
      height = plotHeight,
      res = plotDPI
    )
    runPlotMinfiEwasWater(
      draw_fun = function() drawPlotObjectMinfiEwasWater(residual_significance_plot),
      display = display,
      file = significance_file,
      width = plotWidth,
      height = plotHeight,
      res = plotDPI
    )

    plot_list[[phenotype]] <- list(
      qqplot = qq_plot,
      residualSD = residual_plot,
      residualSignificance = residual_significance_plot
    )
    saved_files[[phenotype]] <- list(
      qqplot = qq_file,
      residualSD = residual_file,
      residualSignificance = significance_file
    )
  }

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      paste("Diagnostic plots generated for phenotypes:", length(plot_list)),
      if (is.null(outputDir)) {
        "Diagnostic plots were returned in memory only."
      } else {
        paste("Diagnostic plots saved to:    ", outputDir)
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      plots = plot_list,
      inflationFactors = inflation_factors,
      files = saved_files
    ),
    class = "dnaEPICO_methylationLME_diagnostic_plots"
  )
}

#' Annotate longitudinal mixed-effects summary tables with array annotation metadata
#'
#' @param modelSummaries Object returned by
#'   `summarizeMethylationLMEModels()` or a named list of summary data
#'   frames.
#' @param annotationObject Character package/object name, annotation data frame,
#'   or annotation object understood by `minfi::getAnnotation()`.
#' @param annotationCols Character vector or comma-separated string of annotation
#'   columns to append.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationLME_annotation"`
#'   containing the annotated summary table and any requested annotation columns
#'   that were unavailable in the chosen annotation object.
#'
#' @description
#' Merge phenotype-specific longitudinal summary tables with probe annotation
#' metadata and return a single annotated result table.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' annotation_data <- annotateMethylationLMESummaries(
#'   modelSummaries = ex$modelSummaries,
#'   annotationObject = ex$annotationData,
#'   annotationCols = "Name,chr,pos",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(annotation_data)
#'
#' @export
annotateMethylationLMESummaries <- function(
    modelSummaries,
    annotationObject,
    annotationCols = c(
      "Name",
      "chr",
      "pos",
      "UCSC_RefGene_Group",
      "UCSC_RefGene_Name",
      "Relation_to_Island",
      "GencodeV41_Group"
    ),
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  summary_list <- modelSummaries
  if (!is.null(modelSummaries$summaries)) {
    summary_list <- modelSummaries$summaries
  }

  annotation_cols <- splitOptionMinfiEwasWater(annotationCols, sep = ",")
  annotation_df <- coerceAnnotationDataMethylationGLM(annotationObject)

  cleaned_summaries <- lapply(
    names(summary_list),
    function(model_name) {
      df <- summary_list[[model_name]]
      if (is.null(df) || nrow(df) == 0L) {
        return(NULL)
      }
      if (!all(c("CpG", "Interaction.Term", "P.value") %in% colnames(df))) {
        return(NULL)
      }

      df_split <- split(df, df$Interaction.Term)
      model_tables <- lapply(
        names(df_split),
        function(term) {
          sub_df <- df_split[[term]][, c("CpG", "P.value"), drop = FALSE]
          clean_term <- gsub("`", "", term, fixed = TRUE)
          interaction_suffix <- clean_term
          if (startsWith(clean_term, paste0(model_name, "."))) {
            interaction_suffix <- sub(
              paste0("^", escapeRegexMethylationGLM(model_name), "\\."),
              "",
              clean_term
            )
          }
          p_col <- paste0(model_name, "_", interaction_suffix, "_P.Value")
          colnames(sub_df)[2] <- p_col
          sub_df
        }
      )

      Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), model_tables)
    }
  )
  cleaned_summaries <- Filter(Negate(is.null), cleaned_summaries)

  if (length(cleaned_summaries) == 0L) {
    merged_summary <- data.frame(CpG = character(0))
  } else if (length(cleaned_summaries) == 1L) {
    merged_summary <- cleaned_summaries[[1L]]
  } else {
    merged_summary <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), cleaned_summaries)
  }

  available_annotation_cols <- intersect(annotation_cols, colnames(annotation_df))
  missing_annotation_cols <- setdiff(annotation_cols, colnames(annotation_df))
  annotated_results <- merge(
    merged_summary,
    annotation_df[, c("CpG", available_annotation_cols), drop = FALSE],
    by = "CpG",
    all.x = TRUE
  )
  if ("CpG" %in% colnames(annotated_results)) {
    colnames(annotated_results)[colnames(annotated_results) == "CpG"] <- "IlmnID"
  }

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      paste("Annotated CpG rows:          ", nrow(annotated_results)),
      paste("Annotation columns used:      ", paste(available_annotation_cols, collapse = ", ")),
      if (length(missing_annotation_cols) == 0L) {
        "Missing annotation columns:   none"
      } else {
        paste("Missing annotation columns:   ", paste(missing_annotation_cols, collapse = ", "))
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      data = annotated_results,
      annotationColumnsUsed = available_annotation_cols,
      missingAnnotationCols = missing_annotation_cols
    ),
    class = "dnaEPICO_methylationLME_annotation"
  )
}

#' Write optional disk outputs for longitudinal mixed-effects analyses
#'
#' @param modelResults Object returned by `fitMethylationLMEModels()`.
#' @param modelSummaries Object returned by `summarizeMethylationLMEModels()`.
#' @param annotatedResults Object returned by
#'   `annotateMethylationLMESummaries()` or a compatible data frame.
#' @param significantInteractions Object returned by
#'   `collectSignificantInteractionsMethylationLME()` or `NULL`.
#' @param outputRData Character. Directory used for serialized model and summary
#'   outputs.
#' @param summaryTxtDir Character. Directory used for tab-delimited summary
#'   tables.
#' @param significantInteractionDir Character. Directory used for significant
#'   interaction coefficient tables.
#' @param annotatedLMEOut Character. Directory used for the annotated summary
#'   XLSX workbook.
#' @param saveTxtSummaries Logical. If `TRUE`, write tab-delimited summary tables.
#' @param saveSignificantInteractions Logical. If `TRUE`, write significant
#'   interaction coefficient tables.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationLME_paths"`
#'   containing the paths of the files written to disk.
#'
#' @description
#' Write optional serialized outputs, summary tables, significant interaction
#' tables, and annotated results from the longitudinal mixed-effects workflow.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationLMEStateDnaEpico()
#' annotation_data <- annotateMethylationLMESummaries(
#'   modelSummaries = ex$modelSummaries,
#'   annotationObject = ex$annotationData,
#'   annotationCols = "Name,chr,pos",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' significant_hits <- collectSignificantInteractionsMethylationLME(
#'   modelResults = ex$modelResults,
#'   pvalThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' output_paths <- writeMethylationLMEOutputs(
#'   modelResults = ex$modelResults,
#'   modelSummaries = ex$modelSummaries,
#'   annotatedResults = annotation_data,
#'   significantInteractions = significant_hits,
#'   outputRData = file.path(ex$tempDir, "models"),
#'   summaryTxtDir = file.path(ex$tempDir, "summary"),
#'   significantInteractionDir = file.path(ex$tempDir, "significant"),
#'   annotatedLMEOut = file.path(ex$tempDir, "annotated"),
#'   saveTxtSummaries = TRUE,
#'   saveSignificantInteractions = TRUE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(output_paths)
#'
#' @export
writeMethylationLMEOutputs <- function(
    modelResults,
    modelSummaries,
    annotatedResults,
    significantInteractions = NULL,
    outputRData,
    summaryTxtDir,
    significantInteractionDir,
    annotatedLMEOut,
    saveTxtSummaries = TRUE,
    saveSignificantInteractions = FALSE,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationLME.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  dir.create(outputRData, recursive = TRUE, showWarnings = FALSE)
  dir.create(annotatedLMEOut, recursive = TRUE, showWarnings = FALSE)

  model_files <- stats::setNames(character(length(modelResults$fits)), names(modelResults$fits))
  summary_files <- stats::setNames(character(length(modelSummaries$summaries)), names(modelSummaries$summaries))
  summary_txt_files <- list()
  significant_files <- list()

  for (phenotype in names(modelResults$fits)) {
    model_file <- file.path(outputRData, paste0(phenotype, "LME.rds"))
    summary_file <- file.path(outputRData, paste0(phenotype, "SummaryLME.rds"))
    saveRDS(modelResults$fits[[phenotype]], file = model_file)
    saveRDS(modelSummaries$summaries[[phenotype]], file = summary_file)
    model_files[[phenotype]] <- model_file
    summary_files[[phenotype]] <- summary_file
  }

  if (isTRUE(saveTxtSummaries)) {
    dir.create(summaryTxtDir, recursive = TRUE, showWarnings = FALSE)
    for (phenotype in names(modelSummaries$summaries)) {
      summary_df <- modelSummaries$summaries[[phenotype]]
      if (is.null(summary_df) || nrow(summary_df) == 0L) {
        next
      }

      if ("P.value" %in% colnames(summary_df)) {
        summary_df <- summary_df[order(summary_df$P.value), , drop = FALSE]
      }

      output_file <- file.path(summaryTxtDir, paste0(phenotype, "SummaryLME.txt"))
      utils::write.table(summary_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      summary_txt_files[[phenotype]] <- output_file
    }
  }

  if (isTRUE(saveSignificantInteractions) && !is.null(significantInteractions)) {
    dir.create(significantInteractionDir, recursive = TRUE, showWarnings = FALSE)
    for (phenotype in names(significantInteractions)) {
      phenotype_hits <- significantInteractions[[phenotype]]
      if (length(phenotype_hits) == 0L) {
        next
      }

      significant_files[[phenotype]] <- character(0)
      phenotype_dir <- file.path(significantInteractionDir, phenotype)
      dir.create(phenotype_dir, recursive = TRUE, showWarnings = FALSE)

      for (cpg in names(phenotype_hits)) {
        cpg_dir <- file.path(phenotype_dir, cpg)
        dir.create(cpg_dir, recursive = TRUE, showWarnings = FALSE)
        output_file <- file.path(cpg_dir, paste0(cpg, ".txt"))
        utils::write.table(phenotype_hits[[cpg]], file = output_file, sep = "\t", quote = FALSE)
        significant_files[[phenotype]] <- c(significant_files[[phenotype]], output_file)
      }
    }
  }

  annotated_df <- annotatedResults
  if (!is.null(annotatedResults$data)) {
    annotated_df <- annotatedResults$data
  }
  annotated_file <- file.path(annotatedLMEOut, "annotatedLME.xlsx")
  dictionary <- buildAnnotatedWorkbookDictionaryMethylationGLM(
    columns = colnames(annotated_df),
    modelDescription = "Pvalue from LME model",
    formulaText = modelResults$formulas,
    modelLabel = "LME",
    responseLabel = inferMethylationValueLabelMethylationGLM(modelResults)
  )
  lme_engine <- modelResults$settings$lmeEngine
  if (is.null(lme_engine) || !nzchar(lme_engine)) {
    requested_libraries <- tolower(as.character(modelResults$settings$lmeLibs))
    lme_engine <- if ("nlme" %in% requested_libraries) "nlme" else "lme4"
  }
  correlation_structure <- modelResults$settings$correlationStructure
  if (is.null(correlation_structure) || !nzchar(correlation_structure)) {
    correlation_structure <- "none"
  }
  correlation_var <- modelResults$settings$correlationVar
  if (is.null(correlation_var) || !nzchar(correlation_var)) {
    correlation_var <- "None"
  }
  interaction_term <- modelResults$settings$interactionTerm
  if (is.null(interaction_term) || !nzchar(interaction_term)) {
    interaction_term <- "None"
  }
  metadata <- data.frame(
    Key = c(
      "analysis",
      "backend",
      "fitting_function",
      "libraries",
      "correlation_structure",
      "correlation_variable",
      "interaction_term",
      "response_label",
      "phenotypes",
      "cpg_count",
      "annotation_columns",
      "missing_annotation_columns"
    ),
    Value = c(
      "methylationLME",
      lme_engine,
      if (identical(lme_engine, "nlme")) "nlme::lme" else "lmerTest::lmer",
      paste(modelResults$settings$lmeLibs, collapse = ","),
      correlation_structure,
      correlation_var,
      interaction_term,
      inferMethylationValueLabelMethylationGLM(modelResults),
      paste(modelResults$phenotypes, collapse = ","),
      as.character(nrow(annotated_df)),
      if (!is.null(annotatedResults$annotationColumnsUsed)) {
        paste(annotatedResults$annotationColumnsUsed, collapse = ",")
      } else {
        ""
      },
      if (!is.null(annotatedResults$missingAnnotationCols)) {
        paste(annotatedResults$missingAnnotationCols, collapse = ",")
      } else {
        ""
      }
    ),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  writeAnnotatedWorkbookMethylationGLM(
    annotated_df = annotated_df,
    file = annotated_file,
    resultSheet = "annotatedLME",
    dictionary = dictionary,
    metadata = metadata
  )

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      paste("Serialized model outputs:    ", length(model_files)),
      paste("Serialized summary outputs:  ", length(summary_files)),
      paste("Annotated results file:      ", annotated_file),
      if (isTRUE(saveTxtSummaries)) {
        paste("Summary text files written:   ", length(summary_txt_files))
      } else {
        "Summary text files written:   0"
      },
      if (isTRUE(saveSignificantInteractions)) {
        paste("Significant interaction files:", sum(vapply(significant_files, length, integer(1))))
      } else {
        "Significant interaction files: 0"
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      modelFiles = model_files,
      summaryFiles = summary_files,
      summaryTxtFiles = summary_txt_files,
      significantInteractionFiles = significant_files,
      annotatedLME = annotated_file
    ),
    class = "dnaEPICO_methylationLME_paths"
  )
}
