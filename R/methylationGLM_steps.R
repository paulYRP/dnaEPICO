#' Normalize optional numeric inputs for methylationGLM helpers
#'
#' @param value Numeric, character, or `NA` input.
#'
#' @return A numeric scalar or `NA_real_`.
#'
#' @description
#' Internal helper that standardizes optional numeric arguments such as
#' `cpgLimit` and `summaryPval`.
#'
#' @keywords internal
#' @noRd
normalizeOptionalNumericMethylationGLM <- function(value) {
  if (is.null(value) || length(value) == 0L) {
    return(NA_real_)
  }

  if (length(value) == 1L && is.character(value)) {
    trimmed <- trimws(value)
    if (!nzchar(trimmed) || tolower(trimmed) %in% c("na", "null")) {
      return(NA_real_)
    }
    return(as.numeric(trimmed))
  }

  as.numeric(value[[1L]])
}

#' Normalize optional chunk-size inputs for methylationGLM helpers
#'
#' @param chunkSize Integer, character, or `NULL` input.
#'
#' @return An integer scalar or `NULL`.
#'
#' @description
#' Internal helper that standardizes optional chunk-size arguments used while
#' extracting CpG-level summaries.
#'
#' @keywords internal
#' @noRd
normalizeChunkSizeMethylationGLM <- function(chunkSize) {
  if (is.null(chunkSize) || length(chunkSize) == 0L) {
    return(NULL)
  }

  if (length(chunkSize) == 1L && is.character(chunkSize)) {
    trimmed <- trimws(chunkSize)
    if (!nzchar(trimmed) || tolower(trimmed) %in% c("na", "null")) {
      return(NULL)
    }
    return(as.integer(trimmed))
  }

  as.integer(chunkSize[[1L]])
}

#' Escape regular-expression metacharacters for methylationGLM helpers
#'
#' @param x Character scalar to escape.
#'
#' @return Character scalar safe to use in `grep()` or `grepl()`.
#'
#' @description
#' Internal helper that escapes user-supplied model terms before matching
#' coefficient names.
#'
#' @keywords internal
#' @noRd
escapeRegexMethylationGLM <- function(x) {
  gsub("([][{}()+*^$|\\?.])", "\\\\\\1", x)
}

#' Backtick variable names for methylationGLM formulas
#'
#' @param x Character vector of variable names.
#'
#' @return Character vector with names wrapped in backticks.
#'
#' @description
#' Internal helper that keeps formula construction robust when variable names
#' contain punctuation.
#'
#' @keywords internal
#' @noRd
quoteNamesMethylationGLM <- function(x) {
  paste0("`", gsub("`", "", x, fixed = TRUE), "`")
}

#' Parse phenotype-to-PRS mappings for methylationGLM helpers
#'
#' @param prsMap Character scalar or vector describing mappings in the form
#'   `"Phenotype:PRS"`.
#'
#' @return A named character vector.
#'
#' @description
#' Internal helper that parses the optional `prsMap` argument.
#'
#' @keywords internal
#' @noRd
parsePrsMapMethylationGLM <- function(prsMap = NULL) {
  if (is.null(prsMap) || length(prsMap) == 0L) {
    return(stats::setNames(character(0), character(0)))
  }

  map_entries <- splitOptionMinfiEwasWater(prsMap, sep = ",")
  if (length(map_entries) == 0L) {
    return(stats::setNames(character(0), character(0)))
  }

  pieces <- strsplit(map_entries, split = ":", fixed = TRUE)
  invalid <- !vapply(pieces, function(x) length(x) >= 2L, logical(1))
  if (any(invalid)) {
    stop(
      "Each prsMap entry must follow the format 'Phenotype:PRS'. Invalid entries: ",
      paste(map_entries[invalid], collapse = ", "),
      call. = FALSE
    )
  }

  keys <- vapply(pieces, function(x) trimws(x[[1L]]), character(1))
  values <- vapply(pieces, function(x) trimws(x[[2L]]), character(1))
  stats::setNames(values, keys)
}

#' Find coefficient rows for a phenotype or interaction term
#'
#' @param coefNames Character vector of model coefficient names.
#' @param variable Character. Phenotype variable of interest.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return Character vector of matching coefficient names.
#'
#' @description
#' Internal helper that locates the phenotype main-effect or interaction rows in
#' a CpG-level coefficient table.
#'
#' @keywords internal
#' @noRd
findCoefficientRowsMethylationGLM <- function(
    coefNames,
    variable,
    interactionTerm = NULL
) {
  normalized_names <- gsub("`", "", coefNames, fixed = TRUE)
  variable_pattern <- escapeRegexMethylationGLM(variable)

  if (!is.null(interactionTerm) && nzchar(interactionTerm)) {
    interaction_pattern <- escapeRegexMethylationGLM(interactionTerm)
    matches <- grepl(
      paste0("^", variable_pattern, ".*:", interaction_pattern),
      normalized_names
    )

    if (!any(matches)) {
      matches <- grepl(paste0("^", variable_pattern), normalized_names)
    }
  } else {
    matches <- grepl(paste0("^", variable_pattern), normalized_names)
  }

  coefNames[matches]
}

#' Build a GLM formula string for methylationGLM helpers
#'
#' @param phenotype Character. Phenotype variable of interest.
#' @param covariates Character vector of covariate variables.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#'
#' @return Character scalar containing a formula string.
#'
#' @description
#' Internal helper that builds the one-CpG Gaussian GLM formula.
#'
#' @keywords internal
#' @noRd
buildFormulaMethylationGLM <- function(
    phenotype,
    covariates = character(0),
    interactionTerm = NULL,
    responseVar = "beta"
) {
  quoted_phenotype <- quoteNamesMethylationGLM(phenotype)
  quoted_covariates <- quoteNamesMethylationGLM(covariates)

  if (!is.null(interactionTerm) && nzchar(interactionTerm)) {
    quoted_interaction <- quoteNamesMethylationGLM(interactionTerm)
    interaction_part <- paste(quoted_phenotype, quoted_interaction, sep = " * ")
    fixed_terms <- setdiff(covariates, interactionTerm)
    quoted_fixed_terms <- quoteNamesMethylationGLM(fixed_terms)
    terms <- c(interaction_part, quoted_fixed_terms)
  } else {
    terms <- c(quoted_phenotype, quoted_covariates)
  }

  terms <- unique(terms[nzchar(terms)])
  if (length(terms) == 0L) {
    stop("At least one phenotype or covariate term is required.", call. = FALSE)
  }

  paste(responseVar, "~", paste(terms, collapse = " + "))
}

#' Fit a single CpG-level Gaussian GLM for methylationGLM helpers
#'
#' @param cpg Character. CpG column name.
#' @param data Data frame containing phenotype and beta columns.
#' @param modelVars Character vector of variables to retain from `data`.
#' @param formulaText Character scalar containing the model formula.
#' @param factorVars Character vector of variables to coerce to factors.
#'
#' @return A list with coefficient, fitted-value, and residual components, or an
#'   error object when fitting fails.
#'
#' @description
#' Internal helper that fits one Gaussian GLM for a single CpG site.
#'
#' @keywords internal
#' @noRd
fitCpGModelMethylationGLM <- function(
    cpg,
    data,
    modelVars,
    formulaText,
    factorVars = character(0),
    responseVar = "beta"
) {
  tryCatch(
    {
      model_data <- data[, modelVars, drop = FALSE]
      model_data[[responseVar]] <- as.numeric(data[[cpg]])

      for (var in intersect(factorVars, colnames(model_data))) {
        model_data[[var]] <- as.factor(model_data[[var]])
      }

      fit <- glm2::glm2(
        formula = stats::as.formula(formulaText),
        data = model_data,
        family = stats::gaussian(),
        na.action = stats::na.exclude
      )

      list(
        coef = summary(fit)$coefficients,
        residuals = stats::residuals(fit),
        fitted = stats::fitted(fit)
      )
    },
    error = function(error) {
      structure(
        list(error = conditionMessage(error)),
        class = "dnaEPICO_methylationGLM_fit_error"
      )
    }
  )
}

#' Summarize a single CpG-level GLM fit
#'
#' @param cpg Character. CpG identifier.
#' @param modelObj List returned by `fitCpGModelMethylationGLM()`.
#' @param variable Character. Phenotype variable of interest.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#' @param includeResidualSD Logical. If `TRUE`, append residual standard
#'   deviations to the output.
#'
#' @return A data frame for the requested CpG, or `NULL`.
#'
#' @description
#' Internal helper that extracts phenotype-specific coefficient rows from one
#' CpG-level GLM fit.
#'
#' @keywords internal
#' @noRd
summarizeCpGFitMethylationGLM <- function(
    cpg,
    modelObj,
    variable,
    interactionTerm = NULL,
    includeResidualSD = TRUE
) {
  if (is.null(modelObj) || inherits(modelObj, "dnaEPICO_methylationGLM_fit_error")) {
    return(NULL)
  }

  coef_table <- modelObj$coef
  if (is.null(coef_table)) {
    return(NULL)
  }

  matched_rows <- findCoefficientRowsMethylationGLM(
    coefNames = rownames(coef_table),
    variable = variable,
    interactionTerm = interactionTerm
  )
  if (length(matched_rows) == 0L) {
    return(NULL)
  }

  summary_df <- as.data.frame(coef_table[matched_rows, , drop = FALSE])
  summary_df$CpG <- cpg
  summary_df$Coefficient <- rownames(summary_df)

  if (isTRUE(includeResidualSD) && !is.null(modelObj$residuals)) {
    summary_df$ResidualSD <- stats::sd(modelObj$residuals, na.rm = TRUE)
  }

  summary_df
}

resolveParallelBackendMethylationModels <- function(nCores) {
  n_cores <- max(1L, as.integer(nCores))
  requested_backend <- tolower(Sys.getenv("DNAEPICO_PARALLEL_BACKEND", "auto"))
  if (!(requested_backend %in% c("auto", "fork", "psock", "serial"))) {
    requested_backend <- "auto"
  }

  if (n_cores <= 1L || identical(requested_backend, "serial")) {
    return("serial")
  }

  fork_supported <- !identical(.Platform$OS.type, "windows")
  if (identical(requested_backend, "fork")) {
    return(if (fork_supported) "fork" else "psock")
  }
  if (identical(requested_backend, "psock")) {
    return("psock")
  }

  if (fork_supported) {
    return("fork")
  }

  "psock"
}

chunkCpGColumnsMethylationModels <- function(
    cpgColumns,
    nCores = 1L,
    batchesPerCore = 8L
) {
  if (length(cpgColumns) == 0L) {
    return(list())
  }

  n_cores <- max(1L, as.integer(nCores))
  batches_per_core <- max(1L, as.integer(batchesPerCore))
  target_batches <- min(length(cpgColumns), max(1L, n_cores * batches_per_core))
  chunk_size <- max(1L, ceiling(length(cpgColumns) / target_batches))

  split(cpgColumns, ceiling(seq_along(cpgColumns) / chunk_size))
}

validateWorkerPackagesMethylationModels <- function(
    libPath = NULL,
    packages = character(0)
) {
  if (!is.null(libPath)) {
    .libPaths(unique(c(libPath, .libPaths())))
  }

  packages <- packages[!is.na(packages) & nzchar(packages)]
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Failed to load package: ", pkg, call. = FALSE)
    }
  }

  invisible(TRUE)
}

makePsockClusterMethylationModels <- function(clusterSize) {
  tryCatch(
    parallel::makeCluster(
      clusterSize,
      type = "PSOCK",
      useXDR = FALSE,
      methods = FALSE
    ),
    error = function(error) {
      parallel::makeCluster(clusterSize, type = "PSOCK")
    }
  )
}

combineFitBatchResultsMethylationModels <- function(
    batchResults,
    cpgColumns
) {
  fit_chunks <- lapply(batchResults, function(result) result$fits)
  fit_chunks <- Filter(function(x) !is.null(x) && length(x) > 0L, fit_chunks)
  fit_list <- if (length(fit_chunks) == 0L) {
    list()
  } else {
    out <- list()
    for (chunk in fit_chunks) {
      for (cpg in names(chunk)) {
        out[[cpg]] <- chunk[[cpg]]
      }
    }
    out
  }
  fit_list <- fit_list[cpgColumns[cpgColumns %in% names(fit_list)]]

  summary_chunks <- lapply(batchResults, function(result) result$summaries)
  summary_chunks <- Filter(function(x) !is.null(x) && nrow(x) > 0L, summary_chunks)
  summary_df <- if (length(summary_chunks) == 0L) {
    data.frame()
  } else {
    out <- do.call(rbind, summary_chunks)
    rownames(out) <- NULL
    out
  }

  list(fits = fit_list, summaries = summary_df)
}

optionalTermMatchesMethylationModels <- function(
    requested,
    cached
) {
  normalize <- function(value) {
    if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) {
      return("")
    }
    value <- as.character(value[[1L]])
    if (!nzchar(value)) {
      return("")
    }
    value
  }

  identical(normalize(requested), normalize(cached))
}

filterSummaryByPvalueMethylationGLM <- function(
    summaryDf,
    pValueFilter,
    includeResidualSD = TRUE
) {
  summary_df <- summaryDf
  if (is.null(summary_df) || nrow(summary_df) == 0L) {
    return(data.frame())
  }

  ordered_columns <- c(
    "CpG",
    "Coefficient",
    "Estimate",
    "Std. Error",
    "t value",
    "Pr(>|t|)",
    if (isTRUE(includeResidualSD)) "ResidualSD"
  )
  ordered_columns <- intersect(ordered_columns, colnames(summary_df))
  summary_df <- summary_df[, ordered_columns, drop = FALSE]

  if (nrow(summary_df) > 0L && !is.na(pValueFilter)) {
    summary_df <- summary_df[summary_df[["Pr(>|t|)"]] < pValueFilter, , drop = FALSE]
  }
  rownames(summary_df) <- NULL

  summary_df
}

fitCpGModelMethylationGLMPrepared <- function(
    cpg,
    cpgValues,
    modelData,
    formulaText,
    responseVar = "beta"
) {
  tryCatch(
    {
      model_data <- modelData
      model_data[[responseVar]] <- as.numeric(cpgValues)

      fit <- glm2::glm2(
        formula = stats::as.formula(formulaText),
        data = model_data,
        family = stats::gaussian(),
        na.action = stats::na.exclude
      )

      list(
        coef = summary(fit)$coefficients,
        residuals = stats::residuals(fit),
        fitted = stats::fitted(fit)
      )
    },
    error = function(error) {
      structure(
        list(error = conditionMessage(error)),
        class = "dnaEPICO_methylationGLM_fit_error"
      )
    }
  )
}

fitMethylationGLMBatch <- function(
    cpgBatch,
    data,
    modelData,
    formulaText,
    phenotype,
    interactionTerm = NULL,
    responseVar = "beta"
) {
  fits <- vector("list", length(cpgBatch))
  names(fits) <- cpgBatch
  summaries <- vector("list", length(cpgBatch))
  names(summaries) <- cpgBatch

  for (cpg in cpgBatch) {
    model_obj <- fitCpGModelMethylationGLMPrepared(
      cpg = cpg,
      cpgValues = data[[cpg]],
      modelData = modelData,
      formulaText = formulaText,
      responseVar = responseVar
    )
    fits[[cpg]] <- model_obj
    summaries[[cpg]] <- summarizeCpGFitMethylationGLM(
      cpg = cpg,
      modelObj = model_obj,
      variable = phenotype,
      interactionTerm = interactionTerm,
      includeResidualSD = TRUE
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

#' Create a distribution plot used by methylationGLM helpers
#'
#' @param values Vector of phenotype or covariate values to plot.
#' @param variable Character. Variable label.
#' @param type Character. Either `"hist"` or `"bar"`.
#' @param fill Character. Fill colour used for the bars.
#'
#' @return A `ggplot2` object.
#'
#' @description
#' Internal helper that builds exploratory distribution plots.
#'
#' @keywords internal
#' @noRd
createDistributionPlotMethylationGLM <- function(
    values,
    variable,
    type = c("hist", "bar"),
    fill = "steelblue"
) {
  type <- match.arg(type)

  if (identical(type, "hist")) {
    plot_data <- data.frame(value = as.numeric(values))
    return(
      ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
        ggplot2::geom_histogram(bins = 30, fill = fill, color = "white") +
        ggplot2::labs(
          title = paste("Distribution of", variable),
          x = variable,
          y = "Frequency"
        ) +
        ggplot2::theme_minimal()
    )
  }

  plot_data <- data.frame(value = as.factor(values))
  ggplot2::ggplot(plot_data, ggplot2::aes(x = value)) +
    ggplot2::geom_bar(fill = fill) +
    ggplot2::labs(
      title = paste("Distribution of", variable),
      x = variable,
      y = "Count"
    ) +
    ggplot2::theme_minimal()
}

#' Resolve an annotation object for methylationGLM helpers
#'
#' @param annotationObject Character package/object name, annotation data frame,
#'   or annotation object understood by `minfi::getAnnotation()`.
#'
#' @return An object accepted by `minfi::getAnnotation()`.
#'
#' @description
#' Internal helper that loads or validates the annotation resource used to
#' annotate CpG-level summary tables.
#'
#' @keywords internal
#' @noRd
resolveAnnotationObjectMethylationGLM <- function(annotationObject) {
  if (!is.character(annotationObject) || length(annotationObject) != 1L) {
    return(annotationObject)
  }

  if (requireNamespace(annotationObject, quietly = TRUE)) {
    annotation_lookup <- suppressPackageStartupMessages(
      tryCatch(
        minfi::getAnnotation(annotationObject),
        error = function(e) NULL
      )
    )
    if (!is.null(annotation_lookup)) {
      return(annotationObject)
    }

    annotation_namespace <- asNamespace(annotationObject)
    if (exists(annotationObject, envir = annotation_namespace, inherits = FALSE)) {
      return(get(annotationObject, envir = annotation_namespace, inherits = FALSE))
    }
  }

  loaded_namespaces <- loadedNamespaces()
  for (namespace_name in loaded_namespaces) {
    namespace_env <- asNamespace(namespace_name)
    if (exists(annotationObject, envir = namespace_env, inherits = FALSE)) {
      return(get(annotationObject, envir = namespace_env, inherits = FALSE))
    }
  }

  if (exists(annotationObject, inherits = TRUE)) {
    return(get(annotationObject, inherits = TRUE))
  }

  stop(
    "Annotation package or object was not found: ",
    annotationObject,
    call. = FALSE
  )
}

coerceAnnotationDataMethylationGLM <- function(annotationObject) {
  if (is.data.frame(annotationObject)) {
    annotation_df <- annotationObject
    if (!("CpG" %in% colnames(annotation_df))) {
      if ("IlmnID" %in% colnames(annotation_df)) {
        annotation_df$CpG <- annotation_df$IlmnID
        return(annotation_df)
      }
      if (is.null(rownames(annotation_df))) {
        stop(
          "annotationObject data frames must include a CpG column or row names.",
          call. = FALSE
        )
      }
      annotation_df$CpG <- rownames(annotation_df)
    }
    return(annotation_df)
  }

  annotation_source <- resolveAnnotationObjectMethylationGLM(annotationObject)
  annotation_df <- as.data.frame(minfi::getAnnotation(annotation_source))
  annotation_df$CpG <- rownames(annotation_df)
  annotation_df
}

buildAnnotatedWorkbookDictionaryMethylationGLM <- function(
    columns,
    modelDescription,
    formulaText,
    modelLabel,
    responseLabel
) {
  pvalue_columns <- grepl("P\\.Value$|P\\.value$", columns)
  default_formula <- paste(responseLabel, "~ formula unavailable")
  formula_values <- as.character(formulaText)
  formula_values <- formula_values[!is.na(formula_values) & nzchar(formula_values)]

  select_formula <- function(column) {
    if (length(formula_values) == 0L) {
      return(default_formula)
    }

    formula_names <- names(formula_values)
    if (!is.null(formula_names) && any(nzchar(formula_names))) {
      for (formula_name in formula_names[nzchar(formula_names)]) {
        if (grepl(paste0("^", escapeRegexMethylationGLM(formula_name)), column)) {
          return(unname(formula_values[[formula_name]]))
        }
      }
    }

    if (length(formula_values) == 1L) {
      return(unname(formula_values[[1L]]))
    }

    paste(unique(unname(formula_values)), collapse = "; ")
  }

  display_formula <- function(column) {
    formula_text <- select_formula(column)
    formula_text <- sub(
      "^\\s*[^~]+\\s*~",
      paste(responseLabel, "~"),
      formula_text,
      ignore.case = TRUE
    )
    paste0(modelLabel, ": ", formula_text)
  }

  descriptions <- ifelse(
    pvalue_columns,
    modelDescription,
    ifelse(
      columns %in% c("IlmnID", "CpG", "Name"),
      "CpG probe identifier",
      "Genomic annotation or supporting result column"
    )
  )
  formulas <- rep("", length(columns))
  formulas[pvalue_columns] <- vapply(
    columns[pvalue_columns],
    display_formula,
    character(1)
  )

  data.frame(
    Column = columns,
    Description = descriptions,
    Formula = formulas,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

writeAnnotatedWorkbookMethylationGLM <- function(
    annotated_df,
    file,
    resultSheet,
    dictionary
) {
  workbook <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(workbook, resultSheet)
  openxlsx::writeData(workbook, sheet = resultSheet, x = annotated_df)
  openxlsx::addWorksheet(workbook, "dictionary")
  openxlsx::writeData(workbook, sheet = "dictionary", x = dictionary)
  openxlsx::saveWorkbook(workbook, file = file, overwrite = TRUE)

  invisible(file)
}

inferMethylationValueLabelMethylationGLM <- function(modelResults) {
  if (!is.null(modelResults$responseLabel) && nzchar(modelResults$responseLabel)) {
    return(modelResults$responseLabel)
  }

  if (!is.null(modelResults$settings$methylationScale)) {
    return(methylationScaleResponseLabelDnaEpico(modelResults$settings$methylationScale))
  }

  if (is.null(modelResults$fits) || length(modelResults$fits) == 0L) {
    return("Beta values")
  }

  found_response <- FALSE
  for (fit_group in modelResults$fits) {
    if (!is.list(fit_group)) {
      next
    }

    for (fit_object in fit_group) {
      if (
        !is.list(fit_object) ||
          inherits(fit_object, "dnaEPICO_methylationGLM_fit_error") ||
          inherits(fit_object, "dnaEPICO_methylationLME_fit_error") ||
          is.null(fit_object$fitted) ||
          is.null(fit_object$residuals)
      ) {
        next
      }

      fitted_values <- fit_object$fitted
      residual_values <- fit_object$residuals
      if (
        !is.numeric(fitted_values) ||
          !is.numeric(residual_values) ||
          length(fitted_values) != length(residual_values)
      ) {
        next
      }

      response_values <- fitted_values + residual_values
      response_values <- response_values[is.finite(response_values)]
      if (length(response_values) == 0L) {
        next
      }

      found_response <- TRUE
      if (any(response_values < 0 | response_values > 1, na.rm = TRUE)) {
        return("M-values")
      }
    }
  }

  if (isTRUE(found_response)) {
    return("Beta values")
  }

  "Beta values"
}

#' Prepare phenotype-plus-methylation data for one-timepoint GLM analyses
#'
#' @param inputPheno Character. Path to the merged phenotype-plus-methylation object
#'   created by `preprocessingPheno()`.
#' @param phenotypes Character vector or comma-separated string of phenotype
#'   variables to model.
#' @param covariates Character vector or comma-separated string of covariate
#'   variables to adjust for.
#' @param factorVars Character vector or comma-separated string of variables that
#'   should be converted to factors before modeling.
#' @param cpgPrefix Character. Prefix used to identify methylation columns.
#' @param cpgLimit Integer or `NA`. Maximum number of CpGs to retain. `NA`
#'   keeps all matching CpGs.
#' @param methylationScale Character. Methylation metric represented by the CpG
#'   columns. One of `"beta"`, `"m"`, or `"cn"`.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#' @param prsMap Character vector or comma-separated string of phenotype-to-PRS
#'   mappings in the form `"Phenotype:PRS"`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationGLM_data"` containing the
#'   prepared analysis data, parsed variable selections, CpG columns, and
#'   exploratory summaries.
#'
#' @description
#' Load the merged phenotype-plus-methylation input object, validate the requested
#' modeling variables, convert selected variables to factors, and return a
#' single in-memory object for downstream helpers.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' prepared_data <- prepareMethylationGLMData(
#'   inputPheno = ex$inputPath,
#'   phenotypes = "status",
#'   covariates = "sex,age",
#'   factorVars = "status,sex",
#'   cpgLimit = 2,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(prepared_data)
#'
#' @export
prepareMethylationGLMData <- function(
    inputPheno,
    phenotypes,
    covariates,
    factorVars,
    cpgPrefix = "cg",
    cpgLimit = NA,
    methylationScale = "beta",
    interactionTerm = NULL,
    prsMap = NULL,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationGLM.txt"
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
  analysis_data <- loadSavedObjectPreprocessingPheno(inputPheno, preferred_name = "phenoBT1")

  if (!is.data.frame(analysis_data)) {
    analysis_data <- as.data.frame(analysis_data, stringsAsFactors = FALSE)
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

  for (var in intersect(factor_list, colnames(analysis_data))) {
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

  requested_columns <- unique(c(phenotype_list, covariate_list))
  missing_counts <- vapply(
    requested_columns,
    function(column_name) sum(is.na(analysis_data[[column_name]])),
    integer(1)
  )
  variable_summary <- summary(analysis_data[, requested_columns, drop = FALSE])
  interaction_table <- NULL
  if (!is.null(resolved_interaction) && nzchar(resolved_interaction)) {
    interaction_table <- table(analysis_data[[resolved_interaction]], useNA = "ifany")
  }

  log_lines <- c(
    "=======================================================================",
    paste("Loaded phenotype + methylation data from:", inputPheno),
    paste("Merged modeling object:     ", methylation_prefix, "*"),
    paste("Data dimensions:             ", paste(dim(analysis_data), collapse = " x ")),
    paste("Phenotypes:                  ", paste(phenotype_list, collapse = ", ")),
    paste("Covariates:                  ", paste(covariate_list, collapse = ", ")),
    paste("Factor variables:            ", paste(factor_list, collapse = ", ")),
    paste("CpG columns retained:        ", length(cpg_columns)),
    "Missing summary:",
    paste(names(missing_counts), missing_counts, sep = ": ", collapse = "; "),
    "Summary statistics:",
    previewLinesMinfiEwasWater(variable_summary)
  )
  if (!is.null(interaction_table)) {
    log_lines <- c(
      log_lines,
      paste("Interaction table for", resolved_interaction, ":"),
      previewLinesMinfiEwasWater(interaction_table)
    )
  }
  log_lines <- c(log_lines, "=======================================================================")
  emitLogMinfiEwasWater(log_lines, verbose = verbose, log_path = log_path)

  structure(
    list(
      data = analysis_data,
      phenotypes = phenotype_list,
      covariates = covariate_list,
      factorVars = factor_list,
      cpgColumns = cpg_columns,
      cpgPrefix = cpgPrefix,
      cpgLimit = cpg_limit,
      methylationScale = methylation_scale,
      responseLabel = methylation_label,
      methylationObjectPrefix = methylation_prefix,
      internalResponseColumn = response_column,
      prsMap = prs_map,
      interactionTerm = resolved_interaction,
      requestedInteractionTerm = interactionTerm,
      missingCounts = missing_counts,
      variableSummary = variable_summary,
      interactionTable = interaction_table
    ),
    class = "dnaEPICO_methylationGLM_data"
  )
}
#' Plot phenotype and covariate distributions for one-timepoint GLM analyses
#'
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
#' @param plotWidth Integer. TIFF width in pixels when plots are written to disk.
#' @param plotHeight Integer. TIFF height in pixels when plots are written to
#'   disk.
#' @param plotDPI Integer. TIFF resolution in DPI when plots are written to
#'   disk.
#' @param outputDir Character or `NULL`. Directory used for TIFF files. When
#'   `NULL`, plots are returned in memory only.
#' @param display Logical. If `TRUE`, draw plots on the active graphics device.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationGLM_distribution_plots"`
#'   containing the generated `ggplot2` objects and any saved TIFF file paths.
#'
#' @description
#' Create phenotype, factor-variable, and numeric-covariate distribution plots
#' from the object returned by `prepareMethylationGLMData()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' distribution_plots <- plotMethylationGLMDistributions(
#'   preparedData = ex$preparedData,
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(distribution_plots)
#'
#' @export
plotMethylationGLMDistributions <- function(
    preparedData,
    plotWidth = 2000L,
    plotHeight = 1000L,
    plotDPI = 150L,
    outputDir = NULL,
    display = FALSE,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationGLM.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  analysis_data <- preparedData$data
  phenotype_plots <- list()
  factor_plots <- list()
  covariate_plots <- list()
  saved_files <- list(phenotypes = list(), factors = list(), covariates = list())

  for (var in preparedData$phenotypes) {
    if (!var %in% colnames(analysis_data)) {
      next
    }

    plot_type <- if (is.numeric(analysis_data[[var]])) "hist" else "bar"
    plot_object <- createDistributionPlotMethylationGLM(
      values = analysis_data[[var]],
      variable = var,
      type = plot_type,
      fill = "steelblue"
    )
    phenotype_plots[[var]] <- plot_object

    file_path <- NULL
    if (!is.null(outputDir)) {
      file_path <- file.path(
        outputDir,
        paste0(if (identical(plot_type, "hist")) "hist_" else "bar_", var, ".tiff")
      )
    }

    runPlotMinfiEwasWater(
      draw_fun = function() drawPlotObjectMinfiEwasWater(plot_object),
      display = display,
      file = file_path,
      width = plotWidth,
      height = plotHeight,
      res = plotDPI
    )

    if (!is.null(file_path)) {
      saved_files$phenotypes[[var]] <- file_path
    }
  }

  factor_vars <- intersect(preparedData$factorVars, colnames(analysis_data))
  for (var in factor_vars) {
    plot_object <- createDistributionPlotMethylationGLM(
      values = analysis_data[[var]],
      variable = var,
      type = "bar",
      fill = "darkorange"
    )
    factor_plots[[var]] <- plot_object

    file_path <- NULL
    if (!is.null(outputDir)) {
      file_path <- file.path(outputDir, paste0("bar_", var, ".tiff"))
    }

    runPlotMinfiEwasWater(
      draw_fun = function() drawPlotObjectMinfiEwasWater(plot_object),
      display = display,
      file = file_path,
      width = plotWidth,
      height = plotHeight,
      res = plotDPI
    )

    if (!is.null(file_path)) {
      saved_files$factors[[var]] <- file_path
    }
  }

  numeric_covariates <- setdiff(preparedData$covariates, preparedData$factorVars)
  numeric_covariates <- intersect(numeric_covariates, colnames(analysis_data))
  for (var in numeric_covariates) {
    plot_object <- createDistributionPlotMethylationGLM(
      values = analysis_data[[var]],
      variable = var,
      type = "hist",
      fill = "darkgreen"
    )
    covariate_plots[[var]] <- plot_object

    file_path <- NULL
    if (!is.null(outputDir)) {
      file_path <- file.path(outputDir, paste0("hist_", var, ".tiff"))
    }

    runPlotMinfiEwasWater(
      draw_fun = function() drawPlotObjectMinfiEwasWater(plot_object),
      display = display,
      file = file_path,
      width = plotWidth,
      height = plotHeight,
      res = plotDPI
    )

    if (!is.null(file_path)) {
      saved_files$covariates[[var]] <- file_path
    }
  }

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      paste("Phenotype distribution plots: ", length(phenotype_plots)),
      paste("Factor distribution plots:    ", length(factor_plots)),
      paste("Numeric covariate plots:      ", length(covariate_plots)),
      if (is.null(outputDir)) {
        "Distribution plots were returned in memory only."
      } else {
        paste("Distribution plots saved to:  ", outputDir)
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      phenotypes = phenotype_plots,
      factors = factor_plots,
      covariates = covariate_plots,
      files = saved_files
    ),
    class = "dnaEPICO_methylationGLM_distribution_plots"
  )
}

#' Fit CpG-wise Gaussian GLMs for one-timepoint methylation analyses
#'
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
#' @param nCores Integer. Number of worker processes to use.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded
#'   to worker processes.
#' @param glmLibs Character vector or comma-separated string of package names to
#'   check on worker processes. The default is `"glm2"`.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationGLM_models"` containing
#'   fitted model lists, model formulas, and counts of failed CpG fits.
#'
#' @description
#' Fit one Gaussian GLM per CpG for each phenotype requested in the object
#' returned by `prepareMethylationGLMData()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' model_results <- fitMethylationGLMModels(
#'   preparedData = ex$preparedData,
#'   nCores = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(model_results$fits)
#'
#' @export
fitMethylationGLMModels <- function(
    preparedData,
    nCores = 1L,
    libPath = NULL,
    glmLibs = "glm2",
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationGLM.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)

  if (is.null(libPath)) {
    libPath <- .libPaths()
  }

  glm_lib_list <- splitOptionMinfiEwasWater(glmLibs, sep = ",")
  if (length(glm_lib_list) == 0L) {
    glm_lib_list <- "glm2"
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
    model_vars <- unique(c(phenotype, covariates, preparedData$interactionTerm))
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

    formula_text <- buildFormulaMethylationGLM(
      phenotype = phenotype,
      covariates = covariates,
      interactionTerm = preparedData$interactionTerm,
      responseVar = preparedData$internalResponseColumn
    )

    base_model_data <- analysis_data[, model_vars, drop = FALSE]
    factor_vars <- preparedData$factorVars
    for (var in intersect(factor_vars, colnames(base_model_data))) {
      base_model_data[[var]] <- as.factor(base_model_data[[var]])
    }
    cpg_batches <- chunkCpGColumnsMethylationModels(
      cpgColumns = cpg_columns,
      nCores = n_cores,
      batchesPerCore = 8L
    )
    batch_worker <- fitMethylationGLMBatch
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
              packages = glm_lib_list
            )
            batch_worker(
              cpgBatch = batch,
              data = analysis_data,
              modelData = base_model_data,
              formulaText = formula_text,
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
                "response_var",
                "batch_worker",
                "libPath",
                "glm_lib_list",
                "validateWorkerPackagesMethylationModels",
                "fitCpGModelMethylationGLMPrepared",
                "summarizeCpGFitMethylationGLM",
                "findCoefficientRowsMethylationGLM",
                "escapeRegexMethylationGLM"
              ),
              envir = environment()
            )

            parallel::clusterEvalQ(
              cl,
              validateWorkerPackagesMethylationModels(
                libPath = libPath,
                packages = glm_lib_list
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
        packages = glm_lib_list
      )
      batch_results <- lapply(
        cpg_batches,
        function(batch) {
          batch_worker(
            cpgBatch = batch,
            data = analysis_data,
            modelData = base_model_data,
            formulaText = formula_text,
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
    phenotype_summary_cache <- filterSummaryByPvalueMethylationGLM(
      summaryDf = combined_results$summaries,
      pValueFilter = NA_real_,
      includeResidualSD = TRUE
    )
    failures <- vapply(
      fit_list,
      function(x) inherits(x, "dnaEPICO_methylationGLM_fit_error"),
      logical(1)
    )

    fits[[phenotype]] <- fit_list
    summary_cache[[phenotype]] <- phenotype_summary_cache
    formulas[[phenotype]] <- formula_text
    failure_counts[[phenotype]] <- sum(failures)

    emitLogMinfiEwasWater(
      c(
        "=======================================================================",
        paste("Fitted phenotype:            ", phenotype),
        paste("Formula:                     ", formula_text),
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
        glmLibs = glm_lib_list,
        methylationScale = preparedData$methylationScale,
        methylationObjectPrefix = preparedData$methylationObjectPrefix,
        responseLabel = preparedData$responseLabel,
        internalResponseColumn = preparedData$internalResponseColumn,
        interactionTerm = preparedData$interactionTerm
      ),
      responseLabel = preparedData$responseLabel
    ),
    class = "dnaEPICO_methylationGLM_models"
  )
}

#' Summarize CpG-wise Gaussian GLM fits for one-timepoint analyses
#'
#' @param modelResults Object returned by `fitMethylationGLMModels()`.
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
#' @param summaryResidualSD Logical. If `TRUE`, add residual standard deviations
#'   to each CpG summary row.
#' @param summaryPval Numeric or `NA`. Optional p-value filter applied to the
#'   returned summary tables. `NA` keeps all rows.
#' @param nCores Integer. Number of worker processes to use while extracting
#'   summary rows.
#' @param libPath Character vector or `NULL`. Optional library paths forwarded
#'   to worker processes.
#' @param glmLibs Character vector or comma-separated string of package names to
#'   check on worker processes. The default is `"glm2"`.
#' @param chunkSize Integer or `NULL`. Number of CpGs to process per parallel
#'   chunk. `NULL` chooses a value automatically.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationGLM_summaries"`
#'   containing one CpG-level summary data frame per phenotype.
#'
#' @description
#' Extract phenotype-specific CpG coefficient tables from the fitted model
#' object returned by `fitMethylationGLMModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' summary_results <- summarizeMethylationGLMModels(
#'   modelResults = ex$modelResults,
#'   preparedData = ex$preparedData,
#'   summaryResidualSD = TRUE,
#'   summaryPval = NA,
#'   nCores = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(summary_results$summaries)
#'
#' @export
summarizeMethylationGLMModels <- function(
    modelResults,
    preparedData,
    summaryResidualSD = TRUE,
    summaryPval = NA,
    nCores = 1L,
    libPath = NULL,
    glmLibs = "glm2",
    chunkSize = NULL,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationGLM.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)

  if (is.null(libPath)) {
    libPath <- .libPaths()
  }

  glm_lib_list <- splitOptionMinfiEwasWater(glmLibs, sep = ",")
  if (length(glm_lib_list) == 0L) {
    glm_lib_list <- "glm2"
  }

  p_value_filter <- normalizeOptionalNumericMethylationGLM(summaryPval)
  chunk_size <- normalizeChunkSizeMethylationGLM(chunkSize)
  n_cores <- max(1L, as.integer(nCores))
  summaries <- list()

  for (phenotype in names(modelResults$fits)) {
    if (
      !is.null(modelResults$summaryCache) &&
        !is.null(modelResults$summaryCache[[phenotype]])
    ) {
      summary_df <- filterSummaryByPvalueMethylationGLM(
        summaryDf = modelResults$summaryCache[[phenotype]],
        pValueFilter = p_value_filter,
        includeResidualSD = isTRUE(summaryResidualSD)
      )
      summaries[[phenotype]] <- summary_df

      emitLogMinfiEwasWater(
        c(
          "=======================================================================",
          paste("Summarized phenotype:        ", phenotype),
          paste("CpG summary rows returned:   ", nrow(summary_df)),
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

    summary_worker <- summarizeCpGFitMethylationGLM
    resolved_interaction <- preparedData$interactionTerm
    include_residual_sd <- isTRUE(summaryResidualSD)

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
          "include_residual_sd",
          "summary_worker",
          "libPath",
          "glm_lib_list"
        ),
        envir = environment()
      )

      parallel::clusterEvalQ(
        cl,
        {
          if (!is.null(libPath)) {
            .libPaths(unique(c(libPath, .libPaths())))
          }

          for (pkg in glm_lib_list) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
              stop("Failed to load package: ", pkg, call. = FALSE)
            }
          }

          NULL
        }
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
                variable = phenotype,
                interactionTerm = resolved_interaction,
                includeResidualSD = include_residual_sd
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
                variable = phenotype,
                interactionTerm = resolved_interaction,
                includeResidualSD = include_residual_sd
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
      summary_df <- summary_df[, c(
        "CpG",
        "Coefficient",
        "Estimate",
        "Std. Error",
        "t value",
        "Pr(>|t|)",
        if (isTRUE(summaryResidualSD)) "ResidualSD"
      ), drop = FALSE]
      rownames(summary_df) <- NULL
    }

    if (nrow(summary_df) > 0L && !is.na(p_value_filter)) {
      summary_df <- summary_df[summary_df[["Pr(>|t|)"]] < p_value_filter, , drop = FALSE]
      rownames(summary_df) <- NULL
    }

    summaries[[phenotype]] <- summary_df

    emitLogMinfiEwasWater(
      c(
        "=======================================================================",
        paste("Summarized phenotype:        ", phenotype),
        paste("CpG summary rows returned:   ", nrow(summary_df)),
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
        summaryResidualSD = isTRUE(summaryResidualSD),
        summaryPval = p_value_filter,
        chunkSize = chunk_size
      )
    ),
    class = "dnaEPICO_methylationGLM_summaries"
  )
}

#' Collect significant CpG coefficient tables from fitted one-timepoint GLMs
#'
#' @param modelResults Object returned by `fitMethylationGLMModels()`.
#' @param pvalThreshold Numeric. Threshold applied to phenotype main-effect or
#'   interaction p-values.
#' @param interactionTerm Character or `NULL`. Optional interaction term.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationGLM_significant_cpgs"`.
#'
#' @description
#' Collect the raw coefficient tables for CpGs whose phenotype main effect or
#' interaction p-value passes the requested threshold.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' significant_cpgs <- collectSignificantCpGsMethylationGLM(
#'   modelResults = ex$modelResults,
#'   pvalThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(significant_cpgs)
#'
#' @export
collectSignificantCpGsMethylationGLM <- function(
    modelResults,
    pvalThreshold = 0.05,
    interactionTerm = NULL,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationGLM.txt"
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
        hit_cpgs <- unique(cached_summary$CpG[cached_summary[["Pr(>|t|)"]] < threshold])
        hit_cpgs <- hit_cpgs[!is.na(hit_cpgs) & hit_cpgs %in% names(fit_list)]
        for (cpg in hit_cpgs) {
          model_obj <- fit_list[[cpg]]
          if (is.null(model_obj) || inherits(model_obj, "dnaEPICO_methylationGLM_fit_error")) {
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
      if (is.null(model_obj) || inherits(model_obj, "dnaEPICO_methylationGLM_fit_error")) {
        next
      }

      coef_table <- model_obj$coef
      if (is.null(coef_table)) {
        next
      }

      matched_rows <- findCoefficientRowsMethylationGLM(
        coefNames = rownames(coef_table),
        variable = phenotype,
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
      paste("Significant CpGs retained at p <", threshold, ":"),
      paste(names(hit_counts), hit_counts, sep = ": ", collapse = "; "),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(retained, class = "dnaEPICO_methylationGLM_significant_cpgs")
}
#' Plot diagnostic summaries for one-timepoint methylation GLMs
#'
#' @param modelSummaries Object returned by `summarizeMethylationGLMModels()`.
#' @param preparedData Object returned by `prepareMethylationGLMData()`.
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
#' @return A list with class `"dnaEPICO_methylationGLM_diagnostic_plots"`
#'   containing the generated `ggplot2` objects, genomic inflation factors, and
#'   any saved TIFF file paths.
#'
#' @description
#' Create Q-Q and residual-diagnostic plots from the CpG summary tables returned
#' by `summarizeMethylationGLMModels()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' diagnostic_plots <- plotMethylationGLMDiagnostics(
#'   modelSummaries = ex$modelSummaries,
#'   preparedData = ex$preparedData,
#'   display = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(diagnostic_plots$plots)
#'
#' @export
plotMethylationGLMDiagnostics <- function(
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
    log_file = "log_methylationGLM.txt"
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

    summary_df$FDR <- stats::p.adjust(summary_df[["Pr(>|t|)"]], method = padjmethod)
    pvals <- summary_df[["Pr(>|t|)"]]
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
    residual_plot <- NULL
    residual_significance_plot <- NULL
    residual_file <- NULL
    significance_file <- NULL

    if ("ResidualSD" %in% colnames(summary_df)) {
      residual_plot <- ggplot2::ggplot(
        summary_df,
        ggplot2::aes(x = log2meanBeta, y = ResidualSD)
      ) +
        ggplot2::geom_point(alpha = 0.6, color = "black") +
        ggplot2::labs(
          title = "plotSA-style: Residual SD vs Average Beta",
          x = "log2(Average Beta)",
          y = "Residual Standard Deviation"
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
        ggplot2::aes(x = -log10(`Pr(>|t|)`), y = ResidualSD, color = FDR < fdrThreshold)
      ) +
        ggplot2::geom_point(alpha = 0.6) +
        ggrepel::geom_text_repel(
          data = subset(summary_df, FDR < fdrThreshold),
          ggplot2::aes(label = CpG)
        ) +
        ggplot2::scale_color_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick")) +
        ggplot2::labs(
          title = paste("Residual SD vs Significance for", phenotype),
          x = "-log10(p-value)",
          y = "Residual SD",
          color = paste("FDR <", fdrThreshold)
        ) +
        ggplot2::theme_minimal()
    }

    qq_file <- NULL
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

    if (!is.null(residual_plot)) {
      runPlotMinfiEwasWater(
        draw_fun = function() drawPlotObjectMinfiEwasWater(residual_plot),
        display = display,
        file = residual_file,
        width = plotWidth,
        height = plotHeight,
        res = plotDPI
      )
    }

    if (!is.null(residual_significance_plot)) {
      runPlotMinfiEwasWater(
        draw_fun = function() drawPlotObjectMinfiEwasWater(residual_significance_plot),
        display = display,
        file = significance_file,
        width = plotWidth,
        height = plotHeight,
        res = plotDPI
      )
    }

    plot_list[[phenotype]] <- list(
      qqplot = qq_plot,
      residualSD = residual_plot,
      residualSignificance = residual_significance_plot
    )
    saved_files[[phenotype]] <- list(
      qqplot = qq_file,
      residualSD = if (!is.null(residual_plot)) residual_file else NULL,
      residualSignificance = if (!is.null(residual_significance_plot)) significance_file else NULL
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
    class = "dnaEPICO_methylationGLM_diagnostic_plots"
  )
}

#' Annotate one-timepoint GLM summary tables with array annotation metadata
#'
#' @param modelSummaries Object returned by `summarizeMethylationGLMModels()`
#'   or a named list of CpG summary data frames.
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
#' @return A list with class `"dnaEPICO_methylationGLM_annotation"`
#'   containing the annotated summary table and any requested annotation columns
#'   that were unavailable in the chosen annotation object.
#'
#' @description
#' Merge phenotype-specific CpG summary tables with probe annotation metadata and
#' return a single annotated result table.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' annotation_data <- annotateMethylationGLMSummaries(
#'   modelSummaries = ex$modelSummaries,
#'   annotationObject = ex$annotationData,
#'   annotationCols = "Name,chr,pos",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(annotation_data)
#'
#' @export
annotateMethylationGLMSummaries <- function(
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
    log_file = "log_methylationGLM.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  summary_list <- modelSummaries
  if (!is.null(modelSummaries$summaries)) {
    summary_list <- modelSummaries$summaries
  }

  annotation_cols <- splitOptionMinfiEwasWater(annotationCols, sep = ",")
  annotation_df <- coerceAnnotationDataMethylationGLM(annotationObject)

  merged_summary_list <- lapply(
    names(summary_list),
    function(phenotype) {
      summary_df <- summary_list[[phenotype]]
      if (is.null(summary_df) || nrow(summary_df) == 0L) {
        return(NULL)
      }

      coefficient_names <- unique(summary_df$Coefficient)
      coefficient_tables <- lapply(
        coefficient_names,
        function(coefficient_name) {
          sub_df <- summary_df[
            summary_df$Coefficient == coefficient_name,
            c("CpG", "Pr(>|t|)"),
            drop = FALSE
          ]
          clean_name <- gsub("`", "", coefficient_name, fixed = TRUE)
          colnames(sub_df)[2] <- paste0(clean_name, "P.Value")
          sub_df
        }
      )

      if (length(coefficient_tables) == 0L) {
        return(NULL)
      }

      Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), coefficient_tables)
    }
  )
  merged_summary_list <- Filter(Negate(is.null), merged_summary_list)

  if (length(merged_summary_list) == 0L) {
    merged_summary <- data.frame(CpG = character(0))
  } else if (length(merged_summary_list) == 1L) {
    merged_summary <- merged_summary_list[[1L]]
  } else {
    merged_summary <- Reduce(function(x, y) merge(x, y, by = "CpG", all = TRUE), merged_summary_list)
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
    class = "dnaEPICO_methylationGLM_annotation"
  )
}

#' Write optional disk outputs for one-timepoint GLM analyses
#'
#' @param modelResults Object returned by `fitMethylationGLMModels()`.
#' @param modelSummaries Object returned by `summarizeMethylationGLMModels()`.
#' @param annotatedResults Object returned by
#'   `annotateMethylationGLMSummaries()` or a compatible data frame.
#' @param significantCpGs Object returned by
#'   `collectSignificantCpGsMethylationGLM()` or `NULL`.
#' @param outputRData Character. Directory used for serialized model and summary
#'   outputs.
#' @param summaryTxtDir Character. Directory used for tab-delimited summary
#'   tables.
#' @param significantCpGDir Character. Directory used for significant-CpG
#'   coefficient tables.
#' @param annotatedGLMOut Character. Directory used for the annotated summary
#'   XLSX workbook.
#' @param saveTxtSummaries Logical. If `TRUE`, write tab-delimited summary
#'   tables.
#' @param saveSignificantCpGs Logical. If `TRUE`, write significant-CpG
#'   coefficient tables.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same messages to a log file.
#' @param log_dir Character or `NULL`. Directory used for the log file when
#'   `logs = TRUE`.
#' @param log_file Character. File name used when `logs = TRUE`.
#'
#' @return A list with class `"dnaEPICO_methylationGLM_paths"` containing
#'   the paths of the files written to disk.
#'
#' @description
#' Write optional serialized outputs, summary tables, significant-CpG tables,
#' and annotated results from the one-timepoint GLM workflow.
#'
#' @examples
#' ex <- dnaEPICO:::exampleMethylationGLMStateDnaEpico()
#' annotation_data <- annotateMethylationGLMSummaries(
#'   modelSummaries = ex$modelSummaries,
#'   annotationObject = ex$annotationData,
#'   annotationCols = "Name,chr,pos",
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' significant_cpgs <- collectSignificantCpGsMethylationGLM(
#'   modelResults = ex$modelResults,
#'   pvalThreshold = 1,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' output_paths <- writeMethylationGLMOutputs(
#'   modelResults = ex$modelResults,
#'   modelSummaries = ex$modelSummaries,
#'   annotatedResults = annotation_data,
#'   significantCpGs = significant_cpgs,
#'   outputRData = file.path(ex$tempDir, "models"),
#'   summaryTxtDir = file.path(ex$tempDir, "summary"),
#'   significantCpGDir = file.path(ex$tempDir, "significant"),
#'   annotatedGLMOut = file.path(ex$tempDir, "annotated"),
#'   saveTxtSummaries = TRUE,
#'   saveSignificantCpGs = TRUE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' names(output_paths)
#'
#' @export
writeMethylationGLMOutputs <- function(
    modelResults,
    modelSummaries,
    annotatedResults,
    significantCpGs = NULL,
    outputRData,
    summaryTxtDir,
    significantCpGDir,
    annotatedGLMOut,
    saveTxtSummaries = TRUE,
    saveSignificantCpGs = FALSE,
    verbose = FALSE,
    logs = FALSE,
    log_dir = NULL,
    log_file = "log_methylationGLM.txt"
) {
  log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir, log_file = log_file)
  dir.create(outputRData, recursive = TRUE, showWarnings = FALSE)
  dir.create(annotatedGLMOut, recursive = TRUE, showWarnings = FALSE)

  model_files <- stats::setNames(character(length(modelResults$fits)), names(modelResults$fits))
  summary_files <- stats::setNames(character(length(modelSummaries$summaries)), names(modelSummaries$summaries))
  summary_txt_files <- list()
  significant_files <- list()

  for (phenotype in names(modelResults$fits)) {
    model_file <- file.path(outputRData, paste0(phenotype, "GLM.rds"))
    summary_file <- file.path(outputRData, paste0(phenotype, "SummaryGLM.rds"))
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

      if ("Pr(>|t|)" %in% colnames(summary_df)) {
        summary_df <- summary_df[order(summary_df[["Pr(>|t|)"]]), , drop = FALSE]
      }

      output_file <- file.path(summaryTxtDir, paste0(phenotype, "SummaryGLM.txt"))
      utils::write.table(summary_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
      summary_txt_files[[phenotype]] <- output_file
    }
  }

  if (isTRUE(saveSignificantCpGs) && !is.null(significantCpGs)) {
    dir.create(significantCpGDir, recursive = TRUE, showWarnings = FALSE)
    for (phenotype in names(significantCpGs)) {
      phenotype_hits <- significantCpGs[[phenotype]]
      if (length(phenotype_hits) == 0L) {
        next
      }

      significant_files[[phenotype]] <- character(0)
      phenotype_dir <- file.path(significantCpGDir, phenotype)
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
  annotated_file <- file.path(annotatedGLMOut, "annotatedGLM.xlsx")
  dictionary <- buildAnnotatedWorkbookDictionaryMethylationGLM(
    columns = colnames(annotated_df),
    modelDescription = "Pvalue from GLM model",
    formulaText = modelResults$formulas,
    modelLabel = "GLM",
    responseLabel = inferMethylationValueLabelMethylationGLM(modelResults)
  )
  writeAnnotatedWorkbookMethylationGLM(
    annotated_df = annotated_df,
    file = annotated_file,
    resultSheet = "annotatedGLM",
    dictionary = dictionary
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
      if (isTRUE(saveSignificantCpGs)) {
        paste("Significant CpG text files:  ", sum(vapply(significant_files, length, integer(1))))
      } else {
        "Significant CpG text files:  0"
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
      significantCpGFiles = significant_files,
      annotatedGLM = annotated_file
    ),
    class = "dnaEPICO_methylationGLM_paths"
  )
}
