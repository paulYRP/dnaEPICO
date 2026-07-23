emptyModelMessagesDnaEpico <- function() {
    data.frame(
        Phenotype = character(0), CpG = character(0),
        Model.Message = character(0), P.Value.Available = logical(0),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

emptyFitFailuresMethylationModels <- function() {
    data.frame(
        Phenotype = character(0), CpG = character(0),
        Status = character(0), Error = character(0),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

parseFiniteNumericMethylationModels <- function(value) {
    if (length(value) != 1L || is.na(value)) {
        return(NA_real_)
    }
    value <- as.character(value)
    numeric_pattern <- paste0(
        "^[[:space:]]*[+-]?(",
        "[0-9]+([.][0-9]*)?|[.][0-9]+",
        ")([eE][+-]?[0-9]+)?[[:space:]]*$"
    )
    if (!grepl(numeric_pattern, value)) {
        return(NA_real_)
    }
    as.numeric(value)
}

inputIdentityMethylationModels <- function(path) {
    normalized <- normalizePath(path, winslash = "/", mustWork = TRUE)
    info <- file.info(normalized)
    list(
        path = normalized,
        size = unname(as.numeric(info$size[[1L]])),
        modified = unname(as.numeric(info$mtime[[1L]])),
        md5 = unname(as.character(tools::md5sum(normalized)[[1L]]))
    )
}

cpgResponseMethylationModels <- function(data, cpg) {
    if (is.matrix(data)) {
        return(data[, cpg, drop = TRUE])
    }
    data[[cpg]]
}

compactCoefficientResultsMethylationModels <- function(
    fits, cpgOrder,
    includeResidualSD = FALSE
) {
    successful <- fits[!vapply(fits, function(x) {
        is.null(x) || inherits(x, c(
            "dnaEPICO_methylationGLM_fit_error",
            "dnaEPICO_methylationLME_fit_error"
        ))
    }, logical(1))]
    coefficient_names <- unique(unlist(lapply(
        successful,
        function(x) rownames(x$coef)
    ), use.names = FALSE))
    coefficient_names <- coefficient_names[
        !is.na(coefficient_names) & nzchar(coefficient_names)
    ]

    make_matrix <- function() {
        matrix(
            NA_real_, nrow = length(cpgOrder),
            ncol = length(coefficient_names),
            dimnames = list(cpgOrder, coefficient_names)
        )
    }
    estimate <- make_matrix()
    std_error <- make_matrix()
    degrees_freedom <- make_matrix()
    statistic <- make_matrix()
    p_value <- make_matrix()
    residual_sd <- if (isTRUE(includeResidualSD)) {
        stats::setNames(rep(NA_real_, length(cpgOrder)), cpgOrder)
    } else {
        NULL
    }
    coefficient_terms <- stats::setNames(
        rep(NA_character_, length(coefficient_names)),
        coefficient_names
    )

    for (cpg in intersect(cpgOrder, names(successful))) {
        fit <- successful[[cpg]]
        coef_table <- fit$coef
        if (is.null(coef_table) || nrow(coef_table) == 0L) {
            next
        }
        coefficient_rows <- intersect(rownames(coef_table), coefficient_names)
        copy_column <- function(target, source) {
            if (!(source %in% colnames(coef_table))) {
                return(target)
            }
            target[cpg, coefficient_rows] <-
                as.numeric(coef_table[coefficient_rows, source, drop = TRUE])
            target
        }
        estimate <- copy_column(estimate, "Estimate")
        std_error <- copy_column(std_error, "Std. Error")
        degrees_freedom <- copy_column(degrees_freedom, "df")
        statistic <- copy_column(statistic, "t value")
        p_value <- copy_column(p_value, "Pr(>|t|)")
        if (!is.null(residual_sd) && !is.null(fit$residualSD)) {
            residual_sd[[cpg]] <- as.numeric(fit$residualSD[[1L]])
        }
        if (!is.null(fit$coefficientTerms)) {
            mapped <- intersect(
                names(fit$coefficientTerms),
                names(coefficient_terms)
            )
            coefficient_terms[mapped] <-
                as.character(fit$coefficientTerms[mapped])
        }
    }

    if (length(degrees_freedom) == 0L ||
        !any(is.finite(degrees_freedom))) {
        degrees_freedom <- NULL
    }
    list(
        cpgOrder = cpgOrder,
        coefficientNames = coefficient_names,
        coefficientTerms = coefficient_terms,
        estimate = estimate, stdError = std_error,
        df = degrees_freedom, statistic = statistic,
        pValue = p_value, residualSD = residual_sd
    )
}

combineCompactCoefficientResultsMethylationModels <- function(
    batchResults,
    cpgOrder
) {
    compact_batches <- lapply(batchResults, function(x) x$coefficientResults)
    compact_batches <- Filter(Negate(is.null), compact_batches)
    coefficient_names <- unique(unlist(lapply(
        compact_batches,
        function(x) x$coefficientNames
    ), use.names = FALSE))
    coefficient_names <- coefficient_names[
        !is.na(coefficient_names) & nzchar(coefficient_names)
    ]

    make_matrix <- function() {
        matrix(
            NA_real_, nrow = length(cpgOrder),
            ncol = length(coefficient_names),
            dimnames = list(cpgOrder, coefficient_names)
        )
    }
    combined <- list(
        cpgOrder = cpgOrder,
        coefficientNames = coefficient_names,
        coefficientTerms = stats::setNames(
            rep(NA_character_, length(coefficient_names)),
            coefficient_names
        ),
        estimate = make_matrix(), stdError = make_matrix(),
        df = make_matrix(), statistic = make_matrix(),
        pValue = make_matrix(),
        residualSD = stats::setNames(
            rep(NA_real_, length(cpgOrder)),
            cpgOrder
        )
    )

    for (batch in compact_batches) {
        rows <- intersect(batch$cpgOrder, cpgOrder)
        columns <- intersect(batch$coefficientNames, coefficient_names)
        if (length(rows) > 0L && length(columns) > 0L) {
            for (field in c(
                "estimate", "stdError", "statistic",
                "pValue"
            )) {
                combined[[field]][rows, columns] <-
                    batch[[field]][rows, columns, drop = FALSE]
            }
            if (!is.null(batch$df)) {
                combined$df[rows, columns] <-
                    batch$df[rows, columns, drop = FALSE]
            }
        }
        if (!is.null(batch$residualSD)) {
            combined$residualSD[rows] <- batch$residualSD[rows]
        }
        if (!is.null(batch$coefficientTerms)) {
            mapped <- intersect(
                names(batch$coefficientTerms),
                names(combined$coefficientTerms)
            )
            missing_terms <- is.na(combined$coefficientTerms[mapped])
            mapped <- mapped[missing_terms]
            combined$coefficientTerms[mapped] <-
                as.character(batch$coefficientTerms[mapped])
        }
    }

    if (!any(is.finite(combined$df))) {
        combined$df <- NULL
    }
    if (!any(is.finite(combined$residualSD))) {
        combined$residualSD <- NULL
    }
    combined
}

coefficientTableFromCompactMethylationModels <- function(
    coefficientResults,
    cpg
) {
    if (is.null(coefficientResults) ||
        !(cpg %in% coefficientResults$cpgOrder)) {
        return(NULL)
    }
    coefficient_names <- coefficientResults$coefficientNames
    if (length(coefficient_names) == 0L) {
        return(NULL)
    }
    table <- data.frame(
        Estimate = as.numeric(
            coefficientResults$estimate[cpg, coefficient_names]
        ),
        `Std. Error` = as.numeric(
            coefficientResults$stdError[cpg, coefficient_names]
        ),
        check.names = FALSE
    )
    if (!is.null(coefficientResults$df)) {
        table$df <- as.numeric(
            coefficientResults$df[cpg, coefficient_names]
        )
    }
    table[["t value"]] <- as.numeric(
        coefficientResults$statistic[cpg, coefficient_names]
    )
    table[["Pr(>|t|)"]] <- as.numeric(
        coefficientResults$pValue[cpg, coefficient_names]
    )
    rownames(table) <- coefficient_names
    table
}

collectBatchModelMessagesMethylationModels <- function(
    fits,
    phenotype
) {
    if (length(fits) == 0L) {
        return(emptyModelMessagesDnaEpico())
    }
    out <- data.frame(
        Phenotype = rep(phenotype, length(fits)),
        CpG = names(fits),
        Model.Message = vapply(
            fits, modelMessageDnaEpico, character(1)
        ),
        P.Value.Available = vapply(
            fits, function(x) isTRUE(x$pValueAvailable), logical(1)
        ),
        stringsAsFactors = FALSE, check.names = FALSE
    )
    rownames(out) <- NULL
    out
}

collectBatchFitFailuresMethylationModels <- function(
    fits,
    phenotype,
    errorClass
) {
    failed <- vapply(fits, inherits, logical(1), what = errorClass)
    if (!any(failed)) {
        return(emptyFitFailuresMethylationModels())
    }
    failed_fits <- fits[failed]
    data.frame(
        Phenotype = rep(phenotype, length(failed_fits)),
        CpG = names(failed_fits), Status = "failed",
        Error = vapply(
            failed_fits, function(x) as.character(x$error), character(1)
        ),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

combineBatchTablesMethylationModels <- function(
    batchResults,
    field,
    empty
) {
    tables <- lapply(batchResults, function(x) x[[field]])
    tables <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, tables)
    if (length(tables) == 0L) {
        return(empty)
    }
    out <- do.call(rbind, tables)
    rownames(out) <- NULL
    out
}

packageVersionsMethylationModels <- function(packages) {
    packages <- unique(c("dnaEPICO", as.character(packages)))
    stats::setNames(vapply(packages, function(package) {
        if (!requireNamespace(package, quietly = TRUE)) {
            return(NA_character_)
        }
        as.character(utils::packageVersion(package))
    }, character(1)), packages)
}

factorSpecificationMethylationModels <- function(data, factorVars) {
    variables <- intersect(factorVars, colnames(data))
    lapply(stats::setNames(variables, variables), function(variable) {
        value <- data[[variable]]
        list(
            levels = levels(value),
            ordered = is.ordered(value),
            contrasts = if (is.factor(value)) {
                unclass(stats::contrasts(value))
            } else {
                NULL
            }
        )
    })
}

buildPhenotypeSignatureMethylationModels <- function(
    analysis, engine, phenotype, formulaText, preparedData,
    modelSettings = list(), packages = character(0)
) {
    list(
        schemaVersion = 1L,
        analysis = analysis, engine = engine,
        phenotype = phenotype, formula = formulaText,
        inputIdentity = preparedData$inputIdentity,
        dimensions = dim(preparedData$data),
        cpgOrder = preparedData$cpgColumns,
        cpgPrefix = preparedData$cpgPrefix,
        cpgLimit = preparedData$cpgLimit,
        methylationScale = preparedData$methylationScale,
        internalResponseColumn = preparedData$internalResponseColumn,
        covariates = preparedData$covariates,
        phenotypePrs = if (phenotype %in% names(preparedData$prsMap)) {
            unname(preparedData$prsMap[[phenotype]])
        } else {
            character(0)
        },
        interactionTerm = preparedData$interactionTerm,
        factors = factorSpecificationMethylationModels(
            preparedData$data,
            preparedData$factorVars
        ),
        scaleVars = preparedData$scaleVars,
        scalingMetadata = preparedData$scalingMetadata,
        modelSettings = modelSettings,
        packageVersions = packageVersionsMethylationModels(packages),
        rVersion = paste(
            R.version$major,
            R.version$minor,
            sep = "."
        )
    )
}

phenotypeSummaryPathMethylationModels <- function(
    outputDir,
    phenotype,
    analysis
) {
    suffix <- if (identical(analysis, "glm")) "SummaryGLM.rds" else
        "SummaryLME.rds"
    file.path(outputDir, paste0(phenotype, suffix))
}

assemblePhenotypeSummaryMethylationModels <- function(
    analysis, engine, phenotype, signature, cpgOrder,
    coefficientResults, targetSummary, omnibusTests,
    modelMessages, fitFailures, failureCount, failureReasons,
    formulaText, settings
) {
    structure(list(
        formatVersion = 1L, complete = TRUE,
        analysis = analysis, engine = engine, phenotype = phenotype,
        signature = signature, cpgOrder = cpgOrder,
        completion = stats::setNames(
            rep(TRUE, length(cpgOrder)),
            cpgOrder
        ),
        coefficientResults = coefficientResults,
        targetSummary = targetSummary,
        omnibusTests = omnibusTests,
        modelMessages = modelMessages,
        fitFailures = fitFailures,
        failureCount = as.integer(failureCount),
        failureReasons = failureReasons,
        formula = formulaText, settings = settings,
        versions = signature$packageVersions
    ), class = "dnaEPICO_methylation_phenotype_summary")
}

validatePhenotypeSummaryMethylationModels <- function(
    object,
    expectedSignature
) {
    required <- c(
        "formatVersion", "complete", "analysis", "engine", "phenotype",
        "signature", "cpgOrder", "completion", "coefficientResults",
        "targetSummary", "modelMessages", "fitFailures", "failureCount",
        "failureReasons", "formula", "settings", "versions"
    )
    if (!is.list(object) || !all(required %in% names(object))) {
        return(list(valid = FALSE, reason = "required fields are missing"))
    }
    if (!identical(object$formatVersion, 1L) ||
        !isTRUE(object$complete)) {
        return(list(valid = FALSE, reason = "summary is incomplete"))
    }
    if (!identical(object$signature, expectedSignature)) {
        return(list(
            valid = FALSE,
            reason = "the input file or model configuration changed"
        ))
    }
    expected_cpgs <- expectedSignature$cpgOrder
    if (!identical(object$cpgOrder, expected_cpgs) ||
        !identical(names(object$completion), expected_cpgs) ||
        length(object$completion) != length(expected_cpgs) ||
        !all(object$completion)) {
        return(list(valid = FALSE, reason = "CpG completion is incomplete"))
    }
    compact <- object$coefficientResults
    if (!is.list(compact) ||
        !identical(compact$cpgOrder, expected_cpgs)) {
        return(list(valid = FALSE, reason = "coefficient results are invalid"))
    }
    if (!is.data.frame(object$targetSummary) ||
        !is.data.frame(object$modelMessages) ||
        nrow(object$modelMessages) != length(expected_cpgs) ||
        !identical(as.character(object$modelMessages$CpG), expected_cpgs) ||
        !is.data.frame(object$fitFailures)) {
        return(list(valid = FALSE, reason = "result tables are invalid"))
    }
    list(valid = TRUE, reason = "valid")
}

loadPhenotypeSummaryMethylationModels <- function(
    path,
    expectedSignature
) {
    if (is.null(path) || !file.exists(path)) {
        return(list(object = NULL, reason = "summary file is absent"))
    }
    object <- tryCatch(readRDS(path), error = function(error) error)
    if (inherits(object, "error")) {
        return(list(
            object = NULL,
            reason = paste("summary cannot be read:", conditionMessage(object))
        ))
    }
    validation <- validatePhenotypeSummaryMethylationModels(
        object,
        expectedSignature
    )
    if (!isTRUE(validation$valid)) {
        return(list(object = NULL, reason = validation$reason))
    }
    list(object = object, reason = "valid")
}

savePhenotypeSummaryMethylationModels <- function(
    object,
    path
) {
    validation <- validatePhenotypeSummaryMethylationModels(
        object,
        object$signature
    )
    if (!isTRUE(validation$valid)) {
        stop(
            "The phenotype summary is invalid: ",
            validation$reason, call. = FALSE
        )
    }
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(
        pattern = paste0(".", basename(path), "-"),
        tmpdir = dirname(path), fileext = ".tmp"
    )
    backup <- tempfile(
        pattern = paste0(".", basename(path), "-"),
        tmpdir = dirname(path), fileext = ".backup"
    )
    on.exit(unlink(temporary), add = TRUE)
    on.exit(unlink(backup), add = TRUE)
    saveRDS(object, file = temporary, compress = TRUE)
    temporary_info <- file.info(temporary)
    if (!file.exists(temporary) || is.na(temporary_info$size[[1L]]) ||
        temporary_info$size[[1L]] <= 0) {
        stop("The phenotype summary could not be written: ", path,
            call. = FALSE
        )
    }
    had_existing <- file.exists(path)
    if (had_existing && !file.rename(path, backup)) {
        stop("Could not preserve the existing phenotype summary: ",
            path, call. = FALSE
        )
    }
    if (!file.rename(temporary, path)) {
        if (had_existing && file.exists(backup)) {
            file.rename(backup, path)
        }
        stop("Could not finalize phenotype summary: ", path, call. = FALSE)
    }
    if (had_existing && file.exists(backup)) {
        unlink(backup)
    }
    invisible(path)
}

availableMemoryMbMethylationModels <- function() {
    configured <- Sys.getenv("DNAEPICO_AVAILABLE_MEMORY_MB", unset = "")
    if (nzchar(configured)) {
        value <- parseFiniteNumericMethylationModels(configured)
        if (length(value) == 1L && is.finite(value) && value > 0) {
            return(value)
        }
    }
    candidates <- numeric(0)
    if (file.exists("/proc/meminfo")) {
        info <- tryCatch(
            readLines("/proc/meminfo", warn = FALSE),
            error = function(error) character(0)
        )
        line <- grep("^MemAvailable:", info, value = TRUE)
        if (length(line) == 1L) {
            kb <- parseFiniteNumericMethylationModels(gsub(
                "[^0-9]", "",
                line
            ))
            if (is.finite(kb)) {
                candidates <- c(candidates, kb / 1024)
            }
        }
    }
    cgroup_limit <- "/sys/fs/cgroup/memory.max"
    cgroup_used <- "/sys/fs/cgroup/memory.current"
    if (file.exists(cgroup_limit) && file.exists(cgroup_used)) {
        limit <- parseFiniteNumericMethylationModels(readLines(
            cgroup_limit,
            n = 1L, warn = FALSE
        ))
        used <- parseFiniteNumericMethylationModels(readLines(
            cgroup_used,
            n = 1L, warn = FALSE
        ))
        if (is.finite(limit) && is.finite(used) && limit > used) {
            candidates <- c(candidates, (limit - used) / 1024^2)
        }
    }
    slurm_node <- parseFiniteNumericMethylationModels(Sys.getenv(
        "SLURM_MEM_PER_NODE",
        unset = ""
    ))
    if (is.finite(slurm_node) && slurm_node > 0) {
        candidates <- c(candidates, slurm_node)
    }
    slurm_per_cpu <- parseFiniteNumericMethylationModels(Sys.getenv(
        "SLURM_MEM_PER_CPU",
        unset = ""
    ))
    slurm_cpus <- parseFiniteNumericMethylationModels(Sys.getenv(
        "SLURM_CPUS_PER_TASK",
        unset = ""
    ))
    if (is.finite(slurm_per_cpu) && slurm_per_cpu > 0 &&
        is.finite(slurm_cpus) && slurm_cpus > 0) {
        candidates <- c(candidates, slurm_per_cpu * slurm_cpus)
    }
    candidates <- candidates[is.finite(candidates) & candidates > 0]
    if (length(candidates) == 0L) {
        return(NA_real_)
    }
    min(candidates)
}

estimateWorkerMemoryMethylationModels <- function(
    engine,
    modelData,
    nSamples
) {
    model_mb <- as.numeric(utils::object.size(modelData)) / 1024^2
    engine_mb <- switch(
        tolower(engine),
        lme4 = 384,
        nlme = 256,
        glm2 = 128,
        384
    )
    max(engine_mb, model_mb * 4 + as.numeric(nSamples) * 0.02)
}

memoryWorkerCapMethylationModels <- function(
    engine,
    analysisData,
    modelData,
    requestedWorkers
) {
    available_mb <- availableMemoryMbMethylationModels()
    worker_mb <- estimateWorkerMemoryMethylationModels(
        engine,
        modelData,
        nrow(analysisData)
    )
    input_mb <- as.numeric(utils::object.size(analysisData)) / 1024^2
    reserve_mb <- max(2048, input_mb * 1.5)
    if (!is.finite(available_mb)) {
        return(list(
            cap = requestedWorkers, availableMB = NA_real_,
            reserveMB = reserve_mb, estimatedWorkerMB = worker_mb,
            reason = "available memory could not be detected"
        ))
    }
    usable_mb <- max(0, available_mb - reserve_mb)
    cap <- max(1L, floor(usable_mb / worker_mb))
    list(
        cap = min(as.integer(requestedWorkers), as.integer(cap)),
        availableMB = available_mb, reserveMB = reserve_mb,
        estimatedWorkerMB = worker_mb,
        reason = "worker count capped by detected available memory"
    )
}

measurePilotMemoryMethylationModels <- function(fitFunction) {
    before <- gc(reset = TRUE)
    value <- fitFunction()
    after <- gc()
    incremental_mb <- max(
        1,
        sum(after[, 6L], na.rm = TRUE) -
            sum(before[, 2L], na.rm = TRUE)
    )
    list(value = value, incrementalMB = incremental_mb)
}

refineParallelPlanWithPilotMethylationModels <- function(
    plan,
    pilotMemoryMB
) {
    empirical_worker_mb <- max(
        plan$estimatedWorkerMemoryMB,
        as.numeric(pilotMemoryMB) * 1.5
    )
    plan$pilotMemoryMB <- as.numeric(pilotMemoryMB)
    plan$estimatedWorkerMemoryMB <- empirical_worker_mb
    if (is.finite(plan$availableMemoryMB)) {
        usable_mb <- max(
            0,
            plan$availableMemoryMB - plan$reservedMemoryMB
        )
        plan$memoryWorkerCap <- max(
            1L,
            as.integer(floor(usable_mb / empirical_worker_mb))
        )
        plan$workerCount <- min(
            plan$workerCount,
            plan$memoryWorkerCap
        )
        if (plan$workerCount <= 1L) {
            plan$backend <- "serial"
            plan$reason <- "one worker after empirical memory cap"
        }
    }
    plan
}
