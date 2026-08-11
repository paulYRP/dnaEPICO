# Shared model-visualisation helpers ----------------------------------------

genomeWideThresholdDnaEpico <- function() {
    9e-8
}

dnaEpicoModelPlotTheme <- function() {
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major = ggplot2::element_line(
        colour = "#E7EDF2", linewidth = 0.35
        ),
        axis.text = ggplot2::element_text(colour = "#243746"),
        axis.title = ggplot2::element_text(colour = "#17324D"),
        legend.position = "bottom",
        plot.margin = ggplot2::margin(12, 18, 12, 12)
    )
}

adaptivePointStyleDnaEpico <- function(n) {
    n <- max(0L, as.integer(n))
    if (n <= 30L) {
    return(list(size = 2.4, alpha = 0.82))
    }
    if (n <= 500L) {
    return(list(size = 1.7, alpha = 0.62))
    }
    if (n <= 5000L) {
    return(list(size = 1.05, alpha = 0.42))
    }
    if (n <= 100000L) {
    return(list(size = 0.55, alpha = 0.28))
    }
    list(size = 0.28, alpha = 0.2)
}

deterministicPlotRowsDnaEpico <- function(
    n, maximum = 50000L, priority = integer()
) {
    n <- max(0L, as.integer(n))
    maximum <- max(1L, as.integer(maximum))
    priority <- unique(as.integer(priority))
    priority <- priority[is.finite(priority) & priority >= 1L & priority <= n]
    if (n <= maximum) {
    return(seq_len(n))
    }
    background_maximum <- max(0L, maximum - length(priority))
    background <- if (background_maximum > 0L) {
    unique(as.integer(round(seq(1L, n, length.out = background_maximum))))
    } else {
    integer()
    }
    sort(unique(c(priority, background)))
}

adaptiveFigureDimensionDnaEpico <- function(
    base, items, pixelsPerItem = 90L, maximum = 4800L
) {
    max(as.integer(base), min(
    as.integer(maximum),
    as.integer(items) * as.integer(pixelsPerItem)
    ))
}

isContinuousModelVariableDnaEpico <- function(values) {
    is.numeric(values) && length(unique(values[is.finite(values)])) > 6L
}

continuousModelAssociationPlotDnaEpico <- function(
    phenotypeValues, variableValues, phenotype, variable, maximumPoints
) {
    data <- data.frame(
    x = coerceNumericDnaEpico(variableValues),
    y = coerceNumericDnaEpico(phenotypeValues)
    )
    data <- data[is.finite(data$x) & is.finite(data$y), , drop = FALSE]
    if (!nrow(data)) {
    return(NULL)
    }
    plot <- ggplot2::ggplot(data, ggplot2::aes(x = x, y = y))
    if (nrow(data) > maximumPoints) {
    plot <- plot + ggplot2::geom_bin_2d(bins = 45) +
        ggplot2::scale_fill_gradient(
        low = "#D8EFF3", high = "#176B87", name = "Observations"
        )
    } else {
    style <- adaptivePointStyleDnaEpico(nrow(data))
    plot <- plot + ggplot2::geom_point(
        size = style$size, alpha = style$alpha, colour = "#176B87"
    )
    }
    plot + ggplot2::geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    colour = "#B42318", fill = "#F4B8B3", linewidth = 0.8
    ) + ggplot2::labs(
    title = NULL, x = variable, y = phenotype
    ) + dnaEpicoModelPlotTheme()
}

prepareMixedAssociationDataDnaEpico <- function(
    phenotypeValues, variableValues, phenotypeNumeric, maximumPoints
) {
    numeric_values <- if (phenotypeNumeric) {
        coerceNumericDnaEpico(phenotypeValues)
    } else {
        coerceNumericDnaEpico(variableValues)
    }
    group_values <- if (phenotypeNumeric) variableValues else phenotypeValues
    data <- data.frame(
        group = as.character(group_values), value = numeric_values,
        stringsAsFactors = FALSE
    )
    missing <- is.na(data$group) | !nzchar(trimws(data$group))
    data$group[missing] <- "Missing"
    data <- data[is.finite(data$value), , drop = FALSE]
    if (!nrow(data)) return(NULL)
    group_order <- names(sort(table(data$group), decreasing = TRUE))
    data$group <- factor(data$group, levels = group_order)
    rows <- deterministicPlotRowsDnaEpico(
        nrow(data), maximum = min(2500L, maximumPoints)
    )
    point_data <- data[rows, , drop = FALSE]
    sizes <- table(data$group)
    violin_data <- data[
        data$group %in% names(sizes[sizes >= 3L]), , drop = FALSE
    ]
    list(
        data = data, groupOrder = group_order, pointData = point_data,
        violinData = violin_data,
        style = adaptivePointStyleDnaEpico(nrow(point_data))
    )
}

mixedAssociationAxisThemeDnaEpico <- function(groupOrder) {
    angle <- if (length(groupOrder) > 5L) 40 else 0
    ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = angle, hjust = 1)
    )
}

mixedModelAssociationPlotDnaEpico <- function(
    phenotypeValues, variableValues, phenotype, variable,
    phenotypeNumeric, maximumPoints
) {
    prepared <- prepareMixedAssociationDataDnaEpico(
        phenotypeValues, variableValues, phenotypeNumeric, maximumPoints
    )
    if (is.null(prepared)) return(NULL)
    x_label <- if (phenotypeNumeric) variable else phenotype
    y_label <- if (phenotypeNumeric) phenotype else variable
    ggplot2::ggplot(
        prepared$data, ggplot2::aes(x = group, y = value)
    ) +
        ggplot2::geom_violin(
            data = prepared$violinData,
            fill = "#D8EFF3", colour = "#176B87",
            alpha = 0.65, scale = "width", trim = TRUE
        ) +
        ggplot2::geom_boxplot(
            width = 0.18, outlier.shape = NA,
            fill = "white", colour = "#17324D", alpha = 0.9
        ) +
        ggplot2::geom_jitter(
            data = prepared$pointData, width = 0.08, height = 0,
            size = min(prepared$style$size, 1.2),
            alpha = min(prepared$style$alpha, 0.38), colour = "#243746"
        ) +
        ggplot2::labs(title = NULL, x = x_label, y = y_label) +
        dnaEpicoModelPlotTheme() +
        mixedAssociationAxisThemeDnaEpico(prepared$groupOrder)
}

categoricalModelAssociationPlotDnaEpico <- function(
    phenotypeValues, variableValues, phenotype, variable
) {
    data <- data.frame(
    phenotype = as.character(phenotypeValues),
    variable = as.character(variableValues), stringsAsFactors = FALSE
    )
    for (column in names(data)) {
    missing <- is.na(data[[column]]) | !nzchar(trimws(data[[column]]))
    data[[column]][missing] <- "Missing"
    }
    counts <- as.data.frame(table(data), stringsAsFactors = FALSE)
    counts <- counts[counts$Freq > 0L, , drop = FALSE]
    if (!nrow(counts)) {
    return(NULL)
    }
    totals <- stats::ave(counts$Freq, counts$phenotype, FUN = sum)
    counts$percentage <- 100 * counts$Freq / totals
    counts$label <- sprintf("%s\n%.1f%%", counts$Freq, counts$percentage)
    ggplot2::ggplot(counts, ggplot2::aes(
    x = variable, y = phenotype, fill = percentage
    )) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.6) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5) +
    ggplot2::scale_fill_gradient(
        low = "#EAF4F6", high = "#176B87", name = "Row %"
    ) +
    ggplot2::labs(title = NULL, x = variable, y = phenotype) +
    dnaEpicoModelPlotTheme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 40, hjust = 1
    ))
}

createModelAssociationPlotDnaEpico <- function(
    data, phenotype, variable, maximumPoints = 5000L
) {
    if (identical(phenotype, variable) ||
    !all(c(phenotype, variable) %in% names(data))) {
    return(NULL)
    }
    phenotype_values <- data[[phenotype]]
    variable_values <- data[[variable]]
    phenotype_numeric <- isContinuousModelVariableDnaEpico(phenotype_values)
    variable_numeric <- isContinuousModelVariableDnaEpico(variable_values)
    if (phenotype_numeric && variable_numeric) {
    return(continuousModelAssociationPlotDnaEpico(
        phenotype_values, variable_values, phenotype, variable,
        maximumPoints
    ))
    }
    if (xor(phenotype_numeric, variable_numeric)) {
    return(mixedModelAssociationPlotDnaEpico(
        phenotype_values, variable_values, phenotype, variable,
        phenotype_numeric, maximumPoints
    ))
    }
    categoricalModelAssociationPlotDnaEpico(
    phenotype_values, variable_values, phenotype, variable
    )
}

modelVariablesDnaEpico <- function(preparedData, includeLongitudinal = FALSE) {
    variables <- unique(c(
    preparedData$phenotypes, preparedData$covariates,
    preparedData$factorVars, preparedData$interactionTerm,
    if (isTRUE(includeLongitudinal)) preparedData$timeVar else NULL
    ))
    variables[!is.na(variables) & nzchar(variables) &
    variables %in% names(preparedData$data)]
}

modelMissingnessPlotDnaEpico <- function(data, variables) {
    missing <- data.frame(
    variable = variables,
    missing = vapply(data, function(value) {
        sum(is.na(value) | (is.character(value) & !nzchar(trimws(value))))
    }, numeric(1)), stringsAsFactors = FALSE
    )
    missing$percentage <- if (nrow(data)) {
    100 * missing$missing / nrow(data)
    } else {
    0
    }
    missing$variable <- stats::reorder(missing$variable, missing$percentage)
    ggplot2::ggplot(missing, ggplot2::aes(x = variable, y = percentage)) +
    ggplot2::geom_col(fill = "#176B87", alpha = 0.85) +
    ggplot2::geom_text(
        ggplot2::aes(label = sprintf(
        "%s (%.1f%%)", format(missing, big.mark = ","), percentage
        )),
        hjust = -0.08, size = 3.5
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_y_continuous(
        limits = c(0, max(5, missing$percentage * 1.18)),
        expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
        title = NULL, x = "Model variable",
        y = "Missing observations (%)"
    ) +
    dnaEpicoModelPlotTheme()
}

modelCorrelationDataDnaEpico <- function(data, variables) {
    numeric_data <- data[, variables, drop = FALSE]
    correlation <- stats::cor(
    numeric_data,
    use = "pairwise.complete.obs", method = "spearman"
    )
    correlation[is.finite(correlation) & abs(correlation) < 5e-4] <- 0
    complete_n <- outer(variables, variables, Vectorize(function(x, y) {
    sum(stats::complete.cases(numeric_data[, c(x, y), drop = FALSE]))
    }))
    dimnames(complete_n) <- list(variables, variables)
    result <- as.data.frame(as.table(correlation), stringsAsFactors = FALSE)
    names(result) <- c("variable1", "variable2", "correlation")
    result$completeN <- as.vector(complete_n)
    result$index1 <- match(result$variable1, variables)
    result$index2 <- match(result$variable2, variables)
    result <- result[result$index1 >= result$index2, , drop = FALSE]
    result$label <- ifelse(
    result$index1 == result$index2,
    paste0("n=", format(result$completeN, big.mark = ",")),
    sprintf(
        "rho=%.2f\nn=%s", result$correlation,
        format(result$completeN, big.mark = ",")
    )
    )
    result
}

modelCorrelationPlotDnaEpico <- function(data) {
    ggplot2::ggplot(data, ggplot2::aes(
    x = variable1, y = variable2, fill = correlation
    )) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.4) +
    ggplot2::geom_text(
        ggplot2::aes(label = label),
        size = 3.1, lineheight = 0.92
    ) +
    ggplot2::scale_fill_gradient2(
        low = "#B42318", mid = "white", high = "#176B87",
        midpoint = 0, limits = c(-1, 1)
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
        title = NULL, x = NULL, y = NULL, fill = "Spearman rho"
    ) +
    dnaEpicoModelPlotTheme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
        angle = 45, hjust = 1
    ))
}

saveModelDesignPlotDnaEpico <- function(
    plot, file, display, width, height, resolution
) {
    runPlotMinfiEwasWater(
    draw_fun = function() drawPlotObjectMinfiEwasWater(plot),
    display = display, file = file, width = width,
    height = height, res = resolution
    )
    invisible(NULL)
}

plotModelAssociationsDnaEpico <- function(
    data, preparedData, variables, analysis, outputDir,
    plotWidth, plotHeight, plotDPI, display
) {
    association_vars <- unique(c(
    preparedData$covariates, preparedData$factorVars,
    preparedData$interactionTerm,
    if (identical(analysis, "LME")) preparedData$timeVar else NULL
    ))
    association_vars <- association_vars[
    !is.na(association_vars) & nzchar(association_vars) &
        association_vars %in% variables
    ]
    plots <- files <- list()
    for (phenotype in intersect(preparedData$phenotypes, variables)) {
    for (variable in setdiff(association_vars, phenotype)) {
        plot <- createModelAssociationPlotDnaEpico(
        data, phenotype, variable
        )
        if (is.null(plot)) {
        next
        }
        key <- paste0("association_", phenotype, "_by_", variable)
        filename <- paste0(
        "association_", safeFigureComponentDnaEpico(phenotype),
        "_by_", safeFigureComponentDnaEpico(variable), ".tiff"
        )
        file <- if (is.null(outputDir)) NULL else file.path(outputDir, filename)
        levels <- if (is.numeric(data[[variable]])) {
        1L
        } else {
        length(unique(data[[variable]][!is.na(data[[variable]])]))
        }
        saveModelDesignPlotDnaEpico(
        plot, file, display,
        adaptiveFigureDimensionDnaEpico(
            plotWidth, levels,
            pixelsPerItem = 170L
        ), plotHeight, plotDPI
        )
        plots[[key]] <- plot
        files[[key]] <- file
    }
    }
    list(plots = plots, files = files)
}

plotModelDesignDnaEpico <- function(preparedData, analysis = c("GLM",
    "LME"), outputDir = NULL, plotWidth = 2000L, plotHeight = 1000L,
    plotDPI = 150L, display = FALSE) {
    analysis <- match.arg(analysis)
    variables <- modelVariablesDnaEpico(preparedData, includeLongitudinal =
        identical(analysis,
        "LME"))
    data <- preparedData$data[, variables, drop = FALSE]
    plots <- files <- list()
    if (!length(variables)) {
        return(structure(list(plots = plots, files = files, variables =
            variables),
            class = "dnaEPICO_model_design_plots"))
    }
    missing_plot <- modelMissingnessPlotDnaEpico(data, variables)
    missing_file <- if (is.null(outputDir)) {
        NULL } else {
        file.path(outputDir, "modelVariables_missingnessPercentage.tiff")
    }
    saveModelDesignPlotDnaEpico(missing_plot, missing_file, display,
        plotWidth, adaptiveFigureDimensionDnaEpico(plotHeight,
            length(variables), pixelsPerItem = 105L), plotDPI)
    plots$missingness <- missing_plot
    files$missingness <- missing_file
    numeric_vars <- setdiff(variables[vapply(data, is.numeric,
        logical(1))], preparedData$factorVars)
    if (length(numeric_vars) >= 2L) {
        correlation_plot <- modelCorrelationPlotDnaEpico(
            modelCorrelationDataDnaEpico(data,
            numeric_vars))
        correlation_file <- if (is.null(outputDir)) {
            NULL } else {
            file.path(outputDir, "modelVariables_spearmanCorrelation.tiff")
        }
        saveModelDesignPlotDnaEpico(correlation_plot, correlation_file,
            display, adaptiveFigureDimensionDnaEpico(plotWidth,
                length(numeric_vars), pixelsPerItem = 145L),
            adaptiveFigureDimensionDnaEpico(plotHeight,
                length(numeric_vars), pixelsPerItem = 125L), plotDPI)
        plots$correlation <- correlation_plot
        files$correlation <- correlation_file
    }
    associations <- plotModelAssociationsDnaEpico(data, preparedData,
        variables, analysis, outputDir, plotWidth, plotHeight,
        plotDPI, display)
    plots <- c(plots, associations$plots)
    files <- c(files, associations$files)
    structure(list(plots = plots, files = files, variables = variables),
        class = "dnaEPICO_model_design_plots")
}

plotLmeVariableDistributionsDnaEpico <- function(
    data, variables, phenotypes, outputDir, display,
    width, height, resolution
) {
    plots <- files <- list()
    for (variable in variables) {
    type <- if (is.numeric(data[[variable]])) "hist" else "bar"
    plot <- createDistributionPlotMethylationGLM(
        values = data[[variable]], variable = variable, type = type,
        fill = if (variable %in% phenotypes) "#2C7FB8" else "#008A83"
    )
    suffix <- if (identical(type, "hist")) "continuous" else "categorical"
    filename <- paste0(
        "distribution_", safeFigureComponentDnaEpico(variable),
        "_", suffix, ".tiff"
    )
    file <- if (is.null(outputDir)) NULL else file.path(outputDir, filename)
    saveModelDesignPlotDnaEpico(
        plot, file, display, width, height, resolution
    )
    plots[[variable]] <- plot
    files[[variable]] <- file
    }
    list(plots = plots, files = files)
}

lmeObservationCountPlotDnaEpico <- function(person) {
    counts <- as.data.frame(table(person), stringsAsFactors = FALSE)
    names(counts) <- c("person", "observations")
    ggplot2::ggplot(counts, ggplot2::aes(x = observations)) +
    ggplot2::geom_histogram(
        binwidth = 1, boundary = 0.5, fill = "#176B87",
        colour = "white"
    ) +
    ggplot2::labs(
        title = NULL, x = "Observations per participant",
        y = "Participants"
    ) +
    dnaEpicoModelPlotTheme()
}

lmeTrajectorySummaryDnaEpico <- function(values) {
    observed <- values[is.finite(values)]
    centre <- if (length(observed)) mean(observed) else NA_real_
    error <- if (length(observed) > 1L) {
    1.96 * stats::sd(observed) / sqrt(length(observed))
    } else {
    0
    }
    data.frame(y = centre, ymin = centre - error, ymax = centre + error)
}

lmeTrajectoryPlotDnaEpico <- function(data, displayed, timeVar, phenotype) {
    ggplot2::ggplot(displayed, ggplot2::aes(
    x = time, y = value, group = person
    )) +
    ggplot2::geom_line(alpha = 0.18, colour = "#64748B") +
    ggplot2::geom_point(alpha = 0.3, colour = "#176B87", size = 1.1) +
    ggplot2::stat_summary(
        data = data, ggplot2::aes(group = 1), fun = mean,
        geom = "line", colour = "#B42318", linewidth = 1.15
    ) +
    ggplot2::stat_summary(
        data = data, ggplot2::aes(group = 1),
        fun.data = lmeTrajectorySummaryDnaEpico,
        geom = "errorbar", colour = "#B42318", width = 0.08
    ) +
    ggplot2::labs(title = NULL, x = timeVar, y = phenotype) +
    dnaEpicoModelPlotTheme()
}

preparePairedLmeDataDnaEpico <- function(data, timeValues) {
    observed_times <- unique(timeValues[!is.na(timeValues)])
    if (length(observed_times) != 2L) {
    return(NULL)
    }
    time_order <- if (is.factor(timeValues)) {
    intersect(levels(timeValues), as.character(observed_times))
    } else {
    sort(observed_times)
    }
    data$time <- factor(data$time, levels = time_order)
    data <- data[order(data$person, data$time), , drop = FALSE]
    counts <- table(data$person, data$time)
    people <- rownames(counts)[
    rowSums(counts > 0L) == 2L & rowSums(counts) == 2L
    ]
    data[data$person %in% people, , drop = FALSE]
}

lmePairedChangePlotDnaEpico <- function(data, displayed, timeVar, phenotype) {
    ggplot2::ggplot(displayed, ggplot2::aes(
    x = time, y = value, group = person
    )) +
    ggplot2::geom_line(
        colour = "#64748B", alpha = 0.28, linewidth = 0.45
    ) +
    ggplot2::geom_point(
        colour = "#176B87", alpha = 0.55, size = 1.2
    ) +
    ggplot2::stat_summary(
        data = data, ggplot2::aes(group = 1), fun = mean,
        geom = "line", colour = "#B42318", linewidth = 1.25
    ) +
    ggplot2::stat_summary(
        data = data, ggplot2::aes(group = 1), fun = mean,
        geom = "point", colour = "#B42318", size = 2.8
    ) +
    ggplot2::labs(title = NULL, x = timeVar, y = phenotype) +
    dnaEpicoModelPlotTheme()
}

selectDisplayedPeopleDnaEpico <- function(data, maximum) {
    people <- unique(data$person)
    if (length(people) > maximum) {
    people <- people[deterministicPlotRowsDnaEpico(
        length(people),
        maximum = maximum
    )]
    }
    data[data$person %in% people, , drop = FALSE]
}

plotLmePhenotypeTrajectoryDnaEpico <- function(data,
    person, timeValues, phenotype, timeVar, outputDir,
    display, width, height, resolution) {
    complete <- data.frame(person = person, time = timeValues,
        value = data[[phenotype]], stringsAsFactors = FALSE)
    complete <- complete[stats::complete.cases(complete) &
        is.finite(complete$value), , drop = FALSE]
    trajectory <- lmeTrajectoryPlotDnaEpico(complete,
        selectDisplayedPeopleDnaEpico(complete, 150L),
        timeVar, phenotype)
    trajectory_file <- if (is.null(outputDir)) {
        NULL
    }
    else {
        file.path(outputDir, paste0("longitudinalTrajectory_",
            safeFigureComponentDnaEpico(phenotype), "_by_",
            safeFigureComponentDnaEpico(timeVar), ".tiff"))
    }
    saveModelDesignPlotDnaEpico(trajectory, trajectory_file,
        display, width, height, resolution)
    trajectory_key <- paste0("trajectory_", phenotype)
    plots <- stats::setNames(list(trajectory), trajectory_key)
    files <- list()
    if (!is.null(trajectory_file)) {
        files[[trajectory_key]] <- trajectory_file
    }
    paired <- preparePairedLmeDataDnaEpico(complete,
        timeValues)
    if (!is.null(paired) && nrow(paired)) {
        change <- lmePairedChangePlotDnaEpico(paired,
            selectDisplayedPeopleDnaEpico(paired, 200L),
            timeVar, phenotype)
        change_file <- if (is.null(outputDir)) {
            NULL
        }
        else {
            file.path(outputDir, paste0("pairedChange_",
                safeFigureComponentDnaEpico(phenotype),
                "_by_", safeFigureComponentDnaEpico(timeVar),
                ".tiff"))
        }
        saveModelDesignPlotDnaEpico(change, change_file,
            display, width, height, resolution)
        plots[[paste0("change_", phenotype)]] <- change
        if (!is.null(change_file)) {
            files[[paste0("change_", phenotype)]] <- change_file
        }
    }
    list(plots = plots, files = files)
}

addLmeObservationDistributionDnaEpico <- function(
    result, person, personVar, outputDir, display,
    plotWidth, plotHeight, plotDPI
) {
    plot <- lmeObservationCountPlotDnaEpico(person)
    file <- if (is.null(outputDir)) {
        NULL
    } else {
        file.path(
            outputDir,
            paste0(
                "participantObservationCount_",
                safeFigureComponentDnaEpico(personVar), ".tiff"
            )
        )
    }
    saveModelDesignPlotDnaEpico(
        plot, file, display, plotWidth, plotHeight, plotDPI
    )
    result$plots$observations <- plot
    if (!is.null(file)) result$files$observations <- file
    result
}

addLmeTimeDistributionDnaEpico <- function(
    result, timeValues, timeVar, outputDir, display,
    plotWidth, plotHeight, plotDPI
) {
    time_type <- if (is.numeric(timeValues)) "hist" else "bar"
    plot <- createDistributionPlotMethylationGLM(
        timeValues, timeVar, type = time_type, fill = "#7A5195"
    )
    file <- if (is.null(outputDir)) {
        NULL
    } else {
        file.path(
            outputDir,
            paste0(
                "timepointDistribution_",
                safeFigureComponentDnaEpico(timeVar), ".tiff"
            )
        )
    }
    saveModelDesignPlotDnaEpico(
        plot, file, display, plotWidth, plotHeight, plotDPI
    )
    result$plots$timeDistribution <- plot
    if (!is.null(file)) result$files$timeDistribution <- file
    result
}

addLmeTrajectoryDistributionsDnaEpico <- function(
    result, data, person, timeValues, preparedData,
    outputDir, display, plotWidth, plotHeight, plotDPI
) {
    numeric_phenotypes <- preparedData$phenotypes[vapply(
        data[preparedData$phenotypes], is.numeric, logical(1)
    )]
    for (phenotype in numeric_phenotypes) {
        item <- plotLmePhenotypeTrajectoryDnaEpico(
            data, person, timeValues, phenotype, preparedData$timeVar,
            outputDir, display, plotWidth, plotHeight, plotDPI
        )
        result$plots <- c(result$plots, item$plots)
        result$files <- c(result$files, item$files)
    }
    result
}

plotMethylationLMEDistributions <- function(
    preparedData, plotWidth = 2000L, plotHeight = 1000L,
    plotDPI = 150L, outputDir = NULL, display = FALSE
) {
    data <- preparedData$data
    variables <- modelVariablesDnaEpico(
        preparedData, includeLongitudinal = TRUE
    )
    result <- plotLmeVariableDistributionsDnaEpico(
        data, setdiff(variables, preparedData$timeVar),
        preparedData$phenotypes, outputDir, display,
        plotWidth, plotHeight, plotDPI
    )
    person <- as.character(data[[preparedData$personVar]])
    result <- addLmeObservationDistributionDnaEpico(
        result, person, preparedData$personVar, outputDir, display,
        plotWidth, plotHeight, plotDPI
    )
    time_values <- data[[preparedData$timeVar]]
    result <- addLmeTimeDistributionDnaEpico(
        result, time_values, preparedData$timeVar, outputDir, display,
        plotWidth, plotHeight, plotDPI
    )
    result <- addLmeTrajectoryDistributionsDnaEpico(
        result, data, person, time_values, preparedData,
        outputDir, display, plotWidth, plotHeight, plotDPI
    )
    structure(result, class = "dnaEPICO_methylationLME_distribution_plots")
}

safeFigureComponentDnaEpico <- function(value, fallback = "value") {
    value <- trimws(as.character(value[[1L]]))
    value <- gsub("[<>:\"/\\\\|?*]+", "_", value)
    value <- gsub("[[:space:]]+", "_", value)
    value <- gsub("_+", "_", value)
    value <- gsub("^[._ ]+|[. _]+$", "", value)
    if (nzchar(value)) value else fallback
}

identifyAnnotatedPvalueColumnsDnaEpico <- function(data, dictionary = NULL) {
    columns <- names(data)
    pvalue_columns <- columns[
    grepl("P\\.Value$|P\\.value$", columns) &
        !grepl("Adjusted\\.P\\.Value$|Adjusted\\.P\\.value$", columns)
    ]

    if (is.data.frame(dictionary) &&
    all(c("Column", "Description") %in% names(dictionary))) {
    dictionary_columns <- as.character(dictionary$Column[
        grepl("p-value", dictionary$Description, ignore.case = TRUE) &
        !grepl("adjust", dictionary$Description, ignore.case = TRUE)
    ])
    documented <- intersect(pvalue_columns, dictionary_columns)
    if (length(documented)) {
        pvalue_columns <- documented
    }
    }

    pvalue_columns
}

normalizeAutosomeDnaEpico <- function(chromosome) {
    chromosome <- toupper(trimws(as.character(chromosome)))
    chromosome <- sub("^CHR", "", chromosome)
    coerceIntegerDnaEpico(chromosome)
}

prepareManhattanDataDnaEpico <- function(annotatedData, pValueColumn,
    cpgColumn = "IlmnID", chromosomeColumn = "chr", positionColumn = "pos",
    chromosomeGap = 5e+06) {
    required <- c(cpgColumn, chromosomeColumn, positionColumn, pValueColumn)
    missing <- setdiff(required, names(annotatedData))
    if (length(missing)) {
        missing_text <- paste(missing, collapse = ", ")
        stop(sprintf("Manhattan input is missing column(s): %s",
            missing_text), call. = FALSE) }
    chromosome <- normalizeAutosomeDnaEpico(annotatedData[[chromosomeColumn]])
    position <- coerceNumericDnaEpico(annotatedData[[positionColumn]])
    pvalue <- coerceNumericDnaEpico(annotatedData[[pValueColumn]])
    valid <- chromosome %in% seq_len(22L) & is.finite(position) &
        position >= 0 & is.finite(pvalue) & pvalue >= 0 & pvalue <=
        1
    data <- data.frame(CpG = as.character(annotatedData[[cpgColumn]][valid]),
        chromosome = chromosome[valid], position = position[valid],
        pvalue = pvalue[valid], stringsAsFactors = FALSE)
    if (!nrow(data)) {
        return(list(data = data, chromosomeSummary = data.frame(),
            tested = 0L, excluded = length(valid), zeroPvalues = 0L))
    }; data$pvalue[data$pvalue == 0] <- .Machine$double.xmin
    data <- data[order(data$chromosome, data$position, data$CpG),
        , drop = FALSE]
    chromosome_max <- vapply(split(data$position, data$chromosome),
        max, numeric(1), na.rm = TRUE)
    present <- as.integer(names(chromosome_max))
    chromosome_max <- chromosome_max[order(present)]
    present <- as.integer(names(chromosome_max))
    offsets <- numeric(length(chromosome_max))
    if (length(offsets) > 1L) {
        offsets[-1L] <- cumsum(chromosome_max[-length(chromosome_max)] +
            chromosomeGap) }; names(offsets) <- names(chromosome_max)
    data$cumulativePosition <- data$position + offsets[as.character(
        data$chromosome)]; data$minusLog10P <- -log10(data$pvalue)
    data$displayCpG <- sub("_(TC|BC)[0-9]+$", "", data$CpG, ignore.case = TRUE)
    chromosome_summary <- data.frame(chromosome = present, start = unname(
        offsets[as.character(present)]),
        end = unname(offsets[as.character(present)] + chromosome_max),
        minimum = vapply(split(data$cumulativePosition, data$chromosome),
            min, numeric(1), na.rm = TRUE), maximum = vapply(split(
            data$cumulativePosition,
            data$chromosome), max, numeric(1), na.rm = TRUE),
            stringsAsFactors = FALSE)
    chromosome_summary$midpoint <- (chromosome_summary$start +
        chromosome_summary$end)/2
    list(data = data, chromosomeSummary = chromosome_summary, tested = nrow(
        data),
        excluded = sum(!valid), zeroPvalues = sum(pvalue[valid] ==
            0, na.rm = TRUE)) }

selectManhattanLabelsDnaEpico <- function(
    data, suggestiveThreshold = 1e-5,
    genomeWideThreshold = genomeWideThresholdDnaEpico(), maximum = 20L
) {
    hits <- data[data$pvalue < genomeWideThreshold, , drop = FALSE]
    if (!nrow(hits)) {
    hits <- data[data$pvalue < suggestiveThreshold, , drop = FALSE]
    }
    hits <- hits[order(hits$pvalue, hits$CpG), , drop = FALSE]
    utils::head(hits, max(0L, as.integer(maximum)))
}

prepareManhattanDisplayDnaEpico <- function(
    prepared, suggestiveThreshold, genomeWideThreshold,
    maximumLabels, maximumBackgroundPoints
) {
    data <- prepared$data
    labels <- selectManhattanLabelsDnaEpico(
    data, suggestiveThreshold, genomeWideThreshold, maximumLabels
    )
    priority <- union(
    which(data$pvalue < suggestiveThreshold),
    which(data$pvalue < genomeWideThreshold)
    )
    rows <- deterministicPlotRowsDnaEpico(
    nrow(data), maximumBackgroundPoints, priority
    )
    display <- data[rows, , drop = FALSE]
    levels <- c("Genome-wide", "Suggestive", "Background")
    display$significance <- "Background"
    display$significance[display$pvalue < suggestiveThreshold] <- "Suggestive"
    display$significance[display$pvalue < genomeWideThreshold] <- "Genome-wide"
    display$significance <- factor(display$significance, levels = levels)
    labels$significance <- factor(ifelse(
    labels$pvalue < genomeWideThreshold, "Genome-wide", "Suggestive"
    ), levels = levels)
    list(data = data, display = display, labels = labels, levels = levels)
}

manhattanChromosomeColoursDnaEpico <- function() {
    stats::setNames(
    rep(RColorBrewer::brewer.pal(10L, "Paired"), length.out = 22L),
    as.character(seq_len(22L))
    )
}

manhattanSignificanceColoursDnaEpico <- function() {
    c(
    `Genome-wide` = "#C83E4D", Suggestive = "#E59B24",
    Background = "#56616C"
    )
}

radialManhattanConfigDnaEpico <- function(
    data, labels, chromosomes, suggestiveThreshold, genomeWideThreshold,
    levels
) {
    radial_limit <- ceiling(max(data$minusLog10P))
    data$radius <- radial_limit - data$minusLog10P
    labels$radius <- radial_limit - labels$minusLog10P
    breaks <- pretty(c(0, radial_limit), n = 5L)
    breaks <- breaks[breaks >= 0 & breaks <= radial_limit]
    grid <- data.frame(
    minusLog10P = breaks, radius = radial_limit - breaks,
    cumulativePosition = 0
    )
    thresholds <- data.frame(
    threshold = factor(c("Genome-wide", "Suggestive"), levels = levels),
    radius = radial_limit - -log10(c(
        genomeWideThreshold, suggestiveThreshold
    ))
    )
    count <- nrow(labels)
    inner <- if (count > 50L) -7 else if (count > 15L) -3.5 else -1.8
    label_size <- if (count > 50L) 2.1 else if (count > 15L) 3 else 3.5
    padding <- if (count > 50L) 0.08 else 0.35
    ring <- c(
    inner = radial_limit + 0.35, outer = radial_limit + 1.15,
    label = radial_limit + 1.85
    )
    thresholds <- thresholds[
    thresholds$radius >= inner & thresholds$radius <= ring[["label"]] + 0.35, ,
    drop = FALSE
    ]
    maximum_position <- max(chromosomes$end)
    chromosomes$angle <- 90 - 360 * chromosomes$midpoint / maximum_position
    adjust <- chromosomes$angle < -90
    chromosomes$angle[adjust] <- chromosomes$angle[adjust] + 180
    chromosomes$labelOffset <- ifelse(chromosomes$chromosome >= 19L, 0.25, 0)
    list(
    data = data, labels = labels, grid = grid, thresholds = thresholds,
    chromosomes = chromosomes, radialLimit = radial_limit,
    innerLimit = inner, labelSize = label_size,
    padding = padding, ring = ring, maximumPosition = maximum_position
    )
}

radialManhattanBaseDnaEpico <- function(config) {
    ggplot2::ggplot(config$data, ggplot2::aes(
    x = cumulativePosition, y = radius
    )) +
    ggplot2::geom_hline(
        yintercept = config$grid$radius,
        colour = "#D7D7D7", linewidth = 0.35
    ) +
    ggplot2::geom_hline(
        data = config$thresholds,
        ggplot2::aes(yintercept = radius, colour = threshold),
        inherit.aes = FALSE, linewidth = 0.55, show.legend = FALSE
    ) +
    ggplot2::annotate(
        "segment",
        x = 0, xend = 0, y = 0,
        yend = config$radialLimit, colour = "#4D545B", linewidth = 0.45
    ) +
    ggplot2::geom_text(
        data = config$grid,
        ggplot2::aes(
        x = cumulativePosition, y = radius, label = minusLog10P
        ), inherit.aes = FALSE, hjust = -0.25, vjust = -0.25,
        size = 3.5, colour = "#303030"
    )
}

addRadialManhattanPointsDnaEpico <- function(
    plot, config, chromosomeColours, significanceColours
) {
    plot + ggplot2::geom_rect(
    data = config$chromosomes,
    ggplot2::aes(
        xmin = start, xmax = end, ymin = config$ring[["inner"]],
        ymax = config$ring[["outer"]], fill = factor(chromosome)
    ), inherit.aes = FALSE, colour = "white", linewidth = 0.8
    ) +
    ggplot2::geom_point(
        ggplot2::aes(
        fill = factor(chromosome), colour = significance,
        size = significance
        ),
        shape = 21, alpha = 0.78, stroke = 0.55,
        show.legend = FALSE
    ) +
    ggplot2::geom_text(
        data = config$chromosomes,
        ggplot2::aes(
        x = midpoint, y = config$ring[["label"]] + labelOffset,
        label = paste0("chr", chromosome), angle = angle
        ), inherit.aes = FALSE, fontface = "bold", size = 3.5,
        colour = "#111111"
    ) +
    ggplot2::scale_fill_manual(values = chromosomeColours, guide = "none") +
    ggplot2::scale_colour_manual(
        values = significanceColours, guide = "none"
    ) +
    ggplot2::scale_size_manual(values = c(
        `Genome-wide` = 2.8, Suggestive = 1.15, Background = 0.28
    ), guide = "none")
}

finishRadialManhattanPlotDnaEpico <- function(plot, config) {
    plot <- plot + ggplot2::scale_x_continuous(
    limits = c(0, config$maximumPosition), expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
        limits = c(
        config$innerLimit, config$ring[["label"]] + 0.35
        ), expand = c(0, 0)
    ) + ggplot2::coord_polar(start = 0, clip = "off") +
    ggplot2::theme_void(base_size = 13) +
    ggplot2::theme(
        legend.position = "none",
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(5, 25, 5, 25)
    )
    if (nrow(config$labels)) {
    plot <- plot + ggrepel::geom_text_repel(
        data = config$labels,
        ggplot2::aes(
        x = cumulativePosition, y = radius, label = displayCpG
        ), inherit.aes = FALSE, size = config$labelSize,
        colour = "#303030", segment.colour = "#7A838B",
        segment.linetype = "dashed", segment.size = 0.35,
        box.padding = config$padding, point.padding = config$padding,
        min.segment.length = 0, max.overlaps = Inf,
        force = 1.8, seed = 2026, show.legend = FALSE
    )
    }
    plot
}

linearManhattanBaseDnaEpico <- function(
    data, chromosomes, suggestiveThreshold, genomeWideThreshold,
    chromosomeColours, significanceColours
) {
    boundaries <- if (nrow(chromosomes) > 1L) {
    (utils::head(chromosomes$end, -1L) +
        utils::tail(chromosomes$start, -1L)) / 2
    } else {
    numeric(0)
    }
    ggplot2::ggplot(data, ggplot2::aes(
    x = cumulativePosition, y = minusLog10P
    )) +
    ggplot2::geom_vline(
        xintercept = boundaries, colour = "#E2E5E8", linewidth = 0.25
    ) +
    ggplot2::geom_hline(
        yintercept = -log10(suggestiveThreshold),
        colour = "#E59B24", linewidth = 0.65
    ) +
    ggplot2::geom_hline(
        yintercept = -log10(genomeWideThreshold),
        colour = "#C83E4D", linewidth = 0.65
    ) +
    ggplot2::geom_point(
        ggplot2::aes(
        fill = factor(chromosome), colour = significance,
        size = significance
        ),
        shape = 21, alpha = 0.78, stroke = 0.45,
        show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = chromosomeColours, guide = "none") +
    ggplot2::scale_colour_manual(
        values = significanceColours, guide = "none"
    ) +
    ggplot2::scale_size_manual(values = c(
        `Genome-wide` = 3, Suggestive = 1.2, Background = 0.3
    ), guide = "none")
}

finishLinearManhattanPlotDnaEpico <- function(
    plot, labels, data, chromosomes, suggestiveThreshold, genomeWideThreshold
) {
    plot <- plot + ggplot2::scale_x_continuous(
    breaks = chromosomes$midpoint,
    labels = paste0("chr", chromosomes$chromosome),
    expand = ggplot2::expansion(mult = c(0.005, 0.005))
    ) + ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0.01, 0.16))
    ) + ggplot2::labs(
    title = NULL, x = NULL, y = expression(-log[10](italic(p)))
    ) + ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
        legend.position = "none",
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(
        colour = "#E5E7E9", linewidth = 0.3
        ),
        axis.text.x = ggplot2::element_text(
        colour = "#252525", angle = 45, hjust = 1, size = 11
        ),
        axis.text.y = ggplot2::element_text(colour = "#252525", size = 11),
        axis.title = ggplot2::element_text(colour = "#252525", size = 13),
        plot.background = ggplot2::element_rect(fill = "white", colour = NA),
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.margin = ggplot2::margin(12, 20, 12, 12)
    )
    if (nrow(labels)) {
    limits <- if (any(data$pvalue < genomeWideThreshold)) {
        c(-log10(genomeWideThreshold) + 0.15, Inf)
    } else {
        c(
        -log10(suggestiveThreshold) + 0.15,
        -log10(genomeWideThreshold) - 0.15
        )
    }
    plot <- plot + ggrepel::geom_text_repel(
        data = labels,
        ggplot2::aes(
        x = cumulativePosition, y = minusLog10P, label = displayCpG
        ), inherit.aes = FALSE, size = 4.2, colour = "#252525",
        segment.colour = "#C83E4D", segment.size = 0.4,
        box.padding = 0.45, point.padding = 0.3,
        min.segment.length = 0, max.overlaps = Inf,
        force = 1.8, ylim = limits, seed = 2026, show.legend = FALSE
    )
    }
    plot
}

createManhattanPlotDnaEpico <- function(
    prepared, version = c("v1", "v2"),
    suggestiveThreshold = 1e-5,
    genomeWideThreshold = genomeWideThresholdDnaEpico(),
    maximumLabels = 20L, maximumBackgroundPoints = 200000L
) {
    version <- match.arg(version)
    if (!nrow(prepared$data)) {
    return(NULL)
    }
    display <- prepareManhattanDisplayDnaEpico(
    prepared, suggestiveThreshold, genomeWideThreshold,
    maximumLabels, maximumBackgroundPoints
    )
    chromosome_colours <- manhattanChromosomeColoursDnaEpico()
    significance_colours <- manhattanSignificanceColoursDnaEpico()
    chromosomes <- prepared$chromosomeSummary
    if (identical(version, "v1")) {
    config <- radialManhattanConfigDnaEpico(
        display$display, display$labels, chromosomes,
        suggestiveThreshold, genomeWideThreshold, display$levels
    )
    plot <- radialManhattanBaseDnaEpico(config)
    plot <- addRadialManhattanPointsDnaEpico(
        plot, config, chromosome_colours, significance_colours
    )
    return(finishRadialManhattanPlotDnaEpico(plot, config))
    }
    plot <- linearManhattanBaseDnaEpico(
    display$display, chromosomes,
    suggestiveThreshold, genomeWideThreshold,
    chromosome_colours, significance_colours
    )
    finishLinearManhattanPlotDnaEpico(
    plot, display$labels, display$data, chromosomes,
    suggestiveThreshold, genomeWideThreshold
    )
}

plotAnnotatedManhattanColumnDnaEpico <- function(
    data, pValueColumn, analysis, outputDir, display,
    plotWidth, plotHeight, plotDPI,
    suggestiveThreshold, genomeWideThreshold, maximumLabels
) {
    prepared <- tryCatch(
    prepareManhattanDataDnaEpico(data, pValueColumn),
    error = function(error) NULL
    )
    if (is.null(prepared) || !nrow(prepared$data)) {
    return(NULL)
    }
    plots <- files <- list()
    rows <- list()
    for (version in c("v1", "v2")) {
    plot <- createManhattanPlotDnaEpico(
        prepared, version, suggestiveThreshold,
        genomeWideThreshold, maximumLabels
    )
    filename <- paste0(
        "manhattan_", safeFigureComponentDnaEpico(pValueColumn),
        "_", version, ".tiff"
    )
    file <- if (is.null(outputDir)) NULL else file.path(outputDir, filename)
    square <- identical(version, "v1")
    dimension <- max(plotWidth, plotHeight)
    saveModelDesignPlotDnaEpico(
        plot, file, display,
        if (square) dimension else plotWidth,
        if (square) dimension else plotHeight, plotDPI
    )
    plots[[version]] <- plot
    if (!is.null(file)) files[[version]] <- file
    rows[[version]] <- data.frame(
        analysis = analysis, plotType = "manhattan",
        pValueColumn = pValueColumn, version = version,
        filename = filename, testedCpGs = prepared$tested,
        excludedCpGs = prepared$excluded,
        suggestiveThreshold = suggestiveThreshold,
        genomeWideThreshold = genomeWideThreshold,
        stringsAsFactors = FALSE, check.names = FALSE
    )
    }
    list(plots = plots, files = files, manifest = do.call(rbind, rows))
}

plotAnnotatedManhattanDnaEpico <- function(
    annotatedResults, analysis = c("GLM", "LME"), outputDir = NULL,
    plotWidth = 2000L, plotHeight = 1000L, plotDPI = 150L,
    display = FALSE, suggestiveThreshold = 1e-5,
    genomeWideThreshold = genomeWideThresholdDnaEpico(), maximumLabels = 20L
) {
    analysis <- match.arg(analysis)
    data <- if (is.data.frame(annotatedResults)) {
    annotatedResults
    } else {
    annotatedResults$data
    }
    columns <- identifyAnnotatedPvalueColumnsDnaEpico(data)
    results <- lapply(columns, function(column) {
    plotAnnotatedManhattanColumnDnaEpico(
        data, column, analysis, outputDir, display,
        plotWidth, plotHeight, plotDPI,
        suggestiveThreshold, genomeWideThreshold, maximumLabels
    )
    })
    names(results) <- columns
    results <- Filter(Negate(is.null), results)
    manifest <- if (length(results)) {
    do.call(rbind, lapply(results, function(result) result$manifest))
    } else {
    data.frame()
    }
    structure(list(
    plots = lapply(results, function(result) result$plots),
    files = lapply(results, function(result) result$files),
    manifest = manifest, pValueColumns = columns
    ), class = "dnaEPICO_manhattan_plots")
}

# Model-level Venn helpers --------------------------------------------------

assignVennGenesDnaEpico <- function(data, columns) {
    genes <- rep(NA_character_, nrow(data))
    for (column in intersect(columns, names(data))) {
    source <- trimws(as.character(data[[column]]))
    source[is.na(source) | !nzchar(source)] <- NA_character_
    missing <- is.na(genes)
    genes[missing] <- source[missing]
    }
    genes
}

vennRegionLabelsDnaEpico <- function(sets, geneMap) {
    venn <- ggVennDiagram::Venn(sets)
    regions <- ggVennDiagram::venn_regionlabel(
    ggVennDiagram::process_data(venn)
    )
    regions$geneN <- vapply(regions$item, function(items) {
    genes <- unname(geneMap[as.character(items)])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    tokens <- trimws(unlist(strsplit(genes, ";", fixed = TRUE),
        use.names = FALSE
    ))
    length(unique(tokens[nzchar(tokens)]))
    }, integer(1))
    regions$displayLabel <- paste0(
    format(regions$count,
        big.mark = ",", scientific = FALSE,
        trim = TRUE
    ),
    " (", format(regions$geneN,
        big.mark = ",", scientific = FALSE,
        trim = TRUE
    ),
    ")"
    )
    regions
}

createVennDPlotDnaEpico <- function(sets, geneMap) {
    colours <- RColorBrewer::brewer.pal(max(3L, length(sets)), "Set2")
    colours <- colours[seq_along(sets)]
    names(colours) <- names(sets)
    category_names <- paste0(
    names(sets), "\n",
    format(lengths(sets),
        big.mark = ",", scientific = FALSE,
        trim = TRUE
    )
    )
    regions <- vennRegionLabelsDnaEpico(sets, geneMap)
    plot <- ggVennDiagram::ggVennDiagram(
    sets,
    category.names = category_names,
    set_color = colours, label = "none", edge_size = 1.05,
    set_size = 4.2
    ) +
    ggplot2::scale_fill_gradient(
        low = "#F4F7FA", high = "#176B87", trans = "log1p"
    ) +
    ggplot2::geom_text(
        data = regions,
        ggplot2::aes(x = X, y = Y, label = displayLabel),
        inherit.aes = FALSE, colour = "#111827", size = 4
    ) +
    ggplot2::labs(title = NULL) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
        plot.title = ggplot2::element_blank(),
        legend.position = "none",
        plot.margin = ggplot2::margin(18, 80, 18, 80)
    )
    plot$coordinates$clip <- "off"
    plot
}

createIntersectionPlotDnaEpico <- function(sets, geneMap,
    maximumIntersections = 25L) {
    identifiers <- sort(unique(unlist(sets, use.names = FALSE)))
    if (!length(identifiers)) { return(NULL) }
    membership <- vapply(sets, function(set) {
        identifiers %in% set
    }, logical(length(identifiers)))
    if (!is.matrix(membership)) {
        membership <- matrix(membership, ncol = length(sets))
    }; colnames(membership) <- names(sets)
    keys <- apply(membership, 1L, function(row) {
        paste0(as.integer(row), collapse = "")
    })
    counts <- sort(table(keys), decreasing = TRUE)
    counts <- utils::head(counts, max(1L, as.integer(maximumIntersections)))
    selected_keys <- names(counts)
    intersection_ids <- lapply(selected_keys, function(key) {
        identifiers[keys == key] })
    gene_counts <- vapply(intersection_ids, function(ids) {
        genes <- unname(geneMap[ids])
        tokens <- trimws(unlist(strsplit(genes[!is.na(genes) &
            nzchar(genes)], ";", fixed = TRUE), use.names = FALSE))
        length(unique(tokens[nzchar(tokens)]))
    }, integer(1))
    intersection <- paste0("I", seq_along(selected_keys))
    bars <- data.frame(intersection = intersection, count = as.integer(counts),
        genes = gene_counts, panel = factor("Intersection size",
            levels = c("Intersection size", "Set membership")),
        stringsAsFactors = FALSE)
    bars$label <- paste0(format(bars$count, big.mark = ","),
        " (", format(bars$genes, big.mark = ","), ")")
    membership_labels <- vapply(selected_keys, function(key) {
        membership <- strsplit(key, "", fixed = TRUE)[[1L]] ==
            "1"
        paste(names(sets)[membership], collapse = " + ")
    }, character(1))
    intersection_labels <- paste0(intersection, ": ", membership_labels)
    bars$intersectionLabel <- factor(intersection_labels, levels = rev(
        intersection_labels))
    ggplot2::ggplot(bars, ggplot2::aes(x = count, y = intersectionLabel)) +
        ggplot2::geom_col(fill = "#176B87", alpha = 0.88) + ggplot2::geom_text(
            ggplot2::aes(label = label),
        hjust = -0.15, size = 3.2) + ggplot2::scale_x_continuous(limits = c(0,
        max(bars$count) * 1.16), expand = ggplot2::expansion(mult = c(0,
        0.01))) + ggplot2::labs(title = NULL, x = "CpGs", y =
            "Intersection membership") +
        dnaEpicoModelPlotTheme() + ggplot2::theme(axis.text.y =
            ggplot2::element_text(size = 9),
        axis.title.y = ggplot2::element_text(size = 10), plot.margin =
            ggplot2::margin(12, 32, 12, 12)) }


normalizeModelVennOptionDnaEpico <- function(value, argument) {
    if (is.null(value) || !length(value) || all(is.na(value))) {
    return(character(0))
    }
    values <- splitOptionMinfiEwasWater(value, sep = ",")
    values <- trimws(as.character(values))
    values <- values[!is.na(values) & nzchar(values) &
    tolower(values) != "null"]
    if (anyDuplicated(values)) {
    stop(argument, " must contain unique values.", call. = FALSE)
    }
    values
}

modelVennSummaryListDnaEpico <- function(modelSummaries) {
    if (!is.null(modelSummaries$diagnosticSummaries)) {
    modelSummaries$diagnosticSummaries
    } else if (!is.null(modelSummaries$summaries)) {
    modelSummaries$summaries
    } else {
    modelSummaries
    }
}

modelVennSettingValueDnaEpico <- function(artifacts, field) {
    values <- lapply(artifacts, function(artifact) artifact$settings[[field]])
    serialized <- vapply(values, function(value) {
    paste(utils::capture.output(dput(value)), collapse = "")
    }, character(1))
    if (length(unique(serialized)) > 1L) {
    stop("Compact phenotype summaries disagree about ", field, ".",
        call. = FALSE
    )
    }
    values[[1L]]
}

modelVennFactorLevelsDnaEpico <- function(artifacts) {
    factor_levels <- list()
    for (artifact in artifacts) {
    levels <- artifact$settings$factorLevels
    if (is.list(levels) && length(levels)) {
        factor_levels[names(levels)] <- levels
    }
    }
    factor_levels
}

loadModelVennSummariesDnaEpico <- function(summaryDir,
    analysis) {
    analysis <- match.arg(toupper(analysis),
        c("GLM", "LME"))
    suffix <- if (identical(analysis, "GLM")) {
        "SummaryGLM\\.rds$"
    }
    else {
        "SummaryLME\\.rds$"
    }
    files <- list.files(summaryDir, pattern = suffix,
        full.names = TRUE)
    if (!length(files)) {
        stop("No compact ", analysis, " phenotype summaries were found in: ",
            summaryDir, call. = FALSE)
    }
    artifacts <- lapply(files, readRDS)
    valid <- vapply(artifacts, function(artifact) {
        is.list(artifact) && isTRUE(artifact$complete) &&
            identical(toupper(artifact$analysis),
                analysis) && is.character(artifact$phenotype) &&
            length(artifact$phenotype) ==
                1L && is.data.frame(artifact$targetSummary) &&
            is.list(artifact$settings)
    }, logical(1))
    if (!all(valid)) {
        stop("One or more compact ", analysis,
            " phenotype summaries are invalid.",
            call. = FALSE)
    }
    phenotypes <- vapply(artifacts, `[[`,
        character(1), "phenotype")
    if (anyDuplicated(phenotypes)) {
        stop("Compact phenotype summary names must be unique.",
            call. = FALSE)
    }
    names(artifacts) <- phenotypes
    factor_levels <- modelVennFactorLevelsDnaEpico(artifacts)
    structure(list(summaries = lapply(artifacts,
        `[[`, "targetSummary"), diagnosticSummaries = lapply(artifacts,
        `[[`, "targetSummary"), omnibusTests = lapply(artifacts,
        `[[`, "omnibusTests"), phenotypes = phenotypes,
        settings = list(interactionTerm = modelVennSettingValueDnaEpico(
            artifacts,
            "interactionTerm"), omnibusTest = isTRUE(
            modelVennSettingValueDnaEpico(artifacts,
            "omnibusTest")), factorLevels = factor_levels)),
        class = paste0("dnaEPICO_methylation",
            analysis, "_summaries"))
}

duplicatedVennCoefficientTermsDnaEpico <- function(summaries, analysis) {
    if (!identical(analysis, "GLM")) {
    return(character(0))
    }
    occurrences <- unlist(lapply(summaries, function(table) {
    if (!is.data.frame(table) || !"Coefficient" %in% names(table)) {
        return(character(0))
    }
    unique(as.character(table$Coefficient))
    }), use.names = FALSE)
    unique(occurrences[duplicated(occurrences)])
}

coefficientVennColumnDnaEpico <- function(
    term, phenotype, analysis, duplicatedTerms
) {
    clean_term <- gsub("`", "", term, fixed = TRUE)
    if (identical(analysis, "GLM")) {
    if (term %in% duplicatedTerms) {
        return(paste0(phenotype, "_", clean_term, "_P.Value"))
    }
    return(paste0(clean_term, "P.Value"))
    }
    suffix <- clean_term
    if (startsWith(clean_term, paste0(phenotype, "."))) {
    suffix <- sub(paste0(
        "^", escapeRegexMethylationGLM(phenotype), "\\."
    ), "", clean_term)
    }
    paste0(phenotype, "_", suffix, "_P.Value")
}

coefficientVennMappingRowDnaEpico <- function(
    summary, phenotype, analysis, duplicatedTerms
) {
    term_column <- if (identical(analysis, "GLM")) {
    "Coefficient"
    } else {
    "Interaction.Term"
    }
    if (!is.data.frame(summary) || !term_column %in% names(summary)) {
    stop(
        "No coefficient summaries were available for Venn phenotype '",
        phenotype, "'.",
        call. = FALSE
    )
    }
    terms <- unique(as.character(summary[[term_column]]))
    terms <- terms[!is.na(terms) & nzchar(terms)]
    if (!length(terms)) {
    stop(
        "No coefficient terms were available for Venn phenotype '",
        phenotype, "'.",
        call. = FALSE
    )
    }
    columns <- vapply(terms, coefficientVennColumnDnaEpico, character(1),
    phenotype = phenotype, analysis = analysis,
    duplicatedTerms = duplicatedTerms
    )
    data.frame(
    phenotype = phenotype, term = terms, pValueColumn = columns,
    originalLabel = sub("_?P\\.Value$", "", columns),
    stringsAsFactors = FALSE, check.names = FALSE
    )
}

resolveCoefficientVennColumnsDnaEpico <- function(
    annotatedData, modelSummaries, phenotypes, analysis = c("GLM", "LME")
) {
    analysis <- match.arg(analysis)
    summaries <- modelVennSummaryListDnaEpico(modelSummaries)
    missing_phenotypes <- setdiff(phenotypes, names(summaries))
    if (length(missing_phenotypes)) {
    missing_phenotypes_text <- paste(missing_phenotypes, collapse = ", ")
    stop(
        sprintf(
        "Venn phenotype(s) were not modelled: %s",
        missing_phenotypes_text
        ),
        call. = FALSE
    )
    }

    duplicated_terms <- duplicatedVennCoefficientTermsDnaEpico(
    summaries, analysis
    )

    rows <- lapply(phenotypes, function(phenotype) {
    coefficientVennMappingRowDnaEpico(
        summaries[[phenotype]], phenotype, analysis, duplicated_terms
    )
    })
    mapping <- do.call(rbind, rows)
    mapping <- mapping[!duplicated(mapping$pValueColumn), , drop = FALSE]
    rownames(mapping) <- NULL
    missing_columns <- setdiff(mapping$pValueColumn, names(annotatedData))
    if (length(missing_columns)) {
    missing_columns_text <- paste(missing_columns, collapse = ", ")
    message_prefix <- paste(
        "Resolved Venn coefficient P-value column(s) were not found",
        "in the annotated results"
    )
    stop(
        sprintf(
        "%s: %s",
        message_prefix, missing_columns_text
        ),
        call. = FALSE
    )
    }
    mapping
}

resolveOmnibusVennColumnsDnaEpico <- function(annotatedData,
    modelSummaries, phenotypes) {
    if (!isTRUE(modelSummaries$settings$omnibusTest)) {
        stop("Omnibus Venn phenotypes require omnibusTest = TRUE.",
            call. = FALSE)
    }
    model_phenotypes <- modelSummaries$phenotypes
    if (is.null(model_phenotypes)) {
        model_phenotypes <- names(modelVennSummaryListDnaEpico(modelSummaries))
    }
    missing_phenotypes <- setdiff(phenotypes,
        model_phenotypes)
    if (length(missing_phenotypes)) {
        missing_phenotypes_text <- paste(missing_phenotypes,
            collapse = ", ")
        stop(sprintf("Omnibus Venn phenotype(s) were not modelled: %s",
            missing_phenotypes_text), call. = FALSE)
    }
    interaction <- modelSummaries$settings$interactionTerm
    prefixes <- vapply(phenotypes, function(phenotype) {
        paste(c(phenotype, if (!is.null(interaction) &&
            nzchar(interaction)) {
            interaction
        } else {
            NULL
        }), collapse = "_")
    }, character(1))
    columns <- paste0(gsub("`", "", prefixes,
        fixed = TRUE), "_Omnibus_P.Value")
    missing_columns <- setdiff(columns, names(annotatedData))
    if (length(missing_columns)) {
        missing_columns_text <- paste(missing_columns,
            collapse = ", ")
        message_prefix <- paste(
            "Resolved omnibus Venn P-value column(s) were not found in",
            "the annotated results")
        stop(sprintf("%s: %s", message_prefix,
            missing_columns_text), call. = FALSE)
    }
    data.frame(phenotype = phenotypes, term = prefixes,
        pValueColumn = columns, originalLabel = phenotypes,
        stringsAsFactors = FALSE, check.names = FALSE)
}

applyModelVennLabelsDnaEpico <- function(mapping, labels, argument) {
    labels <- normalizeModelVennOptionDnaEpico(labels, argument)
    if (!length(labels)) {
    mapping$displayLabel <- mapping$originalLabel
    mapping$labelSource <- "automatic"
    return(mapping)
    }
    if (length(labels) != nrow(mapping)) {
    original_labels_text <- paste(mapping$originalLabel, collapse = ", ")
    stop(
        sprintf(
        "%s must contain exactly %s labels in this %s: %s",
        argument, nrow(mapping), "resolved coefficient order",
        original_labels_text
        ),
        call. = FALSE
    )
    }
    if (anyDuplicated(tolower(labels))) {
    stop(argument, " labels must be unique without regard to letter case.",
        call. = FALSE
    )
    }
    mapping$displayLabel <- labels
    mapping$labelSource <- "user"
    mapping
}

validateResolvedModelVennSetsDnaEpico <- function(mapping, type) {
    n_sets <- nrow(mapping)
    if (n_sets < 2L) {
    stop(type, " Venn analysis resolved fewer than two P-value columns.",
        call. = FALSE
    )
    }
    if (n_sets > 7L) {
    p_value_columns_text <- paste(mapping$pValueColumn, collapse = ", ")
    stop(
        sprintf(
        "%s %s: %s. %s",
        type, "Venn analysis resolved more than seven P-value columns",
        p_value_columns_text, "No columns were omitted."
        ),
        call. = FALSE
    )
    }
    invisible(mapping)
}

modelVennRegionSpecificationDnaEpico <- function(labels) {
    n_sets <- length(labels)
    rows <- list()
    row_index <- 1L
    for (size in seq_len(n_sets)) {
    combinations <- utils::combn(seq_len(n_sets), size, simplify = FALSE)
    for (indices in combinations) {
        key <- paste0(as.integer(seq_len(n_sets) %in% indices),
        collapse = ""
        )
        label <- paste(labels[indices], collapse = "_AND_")
        if (size < n_sets) {
        label <- paste0(label, "_ONLY")
        }
        rows[[row_index]] <- data.frame(
        key = key,
        column = safeFigureComponentDnaEpico(label, "Venn_region"),
        description = paste(
            "CpG significant for", paste(labels[indices],
            collapse = ", "
            ),
            if (size < n_sets) "only" else "in all selected sets"
        ), stringsAsFactors = FALSE, check.names = FALSE
        )
        row_index <- row_index + 1L
    }
    }
    result <- do.call(rbind, rows)
    if (anyDuplicated(tolower(result$column))) {
    stop(
        "Venn display labels create duplicate workbook region columns ",
        "after filename-safe normalization.",
        call. = FALSE
    )
    }
    result
}

buildModelVennWorksheetDnaEpico <- function(data,
    mapping, threshold) {
    values <- lapply(mapping$pValueColumn,
        function(column) {
            coerceNumericDnaEpico(data[[column]])
        })
    significant <- do.call(cbind, lapply(values,
        function(value) {
            is.finite(value) & value < threshold
        }))
    if (!is.matrix(significant)) {
        significant <- matrix(significant,
            ncol = nrow(mapping))
    }
    colnames(significant) <- mapping$displayLabel
    keep <- rowSums(significant) > 0L
    specifications <- modelVennRegionSpecificationDnaEpico(mapping$displayLabel)
    keys <- apply(significant[keep, , drop = FALSE],
        1L, function(row) {
            paste0(as.integer(row), collapse = "")
        })
    region_index <- match(keys, specifications$key)
    membership <- matrix(NA_character_, nrow = sum(keep),
        ncol = nrow(specifications), dimnames = list(NULL,
            specifications$column))
    identifiers <- as.character(data$IlmnID[keep])
    if (length(identifiers)) {
        membership[cbind(seq_along(identifiers),
            region_index)] <- identifiers
    }
    all_pvalues <- names(data)[grepl("P\\.Value$|P\\.value$",
        names(data))]
    omnibus_details <- names(data)[grepl("_Omnibus_",
        names(data), fixed = TRUE)]
    model_messages <- names(data)[grepl("_Model\\.Message$",
        names(data))]
    supporting <- setdiff(names(data), c("IlmnID",
        "CpG", all_pvalues, omnibus_details,
        model_messages))
    selected <- data[keep, c("IlmnID", mapping$pValueColumn,
        supporting), drop = FALSE]
    result <- cbind(as.data.frame(membership,
        check.names = FALSE, stringsAsFactors = FALSE),
        selected)
    list(data = result, regions = specifications)
}

modelVennMetadataRowsDnaEpico <- function(configurations, modelSummaries) {
    rows <- list()
    row_index <- 1L
    for (type in names(configurations)) {
    mapping <- configurations[[type]]
    for (index in seq_len(nrow(mapping))) {
        phenotype <- mapping$phenotype[[index]]
        factor_levels <- modelSummaries$settings$factorLevels[[phenotype]]
        reference <- if (length(factor_levels)) {
        factor_levels[[1L]]
        } else {
        "None"
        }
        prefix <- paste0("vennD.", type, ".label.", index, ".")
        entries <- c(
        phenotype = phenotype, term = mapping$term[[index]],
        p_value_column = mapping$pValueColumn[[index]],
        original = mapping$originalLabel[[index]],
        display = mapping$displayLabel[[index]],
        source = mapping$labelSource[[index]], reference = reference
        )
        for (field in names(entries)) {
        rows[[row_index]] <- data.frame(
            Key = paste0(prefix, field), Value = entries[[field]],
            stringsAsFactors = FALSE, check.names = FALSE
        )
        row_index <- row_index + 1L
        }
    }
    }
    if (!length(rows)) data.frame() else do.call(rbind, rows)
}

emptyModelVennResultDnaEpico <- function() {
    structure(list(
    plots = list(), files = data.frame(), sheets = list(),
    dictionaryRows = data.frame(), metadataRows = data.frame(),
    mappings = list()
    ), class = "dnaEPICO_vennD_plots")
}

removeOldModelVennFiguresDnaEpico <- function(outputDir, analysis) {
    if (is.null(outputDir) || !dir.exists(outputDir)) {
    return(invisible(NULL))
    }
    old_files <- list.files(
    outputDir,
    pattern = paste0("^(vennD|intersection)_", analysis, "_.*\\.tiff$"),
    full.names = TRUE, ignore.case = TRUE
    )
    if (length(old_files)) {
    status <- unlink(old_files, force = TRUE)
    if (status != 0L || any(file.exists(old_files))) {
        stop(
        "Could not replace existing ", analysis,
        " Venn figures. Close open figure files and retry.",
        call. = FALSE
        )
    }
    }
    invisible(NULL)
}

resolveModelVennAnnotationDnaEpico <- function(annotatedResults) {
    data <- if (is.data.frame(annotatedResults)) {
    annotatedResults
    } else {
    annotatedResults$data
    }
    if (!is.data.frame(data) || !"IlmnID" %in% names(data)) {
    stop("Annotated Venn results must contain an IlmnID column.",
        call. = FALSE
    )
    }
    direct <- grep(
    "^GencodeV[0-9]+_RefGene_Name$", names(data),
    value = TRUE
    )
    nearest <- grep(
    "^GencodeV[0-9]+_NonAnnotated_RefGene_Name$",
    names(data),
    value = TRUE
    )
    if (length(direct) != 1L || length(nearest) != 1L) {
    stop(
        "Venn figures require one release-matched pair of GENCODE ",
        "direct and nearest gene-name columns.",
        call. = FALSE
    )
    }
    release <- sub("^GencodeV([0-9]+)_RefGene_Name$", "\\1", direct)
    expected_nearest <- paste0(
    "GencodeV", release, "_NonAnnotated_RefGene_Name"
    )
    required <- c("UCSC_RefGene_Name", direct, expected_nearest)
    missing <- setdiff(required, names(data))
    if (length(missing)) {
    requirement <- paste(
        "Venn figures require UCSC and release-matched GENCODE gene",
        "annotation columns"
    )
    missing_text <- paste(missing, collapse = ", ")
    stop(sprintf(
        "%s. Missing: %s",
        requirement, missing_text
    ), call. = FALSE)
    }
    annotations <- list(UCSC = "UCSC_RefGene_Name")
    annotations[[paste0("GENCODEv", release)]] <- c(direct, expected_nearest)
    list(data = data, release = release, annotations = annotations)
}

resolveModelVennConfigurationsDnaEpico <- function(
    data, modelSummaries, analysis,
    coefficientPhenotypes, coefficientLabels,
    omnibusPhenotypes, omnibusLabels
) {
    configurations <- list()
    if (length(coefficientPhenotypes)) {
    mapping <- resolveCoefficientVennColumnsDnaEpico(
        data, modelSummaries, coefficientPhenotypes, analysis
    )
    mapping <- applyModelVennLabelsDnaEpico(
        mapping, coefficientLabels, "vennDLabels"
    )
    validateResolvedModelVennSetsDnaEpico(mapping, "Coefficient")
    configurations$coefficient <- mapping
    }
    if (length(omnibusPhenotypes)) {
    mapping <- resolveOmnibusVennColumnsDnaEpico(
        data, modelSummaries, omnibusPhenotypes
    )
    mapping <- applyModelVennLabelsDnaEpico(
        mapping, omnibusLabels, "vennDOmnibusLabels"
    )
    validateResolvedModelVennSetsDnaEpico(mapping, "Omnibus")
    configurations$omnibus <- mapping
    }
    configurations
}

modelVennThresholdsDnaEpico <- function() {
    c(
    nominal = 0.05, suggestive = 1e-5,
    genomeWide = genomeWideThresholdDnaEpico()
    )
}

modelVennSheetNamesDnaEpico <- function() {
    list(
    coefficient = c(
        nominal = "vennDNOM", suggestive = "vennDSUGG",
        genomeWide = "vennDGW"
    ),
    omnibus = c(
        nominal = "vennDOmnibusNOM",
        suggestive = "vennDOmnibusSUGG",
        genomeWide = "vennDOmnibusGW"
    )
    )
}

buildModelVennWorksheetsDnaEpico <- function(
    data, configurations, thresholds, sheetNames
) {
    sheets <- list()
    dictionary <- list()
    index <- 1L
    for (type in names(configurations)) {
    mapping <- configurations[[type]]
    for (threshold_name in names(thresholds)) {
        threshold <- thresholds[[threshold_name]]
        worksheet <- buildModelVennWorksheetDnaEpico(
        data, mapping, threshold
        )
        sheet <- unname(sheetNames[[type]][[threshold_name]])
        sheets[[sheet]] <- worksheet$data
        note <- paste0(
        "The applicable threshold for ", sheet, " is p < ",
        format(threshold, scientific = TRUE), "."
        )
        dictionary[[index]] <- data.frame(
        Column = worksheet$regions$column,
        Description = paste(worksheet$regions$description, note),
        Formula = "", stringsAsFactors = FALSE,
        check.names = FALSE
        )
        index <- index + 1L
    }
    }
    dictionary <- if (length(dictionary)) {
    unique(do.call(rbind, dictionary))
    } else {
    data.frame()
    }
    list(sheets = sheets, dictionary = dictionary)
}

drawModelVennFigureDnaEpico <- function(
    plot, file, display, width, height, resolution
) {
    runPlotMinfiEwasWater(
    draw_fun = function() drawPlotObjectMinfiEwasWater(plot),
    display = display, file = file, width = width,
    height = height, res = resolution
    )
    invisible(NULL)
}

modelVennFileRowDnaEpico <- function(
    analysis, type, figureType, mapping, threshold, annotation,
    filename, path
) {
    data.frame(
    analysis = analysis, type = type, figureType = figureType,
    pValueColumns = paste(mapping$pValueColumn, collapse = ", "),
    labels = paste(mapping$displayLabel, collapse = ", "),
    threshold = threshold, annotation = annotation,
    filename = filename,
    path = if (is.null(path)) NA_character_ else path,
    stringsAsFactors = FALSE, check.names = FALSE
    )
}

baseModelVennPlotPairDnaEpico <- function(
    data, ids, mapping, geneMap, context
) {
    values <- lapply(mapping$pValueColumn, function(column) {
        coerceNumericDnaEpico(data[[column]])
    })
    sets <- lapply(values, function(value) {
        unique(ids[is.finite(value) & value < context$threshold])
    })
    names(sets) <- mapping$displayLabel
    type_tag <- if (identical(context$type, "omnibus")) "_omnibus" else ""
    filename <- paste0(
        "vennD_", context$analysis, "_", context$labelComponent,
        type_tag, "_", context$thresholdName, "_",
        context$annotationName, ".tiff"
    )
    file <- if (is.null(context$outputDir)) {
        NULL
    } else {
        file.path(context$outputDir, filename)
    }
    plot <- createVennDPlotDnaEpico(sets, geneMap)
    drawModelVennFigureDnaEpico(
        plot, file, context$display, context$width,
        context$height, context$resolution
    )
    key <- paste(
        context$type, context$thresholdName,
        context$annotationName, sep = ":"
    )
    row <- modelVennFileRowDnaEpico(
        context$analysis, context$type, "vennD", mapping,
        context$thresholdName, context$annotationName, filename, file
    )
    list(
        sets = sets, typeTag = type_tag,
        plots = stats::setNames(list(plot), key), rows = list(row)
    )
}

addModelVennIntersectionDnaEpico <- function(result, mapping, geneMap,
                                                context) {
    if (nrow(mapping) < 4L) return(result)
    plot <- createIntersectionPlotDnaEpico(result$sets, geneMap)
    if (is.null(plot)) return(result)
    filename <- paste0(
        "intersection_", context$analysis, "_", context$labelComponent,
        result$typeTag, "_", context$thresholdName, "_",
        context$annotationName, ".tiff"
    )
    file <- if (is.null(context$outputDir)) {
        NULL
    } else {
        file.path(context$outputDir, filename)
    }
    drawModelVennFigureDnaEpico(
        plot, file, context$display, max(context$width, 2200L),
        max(context$height, 1800L), context$resolution
    )
    key <- paste(
        "intersection", context$type, context$thresholdName,
        context$annotationName, sep = ":"
    )
    result$plots[[key]] <- plot
    result$rows[[2L]] <- modelVennFileRowDnaEpico(
        context$analysis, context$type, "intersection", mapping,
        context$thresholdName, context$annotationName, filename, file
    )
    result
}

buildModelVennPlotPairDnaEpico <- function(
    data, ids, mapping, geneMap, analysis, type, thresholdName,
    threshold, annotationName, labelComponent, outputDir, display,
    width, height, resolution
) {
    context <- list(
        analysis = analysis, type = type, thresholdName = thresholdName,
        threshold = threshold, annotationName = annotationName,
        labelComponent = labelComponent, outputDir = outputDir,
        display = display, width = width, height = height,
        resolution = resolution
    )
    result <- baseModelVennPlotPairDnaEpico(
        data, ids, mapping, geneMap, context
    )
    result <- addModelVennIntersectionDnaEpico(
        result, mapping, geneMap, context
    )
    result[c("plots", "rows")]
}

buildModelVennPlotsDnaEpico <- function(
    data, configurations, annotations, thresholds, analysis,
    outputDir, display, width, height, resolution
) {
    plots <- rows <- list()
    ids <- as.character(data$IlmnID)
    for (type in names(configurations)) {
    mapping <- configurations[[type]]
    label_component <- paste(vapply(
        mapping$displayLabel, safeFigureComponentDnaEpico, character(1)
    ), collapse = "_")
    for (annotation_name in names(annotations)) {
        gene_map <- stats::setNames(assignVennGenesDnaEpico(
        data, annotations[[annotation_name]]
        ), ids)
        for (threshold_name in c("nominal", "genomeWide")) {
        result <- buildModelVennPlotPairDnaEpico(
            data, ids, mapping, gene_map, analysis, type,
            threshold_name, thresholds[[threshold_name]],
            annotation_name, label_component,
            outputDir, display, width, height, resolution
        )
        plots <- c(plots, result$plots)
        rows <- c(rows, result$rows)
        }
    }
    }
    files <- if (length(rows)) do.call(rbind, rows) else data.frame()
    list(plots = plots, files = files)
}

generateModelVennDDnaEpico <- function(annotatedResults, modelSummaries,
    analysis = c("GLM", "LME"), vennDPhenotypes = NULL, vennDLabels = NULL,
    vennDOmnibusPhenotypes = NULL, vennDOmnibusLabels = NULL, outputDir = NULL,
    plotWidth = 1800L, plotHeight = 1600L, plotDPI = 180L, display = FALSE,
    verbose = FALSE, logs = FALSE, log_dir = NULL, log_file = NULL) {
    analysis <- match.arg(analysis)
    coefficient_phenotypes <- normalizeModelVennOptionDnaEpico(vennDPhenotypes,
        "vennDPhenotypes")
    omnibus_phenotypes <- normalizeModelVennOptionDnaEpico(
        vennDOmnibusPhenotypes,
        "vennDOmnibusPhenotypes")
    removeOldModelVennFiguresDnaEpico(outputDir, analysis)
    if (!length(coefficient_phenotypes) && !length(omnibus_phenotypes)) {
        return(emptyModelVennResultDnaEpico()) }
    if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
        stop("Package 'ggVennDiagram' is required when Venn phenotypes are ",
            "supplied.", call. = FALSE) }
    annotation <- resolveModelVennAnnotationDnaEpico(annotatedResults)
    configurations <- resolveModelVennConfigurationsDnaEpico(annotation$data,
        modelSummaries, analysis, coefficient_phenotypes, vennDLabels,
        omnibus_phenotypes, vennDOmnibusLabels)
    log_path <- resolveLogPathMinfiEwasWater(logs = logs, log_dir = log_dir,
        log_file = log_file)
    emitLogMinfiEwasWater(c(
        "============================================================",
        paste("Starting", analysis, "model Venn analysis"), paste(
            "Venn result types:",
            paste(names(configurations), collapse = ", ")), paste(
            "Venn output directory:",
            if (is.null(outputDir)) {
                "In memory only"
            } else { outputDir
            })), verbose = verbose, log_path = log_path)
    thresholds <- modelVennThresholdsDnaEpico()
    worksheets <- buildModelVennWorksheetsDnaEpico(annotation$data,
        configurations, thresholds, modelVennSheetNamesDnaEpico())
    figures <- buildModelVennPlotsDnaEpico(annotation$data, configurations,
        annotation$annotations, thresholds, analysis, outputDir,
        display, plotWidth, plotHeight, plotDPI)
    metadata <- modelVennMetadataRowsDnaEpico(configurations, modelSummaries)
    emitLogMinfiEwasWater(c(paste("Completed", analysis, "model Venn analysis"),
        paste("Generated Venn figures:", nrow(figures$files)),
        paste("Generated Venn worksheets:", paste(names(worksheets$sheets),
            collapse = ", ")),
            "============================================================"),
        verbose = verbose, log_path = log_path)
    structure(list(plots = figures$plots, files = figures$files,
        sheets = worksheets$sheets, dictionaryRows = worksheets$dictionary,
        metadataRows = metadata, mappings = configurations), class =
            "dnaEPICO_vennD_plots") }
