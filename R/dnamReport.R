#' Normalize a report path for cross-platform rendering
#'
#' @param path Character. Path to normalize.
#'
#' @return Character scalar containing a normalized path.
#'
#' @keywords internal
#' @noRd
normalizePathDnamReport <- function(path) {
    gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = FALSE))
}

#' Normalize requested model sections for dashboard reports
#'
#' @param modelSections Character vector containing any of `'glm'` and `'lme'`.
#'
#' @return A unique character vector in canonical report order.
#'
#' @keywords internal
#' @noRd
normalizeModelSectionsDnamReport <- function(modelSections) {
    if (!is.character(modelSections)) {
        stop("modelSections must be a character vector.", call. = FALSE)
    }
    if (length(modelSections) == 0L) {
        return(character(0))
    }
    if (anyNA(modelSections) || any(!nzchar(trimws(modelSections)))) {
        stop("modelSections cannot contain missing or empty values.",
            call. = FALSE
        )
    }

    allowed <- c("glm", "lme")
    normalized <- unique(tolower(trimws(modelSections)))
    invalid <- setdiff(normalized, allowed)
    if (length(invalid) > 0L) {
        stop("Unsupported modelSections value(s): ",
            paste(invalid, collapse = ", "), ". Expected any of: ",
            paste(allowed, collapse = ", "), ".",
            call. = FALSE
        )
    }
    allowed[allowed %in% normalized]
}

#' Format a report print path
#'
#' @param path Character path to format.
#' @param missing Character value used when the path is unavailable.
#'
#' @return Character scalar for console printing.
#'
#' @keywords internal
#' @noRd
formatPrintPathDnamReport <- function(path, missing = "not written") {
    if (is.null(path) || !length(path) || is.na(path[[1]]) ||
        !nzchar(path[[1]])) {
        return(missing)
    }

    normalizePathDnamReport(path[[1]])
}

#' Collect report figure inventory from a directory
#'
#' @param directory Character. Directory containing report figures.
#' @param patterns Character vector of regular expressions used to find matching
#'   files.
#' @param label Character. Short label used in report logs.
#' @param recursive Logical. If `TRUE`, search subdirectories as well.
#'
#' @return A list containing the normalized directory, whether it exists, the
#'   matching files, the number of matches, and the label.
#'
#' @keywords internal
#' @noRd
collectFigureInventoryDnamReport <- function(
    directory, patterns,
    label, recursive = TRUE
) {
    normalized_dir <- normalizePathDnamReport(directory)
    files <- character(0)
    directory_exists <- dir.exists(directory)

    if (directory_exists) {
        file_sets <- lapply(patterns, function(pattern) {
            list.files(directory,
                pattern = pattern, full.names = TRUE,
                recursive = recursive, ignore.case = TRUE
            )
        })
        files <- sort(unique(unlist(file_sets, use.names = FALSE)))
        files <- vapply(files, normalizePathDnamReport, character(1))
    }

    list(
        label = label, directory = normalized_dir, exists = directory_exists,
        files = files, count = length(files)
    )
}

#' Prepare inputs for a DNA methylation report
#'
#' @param outputDir Character. Directory where the report project is written.
#' @param qcDir Character. Directory containing ENmix quality-control figures.
#' @param preprocessingDir Character. Directory containing preprocessing
#'   quality-control figures.
#' @param postprocessingDir Character. Directory containing postprocessing
#' metric
#'   figures.
#' @param svaDir Character. Directory containing SVA or batch-effect figures.
#' @param glmDir Character. Directory containing GLM figures.
#' @param lmeDir Character. Directory containing LME figures.
#' @param figDir Character. Directory used for generated report figure assets.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write progress messages to
#'   `file.path(logDir, 'log_dnamReport.txt')`.
#' @param logDir Character. Directory for optional log files.
#'
#' @return A list with class `'dnaEPICO_dnamReport_prepared'`.
#'
#' @examples
#' report_root <- file.path(tempdir(), "dnaepico-report-inputs")
#' prepared <- prepareDnamReportInputs(
#'     outputDir = file.path(report_root, "reports"),
#'     qcDir = file.path(
#'         report_root,
#'         "figures",
#'         "preprocessingMinfiEwasWater",
#'         "enmix"
#'     ),
#'     preprocessingDir = file.path(
#'         report_root,
#'         "figures",
#'         "preprocessingMinfiEwasWater",
#'         "qc"
#'     ),
#'     postprocessingDir = file.path(
#'         report_root,
#'         "figures",
#'         "preprocessingMinfiEwasWater",
#'         "metrics"
#'     ),
#'     svaDir = file.path(report_root, "figures", "svaEnmix")
#' )
#' inherits(prepared, "dnaEPICO_dnamReport_prepared")
#'
#' @export
prepareDnamReportInputs <- function(
    outputDir = "reports", qcDir = file.path(
        "figures",
        "preprocessingMinfiEwasWater", "enmix"
    ), preprocessingDir = file.path(
        "figures",
        "preprocessingMinfiEwasWater", "qc"
    ), postprocessingDir = file.path(
        "figures",
        "preprocessingMinfiEwasWater", "metrics"
    ), svaDir = file.path(
        "figures",
        "svaEnmix"
    ), glmDir = file.path("figures", "methylationGLM"),
    lmeDir = file.path("figures", "methylationLME"), figDir = file.path(
        outputDir,
        "assets", "figures"
    ), verbose = FALSE, logs = FALSE,
    logDir = outputDir
) {
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = logDir,
        log_file = "log_dnamReport.txt"
    )

    figure_inventory <- list(
        qc = collectFigureInventoryDnamReport(
            directory = qcDir,
            patterns = c(
                "\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$",
                "\\.tif$", "\\.tiff$"
            ), label = "ENmix QC"
        ), preprocessing = collectFigureInventoryDnamReport(
            directory = preprocessingDir,
            patterns = c(
                "\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$",
                "\\.tif$", "\\.tiff$"
            ), label = "Quality control"
        ),
        postprocessing = collectFigureInventoryDnamReport(
            directory = postprocessingDir,
            patterns = c(
                "\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$",
                "\\.tif$", "\\.tiff$"
            ), label = "Postprocessing"
        ),
        sva = collectFigureInventoryDnamReport(
            directory = svaDir,
            patterns = c(
                "\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$",
                "\\.tif$", "\\.tiff$"
            ), label = "SVA"
        ), glm = collectFigureInventoryDnamReport(
            directory = glmDir,
            patterns = c(
                "\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$",
                "\\.tif$", "\\.tiff$"
            ), label = "GLM"
        ), lme = collectFigureInventoryDnamReport(
            directory = lmeDir,
            patterns = c(
                "\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$",
                "\\.tif$", "\\.tiff$"
            ), label = "LME"
        )
    )

    missing_directories <- vapply(
        figure_inventory, function(section) !isTRUE(section$exists),
        logical(1)
    )
    missing_directories <- names(missing_directories)[missing_directories]

    emitLogMinfiEwasWater(
        c(
            "============================================================",
            "Prepared DNA methylation dashboard report inputs", paste(
                "Output directory:",
                normalizePathDnamReport(outputDir)
            ), paste(
                "Figure directory:",
                normalizePathDnamReport(figDir)
            ), "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    output_file <- normalizePathDnamReport(file.path(
        outputDir,
        "docs", "index.html"
    ))

    structure(list(
        output = basename(output_file),
            outputDir = normalizePathDnamReport(outputDir),
        outputFile = output_file, figDir = normalizePathDnamReport(figDir),
        figureInventory = figure_inventory,
            missingFigureDirectories = missing_directories,
        logFile = log_path
    ), class = "dnaEPICO_dnamReport_prepared")
}

#' Render a prepared DNA methylation report
#'
#' @param preparedReport Object returned by `prepareDnamReportInputs()`.
#' @param verbose Logical. If `TRUE`, emit progress messages.
#' @param logs Logical. If `TRUE`, write progress messages to a log file.
#' @param logDir Character or `NULL`. Directory for optional log files.
#' @param clean Logical. Retained for backwards compatibility.
#'
#' @return A list with class `'dnaEPICO_dnamReport_render'`.
#'
#' @examples
#' report_root <- file.path(tempdir(), "dnaepico-render-example")
#' prepared <- prepareDnamReportInputs(
#'     outputDir = file.path(report_root, "reports")
#' )
#' rendered <- renderDnamReport(prepared)
#' rendered$status
#'
#' @export
renderDnamReport <- function(
    preparedReport, verbose = FALSE,
    logs = FALSE, logDir = NULL, clean = TRUE
) {
    if (!inherits(preparedReport, "dnaEPICO_dnamReport_prepared")) {
        stop("preparedReport must be an object returned by prepareDnamReportInputs().",
            call. = FALSE
        )
    }

    log_path <- resolveLogPathMinfiEwasWater(logs = logs,
        log_dir = if (is.null(logDir)) {
        dirname(preparedReport$outputFile)
    } else {
        logDir
    }, log_file = "log_dnamReport.txt")

    emitLogMinfiEwasWater(
        c(
            "============================================================",
            "renderDnamReport() is retained for compatibility.",
            "Use dnamReport() to generate the Quarto dashboard report.",
            "============================================================"
        ),
        verbose = verbose, log_path = log_path
    )

    structure(list(
        preparedReport = preparedReport, status = "skipped",
        renderedFile = NULL,
            errorMessage = "Use dnamReport() to generate the Quarto dashboard report.",
        logFile = log_path
    ), class = "dnaEPICO_dnamReport_render")
}

#' Generate a DNA methylation dashboard report
#'
#' @param outputDir Character. Directory where the Quarto project is written.
#' @param phenoTab Character or `NULL`. CSV file shown in the Data tab.
#'   When `NULL`, the path is inferred from the Makefile output layout.
#' @param enmixTab Character. Directory containing ENmix quality-control
#' figures.
#' @param qcTab Character. Directory containing Quality Control figures.
#' @param svaTab Character. Directory containing Batch Effect or SVA figures.
#' @param metricTab Character. Directory containing Metrics figures.
#' @param glmTab Character or `NULL`. XLSX workbook shown in the GLM Analysis
#' tab.
#'   When `NULL`, the path is inferred from the Makefile output layout. Report
#'   sidecars produced by `methylationGLM()` are read from this report project's
#'   `assets/results/glm_results` directory when available.
#' @param lmeTab Character or `NULL`. XLSX workbook shown in the LME Analysis
#' tab.
#'   When `NULL`, the path is inferred from the Makefile output layout. Report
#'   sidecars produced by `methylationLME()` are read from this report project's
#'   `assets/results/lme_results` directory when available.
#' @param modelSections Character vector containing any of `'glm'` and `'lme'`.
#'   The report includes the corresponding model pages, logs, and summary
#'   sections.
#'   Use `character(0)` for a preprocessing-only report. The default preserves
#'   the complete GLM-and-LME report.
#' @param logTab Character. Directory containing workflow logs shown in the Logs
#'   tab.
#' @param verbose Logical. If `TRUE`, emit progress messages.
#' @param logs Logical. If `TRUE`, write a report log.
#' @param projectName Character. Name used for the generated Quarto project.
#' @param detPPath Character or `NULL`. RData file containing the detection
#'   P-value matrix object `detP`, used to build the quality-control tables.
#'   When `NULL`, the path is inferred from the Makefile output layout.
#' @param detPThreshold Numeric. Detection P-value threshold used when
#'   summarising the `detP` matrix.
#' @param cpgDetectionPath Character or `NULL`. Optional fallback CpG detection
#'   summary CSV.
#' @param sampleDetectionPath Character or `NULL`. Optional fallback sample
#'   detection summary CSV.
#' @param logoPath Character. Path to the navbar logo. Defaults to the packaged
#'   `inst/extdata/dnaEPICO.svg` asset.
#' @param imagePattern Character. Regular expression used to identify image
#'   files inside the section directories.
#' @param recursive Logical. If `TRUE`, search section directories recursively.
#'
#' @return A list with class `'dnaEPICO_dnamReport'`.
#'
#' @examples
#' report_root <- file.path(tempdir(), "dnaepico-dnam-report")
#' pheno_file <- file.path(
#'     report_root,
#'     "data",
#'     "model1",
#'     "preprocessingMinfiEwasWater",
#'     "phenoLC.csv"
#' )
#' dir.create(dirname(pheno_file), recursive = TRUE, showWarnings = FALSE)
#' utils::write.csv(
#'     data.frame(
#'         UID = c("sample1", "sample2"),
#'         Timepoint = c(1, 2),
#'         Sex = c("F", "M")
#'     ),
#'     pheno_file,
#'     row.names = FALSE
#' )
#'
#' result <- dnamReport(
#'     outputDir = file.path(report_root, "reports", "model1"),
#'     phenoTab = pheno_file,
#'     enmixTab = file.path(
#'         report_root,
#'         "figures",
#'         "model1",
#'         "preprocessingMinfiEwasWater",
#'         "enmix"
#'     ),
#'     qcTab = file.path(
#'         report_root,
#'         "figures",
#'         "model1",
#'         "preprocessingMinfiEwasWater",
#'         "qc"
#'     ),
#'     svaTab = file.path(report_root, "figures", "model1", "svaEnmix"),
#'     metricTab = file.path(
#'         report_root,
#'         "figures",
#'         "model1",
#'         "preprocessingMinfiEwasWater",
#'         "metrics"
#'     ),
#'     logTab = file.path(report_root, "logs", "model1")
#' )
#' result$status
#'
#' @export
dnamReport <- function(
    outputDir = "reports", phenoTab = NULL,
    enmixTab = file.path(
        "figures", "preprocessingMinfiEwasWater",
        "enmix"
    ), qcTab = file.path(
        "figures", "preprocessingMinfiEwasWater",
        "qc"
    ), svaTab = file.path("figures", "svaEnmix"), metricTab = file.path(
        "figures",
        "preprocessingMinfiEwasWater", "metrics"
    ), glmTab = NULL,
    lmeTab = NULL, modelSections = c("glm", "lme"),
    logTab = outputDir, verbose = FALSE, logs = FALSE,
    projectName = "dnaEPICO", detPPath = NULL, detPThreshold = 0.01,
    cpgDetectionPath = NULL, sampleDetectionPath = NULL,
        logoPath = system.file("extdata",
        "dnaEPICO.svg",
        package = "dnaEPICO"
    ), imagePattern = "\\.(png|jpg|jpeg|gif|webp|svg|tif|tiff)$",
    recursive = TRUE
) {
    old_options <- options(stringsAsFactors = FALSE)
    on.exit(options(old_options), add = TRUE)
    detPThreshold <- validateProbabilityDnaEpico(
        detPThreshold,
        "detPThreshold"
    )
    model_sections <- normalizeModelSectionsDnamReport(modelSections)
    include_glm <- "glm" %in% model_sections
    include_lme <- "lme" %in% model_sections

    is_absolute_path <- function(path) {
        grepl("^[A-Za-z]:[/\\\\]|^/", path)
    }

    infer_workflow_root <- function(output_dir) {
        normalized_output <- if (is_absolute_path(output_dir)) {
            normalizePath(output_dir, winslash = "/", mustWork = FALSE)
        } else {
            normalizePath(file.path(getwd(), output_dir),
                winslash = "/",
                mustWork = FALSE
            )
        }

        if (identical(basename(dirname(normalized_output)), "reports")) {
            return(dirname(dirname(normalized_output)))
        }
        if (identical(basename(normalized_output), "reports")) {
            return(dirname(normalized_output))
        }

        dirname(normalized_output)
    }

    root_dir <- infer_workflow_root(outputDir)
    resolve_report_path <- function(path, base = root_dir) {
        if (is.null(path) || !nzchar(path)) {
            return("")
        }
        if (is_absolute_path(path)) {
            return(normalizePath(path, winslash = "/", mustWork = FALSE))
        }
        normalizePath(file.path(base, path),
            winslash = "/",
            mustWork = FALSE
        )
    }
    project_name <- projectName
    qc_dir <- resolve_report_path(enmixTab)
    preprocessing_dir <- resolve_report_path(qcTab)
    postprocessing_dir <- resolve_report_path(metricTab)
    sva_dir <- resolve_report_path(svaTab)

    infer_model_from_path <- function(path,
        step = "preprocessingMinfiEwasWater") {
        pieces <- strsplit(gsub("\\\\", "/", normalizePath(path,
            winslash = "/", mustWork = FALSE
        )), "/", fixed = TRUE)[[1]]
        step_idx <- which(pieces == step)
        figures_idx <- which(pieces == "figures")
        if (length(step_idx) && length(figures_idx)) {
            usable <- figures_idx[figures_idx < step_idx[[1]]]
            if (length(usable) && step_idx[[1]] - usable[[length(usable)]] >=
                2L) {
                return(pieces[[step_idx[[1]] - 1L]])
            }
        }
        basename(normalizePath(outputDir, winslash = "/", mustWork = FALSE))
    }

    model_name <- infer_model_from_path(qc_dir)
    logs_dir <- resolve_report_path(logTab)

    pheno_file <- if (is.null(phenoTab) || !nzchar(phenoTab)) {
        resolve_report_path(file.path(
            "data", model_name, "preprocessingMinfiEwasWater",
            "phenoLC.csv"
        ))
    } else {
        resolve_report_path(phenoTab)
    }
    detp_path <- if (is.null(detPPath) || !nzchar(detPPath)) {
        resolve_report_path(file.path(
            "rData", model_name, "preprocessingMinfiEwasWater",
            "qc", "detP_RGSet.RData"
        ))
    } else {
        resolve_report_path(detPPath)
    }
    cpg_detection_path <- if (is.null(cpgDetectionPath) ||
        !nzchar(cpgDetectionPath)) {
        resolve_report_path(file.path(
            "data", model_name, "preprocessingMinfiEwasWater",
            "cpgD.csv"
        ))
    } else {
        resolve_report_path(cpgDetectionPath)
    }
    sample_detection_path <- if (is.null(sampleDetectionPath) ||
        !nzchar(sampleDetectionPath)) {
        resolve_report_path(file.path(
            "data", model_name, "preprocessingMinfiEwasWater",
            "sampleD.csv"
        ))
    } else {
        resolve_report_path(sampleDetectionPath)
    }
    glm_table_path <- if (is.null(glmTab) || !nzchar(glmTab)) {
        resolve_report_path(file.path(
            "data", model_name, "methylationGLM",
            "annotatedGLM.xlsx"
        ))
    } else {
        resolve_report_path(glmTab)
    }
    lme_table_path <- if (is.null(lmeTab) || !nzchar(lmeTab)) {
        resolve_report_path(file.path(
            "data", model_name, "methylationLME",
            "annotatedLME.xlsx"
        ))
    } else {
        resolve_report_path(lmeTab)
    }
    glm_dir <- resolve_report_path(file.path(
        "figures", model_name,
        "methylationGLM"
    ))
    lme_dir <- resolve_report_path(file.path(
        "figures", model_name,
        "methylationLME"
    ))
    logo_candidates <- c(logoPath, system.file("extdata", "dnaEPICO.svg",
        package = "dnaEPICO"
    ), file.path(
        root_dir, "inst", "extdata",
        "dnaEPICO.svg"
    ), file.path(root_dir, "dnaEPICO.svg"))
    logo_candidates <- logo_candidates[nzchar(logo_candidates)]
    logo_source_path <- logo_candidates[file.exists(logo_candidates)][1]
    if (is.na(logo_source_path)) {
        logo_source_path <- ""
    }
    viewer_js_candidates <- c(system.file("extdata", "cpg-viewer.js",
        package = "dnaEPICO"
    ), file.path(
        getwd(), "inst", "extdata",
        "cpg-viewer.js"
    ), file.path(
        root_dir, "inst", "extdata",
        "cpg-viewer.js"
    ))
    viewer_js_candidates <- viewer_js_candidates[nzchar(viewer_js_candidates)]
    viewer_js_source_path <-
        viewer_js_candidates[file.exists(viewer_js_candidates)][1]
    if (is.na(viewer_js_source_path)) {
        viewer_js_source_path <- ""
    }

    project_dir <- resolve_report_path(outputDir)
    assets_dir <- file.path(project_dir, "assets")
    assets_figures_dir <- file.path(assets_dir, "figures")
    assets_logs_dir <- file.path(assets_dir, "logs")
    fig_dir <- assets_figures_dir

    image_pattern <- imagePattern
    magick_available <- requireNamespace("magick", quietly = TRUE)
    log_path <- resolveLogPathMinfiEwasWater(
        logs = logs, log_dir = logs_dir,
        log_file = "log_dnamReport.txt"
    )

    emitLogMinfiEwasWater(
        c(
            "============================================================",
            paste("Starting DNA Methylation Report Step:", project_name),
            paste("Project root:", root_dir), paste(
                "Output project:",
                project_dir
            ), "============================================================"
        ),
        verbose = FALSE, log_path = log_path
    )

    slash <- function(path) {
        gsub("\\\\", "/", path)
    }

    write_utf8 <- function(path, lines) {
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        con <- file(path, open = "w", encoding = "UTF-8")
        on.exit(close(con), add = TRUE)
        writeLines(lines, con = con, useBytes = TRUE)
    }

    find_quarto <- function() {
        env_candidates <- c(Sys.getenv("QUARTO_BIN"), Sys.getenv("QUARTO"))
        env_candidates <- env_candidates[nzchar(env_candidates)]
        env_candidates <- normalizePath(env_candidates,
            winslash = "/",
            mustWork = FALSE
        )
        env_candidates <- env_candidates[file.exists(env_candidates)]
        if (length(env_candidates)) {
            return(env_candidates[[1]])
        }

        quarto_bin <- Sys.which("quarto")
        if (nzchar(quarto_bin)) {
            return(unname(quarto_bin))
        }

        candidates <- c(
            "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe",
            "C:/Program Files/Quarto/bin/quarto.exe", file.path(
                Sys.getenv("LOCALAPPDATA"),
                "Programs", "Quarto", "bin", "quarto.exe"
            ), file.path(
                Sys.getenv("LOCALAPPDATA"),
                "quarto", "bin", "quarto.exe"
            )
        )
        candidates <- normalizePath(candidates,
            winslash = "/",
            mustWork = FALSE
        )
        candidates <- candidates[file.exists(candidates)]
        if (length(candidates)) {
            return(candidates[[1]])
        }

        ""
    }

    html_escape <- function(text) {
        text <- gsub("&", "&amp;", text, fixed = TRUE)
        text <- gsub("<", "&lt;", text, fixed = TRUE)
        text <- gsub(">", "&gt;", text, fixed = TRUE)
        text <- gsub("\"", "&quot;", text, fixed = TRUE)
        text
    }

    r_string <- function(text) {
        text <- gsub("\\\\", "\\\\\\\\", text)
        text <- gsub("\"", "\\\"", text, fixed = TRUE)
        paste0("\"", text, "\"")
    }

    slugify <- function(text) {
        text <- tools::file_path_sans_ext(basename(text))
        text <- tolower(text)
        text <- gsub("[^a-z0-9]+", "-", text)
        text <- gsub("(^-+|-+$)", "", text)
        if (!nzchar(text)) {
            text <- "figure"
        }
        text
    }

    pretty_label <- function(path) {
        label <- tools::file_path_sans_ext(basename(path))
        label <- gsub("([a-z0-9])([A-Z])", "\\1 \\2", label,
            perl = TRUE
        )
        label <- gsub("[_\\-]+", " ", label)
        label <- gsub("\\s+", " ", label)
        label <- trimws(label)
        words <- strsplit(tolower(label), "\\s+", perl = TRUE)[[1]]
        label <- paste0(toupper(substr(words, 1L, 1L)), substr(
            words,
            2L, nchar(words)
        ))
        label <- paste(label, collapse = " ")
        label <- gsub("\\bQc\\b", "QC", label)
        label <- gsub("\\bPca\\b", "PCA", label)
        label <- gsub("\\bSva\\b", "SVA", label)
        label <- gsub("\\bEnmix\\b", "ENmix", label)
        label
    }

    callout_lines <- function(text, type = "note") {
        c(sprintf("::: {.callout-%s}", type), text, ":::", "")
    }

    copy_log_asset <- function(src_path, output_name) {
        dir.create(assets_logs_dir, recursive = TRUE, showWarnings = FALSE)
        dest_path <- file.path(assets_logs_dir, output_name)

        if (file.exists(src_path)) {
            ok <- file.copy(src_path, dest_path, overwrite = TRUE)
            if (!ok) {
                stop("Failed to copy log file: ", src_path)
            }
        }

        list(
            exists = file.exists(src_path), source_path = slash(src_path),
            asset_path = slash(file.path("assets", "logs", output_name))
        )
    }

    copy_first_existing_log_asset <- function(src_paths, output_name) {
        existing <- src_paths[file.exists(src_paths)]
        src_path <- if (length(existing) > 0L) {
            existing[[1L]]
        } else {
            src_paths[[1L]]
        }
        copy_log_asset(src_path, output_name)
    }

    copy_figure_assets <- function(src_dir, asset_subdir) {
        dir.create(file.path(assets_figures_dir, asset_subdir),
            recursive = TRUE, showWarnings = FALSE
        )

        if (!dir.exists(src_dir)) {
            return(data.frame(
                title = character(), original_name = character(),
                asset_path = character(), browser_ready = logical(),
                converted = logical(), stringsAsFactors = FALSE
            ))
        }

        files <- list.files(src_dir,
            pattern = image_pattern,
            full.names = TRUE, recursive = recursive, ignore.case = TRUE
        )
        files <- sort(files)

        if (!length(files)) {
            return(data.frame(
                title = character(), original_name = character(),
                asset_path = character(), browser_ready = logical(),
                converted = logical(), stringsAsFactors = FALSE
            ))
        }

        rows <- vector("list", length(files))

        for (idx in seq_along(files)) {
            src_path <- files[[idx]]
            src_ext <- tolower(tools::file_ext(src_path))
            base_slug <- sprintf("%02d-%s", idx, slugify(src_path))
            converted <- FALSE
            browser_ready <- TRUE

            if (src_ext %in% c("tif", "tiff") && magick_available) {
                dest_name <- paste0(base_slug, ".png")
                dest_path <- file.path(
                    assets_figures_dir, asset_subdir,
                    dest_name
                )
                converted <- tryCatch(
                    {
                        img <- magick::image_read(src_path)
                        magick::image_write(img,
                            path = dest_path,
                            format = "png"
                        )
                        TRUE
                    },
                    error = function(e) FALSE
                )

                if (!isTRUE(converted)) {
                    dest_name <- paste0(base_slug, ".", src_ext)
                    dest_path <- file.path(
                        assets_figures_dir,
                        asset_subdir, dest_name
                    )
                    ok <- file.copy(src_path, dest_path, overwrite = TRUE)
                    if (!ok) {
                        stop("Failed to copy figure: ", src_path)
                    }
                    browser_ready <- FALSE
                }
            } else {
                dest_name <- paste0(base_slug, ".", src_ext)
                dest_path <- file.path(
                    assets_figures_dir, asset_subdir,
                    dest_name
                )
                ok <- file.copy(src_path, dest_path, overwrite = TRUE)
                if (!ok) {
                    stop("Failed to copy figure: ", src_path)
                }
                browser_ready <- !(src_ext %in% c("tif", "tiff"))
            }

            rows[[idx]] <- data.frame(
                title = pretty_label(src_path),
                original_name = basename(src_path),
                    asset_path = slash(file.path(
                    "assets",
                    "figures", asset_subdir, dest_name
                )), browser_ready = browser_ready,
                converted = converted, stringsAsFactors = FALSE
            )
        }

        do.call(rbind, rows)
    }

    build_figure_cards <- function(items, empty_message,
        title_prefix = "Figure",
                                    show_filename = FALSE,
                                        figure_titles = NULL) {
        if (!nrow(items)) {
            return(callout_lines(empty_message, type = "warning"))
        }

        lines <- character()

        if (any(!items$browser_ready)) {
            lines <- c(lines, callout_lines(
                paste(
                    "Some TIFF figures were copied without PNG conversion.",
                    "Install the `magick` package and rerun `dnamReport()` if you want those figures to render directly in the browser."
                ),
                type = "warning"
            ))
        }

        for (idx in seq_len(nrow(items))) {
            row <- items[idx, , drop = FALSE]
            caption_title <- if (!is.null(figure_titles) &&
                length(figure_titles) >=
                idx && nzchar(figure_titles[[idx]])) {
                figure_titles[[idx]]
            } else if ("title" %in% names(row) && nzchar(row$title[[1]])) {
                row$title[[1]]
            } else {
                sprintf("%s %s", title_prefix, idx)
            }
            lines <- c(lines, sprintf(
                "::: {.card title=\"%s\"}",
                html_escape(caption_title)
            ), "")

            if (isTRUE(row$browser_ready)) {
                lines <- c(
                    lines, "<figure class=\"qpasst-figure-card\">",
                    sprintf(
                        "  <img src=\"%s\" alt=\"%s\" loading=\"lazy\" />",
                        row$asset_path, html_escape(row$title)
                    )
                )
                if (isTRUE(show_filename)) {
                    lines <- c(lines, sprintf(
                        "  <figcaption><code>%s</code></figcaption>",
                        html_escape(row$original_name)
                    ))
                }
                lines <- c(lines, "</figure>", "")
            } else {
                download_label <- if (isTRUE(show_filename)) {
                    row$original_name
                } else {
                    caption_title
                }
                lines <- c(
                    lines, sprintf(
                        "<p class=\"qpasst-download\">Browser preview is unavailable for this TIFF file. <a href=\"%s\">Download %s</a>.</p>",
                        row$asset_path, html_escape(download_label)
                    ),
                    ""
                )
            }
            lines <- c(lines, ":::", "")
        }

        lines
    }

    js_quote_result_values <- function(values) {
        missing <- is.na(values)
        values <- enc2utf8(as.character(values))
        values <- gsub("\\", "\\\\", values, fixed = TRUE)
        values <- gsub("\"", "\\\"", values, fixed = TRUE)
        values <- gsub("\r", "\\r", values, fixed = TRUE)
        values <- gsub("\n", "\\n", values, fixed = TRUE)
        values <- gsub("\t", "\\t", values, fixed = TRUE)
        values <- gsub("<", "\\u003c", values, fixed = TRUE)
        values <- gsub(">", "\\u003e", values, fixed = TRUE)
        values <- gsub("&", "\\u0026", values, fixed = TRUE)
        values <- paste0("\"", values, "\"")
        values[missing] <- "null"
        values
    }

    js_result_array <- function(values) {
        paste0(
            "[", paste(js_quote_result_values(values), collapse = ","),
            "]"
        )
    }

    read_optional_workbook_sheet <- function(path, sheets, sheet) {
        if (!(sheet %in% sheets)) {
            return(data.frame())
        }
        tryCatch(openxlsx::read.xlsx(path, sheet = sheet, check.names = FALSE),
            error = function(e) data.frame()
        )
    }

    metadata_frame_to_list <- function(metadata) {
        if (!is.data.frame(metadata) || !all(c("Key", "Value") %in%
            names(metadata))) {
            return(list())
        }
        values <- as.character(metadata$Value)
        names(values) <- as.character(metadata$Key)
        as.list(values[nzchar(names(values))])
    }

    split_report_metadata_values <- function(value) {
        if (is.null(value) || !length(value) || is.na(value[[1L]]) ||
            !nzchar(value[[1L]])) {
            return(character())
        }
        values <- trimws(unlist(strsplit(as.character(value[[1L]]),
            ",",
            fixed = TRUE
        ), use.names = FALSE))
        unique(values[nzchar(values)])
    }

    extract_report_formula_phenotype <- function(formula) {
        formula <- as.character(formula[[1L]])
        rhs <- trimws(sub("^[^~]*~", "", formula))
        quoted_term <- regmatches(rhs, regexpr("`[^`]+`", rhs,
            perl = TRUE
        ))
        if (length(quoted_term) && nzchar(quoted_term) && !identical(
            quoted_term,
            character(0)
        )) {
            return(substring(quoted_term, 2L, nchar(quoted_term) -
                1L))
        }
        trimws(sub("\\s*(?:\\+|\\*|:|\\|).*", "", rhs, perl = TRUE))
    }

    resolve_report_formula_records <- function(dictionary, metadata = list()) {
        empty_records <- data.frame(
            phenotype = character(),
            result_column = character(), formula = character(),
            stringsAsFactors = FALSE, check.names = FALSE
        )
        if (!is.data.frame(dictionary) || !("Formula" %in% names(dictionary))) {
            return(empty_records)
        }

        keep <- !is.na(dictionary$Formula) & nzchar(trimws(as.character(dictionary$Formula)))
        if (!any(keep)) {
            return(empty_records)
        }
        formula_rows <- dictionary[keep, , drop = FALSE]
        formulas <- as.character(formula_rows$Formula)
        result_columns <- if ("Column" %in% names(formula_rows)) {
            as.character(formula_rows$Column)
        } else {
            rep("", length(formulas))
        }
        explicit_phenotypes <- if ("Phenotype" %in% names(formula_rows)) {
            as.character(formula_rows$Phenotype)
        } else {
            rep("", length(formulas))
        }
        known_phenotypes <- split_report_metadata_values(metadata$phenotypes)

        phenotype_labels <- vapply(seq_along(formulas), function(index) {
            if (!is.na(explicit_phenotypes[[index]]) &&
                nzchar(trimws(explicit_phenotypes[[index]]))) {
                return(trimws(explicit_phenotypes[[index]]))
            }
            if (length(known_phenotypes)) {
                matches <- known_phenotypes[vapply(
                    known_phenotypes,
                    function(phenotype) {
                        grepl(paste0("`", phenotype, "`"), formulas[[index]],
                            fixed = TRUE
                        ) || startsWith(
                            result_columns[[index]],
                            phenotype
                        )
                    }, logical(1)
                )]
                if (length(matches)) {
                    return(matches[[1L]])
                }
            }
            extracted <- extract_report_formula_phenotype(formulas[[index]])
            if (nzchar(extracted)) {
                extracted
            } else {
                result_columns[[index]]
            }
        }, character(1))

        records <- data.frame(
            phenotype = phenotype_labels, result_column = result_columns,
            formula = formulas, stringsAsFactors = FALSE, check.names = FALSE
        )
        duplicate_key <- paste(records$phenotype, records$formula,
            sep = "\r"
        )
        records <- records[!duplicated(duplicate_key), , drop = FALSE]
        rownames(records) <- NULL
        records
    }

    last_log_field <- function(path, label) {
        if (!file.exists(path)) {
            return("")
        }
        lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
        matches <- grep(paste0("^", label, ":"), trimws(lines),
            value = TRUE
        )
        if (!length(matches)) {
            return("")
        }
        trimws(sub("^[^:]+:", "", matches[[length(matches)]]))
    }

    resolve_lme_report_metadata <- function(metadata, dictionary,
                                            log_path) {
        values <- metadata_frame_to_list(metadata)
        fallback <- list(libraries = last_log_field(
            log_path,
            "LME libraries"
        ), correlation_structure = last_log_field(
            log_path,
            "Correlation structure"
        ), correlation_variable = last_log_field(
            log_path,
            "Correlation variable"
        ), interaction_term = last_log_field(
            log_path,
            "Interaction term"
        ), annotation_columns = last_log_field(
            log_path,
            "Annotation columns used"
        ), missing_annotation_columns = last_log_field(
            log_path,
            "Missing annotation columns"
        ))
        for (key in names(fallback)) {
            if (is.null(values[[key]]) || !nzchar(values[[key]])) {
                values[[key]] <- fallback[[key]]
            }
        }

        if (is.null(values$backend) || !nzchar(values$backend)) {
            values$backend <- if (grepl("nlme", values$libraries,
                ignore.case = TRUE
            )) {
                "nlme"
            } else if (grepl("lme4|lmerTest", values$libraries,
                ignore.case = TRUE
            )) {
                "lme4"
            } else {
                ""
            }
        }
        if (is.null(values$fitting_function) ||
            !nzchar(values$fitting_function)) {
            values$fitting_function <- if (identical(
                values$backend,
                "nlme"
            )) {
                "nlme::lme"
            } else if (identical(values$backend, "lme4")) {
                "lmerTest::lmer"
            } else {
                ""
            }
        }
        if (is.null(values$correlation_structure) ||
            !nzchar(values$correlation_structure)) {
            values$correlation_structure <- "none"
        }
        if (is.null(values$correlation_variable) ||
            !nzchar(values$correlation_variable)) {
            values$correlation_variable <- "None"
        }
        if (is.null(values$interaction_term) ||
            !nzchar(values$interaction_term)) {
            values$interaction_term <- "None"
        }
        values
    }

    same_report_path <- function(source, target) {
        source <- normalizePath(source, winslash = "/", mustWork = FALSE)
        target <- normalizePath(target, winslash = "/", mustWork = FALSE)
        if (identical(.Platform$OS.type, "windows")) {
            identical(tolower(source), tolower(target))
        } else {
            identical(source, target)
        }
    }

    prepare_table_viewer_directory <- function(var_prefix) {
        result_dir <- file.path(assets_dir, "results", var_prefix)
        dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
        stale_chunks <- list.files(result_dir,
            pattern = "^chunk-[0-9]+\\.js$",
            full.names = TRUE
        )
        if (length(stale_chunks)) {
            unlink(stale_chunks, force = TRUE)
        }
        result_dir
    }

    write_table_viewer_chunk <- function(chunk, chunk_number,
                                            result_dir, var_prefix, id_column) {
        encoded_columns <- lapply(chunk, js_quote_result_values)
        encoded_rows <- if (nrow(chunk)) {
            do.call(paste, c(encoded_columns, sep = ","))
        } else {
            character()
        }
        chunk_key <- paste0(var_prefix, ":", chunk_number)
        write_utf8(file.path(result_dir, sprintf(
            "chunk-%04d.js",
            chunk_number
        )), paste0(
            "window.dnaEPICOResultChunks=window.dnaEPICOResultChunks||{};",
            "window.dnaEPICOResultChunks[", js_quote_result_values(chunk_key),
            "]={rows:[", paste0("[", encoded_rows, "]", collapse = ","),
            "]};"
        ))
        ids <- if (nzchar(id_column)) {
            as.character(chunk[[id_column]])
        } else {
            rep("", nrow(chunk))
        }
        ids <- ids[!is.na(ids)]
        data.frame(number = chunk_number, first_id = if (length(ids)) {
            ids[[1L]]
        } else {
            ""
        }, last_id = if (length(ids)) {
            ids[[length(ids)]]
        } else {
            ""
        }, stringsAsFactors = FALSE)
    }

    write_table_viewer_manifest <- function(columns, n_rows,
                                            chunk_size, chunk_index, var_prefix,
                                                id_column, downloads,
                                            item_singular, item_plural) {
        chunk_objects <- if (nrow(chunk_index)) {
            vapply(seq_len(nrow(chunk_index)), function(i) {
                paste0(
                    "{\"number\":", chunk_index$number[[i]],
                    ",\"firstId\":",
                        js_quote_result_values(chunk_index$first_id[[i]]),
                    ",\"lastId\":",
                        js_quote_result_values(chunk_index$last_id[[i]]),
                    "}"
                )
            }, character(1))
        } else {
            character()
        }
        download_objects <- if (length(downloads)) {
            vapply(downloads, function(download) {
                paste0(
                    "{\"href\":", js_quote_result_values(download$href),
                    ",\"label\":", js_quote_result_values(download$label),
                    "}"
                )
            }, character(1))
        } else {
            character()
        }
        manifest <- paste0(
            "window.dnaEPICOResultManifests=window.dnaEPICOResultManifests||{};",
            "window.dnaEPICOResultManifests[",
                js_quote_result_values(var_prefix),
            "]={", "\"key\":", js_quote_result_values(var_prefix),
            ",", "\"columns\":", js_result_array(columns), ",",
            "\"totalRows\":", n_rows, ",", "\"chunkSize\":",
            chunk_size, ",", "\"maxCachedChunks\":4,", "\"idColumn\":",
            js_quote_result_values(id_column), ",", "\"basePath\":",
            js_quote_result_values(paste0(
                "assets/results/",
                var_prefix
            )), ",", "\"itemSingular\":", js_quote_result_values(item_singular),
            ",", "\"itemPlural\":", js_quote_result_values(item_plural),
            ",", "\"downloads\":[", paste(download_objects, collapse = ","),
            "],", "\"chunks\":[", paste(chunk_objects, collapse = ","),
            "]};"
        )
        result_dir <- file.path(assets_dir, "results", var_prefix)
        write_utf8(file.path(result_dir, "manifest.js"), manifest)
    }

    table_viewer_asset_result <- function(columns, n_rows, var_prefix,
                                            id_column, item_singular,
                                                item_plural, search_label,
                                            search_placeholder) {
        list(
            key = var_prefix, columns = columns, n_rows = n_rows,
            n_cols = length(columns), id_column = id_column,
            item_singular = item_singular, item_plural = item_plural,
            search_label = search_label,
                search_placeholder = search_placeholder,
            manifest_path = paste0(
                "assets/results/", var_prefix,
                "/manifest.js"
            )
        )
    }

    write_table_viewer_assets <- function(table_data, var_prefix,
                                            id_column, downloads,
                                                chunk_size = 5000L,
                                                item_singular = "CpG",
                                            item_plural = "CpGs",
                                                search_label = "Find CpG",
                                                search_placeholder = "e.g. cg00000029") {
        table_data <- as.data.frame(table_data,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        chunk_size <- max(100L, as.integer(chunk_size))
        result_dir <- prepare_table_viewer_directory(var_prefix)

        n_rows <- nrow(table_data)
        n_chunks <- if (n_rows) {
            ceiling(n_rows / chunk_size)
        } else {
            0L
        }
        chunk_rows <- vector("list", n_chunks)
        if (n_chunks) {
            for (chunk_idx in seq_len(n_chunks)) {
                start <- (chunk_idx - 1L) * chunk_size + 1L
                end <- min(chunk_idx * chunk_size, n_rows)
                chunk <- table_data[start:end, , drop = FALSE]
                chunk_rows[[chunk_idx]] <- write_table_viewer_chunk(
                    chunk = chunk,
                    chunk_number = chunk_idx, result_dir = result_dir,
                    var_prefix = var_prefix, id_column = id_column
                )
            }
        }
        chunk_index <- if (length(chunk_rows)) {
            do.call(rbind, chunk_rows)
        } else {
            data.frame()
        }
        write_table_viewer_manifest(
            columns = names(table_data),
            n_rows = n_rows, chunk_size = chunk_size, chunk_index = chunk_index,
            var_prefix = var_prefix, id_column = id_column,
                downloads = downloads,
            item_singular = item_singular, item_plural = item_plural
        )

        table_viewer_asset_result(
            columns = names(table_data),
            n_rows = n_rows, var_prefix = var_prefix, id_column = id_column,
            item_singular = item_singular, item_plural = item_plural,
            search_label = search_label, search_placeholder = search_placeholder
        )
    }

    write_delimited_table_viewer_assets <- function(table_path,
                                                    var_prefix, id_column,
                                                        downloads,
                                                        expected_rows,
                                                        expected_columns,
                                                    chunk_size = 5000L,
                                                        item_singular = "CpG",
                                                        item_plural = "CpGs",
                                                    search_label = "Find CpG",
                                                        search_placeholder = "e.g. cg00000029") {
        chunk_size <- max(100L, as.integer(chunk_size))
        result_dir <- prepare_table_viewer_directory(var_prefix)
        chunk_state <- new.env(parent = emptyenv())
        chunk_state$rows <- list()
        stream_result <- streamReportTableDnaEpico(
            tableFile = table_path,
            chunkSize = chunk_size, expectedRows = expected_rows,
            expectedColumns = expected_columns, chunkHandler = function(chunk,
                                                                        chunk_number) {
                chunk_state$rows[[chunk_number]] <- write_table_viewer_chunk(
                    chunk = chunk,
                    chunk_number = chunk_number, result_dir = result_dir,
                    var_prefix = var_prefix, id_column = id_column
                )
            }
        )
        chunk_rows <- chunk_state$rows
        columns <- stream_result$columns
        n_rows <- stream_result$nRows
        if (!(id_column %in% columns)) {
            stop("The report table identifier column was not found.",
                call. = FALSE
            )
        }
        chunk_index <- if (length(chunk_rows)) {
            do.call(rbind, chunk_rows)
        } else {
            data.frame()
        }
        write_table_viewer_manifest(
            columns = columns, n_rows = n_rows,
            chunk_size = chunk_size, chunk_index = chunk_index,
            var_prefix = var_prefix, id_column = id_column,
                downloads = downloads,
            item_singular = item_singular, item_plural = item_plural
        )

        table_viewer_asset_result(
            columns = columns, n_rows = n_rows,
            var_prefix = var_prefix, id_column = id_column,
                item_singular = item_singular,
            item_plural = item_plural, search_label = search_label,
            search_placeholder = search_placeholder
        )
    }

    prepare_xlsx_table_assets <- function(data_path, sheet, var_prefix,
                                            analysis = c("glm", "lme"),
                                                chunk_size = 5000L,
                                                model_log_path = "") {
        analysis <- match.arg(analysis)
        source_data_path <- if (grepl(
            "^[A-Za-z]:[/\\\\]|^/",
            data_path
        )) {
            data_path
        } else {
            file.path(root_dir, data_path)
        }
        if (!file.exists(source_data_path)) {
            return(list(
                ok = FALSE, error = paste0(
                    "Data file not found at `",
                    slash(source_data_path), "`."
                ), source_path = source_data_path,
                analysis = analysis, metadata = list()
            ))
        }

        result_dir <- file.path(assets_dir, "results", var_prefix)
        dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
        xlsx_name <- basename(source_data_path)
        legacy_xlsx_target <- file.path(result_dir, xlsx_name)
        if (!same_report_path(source_data_path, legacy_xlsx_target) &&
            file.exists(legacy_xlsx_target)) {
            unlink(legacy_xlsx_target, force = TRUE)
        }

        sidecar <- resolveReportTableSidecarDnaEpico(
            workbookFile = source_data_path,
            sidecarDir = result_dir,
            sheet = sheet
        )
        if (isTRUE(sidecar$ok)) {
            dictionary <- tryCatch(utils::read.delim(sidecar$dictionary,
                check.names = FALSE, stringsAsFactors = FALSE,
                quote = "\"", comment.char = ""
            ), error = function(error) error)
            workbook_metadata <- if (file.exists(sidecar$workbookMetadata)) {
                tryCatch(utils::read.delim(sidecar$workbookMetadata,
                    check.names = FALSE, stringsAsFactors = FALSE,
                    quote = "\"", comment.char = ""
                ), error = function(error) error)
            } else {
                data.frame()
            }
            if (inherits(dictionary, "error") || inherits(
                workbook_metadata,
                "error"
            )) {
                sidecar$ok <- FALSE
            }
        }
        if (!isTRUE(sidecar$ok)) {
            sheets <- tryCatch(openxlsx::getSheetNames(source_data_path),
                error = function(e) character()
            )
            dictionary <- read_optional_workbook_sheet(
                source_data_path,
                sheets, "dictionary"
            )
            workbook_metadata <- read_optional_workbook_sheet(
                source_data_path,
                sheets, "metadata"
            )
        }
        metadata <- if (identical(analysis, "lme")) {
            resolve_lme_report_metadata(
                workbook_metadata, dictionary,
                model_log_path
            )
        } else {
            metadata_frame_to_list(workbook_metadata)
        }
        metadata$formula_records <- resolve_report_formula_records(
            dictionary,
            metadata
        )
        metadata$formulas <- unique(metadata$formula_records$formula)

        text_name <- paste0(sheet, ".tsv.gz")
        text_target <- file.path(result_dir, text_name)
        downloads <- list(list(href = paste0(
                "assets/results/", var_prefix,
                "/", text_name
            ), label = "Download all results (TSV.gz)")
        )

        viewer_assets <- NULL
        source_mode <- "xlsx_fallback"
        if (isTRUE(sidecar$ok)) {
            sidecar_staged <- same_report_path(
                sidecar$table,
                text_target
            ) || isTRUE(file.copy(sidecar$table,
                text_target,
                overwrite = TRUE
            ))
            if (sidecar_staged && file.exists(text_target)) {
                viewer_assets <- tryCatch(
                    write_delimited_table_viewer_assets(
                        table_path = text_target,
                        var_prefix = var_prefix, id_column = sidecar$idColumn,
                        downloads = downloads, expected_rows = sidecar$rows,
                        expected_columns = sidecar$columns,
                            chunk_size = chunk_size
                    ),
                    error = function(error) NULL
                )
                if (!is.null(viewer_assets)) {
                    source_mode <- "streamed_sidecar"
                }
            }
        }

        if (is.null(viewer_assets)) {
            table_data <- tryCatch(openxlsx::read.xlsx(source_data_path,
                sheet = sheet, check.names = FALSE
            ), error = function(e) e)
            if (inherits(table_data, "error")) {
                return(list(
                    ok = FALSE, error = conditionMessage(table_data),
                    source_path = source_data_path, analysis = analysis,
                    metadata = list()
                ))
            }
            report_table <- sortReportTableDnaEpico(table_data)
            table_data <- report_table$data
            data.table::fwrite(table_data,
                file = text_target,
                sep = "\t", quote = TRUE, na = "", compress = "gzip",
                showProgress = FALSE
            )
            viewer_assets <- write_table_viewer_assets(
                table_data = table_data,
                var_prefix = var_prefix, id_column = report_table$idColumn,
                downloads = downloads, chunk_size = chunk_size
            )
        }

        c(list(
            ok = TRUE, error = NULL, source_path = source_data_path,
            analysis = analysis, dictionary = dictionary, metadata = metadata,
            source_mode = source_mode, text_path = paste0(
                "assets/results/",
                var_prefix, "/", text_name
            )
        ), viewer_assets)
    }

    prepare_csv_table_assets <- function(data_path, var_prefix,
                                            front_columns = character(),
                                                chunk_size = 5000L) {
        source_data_path <- if (grepl(
            "^[A-Za-z]:[/\\\\]|^/",
            data_path
        )) {
            data_path
        } else {
            file.path(root_dir, data_path)
        }
        if (!file.exists(source_data_path)) {
            return(list(
                ok = FALSE, error = paste0(
                    "Data file not found at `",
                    slash(source_data_path), "`."
                ), source_path = source_data_path,
                metadata = list()
            ))
        }

        table_data <- tryCatch(utils::read.csv(source_data_path,
            check.names = FALSE
        ), error = function(e) e)
        if (inherits(table_data, "error")) {
            return(list(
                ok = FALSE, error = conditionMessage(table_data),
                source_path = source_data_path, metadata = list()
            ))
        }
        table_data <- as.data.frame(table_data,
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
        front_columns <- front_columns[nzchar(front_columns) &
            front_columns %in% names(table_data)]
        table_data <- table_data[, c(front_columns, setdiff(
            names(table_data),
            front_columns
        )), drop = FALSE]
        id_column <- if (length(front_columns)) {
            front_columns[[1L]]
        } else if (ncol(table_data)) {
            names(table_data)[[1L]]
        } else {
            ""
        }
        sort_columns <- unique(c(id_column, front_columns[-1L]))
        sort_columns <- sort_columns[nzchar(sort_columns) & sort_columns %in%
            names(table_data)]
        if (length(sort_columns) && nrow(table_data) > 1L) {
            sort_values <- lapply(table_data[sort_columns], as.character)
            table_data <- table_data[do.call(order, c(
                sort_values,
                list(na.last = TRUE)
            )), , drop = FALSE]
            rownames(table_data) <- NULL
        }

        result_dir <- file.path(assets_dir, "results", var_prefix)
        dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
        csv_name <- basename(source_data_path)
        csv_target <- file.path(result_dir, csv_name)
        staged <- same_report_path(source_data_path, csv_target) ||
            isTRUE(file.copy(source_data_path, csv_target, overwrite = TRUE))
        if (!staged || !file.exists(csv_target)) {
            return(list(
                ok = FALSE, error = paste0(
                    "Could not copy the data file into the report assets from `",
                    slash(source_data_path), "`."
                ), source_path = source_data_path,
                metadata = list()
            ))
        }

        search_label <- if (nzchar(id_column)) {
            paste("Find", id_column)
        } else {
            "Find row"
        }
        first_id <- if (nzchar(id_column) && nrow(table_data)) {
            as.character(table_data[[id_column]][[1L]])
        } else {
            ""
        }
        viewer_assets <- write_table_viewer_assets(
            table_data = table_data,
            var_prefix = var_prefix, id_column = id_column,
                downloads = list(list(href = paste0(
                "assets/results/",
                var_prefix, "/", csv_name
            ), label = "Download data (CSV)")),
            chunk_size = chunk_size, item_singular = "record",
            item_plural = "records", search_label = search_label,
            search_placeholder = if (nzchar(first_id)) {
                paste("e.g.", first_id)
            } else {
                "Enter an identifier"
            }
        )

        c(list(
            ok = TRUE, error = NULL, source_path = source_data_path,
            metadata = list(), csv_path = paste0(
                "assets/results/",
                var_prefix, "/", csv_name
            )
        ), viewer_assets)
    }

    build_result_table_section <- function(title, table_assets) {
        if (!isTRUE(table_assets$ok)) {
            return(c(
                sprintf("::: {.card title=\"%s\"}", html_escape(title)),
                "", callout_lines(table_assets$error, type = "warning"),
                ":::", ""
            ))
        }

        search_label <- if (!is.null(table_assets$search_label) &&
            nzchar(table_assets$search_label)) {
            table_assets$search_label
        } else {
            "Find identifier"
        }
        search_placeholder <- if (!is.null(table_assets$search_placeholder) &&
            nzchar(table_assets$search_placeholder)) {
            table_assets$search_placeholder
        } else {
            "Enter an identifier"
        }
        pagination_lines <- function(position) {
            c(
                sprintf(
                    "  <div class=\"dnaepico-viewer-pagination dnaepico-viewer-pagination-%s\">",
                    position
                ), "    <button type=\"button\" data-role=\"first\" title=\"First page\" aria-label=\"First page\">&laquo;</button>",
                "    <button type=\"button\" data-role=\"previous\" title=\"Previous page\" aria-label=\"Previous page\">&lsaquo;</button>",
                "    <label>Page <input data-role=\"page-number\" type=\"number\" min=\"1\" value=\"1\" /></label>",
                "    <span data-role=\"page-count\"></span>",
                "    <button type=\"button\" data-role=\"next\" title=\"Next page\" aria-label=\"Next page\">&rsaquo;</button>",
                "    <button type=\"button\" data-role=\"last\" title=\"Last page\" aria-label=\"Last page\">&raquo;</button>",
                "  </div>"
            )
        }

        c(
            sprintf("::: {.card title=\"%s\"}", html_escape(title)),
            "", sprintf(
                "<div class=\"dnaepico-cpg-viewer\" data-result-key=\"%s\">",
                html_escape(table_assets$key)
            ), "", "```{=html}",
            "  <div class=\"dnaepico-viewer-toolbar\">",
                "    <label>Rows per page <select data-role=\"page-size\"><option>25</option><option>50</option><option>100</option></select></label>",
            sprintf(
                "    <label>%s <input data-role=\"cpg-search\" type=\"search\" placeholder=\"%s\" /></label>",
                html_escape(search_label), html_escape(search_placeholder)
            ),
            "    <button type=\"button\" data-role=\"find-cpg\">Find</button>",
            "  </div>", "  <div class=\"dnaepico-viewer-filter\">",
            "    <label>Filter column <select data-role=\"filter-column\"><option value=\"\">Choose a column</option></select></label>",
            "    <label>Condition <select data-role=\"filter-operator\"><option value=\"contains\">contains</option><option value=\"equals\">equals</option><option value=\"lt\">&lt;</option><option value=\"lte\">&le;</option><option value=\"gt\">&gt;</option><option value=\"gte\">&ge;</option></select></label>",
            "    <label>Value <input data-role=\"filter-value\" type=\"text\" placeholder=\"e.g. 0.05\" /></label>",
            "    <button type=\"button\" data-role=\"apply-filter\">Apply filter</button>",
            "    <button type=\"button\" data-role=\"clear-filter\" disabled>Clear</button>",
            "    <span data-role=\"filter-summary\" aria-live=\"polite\"></span>",
            "  </div>",
                "  <div class=\"dnaepico-viewer-downloads\" data-role=\"downloads\"></div>",
            "  <div class=\"dnaepico-viewer-status\" data-role=\"status\" aria-live=\"polite\">Loading results&hellip;</div>",
            pagination_lines("top"),
                "  <div class=\"dnaepico-viewer-table-wrap\">",
            "    <table class=\"table table-striped table-sm dnaepico-viewer-table\">",
            "      <thead data-role=\"head\"></thead>",
                "      <tbody data-role=\"body\"></tbody>",
            "    </table>", "  </div>", pagination_lines("bottom"),
            "  <noscript>JavaScript is required for paged browsing. Use a complete-file download instead.</noscript>",
            "```", "", "</div>", sprintf(
                "<script src=\"%s\"></script>",
                table_assets$manifest_path
            ), "<script src=\"assets/cpg-viewer.js\"></script>",
            ":::", ""
        )
    }

    build_data_frame_table_section <- function(title, data, empty_message,
                                                preview_rows = 25L) {
        if (!is.data.frame(data) || !nrow(data)) {
            return(c(
                sprintf("::: {.card title=\"%s\"}", html_escape(title)),
                "", callout_lines(empty_message, type = "warning"),
                ":::", ""
            ))
        }

        if (is.finite(preview_rows) && nrow(data) > as.integer(preview_rows)) {
            data <- utils::head(data, as.integer(preview_rows))
        }

        header_cells <- paste0(
            "    <th>", html_escape(names(data)),
            "</th>"
        )
        body_rows <- unlist(lapply(seq_len(nrow(data)), function(i) {
            row_values <- vapply(data[i, , drop = FALSE], function(value) {
                if (is.na(value)) {
                    ""
                } else {
                    as.character(value)
                }
            }, character(1))
            c("  <tr>", paste0(
                "    <td>", html_escape(row_values),
                "</td>"
            ), "  </tr>")
        }), use.names = FALSE)

        c(
            sprintf("::: {.card title=\"%s\"}", html_escape(title)),
            "", "<div class=\"table-responsive\">",
                "  <table class=\"table table-striped table-sm\">",
            "  <thead>", "  <tr>", header_cells, "  </tr>", "  </thead>",
            "  <tbody>", body_rows, "  </tbody>", "  </table>",
            "</div>", ":::", ""
        )
    }

    format_count <- function(value) {
        if (!length(value) || is.na(value)) {
            return("not available")
        }
        format(as.integer(round(value)),
            big.mark = ",", scientific = FALSE,
            trim = TRUE
        )
    }

    format_model_formula <- function(value) {
        value <- paste(as.character(value), collapse = "")
        value <- gsub("`", "", value, fixed = TRUE)
        value <- gsub("\\s+", " ", value, perl = TRUE)
        value <- gsub("\\s+\\)", ")", value, perl = TRUE)
        trimws(value)
    }

    build_model_formula_notes <- function(formula_records, model_label) {
        if (!is.data.frame(formula_records) || !nrow(formula_records)) {
            return(character())
        }

        display_formulas <- vapply(
            formula_records$formula, format_model_formula,
            character(1)
        )
        if (identical(tolower(model_label), "nlme")) {
            display_formulas <- sub("^LME:", "nlme:", display_formulas,
                ignore.case = TRUE
            )
        }
        if (nrow(formula_records) == 1L) {
            return(paste0(
                "The recorded model formula is <code class=\"dnaepico-model-formula\">",
                html_escape(display_formulas[[1L]]), "</code>."
            ))
        }

        formula_items <- vapply(
            seq_len(nrow(formula_records)),
            function(index) {
                phenotype <- formula_records$phenotype[[index]]
                if (is.na(phenotype) || !nzchar(trimws(phenotype))) {
                    phenotype <- paste("Model", index)
                }
                paste0(
                    "<div class=\"dnaepico-model-formula-item\">",
                    "<strong>", html_escape(phenotype), "</strong>",
                    "<code class=\"dnaepico-model-formula\">",
                    html_escape(display_formulas[[index]]), "</code>",
                    "</div>"
                )
            }, character(1)
        )

        c(sprintf(
            "%d phenotype-specific models were fitted.",
            nrow(formula_records)
        ), paste0(
            "<details class=\"dnaepico-model-formulas\">",
            "<summary>View recorded model formulas (", nrow(formula_records),
            ")</summary>", "<div class=\"dnaepico-model-formula-list\">",
            paste(formula_items, collapse = ""), "</div>", "</details>"
        ))
    }

    format_decimal <- function(value, digits = 2L) {
        if (!length(value) || is.na(value)) {
            return("not available")
        }
        format(round(as.numeric(value), digits),
            nsmall = digits,
            scientific = FALSE, trim = TRUE
        )
    }

    plural <- function(n, singular, plural_form = paste0(
                            singular,
                            "s"
                        )) {
        if (length(n) && !is.na(n) && as.integer(n) == 1L) {
            singular
        } else {
            plural_form
        }
    }

    is_numeric_like <- function(values) {
        values <- as.character(values)
        grepl("^\\s*[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?\\s*$",
            values,
            perl = TRUE
        )
    }

    sort_values <- function(values) {
        values <- unique(as.character(values[!is.na(values)]))
        values <- values[nzchar(trimws(values))]
        if (!length(values)) {
            return(character())
        }
        if (all(is_numeric_like(values))) {
            numeric_values <- as.numeric(values)
            values[order(numeric_values)]
        } else {
            sort(values)
        }
    }

    collapse_values <- function(values) {
        values <- sort_values(values)
        if (!length(values)) {
            return("not available")
        }
        if (length(values) == 1L) {
            return(values[[1]])
        }
        if (length(values) == 2L) {
            return(paste(values, collapse = " and "))
        }
        paste0(
            paste(values[-length(values)], collapse = ", "),
            ", and ", values[[length(values)]]
        )
    }

    is_blank <- function(values) {
        is.na(values) | !nzchar(trimws(as.character(values)))
    }

    safe_read_table_file <- function(path, sheet = 1) {
        if (!file.exists(path)) {
            return(NULL)
        }
        extension <- tolower(sub("^.*\\.([^.]+)$", "\\1", basename(path)))
        tryCatch(if (extension %in% c("xlsx", "xlsm")) {
            openxlsx::read.xlsx(path, sheet = sheet, check.names = FALSE)
        } else {
            utils::read.csv(path, check.names = FALSE)
        }, error = function(e) NULL)
    }

    read_detection_tables <- function(detp_path, threshold, cpg_path = "",
                                        sample_path = "") {
        result <- list(
            source = "detP", path = slash(detp_path),
            exists = FALSE, threshold = threshold, cpg = NULL,
            sample = NULL, error = NULL
        )

        if (file.exists(detp_path)) {
            loaded <- tryCatch(
                {
                    detp_env <- new.env(parent = emptyenv())
                    object_names <- load(detp_path, envir = detp_env)
                    detp_object <- if ("detP" %in% object_names) {
                        detp_env$detP
                    } else {
                        matrix_objects <- object_names[vapply(
                            object_names,
                            function(name) {
                                obj <- get(name, envir = detp_env)
                                is.matrix(obj) || is.data.frame(obj)
                            }, logical(1)
                        )]
                        if (length(matrix_objects)) {
                            get(matrix_objects[[1]], envir = detp_env)
                        } else {
                            NULL
                        }
                    }

                    if (is.null(detp_object)) {
                        stop("No `detP` matrix-like object was found.",
                            call. = FALSE
                        )
                    }

                    detp_matrix <- as.matrix(detp_object)
                    storage.mode(detp_matrix) <- "numeric"
                    if (is.null(colnames(detp_matrix))) {
                        colnames(detp_matrix) <- sprintf(
                            "Sample_%s",
                            seq_len(ncol(detp_matrix))
                        )
                    }

                    assessed <- is.finite(detp_matrix) & detp_matrix >=
                        0 & detp_matrix <= 1
                    detected <- assessed & detp_matrix < threshold
                    assessed_by_cpg <- rowSums(assessed)
                    detected_by_cpg <- rowSums(detected)
                    cpg <- data.frame(
                        metric = c(
                            "Total CpGs assessed",
                            "CpGs detected in at least one sample",
                                "CpGs detected in all samples",
                            "CpGs never detected",
                                "CpGs with no valid detection P values"
                        ),
                        nCpGs = c(nrow(detp_matrix), sum(detected_by_cpg >=
                            1L), sum(assessed_by_cpg == ncol(detp_matrix) &
                            detected_by_cpg == ncol(detp_matrix)),
                                sum(assessed_by_cpg >
                            0L & detected_by_cpg == 0L), sum(assessed_by_cpg ==
                            0L)), stringsAsFactors = FALSE
                    )
                    assessed_by_sample <- colSums(assessed)
                    detected_by_sample <- colSums(detected)
                    percent_detected <- rep(NA_real_,
                        length(assessed_by_sample))
                    has_assessed_values <- assessed_by_sample > 0L
                    percent_detected[has_assessed_values] <- 100 *
                        detected_by_sample[has_assessed_values] / assessed_by_sample[has_assessed_values]
                    sample <- data.frame(
                        UID = colnames(detp_matrix),
                        nAssessed = assessed_by_sample,
                            nDetected = detected_by_sample,
                        pDetected = round(percent_detected, 2),
                            stringsAsFactors = FALSE
                    )
                    list(cpg = cpg, sample = sample)
                },
                error = function(e) e
            )

            if (inherits(loaded, "error")) {
                result$error <- conditionMessage(loaded)
            } else {
                result$exists <- TRUE
                result$cpg <- loaded$cpg
                result$sample <- loaded$sample
                return(result)
            }
        }

        cpg <- safe_read_table_file(cpg_path)
        sample <- safe_read_table_file(sample_path)
        if (!is.null(cpg) || !is.null(sample)) {
            result$source <- "csv"
            result$exists <- TRUE
            result$cpg <- cpg
            result$sample <- sample
        }

        result
    }

    detection_table_warning <- function(detection_tables) {
        if (!is.null(detection_tables$error)) {
            return(paste0(
                "Unable to build detection tables from `",
                detection_tables$path, "`: ", detection_tables$error
            ))
        }
        paste0(
            "Detection P-value file not found at `", detection_tables$path,
            "`."
        )
    }

    pick_column <- function(data, candidates) {
        if (is.null(data) || !ncol(data)) {
            return(NULL)
        }
        if (any(candidates %in% names(data))) {
            return(candidates[candidates %in% names(data)][[1]])
        }
        lower_names <- tolower(names(data))
        for (candidate in candidates) {
            idx <- match(tolower(candidate), lower_names)
            if (!is.na(idx)) {
                return(names(data)[[idx]])
            }
        }
        NULL
    }

    pick_participant_column <- function(data) {
        exact_col <- pick_column(data, c(
            "Person_ID", "PersonID",
            "person_id", "Person", "person", "Participant_ID",
            "ParticipantID", "participant_id", "Participant",
            "Subject_ID", "SubjectID", "Individual_ID", "IndividualID",
            "UID", "ID"
        ))
        if (!is.null(exact_col)) {
            return(exact_col)
        }

        id_cols <- grep("id", names(data),
            ignore.case = TRUE,
            value = TRUE
        )
        technical_id <-
            grepl("sentrix|sample|array|slide|position|time|visit|plate|well|pool|probe",
            id_cols,
            ignore.case = TRUE
        )
        id_cols <- id_cols[!technical_id]
        if (length(id_cols)) {
            return(id_cols[[1]])
        }

        NULL
    }

    pick_timepoint_column <- function(data) {
        pick_column(data, c(
            "Timepoint", "Tiempoint", "Time_Point",
            "Timepoint_ID", "TimepointID", "Visit", "VisitID",
            "Visit_ID"
        ))
    }

    summarize_dataset <- function(data_path) {
        data <- safe_read_table_file(data_path)
        summary <- list(
            path = slash(data_path), exists = !is.null(data),
            n_rows = NA_integer_, n_cols = NA_integer_, participant_col = NULL,
            n_participants = NA_integer_, timepoint_col = NULL,
            timepoints = character(), n_timepoints = NA_integer_,
            complete_timepoint_participants = NA_integer_,
                repeated_participant_timepoint_rows = NA_integer_,
            survey_columns = character(),
                participants_missing_survey = NA_integer_
        )

        if (is.null(data)) {
            return(summary)
        }

        summary$n_rows <- nrow(data)
        summary$n_cols <- ncol(data)
        participant_col <- pick_participant_column(data)
        timepoint_col <- pick_timepoint_column(data)
        summary$participant_col <- participant_col
        summary$timepoint_col <- timepoint_col

        if (!is.null(participant_col)) {
            participants <- as.character(data[[participant_col]])
            participants <- participants[!is_blank(participants)]
            summary$n_participants <- length(unique(participants))
        }

        if (!is.null(timepoint_col)) {
            timepoints <- as.character(data[[timepoint_col]])
            summary$timepoints <- sort_values(timepoints)
            summary$n_timepoints <- length(summary$timepoints)
        }

        if (!is.null(participant_col) && !is.null(timepoint_col)) {
            participants <- as.character(data[[participant_col]])
            timepoints <- as.character(data[[timepoint_col]])
            complete_rows <- !is_blank(participants) & !is_blank(timepoints)
            keys <- paste(participants[complete_rows],
                timepoints[complete_rows],
                sep = "\r"
            )
            summary$repeated_participant_timepoint_rows <- length(keys) -
                length(unique(keys))

            split_timepoints <- split(
                timepoints[complete_rows],
                participants[complete_rows]
            )
            summary$complete_timepoint_participants <- sum(vapply(
                split_timepoints,
                function(values) {
                    all(summary$timepoints %in%
                        unique(as.character(values)))
                }, logical(1)
            ))
        }

        survey_columns <- grep("^(MHC_|BDSST|BRS_|SS_|WHO_)",
            names(data),
            value = TRUE
        )
        summary$survey_columns <- survey_columns
        if (!is.null(participant_col) && length(survey_columns)) {
            participants <- as.character(data[[participant_col]])
            unique_participants <- sort_values(participants)
            missing_survey <- vapply(unique_participants,
                function(participant) {
                rows <- !is_blank(participants) & participants ==
                    participant
                values <- unlist(data[rows, survey_columns, drop = FALSE],
                    use.names = FALSE
                )
                all(is_blank(values))
            }, logical(1))
            summary$participants_missing_survey <- sum(missing_survey)
        }

        summary
    }

    summarize_cpg_detection <- function(data) {
        summary <- list(exists = is.data.frame(data))
        if (is.null(data) || !all(c("metric", "nCpGs") %in% names(data))) {
            summary$total_assessed <- NA_integer_
            summary$detected_any <- NA_integer_
            summary$detected_all <- NA_integer_
            summary$never_detected <- NA_integer_
            return(summary)
        }

        get_metric <- function(pattern) {
            idx <- grep(pattern, data$metric, ignore.case = TRUE)
            if (length(idx)) {
                as.numeric(data$nCpGs[[idx[[1]]]])
            } else {
                NA_real_
            }
        }

        summary$total_assessed <- get_metric("^Total CpGs assessed$")
        summary$detected_any <- get_metric("at least one sample")
        summary$detected_all <- get_metric("all samples")
        summary$never_detected <- get_metric("never detected")
        summary
    }

    summarize_sample_detection <- function(data) {
        summary <- list(exists = is.data.frame(data))
        if (is.null(data) || !"pDetected" %in% names(data)) {
            summary$n_samples <- NA_integer_
            summary$min_p_detected <- NA_real_
            summary$median_p_detected <- NA_real_
            summary$max_p_detected <- NA_real_
            return(summary)
        }

        p_values <- as.character(data$pDetected)
        numeric_idx <- is_numeric_like(p_values)
        p_detected <- rep(NA_real_, length(p_values))
        p_detected[numeric_idx] <- as.numeric(p_values[numeric_idx])
        summary$n_samples <- nrow(data)
        if (any(!is.na(p_detected))) {
            summary$min_p_detected <- min(p_detected, na.rm = TRUE)
            summary$median_p_detected <- stats::median(p_detected,
                na.rm = TRUE
            )
            summary$max_p_detected <- max(p_detected, na.rm = TRUE)
        } else {
            summary$min_p_detected <- NA_real_
            summary$median_p_detected <- NA_real_
            summary$max_p_detected <- NA_real_
        }
        summary
    }

    extract_sva_log_summary <- function(log_path) {
        summary <- list(k = NA_integer_, percent_variation = NA_real_)
        if (!file.exists(log_path)) {
            return(summary)
        }

        text <- paste(readLines(log_path, warn = FALSE, encoding = "UTF-8"),
            collapse = " "
        )
        match <- regexec("([0-9]+)\\s+surrogate variables explain\\s+([0-9.]+)\\s*%",
            text,
            perl = TRUE
        )
        values <- regmatches(text, match)[[1]]
        if (length(values) == 3L) {
            summary$k <- as.integer(values[[2]])
            summary$percent_variation <- as.numeric(values[[3]])
        }
        summary
    }

    summarize_sva <- function(pheno_path, log_path) {
        pheno_data <- safe_read_table_file(pheno_path)
        log_summary <- extract_sva_log_summary(log_path)
        summary <- list(
            path = slash(pheno_path), exists = !is.null(pheno_data),
            n_samples = NA_integer_, n_surrogate_variables = NA_integer_,
            surrogate_variables = character(), log_k = log_summary$k,
            percent_variation = log_summary$percent_variation,
            sentrix_id_levels = NA_integer_,
                sentrix_position_levels = NA_integer_
        )

        if (!is.na(log_summary$k)) {
            summary$n_surrogate_variables <- log_summary$k
            summary$surrogate_variables <- paste0("SV", seq_len(log_summary$k))
        }

        if (!is.null(pheno_data)) {
            summary$n_samples <- nrow(pheno_data)
            sentrix_id_col <- pick_column(pheno_data, c(
                "Sentrix_ID",
                "SentrixID", "Sentrix Id"
            ))
            sentrix_position_col <- pick_column(pheno_data, c(
                "Sentrix_Position",
                "SentrixPosition", "Sentrix Position"
            ))
            if (!is.null(sentrix_id_col)) {
                summary$sentrix_id_levels <-
                    length(sort_values(pheno_data[[sentrix_id_col]]))
            }
            if (!is.null(sentrix_position_col)) {
                summary$sentrix_position_levels <-
                    length(sort_values(pheno_data[[sentrix_position_col]]))
            }
        }

        summary
    }

    summarize_logs <- function(log_assets) {
        labels <- c(
            methylation = "Methylation Analysis", data = "Data Preparation",
            batch = "Batch Effect", glm = "GLM", lme = "LME"
        )
        rows <- lapply(names(log_assets), function(name) {
            asset <- log_assets[[name]]
            line_count <- if (isTRUE(asset$exists)) {
                length(readLines(asset$source_path,
                    warn = FALSE,
                    encoding = "UTF-8"
                ))
            } else {
                NA_integer_
            }
            data.frame(
                key = name, label = labels[[name]],
                    exists = isTRUE(asset$exists),
                source_path = asset$source_path, line_count = line_count,
                stringsAsFactors = FALSE
            )
        })
        rows <- do.call(rbind, rows)
        list(
            rows = rows, total = nrow(rows), found = sum(rows$exists),
            total_lines = sum(rows$line_count, na.rm = TRUE)
        )
    }

    make_data_notes <- function(summary) {
        if (!isTRUE(summary$exists)) {
            return("The source CSV was not found.")
        }

        notes <- c(sprintf(
            "The table has %s rows and %s columns.",
            format_count(summary$n_rows), format_count(summary$n_cols)
        ))

        if (!is.null(summary$participant_col) &&
            !is.null(summary$timepoint_col)) {
            notes <- c(notes, sprintf(
                "Using `%s` as the participant identifier and `%s` as the timepoint column, the file contains %s unique %s across timepoint values %s.",
                summary$participant_col, summary$timepoint_col,
                format_count(summary$n_participants), plural(
                    summary$n_participants,
                    "participant"
                ), collapse_values(summary$timepoints)
            ))
        } else if (!is.null(summary$participant_col)) {
            notes <- c(notes, sprintf(
                "Using `%s` as the participant identifier, the file contains %s unique %s. No Timepoint column was detected.",
                summary$participant_col, format_count(summary$n_participants),
                plural(summary$n_participants, "participant")
            ))
        } else {
            notes <- c(notes,
                "No biological participant identifier column was detected. Technical identifiers such as sample, Sentrix, array, and slide IDs are not treated as participant IDs.")
        }

        notes
    }

    describe_figure_title <- function(title) {
        title <- sub("^Figure\\s+[0-9]+:\\s*", "", title)
        title <- sub("\\.$", "", title)
        first <- substr(title, 1L, 1L)
        rest <- substr(title, 2L, nchar(title))
        paste0(tolower(first), rest)
    }

    highlight_report_label <- function(title) {
        sub("^((?:Figure|Table)\\s+[0-9]+)([:.])", "`\\1`\\2",
            title,
            perl = TRUE
        )
    }

    make_figure_notes <- function(items, figure_titles, section_label,
                                    figure_descriptions = NULL) {
        n_items <- nrow(items)
        if (!n_items) {
            return(sprintf(
                "No supported figures were found for the %s tab.",
                section_label
            ))
        }

        vapply(seq_len(n_items), function(idx) {
            title <- if (length(figure_titles) >= idx &&
                nzchar(figure_titles[[idx]])) {
                figure_titles[[idx]]
            } else if ("title" %in% names(items) &&
                nzchar(items$title[[idx]])) {
                items$title[[idx]]
            } else {
                sprintf("Figure %s", idx)
            }
            highlighted_title <- highlight_report_label(title)
            description <- if (!is.null(figure_descriptions) &&
                length(figure_descriptions) >= idx &&
                    nzchar(figure_descriptions[[idx]])) {
                figure_descriptions[[idx]]
            } else {
                describe_figure_title(title)
            }
            description <- sub("\\.$", "", description)
            if (grepl("^(presents|shows|summarises|compares)\\b",
                description,
                ignore.case = TRUE
            )) {
                sprintf("%s %s.", highlighted_title, description)
            } else {
                sprintf(
                    "%s presents %s.", highlighted_title,
                    description
                )
            }
        }, character(1))
    }

    make_metrics_notes <- function(items, figure_titles, figure_descriptions,
                                    data_summary) {
        make_figure_notes(items, figure_titles, "Metrics", figure_descriptions)
    }

    make_quality_control_notes <- function(items, figure_titles,
                                            figure_descriptions, cpg_summary,
                                                sample_summary) {
        notes <- make_figure_notes(
            items, figure_titles, "Quality Control",
            figure_descriptions
        )

        if (isTRUE(cpg_summary$exists)) {
            notes <- c(notes, sprintf(
                "`Table 1` reports %s assessed CpGs, %s detected in all samples, and %s never detected.",
                format_count(cpg_summary$total_assessed),
                    format_count(cpg_summary$detected_all),
                format_count(cpg_summary$never_detected)
            ))
        }

        if (isTRUE(sample_summary$exists)) {
            notes <- c(notes, sprintf(
                "`Table 2` shows sample detection percentages ranging from %s%% to %s%% across %s %s.",
                format_decimal(sample_summary$min_p_detected),
                format_decimal(sample_summary$max_p_detected),
                format_count(sample_summary$n_samples), plural(
                    sample_summary$n_samples,
                    "sample"
                )
            ))
        }

        notes
    }

    make_batch_effect_notes <- function(items, figure_titles,
                                        figure_descriptions, sva_summary) {
        notes <- make_figure_notes(
            items, figure_titles, "Batch Effect",
            figure_descriptions
        )

        if (!is.na(sva_summary$log_k) &&
            !is.na(sva_summary$percent_variation)) {
            notes <- c(notes, sprintf(
                "The SVA log reports %s surrogate %s explaining %s%% of data variation.",
                format_count(sva_summary$log_k), plural(
                    sva_summary$log_k,
                    "variable"
                ), format_decimal(sva_summary$percent_variation)
            ))
        }

        notes
    }

    make_logs_notes <- function(log_summary, lme_label = "LME Analysis") {
        descriptions <- c(
            methylation = "displays methylation preprocessing, including IDAT loading, normalisation, filtering, and cell composition estimation",
            data = "displays phenotype preparation, timepoint splitting, and methylation matrix export steps",
            batch = "displays the hidden-effect and surrogate-variable analysis workflow",
            glm = "displays the generalised linear model workflow and CpG annotation steps",
            lme = "displays the linear mixed-effects model workflow and CpG annotation steps"
        )
        labels <- c(
            methylation = "Methylation Analysis", data = "Data Preparation",
            batch = "Batch Effect", glm = "GLM Analysis", lme = lme_label
        )

        notes <- c(
            "The Logs tab displays workflow log files generated by the pipeline, with one section per analysis stage.",
            vapply(log_summary$rows$key, function(key) {
                status <-
                    if (isTRUE(log_summary$rows$exists[log_summary$rows$key ==
                    key])) {
                    "is available"
                } else {
                    "was not found"
                }
                sprintf(
                    "`%s` %s and %s.", labels[[key]], status,
                    descriptions[[key]]
                )
            }, character(1))
        )

        notes
    }

    html_paragraph <- function(text) {
        sprintf("<p>%s</p>", html_escape(text))
    }

    html_section <- function(title, paragraphs) {
        c(sprintf("<h3>%s</h3>", html_escape(title)),
            html_paragraph(paragraphs))
    }

    strip_inline_markdown <- function(text) {
        text <- gsub("`([^`]+)`", "\\1", text, perl = TRUE)
        text <- gsub("<[^>]+>", "", text, perl = TRUE)
        text <- gsub("\\s+", " ", text)
        trimws(text)
    }

    sentence_case <- function(text) {
        if (!nzchar(text)) {
            return(text)
        }
        paste0(toupper(substr(text, 1L, 1L)), substr(
            text, 2L,
            nchar(text)
        ))
    }

    tab_report_paragraph <- function(notes, fallback) {
        notes <- notes[!grepl("dnaepico-model-formulas", notes,
            fixed = TRUE
        )]
        notes <- strip_inline_markdown(notes)
        notes <- notes[nzchar(notes)]
        if (!length(notes)) {
            return(sprintf("%s.", sentence_case(fallback)))
        }

        paste(notes, collapse = " ")
    }

    build_report_page <- function(project_name, project_dir,
                                    data_notes, enmix_notes,
                                        quality_control_notes,
                                        batch_effect_notes,
                                    metrics_notes, glm_notes, lme_notes,
                                        logs_notes,
                                        lme_label = "LME Analysis",
                                        model_sections = c("glm", "lme")) {
        overview <-
            sprintf("This report is generated automatically from the available datasets, figures, tables, and workflow logs for the project.")

        data_paragraph <- tab_report_paragraph(data_notes,
            "displays the phenotype data preview and detected participant/timepoint fields")
        enmix_paragraph <- tab_report_paragraph(
            enmix_notes,
            "displays ENmix control plots and interpretation notes"
        )
        quality_control_paragraph <- tab_report_paragraph(
            quality_control_notes,
            "displays quality-control figures and detection summary tables"
        )
        batch_effect_paragraph <- tab_report_paragraph(
            batch_effect_notes,
            "displays SVA and batch-effect figures"
        )
        metrics_paragraph <- tab_report_paragraph(
            metrics_notes,
            "displays post-filtering methylation metric figures"
        )
        glm_paragraph <- tab_report_paragraph(glm_notes,
            "displays the annotated generalised linear model results table")
        lme_paragraph <- tab_report_paragraph(lme_notes,
            "displays the annotated linear mixed-effects model results table")
        logs_paragraph <- tab_report_paragraph(logs_notes,
            "displays workflow log files for each analysis stage")
        model_report_sections <- c(
            if ("glm" %in% model_sections) {
                html_section("GLM", glm_paragraph)
            },
            if ("lme" %in% model_sections) {
                html_section(sub(" Analysis$", "", lme_label), lme_paragraph)
            }
        )

        c(
            "---", "title: \"Report\"", "format:", "  html:", "    css:",
            "      - assets/qpasst.css", "page-layout: full",
            "body-classes: qpasst-report-page", "---", "",
                "<div class=\"qpasst-report-viewport\">",
            "<div class=\"qpasst-report-card\">",
                "<div class=\"qpasst-report-sheet\">",
            sprintf("<h2>%s Report</h2>", html_escape(project_name)),
            html_paragraph(overview), html_section("Data", data_paragraph),
            html_section("ENmix QC", enmix_paragraph), html_section(
                "Quality Control",
                quality_control_paragraph
            ), html_section(
                "Batch Effect",
                batch_effect_paragraph
            ), html_section(
                "Metrics",
                metrics_paragraph
            ), model_report_sections,
            html_section("Logs", logs_paragraph), "</div>", "</div>",
            "</div>"
        )
    }

    sidebar_lines <- function(notes) {
        bullet_lines <- if (length(notes)) {
            paste0("- ", notes)
        } else {
            "- Add notes here."
        }

        c(
            "## {.sidebar width=\"320px\"}", "### Notes", "", bullet_lines,
            ""
        )
    }

    compose_page <- function(title, notes, body_lines, body_classes = NULL) {
        front_matter <- c(
            "---", sprintf("title: \"%s\"", title),
            "format: dashboard"
        )

        if (!is.null(body_classes) && nzchar(body_classes)) {
            front_matter <- c(front_matter, sprintf(
                "body-classes: %s",
                body_classes
            ))
        }

        c(
            front_matter, "---", "", sidebar_lines(notes), "## Column",
            "", body_lines
        )
    }

    prepared_report <- prepareDnamReportInputs(
        outputDir = project_dir,
        qcDir = qc_dir, preprocessingDir = preprocessing_dir,
        postprocessingDir = postprocessing_dir, svaDir = sva_dir,
        glmDir = glm_dir, lmeDir = lme_dir, figDir = fig_dir,
        verbose = FALSE, logs = FALSE, logDir = logs_dir
    )

    enmix_items <- copy_figure_assets(qc_dir, asset_subdir = "enmix-qc")
    metrics_items <- copy_figure_assets(postprocessing_dir,
        asset_subdir = "metrics")
    qc_items <- copy_figure_assets(preprocessing_dir,
        asset_subdir = "quality-control")
    batch_items <- copy_figure_assets(sva_dir, asset_subdir = "batch-effects")

    log_assets <- list(
        methylation = copy_log_asset(file.path(
            logs_dir,
            "log_preprocessingMinfiEwasWater.txt"
        ), "preprocessingMinfiEwasWater.txt"),
        data = copy_log_asset(
            file.path(logs_dir, "log_preprocessingPheno.txt"),
            "preprocessingPheno.txt"
        ), batch = copy_log_asset(file.path(
            logs_dir,
            "log_svaEnmix.txt"
        ), "svaEnmix.txt")
    )
    if (include_glm) {
        log_assets$glm <- copy_first_existing_log_asset(file.path(
            logs_dir,
            "log_methylationGLM.txt"
        ), "methylationGLM.txt")
    }
    if (include_lme) {
        log_assets$lme <- copy_first_existing_log_asset(file.path(
            logs_dir,
            "log_methylationLME.txt"
        ), "methylationLME.txt")
    }

    data_summary <- summarize_dataset(pheno_file)
    detection_tables <- read_detection_tables(
        detp_path = detp_path,
        threshold = detPThreshold, cpg_path = cpg_detection_path,
        sample_path = sample_detection_path
    )
    cpg_detection_summary <- summarize_cpg_detection(detection_tables$cpg)
    sample_detection_summary <-
        summarize_sample_detection(detection_tables$sample)
    sva_summary <- summarize_sva(pheno_path = pheno_file, log_path = file.path(
        logs_dir,
        "log_svaEnmix.txt"
    ))
    log_summary <- summarize_logs(log_assets)

    enmix_figure_titles <- enmix_items$title
    enmix_notes <- c(paste(
        "`ENmix` produces control plots similar to those generated by Illumina's GenomeStudio software.",
        "Infinium controls are not evaluated using absolute intensity values.",
        "Instead, we interpret them based on their expected signal range across the dataset.",
        "This approach is intended to provide greater flexibility and reliability when accounting for natural biological variation."
    ))
    metrics_figure_titles <- c(
        "Figure 1: The Density Distributions of Beta and M-values Across Samples Post-Filtering Methylation Plot by Timepoint",
        "Figure 2: Principal Component Analysis Post-Filtering Methylation Plot by Timepoint"
    )
    metrics_figure_descriptions <- c(
        "beta-value and M-value density distributions across samples after filtering, grouped by timepoint",
        "a principal component analysis plot after filtering, grouped by timepoint"
    )
    quality_control_metadata <- data.frame(
        pattern = c(
            "densityBeta",
            "detection_pvalues", "quality_control", "sexClinical",
            "sexComparison_RawNorm", "sexPrediction"
        ), title = c(
            "Figure 1: Beta Value Density Distribution (MSET)",
            "Figure 2: Detection P-values (RGSET)",
                "Figure 3: Quality Control Plot (MSET)",
            "Figure 4: Clinical Sex Distribution (GSET)",
                "Figure 5: Sex Comparison Before and After Normalisation (MSETF)",
            "Figure 6: Sex Prediction Plot (GSET)"
        ), description = c(
            "presents the distribution of methylation beta values across samples, allowing visual assessment of overall methylation patterns and potential sample outliers",
            "shows the distribution of detection p-values across samples, providing a quality check for probe signal reliability and identifying samples with poor detection performance",
            "presents sample-level quality metrics based on methylated and unmethylated signal intensities, helping to identify low-quality samples or technical outliers",
            "summarises the sex information recorded in the clinical or phenotype data, supporting checks of sample annotation consistency",
            "compares sex-related methylation patterns in raw and normalised data, helping to confirm whether normalisation preserves expected biological differences",
            "presents methylation-based sex prediction across samples, allowing comparison with recorded clinical sex and detection of possible sample swaps or annotation errors"
        ),
        stringsAsFactors = FALSE
    )
    quality_control_figure_titles <- qc_items$title
    quality_control_figure_descriptions <- rep("", nrow(qc_items))
    for (metadata_idx in seq_len(nrow(quality_control_metadata))) {
        item_idx <- grep(quality_control_metadata$pattern[[metadata_idx]],
            qc_items$original_name,
            ignore.case = TRUE
        )
        if (length(item_idx)) {
            quality_control_figure_titles[item_idx] <-
                quality_control_metadata$title[[metadata_idx]]
            quality_control_figure_descriptions[item_idx] <-
                quality_control_metadata$description[[metadata_idx]]
        }
    }
    batch_effect_figure_titles <- c(
        "Figure 1: Distribution of Samples Across the First Two Surrogate Variables by SentrixID",
        "Figure 2: Pairwise Distribution of the First Three Surrogate Variables by SentrixID and SentrixPosition",
        "Figure 3: Distribution of Samples Across the First Two Surrogate Variables by SentrixPosition"
    )
    batch_effect_figure_descriptions <- c(
        "sample distribution across the first two surrogate variables, coloured by Sentrix ID",
        "pairwise distributions of the first three surrogate variables by Sentrix ID and Sentrix position",
        "sample distribution across the first two surrogate variables by Sentrix position"
    )
    data_table_assets <- prepare_csv_table_assets(
        data_path = pheno_file,
        var_prefix = "phenotype_data", front_columns = c(
            if (is.null(data_summary$participant_col)) "" else data_summary$participant_col,
            if (is.null(data_summary$timepoint_col)) "" else data_summary$timepoint_col
        )
    )
    unrequested_table_assets <- function(analysis) {
        list(
            ok = FALSE, error = "The model section was not requested.",
            source_path = "", analysis = analysis, metadata = list(),
            source_mode = NULL
        )
    }
    glm_table_assets <- if (include_glm) {
        prepare_xlsx_table_assets(
            data_path = glm_table_path,
            sheet = "annotatedGLM", var_prefix = "glm_results",
            analysis = "glm",
            model_log_path = file.path(logs_dir, "log_methylationGLM.txt")
        )
    } else {
        unrequested_table_assets("glm")
    }
    lme_table_assets <- if (include_lme) {
        prepare_xlsx_table_assets(
            data_path = lme_table_path,
            sheet = "annotatedLME", var_prefix = "lme_results",
            analysis = "lme",
            model_log_path = file.path(logs_dir, "log_methylationLME.txt")
        )
    } else {
        unrequested_table_assets("lme")
    }
    emitLogMinfiEwasWater(
        c(
            paste(
                "GLM report table source:",
                if (!include_glm) {
                    "not requested"
                } else if (is.null(glm_table_assets$source_mode)) {
                    "unavailable"
                } else {
                    glm_table_assets$source_mode
                }
            ),
            paste("LME report table source:",
                if (!include_lme) {
                    "not requested"
                } else if (is.null(lme_table_assets$source_mode)) {
                    "unavailable"
                } else {
                    lme_table_assets$source_mode
                })
        ),
        verbose = FALSE, log_path = log_path
    )

    lme_backend <- lme_table_assets$metadata$backend
    if (is.null(lme_backend) || !nzchar(lme_backend)) {
        lme_backend <- "lme4"
    }
    lme_backend <- tolower(lme_backend)
    lme_analysis_label <- if (identical(lme_backend, "nlme")) {
        "nlme Analysis"
    } else {
        "LME Analysis"
    }
    lme_interaction_term <- lme_table_assets$metadata$interaction_term
    has_lme_interaction <- !is.null(lme_interaction_term) &&
        nzchar(lme_interaction_term) && !tolower(lme_interaction_term) %in%
        c("none", "null", "na")

    glm_table_title <-
        "Table 1. Generalised Linear Model Results and Genomic Annotation of CpG Sites by Phenotype(s)"
    lme_table_title <- if (has_lme_interaction) {
        sprintf(
            "Table 1. Linear Mixed-Effects Model Results and Genomic Annotation of CpG Sites by Phenotype(s) and %s",
            lme_interaction_term
        )
    } else {
        "Table 1. Linear Mixed-Effects Model Results and Genomic Annotation of CpG Sites by Phenotype(s)"
    }

    describe_available_annotations <- function(columns) {
        notes <- character()
        if (all(c("IlmnID", "Name") %in% columns)) {
            notes <- c(notes,
                "The `IlmnID` and `Name` columns identify each Illumina CpG probe.")
        } else if ("IlmnID" %in% columns) {
            notes <- c(notes,
                "The `IlmnID` column identifies each Illumina CpG probe.")
        }
        if (all(c("chr", "pos") %in% columns)) {
            notes <- c(notes,
                "The `chr` and `pos` columns give the probe's genomic position.")
        }
        if (any(c("UCSC_RefGene_Group", "UCSC_RefGene_Name") %in%
            columns)) {
            notes <- c(notes,
                "Available `UCSC_RefGene` columns provide gene-related annotation.")
        }
        if ("Relation_to_Island" %in% columns) {
            notes <- c(notes,
                "`Relation_to_Island` reports the probe's CpG-island context.")
        }
        if ("GencodeV41_Group" %in% columns) {
            notes <- c(notes,
                "`GencodeV41_Group` provides GENCODE version 41 gene-region annotation.")
        }
        notes
    }

    glm_table_description <- if (isTRUE(glm_table_assets$ok)) {
        glm_p_columns <- grep("P\\.Value$|P\\.value$", glm_table_assets$columns,
            value = TRUE
        )
        c(
            "CpGs were analysed using a generalised linear model fitted with `glm2`.",
            build_model_formula_notes(
                glm_table_assets$metadata$formula_records,
                "GLM"
            ), if (length(glm_p_columns)) {
                paste0(
                    "The result p-value column(s) are `",
                    paste(glm_p_columns, collapse = "`, `"), "`."
                )
            }, describe_available_annotations(glm_table_assets$columns)
        )
    } else {
        "The annotated GLM result workbook was not available when this report was generated."
    }

    lme_table_description <- if (isTRUE(lme_table_assets$ok)) {
        lme_p_columns <- grep("P\\.Value$|P\\.value$", lme_table_assets$columns,
            value = TRUE
        )
        fitting_description <- if (identical(lme_backend, "nlme")) {
            "CpGs were analysed using `nlme::lme()`."
        } else {
            "CpGs were analysed using `lmerTest::lmer()` with the `lme4` mixed-effects framework."
        }
        correlation_structure <- lme_table_assets$metadata$correlation_structure
        correlation_variable <- lme_table_assets$metadata$correlation_variable
        correlation_description <- if (identical(
            lme_backend,
            "nlme"
        ) && !is.null(correlation_structure) && nzchar(correlation_structure) &&
            !tolower(correlation_structure) %in% c(
                "none", "null",
                "na"
            )) {
            sprintf(
                "The `%s` residual correlation structure orders repeated observations using `%s`.",
                correlation_structure, if (is.null(correlation_variable) ||
                    !nzchar(correlation_variable)) {
                    "the configured correlation variable"
                } else {
                    correlation_variable
                }
            )
        } else {
            character()
        }
        interaction_description <- if (has_lme_interaction) {
            sprintf(
                "The reported phenotype terms test their interaction with `%s`.",
                lme_interaction_term
            )
        } else {
            "No phenotype-by-time interaction was fitted; the reported p-value column(s) correspond to phenotype coefficient(s)."
        }
        formula_description <- build_model_formula_notes(
            lme_table_assets$metadata$formula_records,
            if (identical(lme_backend, "nlme")) {
                "nlme"
            } else {
                "LME"
            }
        )
        c(
            fitting_description, correlation_description,
                interaction_description,
            formula_description, if (length(lme_p_columns)) {
                paste0(
                    "The result p-value column(s) are `",
                    paste(lme_p_columns, collapse = "`, `"), "`."
                )
            }, describe_available_annotations(lme_table_assets$columns)
        )
    } else {
        "The annotated LME result workbook was not available when this report was generated."
    }

    data_notes <- make_data_notes(data_summary)
    metrics_notes <- make_metrics_notes(
        metrics_items, metrics_figure_titles,
        metrics_figure_descriptions, data_summary
    )
    quality_control_notes <- make_quality_control_notes(
        qc_items,
        quality_control_figure_titles, quality_control_figure_descriptions,
        cpg_detection_summary, sample_detection_summary
    )
    batch_effect_notes <- make_batch_effect_notes(
        batch_items,
        batch_effect_figure_titles, batch_effect_figure_descriptions,
        sva_summary
    )
    glm_notes <- glm_table_description
    lme_notes <- lme_table_description
    logs_notes <- make_logs_notes(log_summary, lme_label = lme_analysis_label)

    if (file.exists(logo_source_path)) {
        dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
        ok <- file.copy(logo_source_path, file.path(
            assets_dir,
            "dnaEPICO.svg"
        ), overwrite = TRUE)
        if (!ok) {
            warning("Failed to copy navbar logo: ", logo_source_path)
        }
    }
    if (file.exists(viewer_js_source_path)) {
        dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
        ok <- file.copy(viewer_js_source_path, file.path(
            assets_dir,
            "cpg-viewer.js"
        ), overwrite = TRUE)
        if (!ok) {
            warning(
                "Failed to copy CpG viewer JavaScript: ",
                viewer_js_source_path
            )
        }
    } else {
        warning("CpG viewer JavaScript was not found in the installed package.")
    }

    model_navbar <- c(
        if (include_glm) {
            c("      - href: glm.qmd", "        text: \"GLM Analysis\"")
        },
        if (include_lme) {
            c("      - href: lme.qmd", sprintf(
                "        text: \"%s\"",
                lme_analysis_label
            ))
        }
    )
    quarto_yml <- c(
        "project:", "  type: website", "  output-dir: docs",
        "  resources:", "    - assets/", "", "website:", sprintf(
            "  title: \"%s\"",
            project_name
        ), "  navbar:", "    logo: assets/dnaEPICO.svg",
        "    left:", "      - href: index.qmd", "        text: \"Data\"",
        "      - href: enmix-qc.qmd", "        text: \"ENmix QC\"",
        "      - href: quality-control.qmd",
            "        text: \"Quality Control\"",
        "      - href: batch-effect.qmd", "        text: \"Batch Effect\"",
        "      - href: metrics.qmd", "        text: \"Metrics\"",
        model_navbar, "      - href: report.qmd",
        "        text: \"Report\"", "      - href: logs.qmd",
        "        text: \"Logs\"", "    search: true",
            "  page-navigation: false",
        "  bread-crumbs: false", "  reader-mode: false", "",
        "format:", "  dashboard:", "    theme: cosmo", "    css:",
        "      - assets/qpasst.css", "    toc: false", "", "execute:",
        "  warning: false", "  message: false"
    )

    site_css <- c(
        ":root {", "  --qpasst-border: #e5e7eb;", "  --qpasst-bg: #d5d9df;",
        "  --qpasst-navbar-bg: #f4f6f8;", "  --qpasst-card-bg: #ffffff;",
        "  --qpasst-accent: #0d6efd;", "}", "", "body {",
            "  background: #ffffff;",
        "}", "", ".navbar,", ".quarto-navbar,", "#quarto-header .navbar {",
        "  background: var(--qpasst-navbar-bg) !important;",
        "  border-bottom: 1px solid var(--qpasst-border) !important;",
        "  box-shadow: none !important;", "}", "", ".navbar-brand,",
        ".navbar-title {", "  color: #3f444a !important;",
            "  font-weight: 600;",
        "}", "", "#quarto-header .navbar-brand-container {",
        "  display: flex;", "  align-items: center;", "}", "",
        "#quarto-header .navbar-brand.navbar-brand-logo {",
            "  display: flex !important;",
        "  align-items: center;", "  padding-top: 0;", "  padding-bottom: 0;",
        "  margin-right: 0.35rem;", "}", "",
            "#quarto-header .navbar-brand.navbar-brand-logo img,",
        "#quarto-header .navbar-brand.navbar-brand-logo .navbar-logo,",
        "#quarto-header .navbar-brand.navbar-brand-logo .navbar-logo-image {",
        "  display: block;", "  height: 48px !important;", "  width: auto;",
        "  max-height: 48px !important;", "  max-width: 48px !important;",
        "  object-fit: contain;", "}", "",
            "#quarto-header .navbar-brand:not(.navbar-brand-logo) {",
        "  display: none !important;", "}", "", ".navbar-nav .nav-link {",
        "  color: #495057 !important;", "  font-weight: 700;",
        "}", "", ".navbar-nav .nav-link.active,",
            ".navbar-nav .show > .nav-link,",
        ".navbar-nav .nav-link:hover {",
            "  color: var(--qpasst-accent) !important;",
        "}", "", "#quarto-dashboard-header {", "  display: none !important;",
        "}", "", ".quarto-dashboard .bslib-sidebar-layout > .main,",
        ".quarto-dashboard .bslib-sidebar-layout .sidebar-content {",
        "  overflow-y: auto !important;", "}", "",
            ".quarto-dashboard .sidebar-content > .bslib-grid,",
        ".quarto-dashboard .sidebar-content > .bslib-grid > .bslib-grid,",
        ".quarto-dashboard .sidebar-content > .bslib-grid > .bslib-grid > .bslib-grid {",
        "  grid-auto-rows: max-content !important;", "}", "",
        ".quarto-dashboard .bslib-card {", "  min-height: 0 !important;",
        "  margin-bottom: 1.5rem;", "}", "",
            ".quarto-dashboard .bslib-card .card-body {",
        "  height: auto !important;", "  overflow: visible !important;",
        "  padding: 1.5rem !important;", "}", "",
            "body.qpasst-logs-page .sidebar-content > .bslib-grid > .bslib-grid {",
        "  grid-template-rows: minmax(3em, max-content) minmax(3em, max-content) minmax(3em, max-content) !important;",
        "}", "", "body.qpasst-logs-page .panel-tabset.bslib-grid {",
        "  grid-template-columns: minmax(3em, 1fr) !important;",
        "  grid-auto-rows: max-content !important;", "}", "",
        "body.qpasst-logs-page .bslib-card {",
            "  min-height: calc(100vh - 12rem) !important;",
        "}", "", "body.qpasst-logs-page .bslib-card .card-body {",
        "  overflow: auto !important;", "}", "",
            "body.qpasst-logs-page pre.text {",
        "  margin: 0;", "  max-height: calc(100vh - 16rem);",
        "  overflow: auto;", "  white-space: pre;", "}", "",
        "body.qpasst-report-page #title-block-header {",
            "  display: none !important;",
        "}", "", "body.qpasst-report-page #quarto-content,",
        "body.qpasst-report-page main.content,",
            "body.qpasst-report-page #quarto-document-content,",
        "body.qpasst-report-page .page-columns,",
            "body.qpasst-report-page .page-full {",
        "  width: 100% !important;", "  max-width: none !important;",
        "}", "", "body.qpasst-report-page main.content {",
            "  margin: 0 !important;",
        "  padding: 0.75rem !important;", "}", "",
            "body.qpasst-report-page .qpasst-report-viewport {",
        "  width: 100%;", "  min-height: calc(100vh - 7rem);",
        "}", "", "body.qpasst-report-page .qpasst-report-card {",
        "  width: 100%;", "  min-height: calc(100vh - 7rem);",
        "  border: 1px solid var(--qpasst-border);", "  border-radius: 12px;",
        "  background: var(--qpasst-card-bg);",
            "  box-shadow: 0 4px 14px rgba(15, 23, 42, 0.05);",
        "  overflow-y: auto;", "  padding: 1.25rem 1.5rem;",
        "}", "", "body.qpasst-report-page .qpasst-report-sheet {",
        "  width: 100%;", "  min-height: 297mm;", "  margin: 0;",
        "  padding: 0.5rem 0 2rem;", "}", "",
            "body.qpasst-report-page .qpasst-report-sheet h2,",
        "body.qpasst-report-page .qpasst-report-sheet h3 {",
        "  color: #2f3540;", "  margin-top: 0;", "}", "",
            "body.qpasst-report-page .qpasst-report-sheet p {",
        "  font-size: 1.02rem;", "  line-height: 1.75;",
            "  margin-bottom: 1rem;",
        "}", "", ".dashboard .sidebar,", ".quarto-dashboard .sidebar,",
        ".sidebar {", "  background: var(--qpasst-bg) !important;",
        "  border-right: 1px solid var(--qpasst-border) !important;",
        "  overflow-x: hidden !important;", "}", "", ".sidebar h3,",
        ".sidebar .h3 {", "  color: #2f3540 !important;", "}",
        "", ".sidebar p,", ".sidebar li,", ".sidebar code {",
        "  font-size: 1.06rem;", "  line-height: 1.7;",
            "  overflow-wrap: anywhere;",
        "  word-break: break-word;", "}", "",
            "body.qpasst-data-page .sidebar li {",
        "  font-size: 1.08rem;", "  line-height: 1.65;", "}",
        "", ".sidebar li::marker {", "  color: var(--qpasst-accent);",
        "}", "", ".qpasst-gallery {", "  display: grid;",
            "  grid-template-columns: 1fr;",
        "  gap: 1.25rem;", "  margin-top: 0.5rem;", "}", "",
        ".qpasst-gallery-item {", "  margin: 0;",
            "  border: 1px solid var(--qpasst-border);",
        "  border-radius: 12px;", "  background: var(--qpasst-card-bg);",
        "  padding: 1rem;", "  box-shadow: 0 4px 14px rgba(15, 23, 42, 0.05);",
        "}", "", ".qpasst-gallery-item figcaption {",
            "  margin-bottom: 0.75rem;",
        "  line-height: 1.5;", "}", "", ".qpasst-gallery-item img {",
        "  display: block;", "  width: 50%;", "  height: auto;",
        "  border-radius: 8px;", "  background: #fff;", "}",
        "", ".qpasst-gallery-item code {", "  white-space: normal;",
        "  word-break: break-word;", "}", "", ".qpasst-figure-card {",
        "  margin: 0;", "}", "", ".qpasst-figure-card img {",
        "  display: block;", "  width: 100%;", "  height: auto;",
        "  border-radius: 8px;", "  background: #fff;", "}",
        "", ".qpasst-figure-card figcaption {", "  margin-top: 0.75rem;",
        "}", "", ".qpasst-download {", "  margin: 0;",
            "  padding: 0.75rem 0.9rem;",
        "  border-radius: 8px;", "  background: #f8fafc;",
            "  border: 1px dashed var(--qpasst-border);",
        "}", "", ".table {", "  font-size: 0.92rem;", "}", "",
        ".dnaepico-cpg-viewer {", "  display: flex !important;",
        "  flex-direction: column !important;", "  gap: 0.45rem;",
        "  width: 100%;", "  min-height: 0 !important;", "}",
        "", ".dnaepico-viewer-toolbar,", ".dnaepico-viewer-filter,",
        ".dnaepico-viewer-pagination,", ".dnaepico-viewer-downloads {",
        "  display: flex !important;", "  flex: 0 0 auto !important;",
        "  flex-direction: row !important;",
            "  align-items: center !important;",
        "  flex-wrap: wrap !important;", "  gap: 0.45rem;",
            "  min-height: 0 !important;",
        "}", "", ".dnaepico-viewer-toolbar label,",
            ".dnaepico-viewer-filter label,",
        ".dnaepico-viewer-pagination label {",
            "  display: inline-flex !important;",
        "  flex-direction: row !important;",
            "  align-items: center !important;",
        "  gap: 0.4rem;", "  margin: 0;", "}", "", ".dnaepico-viewer-filter {",
        "  padding: 0.4rem 0.5rem;", "  border: 1px solid #e2e8f0;",
        "  border-radius: 0.4rem;", "  background: #f8fafc;",
        "}", "", ".dnaepico-viewer-toolbar input {", "  min-width: 14rem;",
        "}", "", ".dnaepico-viewer-filter input {", "  width: 8rem;",
        "}", "", ".dnaepico-viewer-pagination input {", "  width: 4.5rem;",
        "}", "", ".dnaepico-viewer-toolbar input,",
            ".dnaepico-viewer-toolbar select,",
        ".dnaepico-viewer-filter input,", ".dnaepico-viewer-filter select,",
        ".dnaepico-viewer-pagination input {", "  border: 1px solid #ced4da;",
        "  border-radius: 0.375rem;", "  padding: 0.25rem 0.45rem;",
        "  background: #fff;", "}", "", ".dnaepico-viewer-toolbar button,",
        ".dnaepico-viewer-filter button,",
            ".dnaepico-viewer-pagination button {",
        "  border: 1px solid #adb5bd;", "  border-radius: 0.375rem;",
        "  padding: 0.25rem 0.55rem;", "  background: #fff;",
        "  color: #212529;", "}", "", ".dnaepico-viewer-toolbar button:hover,",
        ".dnaepico-viewer-filter button:hover:not(:disabled),",
        ".dnaepico-viewer-pagination button:hover:not(:disabled) {",
        "  background: #eef3f8;", "}", "",
            ".dnaepico-viewer-filter button:disabled,",
        ".dnaepico-viewer-pagination button:disabled {", "  opacity: 0.5;",
        "}", "", ".dnaepico-viewer-pagination {",
            "  justify-content: center !important;",
        "  font-size: 0.88rem;", "}", "",
            ".dnaepico-viewer-pagination button {",
        "  width: 2rem;", "  height: 1.85rem;", "  padding: 0;",
        "  font-size: 1rem;", "  line-height: 1;", "}", "",
            ".dnaepico-viewer-table-wrap {",
        "  display: block !important;", "  flex: 1 1 auto !important;",
        "  overflow: auto !important;",
            "  border: 1px solid var(--qpasst-border);",
        "  border-radius: 0.4rem;", "  height: clamp(28rem, 62vh, 52rem);",
        "  min-height: 28rem !important;", "  max-height: none;",
        "}", "", ".dnaepico-viewer-table {", "  margin-bottom: 0;",
        "  white-space: nowrap;", "}", "", ".dnaepico-viewer-table thead th {",
        "  position: sticky;", "  top: 0;", "  z-index: 1;",
        "  background: #f4f6f8;", "}", "", ".dnaepico-viewer-status {",
        "  display: block !important;", "  flex: 0 0 auto !important;",
        "  min-height: 0 !important;", "  color: #495057;", "}",
        "", ".dnaepico-viewer-filter [data-role='filter-summary'] {",
        "  color: #495057;", "  font-size: 0.9rem;", "}", "",
        ".sidebar .dnaepico-model-formula {", "  display: inline;",
        "  white-space: normal !important;", "  overflow-wrap: anywhere;",
        "  word-break: normal;", "}", "", ".sidebar .dnaepico-model-formulas {",
        "  margin: 0.35rem 0;", "}", "",
            ".sidebar .dnaepico-model-formulas summary {",
        "  cursor: pointer;", "  font-weight: 600;", "  color: #2f3540;",
        "}", "", ".sidebar .dnaepico-model-formula-list {", "  display: grid;",
        "  gap: 0.65rem;", "  margin-top: 0.65rem;", "}", "",
        ".sidebar .dnaepico-model-formula-item {",
            "  padding: 0.55rem 0.65rem;",
        "  border: 1px solid #e2e8f0;", "  border-radius: 0.4rem;",
        "  background: #f8fafc;", "}", "",
            ".sidebar .dnaepico-model-formula-item strong,",
        ".sidebar .dnaepico-model-formula-item code {", "  display: block;",
        "}", "", ".sidebar .dnaepico-model-formula-item code {",
        "  margin-top: 0.3rem;", "}", "", ".dnaepico-viewer-error {",
        "  color: #b02a37;", "}", "", "pre code {", "  white-space: pre-wrap;",
        "  word-break: break-word;", "}", "", "@media (max-width: 768px) {",
        "  .qpasst-gallery {", "    grid-template-columns: 1fr;",
        "  }", "  .dnaepico-viewer-toolbar input {", "    min-width: 10rem;",
        "  }", "}"
    )

    data_page <- compose_page(
        title = "Data", notes = data_notes,
        body_classes = "qpasst-data-page",
            body_lines = c(build_result_table_section(
            title = "Data Preview",
            table_assets = data_table_assets
        ))
    )

    enmix_page <- compose_page(
        title = "ENmix QC", notes = enmix_notes,
        body_lines = c(build_figure_cards(enmix_items,
            "No supported image files were found for the ENmix QC tab.",
            title_prefix = "Figure", show_filename = FALSE,
                figure_titles = enmix_figure_titles
        ))
    )

    metrics_page <- compose_page(
        title = "Metrics", notes = metrics_notes,
        body_lines = c(build_figure_cards(metrics_items,
            "No supported image files were found for the Metrics tab.",
            title_prefix = "Figure", show_filename = FALSE,
                figure_titles = metrics_figure_titles
        ))
    )

    quality_control_page <- compose_page(
        title = "Quality Control",
        notes = quality_control_notes, body_lines = c(
            build_figure_cards(qc_items,
                "No supported image files were found for the Quality Control tab.",
                title_prefix = "Figure", show_filename = FALSE,
                    figure_titles = quality_control_figure_titles
            ),
            build_data_frame_table_section(
                title = "Table 1: CpG Detection Summary",
                data = detection_tables$cpg,
                    empty_message = detection_table_warning(detection_tables),
                preview_rows = 10L
            ), build_data_frame_table_section(
                title = "Table 2: Sample Detection Summary",
                data = detection_tables$sample,
                    empty_message = detection_table_warning(detection_tables),
                preview_rows = 100L
            )
        )
    )

    batch_effect_page <- compose_page(
        title = "Batch Effect",
        notes = batch_effect_notes,
            body_lines = c(build_figure_cards(batch_items,
            "No supported image files were found for the Batch Effect tab.",
            title_prefix = "Figure", show_filename = FALSE,
                figure_titles = batch_effect_figure_titles
        ))
    )

    glm_page <- compose_page(
        title = "GLM Analysis", notes = glm_notes,
        body_lines = c(build_result_table_section(
            title = glm_table_title,
            table_assets = glm_table_assets
        ))
    )

    lme_page <- compose_page(
        title = lme_analysis_label, notes = lme_notes,
        body_lines = c(build_result_table_section(
            title = lme_table_title,
            table_assets = lme_table_assets
        ))
    )

    log_render_helper <- c(
        "```{r}", "#| echo: false", "#| include: false",
        "render_log_block <- function(path, label) {",
        "  if (!file.exists(path)) {",
        "    cat('::: {.callout-warning}\\n')",
        "    cat(sprintf('Log file not found for `%s`.\\n', label))",
        "    cat(':::\\n')", "    return(invisible(NULL))",
        "  }", "  lines <- readLines(path, warn = FALSE, encoding = 'UTF-8')",
        "  if (!length(lines)) {", "    cat('::: {.callout-note}\\n')",
        "    cat(sprintf('`%s` is empty.\\n', label))", "    cat(':::\\n')",
        "    return(invisible(NULL))", "  }", "  cat('```text\\n')",
        "  cat(paste(lines, collapse = '\\n'))", "  cat('\\n```\\n')",
        "}", "```"
    )
    base_log_cards <- c(
        "", "::: {.card title=\"Methylation Analysis\"}", "",
        "Displays the methylation preprocessing log, including IDAT loading, normalisation, filtering, and cell composition estimation steps.",
        "", "```{r}", "#| echo: false", "#| results: asis",
        sprintf(
            "render_log_block(%s, 'Methylation Analysis')",
            r_string(log_assets$methylation$asset_path)
        ), "```", ":::", "", "::: {.card title=\"Data Preparation\"}",
        "", "Displays the phenotype preparation log, including timepoint splitting and methylation matrix export steps.",
        "", "```{r}", "#| echo: false", "#| results: asis",
        sprintf(
            "render_log_block(%s, 'Data Preparation')",
            r_string(log_assets$data$asset_path)
        ), "```", ":::", "", "::: {.card title=\"Batch Effect\"}",
        "", "Displays the hidden-effect and surrogate-variable analysis log used for batch-effect assessment.",
        "", "```{r}", "#| echo: false", "#| results: asis",
        sprintf("render_log_block(%s, 'Batch Effect')",
            r_string(log_assets$batch$asset_path)),
        "```", ":::"
    )
    glm_log_card <- if (include_glm) {
        c(
            "", "::: {.card title=\"GLM Analysis\"}", "",
            "Displays the generalised linear model log, including phenotype association testing and CpG annotation steps.",
            "", "```{r}", "#| echo: false", "#| results: asis",
            sprintf("render_log_block(%s, 'GLM Analysis')",
                r_string(log_assets$glm$asset_path)),
            "```", ":::"
        )
    } else {
        character(0)
    }
    lme_log_card <- if (include_lme) {
        c(
            "", sprintf("::: {.card title=\"%s\"}",
                html_escape(lme_analysis_label)), "",
            if (has_lme_interaction) {
                sprintf(
                    "Displays the linear mixed-effects model log, including the phenotype interaction with %s and CpG annotation steps.",
                    lme_interaction_term
                )
            } else {
                "Displays the linear mixed-effects model log, including phenotype coefficient testing and CpG annotation steps."
            }, "", "```{r}", "#| echo: false", "#| results: asis",
            sprintf(
                "render_log_block(%s, %s)", r_string(log_assets$lme$asset_path),
                r_string(lme_analysis_label)
            ), "```", ":::"
        )
    } else {
        character(0)
    }
    logs_page <- compose_page(
        title = "Logs", notes = logs_notes,
        body_classes = "qpasst-logs-page",
        body_lines = c(
            log_render_helper, base_log_cards, glm_log_card, lme_log_card
        )
    )

    report_page <- build_report_page(
        project_name = project_name,
        project_dir = project_dir, data_notes = data_notes,
            enmix_notes = enmix_notes,
        quality_control_notes = quality_control_notes,
            batch_effect_notes = batch_effect_notes,
        metrics_notes = metrics_notes, glm_notes = glm_notes,
        lme_notes = lme_notes, logs_notes = logs_notes,
            lme_label = lme_analysis_label,
        model_sections = model_sections
    )

    unlink(c(
        file.path(project_dir, c("data.qmd", "DNAm.html")),
        file.path(project_dir, "docs", "data.html")
    ), force = TRUE)
    unrequested_source_files <- c(
        if (!include_glm) file.path(project_dir, "glm.qmd"),
        if (!include_lme) file.path(project_dir, "lme.qmd")
    )
    unlink(unrequested_source_files, force = TRUE)
    remove_unrequested_rendered_assets <- function() {
        unrequested_files <- c(
            if (!include_glm) {
                c(
                    file.path(project_dir, "docs", "glm.html"),
                    file.path(
                        project_dir, "docs", "assets", "logs",
                        "methylationGLM.txt"
                    )
                )
            },
            if (!include_lme) {
                c(
                    file.path(project_dir, "docs", "lme.html"),
                    file.path(
                        project_dir, "docs", "assets", "logs",
                        "methylationLME.txt"
                    )
                )
            }
        )
        unrequested_directories <- c(
            if (!include_glm) {
                file.path(
                    project_dir, "docs", "assets", "results", "glm_results"
                )
            },
            if (!include_lme) {
                file.path(
                    project_dir, "docs", "assets", "results", "lme_results"
                )
            }
        )
        unlink(unrequested_files, force = TRUE)
        unlink(unrequested_directories, recursive = TRUE, force = TRUE)
        invisible(NULL)
    }
    remove_unrequested_rendered_assets()

    write_utf8(file.path(project_dir, "_quarto.yml"), quarto_yml)
    write_utf8(file.path(assets_dir, "qpasst.css"), site_css)
    write_utf8(file.path(project_dir, "index.qmd"), data_page)
    write_utf8(file.path(project_dir, "enmix-qc.qmd"), enmix_page)
    write_utf8(file.path(project_dir, "metrics.qmd"), metrics_page)
    write_utf8(
        file.path(project_dir, "quality-control.qmd"),
        quality_control_page
    )
    write_utf8(file.path(project_dir, "batch-effect.qmd"), batch_effect_page)
    if (include_glm) {
        write_utf8(file.path(project_dir, "glm.qmd"), glm_page)
    }
    if (include_lme) {
        write_utf8(file.path(project_dir, "lme.qmd"), lme_page)
    }
    write_utf8(file.path(project_dir, "report.qmd"), report_page)
    write_utf8(file.path(project_dir, "logs.qmd"), logs_page)

    quarto_bin <- find_quarto()

    source_files <- file.path(project_dir, c(
        "_quarto.yml", "index.qmd",
        "enmix-qc.qmd", "metrics.qmd", "quality-control.qmd",
        "batch-effect.qmd",
        if (include_glm) "glm.qmd",
        if (include_lme) "lme.qmd",
        "report.qmd",
        "logs.qmd"
    ))

    emitLogMinfiEwasWater(
        c(paste(
            "Generated", length(source_files),
            "report source files in:", project_dir
        )),
        verbose = FALSE,
        log_path = log_path
    )

    render_status <- "skipped"
    rendered_file <- file.path(project_dir, "docs", "index.html")
    error_message <- NULL

    if (nzchar(quarto_bin)) {
        emitLogMinfiEwasWater(c("Rendering report site..."),
            verbose = FALSE, log_path = log_path
        )
        old_wd <- getwd()
        on.exit(setwd(old_wd), add = TRUE)
        setwd(project_dir)
        render_output <- system2(quarto_bin, "render",
            stdout = TRUE,
            stderr = TRUE
        )
        status <- attr(render_output, "status")
        if (is.null(status)) {
            status <- 0L
        }
        if (length(render_output)) {
            emitLogMinfiEwasWater(c("Quarto output:", render_output),
                verbose = FALSE, log_path = log_path
            )
        }
        if (identical(status, 0L)) {
            render_status <- "rendered"
            emitLogMinfiEwasWater(c("Render complete."),
                verbose = FALSE,
                log_path = log_path
            )
        } else {
            render_status <- "failed"
            error_message <- paste(
                "Quarto render failed with status",
                status
            )
            emitLogMinfiEwasWater(error_message,
                verbose = verbose,
                log_path = log_path
            )
        }
    } else {
        render_status <- "failed"
        error_message <- paste(
            "Quarto CLI not found. The project files were created but not rendered.",
            "Install Quarto or set QUARTO_BIN to the full path of the quarto executable."
        )
        emitLogMinfiEwasWater(c(error_message, project_dir),
            verbose = verbose, log_path = log_path
        )
    }
    remove_unrequested_rendered_assets()

    if (!magick_available && any(grepl("\\.tiff?$", c(
        enmix_items$original_name,
        metrics_items$original_name, qc_items$original_name,
        batch_items$original_name
    ), ignore.case = TRUE))) {
        emitLogMinfiEwasWater(
            c(
                "Note: TIFF figures were copied as-is because the `magick` package is not installed.",
                "Install `magick` and rerun this function if you want automatic PNG conversion for browser previews."
            ),
            verbose = FALSE, log_path = log_path
        )
    }

    emitLogMinfiEwasWater(paste("Report path:",
        normalizePathDnamReport(file.path(
        project_dir,
        "docs", "index.html"
    ))), verbose = verbose, log_path = log_path)

    render_result <- structure(list(
        preparedReport = prepared_report,
        status = render_status, renderedFile = if (identical(
            render_status,
            "rendered"
        )) {
            rendered_file
        } else {
            NULL
        }, errorMessage = error_message,
        logFile = log_path
    ), class = "dnaEPICO_dnamReport_render")

    structure(list(
        preparedReport = prepared_report, renderResult = render_result,
        status = render_status, outputFile = if (identical(
            render_status,
            "rendered"
        )) {
            rendered_file
        } else {
            normalizePathDnamReport(file.path(
                project_dir, "docs",
                "index.html"
            ))
        }, projectDir = project_dir, sourceFiles = source_files,
        resultTableSources = c(
            GLM = if (!include_glm) {
                "not_requested"
            } else if (is.null(glm_table_assets$source_mode)) {
                "unavailable"
            } else {
                glm_table_assets$source_mode
            },
            LME = if (!include_lme) {
                "not_requested"
            } else if (is.null(lme_table_assets$source_mode)) {
                "unavailable"
            } else {
                lme_table_assets$source_mode
            }
        ),
        modelSections = model_sections,
        docsDir = file.path(project_dir, "docs"), logoPath = logo_source_path,
        errorMessage = error_message, logFile = log_path
    ), class = "dnaEPICO_dnamReport")
}

#' Print a DNA methylation report result
#'
#' @param x Object returned by [dnamReport()].
#' @param ... Additional arguments ignored.
#'
#' @return Invisibly returns `x`.
#'
#' @export
print.dnaEPICO_dnamReport <- function(x, ...) {
    report_path <- x$outputFile
    if (is.null(report_path) || !length(report_path) ||
        is.na(report_path[[1]]) ||
        !nzchar(report_path[[1]])) {
        report_path <- x$renderResult$renderedFile
    }

    cat("Class type: ", class(x)[[1]], "\n", sep = "")
    cat("Log output path: ", formatPrintPathDnamReport(x$logFile),
        "\n",
        sep = ""
    )
    cat("Report output path: ", formatPrintPathDnamReport(report_path),
        "\n",
        sep = ""
    )

    invisible(x)
}
