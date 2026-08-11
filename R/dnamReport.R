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
    invalid_text <- paste(invalid, collapse = ", ")
    allowed_text <- paste(allowed, collapse = ", ")
    stop(
        sprintf(
        "Unsupported modelSections value(s): %s. Expected any of: %s.",
        invalid_text, allowed_text
        ),
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
#' @param patterns Character vector of regular expressions used to find
#'   matching files.
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
#'   outputDir = file.path(report_root, "reports"),
#'   qcDir = file.path(
#'     report_root,
#'     "figures",
#'     "preprocessingMinfiEwasWater",
#'     "enmix"
#'   ),
#'   preprocessingDir = file.path(
#'     report_root,
#'     "figures",
#'     "preprocessingMinfiEwasWater",
#'     "qc"
#'   ),
#'   postprocessingDir = file.path(
#'     report_root,
#'     "figures",
#'     "preprocessingMinfiEwasWater",
#'     "metrics"
#'   ),
#'   svaDir = file.path(report_root, "figures", "svaEnmix")
#' )
#' inherits(prepared, "dnaEPICO_dnamReport_prepared")
#'
#' @export
prepareDnamReportInputs <- function(outputDir = "reports",
    qcDir = file.path("figures", "preprocessingMinfiEwasWater",
        "enmix"), preprocessingDir = file.path("figures",
        "preprocessingMinfiEwasWater", "qc"),
    postprocessingDir = file.path("figures",
        "preprocessingMinfiEwasWater", "metrics"),
    svaDir = file.path("figures", "svaEnmix"),
    glmDir = file.path("figures", "methylationGLM"),
    lmeDir = file.path("figures", "methylationLME"),
    figDir = file.path(outputDir, "assets", "figures"),
    verbose = FALSE, logs = FALSE, logDir = outputDir) {
    log_path <- resolveLogPathMinfiEwasWater(logs = logs,
        log_dir = logDir, log_file = "log_dnamReport.txt")
    figure_inventory <- list(qc = collectFigureInventoryDnamReport(directory =
        qcDir, patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$",
            "\\.svg$", "\\.tif$", "\\.tiff$"),
        label = "ENmix QC"), preprocessing = collectFigureInventoryDnamReport(
            directory = preprocessingDir,
        patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$",
            "\\.svg$", "\\.tif$", "\\.tiff$"),
        label = "Quality control"), postprocessing =
            collectFigureInventoryDnamReport(directory = postprocessingDir,
        patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$",
            "\\.svg$", "\\.tif$", "\\.tiff$"),
        label = "Postprocessing"), sva = collectFigureInventoryDnamReport(
            directory = svaDir, patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$",
            "\\.svg$", "\\.tif$", "\\.tiff$"),
        label = "SVA"), glm = collectFigureInventoryDnamReport(directory =
            glmDir, patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$",
            "\\.svg$", "\\.tif$", "\\.tiff$"),
        label = "GLM"), lme = collectFigureInventoryDnamReport(directory =
            lmeDir, patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$",
            "\\.svg$", "\\.tif$", "\\.tiff$"), label = "LME"))
    missing_directories <- vapply(figure_inventory,
        function(section) !isTRUE(section$exists), logical(1))
    missing_directories <- names(missing_directories)[missing_directories]
    emitLogMinfiEwasWater(c(
        "============================================================",
        "Prepared DNA methylation dashboard report inputs",
        paste("Output directory:", normalizePathDnamReport(outputDir)),
        paste("Figure directory:", normalizePathDnamReport(figDir)),
        "============================================================"),
        verbose = verbose, log_path = log_path)
    output_file <- normalizePathDnamReport(file.path(outputDir,
        "docs", "index.html")); structure(list(output = basename(output_file),
        outputDir = normalizePathDnamReport(outputDir),
        outputFile = output_file, figDir = normalizePathDnamReport(figDir),
        figureInventory = figure_inventory, missingFigureDirectories =
            missing_directories,
        logFile = log_path), class = "dnaEPICO_dnamReport_prepared") }

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
#'   outputDir = file.path(report_root, "reports")
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
    stop(
        "preparedReport must be an object returned by ",
        "prepareDnamReportInputs().",
        call. = FALSE
    )
    }

    log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = if (is.null(logDir)) {
        dirname(preparedReport$outputFile)
    } else {
        logDir
    }, log_file = "log_dnamReport.txt"
    )

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
    errorMessage = paste(
        "Use dnamReport() to generate the Quarto",
        "dashboard report."
    ),
    logFile = log_path
    ), class = "dnaEPICO_dnamReport_render")
}

.isAbsolutePathDnamReport <- function(path) {
    grepl("^[A-Za-z]:[/\\\\]|^/", path)
}

.inferWorkflowRootDnamReport <- function(outputDir) {
    normalizedOutput <- if (.isAbsolutePathDnamReport(outputDir)) {
    normalizePath(outputDir, winslash = "/", mustWork = FALSE)
    } else {
    normalizePath(
        file.path(getwd(), outputDir),
        winslash = "/",
        mustWork = FALSE
    )
    }

    if (identical(basename(dirname(normalizedOutput)), "reports")) {
    return(dirname(dirname(normalizedOutput)))
    }
    if (identical(basename(normalizedOutput), "reports")) {
    return(dirname(normalizedOutput))
    }

    dirname(normalizedOutput)
}

.slashDnamReport <- function(path) {
    gsub("\\\\", "/", path)
}

.writeUtf8DnamReport <- function(path, lines) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    connection <- file(path, open = "w", encoding = "UTF-8")
    on.exit(close(connection), add = TRUE)
    writeLines(lines, con = connection, useBytes = TRUE)
}

.findQuartoDnamReport <- function() {
    environmentCandidates <- c(Sys.getenv("QUARTO_BIN"), Sys.getenv("QUARTO"))
    environmentCandidates <- environmentCandidates[
    nzchar(environmentCandidates)
    ]
    environmentCandidates <- normalizePath(
    environmentCandidates,
    winslash = "/",
    mustWork = FALSE
    )
    environmentCandidates <- environmentCandidates[
    file.exists(environmentCandidates)
    ]
    if (length(environmentCandidates)) {
    return(environmentCandidates[[1]])
    }

    quartoBinary <- Sys.which("quarto")
    if (nzchar(quartoBinary)) {
    return(unname(quartoBinary))
    }

    candidates <- c(
    "C:/Program Files/RStudio/resources/app/bin/quarto/bin/quarto.exe",
    "C:/Program Files/Quarto/bin/quarto.exe",
    file.path(
        Sys.getenv("LOCALAPPDATA"),
        "Programs", "Quarto", "bin", "quarto.exe"
    ),
    file.path(
        Sys.getenv("LOCALAPPDATA"),
        "quarto", "bin", "quarto.exe"
    )
    )
    candidates <- normalizePath(
    candidates,
    winslash = "/",
    mustWork = FALSE
    )
    candidates <- candidates[file.exists(candidates)]
    if (length(candidates)) {
    return(candidates[[1]])
    }

    ""
}

.htmlEscapeDnamReport <- function(text) {
    text <- gsub("&", "&amp;", text, fixed = TRUE)
    text <- gsub("<", "&lt;", text, fixed = TRUE)
    text <- gsub(">", "&gt;", text, fixed = TRUE)
    text <- gsub("\"", "&quot;", text, fixed = TRUE)
    text
}

.rStringDnamReport <- function(text) {
    text <- gsub("\\\\", "\\\\\\\\", text)
    text <- gsub("\"", "\\\"", text, fixed = TRUE)
    paste0("\"", text, "\"")
}

.slugifyDnamReport <- function(text) {
    text <- tools::file_path_sans_ext(basename(text))
    text <- tolower(text)
    text <- gsub("[^a-z0-9]+", "-", text)
    text <- gsub("(^-+|-+$)", "", text)
    if (!nzchar(text)) {
    text <- "figure"
    }
    text
}

.prettyLabelDnamReport <- function(path) {
    label <- tools::file_path_sans_ext(basename(path))
    label <- gsub("([a-z0-9])([A-Z])", "\\1 \\2", label, perl = TRUE)
    label <- gsub("[_\\-]+", " ", label)
    label <- gsub("\\s+", " ", label)
    label <- trimws(label)
    words <- strsplit(tolower(label), "\\s+", perl = TRUE)[[1]]
    label <- paste0(
    toupper(substr(words, 1L, 1L)),
    substr(words, 2L, nchar(words))
    )
    label <- paste(label, collapse = " ")
    label <- gsub("\\bQc\\b", "QC", label)
    label <- gsub("\\bPca\\b", "PCA", label)
    label <- gsub("\\bSva\\b", "SVA", label)
    label <- gsub("\\bEnmix\\b", "ENmix", label)
    label
}

.calloutLinesDnamReport <- function(text, type = "note") {
    c(sprintf("::: {.callout-%s}", type), text, ":::", "")
}

.dr_resolve_report_path <- function(path, base = root_dir) {
    if (is.null(path) || !nzchar(path)) {
    return("")
    }
    if (.isAbsolutePathDnamReport(path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
    }
    normalizePath(file.path(base, path),
    winslash = "/",
    mustWork = FALSE
    )
}

.dr_infer_model_from_path <- function(path,
                                                    step =
            "preprocessingMinfiEwasWater") {
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

.dr_copy_log_asset <- function(src_path, output_name) {
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

.dr_copy_first_existing_log_asset <- function(src_paths, output_name) {
    existing <- src_paths[file.exists(src_paths)]
    src_path <- if (length(existing) > 0L) {
    existing[[1L]]
    } else {
    src_paths[[1L]]
    }
    copy_log_asset(src_path, output_name)
}

.emptyReportFigureItemsDnaEpico <- function() {
    data.frame(
    title = character(), original_name = character(),
    asset_path = character(), download_path = character(),
    download_name = character(), browser_ready = logical(),
    converted = logical(), stringsAsFactors = FALSE
    )
}

.copyReportFigureFileDnaEpico <- function(
    srcPath, baseSlug, destinationDir, magickAvailable
) {
    extension <- tolower(tools::file_ext(srcPath))
    converted <- FALSE
    browser_ready <- TRUE
    if (extension %in% c("tif", "tiff") && magickAvailable) {
    name <- paste0(baseSlug, ".png")
    path <- file.path(destinationDir, name)
    converted <- tryCatch(
        {
        image <- magick::image_read(srcPath)
        magick::image_write(image, path = path, format = "png")
        TRUE
        },
        error = function(error) FALSE
    )
    if (converted) {
        download <- file.path(
        destinationDir, paste0(baseSlug, "-original.", extension)
        )
    } else {
        name <- paste0(baseSlug, ".", extension)
        path <- file.path(destinationDir, name)
        browser_ready <- FALSE
        download <- path
    }
    } else {
    name <- paste0(baseSlug, ".", extension)
    path <- file.path(destinationDir, name)
    browser_ready <- !(extension %in% c("tif", "tiff"))
    download <- path
    }
    copy_source <- if (converted) srcPath else srcPath
    copy_target <- if (converted) download else path
    if (!file.copy(copy_source, copy_target, overwrite = TRUE)) {
    message_type <- if (converted) "original figure" else "figure"
    stop("Failed to copy ", message_type, ": ", srcPath)
    }
    list(
    name = name, path = path, download = download,
    browserReady = browser_ready, converted = converted
    )
}

.reportFigureItemDnaEpico <- function(
    srcPath, index, assetSubdir, assetsFiguresDir, magickAvailable
) {
    base_slug <- sprintf("%02d-%s", index, .slugifyDnamReport(srcPath))
    copied <- .copyReportFigureFileDnaEpico(
    srcPath, base_slug, file.path(assetsFiguresDir, assetSubdir),
    magickAvailable
    )
    data.frame(
    title = .prettyLabelDnamReport(srcPath),
    original_name = basename(srcPath),
    asset_path = .slashDnamReport(file.path(
        "assets", "figures", assetSubdir, copied$name
    )),
    download_path = .slashDnamReport(file.path(
        "assets", "figures", assetSubdir, basename(copied$download)
    )),
    download_name = basename(srcPath),
    browser_ready = copied$browserReady, converted = copied$converted,
    stringsAsFactors = FALSE
    )
}

.dr_copy_figure_assets <- function(src_dir, asset_subdir) {
    destination_dir <- file.path(assets_figures_dir, asset_subdir)
    dir.create(destination_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(src_dir)) {
    return(.emptyReportFigureItemsDnaEpico())
    }
    files <- sort(list.files(
    src_dir,
    pattern = image_pattern, full.names = TRUE,
    recursive = recursive, ignore.case = TRUE
    ))
    if (!length(files)) {
    return(.emptyReportFigureItemsDnaEpico())
    }
    rows <- lapply(seq_along(files), function(index) {
    .reportFigureItemDnaEpico(
        files[[index]], index, asset_subdir, assets_figures_dir,
        magick_available
    )
    })
    do.call(rbind, rows)
}

.dr_js_quote_result_values <- function(values) {
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

.dr_js_result_array <- function(values) {
    paste0(
    "[", paste(js_quote_result_values(values), collapse = ","),
    "]"
    )
}

.dr_read_optional_workbook_sheet <- function(path, sheets, sheet) {
    if (!(sheet %in% sheets)) {
    return(data.frame())
    }
    tryCatch(openxlsx::read.xlsx(path, sheet = sheet, check.names = FALSE),
    error = function(e) data.frame()
    )
}

.dr_metadata_frame_to_list <- function(metadata) {
    if (!is.data.frame(metadata) || !all(c("Key", "Value") %in%
    names(metadata))) {
    return(list())
    }
    values <- as.character(metadata$Value)
    names(values) <- as.character(metadata$Key)
    base::as.list(values[nzchar(names(values))])
}

.dr_split_report_metadata_values <- function(value) {
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

.dr_extract_report_formula_phenotype <- function(formula) {
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

.dr_resolve_report_formula_records <- function(dictionary, metadata = list()) {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
    .drs_resolve_report_formula_records_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_last_log_field <- function(path, label) {
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

.dr_resolve_lme_report_metadata <- function(metadata,
    dictionary, log_path) {
    values <- metadata_frame_to_list(metadata)
    fallback <- list(libraries = last_log_field(log_path, "LME libraries"),
        correlation_structure = last_log_field(log_path,
            "Correlation structure"),
        correlation_variable = last_log_field(log_path, "Correlation variable"),
        interaction_term = last_log_field(log_path, "Interaction term"),
        annotation_columns = last_log_field(log_path,
            "Annotation columns used"),
        missing_annotation_columns = last_log_field(log_path,
            "Missing annotation columns"))
    for (key in names(fallback)) {
        if (is.null(values[[key]]) || !nzchar(values[[key]])) {
            values[[key]] <- fallback[[key]]
        } }
    if (is.null(values$backend) || !nzchar(values$backend)) {
        values$backend <- if (grepl("nlme", values$libraries,
            ignore.case = TRUE)) {
            "nlme"
        }
        else if (grepl("lme4|lmerTest", values$libraries, ignore.case = TRUE)) {
            "lme4"
        }
        else {
            ""
        } }
    if (is.null(values$fitting_function) || !nzchar(values$fitting_function)) {
        values$fitting_function <- if (identical(values$backend,
            "nlme")) {
            "nlme::lme"
        }
        else if (identical(values$backend, "lme4")) {
            "lmerTest::lmer"
        }
        else {
            ""
        } }
    if (is.null(values$correlation_structure) || !nzchar(
        values$correlation_structure)) {
        values$correlation_structure <- "none"
    }
    if (is.null(values$correlation_variable) || !nzchar(
        values$correlation_variable)) {
        values$correlation_variable <- "None"
    }
    if (is.null(values$interaction_term) || !nzchar(values$interaction_term)) {
        values$interaction_term <- "None"
    }; values
}

.dr_same_report_path <- function(source, target) {
    source <- normalizePath(source, winslash = "/", mustWork = FALSE)
    target <- normalizePath(target, winslash = "/", mustWork = FALSE)
    if (identical(.Platform$OS.type, "windows")) {
    identical(tolower(source), tolower(target))
    } else {
    identical(source, target)
    }
}

.dr_prepare_table_viewer_directory <- function(var_prefix) {
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

.dr_write_table_viewer_chunk <- function(chunk, chunk_number, result_dir,
    var_prefix, id_column) {
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

.dr_write_table_viewer_manifest <- function(columns,
    n_rows, chunk_size, chunk_index, var_prefix,
    id_column, downloads, item_singular,
    item_plural) {
    chunk_objects <- if (nrow(chunk_index)) {
        vapply(seq_len(nrow(chunk_index)),
            function(i) {
                paste0("{\"number\":", chunk_index$number[[i]],
                    ",\"firstId\":", js_quote_result_values(
            chunk_index$first_id[[i]]),
                    ",\"lastId\":", js_quote_result_values(
            chunk_index$last_id[[i]]),
                    "}")
            }, character(1))
    } else {
        character()
    }
    download_objects <- if (length(downloads)) {
        vapply(downloads, function(download) {
            paste0("{\"href\":", js_quote_result_values(download$href),
                ",\"label\":", js_quote_result_values(download$label),
                "}")
        }, character(1))
    }
    else {
        character()
    }
    manifest <- paste0(paste0("window.dnaEPICOResultManifests=",
        "window.dnaEPICOResultManifests||{};"),
        "window.dnaEPICOResultManifests[",
        js_quote_result_values(var_prefix),
        "]={", "\"key\":", js_quote_result_values(var_prefix),
        ",", "\"columns\":", js_result_array(columns),
        ",", "\"totalRows\":", n_rows, ",",
        "\"chunkSize\":", chunk_size, ",",
        "\"maxCachedChunks\":4,", "\"idColumn\":",
        js_quote_result_values(id_column),
        ",", "\"basePath\":", js_quote_result_values(paste0("assets/results/",
            var_prefix)), ",", "\"itemSingular\":",
        js_quote_result_values(item_singular),
        ",", "\"itemPlural\":", js_quote_result_values(item_plural),
        ",", "\"downloads\":[", paste(download_objects,
            collapse = ","), "],", "\"chunks\":[",
        paste(chunk_objects, collapse = ","),
        "]};")
    result_dir <- file.path(assets_dir, "results",
        var_prefix)
    write_utf8(file.path(result_dir, "manifest.js"),
        manifest)
}

.dr_table_viewer_asset_result <- function(columns, n_rows, var_prefix,
                                                        id_column,
            item_singular,
                                                        item_plural,
            search_label,
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

.dr_write_table_viewer_assets <- function(table_data,
    var_prefix, id_column, downloads, chunk_size = 5000L,
    item_singular = "CpG", item_plural = "CpGs",
    search_label = "Find CpG", search_placeholder = "e.g. cg00000029") {
    table_data <- as.data.frame(table_data,
        stringsAsFactors = FALSE, check.names = FALSE)
    chunk_size <- max(100L, as.integer(chunk_size))
    result_dir <- prepare_table_viewer_directory(var_prefix)
    n_rows <- nrow(table_data)
    n_chunks <- if (n_rows) {
        ceiling(n_rows/chunk_size)
    }
    else {
        0L
    }
    chunk_rows <- vector("list", n_chunks)
    if (n_chunks) {
        for (chunk_idx in seq_len(n_chunks)) {
            start <- (chunk_idx - 1L) * chunk_size +
                1L
            end <- min(chunk_idx * chunk_size,
                n_rows)
            chunk <- table_data[start:end,
                , drop = FALSE]
            chunk_rows[[chunk_idx]] <- write_table_viewer_chunk(chunk = chunk,
                chunk_number = chunk_idx,
                result_dir = result_dir,
                var_prefix = var_prefix,
                id_column = id_column)
        }
    }
    chunk_index <- if (length(chunk_rows)) {
        do.call(rbind, chunk_rows)
    }
    else {
        data.frame()
    }
    write_table_viewer_manifest(columns = names(table_data),
        n_rows = n_rows, chunk_size = chunk_size,
        chunk_index = chunk_index, var_prefix = var_prefix,
        id_column = id_column, downloads = downloads,
        item_singular = item_singular, item_plural = item_plural)
    table_viewer_asset_result(columns = names(table_data),
        n_rows = n_rows, var_prefix = var_prefix,
        id_column = id_column, item_singular = item_singular,
        item_plural = item_plural, search_label = search_label,
        search_placeholder = search_placeholder)
}

.dr_write_delimited_table_viewer_assets <- function(table_path,
    var_prefix, id_column, downloads, expected_rows,
    expected_columns, chunk_size = 5000L,
    item_singular = "CpG", item_plural = "CpGs",
    search_label = "Find CpG", search_placeholder = "e.g. cg00000029") {
    chunk_size <- max(100L, as.integer(chunk_size))
    result_dir <- prepare_table_viewer_directory(var_prefix)
    chunk_state <- new.env(parent = emptyenv())
    chunk_state$rows <- list()
    stream_result <- streamReportTableDnaEpico(tableFile = table_path,
        chunkSize = chunk_size, expectedRows = expected_rows,
        expectedColumns = expected_columns,
        chunkHandler = function(chunk, chunk_number) {
            chunk_state$rows[[chunk_number]] <- write_table_viewer_chunk(
            chunk = chunk,
                chunk_number = chunk_number,
                result_dir = result_dir,
                var_prefix = var_prefix,
                id_column = id_column)
        })
    chunk_rows <- chunk_state$rows
    columns <- stream_result$columns
    n_rows <- stream_result$nRows
    if (!(id_column %in% columns)) {
        stop("The report table identifier column was not found.",
            call. = FALSE)
    }
    chunk_index <- if (length(chunk_rows)) {
        do.call(rbind, chunk_rows)
    }
    else {
        data.frame()
    }
    write_table_viewer_manifest(columns = columns,
        n_rows = n_rows, chunk_size = chunk_size,
        chunk_index = chunk_index, var_prefix = var_prefix,
        id_column = id_column, downloads = downloads,
        item_singular = item_singular, item_plural = item_plural)
    table_viewer_asset_result(columns = columns,
        n_rows = n_rows, var_prefix = var_prefix,
        id_column = id_column, item_singular = item_singular,
        item_plural = item_plural, search_label = search_label,
        search_placeholder = search_placeholder)
}

.dr_prepare_xlsx_table_assets <- function(
    data_path, sheet, var_prefix,
    analysis = c("glm", "lme"), chunk_size = 5000L, model_log_path = "",
        additional_downloads = list()
) {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
    .drs_prepare_xlsx_table_assets_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.secondaryWorkbookSheetAssetsDnaEpico <- function(sourcePath, sheet, varPrefix,
    analysis, chunkSize, workbookDownload, assetsDir, slugifyFunction,
    viewerWriter) {
    sheet_data <- tryCatch(openxlsx::read.xlsx(sourcePath, sheet = sheet,
        check.names = FALSE), error = function(error) error)
    if (inherits(sheet_data, "error")) {
        return(list(ok = FALSE, error = conditionMessage(sheet_data), sheet =
            sheet,
            metadata = list()))
    }
    sheet_data <- as.data.frame(sheet_data, stringsAsFactors = FALSE,
        check.names = FALSE)
    sheet_key <- paste0(varPrefix, "_", slugifyFunction(sheet))
    sheet_dir <- file.path(assetsDir, "results", sheet_key)
    dir.create(sheet_dir, recursive = TRUE, showWarnings = FALSE)
    text_name <- paste0(safeFigureComponentDnaEpico(sheet), ".tsv.gz")
    text_path <- file.path(sheet_dir, text_name)
    data.table::fwrite(sheet_data, file = text_path, sep = "\t", quote = TRUE,
        na = "", compress = "gzip", showProgress = FALSE)
    id_column <- if (ncol(sheet_data)) {
        names(sheet_data)[[1L]]
    } else {
        "" }
    first_id <- if (nzchar(id_column) && nrow(sheet_data)) {
        as.character(sheet_data[[id_column]][[1L]])
    } else {
        "" }
    downloads <- c(list(list(href = paste0("assets/results/", sheet_key,
        "/", text_name), label = "Download selected sheet (TSV.gz)")),
        workbookDownload)
    viewer <- viewerWriter(table_data = sheet_data, var_prefix = sheet_key,
        id_column = id_column, downloads = downloads, chunk_size = chunkSize,
        item_singular = "row", item_plural = "rows", search_label = if (nzchar(
            id_column)) {
            paste("Find", id_column)
        }
        else {
            "Find row"
        }, search_placeholder = if (nzchar(first_id)) {
            paste("e.g.", first_id)
        }
        else {
            "Enter a value"
        })
    c(list(ok = TRUE, error = NULL, source_path = sourcePath, analysis =
        analysis,
        metadata = list(), sheet = sheet, source_mode = "xlsx_sheet",
            text_path = paste0("assets/results/",
            sheet_key, "/", text_name)), viewer)
}

.dr_prepare_xlsx_workbook_assets <- function(
    data_path, primary_sheet,
    var_prefix, analysis = c("glm", "lme"), chunk_size = 5000L,
        model_log_path = ""
) {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
    .drs_prepare_xlsx_workbook_assets_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_prepare_csv_table_assets <- function(
    data_path, var_prefix,
    front_columns = character(), chunk_size = 5000L
) {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
    .drs_prepare_csv_table_assets_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_content_description_html <- function(description = character(),
                                                        dynamic_role = NULL) {
    description <- as.character(description)
    description <- description[!is.na(description) &
    nzchar(trimws(description))]
    if (!length(description) &&
    (is.null(dynamic_role) || !nzchar(dynamic_role))) {
    return("")
    }
    description <- gsub("`([^`]+)`", "\\1", description, perl = TRUE)
    description <- trimws(gsub("\\s+", " ", description))
    dynamic_paragraph <- if (is.null(dynamic_role) ||
    !nzchar(dynamic_role)) {
    ""
    } else {
    sprintf(
        "<p data-role=\"%s\"></p>", html_escape(dynamic_role)
    )
    }
    paste0(
    "<div class=\"dnaepico-content-description\">",
    paste0("<p>", html_escape(description), "</p>",
        collapse = ""
    ),
    dynamic_paragraph,
    "</div>"
    )
}

.dr_content_note_block <- function(notes) {
    if (is.null(notes) || !length(notes)) {
    return(character(0))
    }
    notes <- as.character(notes)
    notes <- notes[!is.na(notes) & nzchar(trimws(notes))]
    if (!length(notes)) {
    return(character(0))
    }
    note_lines <- unlist(lapply(notes, function(note) {
    c(note, "")
    }), use.names = FALSE)
    c("::: {.dnaepico-content-notes}", "", note_lines, ":::", "")
}

.resultTablePaginationDnaEpico <- function(position) {
    c(
    sprintf(paste0(
        "  <div class=\"dnaepico-viewer-pagination dnaepico-",
        "viewer-pagination-%s\">"
    ), position),
    paste0(
        "    <button type=\"button\" data-role=\"first\" title=",
        "\"First page\" aria-label=\"First page\">&laquo;</button>"
    ),
    paste0(
        "    <button type=\"button\" data-role=\"previous\" ",
        "title=\"Previous page\" aria-label=\"Previous page\">",
        "&lsaquo;</button>"
    ),
    paste0(
        "    <label>Page <input data-role=\"page-number\" type=",
        "\"number\" min=\"1\" value=\"1\" /></label>"
    ),
    "    <span data-role=\"page-count\"></span>",
    paste0(
        "    <button type=\"button\" data-role=\"next\" title=",
        "\"Next page\" aria-label=\"Next page\">&rsaquo;</button>"
    ),
    paste0(
        "    <button type=\"button\" data-role=\"last\" title=",
        "\"Last page\" aria-label=\"Last page\">&raquo;</button>"
    ),
    "  </div>"
    )
}

.resultTableControlsDnaEpico <- function(searchLabel,
    searchPlaceholder, descriptionHtml) {
    c("```{=html}", "  <div class=\"dnaepico-viewer-controls\">",
        paste0("    <label><span>Rows per page</span><select data-role=",
            "\"page-size\"><option>25</option><option>50</option>",
            "<option>100</option></select></label>"),
        sprintf(paste0("    <label><span>%s</span><input data-role=",
            "\"cpg-search\" type=\"search\" placeholder=\"%s\" /></label>"),
            searchLabel, searchPlaceholder),
        "    <button type=\"button\" data-role=\"find-cpg\">Find</button>",
        paste0("    <label><span>Filter column</span><select data-role=",
            "\"filter-column\"><option value=\"\">Choose a column</option>",
            "</select></label>"), paste0(
            "    <label><span>Condition</span><select data-role=",
            "\"filter-operator\"><option value=\"contains\">contains</option>",
            "<option value=\"equals\">equals</option><option value=\"lt\">",
            "&lt;</option><option value=\"lte\">&le;</option><option value=",
            "\"gt\">&gt;</option><option value=\"gte\">&ge;</option></select>",
            "</label>"), paste0(
            "    <label><span>Value</span><input data-role=",
            "\"filter-value\" type=\"text\" ",
            "placeholder=\"e.g. 0.05\" /></label>"),
        paste0("    <button type=\"button\" data-role=\"apply-filter\">",
            "Apply filter</button>"), paste0(
            "    <button type=\"button\" data-role=\"clear-filter\" ",
            "disabled>Clear</button>"), "  </div>",
        paste0("  <span class=\"dnaepico-viewer-filter-summary\" data-role=",
            "\"filter-summary\" aria-live=\"polite\"></span>"),
        descriptionHtml, paste0(
            "  <div class=\"dnaepico-viewer-downloads\" data-role=",
            "\"downloads\"></div>"), "```",
        "")
}

.resultTableCardDnaEpico <- function(title, key) {
    c(
    sprintf(paste0(
        "::: {.card .dnaepico-result-card title=\"%s\" ",
        "expandable=\"true\" fill=\"false\"}"
    ), title), "",
    sprintf(paste0(
        "::: {.dnaepico-result-content data-result-key=",
        "\"%s\"}"
    ), key), "", "```{=html}",
    paste0(
        "  <div class=\"dnaepico-viewer-status\" data-role=",
        "\"status\" aria-live=\"polite\">Loading results&hellip;</div>"
    ),
    .resultTablePaginationDnaEpico("top"),
    "  <div class=\"dnaepico-viewer-table-wrap\">",
    paste0(
        "    <table class=\"table table-striped table-sm ",
        "dnaepico-viewer-table\">"
    ),
    "      <thead data-role=\"head\"></thead>",
    "      <tbody data-role=\"body\"></tbody>", "    </table>",
    "  </div>", .resultTablePaginationDnaEpico("bottom"),
    paste0(
        "  <noscript>JavaScript is required for paged browsing. ",
        "Use a complete-file download instead.</noscript>"
    ),
    "```", "", ":::", "", ":::", ""
    )
}

.dr_build_result_table_section <- function(
    title, table_assets, description = NULL, sheet = NULL, hidden = FALSE
) {
    if (!isTRUE(table_assets$ok)) {
    return(c(
        sprintf("::: {.card title=\"%s\"}", html_escape(title)), "",
        callout_lines(table_assets$error, type = "warning"), ":::", ""
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
    sheet_class <- if (is.null(sheet)) "" else " .dnaepico-workbook-sheet"
    sheet_attribute <- if (is.null(sheet)) {
    ""
    } else {
    sprintf(" data-sheet=\"%s\"", html_escape(sheet))
    }
    hidden_attribute <- if (isTRUE(hidden)) " hidden=\"hidden\"" else ""
    title <- html_escape(title)
    key <- html_escape(table_assets$key)
    c(
    sprintf(
        paste0(
        "::: {.dnaepico-cpg-viewer%s data-result-key=\"%s\"",
        "%s%s expandable=\"false\" fill=\"false\"}"
        ), sheet_class,
        key, sheet_attribute, hidden_attribute
    ), "",
    .resultTableControlsDnaEpico(
        html_escape(search_label), html_escape(search_placeholder),
        content_description_html(description)
    ),
    .resultTableCardDnaEpico(title, key),
    "```{=html}",
    sprintf("<script src=\"%s\"></script>", table_assets$manifest_path),
    "<script src=\"assets/cpg-viewer.js\"></script>",
    "```", "", ":::", ""
    )
}

.dr_build_data_frame_table_section <- function(title,
    data, empty_message, preview_rows = 25L,
    description = NULL) {
    if (!is.data.frame(data) || !nrow(data)) {
        return(c(sprintf("::: {.card title=\"%s\"}",
            html_escape(title)), "", callout_lines(empty_message,
            type = "warning"), ":::", ""))
    }
    if (is.finite(preview_rows) && nrow(data) >
        as.integer(preview_rows)) {
        data <- utils::head(data, as.integer(preview_rows))
    }
    header_cells <- paste0("    <th scope=\"col\">",
        html_escape(names(data)), "</th>")
    body_rows <- unlist(lapply(seq_len(nrow(data)),
        function(i) {
            row_values <- vapply(data[i,
                , drop = FALSE], function(value) {
                if (is.na(value)) {
                    ""
                }
                else {
                    as.character(value)
                }
            }, character(1))
            c("  <tr>", paste0("    <td>",
                html_escape(row_values),
                "</td>"), "  </tr>")
        }), use.names = FALSE)
    c(content_description_html(description),
        sprintf(paste0("::: {.card .dnaepico-summary-table-card title=\"%s\" ",
            "expandable=\"true\" fill=\"false\"}"),
            html_escape(title)), "", "```{=html}",
        sprintf(paste0("<div class=\"dnaepico-summary-table-wrap\" ",
            "tabindex=\"0\" role=\"region\" ",
            "aria-label=\"Scrollable %s\">"),
            html_escape(title)), paste0(
            "  <table class=\"table table-striped table-sm ",
            "dnaepico-summary-table\">"),
        "  <thead>", "  <tr>", header_cells,
        "  </tr>", "  </thead>", "  <tbody>",
        body_rows, "  </tbody>", "  </table>",
        "</div>", "```", "", ":::", "")
}

.dr_format_count <- function(value) {
    if (!length(value) || is.na(value)) {
    return("not available")
    }
    format(as.integer(round(value)),
    big.mark = ",", scientific = FALSE,
    trim = TRUE
    )
}

.dr_format_model_formula <- function(value) {
    value <- paste(as.character(value), collapse = "")
    value <- gsub("`", "", value, fixed = TRUE)
    value <- gsub("\\s+", " ", value, perl = TRUE)
    value <- gsub("\\s+\\)", ")", value, perl = TRUE)
    trimws(value)
}

.dr_build_model_formula_notes <- function(formula_records,
    model_label) {
    if (!is.data.frame(formula_records) ||
        !nrow(formula_records)) {
        return(character())
    }
    display_formulas <- vapply(formula_records$formula,
        format_model_formula, character(1))
    if (identical(tolower(model_label), "nlme")) {
        display_formulas <- sub("^LME:",
            "nlme:", display_formulas, ignore.case = TRUE)
    }
    if (nrow(formula_records) == 1L) {
        return(paste0(paste0(
            "The recorded model formula is <code class=\"dnaepico-",
            "model-formula\">"), html_escape(display_formulas[[1L]]),
            "</code>."))
    }
    formula_items <- vapply(seq_len(nrow(formula_records)),
        function(index) {
            phenotype <- formula_records$phenotype[[index]]
            if (is.na(phenotype) || !nzchar(trimws(phenotype))) {
                phenotype <- paste("Model",
                    index)
            }
            paste0("<div class=\"dnaepico-model-formula-item\">",
                "<strong>", html_escape(phenotype),
                "</strong>", "<code class=\"dnaepico-model-formula\">",
                html_escape(display_formulas[[index]]),
                "</code>", "</div>")
        }, character(1))
    c(sprintf("%d phenotype-specific models were fitted.",
        nrow(formula_records)), paste0(
            "<details class=\"dnaepico-model-formulas\">",
        "<summary>View recorded model formulas (",
        nrow(formula_records), ")</summary>",
        "<div class=\"dnaepico-model-formula-list\">",
        paste(formula_items, collapse = ""),
        "</div>", "</details>"))
}

.dr_format_decimal <- function(value, digits = 2L) {
    if (!length(value) || is.na(value)) {
    return("not available")
    }
    format(round(as.numeric(value), digits),
    nsmall = digits,
    scientific = FALSE, trim = TRUE
    )
}

.dr_plural <- function(n, singular, plural_form = paste0(
                                        singular,
                                        "s"
                                    )) {
    if (length(n) && !is.na(n) && as.integer(n) == 1L) {
    singular
    } else {
    plural_form
    }
}

.dr_is_numeric_like <- function(values) {
    values <- as.character(values)
    grepl("^\\s*[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?\\s*$",
    values,
    perl = TRUE
    )
}

.dr_sort_values <- function(values) {
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

.dr_collapse_values <- function(values) {
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

.dr_is_blank <- function(values) {
    is.na(values) | !nzchar(trimws(as.character(values)))
}

.dr_safe_read_table_file <- function(path, sheet = 1) {
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

.findDetectionPvalueObjectDnaEpico <- function(environment, objectNames) {
    if ("detP" %in% objectNames) {
    return(environment$detP)
    }
    matrix_names <- objectNames[vapply(objectNames, function(name) {
    object <- get(name, envir = environment)
    is.matrix(object) || is.data.frame(object)
    }, logical(1))]
    if (length(matrix_names)) {
    get(matrix_names[[1L]], envir = environment)
    } else {
    NULL
    }
}

.detectionTablesFromMatrixDnaEpico <- function(detP, threshold) {
    detP <- as.matrix(detP)
    storage.mode(detP) <- "numeric"
    if (is.null(colnames(detP))) {
    colnames(detP) <- sprintf("Sample_%s", seq_len(ncol(detP)))
    }
    assessed <- is.finite(detP) & detP >= 0 & detP <= 1
    detected <- assessed & detP < threshold
    assessed_cpg <- rowSums(assessed)
    detected_cpg <- rowSums(detected)
    cpg <- data.frame(
    metric = c(
        "Total CpGs assessed", "CpGs detected in at least one sample",
        "CpGs detected in all samples", "CpGs never detected",
        "CpGs with no valid detection P values"
    ),
    nCpGs = c(
        nrow(detP), sum(detected_cpg >= 1L),
        sum(assessed_cpg == ncol(detP) & detected_cpg == ncol(detP)),
        sum(assessed_cpg > 0L & detected_cpg == 0L),
        sum(assessed_cpg == 0L)
    ), stringsAsFactors = FALSE
    )
    assessed_sample <- colSums(assessed)
    detected_sample <- colSums(detected)
    percent <- rep(NA_real_, length(assessed_sample))
    available <- assessed_sample > 0L
    percent[available] <- 100 * detected_sample[available] /
    assessed_sample[available]
    sample <- data.frame(
    UID = colnames(detP), nAssessed = assessed_sample,
    nDetected = detected_sample, pDetected = round(percent, 2),
    stringsAsFactors = FALSE
    )
    list(cpg = cpg, sample = sample)
}

.loadDetectionTablesDnaEpico <- function(path, threshold) {
    tryCatch(
    {
        environment <- new.env(parent = emptyenv())
        object_names <- load(path, envir = environment)
        object <- .findDetectionPvalueObjectDnaEpico(
        environment, object_names
        )
        if (is.null(object)) {
        stop("No `detP` matrix-like object was found.", call. = FALSE)
        }
        .detectionTablesFromMatrixDnaEpico(object, threshold)
    },
    error = function(error) error
    )
}

.dr_read_detection_tables <- function(
    detp_path, threshold, cpg_path = "", sample_path = ""
) {
    result <- list(
    source = "detP", path = slash(detp_path), exists = FALSE,
    threshold = threshold, cpg = NULL, sample = NULL, error = NULL
    )
    if (file.exists(detp_path)) {
    loaded <- .loadDetectionTablesDnaEpico(detp_path, threshold)
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

.dr_detection_table_warning <- function(detection_tables) {
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

.dr_pick_column <- function(data, candidates) {
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

.dr_pick_participant_column <- function(data) {
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
    grepl(
        paste0(
        "sentrix|sample|array|slide|position|time|visit|plate",
        "|well|pool|probe"
        ),
        id_cols,
        ignore.case = TRUE
    )
    id_cols <- id_cols[!technical_id]
    if (length(id_cols)) {
    return(id_cols[[1]])
    }

    NULL
}

.dr_pick_timepoint_column <- function(data) {
    pick_column(data, c(
    "Timepoint", "Tiempoint", "Time_Point",
    "Timepoint_ID", "TimepointID", "Visit", "VisitID",
    "Visit_ID"
    ))
}

.dr_summarize_dataset <- function(data_path) {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
    .drs_summarize_dataset_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_summarize_cpg_detection <- function(data) {
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

.dr_summarize_sample_detection <- function(data) {
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

.dr_extract_sva_log_summary <- function(log_path) {
    summary <- list(k = NA_integer_, percent_variation = NA_real_)
    if (!file.exists(log_path)) {
    return(summary)
    }

    text <- paste(readLines(log_path, warn = FALSE, encoding = "UTF-8"),
    collapse = " "
    )
    pattern <- paste0(
    "([0-9]+)\\s+surrogate variables explain\\s+",
    "([0-9.]+)\\s*%"
    )
    match <- regexec(pattern,
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

.dr_summarize_sva <- function(pheno_path, log_path) {
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

.dr_normalize_log_field <- function(value) {
    value <- tolower(trimws(as.character(value)))
    trimws(gsub("[^[:alnum:]]+", " ", value))
}

.dr_parse_log_fields <- function(lines) {
    lines <- enc2utf8(as.character(lines))
    field_lines <- lines[grepl("^[[:space:]]*[^:]+:[[:space:]]*.*$",
    lines,
    perl = TRUE
    )]
    if (!length(field_lines)) {
    return(character())
    }
    keys <- normalize_log_field(sub(":.*$", "", field_lines))
    values <- trimws(sub("^[^:]+:[[:space:]]*", "", field_lines,
    perl = TRUE
    ))
    keep <- nzchar(keys) & nzchar(values)
    keys <- keys[keep]
    values <- values[keep]
    if (!length(keys)) {
    return(character())
    }
    grouped <- split(seq_along(keys), keys)
    vapply(grouped, function(index) {
    values[[utils::tail(index, 1L)]]
    }, character(1))
}

.dr_log_field <- function(fields, labels) {
    if (!length(fields)) {
    return("")
    }
    labels <- normalize_log_field(labels)
    for (label in labels) {
    if (label %in% names(fields) && nzchar(fields[[label]])) {
        return(fields[[label]])
    }
    }
    ""
}

.dr_log_capture <- function(lines, pattern) {
    if (!length(lines)) {
    return("")
    }
    for (line in rev(lines)) {
    match <- regexec(pattern, line, perl = TRUE, ignore.case = TRUE)
    values <- regmatches(line, match)[[1L]]
    if (length(values) >= 2L && nzchar(trimws(values[[2L]]))) {
        return(trimws(values[[2L]]))
    }
    }
    ""
}

.dr_summarize_logs <- function(log_assets) {
    labels <- c(methylation = "Methylation Analysis",
        data = "Data Preparation", batch = "Batch Effect",
        glm = "GLM", lme = "LME")
    stages <- lapply(names(log_assets), function(name) {
        asset <- log_assets[[name]]
        lines <- if (isTRUE(asset$exists)) {
            readLines(asset$source_path,
                warn = FALSE, encoding = "UTF-8")
        }
        else {
            character()
        }
        list(key = name, label = labels[[name]],
            exists = isTRUE(asset$exists),
            lines = lines, fields = parse_log_fields(lines),
            warning_count = sum(grepl("^[[:space:]]*warning(?:[[:space:]]|:)",
                lines, ignore.case = TRUE,
                perl = TRUE)), error_count = sum(grepl(
            "^[[:space:]]*(?:error|fatal)(?:[[:space:]]|:)",
                lines, ignore.case = TRUE,
                perl = TRUE)))
    })
    names(stages) <- names(log_assets)
    rows <- lapply(names(log_assets), function(name) {
        stage <- stages[[name]]
        line_count <- if (isTRUE(stage$exists)) {
            length(stage$lines)
        }
        else {
            NA_integer_
        }
        data.frame(key = name, label = labels[[name]],
            exists = isTRUE(stage$exists),
            source_path = log_assets[[name]]$source_path,
            line_count = line_count, warning_count = stage$warning_count,
            error_count = stage$error_count,
            stringsAsFactors = FALSE)
    })
    rows <- do.call(rbind, rows)
    list(rows = rows, total = nrow(rows),
        found = sum(rows$exists), total_lines = sum(rows$line_count,
            na.rm = TRUE), stages = stages)
}

.dr_make_data_notes <- function(summary) {
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
        paste0(
        "Using `%s` as the participant identifier and `%s` ",
        "as the timepoint column, the file contains %s ",
        "unique %s across timepoint values %s."
        ),
        summary$participant_col, summary$timepoint_col,
        format_count(summary$n_participants), plural(
        summary$n_participants,
        "participant"
        ), collapse_values(summary$timepoints)
    ))
    } else if (!is.null(summary$participant_col)) {
    notes <- c(notes, sprintf(
        paste0(
        "Using `%s` as the participant identifier, the file ",
        "contains %s unique %s. No Timepoint column was ",
        "detected."
        ),
        summary$participant_col, format_count(summary$n_participants),
        plural(summary$n_participants, "participant")
    ))
    } else {
    notes <- c(
        notes,
        paste0(
        "No biological participant identifier column was ",
        "detected. Technical identifiers such as sample, ",
        "Sentrix, array, and slide IDs are not treated as ",
        "participant IDs."
        )
    )
    }

    notes
}

.dr_make_metrics_notes <- function(items, figure_titles, figure_descriptions,
                                                data_summary) {
    if (!is.data.frame(items) || !nrow(items)) {
    return("No methylation metric figures were available.")
    }
    sprintf(
    paste0(
        "%s metric %s available for assessment of signal, ",
        "sample separation, and post-filtering methylation quality."
    ),
    format_count(nrow(items)), plural(nrow(items), "figure")
    )
}

.dr_make_quality_control_notes <- function(items, figure_titles,
                                                        figure_descriptions,
            cpg_summary,
                                                        sample_summary) {
    notes <- if (is.data.frame(items) && nrow(items)) {
    sprintf(
        "%s methylation quality-control %s available.",
        format_count(nrow(items)), plural(nrow(items), "figure")
    )
    } else {
    "No methylation quality-control figures were available."
    }

    if (isTRUE(cpg_summary$exists)) {
    notes <- c(notes, sprintf(
        paste0(
        "`Table 1` reports %s assessed CpGs, %s detected in ",
        "all samples, and %s never detected."
        ),
        format_count(cpg_summary$total_assessed),
        format_count(cpg_summary$detected_all),
        format_count(cpg_summary$never_detected)
    ))
    }

    if (isTRUE(sample_summary$exists)) {
    notes <- c(notes, sprintf(
        paste0(
        "`Table 2` shows sample detection percentages ",
        "ranging from %s%% to %s%% across %s %s."
        ),
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

.dr_make_batch_effect_notes <- function(items, figure_titles,
                                                        figure_descriptions,
            sva_summary) {
    notes <- if (is.data.frame(items) && nrow(items)) {
    sprintf(
        paste0(
        "%s batch-effect %s available for assessing latent ",
        "variation and technical structure."
        ),
        format_count(nrow(items)), plural(nrow(items), "figure")
    )
    } else {
    "No batch-effect figures were available."
    }

    if (!is.na(sva_summary$log_k) &&
    !is.na(sva_summary$percent_variation)) {
    notes <- c(notes, sprintf(
        paste0(
        "The SVA log reports %s surrogate %s explaining %s%% ",
        "of data variation."
        ),
        format_count(sva_summary$log_k), plural(
        sva_summary$log_k,
        "variable"
        ), format_decimal(sva_summary$percent_variation)
    ))
    }

    notes
}

.reportLogSampleCountDnaEpico <- function(lines, capture) {
    count <- capture(lines,
        "loaded with[[:space:]]+([0-9,]+)[[:space:]]+samples")
    if (!nzchar(count)) {
    count <- capture(
        lines,
        "using all[[:space:]]+([0-9,]+)[[:space:]]+samples"
    )
    }
    count
}

.methylationLogDescriptionDnaEpico <- function(
    fields, sampleCount, field,
    sentenceList
) {
    details <- c(
    if (nzchar(sampleCount)) paste(sampleCount, "samples"),
    if (nzchar(field(fields, "Array type"))) {
        paste(field(fields, "Array type"), "arrays")
    },
    if (nzchar(field(fields, "Normalization methods"))) {
        paste(field(fields, "Normalization methods"), "normalisation")
    },
    if (nzchar(field(fields, "Annotation version"))) {
        paste(field(fields, "Annotation version"), "annotation")
    },
    if (nzchar(field(fields, "Detection p-value threshold"))) {
        paste(
        "detection p-value threshold",
        field(fields, "Detection p-value threshold")
        )
    },
    if (nzchar(field(fields, "Reference"))) {
        paste(field(fields, "Reference"), "cell-composition reference")
    }
    )
    paste0(
    "IDAT loading, quality control, filtering, and cell-composition ",
    "estimation were recorded",
    if (length(details)) paste0(" with ", sentenceList(details)) else "",
    "."
    )
}

.dataLogDescriptionDnaEpico <- function(
    fields, sampleCount, field,
    sentenceList
) {
    identifier <- field(fields, c("Identifier column", "SampleID column"))
    timepoints <- field(fields, "Timepoints")
    dimensions <- field(fields, "Beta dimensions")
    details <- c(
    if (nzchar(sampleCount)) paste(sampleCount, "samples"),
    if (nzchar(identifier)) paste("identifier", sprintf("`%s`", identifier)),
    if (nzchar(timepoints)) paste("timepoints", timepoints),
    if (nzchar(dimensions)) paste("beta-matrix dimensions", dimensions)
    )
    paste0(
    "Phenotype preparation, timepoint handling, and methylation-matrix ",
    "creation were recorded",
    if (length(details)) paste0(" for ", sentenceList(details)) else "",
    "."
    )
}

.batchLogDescriptionDnaEpico <- function(
    fields, sampleCount, field,
    sentenceList, pluralFunction
) {
    count <- field(fields, c(
    "Number of surrogate variables (K)",
    "Number of surrogate variables"
    ))
    numeric_count <- coerceNumericDnaEpico(count)[[1L]]
    technical <- field(fields, "Technical terms modelled")
    details <- c(
    if (nzchar(sampleCount)) paste(sampleCount, "samples"),
    if (nzchar(count)) {
        paste(count, if (!is.na(numeric_count)) {
        pluralFunction(numeric_count, "surrogate variable")
        } else {
        "surrogate variables"
        })
    },
    if (nzchar(technical)) paste("technical terms", technical)
    )
    paste0(
    "Surrogate-variable estimation and technical-factor assessment ",
    "were recorded",
    if (length(details)) paste0(" with ", sentenceList(details)) else "",
    "."
    )
}

.glmLogDescriptionDnaEpico <- function(fields, field, sentenceList) {
    values <- list(
    phenotypes = field(fields, "Phenotypes"),
    covariates = field(fields, "Covariates"),
    interaction = field(fields, "Interaction term"),
    retained = field(fields, "CpG columns retained"),
    adjustment = field(fields, "P-value adjustment method")
    )
    details <- c(
    if (nzchar(values$phenotypes)) paste("phenotypes", values$phenotypes),
    if (nzchar(values$covariates)) paste("covariates", values$covariates),
    if (nzchar(values$interaction)) paste("interaction", values$interaction),
    if (nzchar(values$retained)) paste(values$retained, "CpGs tested"),
    if (nzchar(values$adjustment)) paste(values$adjustment, "adjustment")
    )
    paste0(
    "Model fitting, association testing, CpG annotation, and output ",
    "generation were recorded",
    if (length(details)) paste0(" for ", sentenceList(details)) else "",
    "."
    )
}

.lmeLogDescriptionDnaEpico <- function(fields, field, sentenceList) {
    participant <- field(fields, "Participant key")
    outputs <- field(fields, "Outputs")
    details <- c(
    if (nzchar(participant)) paste("participant key", participant),
    if (nzchar(outputs)) paste("outputs", outputs)
    )
    paste0(
    "Longitudinal model fitting, CpG annotation, and output generation ",
    "were recorded",
    if (length(details)) paste0(" with ", sentenceList(details)) else "",
    "."
    )
}

.reportLogStageDescriptionDnaEpico <- function(
    key, stage, field, capture, sentenceList, pluralFunction
) {
    sample_count <- .reportLogSampleCountDnaEpico(stage$lines, capture)
    switch(key,
    methylation = .methylationLogDescriptionDnaEpico(
        stage$fields, sample_count, field, sentenceList
    ),
    data = .dataLogDescriptionDnaEpico(
        stage$fields, sample_count, field, sentenceList
    ),
    batch = .batchLogDescriptionDnaEpico(
        stage$fields, sample_count, field, sentenceList, pluralFunction
    ),
    glm = .glmLogDescriptionDnaEpico(stage$fields, field, sentenceList),
    lme = .lmeLogDescriptionDnaEpico(stage$fields, field, sentenceList)
    )
}

.dr_make_logs_notes <- function(
    log_summary, lme_label = "LME Analysis"
) {
    labels <- c(
    methylation = "Methylation Analysis", data = "Data Preparation",
    batch = "Batch Effect", glm = "GLM Analysis", lme = lme_label
    )
    notes <- vapply(log_summary$rows$key, function(key) {
    stage <- log_summary$stages[[key]]
    if (!isTRUE(stage$exists)) {
        return(sprintf("`%s` was not found.", labels[[key]]))
    }
    description <- .reportLogStageDescriptionDnaEpico(
        key, stage, log_field, log_capture, sentence_list, plural
    )
    sprintf("`%s` is available. %s", labels[[key]], description)
    }, character(1))
    names(notes) <- log_summary$rows$key
    notes
}

.dr_report_number_markup <- function(text) {
    gsub(
    paste0(
        "(?<![[:alnum:]_.])",
        "([0-9]+(?:,[0-9]{3})*(?:[.][0-9]+)?%?)",
        "(?![[:alnum:]_.])"
    ),
    paste0(
        "<span class=\"dnaepico-data-value\">",
        "\\1</span>"
    ),
    text,
    perl = TRUE
    )
}

.dr_report_inline_markup <- function(text) {
    vapply(as.character(text), function(value) {
    parts <- strsplit(value, "`", fixed = TRUE)[[1L]]
    markup <- vapply(seq_along(parts), function(index) {
        escaped <- html_escape(parts[[index]])
        if (index %% 2L == 0L) {
        return(paste0(
            "<span class=\"dnaepico-data-value\">",
            escaped, "</span>"
        ))
        }
        report_number_markup(escaped)
    }, character(1))
    paste0(markup, collapse = "")
    }, character(1))
}

.dr_html_paragraph <- function(text, highlight_values = FALSE) {
    content <- if (isTRUE(highlight_values)) {
    report_inline_markup(text)
    } else {
    html_escape(text)
    }
    sprintf("<p>%s</p>", content)
}

.dr_html_bullet_list <- function(items, labelled = FALSE,
                                                class_name = "") {
    items <- as.character(items)
    items <- items[nzchar(trimws(items))]
    if (!length(items)) {
    return(character())
    }
    list_class <- if (nzchar(class_name)) {
    sprintf(" class=\"%s\"", html_escape(class_name))
    } else {
    ""
    }
    item_markup <- vapply(items, function(item) {
    if (isTRUE(labelled) && grepl(":", item, fixed = TRUE)) {
        label <- trimws(sub(":.*$", "", item))
        value <- trimws(sub("^[^:]+:[[:space:]]*", "", item,
        perl = TRUE
        ))
        return(sprintf(
        "<li><strong>%s:</strong> %s</li>",
        html_escape(label), report_inline_markup(value)
        ))
    }
    sprintf("<li>%s</li>", report_inline_markup(item))
    }, character(1))
    c(sprintf("<ul%s>", list_class), item_markup, "</ul>")
}

.dr_html_section <- function(title, paragraphs = character(),
                                            bullets = character(), labelled =
            FALSE) {
    c(
    sprintf(
        "<section class=\"dnaepico-report-section\"><h3>%s</h3>",
        html_escape(title)
    ),
    html_paragraph(paragraphs, highlight_values = TRUE),
    html_bullet_list(bullets, labelled = labelled),
    "</section>"
    )
}

.dr_prepare_report_text <- function(text) {
    text <- gsub(
    "<code(?:\\s[^>]*)?>(.*?)</code>",
    "`\\1`",
    text,
    perl = TRUE
    )
    text <- gsub("<[^>]+>", "", text, perl = TRUE)
    text <- gsub("&quot;", "\"", text, fixed = TRUE)
    text <- gsub("&lt;", "<", text, fixed = TRUE)
    text <- gsub("&gt;", ">", text, fixed = TRUE)
    text <- gsub("&amp;", "&", text, fixed = TRUE)
    text <- gsub("\\s+", " ", text)
    trimws(text)
}

.dr_sentence_case <- function(text) {
    if (!nzchar(text)) {
    return(text)
    }
    paste0(toupper(substr(text, 1L, 1L)), substr(
    text, 2L,
    nchar(text)
    ))
}

.dr_report_note_text <- function(notes, fallback = character()) {
    notes <- notes[!grepl("dnaepico-model-formulas", notes,
    fixed = TRUE
    )]
    notes <- prepare_report_text(notes)
    notes <- notes[nzchar(notes)]
    if (!length(notes)) {
    if (!length(fallback) || !nzchar(fallback)) {
        return(character())
    }
    return(sprintf("%s.", sentence_case(fallback)))
    }
    unique(notes)
}

.dr_build_report_page <- function(
    project_name, project_dir, data_notes,
    enmix_notes, quality_control_notes, batch_effect_notes, metrics_notes,
        glm_notes,
    lme_notes, logs_notes, glm_visualisation_notes = character(),
        lme_visualisation_notes = character(),
    summary_items = character(), preprocessing_notes = character(),
        overlap_notes = character(),
    observations = character(), lme_label = "LME Analysis", model_sections = c(
    "glm",
    "lme"
    )
) {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
    .drs_build_report_page_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_search_index_markup <- function(title,
    body_lines) {
    body_lines <- enc2utf8(as.character(body_lines))
    keep <- logical(length(body_lines))
    in_script <- FALSE
    in_executable_chunk <- FALSE
    for (index in seq_along(body_lines)) {
        opens_executable_chunk <- grepl(
            "^```\\{(?:r|python|julia|ojs)(?:[ ,}])",
            trimws(body_lines[[index]]), perl = TRUE,
            ignore.case = TRUE)
        closes_chunk <- identical(trimws(body_lines[[index]]),
            "```")
        opens_script <- grepl("<script(?:\\s|>)", body_lines[[index]],
            perl = TRUE, ignore.case = TRUE)
        closes_script <- grepl("</script>", tolower(body_lines[[index]]),
            fixed = TRUE)
        if (opens_executable_chunk) {
            in_executable_chunk <- TRUE
        }
        if (opens_script) {
            in_script <- TRUE
        }
        keep[[index]] <- !in_script && !in_executable_chunk
        if (closes_script) {
            in_script <- FALSE
        }
        if (closes_chunk && in_executable_chunk) {
            in_executable_chunk <- FALSE
        } }
    text <- paste(body_lines[keep], collapse = " ")
    text <- gsub("<[^>]+>", " ", text, perl = TRUE)
    text <- gsub("&nbsp;", " ", text, fixed = TRUE)
    text <- gsub("&amp;", " and ", text, fixed = TRUE)
    text <- gsub("&lt;", " less than ", text, fixed = TRUE)
    text <- gsub("&gt;", " greater than ", text, fixed = TRUE)
    text <- gsub("[`*_#{}|]", " ", text, perl = TRUE)
    text <- gsub("\\s+", " ", text, perl = TRUE)
    text <- trimws(text)
    if (!nzchar(text)) {
        text <- title
    }
    if (nchar(text, type = "chars") > 50000L) {
        text <- substr(text, 1L, 50000L)
    }
    c("```{=html}", sprintf(paste0("<main class=\"dnaepico-search-index ",
        "quarto-include-in-search-index\" hidden>",
        "<h1>%s</h1><p>%s</p></main>"), html_escape(title),
        html_escape(text)), "```", "")
}

.dr_compose_page <- function(title, body_lines, body_classes = NULL) {
    front_matter <- c(
    "---", sprintf("title: \"%s\"", title),
    "format:", "  dashboard:", "    scrolling: true",
    "    orientation: columns"
    )

    if (!is.null(body_classes) && nzchar(body_classes)) {
    front_matter <- c(front_matter, sprintf(
        "body-classes: %s",
        body_classes
    ))
    }

    c(
    front_matter, "---", "", "## Column", "",
    search_index_markup(title, body_lines), body_lines, "",
    "```{=html}",
    "<script src=\"assets/report-interactions.js\"></script>",
    "```", ""
    )
}

.dr_build_workbook_table_section <- function(
    title, workbook_assets,
    model_notes = character()
) {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
    .drs_build_workbook_table_section_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_subset_figure_items <- function(items, pattern) {
    if (!is.data.frame(items) || !nrow(items) ||
    !("original_name" %in% names(items))) {
    return(data.frame())
    }
    items[grepl(pattern, items$original_name, ignore.case = TRUE), ,
    drop = FALSE
    ]
}

.dr_sentence_list <- function(values) {
    values <- trimws(as.character(values))
    values <- values[nzchar(values)]
    if (!length(values)) {
    return("")
    }
    if (length(values) == 1L) {
    return(values[[1L]])
    }
    if (length(values) == 2L) {
    return(paste(values, collapse = " and "))
    }
    paste0(
    paste(utils::head(values, -1L), collapse = ", "),
    ", and ", utils::tail(values, 1L)
    )
}

.dr_pretty_model_term <- function(value) {
    value <- gsub("Treatmentgroup", "Treatment group", value,
    ignore.case = TRUE
    )
    value <- gsub("([a-z0-9])([A-Z])", "\\1 \\2", value,
    perl = TRUE
    )
    value <- gsub("[_-]+", " ", value)
    trimws(gsub("\\s+", " ", value))
}

.dr_venn_effect_label <- function(value) {
    parts <- strsplit(value, "x", fixed = TRUE)[[1L]]
    parts <- parts[nzchar(parts)]
    if (length(parts) >= 2L) {
    parts <- vapply(parts, pretty_model_term, character(1))
    return(sprintf(
        "the interaction between %s",
        sentence_list(parts)
    ))
    }

    term <- pretty_model_term(value)
    if (identical(tolower(term), "time")) {
    return("the Time main effect")
    }
    sprintf("the %s effect", term)
}

.dr_combine_venn_interaction_tokens <- function(tokens) {
    tokens <- tokens[nzchar(tokens)]
    if (!length(tokens)) {
    return(character())
    }
    effects <- character()
    index <- 1L
    while (index <= length(tokens)) {
    if (index + 2L <= length(tokens) &&
        identical(tolower(tokens[[index + 1L]]), "x")) {
        effects <- c(effects, paste0(
        tokens[[index]], "x", tokens[[index + 2L]]
        ))
        index <- index + 3L
    } else {
        effects <- c(effects, tokens[[index]])
        index <- index + 1L
    }
    }
    effects
}

.dr_venn_figure_metadata <- function(filename) {
    .runDrStages(base::as.list.environment(environment(), all.names =
        TRUE),
    .drs_venn_figure_metadata_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_model_figure_title <- function(filename, analysis) {
    stem <- tools::file_path_sans_ext(basename(filename))
    cleaned <- gsub("[_.-]+", " ", stem)
    cleaned <- gsub("\\bGLM\\b|\\bLME\\b", "", cleaned,
    ignore.case = TRUE
    )
    cleaned <- trimws(gsub("\\s+", " ", cleaned))
    if (grepl("^manhattan", stem, ignore.case = TRUE)) {
    version <- if (grepl("_v1$", stem, ignore.case = TRUE)) {
        "Circular"
    } else {
        "Rectangular"
    }
    test <- sub("^manhattan_", "", stem, ignore.case = TRUE)
    test <- sub("_v[12]$", "", test, ignore.case = TRUE)
    return(paste(version, "Manhattan Plot for", gsub("_", " ", test)))
    }
    if (grepl("^vennD", stem, ignore.case = TRUE)) {
    return(venn_figure_metadata(filename)$title)
    }
    if (grepl("^intersection", stem, ignore.case = TRUE)) {
    cleaned <- gsub("([[:alnum:]])x([[:upper:]])", "\\1 by \\2",
        cleaned,
        perl = TRUE
    )
    return(paste("Ranked CpG and Gene Intersections:", cleaned))
    }
    if (grepl("^qqplot", stem, ignore.case = TRUE)) {
    test <- sub("^qqplot_", "", stem, ignore.case = TRUE)
    test <- sub("_coefficientPvalue$", "", test,
        ignore.case = TRUE
    )
    return(paste(analysis, "EWAS QQ Plot for", gsub("_", " ", test)))
    }
    if (grepl("^(hist|bar|distribution)_", stem, ignore.case = TRUE)) {
    variable <- sub(
        "^(hist|bar|distribution)_", "", stem,
        ignore.case = TRUE
    )
    variable <- sub("_(continuous|categorical)$", "", variable,
        ignore.case = TRUE
    )
    return(paste("Distribution of", gsub("_", " ", variable)))
    }
    pretty_label(filename)
}

.dr_model_figure_description <- function(filename, title,
    analysis) {
    .runDrStages(base::as.list.environment(environment(), all.names =
        TRUE),
    .drs_model_figure_description_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_extract_report_count <- function(value) {
    value <- as.character(value)
    if (!length(value) || is.na(value[[1L]]) || !nzchar(value[[1L]])) {
    return(NA_real_)
    }
    matched <- regmatches(value[[1L]], regexpr(
    "[0-9][0-9,]*", value[[1L]],
    perl = TRUE
    ))
    if (!length(matched) || !nzchar(matched)) {
    return(NA_real_)
    }
    as.numeric(gsub(",", "", matched, fixed = TRUE))
}

.dr_stage_sample_count <- function(stage) {
    if (is.null(stage) || !isTRUE(stage$exists)) {
    return(NA_real_)
    }
    field_value <- log_field(stage$fields, c(
    "Samples after filtering", "Number of samples",
    "Sample count", "nSamples"
    ))
    count <- extract_report_count(field_value)
    if (!is.na(count)) {
    return(count)
    }
    patterns <- c(
    "loaded with[[:space:]]+([0-9,]+)[[:space:]]+samples",
    "using all[[:space:]]+([0-9,]+)[[:space:]]+samples",
    "RGSet loaded with[[:space:]]+([0-9,]+)[[:space:]]+samples"
    )
    for (pattern in patterns) {
    count <- extract_report_count(log_capture(stage$lines, pattern))
    if (!is.na(count)) {
        return(count)
    }
    }
    NA_real_
}

.dr_make_report_summary_items <- function(
    data_summary, cpg_summary,
    sample_summary, sva_summary, log_summary, model_sections, figure_items,
        workbook_available
) {
    .runDrStages(base::as.list.environment(environment(), all.names =
        TRUE),
    .drs_make_report_summary_items_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_make_report_preprocessing_notes <- function(log_summary,
    logs_notes) {
    notes <- paste0(
    "Preprocessing converts raw array signals into comparable ",
    "methylation values, removes measurements that do not meet the ",
    "configured quality rules, and prepares the data for modelling."
    )
    for (key in intersect(c("methylation", "data"), names(
    log_summary$stages
    ))) {
    if (isTRUE(log_summary$stages[[key]]$exists) &&
        key %in% names(logs_notes)) {
        notes <- c(notes, logs_notes[[key]])
    }
    }
    unique(notes)
}

.dr_make_report_overlap_notes <- function(glm_items, lme_items) {
    collect_notes <- function(items, analysis) {
    if (!is.data.frame(items) || !nrow(items)) {
        return(character())
    }
    vapply(seq_len(nrow(items)), function(index) {
        filename <- items$original_name[[index]]
        title <- model_figure_title(filename, analysis)
        description <- model_figure_description(
        filename, title, analysis
        )
        sprintf("%s: %s", title, description)
    }, character(1))
    }
    unique(c(
    collect_notes(glm_items, "GLM"),
    collect_notes(lme_items, "LME")
    ))
}

.dr_make_report_observations <- function(
    data_summary, cpg_summary,
    sample_summary, log_summary
) {
    .runDrStages(base::as.list.environment(environment(), all.names =
        TRUE),
    .drs_make_report_observations_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_build_figure_browser <- function(
    items, empty_message, browser_id,
    analysis, title_prefix = "Figure", figure_titles = NULL,
        figure_descriptions = NULL,
    browser_notes = character()
) {
    .runDrStages(base::as.list.environment(environment(), all.names =
        TRUE),
    .drs_build_figure_browser_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.dr_build_figure_sections <- function(sections, empty_message,
    browser_prefix, titles = list(),
                                                    descriptions = list(),
            notes = list()) {
    lines <- c("::: {.panel-tabset}", "")
    for (index in seq_along(sections)) {
    section <- names(sections)[[index]]
    items <- sections[[index]]
    section_titles <- titles[[section]]
    section_descriptions <- descriptions[[section]]
    section_notes <- notes[[section]]
    lines <- c(
        lines, paste0("## ", section), "",
        build_figure_browser(
        items = items, empty_message = empty_message,
        browser_id = paste0(browser_prefix, "-", index),
        analysis = section, figure_titles = section_titles,
        figure_descriptions = section_descriptions,
        browser_notes = section_notes
        )
    )
    }
    c(lines, ":::", "")
}

.dr_build_model_visualisation_tabs <- function(items,
    venn, analysis) {
    prefix <- tolower(analysis)
    sections <- list(`Model Variables` = subset_figure_items(items,
        "^(hist_|bar_|distribution_)"), `Model Design` = subset_figure_items(
            items,
        paste0("^(missingness_|correlation_|modelVariables_|observat",
            "ions_|participantObservationCount_|timeDistribution_",
            "|timepointDistribution_|association_)")),
        `EWAS QQ Plots` = subset_figure_items(items,
            "^qqplot_"), `Model Diagnostics` = subset_figure_items(items,
            "^(residual|standardError|scaleLocation|influence)"),
        `Effect Summaries` = subset_figure_items(items,
            "^(volcano_|effectForest_)"))
    if (identical(analysis, "LME")) {
        sections[["Methylation Trajectories"]] <- subset_figure_items(items,
            paste0("^(trajectory_|change_|longitudinalTrajectory_|paired",
                "Change_)"))
    }
    sections[["Manhattan Plots"]] <- subset_figure_items(items,
        "^manhattan_")
    if (is.data.frame(venn) && nrow(venn)) {
        sections[["Venn Diagrams"]] <- venn
    }
    sections <- sections[vapply(sections,
        function(section_items) {
            is.data.frame(section_items) &&
                nrow(section_items) > 0L
        }, logical(1))]
    if (!length(sections)) {
        return(callout_lines(paste("No",
            analysis, "visualisation files were available."),
            type = "note"))
    }
    lines <- c("::: {.panel-tabset}", "")
    for (section_index in seq_along(sections)) {
        section <- names(sections)[[section_index]]
        section_items <- sections[[section_index]]
        lines <- c(lines, paste0("## ", section),
            "", build_figure_browser(section_items,
                paste("No", tolower(section),
                    "were available."), browser_id = paste0("dnaepico-",
                    prefix, "-figures-", section_index),
                analysis = analysis))
    }
    c(lines, ":::", "")
}

.dr_unrequested_table_assets <- function(analysis) {
    list(
    ok = FALSE, error = "The model section was not requested.",
    source_path = "", analysis = analysis, metadata = list(),
    source_mode = NULL
    )
}

.dr_describe_available_annotations <- function(columns) {
    notes <- character()
    if (all(c("IlmnID", "Name") %in% columns)) {
    notes <- c(notes, paste0(
        "The `IlmnID` and `Name` columns identify each ",
        "Illumina CpG probe."
    ))
    } else if ("IlmnID" %in% columns) {
    notes <- c(notes, "The `IlmnID` column identifies each Illumina CpG probe.")
    }
    if (all(c("chr", "pos") %in% columns)) {
    notes <- c(notes, paste0(
        "The `chr` and `pos` columns give the probe's ",
        "genomic position."
    )) }
    if (any(c("UCSC_RefGene_Group", "UCSC_RefGene_Name") %in%
    columns)) {
    notes <- c(notes, paste0(
        "Available `UCSC_RefGene` columns provide gene-",
        "related annotation."
    ))
    }
    if ("Relation_to_Island" %in% columns) {
    notes <- c(notes,
        "`Relation_to_Island` reports the probe's CpG-island context.")
    }
    if ("GencodeV41_Group" %in% columns) {
    notes <- c(notes, paste0(
        "`GencodeV41_Group` provides GENCODE version 41 gene-",
        "region annotation."
    ))
    }
    gencode_releases <- unique(sub(
    "^GencodeV([0-9]+)_.*$",
    "\\1", grep("^GencodeV[0-9]+_RefGene_Name$",
        columns,
        value = TRUE
    )
    ))
    if (length(gencode_releases)) {
    notes <- c(notes, sprintf(paste0(
        "GENCODE release %s direct gene-body and nearest-TSS ",
        "annotations are available in the release-labelled ",
        "columns."
    ), paste(gencode_releases,
        collapse = ", "
    )))
    }
    notes
}

.dr_log_status_note <- function(label) {
    logs_notes[grepl(paste0("`", label, "`"), logs_notes, fixed = TRUE)]
}

.dr_remove_unrequested_rendered_assets <- function() {
    unrequested_files <- c(if (!include_glm) {
        c(file.path(project_dir, "docs",
            "glm.html"), file.path(project_dir,
            "docs", "glm-visualisations.html"),
            file.path(project_dir, "docs",
                "assets", "logs", "methylationGLM.txt"))
    }, if (!include_lme) {
        c(file.path(project_dir, "docs",
            "lme.html"), file.path(project_dir,
            "docs", "lme-visualisations.html"),
            file.path(project_dir, "docs",
                "assets", "logs", "methylationLME.txt"))
    })
    unrequested_directories <- c(if (!include_glm) {
        c(file.path(project_dir, "docs",
            "assets", "results", "glm_results"),
            file.path(project_dir, "docs",
                "assets", "figures", "glm-visualisations"))
    }, if (!include_lme) {
        c(file.path(project_dir, "docs",
            "assets", "results", "lme_results"),
            file.path(project_dir, "docs",
                "assets", "figures", "lme-visualisations"))
    })
    unlink(unrequested_files, force = TRUE)
    unlink(unrequested_directories, recursive = TRUE,
        force = TRUE)
    invisible(NULL)
}

.drx0072332_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    overview <- paste0("This report summarises the study data, methylation ",
        "preprocessing, quality-control assessment, batch-effect ",
            "evaluation, and statistical models.")
    data_paragraphs <- report_note_text(data_notes, paste0(
        "displays the phenotype data preview and detected ",
        "participant/timepoint fields"))
    preprocessing_paragraphs <- report_note_text(preprocessing_notes, paste0(
        "Methylation preprocessing details were not available ",
        "when this report was generated"))
    quality_control_paragraphs <- report_note_text(c(paste0(
        "Quality control checks whether the methylation ",
        "measurements are technically reliable before the ",
            "biological results are interpreted."),
        enmix_notes, quality_control_notes), paste0(
            "displays ENmix control plots, methylation quality-",
        "control figures, and detection summary tables"))
    batch_effect_paragraphs <- report_note_text(c(paste0(
        "Batch effects are differences associated with ",
        "processing or other technical factors rather than the ",
            "biological questions being studied."),
        batch_effect_notes), "displays SVA and batch-effect figures")
    metrics_paragraphs <- report_note_text(c(paste0(
        "Beta values describe the proportion of methylation on ",
        "a scale from 0 to 1, while M-values are transformed ",
            "values commonly used for statistical modelling."),
        metrics_notes), "displays post-filtering methylation metric figures")
    glm_paragraphs <- report_note_text(c(paste0(
        "A generalised linear model tests whether each CpG is ",
        "associated with the selected phenotype while accounting ",
            "for the configured covariates."),
        glm_notes),
            "displays the annotated generalised linear model results table")
    environment()
}

.drx0072332_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    lme_paragraphs <- report_note_text(c(paste0(
        "A linear mixed-effects model tests CpG associations while ",
        "accounting for repeated observations from the same ", "participant."),
        lme_notes),
            "displays the annotated linear mixed-effects model results table")
    glm_visualisation_bullets <- report_note_text(glm_visualisation_notes,
        "displays GLM model, diagnostic, and genomic visualisations")
    lme_visualisation_bullets <- report_note_text(lme_visualisation_notes,
        "displays LME longitudinal, diagnostic, and genomic visualisations")
    logs_bullets <- report_note_text(logs_notes,
        "displays workflow log files for each analysis stage")
    observation_bullets <- report_note_text(observations,
        "No recorded workflow result requires specific attention")
    model_report_sections <- c(if ("glm" %in% model_sections) {
        html_section("Generalised linear model analysis", glm_paragraphs,
            bullets = glm_visualisation_bullets)
    }, if ("lme" %in% model_sections) {
        html_section("Linear mixed-effects analysis", lme_paragraphs, bullets =
            lme_visualisation_bullets)
    })
    report_sections <- c(html_section("Summary", bullets = summary_items,
        labelled = TRUE),
        html_section("Data and study design", data_paragraphs), html_section(
            "Methylation preprocessing",
            preprocessing_paragraphs), html_section("Quality control",
            quality_control_paragraphs),
        html_section("Methylation metrics", metrics_paragraphs), html_section(
            "Batch-effect assessment",
            batch_effect_paragraphs), model_report_sections, if (length(
            overlap_notes)) {
            html_section("CpG and gene-overlap figures", bullets =
            report_note_text(overlap_notes))
        }, html_section("Logs", paste0(
            "The full technical records remain available in the ",
            "Logs tab."), bullets = logs_bullets), html_section("Observations",
            bullets = observation_bullets))
    environment()
}

.drx0072332_part_03DnaEpico <- function() {
    .installDrHelpers(environment())
    .dnamReportFunctionResult <- environment()
    environment()
}

.drx0072332_stagesDnaEpico <- c(
    .drx0072332_part_01DnaEpico,
    .drx0072332_part_02DnaEpico,
    .drx0072332_part_03DnaEpico
)

.drs_build_report_page_part_01DnaEpico <- function() {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
        .drx0072332_stagesDnaEpico,
        parent = environment(sys.function()))
}

.drs_build_report_page_part_02DnaEpico <-
    function() {
    .installDrHelpers(environment())
    .dnamReportFunctionResult <- c(
    "---", "title: \"Report\"", "format:", "  dashboard:",
    "    scrolling: true", "    orientation: columns",
        "body-classes: qpasst-report-page",
    "---", "", "## Column", "", search_index_markup("Report", c(
        overview,
        summary_items, data_paragraphs, preprocessing_paragraphs,
            quality_control_paragraphs,
        batch_effect_paragraphs, metrics_paragraphs, glm_paragraphs,
            glm_visualisation_bullets,
        lme_paragraphs, lme_visualisation_bullets, overlap_notes, logs_bullets,
        observation_bullets
    )), "```{=html}", "<div class=\"qpasst-report-viewport\">",
    "<div class=\"qpasst-report-card\">", "<div class=\"qpasst-report-sheet\">",
    "<h2>Analysis Report</h2>", html_paragraph(overview), report_sections,
    "</div>", "</div>", "</div>", "```", "", "```{=html}",
        "<script src=\"assets/report-interactions.js\"></script>",
    "```", ""
    )
    environment()
}

.drs_build_report_page_stagesDnaEpico <- c(
    .drs_build_report_page_part_01DnaEpico,
    .drs_build_report_page_part_02DnaEpico
)

.drs_make_report_summary_items_part_01DnaEpico <-
    function() {
    .installDrHelpers(environment())
    items <- character()
    if (isTRUE(data_summary$exists)) {
        data_parts <- sprintf("%s records across %s variables", format_count(
            data_summary$n_rows),
            format_count(data_summary$n_cols))
        if (!is.na(data_summary$n_participants)) {
            data_parts <- paste0(data_parts, ", representing ", format_count(
            data_summary$n_participants),
                " ", plural(data_summary$n_participants, "participant"))
        }
        items <- c(items, sprintf("Study data: %s.", data_parts))
    }
    if (!is.na(data_summary$n_timepoints) && data_summary$n_timepoints > 0L) {
        design <- if (data_summary$n_timepoints > 1L) {
            "Longitudinal"
        }
        else {
            "Single-timepoint"
        }
        items <- c(items, sprintf("Study design: %s data covering %s %s (%s).",
            design, format_count(data_summary$n_timepoints), plural(
            data_summary$n_timepoints,
                "timepoint"), collapse_values(data_summary$timepoints)))
    }
    methylation_stage <- log_summary$stages$methylation
    assay <- if (!is.null(methylation_stage)) {
        log_field(methylation_stage$fields, "Array type")
    } else {
        "" }
    annotation <- if (!is.null(methylation_stage)) {
        log_field(methylation_stage$fields, "Annotation version")
    } else {
        "" }
    if (nzchar(assay) || nzchar(annotation)) {
        assay_details <- sentence_list(c(if (nzchar(assay)) assay, if (nzchar(
            annotation)) {
            paste(annotation, "annotation")
        }))
        items <- c(items, sprintf("Methylation assay: %s.", assay_details))
    }
    normalization <- if (!is.null(methylation_stage)) {
        log_field(methylation_stage$fields, "Normalization methods")
    } else {
        ""
    }
    environment()
}

.drs_make_report_summary_items_part_02DnaEpico <-
    function() {
    .runDrStages(base::as.list.environment(environment(), all.names =
        TRUE),
    .drs2_make_report_summary_items_part_02DnaEpico_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.drs_make_report_summary_items_part_03DnaEpico <-
    function() {
    .installDrHelpers(environment())
    figure_count <- sum(vapply(figure_items, function(items) {
    if (is.data.frame(items)) nrow(items) else 0L
    }, integer(1)))
    workbook_count <- sum(as.logical(workbook_available), na.rm = TRUE)
    output_parts <- c(sprintf("%s %s", format_count(figure_count), plural(
    figure_count,
    "figure"
    )), sprintf("%s result %s", format_count(workbook_count), plural(
    workbook_count,
    "workbook"
    )))
    items <- c(items, sprintf("Available report content: %s.", sentence_list(
        output_parts)))
    .dnamReportFunctionResult <- unique(items)
    environment()
}

.drs_make_report_summary_items_stagesDnaEpico <- c(
    .drs_make_report_summary_items_part_01DnaEpico,
    .drs_make_report_summary_items_part_02DnaEpico,
    .drs_make_report_summary_items_part_03DnaEpico
)

.drs_build_figure_browser_part_01DnaEpico <-
    function() {
    .runDrStages(base::as.list.environment(environment(), all.names =
        TRUE),
    .drs2_build_figure_browser_part_01DnaEpico_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.drs_build_figure_browser_part_02DnaEpico <-
    function() { .installDrHelpers(environment())
    .dnamReportFunctionResult <- c("```{=html}", sprintf(paste0(
            "<div id=\"%s-controls\" class=\"dnaepico-figure-browser\" ",
            "tabindex=\"0\" role=\"region\" ",
            "aria-label=\"Interactive figure viewer controls\">"),
            browser_id), "<div class=\"dnaepico-figure-toolbar\">",
        paste0("<button type=\"button\" data-role=\"figure-prev\" aria-",
            "label=\"Previous figure\">Previous</button>"),
        paste0("<label><span>Figure</span><select data-role=",
            "\"figure-select\">"), options, "</select></label>", paste0(
            "<button type=\"button\" data-role=\"figure-next\" aria-",
            "label=\"Next figure\">Next</button>"),
        sprintf(paste0("<span class=\"dnaepico-figure-count\" data-role=",
            "\"figure-count\" aria-live=\"polite\">",
            "1 of %d figures</span>"), nrow(items)),
        "</div>", content_description_html(browser_notes,
            dynamic_role = "figure-description"),
        "<div class=\"dnaepico-figure-actions\">",
        paste0("<a class=\"btn btn-sm btn-outline-primary\" data-role=",
            "\"figure-download\">Download original figure</a>"),
        "</div>", paste0("<script type=\"application/json\" data-role=",
            "\"figure-data\">", figure_data, "</script>"), "</div>", "```",
        "", paste0("::: {.card .dnaepico-selected-figure ",
            "title=\"Selected figure\" expandable=\"true\" ",
            "fill=\"false\"}"), "", sprintf(paste0(
            "::: {.dnaepico-figure-content ",
            "data-browser-id=\"%s\"}"), browser_id),
        "", "```{=html}", "<h3 data-role=\"figure-title\"></h3>",
        sprintf(paste0("<p id=\"%s-zoom-help\" class=",
            "\"dnaepico-figure-zoom-help\">",
            "Use Control or Command plus the mouse wheel to zoom. ",
            "Use the normal mouse wheel to scroll vertically, ",
            "Shift plus the mouse wheel to scroll horizontally, ",
            "and the plus, minus, or zero key to adjust or reset zoom.",
            "</p>"), browser_id), paste0(
            "<div class=\"dnaepico-figure-canvas\" data-role=",
            "\"figure-canvas\" tabindex=\"0\" role=\"region\" ",
            "aria-label=\"Interactive figure\" aria-describedby=\"",
            browser_id, "-zoom-help\">"),
        paste0("<div class=\"dnaepico-figure-stage\" data-role=",
            "\"figure-stage\">"), paste0(
            "<img data-role=\"figure-image\" alt=\"\" loading=\"lazy\" ",
            "draggable=\"false\" />"), "</div>", "</div>", paste0(
            "<span class=\"dnaepico-figure-zoom-status\" data-role=",
            "\"figure-zoom-status\" aria-live=\"polite\"></span>"),
        "<p data-role=\"figure-fallback\" hidden></p>",
        "```", "", ":::", "", ":::", "",
        "```{=html}", "<script src=\"assets/figure-viewer.js\"></script>",
        "```", ""); environment() }

.drs_build_figure_browser_stagesDnaEpico <- c(
    .drs_build_figure_browser_part_01DnaEpico,
    .drs_build_figure_browser_part_02DnaEpico
)

.drs_prepare_xlsx_table_assets_part_01DnaEpico <-
    function() {
    .installDrHelpers(environment())
    analysis <- match.arg(analysis, c("glm", "lme"))
    source_data_path <- if (grepl("^[A-Za-z]:[/\\\\]|^/", data_path)) {
        data_path } else {
        file.path(root_dir, data_path)
    }
    if (!file.exists(source_data_path)) {
        {
            .dnamReportReturnValue <- list(ok = FALSE, error = paste0(
            "Data file not found at `",
                slash(source_data_path), "`."), source_path = source_data_path,
                analysis = analysis, metadata = list())
            .dnamReportDidReturn <- TRUE
            return(environment())
        } }
    result_dir <- file.path(assets_dir, "results", var_prefix)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    xlsx_name <- basename(source_data_path)
    legacy_xlsx_target <- file.path(result_dir, xlsx_name)
    if (!same_report_path(source_data_path, legacy_xlsx_target) && file.exists(
        legacy_xlsx_target)) {
        unlink(legacy_xlsx_target, force = TRUE)
    }
    sidecar <- resolveReportTableSidecarDnaEpico(workbookFile =
        source_data_path,
        sidecarDir = result_dir, sheet = sheet)
    if (isTRUE(sidecar$ok)) {
        dictionary <- tryCatch(utils::read.delim(sidecar$dictionary,
            check.names = FALSE,
            stringsAsFactors = FALSE, quote = "\"", comment.char = ""), error =
            function(error) error)
        workbook_metadata <- if (file.exists(sidecar$workbookMetadata)) {
            tryCatch(utils::read.delim(sidecar$workbookMetadata, check.names =
            FALSE,
                stringsAsFactors = FALSE, quote = "\"", comment.char = ""),
            error = function(error) error)
        } else { data.frame() }
        if (inherits(dictionary, "error") || inherits(workbook_metadata,
            "error")) {
            sidecar$ok <- FALSE } }
    if (!isTRUE(sidecar$ok)) {
        sheets <- tryCatch(openxlsx::getSheetNames(source_data_path), error =
            function(e) character())
        dictionary <- read_optional_workbook_sheet(source_data_path, sheets,
            "dictionary")
        workbook_metadata <- read_optional_workbook_sheet(source_data_path,
            sheets, "metadata")
    }; environment() }

.drs_prepare_xlsx_table_assets_part_02DnaEpico <-
    function() { .installDrHelpers(environment())
    metadata <- if (identical(analysis, "lme")) {
        resolve_lme_report_metadata(workbook_metadata, dictionary,
            model_log_path) } else {
        metadata_frame_to_list(workbook_metadata) }
    metadata$formula_records <- resolve_report_formula_records(dictionary,
        metadata)
    metadata$formulas <- unique(metadata$formula_records$formula)
    text_name <- paste0(sheet, ".tsv.gz")
    text_target <- file.path(result_dir, text_name)
    downloads <- c(list(list(href = paste0("assets/results/", var_prefix, "/",
        text_name), label = "Download all results (TSV.gz)")),
            additional_downloads); viewer_assets <- NULL
    source_mode <- "xlsx_fallback"
    if (isTRUE(sidecar$ok)) {
        sidecar_staged <- same_report_path(sidecar$table, text_target) ||
            isTRUE(file.copy(sidecar$table,
            text_target, overwrite = TRUE))
        if (sidecar_staged && file.exists(text_target)) {
            viewer_assets <- tryCatch(write_delimited_table_viewer_assets(
            table_path = text_target,
                var_prefix = var_prefix, id_column = sidecar$idColumn,
            downloads = downloads,
                expected_rows = sidecar$rows, expected_columns =
            sidecar$columns,
                chunk_size = chunk_size), error = function(error) NULL)
            if (!is.null(viewer_assets)) {
                source_mode <- "streamed_sidecar"
            } } }; if (is.null(viewer_assets)) {
        table_data <- tryCatch(openxlsx::read.xlsx(source_data_path, sheet =
            sheet,
            check.names = FALSE), error = function(e) e)
        if (inherits(table_data, "error")) { {
                .dnamReportReturnValue <- list(ok = FALSE, error =
            conditionMessage(table_data),
                    source_path = source_data_path, analysis = analysis,
            metadata = list())
                .dnamReportDidReturn <- TRUE
                return(environment()) } }
        report_table <- sortReportTableDnaEpico(table_data)
        table_data <- report_table$data
        data.table::fwrite(table_data, file = text_target, sep = "\t", quote =
            TRUE,
            na = "", compress = "gzip", showProgress = FALSE)
        viewer_assets <- write_table_viewer_assets(table_data = table_data,
            var_prefix = var_prefix,
            id_column = report_table$idColumn, downloads = downloads,
            chunk_size = chunk_size)
    }; environment() }

.drs_prepare_xlsx_table_assets_part_03DnaEpico <-
    function() {
    .installDrHelpers(environment())
    .dnamReportFunctionResult <- c(list(
    ok = TRUE, error = NULL, source_path = source_data_path,
    analysis = analysis, dictionary = dictionary, metadata = metadata,
        source_mode = source_mode,
    text_path = paste0("assets/results/", var_prefix, "/", text_name)
    ), viewer_assets)
    environment()
}

.drs_prepare_xlsx_table_assets_stagesDnaEpico <- c(
    .drs_prepare_xlsx_table_assets_part_01DnaEpico,
    .drs_prepare_xlsx_table_assets_part_02DnaEpico,
    .drs_prepare_xlsx_table_assets_part_03DnaEpico
)

.drx0103252_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    observations <- character()
    if (!is.na(data_summary$repeated_participant_timepoint_rows) &&
        data_summary$repeated_participant_timepoint_rows >
        0L) {
        observations <- c(observations, sprintf(paste0(
            "%s repeated participant-timepoint rows were detected ",
            "in the Data tab and should be considered when ",
            "interpreting repeated-measures analyses."),
            format_count(data_summary$repeated_participant_timepoint_rows)))
    }
    if (!is.na(data_summary$n_timepoints) && data_summary$n_timepoints > 1L &&
        !is.na(data_summary$complete_timepoint_participants) && !is.na(
            data_summary$n_participants) &&
        data_summary$complete_timepoint_participants <
            data_summary$n_participants) {
        observations <- c(observations, sprintf(paste0(
            "%s of %s participants have records at every observed ",
            "timepoint; the remaining participants have incomplete ",
            "timepoint coverage."),
            format_count(data_summary$complete_timepoint_participants),
            format_count(data_summary$n_participants)))
    }
    stage_counts <- vapply(log_summary$stages, stage_sample_count, numeric(1))
    available_stage_counts <- stage_counts[!is.na(stage_counts)]
    environment()
}

.drx0103252_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    if (length(unique(available_stage_counts)) > 1L) {
        count_details <- vapply(names(available_stage_counts), function(key) {
            sprintf("%s: %s samples", log_summary$stages[[key]]$label,
            format_count(available_stage_counts[[key]]))
        }, character(1))
        observations <- c(observations, sprintf(paste0(
            "Workflow logs record different sample counts (%s). ",
            "These values were not combined."), paste(count_details, collapse =
            "; ")))
    }
    representative_log_count <- if (length(available_stage_counts)) {
        names(sort(table(available_stage_counts), decreasing = TRUE))[[1L]]
    }
    else {
        NA_character_
    }
    representative_log_count <- as.numeric(representative_log_count)
    if (isTRUE(data_summary$exists) && !is.na(representative_log_count) &&
        !is.na(data_summary$n_rows) &&
        data_summary$n_rows != representative_log_count) {
        observations <- c(observations, sprintf(paste0(
            "The Data tab contains %s records, while the workflow ",
            "logs most commonly record %s samples. These values were ",
            "kept separate because the available files do not ",
            "demonstrate that they describe the same report population."),
            format_count(data_summary$n_rows),
            format_count(representative_log_count)))
    }
    environment()
}

.drx0103252_part_03DnaEpico <- function() {
    .installDrHelpers(environment())
    if (isTRUE(sample_summary$exists) && !is.na(sample_summary$n_samples) &&
        !is.na(representative_log_count) && sample_summary$n_samples !=
            representative_log_count) {
        observations <- c(observations, sprintf(paste0(
            "Quality-control Table 2 summarises %s samples, while ",
            "the workflow logs most commonly record %s. The report ",
            "does not assume that these sources are equivalent."),
            format_count(sample_summary$n_samples), format_count(
            representative_log_count)))
    }
    .dnamReportFunctionResult <- environment()
    environment()
}

.drx0103252_stagesDnaEpico <- c(
    .drx0103252_part_01DnaEpico,
    .drx0103252_part_02DnaEpico,
    .drx0103252_part_03DnaEpico
)

.drs_make_report_observations_part_01DnaEpico <- function() {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
        .drx0103252_stagesDnaEpico,
        parent = environment(sys.function()))
}

.drs_make_report_observations_part_02DnaEpico <-
    function() {
    .installDrHelpers(environment())
    if (isTRUE(cpg_summary$exists) && !is.na(cpg_summary$never_detected) &&
        cpg_summary$never_detected >
    0L) {
    observations <- c(observations, sprintf(paste0(
        "Quality-control Table 1 records %s CpGs that were never ",
        "detected in the assessed samples."
    ), format_count(cpg_summary$never_detected)))
    }
    missing_stages <- vapply(log_summary$stages, function(stage) {
    !isTRUE(stage$exists)
    }, logical(1))
    if (any(missing_stages)) {
    missing_labels <- vapply(
        log_summary$stages[missing_stages], function(stage) stage$label,
        character(1)
    )
    observations <- c(observations, sprintf(
        "No log file was available for %s.",
        sentence_list(missing_labels)
    )) }
    issue_stages <- vapply(log_summary$stages, function(stage) {
    stage$warning_count > 0L || stage$error_count > 0L
    }, logical(1))
    if (any(issue_stages)) {
    for (stage in log_summary$stages[issue_stages]) {
        issue_parts <- c(if (stage$warning_count > 0L) {
        sprintf("%s warning %s", stage$warning_count, plural(
            stage$warning_count,
            "entry"
        ))
        }, if (stage$error_count > 0L) {
        sprintf("%s error %s", stage$error_count, plural(
            stage$error_count,
            "entry"
        )) })
        observations <- c(observations, sprintf(
        "The %s log contains %s.",
        stage$label, sentence_list(issue_parts)
        ))
    } }
    if (!length(observations)) {
    observations <- paste0(
        "No recorded workflow result requires specific attention.")
    }
    .dnamReportFunctionResult <- unique(observations)
    environment()
}

.drs_make_report_observations_stagesDnaEpico <- c(
    .drs_make_report_observations_part_01DnaEpico,
    .drs_make_report_observations_part_02DnaEpico
)

.drs_build_workbook_table_section_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    sheet_description <- function(sheet) {
        lower <- tolower(sheet)
        if (grepl("^annotated", lower)) {
            return(paste0("Association results and genomic annotation for the ",
                "selected model. Search, filter, or download the table ",
            "for detailed review.")) }
        if (grepl("metadata", lower)) {
            return(paste0(
            "Analysis settings, model provenance, and the recorded ",
                "meaning of any renamed Venn labels."))
        }; if (grepl("dictionary", lower)) {
            return(paste0(
            "Definitions and formulas for the columns supplied in ",
                "this annotated workbook.")) }
        if (grepl("venn|intersection", lower)) {
            return(paste0(
            "CpG or gene membership used to construct the requested ",
                "model-overlap visualisation.")) }
        paste0("Preview of the selected worksheet from the annotated model ",
            "workbook.") }
    if (is.null(workbook_assets$sheetAssets) || !length(
        workbook_assets$sheetAssets)) { {
            .dnamReportReturnValue <- build_result_table_section(title,
            workbook_assets, description = c(paste0(
            "Search, filter, and download the annotated model ",
                    "results from this table."), model_notes))
            .dnamReportDidReturn <- TRUE
            return(environment()) } }
    selector_id <- paste0("dnaepico-workbook-", safeFigureComponentDnaEpico(
        workbook_assets$analysis))
    options <- vapply(workbook_assets$sheets, function(sheet) {
        sprintf("<option value=\"%s\">%s</option>", html_escape(sheet),
            html_escape(sheet)) }, character(1))
    lines <- c(sprintf(paste0("::: {.card title=\"%s\" expandable=\"false\" ",
        "fill=\"false\"}"), html_escape(title)), "", content_note_block(c(
            paste0("Select a worksheet to inspect association results, ",
        "metadata, Venn membership tables, or column definitions."),
            model_notes)), "```{=html}", sprintf(
            "<div class=\"dnaepico-workbook-selector\" id=\"%s\">",
            selector_id), "<label>Sheet <select data-role=\"workbook-sheet\">",
        options, "</select></label>", if (!is.null(
            workbook_assets$workbookPath) &&
            nzchar(workbook_assets$workbookPath)) {
            sprintf(paste0("<a class=\"btn btn-sm btn-outline-primary\" ",
                "href=\"%s\">", "Download complete workbook (XLSX)</a>"),
            workbook_assets$workbookPath) } else { ""
        }, "</div>", "```", "", ":::", "")
    environment() }

.drs_build_workbook_table_section_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    for (sheet_index in seq_along(workbook_assets$sheets)) {
    sheet <- workbook_assets$sheets[[sheet_index]]
    lines <- c(lines, build_result_table_section(paste0("Sheet: ", sheet),
        workbook_assets$sheetAssets[[sheet]],
        description = sheet_description(sheet),
        sheet = sheet, hidden = sheet_index != 1L
    ))
    }
    lines <- c(
    lines, "```{=html}", "<script>", "(function(){", sprintf(
        "var root=document.getElementById('%s');",
        selector_id
    ), "if(!root){return;}",
        "var select=root.querySelector('[data-role=workbook-sheet]');",
    paste0("var panels=Array.from(document.querySelectorAll(",
        "'.dnaepico-workbook-sheet'));"),
    "function activate(){", "  panels.forEach(function(panel){",
        "    var active=panel.getAttribute('data-sheet')===select.value;",
    "    panel.hidden=!active;", paste0(
        "    if(active){var viewer=panel.querySelector(",
        "'.dnaepico-result-content');if(viewer){",
            "viewer.dispatchEvent(new CustomEvent('dnaepico:",
        "activate-viewer'));}}"
    ), "  });", "}", "select.addEventListener('change',activate);",
    "activate();", paste0(
        "if(document.readyState==='loading'){",
            "document.addEventListener('DOMContentLoaded',",
        "activate,{once:true});}"
    ), "})();", "</script>", "```", ""
    )
    .dnamReportFunctionResult <- lines
    environment()
}

.drs_build_workbook_table_section_stagesDnaEpico <- c(
    .drs_build_workbook_table_section_part_01DnaEpico,
    .drs_build_workbook_table_section_part_02DnaEpico
)

.drs_prepare_csv_table_assets_part_01DnaEpico <-
    function() {
    .runDrStages(base::as.list.environment(environment(), all.names =
        TRUE),
    .drs2_prepare_csv_table_assets_part_01DnaEpico_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.drs_prepare_csv_table_assets_part_02DnaEpico <-
    function() {
    .installDrHelpers(environment())
    csv_target <- file.path(result_dir, csv_name)
    staged <- same_report_path(source_data_path, csv_target) || isTRUE(
        file.copy(source_data_path,
    csv_target,
    overwrite = TRUE
    ))
    if (!staged || !file.exists(csv_target)) {{ .dnamReportReturnValue <- list(
    ok = FALSE, error = paste0(
        "Could not copy the data file into ",
        "the report assets from `", slash(source_data_path), "`."
    ), source_path = source_data_path,
    metadata = list()
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
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
    table_data = table_data, var_prefix = var_prefix,
    id_column = id_column, downloads = list(list(href = paste0(
        "assets/results/",
        var_prefix, "/", csv_name
    ), label = "Download data (CSV)")), chunk_size = chunk_size,
    item_singular = "record", item_plural = "records", search_label =
        search_label,
    search_placeholder = if (nzchar(first_id)) {
        paste("e.g.", first_id)
    } else {
        "Enter an identifier"
    } )
    .dnamReportFunctionResult <- c(list(
    ok = TRUE, error = NULL, source_path = source_data_path,
    metadata = list(), csv_path = paste0(
        "assets/results/", var_prefix, "/",
        csv_name
    )
    ), viewer_assets)
    environment()
}

.drs_prepare_csv_table_assets_stagesDnaEpico <- c(
    .drs_prepare_csv_table_assets_part_01DnaEpico,
    .drs_prepare_csv_table_assets_part_02DnaEpico
)

.drs_model_figure_description_part_01DnaEpico <-
    function() {
    .installDrHelpers(environment())
    original_stem <- tools::file_path_sans_ext(basename(filename))
    stem <- tolower(original_stem)
    if (grepl("^manhattan_", stem)) {
    if (grepl("_v1$", stem)) {{ .dnamReportReturnValue <- paste0(
        "Circular genome-wide view of association p-values by ",
        "chromosome. Alternating chromosome colours support ",
            "genomic localisation without a separate legend."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    {
        .dnamReportReturnValue <- paste0(
        "Genome-wide association p-values ordered by chromosome and ",
        "position, with reference lines and labels highlighting the ",
        "strongest CpG signals."
        )
        .dnamReportDidReturn <- TRUE
        return(environment())
    }
    }
    if (grepl("^vennd", stem)) {{ .dnamReportReturnValue <-
        venn_figure_metadata(filename)$description
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^intersection", stem)) {{ .dnamReportReturnValue <- paste0(
    "Ranked intersection layout showing shared and term-specific ",
    "CpG or gene sets when a Venn diagram becomes crowded."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^qqplot_", stem)) {{ .dnamReportReturnValue <- paste0(
    "Observed against expected -log10 p-values for the selected ",
    analysis, " test. Departures from the diagonal help identify ",
    "association signal or systematic test-statistic inflation."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    environment()
}

.drs_model_figure_description_part_02DnaEpico <-
    function() {
    .installDrHelpers(environment())
    if (grepl("^(hist|bar|distribution)_", stem)) {{ .dnamReportReturnValue <-
        paste0(
    "Distribution of the model variable, including its observed ",
    "range or factor-level balance and missing values where ", "present."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^missingness_", stem)) {{ .dnamReportReturnValue <- paste0(
    "Missing-data pattern across model variables, used to assess ",
    "whether exclusions may be systematic."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^correlation_", stem)) {{ .dnamReportReturnValue <- paste0(
    "Pairwise association structure among numeric model ",
    "variables, supporting review of redundancy and potential ",
    "collinearity."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^(association_|modelvariables_)", stem)) {{
        .dnamReportReturnValue <- paste0(
    "Joint display of variables entering the model, used to ",
    "inspect group balance, separation, and covariate structure."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^(participantobservationcount_|observations_)", stem)) {{
        .dnamReportReturnValue <- paste0(
    "Number of observations per participant, documenting the ",
    "repeated-measures support available to the mixed model."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    environment()
}

.drs_model_figure_description_part_03DnaEpico <-
    function() {
    .installDrHelpers(environment())
    if (grepl("^(timedistribution_|timepointdistribution_)", stem)) {{
        .dnamReportReturnValue <- paste0(
    "Distribution of observations across analysis timepoints, ",
    "used to assess longitudinal balance and sparse visits."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^(trajectory_|longitudinaltrajectory_)", stem)) {{
        .dnamReportReturnValue <- paste0(
    "Participant-level trajectories across time, with the ",
    "population pattern shown to support interpretation of ",
        "within-person change."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^(change_|pairedchange_)", stem)) {{ .dnamReportReturnValue <-
        paste0(
    "Paired within-participant change between observed ",
    "timepoints, separating longitudinal change from ",
        "between-participant variability."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^effectforest_", stem)) {{ .dnamReportReturnValue <- paste0(
    "Top-ranked CpG effect estimates with 95% confidence ",
    "intervals, combining direction, magnitude, and precision."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    if (grepl("^volcano_", stem)) {{ .dnamReportReturnValue <- paste0(
    "Effect magnitude against statistical significance for all ",
    "tested CpGs, highlighting signals that combine both."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    environment()
}

.drs_model_figure_description_part_04DnaEpico <-
    function() {
    .installDrHelpers(environment())
    if (grepl("^(residual|standarderror|scalelocation|influence)", stem)) {{
        .dnamReportReturnValue <- paste0(
    "Model diagnostic used to assess residual behaviour, ",
    "variance stability, precision, or influential observations."
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    .dnamReportFunctionResult <- paste0("Analytical output for ", title,
        ", provided with its original high-resolution download.")
    environment()
}

.drs_model_figure_description_stagesDnaEpico <- c(
    .drs_model_figure_description_part_01DnaEpico,
    .drs_model_figure_description_part_02DnaEpico,
    .drs_model_figure_description_part_03DnaEpico,
    .drs_model_figure_description_part_04DnaEpico
)

.drs_venn_figure_metadata_part_01DnaEpico <-
    function() {
    .runDrStages(base::as.list.environment(environment(), all.names =
        TRUE),
    .drs2_venn_figure_metadata_part_01DnaEpico_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.drs_venn_figure_metadata_part_02DnaEpico <-
    function() {
    .installDrHelpers(environment())
    scope <- if (length(terms) && all(grepl("time", terms, ignore.case =
        TRUE))) {
    "Across Time-related Model Effects"
    } else {
    "Across Model Effects"
    }
    effects <- vapply(terms, venn_effect_label, character(1))
    .dnamReportFunctionResult <- list(title = trimws(paste(
    threshold_title, "CpG and",
    annotation_title, "Gene Overlap", scope
    )), description = sprintf(
    paste0("This Venn diagram compares %s and their %s for %s."),
    threshold_description, annotation_description, if (length(effects)) {
        sentence_list(effects)
    } else {
        "the selected model effects"
    }
    ))
    environment()
}

.drs_venn_figure_metadata_stagesDnaEpico <- c(
    .drs_venn_figure_metadata_part_01DnaEpico,
    .drs_venn_figure_metadata_part_02DnaEpico
)

.drs_summarize_dataset_part_01DnaEpico <-
    function() {
    .installDrHelpers(environment())
    data <- safe_read_table_file(data_path)
    summary <- list(
    path = slash(data_path), exists = !is.null(data), n_rows = NA_integer_,
    n_cols = NA_integer_, participant_col = NULL, n_participants = NA_integer_,
    timepoint_col = NULL, timepoints = character(), n_timepoints = NA_integer_,
    complete_timepoint_participants = NA_integer_,
        repeated_participant_timepoint_rows = NA_integer_,
    survey_columns = character(), participants_missing_survey = NA_integer_
    )
    if (is.null(data)) {{ .dnamReportReturnValue <- summary
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
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
    keys <- paste(participants[complete_rows], timepoints[complete_rows],
        sep = "\r" )
    summary$repeated_participant_timepoint_rows <- length(keys) - length(
        unique(keys))
    split_timepoints <- split(timepoints[complete_rows], participants[
        complete_rows])
    summary$complete_timepoint_participants <- sum(vapply(
        split_timepoints,
        function(values) {
        all(summary$timepoints %in% unique(as.character(values)))
        }, logical(1) )) }
    survey_columns <- grep("^(MHC_|BDSST|BRS_|SS_|WHO_)", names(data), value =
        TRUE)
    summary$survey_columns <- survey_columns
    environment() }

.drs_summarize_dataset_part_02DnaEpico <-
    function() {
    .installDrHelpers(environment())
    if (!is.null(participant_col) && length(survey_columns)) {
    participants <- as.character(data[[participant_col]])
    unique_participants <- sort_values(participants)
    missing_survey <- vapply(unique_participants, function(participant) {
        rows <- !is_blank(participants) & participants == participant
        values <- unlist(data[rows, survey_columns, drop = FALSE], use.names =
            FALSE)
        all(is_blank(values))
    }, logical(1))
    summary$participants_missing_survey <- sum(missing_survey)
    }
    .dnamReportFunctionResult <- summary
    environment()
}

.drs_summarize_dataset_stagesDnaEpico <- c(
    .drs_summarize_dataset_part_01DnaEpico,
    .drs_summarize_dataset_part_02DnaEpico
)

.drs_resolve_report_formula_records_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    empty_records <- data.frame(
    phenotype = character(), result_column = character(),
    formula = character(), stringsAsFactors = FALSE, check.names = FALSE
    )
    if (!is.data.frame(dictionary) || !("Formula" %in% names(dictionary))) {{
        .dnamReportReturnValue <- empty_records
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    keep <- !is.na(dictionary$Formula) & nzchar(trimws(as.character(
        dictionary$Formula)))
    if (!any(keep)) {{ .dnamReportReturnValue <- empty_records
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
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
    environment()
}

.drs_resolve_report_formula_records_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    phenotype_labels <- vapply(seq_along(formulas), function(index) {
    if (!is.na(explicit_phenotypes[[index]]) && nzchar(trimws(
        explicit_phenotypes[[index]]))) {
        return(trimws(explicit_phenotypes[[index]]))
    }
    if (length(known_phenotypes)) {
        matches <- known_phenotypes[vapply(known_phenotypes, function(
            phenotype) {
        grepl(paste0("`", phenotype, "`"), formulas[[index]], fixed = TRUE) ||
            startsWith(result_columns[[index]], phenotype)
        }, logical(1))]
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
    duplicate_key <- paste(records$phenotype, records$formula, sep = "\r")
    records <- records[!duplicated(duplicate_key), , drop = FALSE]
    rownames(records) <- NULL
    .dnamReportFunctionResult <- records
    environment()
}

.drs_resolve_report_formula_records_stagesDnaEpico <- c(
    .drs_resolve_report_formula_records_part_01DnaEpico,
    .drs_resolve_report_formula_records_part_02DnaEpico
)

.drs2_venn_figure_metadata_part_01DnaEpico_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    stem <- tools::file_path_sans_ext(basename(filename))
    if (!grepl("^vennD_", stem, ignore.case = TRUE)) {{
        .dnamReportReturnValue <- NULL
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    remainder <- sub("^vennD_(?:GLM|LME)_", "", stem, ignore.case = TRUE)
    tokens <- strsplit(remainder, "_", fixed = TRUE)[[1L]]
    annotation_idx <- grep("^(?:GENCODEv?[0-9]+|UCSC)$", tokens, ignore.case =
        TRUE)
    threshold_idx <- grep("^(?:genomeWide|nominal|suggestive|combined)$",
        tokens,
    ignore.case = TRUE
    )
    annotation <- if (length(annotation_idx)) {
    tokens[[annotation_idx[[length(annotation_idx)]]]]
    } else {
    ""
    }
    threshold <- if (length(threshold_idx)) {
    tokens[[threshold_idx[[length(threshold_idx)]]]]
    } else {
    ""
    }
    excluded <- unique(c(annotation_idx, threshold_idx))
    terms <- if (length(excluded)) {
    tokens[-excluded]
    } else {
    tokens
    }
    terms <- combine_venn_interaction_tokens(terms)
    environment()
}

.drs2_venn_figure_metadata_part_01DnaEpico_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    annotation_title <- if (grepl("^GENCODE", annotation, ignore.case = TRUE)) {
    version <- sub("^GENCODEv?", "", annotation, ignore.case = TRUE)
    paste("GENCODE v", version, sep = "")
    } else if (identical(toupper(annotation), "UCSC")) {
    "UCSC"
    } else {
    "Annotated"
    }
    annotation_description <- if (identical(annotation_title, "Annotated")) {
    "available gene annotations"
    } else {
    paste(annotation_title, "gene annotations")
    }
    threshold_title <- switch(tolower(threshold),
    genomewide = "Genome-wide",
    nominal = "Nominal-threshold",
    suggestive = "Suggestive-threshold",
    combined = "Combined-threshold",
    ""
    )
    threshold_description <- switch(tolower(threshold),
    genomewide = "genome-wide significant CpGs",
    nominal = "CpGs meeting the nominal threshold",
    suggestive = "CpGs meeting the suggestive threshold",
    combined = "CpGs meeting the selected thresholds",
    "selected CpGs"
    )
    .dnamReportFunctionResult <- environment()
    environment()
}

.drs2_venn_figure_metadata_part_01DnaEpico_stagesDnaEpico <- c(
    .drs2_venn_figure_metadata_part_01DnaEpico_part_01DnaEpico,
    .drs2_venn_figure_metadata_part_01DnaEpico_part_02DnaEpico
)

.drs2_prepare_csv_table_assets_part_01DnaEpico_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    source_data_path <- if (grepl("^[A-Za-z]:[/\\\\]|^/", data_path)) {
    data_path
    } else {
    file.path(root_dir, data_path)
    }
    if (!file.exists(source_data_path)) {{ .dnamReportReturnValue <- list(
    ok = FALSE, error = paste0(
        "Data file not found at `",
        slash(source_data_path), "`."
    ), source_path = source_data_path,
    metadata = list()
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    table_data <- tryCatch(utils::read.csv(source_data_path, check.names =
        FALSE),
    error = function(e) e
    )
    if (inherits(table_data, "error")) {{ .dnamReportReturnValue <- list(
    ok = FALSE, error = conditionMessage(table_data),
    source_path = source_data_path, metadata = list()
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    table_data <- as.data.frame(table_data, stringsAsFactors = FALSE,
        check.names = FALSE)
    front_columns <- front_columns[nzchar(front_columns) & front_columns %in%
    names(table_data)]
    environment()
}

.drs2_prepare_csv_table_assets_part_01DnaEpico_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    table_data <- table_data[, c(front_columns, setdiff(names(table_data),
        front_columns)),
    drop = FALSE
    ]
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
    table_data <- table_data[do.call(order, c(sort_values, list(na.last =
        TRUE))), ,
        drop = FALSE
    ]
    rownames(table_data) <- NULL
    }
    result_dir <- file.path(assets_dir, "results", var_prefix)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    csv_name <- basename(source_data_path)
    .dnamReportFunctionResult <- environment()
    environment()
}

.drs2_prepare_csv_table_assets_part_01DnaEpico_stagesDnaEpico <- c(
    .drs2_prepare_csv_table_assets_part_01DnaEpico_part_01DnaEpico,
    .drs2_prepare_csv_table_assets_part_01DnaEpico_part_02DnaEpico
)

.drs2_make_report_summary_items_part_02DnaEpico_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    detection <- if (!is.null(methylation_stage)) {
    log_field(methylation_stage$fields, "Detection p-value threshold")
    } else {
    ""
    }
    if (nzchar(normalization) || nzchar(detection)) {
    preprocessing <- sentence_list(c(if (nzchar(normalization)) {
        paste(normalization, "normalisation")
    }, if (nzchar(detection)) {
        paste("detection p-value threshold", detection)
    }))
    items <- c(items, sprintf("Preprocessing: %s.", preprocessing))
    }
    if (isTRUE(sample_summary$exists) && !is.na(sample_summary$n_samples)) {
    quality <- sprintf(
        "%s samples are represented in the detection summary",
        format_count(sample_summary$n_samples)
    )
    if (!is.na(sample_summary$min_p_detected) && !is.na(
        sample_summary$max_p_detected)) {
        quality <- paste0(
        quality, ", with detected-CpG percentages from ",
        format_decimal(sample_summary$min_p_detected), "% to ", format_decimal(
            sample_summary$max_p_detected),
        "%"
        )
    }
    items <- c(items, sprintf("Quality control: %s.", quality))
    } else if (isTRUE(cpg_summary$exists) && !is.na(
        cpg_summary$total_assessed)) {
    items <- c(items, sprintf(
        "Quality control: %s CpGs were assessed for detection.",
        format_count(cpg_summary$total_assessed)
    ))
    }
    batch_stage <- log_summary$stages$batch
    batch_k <- sva_summary$log_k
    if (is.na(batch_k) && !is.null(batch_stage)) {
    batch_k <- extract_report_count(log_field(batch_stage$fields, c(
        "Number of surrogate variables (K)",
        "Number of surrogate variables"
    )))
    }
    environment()
}

.drs2_make_report_summary_items_part_02DnaEpico_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    if (!is.na(batch_k)) {
    batch_text <- sprintf(
        "%s surrogate %s estimated to describe latent variation",
        format_count(batch_k), plural(batch_k, "variable")
    )
    if (!is.na(sva_summary$percent_variation)) {
        batch_text <- paste0(
        batch_text, " (", format_decimal(sva_summary$percent_variation),
        "% of recorded variation)"
        )
    }
    items <- c(items, sprintf("Batch-effect assessment: %s.", batch_text))
    }
    model_labels <- c(if ("glm" %in% model_sections) "GLM", if ("lme" %in%
        model_sections) "LME")
    if (length(model_labels)) {
    items <- c(items, sprintf(
        "Statistical analyses: %s result sections are available.",
        sentence_list(model_labels)
    ))
    }
    .dnamReportFunctionResult <- environment()
    environment()
}

.drs2_make_report_summary_items_part_02DnaEpico_stagesDnaEpico <- c(
    .drs2_make_report_summary_items_part_02DnaEpico_part_01DnaEpico,
    .drs2_make_report_summary_items_part_02DnaEpico_part_02DnaEpico
)

.drsg_setup_01DnaEpico_part_01DnaEpico <- function(
    ) {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    slash <- .slashDnamReport
    write_utf8 <- .writeUtf8DnamReport
    find_quarto <- .findQuartoDnamReport
    html_escape <- .htmlEscapeDnamReport
    r_string <- .rStringDnamReport
    slugify <- .slugifyDnamReport
    pretty_label <- .prettyLabelDnamReport
    callout_lines <- .calloutLinesDnamReport
    detPThreshold <- validateProbabilityDnaEpico(detPThreshold, "detPThreshold")
    model_sections <- normalizeModelSectionsDnamReport(modelSections)
    include_glm <- "glm" %in% model_sections
    include_lme <- "lme" %in% model_sections
    root_dir <- .inferWorkflowRootDnamReport(outputDir)
    project_name <- projectName
    qc_dir <- resolve_report_path(enmixTab)
    preprocessing_dir <- resolve_report_path(qcTab)
    postprocessing_dir <- resolve_report_path(metricTab)
    sva_dir <- resolve_report_path(svaTab)
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
    environment()
}

.drsg_setup_01DnaEpico_part_02DnaEpico <- function(
    ) {
    .installDrHelpers(environment())
    cpg_detection_path <- if (is.null(cpgDetectionPath) || !nzchar(
        cpgDetectionPath)) {
    resolve_report_path(file.path(
        "data", model_name, "preprocessingMinfiEwasWater",
        "cpgD.csv"
    ))
    } else {
    resolve_report_path(cpgDetectionPath)
    }
    sample_detection_path <- if (is.null(sampleDetectionPath) || !nzchar(
        sampleDetectionPath)) {
    resolve_report_path(file.path(
        "data", model_name, "preprocessingMinfiEwasWater",
        "sampleD.csv"
    ))
    } else {
    resolve_report_path(sampleDetectionPath)
    }
    .dnamReportFunctionResult <- environment()
    environment()
}

.drsg_setup_01DnaEpico_stagesDnaEpico <- c(
    .drsg_setup_01DnaEpico_part_01DnaEpico,
    .drsg_setup_01DnaEpico_part_02DnaEpico
)

.drsg_setup_02DnaEpico_part_01DnaEpico <- function(
    ) {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    glm_table_path <- if (is.null(glmTab) || !nzchar(glmTab)) {
    resolve_report_path(file.path("data", model_name, "methylationGLM",
        "annotatedGLM.xlsx"))
    } else {
    resolve_report_path(glmTab) }
    lme_table_path <- if (is.null(lmeTab) || !nzchar(lmeTab)) {
    resolve_report_path(file.path("data", model_name, "methylationLME",
        "annotatedLME.xlsx"))
    } else {
    resolve_report_path(lmeTab) }
    glm_dir <- resolve_report_path(file.path("figures", model_name,
        "methylationGLM"))
    lme_dir <- resolve_report_path(file.path("figures", model_name,
        "methylationLME"))
    logo_candidates <- c(
    logoPath, system.file("extdata", "dnaEPICORM.svg", package = "dnaEPICO"),
    file.path(root_dir, "inst", "extdata", "dnaEPICORM.svg"), file.path(
        root_dir,
        "dnaEPICORM.svg" ) )
    logo_candidates <- logo_candidates[nzchar(logo_candidates)]
    logo_source_path <- logo_candidates[file.exists(logo_candidates)][1]
    if (is.na(logo_source_path)) {
    logo_source_path <- "" }
    viewer_js_candidates <- c(
    system.file("extdata", "cpg-viewer.js", package = "dnaEPICO"),
    file.path(getwd(), "inst", "extdata", "cpg-viewer.js"), file.path(
        root_dir,
        "inst", "extdata", "cpg-viewer.js"
    ) )
    viewer_js_candidates <- viewer_js_candidates[nzchar(viewer_js_candidates)]
    viewer_js_source_path <- viewer_js_candidates[file.exists(
        viewer_js_candidates)][1]
    if (is.na(viewer_js_source_path)) {
    viewer_js_source_path <- ""
    }
    figure_viewer_js_candidates <- c(
    system.file("extdata", "figure-viewer.js",
        package = "dnaEPICO"
    ), file.path(getwd(), "inst", "extdata", "figure-viewer.js"),
    file.path(root_dir, "inst", "extdata", "figure-viewer.js")
    )
    figure_viewer_js_candidates <- figure_viewer_js_candidates[nzchar(
        figure_viewer_js_candidates)]
    figure_viewer_js_source_path <- figure_viewer_js_candidates[file.exists(
        figure_viewer_js_candidates)][1]
    environment() }

.drsg_setup_02DnaEpico_part_02DnaEpico <- function(
    ) {
    .installDrHelpers(environment())
    if (is.na(figure_viewer_js_source_path)) {
    figure_viewer_js_source_path <- ""
    }
    report_interactions_js_candidates <- c(
    system.file("extdata", "report-interactions.js",
        package = "dnaEPICO"
    ), file.path(getwd(), "inst", "extdata", "report-interactions.js"),
    file.path(root_dir, "inst", "extdata", "report-interactions.js")
    )
    report_interactions_js_candidates <- report_interactions_js_candidates[
        nzchar(report_interactions_js_candidates)]
    report_interactions_js_source_path <- report_interactions_js_candidates[
        file.exists(report_interactions_js_candidates)][1]
    if (is.na(report_interactions_js_source_path)) {
    report_interactions_js_source_path <- ""
    }
    project_dir <- resolve_report_path(outputDir)
    assets_dir <- file.path(project_dir, "assets")
    .dnamReportFunctionResult <- environment()
    environment()
}

.drsg_setup_02DnaEpico_stagesDnaEpico <- c(
    .drsg_setup_02DnaEpico_part_01DnaEpico,
    .drsg_setup_02DnaEpico_part_02DnaEpico
)

.drs2_build_figure_browser_part_01DnaEpico_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    if (!is.data.frame(items) || !nrow(items)) {{ .dnamReportReturnValue <-
        callout_lines(empty_message, type = "note")
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    titles <- if (!is.null(figure_titles) && length(figure_titles) == nrow(
        items)) {
    as.character(figure_titles)
    } else {
    vapply(items$original_name, model_figure_title, character(1), analysis =
        analysis)
    }
    titles <- sub("^Figure [0-9]+:[[:space:]]*", "", titles)
    descriptions <- if (!is.null(figure_descriptions) && length(
        figure_descriptions) ==
    nrow(items)) {
    as.character(figure_descriptions)
    } else {
    vapply(seq_len(nrow(items)), function(index) {
        model_figure_description(
        items$original_name[[index]], titles[[index]],
        analysis
        )
    }, character(1))
    }
    missing_descriptions <- is.na(descriptions) | !nzchar(descriptions)
    if (any(missing_descriptions)) {
    descriptions[missing_descriptions] <- vapply(
        which(missing_descriptions),
        function(index) {
        model_figure_description(
            items$original_name[[index]], titles[[index]],
            analysis
        )
        }, character(1)
    )
    }
    environment()
}

.drs2_build_figure_browser_part_01DnaEpico_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    options <- vapply(seq_len(nrow(items)), function(index) {
    sprintf(
        "<option value=\"%s\">%s %d: %s</option>", index, title_prefix,
        index, html_escape(titles[[index]])
    )
    }, character(1))
    objects <- vapply(seq_len(nrow(items)), function(index) {
    paste0(
        "{\"title\":", js_quote_result_values(paste0(
        title_prefix, " ",
        index, ": ", titles[[index]]
        )), ",\"previewPath\":", js_quote_result_values(items$asset_path[[
            index]]),
        ",\"downloadPath\":", js_quote_result_values(items$download_path[[
            index]]),
        ",\"downloadName\":", js_quote_result_values(items$download_name[[
            index]]),
        ",\"description\":", js_quote_result_values(descriptions[[index]]),
        ",\"browserReady\":", if (isTRUE(items$browser_ready[[index]])) {
        "true"
        } else {
        "false"
        }, "}"
    )
    }, character(1))
    figure_data <- paste0("[", paste(objects, collapse = ","), "]")
    .dnamReportFunctionResult <- environment()
    environment()
}

.drs2_build_figure_browser_part_01DnaEpico_stagesDnaEpico <- c(
    .drs2_build_figure_browser_part_01DnaEpico_part_01DnaEpico,
    .drs2_build_figure_browser_part_01DnaEpico_part_02DnaEpico
)

.drs_prepare_xlsx_workbook_assets_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    analysis <- match.arg(analysis, c("glm", "lme"))
    source_path <- if (grepl("^[A-Za-z]:[/\\\\]|^/", data_path)) {
    data_path
    } else {
    file.path(root_dir, data_path)
    }
    if (!file.exists(source_path)) {{ .dnamReportReturnValue <-
        prepare_xlsx_table_assets(data_path, primary_sheet,
    var_prefix,
    analysis = analysis, chunk_size = chunk_size, model_log_path =
        model_log_path
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    sheet_names <- tryCatch(openxlsx::getSheetNames(source_path), error =
        function(error) character())
    if (!length(sheet_names)) {{ .dnamReportReturnValue <-
        prepare_xlsx_table_assets(data_path, primary_sheet,
    var_prefix,
    analysis = analysis, chunk_size = chunk_size, model_log_path =
        model_log_path
    )
    .dnamReportDidReturn <- TRUE
    return(environment()) }}
    workbook_dir <- file.path(assets_dir, "workbooks", analysis)
    dir.create(workbook_dir, recursive = TRUE, showWarnings = FALSE)
    workbook_name <- basename(source_path)
    workbook_target <- file.path(workbook_dir, workbook_name)
    workbook_staged <- same_report_path(source_path, workbook_target) ||
        isTRUE(file.copy(source_path,
    workbook_target,
    overwrite = TRUE
    ))
    workbook_href <- if (workbook_staged && file.exists(workbook_target)) {
    paste0("assets/workbooks/", analysis, "/", workbook_name)
    } else {
    ""
    }
    environment()
}

.drs_prepare_xlsx_workbook_assets_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    workbook_download <- if (nzchar(workbook_href)) {
    list(list(href = workbook_href, label =
        "Download complete workbook (XLSX)"))
    } else {
    list()
    }
    ordered_sheets <- c(intersect(primary_sheet, sheet_names), setdiff(
    sheet_names,
    primary_sheet
    ))
    table_assets <- list()
    for (sheet_index in seq_along(ordered_sheets)) {
    sheet <- ordered_sheets[[sheet_index]]
    if (identical(sheet, primary_sheet)) {
        table_assets[[sheet]] <- prepare_xlsx_table_assets(
        data_path = source_path,
        sheet = sheet, var_prefix = var_prefix, analysis = analysis,
        chunk_size = chunk_size, model_log_path = model_log_path,
            additional_downloads = workbook_download
        )
        next
    }
    table_assets[[sheet]] <- .secondaryWorkbookSheetAssetsDnaEpico(
        source_path,
        sheet, var_prefix, analysis, chunk_size, workbook_download, assets_dir,
        slugify, write_table_viewer_assets
    )
    }
    primary_assets <- table_assets[[primary_sheet]]
    if (is.null(primary_assets)) {
    primary_assets <- table_assets[[1L]]
    }
    .dnamReportFunctionResult <- c(primary_assets, list(
    sheets = ordered_sheets,
    sheetAssets = table_assets, workbookPath = workbook_href, primarySheet =
        primary_sheet
    ))
    environment()
}

.drs_prepare_xlsx_workbook_assets_stagesDnaEpico <- c(
    .drs_prepare_xlsx_workbook_assets_part_01DnaEpico,
    .drs_prepare_xlsx_workbook_assets_part_02DnaEpico
)

.runDrStages <- function(arguments, stages, parent) {
    state <- list2env(arguments, parent = parent)
    state$.dnamReportDidReturn <- FALSE
    for (stage_template in stages) {
    stage <- stage_template
    environment(stage) <- state
    frame <- stage()
    if (!is.environment(frame)) {
        return(frame)
    }
    frame_values <- base::as.list.environment(frame, all.names = TRUE)
    base::list2env(frame_values, envir = state)
    if (isTRUE(state$.dnamReportDidReturn)) {
        return(state$.dnamReportReturnValue)
    }
    }
    result <- state$.dnamReportFunctionResult
    if (is.environment(result)) state else result
}

.drHelperBindings <- c(
    resolve_report_path = .dr_resolve_report_path,
    infer_model_from_path = .dr_infer_model_from_path,
    copy_log_asset = .dr_copy_log_asset,
    copy_first_existing_log_asset =
        .dr_copy_first_existing_log_asset,
    copy_figure_assets = .dr_copy_figure_assets,
    js_quote_result_values = .dr_js_quote_result_values,
    js_result_array = .dr_js_result_array,
    read_optional_workbook_sheet =
        .dr_read_optional_workbook_sheet,
    metadata_frame_to_list = .dr_metadata_frame_to_list,
    split_report_metadata_values =
        .dr_split_report_metadata_values,
    extract_report_formula_phenotype =
        .dr_extract_report_formula_phenotype,
    resolve_report_formula_records =
        .dr_resolve_report_formula_records,
    last_log_field = .dr_last_log_field,
    resolve_lme_report_metadata = .dr_resolve_lme_report_metadata,
    same_report_path = .dr_same_report_path,
    prepare_table_viewer_directory =
        .dr_prepare_table_viewer_directory,
    write_table_viewer_chunk = .dr_write_table_viewer_chunk,
    write_table_viewer_manifest = .dr_write_table_viewer_manifest,
    table_viewer_asset_result = .dr_table_viewer_asset_result,
    write_table_viewer_assets = .dr_write_table_viewer_assets,
    write_delimited_table_viewer_assets =
        .dr_write_delimited_table_viewer_assets,
    prepare_xlsx_table_assets = .dr_prepare_xlsx_table_assets,
    prepare_xlsx_workbook_assets =
        .dr_prepare_xlsx_workbook_assets,
    prepare_csv_table_assets = .dr_prepare_csv_table_assets,
    content_description_html = .dr_content_description_html,
    content_note_block = .dr_content_note_block,
    build_result_table_section = .dr_build_result_table_section,
    build_data_frame_table_section =
        .dr_build_data_frame_table_section,
    format_count = .dr_format_count,
    format_model_formula = .dr_format_model_formula,
    build_model_formula_notes = .dr_build_model_formula_notes,
    format_decimal = .dr_format_decimal,
    plural = .dr_plural,
    is_numeric_like = .dr_is_numeric_like,
    sort_values = .dr_sort_values,
    collapse_values = .dr_collapse_values,
    is_blank = .dr_is_blank,
    safe_read_table_file = .dr_safe_read_table_file,
    read_detection_tables = .dr_read_detection_tables,
    detection_table_warning = .dr_detection_table_warning,
    pick_column = .dr_pick_column,
    pick_participant_column = .dr_pick_participant_column,
    pick_timepoint_column = .dr_pick_timepoint_column,
    summarize_dataset = .dr_summarize_dataset,
    summarize_cpg_detection = .dr_summarize_cpg_detection,
    summarize_sample_detection = .dr_summarize_sample_detection,
    extract_sva_log_summary = .dr_extract_sva_log_summary,
    summarize_sva = .dr_summarize_sva,
    normalize_log_field = .dr_normalize_log_field,
    parse_log_fields = .dr_parse_log_fields,
    log_field = .dr_log_field,
    log_capture = .dr_log_capture,
    summarize_logs = .dr_summarize_logs,
    make_data_notes = .dr_make_data_notes,
    make_metrics_notes = .dr_make_metrics_notes,
    make_quality_control_notes = .dr_make_quality_control_notes,
    make_batch_effect_notes = .dr_make_batch_effect_notes,
    make_logs_notes = .dr_make_logs_notes,
    report_number_markup = .dr_report_number_markup,
    report_inline_markup = .dr_report_inline_markup,
    html_paragraph = .dr_html_paragraph,
    html_bullet_list = .dr_html_bullet_list,
    html_section = .dr_html_section,
    prepare_report_text = .dr_prepare_report_text,
    sentence_case = .dr_sentence_case,
    report_note_text = .dr_report_note_text,
    build_report_page = .dr_build_report_page,
    search_index_markup = .dr_search_index_markup,
    compose_page = .dr_compose_page,
    build_workbook_table_section = .dr_build_workbook_table_section,
    subset_figure_items = .dr_subset_figure_items,
    sentence_list = .dr_sentence_list,
    pretty_model_term = .dr_pretty_model_term,
    venn_effect_label = .dr_venn_effect_label,
    combine_venn_interaction_tokens = .dr_combine_venn_interaction_tokens,
    venn_figure_metadata = .dr_venn_figure_metadata,
    model_figure_title = .dr_model_figure_title,
    model_figure_description = .dr_model_figure_description,
    extract_report_count = .dr_extract_report_count,
    stage_sample_count = .dr_stage_sample_count,
    make_report_summary_items = .dr_make_report_summary_items,
    make_report_preprocessing_notes = .dr_make_report_preprocessing_notes,
    make_report_overlap_notes = .dr_make_report_overlap_notes,
    make_report_observations = .dr_make_report_observations,
    build_figure_browser = .dr_build_figure_browser,
    build_figure_sections = .dr_build_figure_sections,
    build_model_visualisation_tabs = .dr_build_model_visualisation_tabs,
    unrequested_table_assets = .dr_unrequested_table_assets,
    describe_available_annotations = .dr_describe_available_annotations,
    log_status_note = .dr_log_status_note,
    remove_unrequested_rendered_assets = .dr_remove_unrequested_rendered_assets
)

.installDrHelpers <- function(state) {
    for (name in names(.drHelperBindings)) {
    helper <- .drHelperBindings[[name]]
    environment(helper) <- state
    assign(name, helper, envir = state)
    }
    invisible(state)
}

.drg_setup_01DnaEpico <- function() {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
    .drsg_setup_01DnaEpico_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.drg_setup_02DnaEpico <- function() {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
    .drsg_setup_02DnaEpico_stagesDnaEpico,
    parent = environment(sys.function())
    )
}

.drg_setup_03DnaEpico <- function() {
    .installDrHelpers(environment())
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
        paste("Starting DNA Methylation Report Step:", project_name), paste(
        "Project root:",
        root_dir
        ), paste("Output project:", project_dir),
            "============================================================"
    ),
    verbose = FALSE, log_path = log_path
    )
    environment()
}

.drg_inventory_01DnaEpico <- function() {
    .installDrHelpers(environment())
    prepared_report <- prepareDnamReportInputs(
    outputDir = project_dir, qcDir = qc_dir,
    preprocessingDir = preprocessing_dir, postprocessingDir =
        postprocessing_dir,
    svaDir = sva_dir, glmDir = glm_dir, lmeDir = lme_dir, figDir = fig_dir,
    verbose = FALSE, logs = FALSE, logDir = logs_dir
    )
    enmix_items <- copy_figure_assets(qc_dir, asset_subdir = "enmix-qc")
    metrics_items <- copy_figure_assets(postprocessing_dir, asset_subdir =
        "metrics")
    qc_items <- copy_figure_assets(preprocessing_dir, asset_subdir =
        "quality-control")
    batch_items <- copy_figure_assets(sva_dir, asset_subdir = "batch-effects")
    glm_items <- if (include_glm) {
    copy_figure_assets(glm_dir, asset_subdir = "glm-visualisations")
    } else {
    data.frame()
    }
    environment()
}

.drg_summaries_01DnaEpico <- function() {
    .installDrHelpers(environment()); lme_items <- if (include_lme) {
        copy_figure_assets(lme_dir, asset_subdir = "lme-visualisations")
    } else { data.frame() }
    log_assets <- list(methylation = copy_log_asset(file.path(logs_dir,
        "log_preprocessingMinfiEwasWater.txt"),
            "preprocessingMinfiEwasWater.txt"),
        data = copy_log_asset(file.path(logs_dir, "log_preprocessingPheno.txt"),
            "preprocessingPheno.txt"), batch = copy_log_asset(file.path(
            logs_dir, "log_svaEnmix.txt"), "svaEnmix.txt"))
    if (include_glm) {
        log_assets$glm <- copy_first_existing_log_asset(file.path(logs_dir,
            "log_methylationGLM.txt"), "methylationGLM.txt")
    }; if (include_lme) {
        log_assets$lme <- copy_first_existing_log_asset(file.path(logs_dir,
            "log_methylationLME.txt"), "methylationLME.txt")
    }; data_summary <- summarize_dataset(pheno_file)
    detection_tables <- read_detection_tables(detp_path = detp_path,
        threshold = detPThreshold, cpg_path = cpg_detection_path,
        sample_path = sample_detection_path)
    cpg_detection_summary <- summarize_cpg_detection(detection_tables$cpg)
    sample_detection_summary <- summarize_sample_detection(
        detection_tables$sample)
    sva_summary <- summarize_sva(pheno_path = pheno_file, log_path = file.path(
        logs_dir, "log_svaEnmix.txt"))
    log_summary <- summarize_logs(log_assets)
    enmix_figure_titles <- enmix_items$title
    enmix_figure_descriptions <- vapply(enmix_figure_titles, function(title) {
        paste0("ENmix control-signal display for ", tolower(sub("[.]$",
            "", title)),
            ", used to assess expected assay-control behaviour and ",
            "identify atypical arrays.")
    }, character(1)); enmix_notes <- c(paste(paste0(
        "`ENmix` produces control plots comparable to ",
        "Illumina GenomeStudio plots."),
            "Infinium controls are interpreted by their expected signal range",
        "rather than absolute intensity to accommodate biological variation."))
    metrics_figure_titles <- metrics_items$title
    metrics_figure_descriptions <- rep("", nrow(metrics_items))
    metrics_metadata <- data.frame(pattern = c("densityBeta&M",
        "examineMDS", "cellComposition"), title = c(
            "Post-filtering Beta- and M-value Distributions by Timepoint",
        "Post-filtering Principal Component Analysis by Timepoint",
        "Estimated Cell-composition Distributions"), description = c(paste0(
            "beta-value and M-value density distributions across ",
        "samples after filtering, grouped by timepoint"), paste0(
            "a principal component analysis plot after filtering,",
        " grouped by timepoint"),
            "estimated cell-type proportions across all retained samples"),
        stringsAsFactors = FALSE); environment() }

.drx0031020_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    for (metadata_idx in seq_len(nrow(metrics_metadata))) {
        item_idx <- grep(metrics_metadata$pattern[[metadata_idx]],
            metrics_items$original_name,
            ignore.case = TRUE)
        if (length(item_idx)) {
            metrics_figure_titles[item_idx] <- metrics_metadata$title[[
            metadata_idx]]
            metrics_figure_descriptions[item_idx] <-
            metrics_metadata$description[[metadata_idx]]
        }
    }
    metrics_figure_titles <- paste0("Figure ", seq_along(
        metrics_figure_titles),
        ": ", sub("^Figure [0-9]+:[[:space:]]*", "", metrics_figure_titles))
    environment()
}

.drx0031020_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    quality_control_metadata <- data.frame(pattern = c("densityBeta",
        "detection_pvalues|detectionPvalue",
        "quality_control", "sexClinical", "sexComparison_RawNorm",
            "sexPrediction",
        "sampleRetention", "probeRetention"), title = c(
            "Beta Value Density Distribution (MSET)",
        "Detection P-values (RGSET)", "Quality Control Plot (MSET)",
            "Clinical Sex Distribution (GSET)",
        "Sex Comparison Before and After Normalisation (MSETF)",
            "Sex Prediction Plot (GSET)",
        "Sample Retention Across Preprocessing Filters",
            "CpG Retention Across Probe Filters"),
        description = c(paste0(
            "methylation beta-value distributions across samples ",
            "for identifying global patterns and sample outliers"), paste0(
            "detection p-value distributions across samples for ",
            "assessing probe-signal reliability and poor ", "detection"),
            paste0("sample-level methylated and unmethylated signal ",
            "intensities for identifying low-quality samples or ",
            "technical outliers"),
            "reported sex distribution for checking sample annotations",
            "sex-related methylation patterns before and after normalisation",
            paste0("methylation-based sex predictions for comparison ",
            "with reported sex and identification of possible ",
                "sample or annotation errors"), paste0(
            "sample counts retained after detection and sex-",
                "concordance filters"), paste0(
            "CpG counts retained after detection, chromosome, ",
                "SNP, and probe-exclusion filters")), stringsAsFactors = FALSE)
    quality_control_figure_titles <- qc_items$title
    quality_control_figure_descriptions <- rep("", nrow(qc_items))
    environment()
}

.drx0031020_part_03DnaEpico <- function() {
    .installDrHelpers(environment())
    for (metadata_idx in seq_len(nrow(quality_control_metadata))) {
        item_idx <- grep(quality_control_metadata$pattern[[metadata_idx]],
            qc_items$original_name,
            ignore.case = TRUE)
        if (length(item_idx)) {
            quality_control_figure_titles[item_idx] <-
            quality_control_metadata$title[[metadata_idx]]
            quality_control_figure_descriptions[item_idx] <-
            quality_control_metadata$description[[metadata_idx]]
        }
    }
    quality_control_figure_titles <- paste0("Figure ", seq_along(
        quality_control_figure_titles),
        ": ", sub("^Figure [0-9]+:[[:space:]]*", "",
            quality_control_figure_titles))
    batch_effect_figure_titles <- batch_items$title
    batch_effect_figure_descriptions <- rep("", nrow(batch_items))
    .dnamReportFunctionResult <- environment()
    environment()
}

.drx0031020_stagesDnaEpico <- c(
    .drx0031020_part_01DnaEpico,
    .drx0031020_part_02DnaEpico,
    .drx0031020_part_03DnaEpico
)

.drg_summaries_02DnaEpico <- function() {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
        .drx0031020_stagesDnaEpico,
        parent = environment(sys.function()))
}

.drg_summaries_03DnaEpico <- function() {
    .installDrHelpers(environment())
    batch_metadata <- data.frame(pattern = c(
    "sva_SentrixID.tiff|surrogateVariables_by_SentrixID.tiff",
    paste0("sva_SentrixIDPosition|surrogateVariableMatrix_by_Sen",
        "trixIDPosition"),
    "sva_SentrixPosition|surrogateVariables_by_SentrixPosition", paste0(
        "sva_technicalFactorAssociations|surrogateVariable_te",
        "chnicalFactorAssociations"
    )
    ), title = c(
    "First Two Surrogate Variables by Sentrix ID",
    "Surrogate-variable Matrix by Sentrix ID and Position",
        "First Two Surrogate Variables by Sentrix Position",
    "Surrogate-variable Associations with Technical Factors"
    ), description = c(paste0(
    "sample distribution across the first two surrogate ",
    "variables, coloured by Sentrix ID"
    ), paste0(
    "pairwise distributions of surrogate variables by ",
    "Sentrix ID and Sentrix position"
    ), paste0(
    "sample distribution across the first two surrogate ",
    "variables by Sentrix position"
    ), paste0(
    "ANOVA p-values describing associations between ",
    "surrogate variables and array technical factors"
    )), stringsAsFactors = FALSE)
    for (metadata_idx in seq_len(nrow(batch_metadata))) {
    item_idx <- grep(batch_metadata$pattern[[metadata_idx]],
        batch_items$original_name,
        ignore.case = TRUE
    )
    if (length(item_idx)) {
        batch_effect_figure_titles[item_idx] <- batch_metadata$title[[
            metadata_idx]]
        batch_effect_figure_descriptions[item_idx] <-
            batch_metadata$description[[metadata_idx]]
    }
    }
    batch_effect_figure_titles <- paste0(
    "Figure ", seq_along(batch_effect_figure_titles),
    ": ", sub("^Figure [0-9]+:[[:space:]]*", "", batch_effect_figure_titles)
    )
    environment()
}

.drg_tables_01DnaEpico <- function() {
    .installDrHelpers(environment())
    data_table_assets <- prepare_csv_table_assets(data_path = pheno_file,
        var_prefix = "phenotype_data", front_columns = c(if (is.null(
            data_summary$participant_col)) {
            ""
        } else {
            data_summary$participant_col
        }, if (is.null(data_summary$timepoint_col)) {
            ""
        } else {
            data_summary$timepoint_col
        }))
    glm_table_assets <- if (include_glm) {
        prepare_xlsx_workbook_assets(data_path = glm_table_path,
            primary_sheet = "annotatedGLM", var_prefix = "glm_results",
            analysis = "glm", model_log_path = file.path(logs_dir,
                "log_methylationGLM.txt"))
    } else {
        unrequested_table_assets("glm")
    }
    lme_table_assets <- if (include_lme) {
        prepare_xlsx_workbook_assets(data_path = lme_table_path,
            primary_sheet = "annotatedLME", var_prefix = "lme_results",
            analysis = "lme", model_log_path = file.path(logs_dir,
                "log_methylationLME.txt"))
    } else {
        unrequested_table_assets("lme")
    }
    emitLogMinfiEwasWater(c(paste("GLM report table source:", if (
        !include_glm) {
        "not requested"
    } else if (is.null(glm_table_assets$source_mode)) {
        "unavailable"
    } else {
        glm_table_assets$source_mode
    }), paste("LME report table source:", if (!include_lme) {
        "not requested"
    } else if (is.null(lme_table_assets$source_mode)) {
        "unavailable"
    } else {
        lme_table_assets$source_mode
    })), verbose = FALSE, log_path = log_path)
    lme_backend <- lme_table_assets$metadata$backend
    if (is.null(lme_backend) || !nzchar(lme_backend)) {
        lme_backend <- "lme4"
    }
    lme_backend <- tolower(lme_backend)
    environment()
}

.drg_tables_02DnaEpico <- function() {
    .installDrHelpers(environment())
    lme_analysis_label <- if (identical(lme_backend, "nlme")) {
    "nlme Analysis"
    } else {
    "LME Analysis"
    }
    lme_interaction_term <- lme_table_assets$metadata$interaction_term
    has_lme_interaction <- !is.null(lme_interaction_term) && nzchar(
        lme_interaction_term) &&
    !tolower(lme_interaction_term) %in% c("none", "null", "na")
    glm_table_title <- paste0(
    "Table 1. Generalised Linear Model Results and CpG ",
    "Annotation by Phenotype"
    )
    lme_table_title <- if (has_lme_interaction) {
    sprintf(paste0(
        "Table 1. Linear Mixed-Effects Model Results and CpG ",
        "Annotation by Phenotype and %s Interaction"
    ), lme_interaction_term)
    } else {
    paste0("Table 1. Linear Mixed-Effects Model Results and CpG ",
        "Annotation by Phenotype")
    }
    glm_table_description <- if (isTRUE(glm_table_assets$ok)) {
    glm_p_columns <- grep("P\\.Value$|P\\.value$", glm_table_assets$columns,
        value = TRUE
    )
    c(
        paste0("CpGs were analysed using a generalised linear model ",
            "fitted with `glm2`."),
        build_model_formula_notes(
        glm_table_assets$metadata$formula_records,
        "GLM"
        ), if (length(glm_p_columns)) {
        paste0("The result p-value column(s) are `", paste(glm_p_columns,
            collapse = "`, `"
        ), "`.")
        }, describe_available_annotations(glm_table_assets$columns)
    )
    } else {
    paste0("The annotated GLM result workbook was not available ",
        "when this report was generated.")
    }
    environment()
}

.lmeCorrelationDescriptionDnaEpico <- function(backend, structure, variable) {
    available <- identical(backend, "nlme") && !is.null(structure) &&
    nzchar(structure) &&
    !(tolower(structure) %in% c("none", "null", "na"))
    if (!available) {
    return(character())
    }
    if (is.null(variable) || !nzchar(variable)) {
    variable <- "the configured correlation variable"
    }
    sprintf(paste0(
    "The `%s` residual correlation structure orders ",
    "repeated observations using `%s`."
    ), structure, variable)
}

.lmeTableDescriptionDnaEpico <- function(
    assets, backend, hasInteraction, interactionTerm,
    formulaNotes, annotationNotes
) {
    if (!isTRUE(assets$ok)) {
    return(paste0(
        "The annotated LME result workbook was not available ",
        "when this report was generated."
    ))
    }
    p_columns <- grep("P\\.Value$|P\\.value$", assets$columns, value = TRUE)
    fitting <- if (identical(backend, "nlme")) {
    "CpGs were analysed using `nlme::lme()`."
    } else {
    paste0(
        "CpGs were analysed using `lmerTest::lmer()` with ",
        "the `lme4` mixed-effects framework."
    )
    }
    correlation <- .lmeCorrelationDescriptionDnaEpico(
    backend, assets$metadata$correlation_structure,
    assets$metadata$correlation_variable
    )
    interaction <- if (hasInteraction) {
    sprintf(paste0(
        "The reported phenotype terms test their interaction ",
        "with `%s`."
    ), interactionTerm)
    } else {
    paste0(
        "No interaction was fitted; the reported p-value ",
        "columns correspond to phenotype coefficients."
    )
    }
    model_label <- if (identical(backend, "nlme")) "nlme" else "LME"
    c(
    fitting, correlation, interaction,
    formulaNotes(assets$metadata$formula_records, model_label),
    if (length(p_columns)) {
        paste0(
        "The result p-value column(s) are `",
        paste(p_columns, collapse = "`, `"), "`."
        )
    },
    annotationNotes(assets$columns)
    )
}

.drg_tables_03DnaEpico <- function() {
    .installDrHelpers(environment())
    lme_table_description <- .lmeTableDescriptionDnaEpico(
    lme_table_assets, lme_backend, has_lme_interaction,
    lme_interaction_term, build_model_formula_notes,
    describe_available_annotations
    )
    environment()
}

.drx0030965_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    data_notes <- make_data_notes(data_summary)
    metrics_notes <- make_metrics_notes(metrics_items, metrics_figure_titles,
        metrics_figure_descriptions, data_summary)
    quality_control_notes <- make_quality_control_notes(qc_items,
        quality_control_figure_titles,
        quality_control_figure_descriptions, cpg_detection_summary,
            sample_detection_summary)
    batch_effect_notes <- make_batch_effect_notes(batch_items,
        batch_effect_figure_titles,
        batch_effect_figure_descriptions, sva_summary)
    glm_notes <- glm_table_description
    lme_notes <- lme_table_description
    logs_notes <- make_logs_notes(log_summary, lme_label = lme_analysis_label)
    glm_venn_items <- subset_figure_items(glm_items,
        "^(vennD|intersection)_GLM_")
    lme_venn_items <- subset_figure_items(lme_items,
        "^(vennD|intersection)_LME_")
    glm_visualisation_notes <- c(sprintf(
        "%s GLM visualisation files are available.",
        format_count(nrow(glm_items))),
            "Select an internal section and figure in the main panel.",
        if (nrow(glm_venn_items)) {
            sprintf("%s requested GLM Venn figures are available.",
            format_count(nrow(glm_venn_items)))
        } else {
            "No GLM Venn output was requested."
        })
    environment()
}

.drx0030965_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    lme_visualisation_notes <- c(sprintf(
        "%s LME visualisation files are available.",
        format_count(nrow(lme_items))), paste0(
            "Longitudinal design, trajectory, association, and ",
        "diagnostic figures are shown separately."), if (nrow(lme_venn_items)) {
        sprintf("%s requested LME Venn figures are available.", format_count(
            nrow(lme_venn_items)))
    } else {
        "No LME Venn output was requested."
    })
    report_summary_items <- make_report_summary_items(data_summary =
        data_summary,
        cpg_summary = cpg_detection_summary, sample_summary =
            sample_detection_summary,
        sva_summary = sva_summary, log_summary = log_summary, model_sections =
            model_sections,
        figure_items = c(list(enmix_items, qc_items, metrics_items,
            batch_items),
            if (include_glm) list(glm_items), if (include_lme) list(
            lme_items)),
        workbook_available = c(if (include_glm) isTRUE(glm_table_assets$ok),
            if (include_lme) isTRUE(lme_table_assets$ok)))
    report_preprocessing_notes <- make_report_preprocessing_notes(log_summary,
        logs_notes)
    report_overlap_notes <- make_report_overlap_notes(glm_venn_items,
        lme_venn_items)
    report_observations <- make_report_observations(data_summary =
        data_summary,
        cpg_summary = cpg_detection_summary, sample_summary =
            sample_detection_summary,
        log_summary = log_summary)
    environment()
}

.drx0030965_part_03DnaEpico <- function() {
    .installDrHelpers(environment())
    if (file.exists(logo_source_path)) {
        dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
        ok <- file.copy(logo_source_path, file.path(assets_dir,
            "dnaEPICORM.svg"),
            overwrite = TRUE)
        if (!ok) {
            warning("Failed to copy navbar logo: ", logo_source_path)
        }
    }
    .dnamReportFunctionResult <- environment()
    environment()
}

.drx0030965_stagesDnaEpico <- c(
    .drx0030965_part_01DnaEpico,
    .drx0030965_part_02DnaEpico,
    .drx0030965_part_03DnaEpico
)

.drg_narrative_01DnaEpico <- function() {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
        .drx0030965_stagesDnaEpico,
        parent = environment(sys.function()))
}

.drg_narrative_02DnaEpico <- function() {
    .installDrHelpers(environment())
    if (file.exists(viewer_js_source_path)) {
    dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(viewer_js_source_path, file.path(assets_dir,
        "cpg-viewer.js"),
        overwrite = TRUE
    )
    if (!ok) {
        warning("Failed to copy CpG viewer JavaScript: ", viewer_js_source_path)
    }
    } else {
    warning("CpG viewer JavaScript was not found in the installed package.")
    }
    if (file.exists(figure_viewer_js_source_path)) {
    dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(figure_viewer_js_source_path, file.path(assets_dir,
        "figure-viewer.js"),
        overwrite = TRUE
    )
    if (!ok) {
        warning("Failed to copy figure viewer JavaScript: ",
            figure_viewer_js_source_path)
    }
    } else {
    warning("Figure viewer JavaScript was not found in the installed package.")
    }
    if (file.exists(report_interactions_js_source_path)) {
    dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
    ok <- file.copy(report_interactions_js_source_path, file.path(
        assets_dir,
        "report-interactions.js"
    ), overwrite = TRUE)
    if (!ok) {
        warning("Failed to copy report interaction JavaScript: ",
            report_interactions_js_source_path)
    }
    } else {
    warning(
        "Report interaction JavaScript was not found in the installed ",
        "package."
    )
    }
    environment()
}

.drx0021616_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    model_navbar <- c(if (include_glm) {
        c("      - href: glm.qmd", "        text: \"GLM Analysis\"",
            "      - href: glm-visualisations.qmd",
            "        text: \"GLM Visualisations\"")
    }, if (include_lme) {
        c("      - href: lme.qmd", "        text: \"LME Analysis\"",
            "      - href: lme-visualisations.qmd",
            "        text: \"LME Visualisations\"")
    })
    environment()
}

.drx0021616_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    quarto_yml <- c("project:", "  type: website", "  output-dir: docs",
        "  resources:",
        "    - assets/", "", "website:", sprintf("  title: \"%s\"",
            project_name),
        "  search:", "    location: navbar", "    type: overlay", "  navbar:",
        "    logo: assets/dnaEPICORM.svg", "    right:",
            "      - href: index.qmd",
        "        text: \"Data\"", "      - href: quality-control.qmd",
            "        text: \"Quality Control\"",
        "      - href: batch-effect.qmd", "        text: \"Batch effect\"",
            "      - href: metrics.qmd",
        "        text: \"Metrics\"", model_navbar, "      - href: report.qmd",
        "        text: \"Report\"", "      - href: logs.qmd",
            "        text: \"Logs\"",
        "    tools:", "      - icon: github",
            "        href: https://github.com/paulYRP/dnaEPICO",
        "        aria-label: dnaEPICO on GitHub", "        target: _blank",
            "      - icon: hexagon",
        "        href: https://bioconductor.org/packages/dnaEPICO/",
            "        aria-label: dnaEPICO on Bioconductor",
        "        target: _blank", "  page-navigation: false",
            "  bread-crumbs: false",
        "  reader-mode: false", "", "format:", "  dashboard:",
            "    theme: cosmo",
        "    scrolling: true", "    orientation: columns",
            "    expandable: false",
        "    css:", "      - assets/qpasst.css", "    toc: false", "",
            "execute:",
        "  warning: false", "  message: false")
    site_css <- siteCssDnamReport()
    quality_figure_notes <- quality_control_notes[!grepl("^`Table [12]`",
        quality_control_notes)]
    cpg_table_notes <- quality_control_notes[grepl("^`Table 1`",
        quality_control_notes)]
    sample_table_notes <- quality_control_notes[grepl("^`Table 2`",
        quality_control_notes)]
    environment()
}

.drx0021616_part_03DnaEpico <- function() {
    .installDrHelpers(environment())
    data_page <- compose_page(title = "Data", body_classes =
        "qpasst-data-page",
        body_lines = c(build_result_table_section(title = "Data Preview",
            table_assets = data_table_assets,
            description = c(data_notes, paste0(
            "Phenotype, participant, and technical variables are ",
                "available to the selected analysis. Search and filters ",
            "operate on the complete table without placing every ",
                "record on the page.")))))
    metrics_page <- compose_page(title = "Metrics", body_lines = c(
        build_figure_browser(metrics_items,
        "No supported image files were found for the Metrics tab.",
            browser_id = "dnaepico-metrics-figures",
        analysis = "Metrics", figure_titles = metrics_figure_titles,
            figure_descriptions = metrics_figure_descriptions,
        browser_notes = metrics_notes)))
    .dnamReportFunctionResult <- environment()
    environment()
}

.drx0021616_stagesDnaEpico <- c(
    .drx0021616_part_01DnaEpico,
    .drx0021616_part_02DnaEpico,
    .drx0021616_part_03DnaEpico
)

.drg_pages_01DnaEpico <- function() {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
        .drx0021616_stagesDnaEpico,
        parent = environment(sys.function()))
}

.drx0021629_part_01DnaEpico <- function() {
    .installDrHelpers(environment())
    .installDrHelpers(environment())
    quality_control_page <- compose_page(title = "Quality Control",
        body_lines = c(build_figure_sections(sections = list(`ENmix QC` =
            enmix_items,
        `Methylation QC` = qc_items), empty_message = paste0(
            "No supported image files were found ",
        "for this ", "quality-control section."), browser_prefix =
            "dnaepico-quality-control",
        titles = list(`ENmix QC` = enmix_figure_titles, `Methylation QC` =
            quality_control_figure_titles),
        descriptions = list(`ENmix QC` = enmix_figure_descriptions,
            `Methylation QC` = quality_control_figure_descriptions),
        notes = list(`ENmix QC` = enmix_notes, `Methylation QC` =
            quality_figure_notes)),
        build_data_frame_table_section(title =
            "Table 1: CpG Detection Summary",
            data = detection_tables$cpg, empty_message =
            detection_table_warning(detection_tables),
            preview_rows = 10L, description = c(cpg_table_notes, paste0(
            "CpG-level detection counts and percentages across ",
                "samples identify consistently detected and ",
            "unreliable probes."))),
        build_data_frame_table_section(title =
            "Table 2: Sample Detection Summary",
            data = detection_tables$sample, empty_message =
            detection_table_warning(detection_tables),
            preview_rows = 100L, description = c(sample_table_notes, paste0(
            "Sample-level detected CpG counts and percentages ",
                "support review of assay completeness and possible ",
            "low-quality samples.")))))
    batch_effect_page <- compose_page(title = "Batch Effect", body_lines = c(
        build_figure_browser(batch_items,
        "No supported image files were found for the Batch Effect tab.",
            browser_id = "dnaepico-batch-effect-figures",
        analysis = "Batch Effect", figure_titles = batch_effect_figure_titles,
        figure_descriptions = batch_effect_figure_descriptions, browser_notes =
            batch_effect_notes)))
    environment()
}

.drx0021629_part_02DnaEpico <- function() {
    .installDrHelpers(environment())
    glm_page <- compose_page(title = "GLM Analysis", body_lines = c(
        build_workbook_table_section(title = glm_table_title,
        workbook_assets = glm_table_assets, model_notes = glm_notes)))
    glm_visualisations_page <- compose_page(title = "GLM Visualisations",
        body_lines = build_model_visualisation_tabs(glm_items,
        glm_venn_items, "GLM"))
    lme_page <- compose_page(title = lme_analysis_label, body_lines = c(
        build_workbook_table_section(title = lme_table_title,
        workbook_assets = lme_table_assets, model_notes = lme_notes)))
    lme_visualisations_page <- compose_page(title = "LME Visualisations",
        body_lines = build_model_visualisation_tabs(lme_items,
        lme_venn_items, "LME"))
    log_render_helper <- c("```{r}", "#| echo: false", "#| include: false",
        "render_log_block <- function(path, label) {",
        "  if (!file.exists(path)) {", "    cat('::: {.callout-warning}\\n')",
        "    cat(sprintf('Log file not found for `%s`.\\n', label))",
            "    cat(':::\\n')",
        "    return(invisible(NULL))", "  }",
            "  lines <- readLines(path, warn = FALSE, encoding = 'UTF-8')",
        "  if (!length(lines)) {", "    cat('::: {.callout-note}\\n')",
            "    cat(sprintf('`%s` is empty.\\n', label))",
        "    cat(':::\\n')", "    return(invisible(NULL))", "  }",
            "  cat('```text\\n')",
        "  cat(paste(lines, collapse = '\\n'))", "  cat('\\n```\\n')", "}",
            "```")
    environment()
}

.drx0021629_part_03DnaEpico <- function() {
    .installDrHelpers(environment())
    base_log_cards <- c("", "::: {.card title=\"Methylation Analysis\"}", "",
        content_note_block(c(paste0(
            "IDAT loading, normalisation, filtering, and cell-",
            "composition estimation."), log_status_note(
            "Methylation Analysis"))),
        "", "```{r}", "#| echo: false", "#| results: asis", sprintf(
            "render_log_block(%s, 'Methylation Analysis')",
            r_string(log_assets$methylation$asset_path)), "```", ":::", "",
            "::: {.card title=\"Data Preparation\"}",
        "", content_note_block(c(paste0(
            "Phenotype preparation, timepoint splitting, and ",
            "methylation-matrix export."), log_status_note(
            "Data Preparation"))),
        "", "```{r}", "#| echo: false", "#| results: asis", sprintf(
            "render_log_block(%s, 'Data Preparation')",
            r_string(log_assets$data$asset_path)), "```", ":::", "",
            "::: {.card title=\"Batch Effect\"}",
        "", content_note_block(c(
            "Surrogate-variable analysis for batch-effect assessment.",
            log_status_note("Batch Effect"))), "", "```{r}", "#| echo: false",
        "#| results: asis", sprintf("render_log_block(%s, 'Batch Effect')",
            r_string(log_assets$batch$asset_path)),
        "```", ":::")
    .dnamReportFunctionResult <- environment()
    environment()
}

.drx0021629_stagesDnaEpico <- c(
    .drx0021629_part_01DnaEpico,
    .drx0021629_part_02DnaEpico,
    .drx0021629_part_03DnaEpico
)

.drg_pages_02DnaEpico <- function() {
    .runDrStages(base::as.list.environment(environment(), all.names = TRUE),
        .drx0021629_stagesDnaEpico,
        parent = environment(sys.function()))
}

.drg_pages_03DnaEpico <- function() {
    .installDrHelpers(environment())
    glm_log_card <- if (include_glm) {
        c("", "::: {.card title=\"GLM Analysis\"}",
            "", content_note_block(c(paste0(
            "Generalised linear model fitting, association ",
                "testing, and CpG annotation."),
                log_status_note("GLM Analysis"))),
            "", "```{r}", "#| echo: false", "#| results: asis",
            sprintf("render_log_block(%s, 'GLM Analysis')",
                r_string(log_assets$glm$asset_path)),
            "```", ":::")
    } else { character(0) }
    lme_log_card <- if (include_lme) {
        c("", sprintf("::: {.card title=\"%s\"}",
            html_escape(lme_analysis_label)), "",
            content_note_block(c(if (has_lme_interaction) {
                sprintf(paste0(
            "Linear mixed-effects model fitting, the phenotype ",
                    "interaction with %s, and CpG annotation."),
                    lme_interaction_term)
            } else {
                paste0("Linear mixed-effects model fitting, phenotype ",
                    "coefficient testing, and CpG annotation.")
            }, log_status_note(lme_analysis_label))),
            "", "```{r}", "#| echo: false", "#| results: asis",
            sprintf("render_log_block(%s, %s)", r_string(
            log_assets$lme$asset_path),
                r_string(lme_analysis_label)), "```",
            ":::")
    } else {
        character(0) }
    logs_page <- compose_page(title = "Logs", body_classes = "qpasst-logs-page",
        body_lines = c(log_render_helper, base_log_cards,
            glm_log_card, lme_log_card))
    report_page <- build_report_page(project_name = project_name,
        project_dir = project_dir, data_notes = data_notes,
        enmix_notes = enmix_notes, quality_control_notes =
            quality_control_notes,
        batch_effect_notes = batch_effect_notes,
        metrics_notes = metrics_notes, glm_notes = glm_notes,
        lme_notes = lme_notes, logs_notes = logs_notes,
        glm_visualisation_notes = glm_visualisation_notes,
        lme_visualisation_notes = lme_visualisation_notes,
        summary_items = report_summary_items, preprocessing_notes =
            report_preprocessing_notes,
        overlap_notes = report_overlap_notes, observations =
            report_observations,
        lme_label = lme_analysis_label, model_sections = model_sections)
    environment() }

.drg_sources_01DnaEpico <- function() {
    .installDrHelpers(environment())
    unlink(c(file.path(project_dir, c("data.qmd", "DNAm.html")),
        file.path(project_dir, "docs", "data.html")), force = TRUE)
    unrequested_source_files <- c(if (!include_glm) {
        file.path(project_dir, c("glm.qmd", "glm-visualisations.qmd"))
    }, if (!include_lme) {
        file.path(project_dir, c("lme.qmd", "lme-visualisations.qmd"))
    })
    unlink(unrequested_source_files, force = TRUE)
    remove_unrequested_rendered_assets()
    write_utf8(file.path(project_dir, "_quarto.yml"), quarto_yml)
    write_utf8(file.path(assets_dir, "qpasst.css"), site_css)
    write_utf8(file.path(project_dir, "index.qmd"), data_page)
    unlink(c(file.path(project_dir, "enmix-qc.qmd"), file.path(project_dir,
        "docs", "enmix-qc.html")), force = TRUE)
    write_utf8(file.path(project_dir, "metrics.qmd"), metrics_page)
    write_utf8(file.path(project_dir, "quality-control.qmd"),
        quality_control_page)
    write_utf8(file.path(project_dir, "batch-effect.qmd"),
        batch_effect_page)
    if (include_glm) {
        write_utf8(file.path(project_dir, "glm.qmd"), glm_page)
        write_utf8(file.path(project_dir, "glm-visualisations.qmd"),
            glm_visualisations_page)
    }
    if (include_lme) {
        write_utf8(file.path(project_dir, "lme.qmd"), lme_page)
        write_utf8(file.path(project_dir, "lme-visualisations.qmd"),
            lme_visualisations_page)
    }
    write_utf8(file.path(project_dir, "report.qmd"), report_page)
    write_utf8(file.path(project_dir, "logs.qmd"), logs_page)
    quarto_bin <- find_quarto()
    source_files <- file.path(project_dir, c("_quarto.yml",
        "index.qmd", "metrics.qmd", "quality-control.qmd",
        "batch-effect.qmd", if (include_glm) "glm.qmd",
        if (include_glm) "glm-visualisations.qmd", if (include_lme) "lme.qmd",
        if (include_lme) "lme-visualisations.qmd", "report.qmd",
        "logs.qmd"))
    emitLogMinfiEwasWater(c(paste("Generated", length(source_files),
        "report source files in:", project_dir)), verbose = FALSE,
        log_path = log_path)
    environment()
}

.drg_render_01DnaEpico <- function() {
    .installDrHelpers(environment())
    render_status <- "skipped"
    rendered_file <- file.path(project_dir, "docs", "index.html")
    error_message <- NULL
    environment()
}

.quartoRenderFilesDnaEpico <- function(paths) {
    unlist(lapply(paths, function(path) {
    if (!file.exists(path)) {
        return(character())
    }
    readLines(path, warn = FALSE, encoding = "UTF-8")
    }), use.names = FALSE)
}

.quartoRenderStateDnaEpico <- function(
    status, renderedFile, projectDir, verbose, logPath
) {
    if (!identical(status, 0L)) {
    error <- paste("Quarto render failed with status", status)
    emitLogMinfiEwasWater(error, verbose = verbose, log_path = logPath)
    return(list(status = "failed", error = error))
    }
    required <- c(renderedFile, file.path(projectDir, "docs", "search.json"))
    missing <- required[!file.exists(required)]
    if (length(missing)) {
    error <- paste(
        "Quarto render did not create required website assets:",
        paste(basename(missing), collapse = ", ")
    )
    emitLogMinfiEwasWater(error, verbose = verbose, log_path = logPath)
    return(list(status = "failed", error = error))
    }
    emitLogMinfiEwasWater(
    c("Render complete."),
    verbose = FALSE, log_path = logPath
    )
    list(status = "rendered", error = NULL)
}

.runQuartoReportDnaEpico <- function(
    quartoBin, projectDir, renderedFile, verbose, logPath
) {
    emitLogMinfiEwasWater(
    c("Rendering report site..."),
    verbose = FALSE, log_path = logPath
    )
    old_directory <- getwd()
    on.exit(setwd(old_directory), add = TRUE)
    setwd(projectDir)
    stdout <- tempfile("dnaEPICO-quarto-stdout-", fileext = ".log")
    stderr <- tempfile("dnaEPICO-quarto-stderr-", fileext = ".log")
    on.exit(unlink(c(stdout, stderr)), add = TRUE)
    status <- system2(quartoBin, "render", stdout = stdout, stderr = stderr)
    output <- .quartoRenderFilesDnaEpico(c(stdout, stderr))
    if (length(output)) {
    emitLogMinfiEwasWater(
        c("Quarto output:", output),
        verbose = FALSE, log_path = logPath
    )
    }
    .quartoRenderStateDnaEpico(
    status, renderedFile, projectDir, verbose, logPath
    )
}

.missingQuartoReportStateDnaEpico <- function(projectDir, verbose, logPath) {
    error <- paste(
    paste0(
        "Quarto CLI not found. The project files were ",
        "created but not rendered."
    ),
    paste0(
        "Install Quarto or set QUARTO_BIN to the full path ",
        "of the quarto executable."
    )
    )
    emitLogMinfiEwasWater(
    c(error, projectDir),
    verbose = verbose, log_path = logPath
    )
    list(status = "failed", error = error)
}

.drg_render_02DnaEpico <- function() {
    .installDrHelpers(environment())
    render <- if (nzchar(quarto_bin)) {
    .runQuartoReportDnaEpico(
        quarto_bin, project_dir, rendered_file, verbose, log_path
    )
    } else {
    .missingQuartoReportStateDnaEpico(project_dir, verbose, log_path)
    }
    render_status <- render$status
    error_message <- render$error
    environment()
}

.drg_render_03DnaEpico <- function() {
    .installDrHelpers(environment())
    remove_unrequested_rendered_assets()
    if (!magick_available && any(grepl("\\.tiff?$", c(enmix_items$original_name,
        metrics_items$original_name, qc_items$original_name,
            batch_items$original_name,
        if (is.data.frame(glm_items) && nrow(glm_items)) {
            glm_items$original_name
        }, if (is.data.frame(lme_items) && nrow(lme_items)) {
            lme_items$original_name
        }), ignore.case = TRUE))) {
        emitLogMinfiEwasWater(c(paste0(
            "Note: TIFF figures were copied as-is because the ",
            "`magick` package is not installed."), paste0(
            "Install `magick` and rerun this function to ",
            "generate PNG browser previews.")), verbose = FALSE,
            log_path = log_path) }
    emitLogMinfiEwasWater(paste("Report path:", normalizePathDnamReport(
        file.path(project_dir,
        "docs", "index.html"))), verbose = verbose, log_path = log_path)
    render_result <- structure(list(preparedReport = prepared_report,
        status = render_status, renderedFile = if (identical(render_status,
            "rendered")) {
            rendered_file
        } else { NULL
        }, errorMessage = error_message, logFile = log_path), class =
            "dnaEPICO_dnamReport_render")
    .dnamReportResult <- structure(list(preparedReport = prepared_report,
        renderResult = render_result, status = render_status, outputFile = if (
            identical(render_status,
            "rendered")) {
            rendered_file } else {
            normalizePathDnamReport(file.path(project_dir, "docs",
                "index.html"))
        }, projectDir = project_dir, sourceFiles = source_files,
        resultTableSources = c(GLM = if (!include_glm) {
            "not_requested"
        } else if (is.null(glm_table_assets$source_mode)) {
            "unavailable" } else {
            glm_table_assets$source_mode
        }, LME = if (!include_lme) {
            "not_requested"
        } else if (is.null(lme_table_assets$source_mode)) {
            "unavailable"
        } else {
            lme_table_assets$source_mode
        }), modelSections = model_sections, docsDir = file.path(project_dir,
            "docs"), logoPath = logo_source_path, errorMessage = error_message,
        logFile = log_path), class = "dnaEPICO_dnamReport")
    environment() }

.drWorkflowStages <- c(
    .drg_setup_01DnaEpico,
    .drg_setup_02DnaEpico,
    .drg_setup_03DnaEpico,
    .drg_inventory_01DnaEpico,
    .drg_summaries_01DnaEpico,
    .drg_summaries_02DnaEpico,
    .drg_summaries_03DnaEpico,
    .drg_tables_01DnaEpico,
    .drg_tables_02DnaEpico,
    .drg_tables_03DnaEpico,
    .drg_narrative_01DnaEpico,
    .drg_narrative_02DnaEpico,
    .drg_pages_01DnaEpico,
    .drg_pages_02DnaEpico,
    .drg_pages_03DnaEpico,
    .drg_sources_01DnaEpico,
    .drg_render_01DnaEpico,
    .drg_render_02DnaEpico,
    .drg_render_03DnaEpico
)

#' Generate a DNA methylation dashboard report
#'
#' @details
#' The Quarto command-line interface is required to render the website. It is
#' not required to install or load `dnaEPICO`, or to use the package's
#' preprocessing and statistical-modeling functions.
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
#'   sidecars produced by `methylationGLM()` are read from this report
#'   project's `assets/results/glm_results` directory when available.
#' @param lmeTab Character or `NULL`. XLSX workbook shown in the LME Analysis
#' tab.
#'   When `NULL`, the path is inferred from the Makefile output layout. Report
#'   sidecars produced by `methylationLME()` are read from this report
#'   project's `assets/results/lme_results` directory when available.
#' @param modelSections Character vector containing any of `'glm'` and `'lme'`.
#'   The report includes the corresponding model pages, logs, and summary
#'   sections.
#'   Use `character(0)` for a preprocessing-only report. The default preserves
#'   the complete GLM-and-LME report.
#' @param logTab Character. Directory containing workflow logs shown in the
#'   Logs tab.
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
#' @param logoPath Character. Path to the navigation-panel logo. Defaults to
#'   the packaged `inst/extdata/dnaEPICORM.svg` asset.
#' @param imagePattern Character. Regular expression used to identify image
#'   files inside the section directories.
#' @param recursive Logical. If `TRUE`, search section directories recursively.
#' @return A list with class `'dnaEPICO_dnamReport'`.
#'
#' @examples
#' report_root <- file.path(tempdir(), "dnaepico-dnam-report")
#' pheno_file <- file.path(
#'   report_root,
#'   "data",
#'   "model1",
#'   "preprocessingMinfiEwasWater",
#'   "phenoLC.csv"
#' )
#' dir.create(dirname(pheno_file), recursive = TRUE, showWarnings = FALSE)
#' utils::write.csv(
#'   data.frame(
#'     UID = c("sample1", "sample2"),
#'     Timepoint = c(1, 2),
#'     Sex = c("F", "M")
#'   ),
#'   pheno_file,
#'   row.names = FALSE
#' )
#'
#' result <- dnamReport(
#'   outputDir = file.path(report_root, "reports", "model1"),
#'   phenoTab = pheno_file,
#'   enmixTab = file.path(
#'     report_root,
#'     "figures",
#'     "model1",
#'     "preprocessingMinfiEwasWater",
#'     "enmix"
#'   ),
#'   qcTab = file.path(
#'     report_root,
#'     "figures",
#'     "model1",
#'     "preprocessingMinfiEwasWater",
#'     "qc"
#'   ),
#'   svaTab = file.path(report_root, "figures", "model1", "svaEnmix"),
#'   metricTab = file.path(
#'     report_root,
#'     "figures",
#'     "model1",
#'     "preprocessingMinfiEwasWater",
#'     "metrics"
#'   ),
#'   logTab = file.path(report_root, "logs", "model1")
#' )
#' result$status
#'
#' @export
dnamReport <- function(
    outputDir = "reports", phenoTab = NULL, enmixTab = file.path(
    "figures",
    "preprocessingMinfiEwasWater", "enmix"
    ), qcTab = file.path(
    "figures", "preprocessingMinfiEwasWater",
    "qc"
    ), svaTab = file.path("figures", "svaEnmix"), metricTab = file.path(
    "figures",
    "preprocessingMinfiEwasWater", "metrics"
    ), glmTab = NULL, lmeTab = NULL,
    modelSections = c("glm", "lme"), logTab = outputDir, verbose = FALSE,
        logs = FALSE,
    projectName = "dnaEPICO", detPPath = NULL, detPThreshold = 0.01,
        cpgDetectionPath = NULL,
    sampleDetectionPath = NULL, logoPath = system.file("extdata",
        "dnaEPICORM.svg",
    package = "dnaEPICO"
    ), imagePattern = "\\.(png|jpg|jpeg|gif|webp|svg|tif|tiff)$",
    recursive = TRUE
) {
    old_options <- options(stringsAsFactors = FALSE)
    on.exit(options(old_options), add = TRUE)
    report_values <- base::as.list.environment(environment(), all.names = TRUE)
    state <- base::list2env(report_values, parent =
        environment(dnamReport))
    .installDrHelpers(state)
    for (stage_template in .drWorkflowStages) {
    stage <- stage_template
    environment(stage) <- state
    frame <- stage()
    frame_values <- base::as.list.environment(frame, all.names = TRUE)
    base::list2env(frame_values, envir = state)
    }
    state$.dnamReportResult
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
