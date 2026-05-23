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
  if (is.null(path) || !length(path) || is.na(path[[1]]) || !nzchar(path[[1]])) {
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
    directory,
    patterns,
    label,
    recursive = TRUE
) {
  normalized_dir <- normalizePathDnamReport(directory)
  files <- character(0)
  directory_exists <- dir.exists(directory)

  if (directory_exists) {
    file_sets <- lapply(
      patterns,
      function(pattern) {
        list.files(
          directory,
          pattern = pattern,
          full.names = TRUE,
          recursive = recursive,
          ignore.case = TRUE
        )
      }
    )
    files <- sort(unique(unlist(file_sets, use.names = FALSE)))
    files <- vapply(files, normalizePathDnamReport, character(1))
  }

  list(
    label = label,
    directory = normalized_dir,
    exists = directory_exists,
    files = files,
    count = length(files)
  )
}

#' Prepare inputs for a DNA methylation report
#'
#' @param outputDir Character. Directory where the report project is written.
#' @param qcDir Character. Directory containing ENmix quality-control figures.
#' @param preprocessingDir Character. Directory containing preprocessing
#'   quality-control figures.
#' @param postprocessingDir Character. Directory containing postprocessing metric
#'   figures.
#' @param svaDir Character. Directory containing SVA or batch-effect figures.
#' @param glmDir Character. Directory containing GLM figures.
#' @param glmmDir Character. Directory containing GLMM figures.
#' @param figDir Character. Directory used for generated report figure assets.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write progress messages to
#'   `file.path(logDir, "log_dnamReport.txt")`.
#' @param logDir Character. Directory for optional log files.
#'
#' @return A list with class `"dnaEPICO_dnamReport_prepared"`.
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
prepareDnamReportInputs <- function(
    outputDir = "reports",
    qcDir = file.path("figures", "preprocessingMinfiEwasWater", "enmix"),
    preprocessingDir = file.path("figures", "preprocessingMinfiEwasWater", "qc"),
    postprocessingDir = file.path("figures", "preprocessingMinfiEwasWater", "metrics"),
    svaDir = file.path("figures", "svaEnmix"),
    glmDir = file.path("figures", "methylationGLM_T1"),
    glmmDir = file.path("figures", "methylationGLMM_T1T2"),
    figDir = file.path(outputDir, "assets", "figures"),
    verbose = FALSE,
    logs = FALSE,
    logDir = outputDir
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = logDir,
    log_file = "log_dnamReport.txt"
  )

  figure_inventory <- list(
    qc = collectFigureInventoryDnamReport(
      directory = qcDir,
      patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$", "\\.tif$", "\\.tiff$"),
      label = "ENmix QC"
    ),
    preprocessing = collectFigureInventoryDnamReport(
      directory = preprocessingDir,
      patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$", "\\.tif$", "\\.tiff$"),
      label = "Quality control"
    ),
    postprocessing = collectFigureInventoryDnamReport(
      directory = postprocessingDir,
      patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$", "\\.tif$", "\\.tiff$"),
      label = "Postprocessing"
    ),
    sva = collectFigureInventoryDnamReport(
      directory = svaDir,
      patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$", "\\.tif$", "\\.tiff$"),
      label = "SVA"
    ),
    glm = collectFigureInventoryDnamReport(
      directory = glmDir,
      patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$", "\\.tif$", "\\.tiff$"),
      label = "GLM"
    ),
    glmm = collectFigureInventoryDnamReport(
      directory = glmmDir,
      patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$", "\\.svg$", "\\.tif$", "\\.tiff$"),
      label = "GLMM"
    )
  )

  missing_directories <- vapply(
    figure_inventory,
    function(section) !isTRUE(section$exists),
    logical(1)
  )
  missing_directories <- names(missing_directories)[missing_directories]

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      "Prepared DNA methylation dashboard report inputs",
      paste("Output directory:", normalizePathDnamReport(outputDir)),
      paste("Figure directory:", normalizePathDnamReport(figDir)),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  output_file <- normalizePathDnamReport(file.path(outputDir, "docs", "index.html"))

  structure(
    list(
      output = basename(output_file),
      outputDir = normalizePathDnamReport(outputDir),
      outputFile = output_file,
      figDir = normalizePathDnamReport(figDir),
      figureInventory = figure_inventory,
      missingFigureDirectories = missing_directories,
      logFile = log_path
    ),
    class = "dnaEPICO_dnamReport_prepared"
  )
}

#' Render a prepared DNA methylation report
#'
#' @param preparedReport Object returned by `prepareDnamReportInputs()`.
#' @param verbose Logical. If `TRUE`, emit progress messages.
#' @param logs Logical. If `TRUE`, write progress messages to a log file.
#' @param logDir Character or `NULL`. Directory for optional log files.
#' @param clean Logical. Retained for backwards compatibility.
#'
#' @return A list with class `"dnaEPICO_dnamReport_render"`.
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
    preparedReport,
    verbose = FALSE,
    logs = FALSE,
    logDir = NULL,
    clean = TRUE
) {
  if (!inherits(preparedReport, "dnaEPICO_dnamReport_prepared")) {
    stop(
      "preparedReport must be an object returned by prepareDnamReportInputs().",
      call. = FALSE
    )
  }

  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = if (is.null(logDir)) dirname(preparedReport$outputFile) else logDir,
    log_file = "log_dnamReport.txt"
  )

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      "renderDnamReport() is retained for compatibility.",
      "Use dnamReport() to generate the Quarto dashboard report.",
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      preparedReport = preparedReport,
      status = "skipped",
      renderedFile = NULL,
      errorMessage = "Use dnamReport() to generate the Quarto dashboard report.",
      logFile = log_path
    ),
    class = "dnaEPICO_dnamReport_render"
  )
}

#' Generate a DNA methylation dashboard report
#'
#' @param outputDir Character. Directory where the Quarto project is written.
#' @param phenoTab Character or `NULL`. CSV file shown in the Data tab.
#'   When `NULL`, the path is inferred from the Makefile output layout.
#' @param enmixTab Character. Directory containing ENmix quality-control figures.
#' @param qcTab Character. Directory containing Quality Control figures.
#' @param svaTab Character. Directory containing Batch Effect or SVA figures.
#' @param metricTab Character. Directory containing Metrics figures.
#' @param glmTab Character or `NULL`. XLSX workbook shown in the GLM Analysis tab.
#'   When `NULL`, the path is inferred from the Makefile output layout.
#' @param lmerTab Character or `NULL`. XLSX workbook shown in the LMER Analysis tab.
#'   When `NULL`, the path is inferred from the Makefile output layout.
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
#' @return A list with class `"dnaEPICO_dnamReport"`.
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
    outputDir = "reports",
    phenoTab = NULL,
    enmixTab = file.path("figures", "preprocessingMinfiEwasWater", "enmix"),
    qcTab = file.path("figures", "preprocessingMinfiEwasWater", "qc"),
    svaTab = file.path("figures", "svaEnmix"),
    metricTab = file.path("figures", "preprocessingMinfiEwasWater", "metrics"),
    glmTab = NULL,
    lmerTab = NULL,
    logTab = outputDir,
    verbose = FALSE,
    logs = FALSE,
    projectName = "dnaEPICO",
    detPPath = NULL,
    detPThreshold = 0.01,
    cpgDetectionPath = NULL,
    sampleDetectionPath = NULL,
    logoPath = system.file("extdata", "dnaEPICO.svg", package = "dnaEPICO"),
    imagePattern = "\\.(png|jpg|jpeg|gif|webp|svg|tif|tiff)$",
    recursive = TRUE
) {
old_options <- options(stringsAsFactors = FALSE)
on.exit(options(old_options), add = TRUE)

is_absolute_path <- function(path) {
  grepl("^[A-Za-z]:[/\\\\]|^/", path)
}

infer_workflow_root <- function(output_dir) {
  normalized_output <- if (is_absolute_path(output_dir)) {
    normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  } else {
    normalizePath(file.path(getwd(), output_dir), winslash = "/", mustWork = FALSE)
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
  normalizePath(file.path(base, path), winslash = "/", mustWork = FALSE)
}

project_name <- projectName
qc_dir <- resolve_report_path(enmixTab)
preprocessing_dir <- resolve_report_path(qcTab)
postprocessing_dir <- resolve_report_path(metricTab)
sva_dir <- resolve_report_path(svaTab)

infer_model_from_path <- function(path, step = "preprocessingMinfiEwasWater") {
  pieces <- strsplit(gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = FALSE)), "/", fixed = TRUE)[[1]]
  step_idx <- which(pieces == step)
  figures_idx <- which(pieces == "figures")
  if (length(step_idx) && length(figures_idx)) {
    usable <- figures_idx[figures_idx < step_idx[[1]]]
    if (length(usable) && step_idx[[1]] - usable[[length(usable)]] >= 2L) {
      return(pieces[[step_idx[[1]] - 1L]])
    }
  }
  basename(normalizePath(outputDir, winslash = "/", mustWork = FALSE))
}

model_name <- infer_model_from_path(qc_dir)
logs_dir <- resolve_report_path(logTab)

pheno_file <- if (is.null(phenoTab) || !nzchar(phenoTab)) {
  resolve_report_path(file.path("data", model_name, "preprocessingMinfiEwasWater", "phenoLC.csv"))
} else {
  resolve_report_path(phenoTab)
}
detp_path <- if (is.null(detPPath) || !nzchar(detPPath)) {
  resolve_report_path(file.path("rData", model_name, "preprocessingMinfiEwasWater", "qc", "detP_RGSet.RData"))
} else {
  resolve_report_path(detPPath)
}
cpg_detection_path <- if (is.null(cpgDetectionPath) || !nzchar(cpgDetectionPath)) {
  resolve_report_path(file.path("data", model_name, "preprocessingMinfiEwasWater", "cpgD.csv"))
} else {
  resolve_report_path(cpgDetectionPath)
}
sample_detection_path <- if (is.null(sampleDetectionPath) || !nzchar(sampleDetectionPath)) {
  resolve_report_path(file.path("data", model_name, "preprocessingMinfiEwasWater", "sampleD.csv"))
} else {
  resolve_report_path(sampleDetectionPath)
}
glm_table_path <- if (is.null(glmTab) || !nzchar(glmTab)) {
  resolve_report_path(file.path("data", model_name, "methylationGLM_T1", "annotatedGLM.xlsx"))
} else {
  resolve_report_path(glmTab)
}
lmer_table_path <- if (is.null(lmerTab) || !nzchar(lmerTab)) {
  resolve_report_path(file.path("data", model_name, "methylationGLMM_T1T2", "annotatedLME.xlsx"))
} else {
  resolve_report_path(lmerTab)
}
glm_dir <- resolve_report_path(file.path("figures", model_name, "methylationGLM_T1"))
glmm_dir <- resolve_report_path(file.path("figures", model_name, "methylationGLMM_T1T2"))
logo_candidates <- c(
  logoPath,
  system.file("extdata", "dnaEPICO.svg", package = "dnaEPICO"),
  file.path(root_dir, "inst", "extdata", "dnaEPICO.svg"),
  file.path(root_dir, "dnaEPICO.svg")
)
logo_candidates <- logo_candidates[nzchar(logo_candidates)]
logo_source_path <- logo_candidates[file.exists(logo_candidates)][1]
if (is.na(logo_source_path)) {
  logo_source_path <- ""
}

project_dir <- resolve_report_path(outputDir)
assets_dir <- file.path(project_dir, "assets")
assets_figures_dir <- file.path(assets_dir, "figures")
assets_logs_dir <- file.path(assets_dir, "logs")
fig_dir <- assets_figures_dir

image_pattern <- imagePattern
magick_available <- requireNamespace("magick", quietly = TRUE)
log_path <- resolveLogPathMinfiEwasWater(
  logs = logs,
  log_dir = logs_dir,
  log_file = "log_dnamReport.txt"
)

emitLogMinfiEwasWater(
  c(
    "=======================================================================",
    paste("Starting DNA Methylation Report Step:", project_name),
    paste("Project root:", root_dir),
    paste("Output project:", project_dir),
    "======================================================================="
  ),
  verbose = FALSE,
  log_path = log_path
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
  env_candidates <- normalizePath(env_candidates, winslash = "/", mustWork = FALSE)
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
    "C:/Program Files/Quarto/bin/quarto.exe",
    file.path(Sys.getenv("LOCALAPPDATA"), "Programs", "Quarto", "bin", "quarto.exe"),
    file.path(Sys.getenv("LOCALAPPDATA"), "quarto", "bin", "quarto.exe")
  )
  candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
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
  text <- gsub('"', "&quot;", text, fixed = TRUE)
  text
}

r_string <- function(text) {
  text <- gsub("\\\\", "\\\\\\\\", text)
  text <- gsub('"', '\\"', text, fixed = TRUE)
  paste0('"', text, '"')
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
    label <- gsub("([a-z0-9])([A-Z])", "\\1 \\2", label, perl = TRUE)
    label <- gsub("[_\\-]+", " ", label)
    label <- gsub("\\s+", " ", label)
    label <- trimws(label)
    words <- strsplit(tolower(label), "\\s+", perl = TRUE)[[1]]
    label <- paste0(toupper(substr(words, 1L, 1L)), substr(words, 2L, nchar(words)))
    label <- paste(label, collapse = " ")
    label <- gsub("\\bQc\\b", "QC", label)
    label <- gsub("\\bPca\\b", "PCA", label)
    label <- gsub("\\bSva\\b", "SVA", label)
    label <- gsub("\\bEnmix\\b", "ENmix", label)
    label
  }

callout_lines <- function(text, type = "note") {
  c(
    sprintf("::: {.callout-%s}", type),
    text,
    ":::",
    ""
  )
}

to_dashboard_data_path <- function(raw_path) {
  if (!nzchar(raw_path)) {
    return("")
  }

  if (grepl("^[A-Za-z]:[/\\\\]|^/", raw_path)) {
    return(slash(raw_path))
  }

  slash(file.path("..", "..", raw_path))
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
    exists = file.exists(src_path),
    source_path = slash(src_path),
    asset_path = slash(file.path("assets", "logs", output_name))
  )
}

copy_figure_assets <- function(src_dir, asset_subdir) {
  dir.create(file.path(assets_figures_dir, asset_subdir), recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(src_dir)) {
    return(
      data.frame(
        title = character(),
        original_name = character(),
        asset_path = character(),
        browser_ready = logical(),
        converted = logical(),
        stringsAsFactors = FALSE
      )
    )
  }

  files <- list.files(
    src_dir,
    pattern = image_pattern,
    full.names = TRUE,
    recursive = recursive,
    ignore.case = TRUE
  )
  files <- sort(files)

  if (!length(files)) {
    return(
      data.frame(
        title = character(),
        original_name = character(),
        asset_path = character(),
        browser_ready = logical(),
        converted = logical(),
        stringsAsFactors = FALSE
      )
    )
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
      dest_path <- file.path(assets_figures_dir, asset_subdir, dest_name)
      converted <- tryCatch(
        {
          img <- magick::image_read(src_path)
          magick::image_write(img, path = dest_path, format = "png")
          TRUE
        },
        error = function(e) FALSE
      )

      if (!isTRUE(converted)) {
        dest_name <- paste0(base_slug, ".", src_ext)
        dest_path <- file.path(assets_figures_dir, asset_subdir, dest_name)
        ok <- file.copy(src_path, dest_path, overwrite = TRUE)
        if (!ok) {
          stop("Failed to copy figure: ", src_path)
        }
        browser_ready <- FALSE
      }
    } else {
      dest_name <- paste0(base_slug, ".", src_ext)
      dest_path <- file.path(assets_figures_dir, asset_subdir, dest_name)
      ok <- file.copy(src_path, dest_path, overwrite = TRUE)
      if (!ok) {
        stop("Failed to copy figure: ", src_path)
      }
      browser_ready <- !(src_ext %in% c("tif", "tiff"))
    }

    rows[[idx]] <- data.frame(
      title = pretty_label(src_path),
      original_name = basename(src_path),
      asset_path = slash(file.path("assets", "figures", asset_subdir, dest_name)),
      browser_ready = browser_ready,
      converted = converted,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

build_gallery_markup <- function(
    items,
    empty_message,
    gallery_class = "qpasst-gallery",
    title_prefix = NULL,
    show_filename = TRUE
) {
  if (!nrow(items)) {
    return(callout_lines(empty_message, type = "warning"))
  }

  lines <- character()

  if (any(!items$browser_ready)) {
    lines <- c(
      lines,
      callout_lines(
        paste(
          "Some TIFF figures were copied without PNG conversion.",
          "Install the `magick` package and rerun `dnamReport()` if you want those figures to render directly in the browser."
        ),
        type = "warning"
      )
    )
  }

  lines <- c(lines, sprintf('<div class="%s">', gallery_class))

  for (idx in seq_len(nrow(items))) {
    row <- items[idx, , drop = FALSE]
    caption_title <- if (is.null(title_prefix)) {
      row$title
    } else {
      sprintf("%s %s", title_prefix, idx)
    }
    caption_html <- if (isTRUE(show_filename)) {
      sprintf(
        "    <figcaption><strong>%s</strong><br /><code>%s</code></figcaption>",
        html_escape(caption_title),
        html_escape(row$original_name)
      )
    } else {
      sprintf(
        "    <figcaption><strong>%s</strong></figcaption>",
        html_escape(caption_title)
      )
    }
    lines <- c(
      lines,
      '  <figure class="qpasst-gallery-item">',
      caption_html
    )

    if (isTRUE(row$browser_ready)) {
      lines <- c(
        lines,
        sprintf(
          '    <img src="%s" alt="%s" loading="lazy" />',
          row$asset_path,
          html_escape(row$title)
        )
      )
    } else {
      download_label <- if (isTRUE(show_filename)) row$original_name else caption_title
      lines <- c(
        lines,
        sprintf(
          '    <p class="qpasst-download">Browser preview is unavailable for this TIFF file. <a href="%s">Download %s</a>.</p>',
          row$asset_path,
          html_escape(download_label)
        )
      )
    }

    lines <- c(lines, "  </figure>")
  }

  c(lines, "</div>", "")
}

build_figure_cards <- function(
    items,
    empty_message,
    title_prefix = "Figure",
    show_filename = FALSE,
    figure_titles = NULL
) {
  if (!nrow(items)) {
    return(callout_lines(empty_message, type = "warning"))
  }

  lines <- character()

  if (any(!items$browser_ready)) {
    lines <- c(
      lines,
      callout_lines(
        paste(
          "Some TIFF figures were copied without PNG conversion.",
          "Install the `magick` package and rerun `dnamReport()` if you want those figures to render directly in the browser."
        ),
        type = "warning"
      )
    )
  }

  for (idx in seq_len(nrow(items))) {
    row <- items[idx, , drop = FALSE]
    caption_title <- if (!is.null(figure_titles) && length(figure_titles) >= idx && nzchar(figure_titles[[idx]])) {
      figure_titles[[idx]]
    } else if ("title" %in% names(row) && nzchar(row$title[[1]])) {
      row$title[[1]]
    } else {
      sprintf("%s %s", title_prefix, idx)
    }
    lines <- c(lines, sprintf("### %s", caption_title), "")

    if (isTRUE(row$browser_ready)) {
      lines <- c(
        lines,
        '<figure class="qpasst-figure-card">',
        sprintf(
          '  <img src="%s" alt="%s" loading="lazy" />',
          row$asset_path,
          html_escape(row$title)
        )
      )
      if (isTRUE(show_filename)) {
        lines <- c(
          lines,
          sprintf("  <figcaption><code>%s</code></figcaption>", html_escape(row$original_name))
        )
      }
      lines <- c(lines, "</figure>", "")
    } else {
      download_label <- if (isTRUE(show_filename)) row$original_name else caption_title
      lines <- c(
        lines,
        sprintf(
          '<p class="qpasst-download">Browser preview is unavailable for this TIFF file. <a href="%s">Download %s</a>.</p>',
          row$asset_path,
          html_escape(download_label)
        ),
        ""
      )
    }
  }

  lines
}

build_xlsx_table_section <- function(
    title,
    data_path,
    sheet,
    page_length = 10L,
    preview_rows = 25L,
    var_prefix = NULL,
    interactive = TRUE
) {
  if (is.null(var_prefix) || !nzchar(var_prefix)) {
    var_prefix <- gsub("-", "_", slugify(title))
  }

  source_data_path <- if (grepl("^[A-Za-z]:[/\\\\]|^/", data_path)) {
    data_path
  } else {
    file.path(root_dir, data_path)
  }

  if (!file.exists(source_data_path)) {
    return(c(
      sprintf("### %s", title),
      "",
      callout_lines(
        paste0("Data file not found at `", slash(source_data_path), "`."),
        type = "warning"
      )
    ))
  }

  table_data <- tryCatch(
    openxlsx::read.xlsx(source_data_path, sheet = sheet, check.names = FALSE),
    error = function(e) e
  )

  if (inherits(table_data, "error")) {
    return(c(
      sprintf("### %s", title),
      "",
      callout_lines(conditionMessage(table_data), type = "warning")
    ))
  }

  if (!isTRUE(interactive)) {
    if (is.finite(preview_rows) && nrow(table_data) > as.integer(preview_rows)) {
      table_data <- utils::head(table_data, as.integer(preview_rows))
    }

    header_cells <- paste0("    <th>", html_escape(names(table_data)), "</th>")
    body_rows <- unlist(
      lapply(seq_len(nrow(table_data)), function(i) {
        row_values <- vapply(
          table_data[i, , drop = FALSE],
          function(value) {
            if (is.na(value)) "" else as.character(value)
          },
          character(1)
        )
        c(
          "  <tr>",
          paste0("    <td>", html_escape(row_values), "</td>"),
          "  </tr>"
        )
      }),
      use.names = FALSE
    )

    return(c(
      sprintf("### %s", title),
      "",
      '<div class="table-responsive">',
      '  <table class="table table-striped table-sm">',
      "  <thead>",
      "  <tr>",
      header_cells,
      "  </tr>",
      "  </thead>",
      "  <tbody>",
      body_rows,
      "  </tbody>",
      "  </table>",
      "</div>",
      ""
    ))
  }

  path_var <- paste0(var_prefix, "_path")
  path_display_var <- paste0(var_prefix, "_path_display")
  data_var <- paste0(var_prefix, "_data")
  error_var <- paste0(var_prefix, "_error")
  sheet_var <- paste0(var_prefix, "_sheet")
  dashboard_data_path <- to_dashboard_data_path(data_path)

  c(
    sprintf("### %s", title),
    "",
    "```{r}",
    "#| echo: false",
    "#| include: false",
    sprintf("%s <- %s", path_var, r_string(dashboard_data_path)),
    sprintf("%s <- %s", path_display_var, r_string(dashboard_data_path)),
    sprintf("%s <- %s", sheet_var, r_string(sheet)),
    sprintf("%s <- NULL", data_var),
    sprintf("%s <- NULL", error_var),
    sprintf("if (!nzchar(%s)) {", path_var),
    sprintf("  %s <- 'No data file path is configured for this table.'", error_var),
    sprintf("} else if (!file.exists(%s)) {", path_var),
    sprintf("  %s <- paste0('Data file not found at `', %s, '`.')", error_var, path_display_var),
    "} else {",
    sprintf("  %s <- tryCatch(", data_var),
    sprintf("    openxlsx::read.xlsx(%s, sheet = %s, check.names = FALSE),", path_var, sheet_var),
    "    error = function(e) {",
    sprintf("      %s <<- conditionMessage(e)", error_var),
    "      NULL",
    "    }",
    "  )",
    "}",
    "```",
    "",
    "```{r}",
    "#| echo: false",
    "#| results: asis",
    sprintf("if (is.null(%s)) {", data_var),
    "  cat('::: {.callout-warning}\\n')",
    sprintf("  cat(%s, '\\n')", error_var),
    "  cat(':::\\n')",
    "}",
    "```",
    "",
    "```{r}",
    "#| echo: false",
    sprintf("if (!is.null(%s)) {", data_var),
    sprintf("  if (%s && requireNamespace('DT', quietly = TRUE)) {", if (isTRUE(interactive)) "TRUE" else "FALSE"),
    "    DT::datatable(",
    sprintf("      %s,", data_var),
    "      filter = 'top',",
    sprintf("      options = list(pageLength = %d, scrollX = TRUE, autoWidth = TRUE),", as.integer(page_length)),
    "      rownames = FALSE",
    "    )",
    "  } else {",
    "    knitr::kable(",
    sprintf("      utils::head(%s, %d),", data_var, as.integer(preview_rows)),
    "      format = 'html',",
    "      table.attr = \"class='table table-striped table-sm'\"",
    "    )",
    "  }",
    "}",
    "```",
    ""
  )
}

build_data_frame_table_section <- function(title, data, empty_message, preview_rows = 25L) {
  if (!is.data.frame(data) || !nrow(data)) {
    return(c(
      sprintf("### %s", title),
      "",
      callout_lines(empty_message, type = "warning")
    ))
  }

  if (is.finite(preview_rows) && nrow(data) > as.integer(preview_rows)) {
    data <- utils::head(data, as.integer(preview_rows))
  }

  header_cells <- paste0("    <th>", html_escape(names(data)), "</th>")
  body_rows <- unlist(
    lapply(seq_len(nrow(data)), function(i) {
      row_values <- vapply(
        data[i, , drop = FALSE],
        function(value) {
          if (is.na(value)) "" else as.character(value)
        },
        character(1)
      )
      c(
        "  <tr>",
        paste0("    <td>", html_escape(row_values), "</td>"),
        "  </tr>"
      )
    }),
    use.names = FALSE
  )

  c(
    sprintf("### %s", title),
    "",
    '<div class="table-responsive">',
    '  <table class="table table-striped table-sm">',
    "  <thead>",
    "  <tr>",
    header_cells,
    "  </tr>",
    "  </thead>",
    "  <tbody>",
    body_rows,
    "  </tbody>",
    "  </table>",
    "</div>",
    ""
  )
}

format_count <- function(value) {
  if (!length(value) || is.na(value)) {
    return("not available")
  }
  format(as.integer(round(value)), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_decimal <- function(value, digits = 2L) {
  if (!length(value) || is.na(value)) {
    return("not available")
  }
  format(round(as.numeric(value), digits), nsmall = digits, scientific = FALSE, trim = TRUE)
}

plural <- function(n, singular, plural_form = paste0(singular, "s")) {
  if (length(n) && !is.na(n) && as.integer(n) == 1L) singular else plural_form
}

is_numeric_like <- function(values) {
  values <- as.character(values)
  grepl(
    "^\\s*[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:[eE][+-]?\\d+)?\\s*$",
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
  paste0(paste(values[-length(values)], collapse = ", "), ", and ", values[[length(values)]])
}

is_blank <- function(values) {
  is.na(values) | !nzchar(trimws(as.character(values)))
}

safe_read_table_file <- function(path, sheet = 1) {
  if (!file.exists(path)) {
    return(NULL)
  }
  extension <- tolower(sub("^.*\\.([^.]+)$", "\\1", basename(path)))
  tryCatch(
    if (extension %in% c("xlsx", "xlsm")) {
      openxlsx::read.xlsx(path, sheet = sheet, check.names = FALSE)
    } else {
      utils::read.csv(path, check.names = FALSE)
    },
    error = function(e) NULL
  )
}

read_detection_tables <- function(detp_path, threshold, cpg_path = "", sample_path = "") {
  result <- list(
    source = "detP",
    path = slash(detp_path),
    exists = FALSE,
    threshold = threshold,
    cpg = NULL,
    sample = NULL,
    error = NULL
  )

  if (file.exists(detp_path)) {
    loaded <- tryCatch({
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
          },
          logical(1)
        )]
        if (length(matrix_objects)) get(matrix_objects[[1]], envir = detp_env) else NULL
      }

      if (is.null(detp_object)) {
        stop("No `detP` matrix-like object was found.", call. = FALSE)
      }

      detp_matrix <- as.matrix(detp_object)
      storage.mode(detp_matrix) <- "numeric"
      if (is.null(colnames(detp_matrix))) {
        colnames(detp_matrix) <- sprintf("Sample_%s", seq_len(ncol(detp_matrix)))
      }

      detected <- detp_matrix < threshold
      cpg <- data.frame(
        metric = c(
          "Total CpGs assessed",
          "CpGs detected in at least one sample",
          "CpGs detected in all samples",
          "CpGs never detected"
        ),
        nCpGs = c(
          nrow(detp_matrix),
          sum(rowSums(detected, na.rm = TRUE) >= 1),
          sum(rowSums(detected, na.rm = TRUE) == ncol(detp_matrix)),
          sum(rowSums(detected, na.rm = TRUE) == 0)
        ),
        stringsAsFactors = FALSE
      )
      sample <- data.frame(
        UID = colnames(detp_matrix),
        nDetected = colSums(detected, na.rm = TRUE),
        pDetected = round(100 * colMeans(detected, na.rm = TRUE), 2),
        stringsAsFactors = FALSE
      )
      list(cpg = cpg, sample = sample)
    }, error = function(e) e)

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
      detection_tables$path,
      "`: ",
      detection_tables$error
    ))
  }
  paste0("Detection P-value file not found at `", detection_tables$path, "`.")
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
  exact_col <- pick_column(
    data,
    c(
      "UID",
      "Participant_ID",
      "ParticipantID",
      "participant_id",
      "Participant",
      "Subject_ID",
      "SubjectID",
      "Individual_ID",
      "IndividualID",
      "ID"
    )
  )
  if (!is.null(exact_col)) {
    return(exact_col)
  }

  id_cols <- grep("id", names(data), ignore.case = TRUE, value = TRUE)
  if (length(id_cols)) {
    return(id_cols[[1]])
  }

  NULL
}

pick_timepoint_column <- function(data) {
  pick_column(
    data,
    c(
      "Timepoint",
      "Tiempoint",
      "Time_Point",
      "Timepoint_ID",
      "TimepointID",
      "Visit",
      "VisitID",
      "Visit_ID"
    )
  )
}

summarize_dataset <- function(data_path) {
  data <- safe_read_table_file(data_path)
  summary <- list(
    path = slash(data_path),
    exists = !is.null(data),
    n_rows = NA_integer_,
    n_cols = NA_integer_,
    participant_col = NULL,
    n_participants = NA_integer_,
    timepoint_col = NULL,
    timepoints = character(),
    n_timepoints = NA_integer_,
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
    keys <- paste(participants[complete_rows], timepoints[complete_rows], sep = "\r")
    summary$repeated_participant_timepoint_rows <- length(keys) - length(unique(keys))

    split_timepoints <- split(timepoints[complete_rows], participants[complete_rows])
    summary$complete_timepoint_participants <- sum(vapply(
      split_timepoints,
      function(values) all(summary$timepoints %in% unique(as.character(values))),
      logical(1)
    ))
  }

  survey_columns <- grep("^(MHC_|BDSST|BRS_|SS_|WHO_)", names(data), value = TRUE)
  summary$survey_columns <- survey_columns
  if (!is.null(participant_col) && length(survey_columns)) {
    participants <- as.character(data[[participant_col]])
    unique_participants <- sort_values(participants)
    missing_survey <- vapply(
      unique_participants,
      function(participant) {
        rows <- !is_blank(participants) & participants == participant
        values <- unlist(data[rows, survey_columns, drop = FALSE], use.names = FALSE)
        all(is_blank(values))
      },
      logical(1)
    )
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
    if (length(idx)) as.numeric(data$nCpGs[[idx[[1]]]]) else NA_real_
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
    summary$median_p_detected <- stats::median(p_detected, na.rm = TRUE)
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

  text <- paste(readLines(log_path, warn = FALSE, encoding = "UTF-8"), collapse = " ")
  match <- regexec("([0-9]+)\\s+surrogate variables explain\\s+([0-9.]+)\\s*%", text, perl = TRUE)
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
    path = slash(pheno_path),
    exists = !is.null(pheno_data),
    n_samples = NA_integer_,
    n_surrogate_variables = NA_integer_,
    surrogate_variables = character(),
    log_k = log_summary$k,
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
    sentrix_id_col <- pick_column(pheno_data, c("Sentrix_ID", "SentrixID", "Sentrix Id"))
    sentrix_position_col <- pick_column(pheno_data, c("Sentrix_Position", "SentrixPosition", "Sentrix Position"))
    if (!is.null(sentrix_id_col)) {
      summary$sentrix_id_levels <- length(sort_values(pheno_data[[sentrix_id_col]]))
    }
    if (!is.null(sentrix_position_col)) {
      summary$sentrix_position_levels <- length(sort_values(pheno_data[[sentrix_position_col]]))
    }
  }

  summary
}

summarize_logs <- function(log_assets) {
  labels <- c(
    methylation = "Methylation Analysis",
    data = "Data Preparation",
    batch = "Batch Effect",
    glm = "GLM",
    lmer = "LMER"
  )
  rows <- lapply(names(log_assets), function(name) {
    asset <- log_assets[[name]]
    line_count <- if (isTRUE(asset$exists)) {
      length(readLines(asset$source_path, warn = FALSE, encoding = "UTF-8"))
    } else {
      NA_integer_
    }
    data.frame(
      key = name,
      label = labels[[name]],
      exists = isTRUE(asset$exists),
      source_path = asset$source_path,
      line_count = line_count,
      stringsAsFactors = FALSE
    )
  })
  rows <- do.call(rbind, rows)
  list(
    rows = rows,
    total = nrow(rows),
    found = sum(rows$exists),
    total_lines = sum(rows$line_count, na.rm = TRUE)
  )
}

make_data_notes <- function(summary) {
  if (!isTRUE(summary$exists)) {
    return("The source CSV was not found.")
  }

  notes <- c(sprintf(
    "The table has %s rows and %s columns.",
    format_count(summary$n_rows),
    format_count(summary$n_cols)
  ))

  if (!is.null(summary$participant_col) && !is.null(summary$timepoint_col)) {
    notes <- c(notes, sprintf(
      "Using `%s` as the participant identifier and `%s` as the timepoint column, the file contains %s unique %s across timepoint values %s.",
      summary$participant_col,
      summary$timepoint_col,
      format_count(summary$n_participants),
      plural(summary$n_participants, "participant"),
      collapse_values(summary$timepoints)
    ))
  } else if (!is.null(summary$participant_col)) {
    notes <- c(notes, sprintf(
      "Using `%s` as the participant identifier, the file contains %s unique %s. No Timepoint column was detected.",
      summary$participant_col,
      format_count(summary$n_participants),
      plural(summary$n_participants, "participant")
    ))
  } else {
    notes <- c(notes, "No participant identifier column was detected. The report looks first for `UID`, then for other columns containing `ID`.")
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
  sub("^((?:Figure|Table)\\s+[0-9]+)([:.])", "`\\1`\\2", title, perl = TRUE)
}

make_figure_notes <- function(items, figure_titles, section_label, figure_descriptions = NULL) {
  n_items <- nrow(items)
  if (!n_items) {
    return(sprintf("No supported figures were found for the %s tab.", section_label))
  }

  vapply(
    seq_len(n_items),
    function(idx) {
      title <- if (length(figure_titles) >= idx && nzchar(figure_titles[[idx]])) {
        figure_titles[[idx]]
      } else if ("title" %in% names(items) && nzchar(items$title[[idx]])) {
        items$title[[idx]]
      } else {
        sprintf("Figure %s", idx)
      }
      highlighted_title <- highlight_report_label(title)
      description <- if (!is.null(figure_descriptions) && length(figure_descriptions) >= idx && nzchar(figure_descriptions[[idx]])) {
        figure_descriptions[[idx]]
      } else {
        describe_figure_title(title)
      }
      description <- sub("\\.$", "", description)
      if (grepl("^(presents|shows|summarises|compares)\\b", description, ignore.case = TRUE)) {
        sprintf("%s %s.", highlighted_title, description)
      } else {
        sprintf("%s presents %s.", highlighted_title, description)
      }
    },
    character(1)
  )
}

make_metrics_notes <- function(items, figure_titles, figure_descriptions, data_summary) {
  make_figure_notes(
    items,
    figure_titles,
    "Metrics",
    figure_descriptions
  )
}

make_quality_control_notes <- function(items, figure_titles, figure_descriptions, cpg_summary, sample_summary) {
  notes <- make_figure_notes(
    items,
    figure_titles,
    "Quality Control",
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
      format_count(sample_summary$n_samples),
      plural(sample_summary$n_samples, "sample")
    ))
  }

  notes
}

make_batch_effect_notes <- function(items, figure_titles, figure_descriptions, sva_summary) {
  notes <- make_figure_notes(
    items,
    figure_titles,
    "Batch Effect",
    figure_descriptions
  )

  if (!is.na(sva_summary$log_k) && !is.na(sva_summary$percent_variation)) {
    notes <- c(notes, sprintf(
      "The SVA log reports %s surrogate %s explaining %s%% of data variation.",
      format_count(sva_summary$log_k),
      plural(sva_summary$log_k, "variable"),
      format_decimal(sva_summary$percent_variation)
    ))
  }

  notes
}

make_logs_notes <- function(log_summary) {
  descriptions <- c(
    methylation = "displays methylation preprocessing, including IDAT loading, normalisation, filtering, and cell composition estimation",
    data = "displays phenotype preparation, timepoint splitting, and methylation matrix export steps",
    batch = "displays the hidden-effect and surrogate-variable analysis workflow",
    glm = "displays the generalised linear model workflow and CpG annotation steps",
    lmer = "displays the linear mixed-effects model workflow and CpG annotation steps"
  )
  labels <- c(
    methylation = "Methylation Analysis",
    data = "Data Preparation",
    batch = "Batch Effect",
    glm = "GLM Analysis",
    lmer = "LMER Analysis"
  )

  notes <- c(
    "The Logs tab displays workflow log files generated by the pipeline, with one section per analysis stage.",
    vapply(
    log_summary$rows$key,
    function(key) {
      status <- if (isTRUE(log_summary$rows$exists[log_summary$rows$key == key])) {
        "is available"
      } else {
        "was not found"
      }
      sprintf("`%s` %s and %s.", labels[[key]], status, descriptions[[key]])
    },
    character(1)
    )
  )

  notes
}

html_paragraph <- function(text) {
  sprintf("<p>%s</p>", html_escape(text))
}

html_section <- function(title, paragraphs) {
  c(sprintf("<h3>%s</h3>", html_escape(title)), html_paragraph(paragraphs))
}

strip_inline_markdown <- function(text) {
  text <- gsub("`([^`]+)`", "\\1", text, perl = TRUE)
  text <- gsub("\\s+", " ", text)
  trimws(text)
}

sentence_case <- function(text) {
  if (!nzchar(text)) {
    return(text)
  }
  paste0(toupper(substr(text, 1L, 1L)), substr(text, 2L, nchar(text)))
}

tab_report_paragraph <- function(notes, fallback) {
  notes <- strip_inline_markdown(notes)
  notes <- notes[nzchar(notes)]
  if (!length(notes)) {
    return(sprintf("%s.", sentence_case(fallback)))
  }

  paste(notes, collapse = " ")
}

build_report_page <- function(
    project_name,
    project_dir,
    data_notes,
    enmix_notes,
    quality_control_notes,
    batch_effect_notes,
    metrics_notes,
    glm_notes,
    lmer_notes,
    logs_notes
) {
  overview <- sprintf(
    "This report is generated automatically from the available datasets, figures, tables, and workflow logs for the project."
  )

  data_paragraph <- tab_report_paragraph(data_notes, "displays the phenotype data preview and detected participant/timepoint fields")
  enmix_paragraph <- tab_report_paragraph(enmix_notes, "displays ENmix control plots and interpretation notes")
  quality_control_paragraph <- tab_report_paragraph(quality_control_notes, "displays quality-control figures and detection summary tables")
  batch_effect_paragraph <- tab_report_paragraph(batch_effect_notes, "displays SVA and batch-effect figures")
  metrics_paragraph <- tab_report_paragraph(metrics_notes, "displays post-filtering methylation metric figures")
  glm_paragraph <- tab_report_paragraph(glm_notes, "displays the annotated generalised linear model results table")
  lmer_paragraph <- tab_report_paragraph(lmer_notes, "displays the annotated linear mixed-effects model results table")
  logs_paragraph <- tab_report_paragraph(logs_notes, "displays workflow log files for each analysis stage")

  c(
    "---",
    'title: "Report"',
    "format:",
    "  html:",
    "    css:",
    "      - assets/qpasst.css",
    "page-layout: full",
    "body-classes: qpasst-report-page",
    "---",
    "",
    '<div class="qpasst-report-viewport">',
    '<div class="qpasst-report-card">',
    '<div class="qpasst-report-sheet">',
    sprintf("<h2>%s Report</h2>", html_escape(project_name)),
    html_paragraph(overview),
    html_section("Data", data_paragraph),
    html_section("ENmix QC", enmix_paragraph),
    html_section("Quality Control", quality_control_paragraph),
    html_section("Batch Effect", batch_effect_paragraph),
    html_section("Metrics", metrics_paragraph),
    html_section("GLM", glm_paragraph),
    html_section("LMER", lmer_paragraph),
    html_section("Logs", logs_paragraph),
    "</div>",
    "</div>",
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
    '## {.sidebar width="320px"}',
    "### Notes",
    "",
    bullet_lines,
    ""
  )
}

compose_page <- function(title, notes, body_lines, body_classes = NULL) {
  front_matter <- c(
    "---",
    sprintf('title: "%s"', title),
    "format: dashboard"
  )

  if (!is.null(body_classes) && nzchar(body_classes)) {
    front_matter <- c(front_matter, sprintf("body-classes: %s", body_classes))
  }

  c(
    front_matter,
    "---",
    "",
    sidebar_lines(notes),
    "## Column",
    "",
    body_lines
  )
}

compose_plain_page <- function(title, body_lines, body_classes = NULL) {
  front_matter <- c(
    "---",
    sprintf('title: "%s"', title),
    "format: dashboard"
  )

  if (!is.null(body_classes) && nzchar(body_classes)) {
    front_matter <- c(front_matter, sprintf("body-classes: %s", body_classes))
  }

  c(
    front_matter,
    "---",
    "",
    "## Column",
    "",
    body_lines
  )
}

inject_card_headers <- function(html_path, headers) {
  if (!file.exists(html_path)) {
    return(invisible(FALSE))
  }

  lines <- readLines(html_path, warn = FALSE, encoding = "UTF-8")
  header_markup <- sprintf('<div class="card-header">%s</div>', headers)

  if (all(header_markup %in% lines)) {
    return(invisible(TRUE))
  }

  card_idx <- grep('^<div class=\"card(?: cell)? bslib-card', lines, perl = TRUE)
  if (length(card_idx) == 0L) {
    return(invisible(FALSE))
  }

  if (length(card_idx) < length(headers)) {
    return(invisible(FALSE))
  }

  offset <- 0L
  for (i in seq_along(headers)) {
    lines <- append(lines, header_markup[[i]], after = card_idx[[i]] + offset)
    offset <- offset + 1L
  }

  write_utf8(html_path, lines)
  invisible(TRUE)
}

post_process_dashboard_titles <- function(
    project_dir,
    enmix_titles = character(),
    metrics_titles = character(),
    quality_control_titles = character(),
    batch_effect_titles = character(),
    glm_table_title = "Table 1. Generalised Linear Model Results and Genomic Annotation of CpG Sites by Phenotype(s)",
    lmer_table_title = "Table 1. Linear Mixed-Effects Model Results and Genomic Annotation of CpG Sites by Phenotype(s) and Timepoint"
) {
rewrite_logo_links <- function(project_dir, href = "./index.html") {
  docs_dir <- file.path(project_dir, "docs")
  if (!dir.exists(docs_dir)) {
    return(invisible(FALSE))
  }

  html_files <- list.files(docs_dir, pattern = "\\.html$", full.names = TRUE)
  if (!length(html_files)) {
    return(invisible(FALSE))
  }

  pattern <- '<a href="\\./index\\.html" class="navbar-brand navbar-brand-logo">'
  replacement <- sprintf('<a href="%s" class="navbar-brand navbar-brand-logo">', href)

  for (html_path in html_files) {
    lines <- readLines(html_path, warn = FALSE, encoding = "UTF-8")
    lines <- gsub(pattern, replacement, lines)
    write_utf8(html_path, lines)
  }

  invisible(TRUE)
}

  docs_dir <- file.path(project_dir, "docs")

  figure_headers <- function(titles) {
    titles <- titles[nzchar(titles)]
    if (!length(titles)) {
      return(character())
    }
    titles <- sub("^Figure\\s+[0-9]+:\\s*", "", titles)
    sprintf("Figure %s: %s", seq_along(titles), titles)
  }

  inject_card_headers(
    file.path(docs_dir, "enmix-qc.html"),
    figure_headers(enmix_titles)
  )

  inject_card_headers(
    file.path(docs_dir, "metrics.html"),
    figure_headers(metrics_titles)
  )

  inject_card_headers(
    file.path(docs_dir, "quality-control.html"),
    c(
      figure_headers(quality_control_titles),
      "Table 1: CpG Detection Summary",
      "Table 2: Sample Detection Summary"
    )
  )

  inject_card_headers(
    file.path(docs_dir, "batch-effect.html"),
    figure_headers(batch_effect_titles)
  )

  inject_card_headers(
    file.path(docs_dir, "glm.html"),
    glm_table_title
  )

  inject_card_headers(
    file.path(docs_dir, "lmer.html"),
    lmer_table_title
  )

  inject_card_headers(
    file.path(docs_dir, "logs.html"),
    c(
      "Methylation Analysis",
      "Data Preparation",
      "Batch Effect",
      "GLM Analysis",
      "LMER Analysis"
    )
  )
  rewrite_logo_links(project_dir, href = "./index.html")
}

prepared_report <- prepareDnamReportInputs(
  outputDir = project_dir,
  qcDir = qc_dir,
  preprocessingDir = preprocessing_dir,
  postprocessingDir = postprocessing_dir,
  svaDir = sva_dir,
  glmDir = glm_dir,
  glmmDir = glmm_dir,
  figDir = fig_dir,
  verbose = FALSE,
  logs = FALSE,
  logDir = logs_dir
)

data_path_for_qmd <- to_dashboard_data_path(pheno_file)

enmix_items <- copy_figure_assets(
  qc_dir,
  asset_subdir = "enmix-qc"
)
metrics_items <- copy_figure_assets(
  postprocessing_dir,
  asset_subdir = "metrics"
)
qc_items <- copy_figure_assets(
  preprocessing_dir,
  asset_subdir = "quality-control"
)
batch_items <- copy_figure_assets(
  sva_dir,
  asset_subdir = "batch-effects"
)

log_assets <- list(
  methylation = copy_log_asset(
    file.path(logs_dir, "log_preprocessingMinfiEwasWater.txt"),
    "preprocessingMinfiEwasWater.txt"
  ),
  data = copy_log_asset(
    file.path(logs_dir, "log_preprocessingPheno.txt"),
    "preprocessingPheno.txt"
  ),
  batch = copy_log_asset(
    file.path(logs_dir, "log_svaEnmix.txt"),
    "svaEnmix.txt"
  ),
  glm = copy_log_asset(
    file.path(logs_dir, "log_methylationGLM_T1.txt"),
    "methylationGLM_T1.txt"
  ),
  lmer = copy_log_asset(
    file.path(logs_dir, "log_methylationGLMM_T1T2.txt"),
    "methylationGLMM_T1T2.txt"
  )
)

data_summary <- summarize_dataset(pheno_file)
detection_tables <- read_detection_tables(
  detp_path = detp_path,
  threshold = detPThreshold,
  cpg_path = cpg_detection_path,
  sample_path = sample_detection_path
)
cpg_detection_summary <- summarize_cpg_detection(detection_tables$cpg)
sample_detection_summary <- summarize_sample_detection(detection_tables$sample)
sva_summary <- summarize_sva(
  pheno_path = pheno_file,
  log_path = file.path(logs_dir, "log_svaEnmix.txt")
)
log_summary <- summarize_logs(log_assets)

enmix_figure_titles <- enmix_items$title
enmix_notes <- c(
  paste(
    "`ENmix` produces control plots similar to those generated by Illumina's GenomeStudio software.",
    "Infinium controls are not evaluated using absolute intensity values.",
    "Instead, we interpret them based on their expected signal range across the dataset.",
    "This approach is intended to provide greater flexibility and reliability when accounting for natural biological variation."
  )
)
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
    "detection_pvalues",
    "quality_control",
    "sexClinical",
    "sexComparison_RawNorm",
    "sexPrediction"
  ),
  title = c(
    "Figure 1: Beta Value Density Distribution (MSET)",
    "Figure 2: Detection P-values (RGSET)",
    "Figure 3: Quality Control Plot (MSET)",
    "Figure 4: Clinical Sex Distribution (GSET)",
    "Figure 5: Sex Comparison Before and After Normalisation (MSETF)",
    "Figure 6: Sex Prediction Plot (GSET)"
  ),
  description = c(
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
  item_idx <- grep(
    quality_control_metadata$pattern[[metadata_idx]],
    qc_items$original_name,
    ignore.case = TRUE
  )
  if (length(item_idx)) {
    quality_control_figure_titles[item_idx] <- quality_control_metadata$title[[metadata_idx]]
    quality_control_figure_descriptions[item_idx] <- quality_control_metadata$description[[metadata_idx]]
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
glm_table_title <- "Table 1. Generalised Linear Model Results and Genomic Annotation of CpG Sites by Phenotype(s)"
glm_table_description <- paste(
  "`Table 1` presents CpG methylation probes analysed using a generalised linear model with the `glm2` package to assess their association with the selected phenotype(s).",
  "Each row represents one CpG site, with statistical results and genomic annotation information.",
  "The IlmnID column gives the unique Illumina probe identifier, while Name provides the probe label used in the analysis or annotation file.",
  "The P.Value column reports the statistical evidence for association between methylation at each CpG site and the phenotype(s), where smaller values indicate stronger evidence before multiple-testing correction.",
  "The chr and pos columns indicate the chromosome and genomic position of the probe.",
  "The UCSC_RefGene_Group and UCSC_RefGene_Name columns describe the gene-related annotation and associated gene name, when available.",
  "The Relation_to_Island column indicates whether the CpG site is located in an OpenSea, Shore, Shelf, or CpG island region.",
  "The GencodeV41_Group column provides additional gene-region annotation based on GENCODE version 41, including transcript or transcription start site information where applicable."
)
lmer_table_title <- "Table 1. Linear Mixed-Effects Model Results and Genomic Annotation of CpG Sites by Phenotype(s) and Timepoint"
lmer_table_description <- paste(
  "`Table 1` presents CpG methylation probes analysed using a linear mixed-effects model with the `lmer` package to assess phenotype-related methylation changes across timepoints.",
  "Each row represents one CpG site, with model results and genomic annotation information.",
  "The IlmnID column gives the unique Illumina probe identifier, while Name provides the probe label used in the analysis or annotation file.",
  "The phenotype:Timepoint P.Value column reports the statistical evidence for an interaction between the selected phenotype(s) and timepoint, where smaller values suggest stronger evidence of time-dependent methylation differences.",
  "The chr and pos columns indicate the chromosome and genomic position of the probe.",
  "The UCSC_RefGene_Group and UCSC_RefGene_Name columns describe the gene-related annotation and associated gene name, when available.",
  "The Relation_to_Island column indicates whether the CpG site is located in an OpenSea, Shore, Shelf, or CpG island region.",
  "The GencodeV41_Group column provides additional gene-region annotation based on GENCODE version 41."
)

data_notes <- make_data_notes(data_summary)
metrics_notes <- make_metrics_notes(
  metrics_items,
  metrics_figure_titles,
  metrics_figure_descriptions,
  data_summary
)
quality_control_notes <- make_quality_control_notes(
  qc_items,
  quality_control_figure_titles,
  quality_control_figure_descriptions,
  cpg_detection_summary,
  sample_detection_summary
)
batch_effect_notes <- make_batch_effect_notes(
  batch_items,
  batch_effect_figure_titles,
  batch_effect_figure_descriptions,
  sva_summary
)
glm_notes <- glm_table_description
lmer_notes <- lmer_table_description
logs_notes <- make_logs_notes(log_summary)

if (file.exists(logo_source_path)) {
  dir.create(assets_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(logo_source_path, file.path(assets_dir, "dnaEPICO.svg"), overwrite = TRUE)
  if (!ok) {
    warning("Failed to copy navbar logo: ", logo_source_path)
  }
}

quarto_yml <- c(
  "project:",
  "  type: website",
  "  output-dir: docs",
  "  resources:",
  "    - assets/",
  "",
  "website:",
  sprintf('  title: "%s"', project_name),
  "  navbar:",
  "    logo: assets/dnaEPICO.svg",
  "    left:",
  "      - href: index.qmd",
  '        text: "Data"',
  "      - href: enmix-qc.qmd",
  '        text: "ENmix QC"',
  "      - href: quality-control.qmd",
  '        text: "Quality Control"',
  "      - href: batch-effect.qmd",
  '        text: "Batch Effect"',
  "      - href: metrics.qmd",
  '        text: "Metrics"',
  "      - href: glm.qmd",
  '        text: "GLM Analysis"',
  "      - href: lmer.qmd",
  '        text: "LMER Analysis"',
  "      - href: report.qmd",
  '        text: "Report"',
  "      - href: logs.qmd",
  '        text: "Logs"',
  "    search: true",
  "  page-navigation: false",
  "  bread-crumbs: false",
  "  reader-mode: false",
  "",
  "format:",
  "  dashboard:",
  "    theme: cosmo",
  "    css:",
  "      - assets/qpasst.css",
  "    toc: false",
  "",
  "execute:",
  "  warning: false",
  "  message: false"
)

site_css <- c(
  ":root {",
  "  --qpasst-border: #e5e7eb;",
  "  --qpasst-bg: #d5d9df;",
  "  --qpasst-navbar-bg: #f4f6f8;",
  "  --qpasst-card-bg: #ffffff;",
  "  --qpasst-accent: #0d6efd;",
  "}",
  "",
  "body {",
  "  background: #ffffff;",
  "}",
  "",
  ".navbar,",
  ".quarto-navbar,",
  "#quarto-header .navbar {",
  "  background: var(--qpasst-navbar-bg) !important;",
  "  border-bottom: 1px solid var(--qpasst-border) !important;",
  "  box-shadow: none !important;",
  "}",
  "",
  ".navbar-brand,",
  ".navbar-title {",
  "  color: #3f444a !important;",
  "  font-weight: 600;",
  "}",
  "",
  "#quarto-header .navbar-brand-container {",
  "  display: flex;",
  "  align-items: center;",
  "}",
  "",
  "#quarto-header .navbar-brand.navbar-brand-logo {",
  "  display: flex !important;",
  "  align-items: center;",
  "  padding-top: 0;",
  "  padding-bottom: 0;",
  "  margin-right: 0.35rem;",
  "}",
  "",
  "#quarto-header .navbar-brand.navbar-brand-logo img,",
  "#quarto-header .navbar-brand.navbar-brand-logo .navbar-logo,",
  "#quarto-header .navbar-brand.navbar-brand-logo .navbar-logo-image {",
  "  display: block;",
  "  height: 48px !important;",
  "  width: auto;",
  "  max-height: 48px !important;",
  "  max-width: 48px !important;",
  "  object-fit: contain;",
  "}",
  "",
  "#quarto-header .navbar-brand:not(.navbar-brand-logo) {",
  "  display: none !important;",
  "}",
  "",
  ".navbar-nav .nav-link {",
  "  color: #495057 !important;",
  "  font-weight: 700;",
  "}",
  "",
  ".navbar-nav .nav-link.active,",
  ".navbar-nav .show > .nav-link,",
  ".navbar-nav .nav-link:hover {",
  "  color: var(--qpasst-accent) !important;",
  "}",
  "",
  "#quarto-dashboard-header {",
  "  display: none !important;",
  "}",
  "",
  ".quarto-dashboard .bslib-sidebar-layout > .main,",
  ".quarto-dashboard .bslib-sidebar-layout .sidebar-content {",
  "  overflow-y: auto !important;",
  "}",
  "",
  ".quarto-dashboard .sidebar-content > .bslib-grid,",
  ".quarto-dashboard .sidebar-content > .bslib-grid > .bslib-grid,",
  ".quarto-dashboard .sidebar-content > .bslib-grid > .bslib-grid > .bslib-grid {",
  "  grid-auto-rows: max-content !important;",
  "}",
  "",
  ".quarto-dashboard .bslib-card {",
  "  min-height: 0 !important;",
  "  margin-bottom: 1.5rem;",
  "}",
  "",
  ".quarto-dashboard .bslib-card .card-body {",
  "  height: auto !important;",
  "  overflow: visible !important;",
  "  padding: 1.5rem !important;",
  "}",
  "",
  "body.qpasst-logs-page .sidebar-content > .bslib-grid > .bslib-grid {",
  "  grid-template-rows: minmax(3em, max-content) minmax(3em, max-content) minmax(3em, max-content) !important;",
  "}",
  "",
  "body.qpasst-logs-page .panel-tabset.bslib-grid {",
  "  grid-template-columns: minmax(3em, 1fr) !important;",
  "  grid-auto-rows: max-content !important;",
  "}",
  "",
  "body.qpasst-logs-page .bslib-card {",
  "  min-height: calc(100vh - 12rem) !important;",
  "}",
  "",
  "body.qpasst-logs-page .bslib-card .card-body {",
  "  overflow: auto !important;",
  "}",
  "",
  "body.qpasst-logs-page pre.text {",
  "  margin: 0;",
  "  max-height: calc(100vh - 16rem);",
  "  overflow: auto;",
  "  white-space: pre;",
  "}",
  "",
  "body.qpasst-report-page #title-block-header {",
  "  display: none !important;",
  "}",
  "",
  "body.qpasst-report-page #quarto-content,",
  "body.qpasst-report-page main.content,",
  "body.qpasst-report-page #quarto-document-content,",
  "body.qpasst-report-page .page-columns,",
  "body.qpasst-report-page .page-full {",
  "  width: 100% !important;",
  "  max-width: none !important;",
  "}",
  "",
  "body.qpasst-report-page main.content {",
  "  margin: 0 !important;",
  "  padding: 0.75rem !important;",
  "}",
  "",
  "body.qpasst-report-page .qpasst-report-viewport {",
  "  width: 100%;",
  "  min-height: calc(100vh - 7rem);",
  "}",
  "",
  "body.qpasst-report-page .qpasst-report-card {",
  "  width: 100%;",
  "  min-height: calc(100vh - 7rem);",
  "  border: 1px solid var(--qpasst-border);",
  "  border-radius: 12px;",
  "  background: var(--qpasst-card-bg);",
  "  box-shadow: 0 4px 14px rgba(15, 23, 42, 0.05);",
  "  overflow-y: auto;",
  "  padding: 1.25rem 1.5rem;",
  "}",
  "",
  "body.qpasst-report-page .qpasst-report-sheet {",
  "  width: 100%;",
  "  min-height: 297mm;",
  "  margin: 0;",
  "  padding: 0.5rem 0 2rem;",
  "}",
  "",
  "body.qpasst-report-page .qpasst-report-sheet h2,",
  "body.qpasst-report-page .qpasst-report-sheet h3 {",
  "  color: #2f3540;",
  "  margin-top: 0;",
  "}",
  "",
  "body.qpasst-report-page .qpasst-report-sheet p {",
  "  font-size: 1.02rem;",
  "  line-height: 1.75;",
  "  margin-bottom: 1rem;",
  "}",
  "",
  ".dashboard .sidebar,",
  ".quarto-dashboard .sidebar,",
  ".sidebar {",
  "  background: var(--qpasst-bg) !important;",
  "  border-right: 1px solid var(--qpasst-border) !important;",
  "  overflow-x: hidden !important;",
  "}",
  "",
  ".sidebar h3,",
  ".sidebar .h3 {",
  "  color: #2f3540 !important;",
  "}",
  "",
  ".sidebar p,",
  ".sidebar li,",
  ".sidebar code {",
  "  font-size: 1.06rem;",
  "  line-height: 1.7;",
  "  overflow-wrap: anywhere;",
  "  word-break: break-word;",
  "}",
  "",
  "body.qpasst-data-page .sidebar li {",
  "  font-size: 1.08rem;",
  "  line-height: 1.65;",
  "}",
  "",
  ".sidebar li::marker {",
  "  color: var(--qpasst-accent);",
  "}",
  "",
  ".qpasst-gallery {",
  "  display: grid;",
  "  grid-template-columns: 1fr;",
  "  gap: 1.25rem;",
  "  margin-top: 0.5rem;",
  "}",
  "",
  ".qpasst-gallery-item {",
  "  margin: 0;",
  "  border: 1px solid var(--qpasst-border);",
  "  border-radius: 12px;",
  "  background: var(--qpasst-card-bg);",
  "  padding: 1rem;",
  "  box-shadow: 0 4px 14px rgba(15, 23, 42, 0.05);",
  "}",
  "",
  ".qpasst-gallery-item figcaption {",
  "  margin-bottom: 0.75rem;",
  "  line-height: 1.5;",
  "}",
  "",
  ".qpasst-gallery-item img {",
  "  display: block;",
  "  width: 50%;",
  "  height: auto;",
  "  border-radius: 8px;",
  "  background: #fff;",
  "}",
  "",
  ".qpasst-gallery-item code {",
  "  white-space: normal;",
  "  word-break: break-word;",
  "}",
  "",
  ".qpasst-figure-card {",
  "  margin: 0;",
  "}",
  "",
  ".qpasst-figure-card img {",
  "  display: block;",
  "  width: 100%;",
  "  height: auto;",
  "  border-radius: 8px;",
  "  background: #fff;",
  "}",
  "",
  ".qpasst-figure-card figcaption {",
  "  margin-top: 0.75rem;",
  "}",
  "",
  ".qpasst-download {",
  "  margin: 0;",
  "  padding: 0.75rem 0.9rem;",
  "  border-radius: 8px;",
  "  background: #f8fafc;",
  "  border: 1px dashed var(--qpasst-border);",
  "}",
  "",
  ".table {",
  "  font-size: 0.92rem;",
  "}",
  "",
  "pre code {",
  "  white-space: pre-wrap;",
  "  word-break: break-word;",
  "}",
  "",
  "@media (max-width: 768px) {",
  "  .qpasst-gallery {",
  "    grid-template-columns: 1fr;",
  "  }",
  "}"
)

data_page <- compose_page(
  title = "Data",
  notes = data_notes,
  body_classes = "qpasst-data-page",
  body_lines = c(
    "### Data Preview",
    "",
    "```{r}",
    "#| echo: false",
    "#| include: false",
    sprintf("data_path <- %s", r_string(data_path_for_qmd)),
    sprintf("data_path_display <- %s", r_string(if (nzchar(data_path_for_qmd)) data_path_for_qmd else "Not detected")),
    sprintf("participant_col <- %s", r_string(if (is.null(data_summary$participant_col)) "" else data_summary$participant_col)),
    sprintf("timepoint_col <- %s", r_string(if (is.null(data_summary$timepoint_col)) "" else data_summary$timepoint_col)),
    "qp_data <- NULL",
    "qp_data_error <- NULL",
    "if (!nzchar(data_path)) {",
    "  qp_data_error <- 'No data file path is configured for this dashboard.'",
    "} else if (!file.exists(data_path)) {",
    "  qp_data_error <- paste0('Data file not found at `', data_path_display, '`.')",
    "} else {",
    "  qp_data <- tryCatch(",
    "    utils::read.csv(data_path, check.names = FALSE),",
    "    error = function(e) {",
    "      qp_data_error <<- conditionMessage(e)",
    "      NULL",
    "    }",
    "  )",
    "  if (!is.null(qp_data)) {",
    "    front_cols <- c(participant_col, timepoint_col)",
    "    front_cols <- front_cols[nzchar(front_cols) & front_cols %in% names(qp_data)]",
    "    qp_data <- qp_data[, c(front_cols, setdiff(names(qp_data), front_cols)), drop = FALSE]",
    "  }",
    "}",
    "```",
    "",
    "```{r}",
    "#| echo: false",
    "#| results: asis",
    "if (is.null(qp_data)) {",
    "  cat('::: {.callout-warning}\\n')",
    "  cat(qp_data_error, '\\n\\n')",
    "  cat('Update the configured CSV path in the `dnamReport()` call or place the file at the expected location before rendering.\\n')",
    "  cat(':::\\n')",
    "}",
    "```",
    "",
    "```{r}",
    "#| echo: false",
    "if (!is.null(qp_data)) {",
    "  if (requireNamespace('DT', quietly = TRUE)) {",
    "    DT::datatable(",
    "      qp_data,",
    "      filter = 'top',",
    "      options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),",
    "      rownames = FALSE",
    "    )",
    "  } else {",
    "    knitr::kable(",
    "      utils::head(qp_data, 25),",
    "      format = 'html',",
    "      table.attr = \"class='table table-striped table-sm'\"",
    "    )",
    "  }",
    "}",
    "```"
  )
)

enmix_page <- compose_page(
  title = "ENmix QC",
  notes = enmix_notes,
  body_lines = c(
    build_figure_cards(
      enmix_items,
      "No supported image files were found for the ENmix QC tab.",
      title_prefix = "Figure",
      show_filename = FALSE,
      figure_titles = enmix_figure_titles
    )
  )
)

metrics_page <- compose_page(
  title = "Metrics",
  notes = metrics_notes,
  body_lines = c(
    build_figure_cards(
      metrics_items,
      "No supported image files were found for the Metrics tab.",
      title_prefix = "Figure",
      show_filename = FALSE,
      figure_titles = metrics_figure_titles
    )
  )
)

quality_control_page <- compose_page(
  title = "Quality Control",
  notes = quality_control_notes,
  body_lines = c(
    build_figure_cards(
      qc_items,
      "No supported image files were found for the Quality Control tab.",
      title_prefix = "Figure",
      show_filename = FALSE,
      figure_titles = quality_control_figure_titles
    ),
    build_data_frame_table_section(
      title = "Table 1: CpG Detection Summary",
      data = detection_tables$cpg,
      empty_message = detection_table_warning(detection_tables),
      preview_rows = 10L
    ),
    build_data_frame_table_section(
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
  body_lines = c(
    build_figure_cards(
      batch_items,
      "No supported image files were found for the Batch Effect tab.",
      title_prefix = "Figure",
      show_filename = FALSE,
      figure_titles = batch_effect_figure_titles
    )
  )
)

glm_page <- compose_page(
  title = "GLM Analysis",
  notes = glm_notes,
  body_lines = c(
    build_xlsx_table_section(
      title = glm_table_title,
      data_path = glm_table_path,
      sheet = "annotatedGLM",
      page_length = 10L,
      preview_rows = 25L,
      var_prefix = "glm_results",
      interactive = TRUE
    )
  )
)

lmer_page <- compose_page(
  title = "LMER Analysis",
  notes = lmer_notes,
  body_lines = c(
    build_xlsx_table_section(
      title = lmer_table_title,
      data_path = lmer_table_path,
      sheet = "annotatedLME",
      page_length = 10L,
      preview_rows = 25L,
      var_prefix = "lmer_results",
      interactive = TRUE
    )
  )
)

logs_page <- compose_page(
  title = "Logs",
  notes = logs_notes,
  body_classes = "qpasst-logs-page",
  body_lines = c(
    "```{r}",
    "#| echo: false",
    "#| include: false",
    "render_log_block <- function(path, label) {",
    "  if (!file.exists(path)) {",
    "    cat('::: {.callout-warning}\\n')",
    "    cat(sprintf('Log file not found for `%s`.\\n', label))",
    "    cat(':::\\n')",
    "    return(invisible(NULL))",
    "  }",
    "  lines <- readLines(path, warn = FALSE, encoding = 'UTF-8')",
    "  if (!length(lines)) {",
    "    cat('::: {.callout-note}\\n')",
    "    cat(sprintf('`%s` is empty.\\n', label))",
    "    cat(':::\\n')",
    "    return(invisible(NULL))",
    "  }",
    "  cat('```text\\n')",
    "  cat(paste(lines, collapse = '\\n'))",
    "  cat('\\n```\\n')",
    "}",
    "```",
    "",
    "### Methylation Analysis",
    "",
    "Displays the methylation preprocessing log, including IDAT loading, normalisation, filtering, and cell composition estimation steps.",
    "",
    "```{r}",
    "#| echo: false",
    "#| results: asis",
    sprintf("render_log_block(%s, 'Methylation Analysis')", r_string(log_assets$methylation$asset_path)),
    "```",
    "",
    "### Data Preparation",
    "",
    "Displays the phenotype preparation log, including timepoint splitting and methylation matrix export steps.",
    "",
    "```{r}",
    "#| echo: false",
    "#| results: asis",
    sprintf("render_log_block(%s, 'Data Preparation')", r_string(log_assets$data$asset_path)),
    "```",
    "",
    "### Batch Effect",
    "",
    "Displays the hidden-effect and surrogate-variable analysis log used for batch-effect assessment.",
    "",
    "```{r}",
    "#| echo: false",
    "#| results: asis",
    sprintf("render_log_block(%s, 'Batch Effect')", r_string(log_assets$batch$asset_path)),
    "```",
    "",
    "### GLM Analysis",
    "",
    "Displays the generalised linear model log, including phenotype association testing and CpG annotation steps.",
    "",
    "```{r}",
    "#| echo: false",
    "#| results: asis",
    sprintf("render_log_block(%s, 'GLM Analysis')", r_string(log_assets$glm$asset_path)),
    "```",
    "",
    "### LMER Analysis",
    "",
    "Displays the linear mixed-effects model log, including timepoint interaction testing and CpG annotation steps.",
    "",
    "```{r}",
    "#| echo: false",
    "#| results: asis",
    sprintf("render_log_block(%s, 'LMER Analysis')", r_string(log_assets$lmer$asset_path)),
    "```"
  )
)

report_page <- build_report_page(
  project_name = project_name,
  project_dir = project_dir,
  data_notes = data_notes,
  enmix_notes = enmix_notes,
  quality_control_notes = quality_control_notes,
  batch_effect_notes = batch_effect_notes,
  metrics_notes = metrics_notes,
  glm_notes = glm_notes,
  lmer_notes = lmer_notes,
  logs_notes = logs_notes
)

unlink(
  c(
    file.path(project_dir, c("data.qmd", "DNAm.html")),
    file.path(project_dir, "docs", "data.html")
  ),
  force = TRUE
)

write_utf8(file.path(project_dir, "_quarto.yml"), quarto_yml)
write_utf8(file.path(assets_dir, "qpasst.css"), site_css)
write_utf8(file.path(project_dir, "index.qmd"), data_page)
write_utf8(file.path(project_dir, "enmix-qc.qmd"), enmix_page)
write_utf8(file.path(project_dir, "metrics.qmd"), metrics_page)
write_utf8(file.path(project_dir, "quality-control.qmd"), quality_control_page)
write_utf8(file.path(project_dir, "batch-effect.qmd"), batch_effect_page)
write_utf8(file.path(project_dir, "glm.qmd"), glm_page)
write_utf8(file.path(project_dir, "lmer.qmd"), lmer_page)
write_utf8(file.path(project_dir, "report.qmd"), report_page)
write_utf8(file.path(project_dir, "logs.qmd"), logs_page)

quarto_bin <- find_quarto()

source_files <- file.path(
  project_dir,
  c(
    "_quarto.yml",
    "index.qmd",
    "enmix-qc.qmd",
    "metrics.qmd",
    "quality-control.qmd",
    "batch-effect.qmd",
    "glm.qmd",
    "lmer.qmd",
    "report.qmd",
    "logs.qmd"
  )
)

emitLogMinfiEwasWater(
  c(
    paste("Generated", length(source_files), "report source files in:", project_dir)
  ),
  verbose = FALSE,
  log_path = log_path
)

render_status <- "skipped"
rendered_file <- file.path(project_dir, "docs", "index.html")
error_message <- NULL

if (nzchar(quarto_bin)) {
  emitLogMinfiEwasWater(
    c(
      "Rendering report site..."
    ),
    verbose = FALSE,
    log_path = log_path
  )
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(project_dir)
  render_output <- system2(quarto_bin, "render", stdout = TRUE, stderr = TRUE)
  status <- attr(render_output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (length(render_output)) {
    emitLogMinfiEwasWater(
      c("Quarto output:", render_output),
      verbose = FALSE,
      log_path = log_path
    )
  }
  if (identical(status, 0L)) {
    render_status <- "rendered"
    post_process_dashboard_titles(
      project_dir,
      enmix_titles = enmix_figure_titles,
      metrics_titles = metrics_figure_titles,
      quality_control_titles = quality_control_figure_titles,
      batch_effect_titles = batch_effect_figure_titles,
      glm_table_title = glm_table_title,
      lmer_table_title = lmer_table_title
    )
    emitLogMinfiEwasWater(
      c(
        "Render complete."
      ),
      verbose = FALSE,
      log_path = log_path
    )
  } else {
    render_status <- "failed"
    error_message <- paste("Quarto render failed with status", status)
    emitLogMinfiEwasWater(
      error_message,
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
  emitLogMinfiEwasWater(
    c(error_message, project_dir),
    verbose = verbose,
    log_path = log_path
  )
}

if (!magick_available && any(grepl("\\.tiff?$", c(enmix_items$original_name, metrics_items$original_name, qc_items$original_name, batch_items$original_name), ignore.case = TRUE))) {
  emitLogMinfiEwasWater(
    c(
      "Note: TIFF figures were copied as-is because the `magick` package is not installed.",
      "Install `magick` and rerun this function if you want automatic PNG conversion for browser previews."
    ),
    verbose = FALSE,
    log_path = log_path
  )
}

emitLogMinfiEwasWater(
  paste("Report path:", normalizePathDnamReport(file.path(project_dir, "docs", "index.html"))),
  verbose = verbose,
  log_path = log_path
)

render_result <- structure(
  list(
    preparedReport = prepared_report,
    status = render_status,
    renderedFile = if (identical(render_status, "rendered")) rendered_file else NULL,
    errorMessage = error_message,
    logFile = log_path
  ),
  class = "dnaEPICO_dnamReport_render"
)

structure(
  list(
    preparedReport = prepared_report,
    renderResult = render_result,
    status = render_status,
    outputFile = if (identical(render_status, "rendered")) {
      rendered_file
    } else {
      normalizePathDnamReport(file.path(project_dir, "docs", "index.html"))
    },
    projectDir = project_dir,
    sourceFiles = source_files,
    docsDir = file.path(project_dir, "docs"),
    logoPath = logo_source_path,
    errorMessage = error_message,
    logFile = log_path
  ),
  class = "dnaEPICO_dnamReport"
)
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
  if (
    is.null(report_path) ||
      !length(report_path) ||
      is.na(report_path[[1]]) ||
      !nzchar(report_path[[1]])
  ) {
    report_path <- x$renderResult$renderedFile
  }

  cat("Class type: ", class(x)[[1]], "\n", sep = "")
  cat("Log output path: ", formatPrintPathDnamReport(x$logFile), "\n", sep = "")
  cat("Report output path: ", formatPrintPathDnamReport(report_path), "\n", sep = "")

  invisible(x)
}
