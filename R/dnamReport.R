#' Normalize a report path for cross-platform rendering
#'
#' @param path Character. Path to normalize.
#'
#' @return Character scalar containing a normalized path.
#'
#' @description
#' Internal helper that normalizes paths for `rmarkdown::render()` on Windows and
#' other platforms without requiring the path to exist yet.
#'
#' @keywords internal
#' @noRd
normalizePathDnamReport <- function(path) {
  gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = FALSE))
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
#' @description
#' Internal helper that inventories report-ready figure files for one report
#' section.
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
          recursive = recursive
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

#' Resolve the packaged DNA methylation report template
#'
#' @param templatePath Character or `NULL`. Optional custom template path.
#'
#' @return Character scalar containing the resolved template path.
#'
#' @description
#' Internal helper that resolves a custom report template path or the packaged
#' `dnamReport.Rmd` template.
#'
#' @keywords internal
#' @noRd
resolveTemplatePathDnamReport <- function(templatePath = NULL) {
  if (!is.null(templatePath) && nzchar(templatePath)) {
    return(normalizePathDnamReport(templatePath))
  }

  packaged_template <- system.file("extdata", "dnamReport.Rmd", package = "dnaEPICO")
  if (!nzchar(packaged_template)) {
    stop(
      "The packaged report template 'dnamReport.Rmd' was not found in inst/extdata.",
      call. = FALSE
    )
  }

  normalizePathDnamReport(packaged_template)
}

#' Prepare inputs for a DNA methylation report
#'
#' @param output Character. Name of the output PDF file.
#' @param outputDir Character. Directory where the report would be written.
#' @param qcDir Character. Directory containing ENmix QC figures, typically `.jpg`
#'   files.
#' @param preprocessingDir Character. Directory containing preprocessing QC
#'   figures, typically `.tiff` files.
#' @param postprocessingDir Character. Directory containing postprocessing metric
#'   figures, typically `.tiff` files.
#' @param svaDir Character. Directory containing SVA figures, typically `.tiff`
#'   files.
#' @param glmDir Character. Directory containing GLM figures, including optional
#'   phenotype subdirectories.
#' @param glmmDir Character. Directory containing GLMM figures, including optional
#'   phenotype subdirectories.
#' @param figDir Character. Directory used by `knitr` for temporary report figure
#'   output while rendering.
#' @param reportTitle Character. Title shown on the rendered report.
#' @param author Character. Author name shown on the rendered report.
#' @param date Character. Date string shown on the rendered report.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`. The
#'   default is `FALSE`.
#' @param logs Logical. If `TRUE`, write the same progress messages to
#'   `file.path(logDir, "log_dnamReport.txt")`.
#' @param logDir Character. Directory for optional log files.
#' @param templatePath Character or `NULL`. Optional path to a custom report
#'   template. When `NULL`, the packaged `dnamReport.Rmd` template is used.
#'
#' @return A list with class `"dnaEPICO_dnamReport_prepared"`.
#' \describe{
#'   \item{output}{Requested output file name.}
#'   \item{outputDir}{Normalized directory where the report would be written.}
#'   \item{outputFile}{Normalized full path to the intended output file.}
#'   \item{figDir}{Normalized directory used by the template for figure copies.}
#'   \item{templatePath}{Resolved path to the packaged or custom report
#'   template.}
#'   \item{params}{Named list of parameters passed to `rmarkdown::render()`.}
#'   \item{figureInventory}{Named list describing the available figures for each
#'   report section.}
#'   \item{missingFigureDirectories}{Character vector of expected figure
#'   directories that were not present at preparation time.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#' See [dnaEPICO_dnamReport_prepared-class] for a class-level overview.
#'
#' @description
#' `prepareDnamReportInputs()` resolves the report template, normalizes the input
#' and output paths, inventories the available figures for each report section,
#' and returns a structured object that can later be passed to
#' `renderDnamReport()`.
#'
#' @examples
#' ex <- dnaEPICO:::exampleDnamReportStateDnaEpico()
#' prepared_report <- prepareDnamReportInputs(
#'   output = "DNAm_Report.pdf",
#'   outputDir = ex$tempDir,
#'   qcDir = file.path(ex$tempDir, "qc"),
#'   preprocessingDir = file.path(ex$tempDir, "preprocessing"),
#'   postprocessingDir = file.path(ex$tempDir, "postprocessing"),
#'   svaDir = file.path(ex$tempDir, "sva"),
#'   glmDir = file.path(ex$tempDir, "glm"),
#'   glmmDir = file.path(ex$tempDir, "glmm"),
#'   figDir = file.path(ex$tempDir, "report-figures"),
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' prepared_report$output
#'
#' @export
prepareDnamReportInputs <- function(
    output = "DNAm_Report.pdf",
    outputDir = "reports",
    qcDir = "figures/preprocessingMinfiEwasWater/enMix",
    preprocessingDir = "figures/preprocessingMinfiEwasWater/qc",
    postprocessingDir = "figures/preprocessingMinfiEwasWater/metrics",
    svaDir = "figures/svaEnmix/sva",
    glmDir = "figures/methylationGLM_T1",
    glmmDir = "figures/methylationGLMM_T1T2",
    figDir = "reports/figures",
    reportTitle = "DNA methylation",
    author = "School of Biomedical Sciences",
    date = format(Sys.Date(), "%B %d, %Y"),
    verbose = FALSE,
    logs = FALSE,
    logDir = outputDir,
    templatePath = NULL
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = logDir,
    log_file = "log_dnamReport.txt"
  )

  resolved_template <- resolveTemplatePathDnamReport(templatePath)
  normalized_output_dir <- normalizePathDnamReport(outputDir)
  normalized_fig_dir <- normalizePathDnamReport(figDir)
  normalized_output_file <- normalizePathDnamReport(file.path(outputDir, output))

  figure_inventory <- list(
    qc = collectFigureInventoryDnamReport(
      directory = qcDir,
      patterns = c("\\.jpg$", "\\.jpeg$", "\\.png$"),
      label = "ENmix QC"
    ),
    preprocessing = collectFigureInventoryDnamReport(
      directory = preprocessingDir,
      patterns = c("\\.tif$", "\\.tiff$", "\\.png$"),
      label = "Preprocessing"
    ),
    postprocessing = collectFigureInventoryDnamReport(
      directory = postprocessingDir,
      patterns = c("\\.tif$", "\\.tiff$", "\\.png$"),
      label = "Postprocessing"
    ),
    sva = collectFigureInventoryDnamReport(
      directory = svaDir,
      patterns = c("\\.tif$", "\\.tiff$", "\\.png$"),
      label = "SVA"
    ),
    glm = collectFigureInventoryDnamReport(
      directory = glmDir,
      patterns = c("\\.tif$", "\\.tiff$", "\\.png$"),
      label = "GLM"
    ),
    glmm = collectFigureInventoryDnamReport(
      directory = glmmDir,
      patterns = c("\\.tif$", "\\.tiff$", "\\.png$"),
      label = "GLMM"
    )
  )

  missing_directories <- vapply(
    figure_inventory,
    function(section) !isTRUE(section$exists),
    logical(1)
  )
  missing_directories <- names(missing_directories)[missing_directories]

  figure_counts <- vapply(
    figure_inventory,
    function(section) section$count,
    integer(1)
  )

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      paste("Prepared report template:       ", resolved_template),
      paste("Normalized output directory:    ", normalized_output_dir),
      paste("Normalized output file:         ", normalized_output_file),
      paste("Normalized figure directory:    ", normalized_fig_dir),
      paste("Figure inventory counts:", paste(names(figure_counts), figure_counts, sep = "=", collapse = "; ")),
      if (length(missing_directories) == 0L) {
        "Missing figure directories:      none"
      } else {
        paste("Missing figure directories:      ", paste(missing_directories, collapse = ", "))
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      output = output,
      outputDir = normalized_output_dir,
      outputFile = normalized_output_file,
      figDir = normalized_fig_dir,
      templatePath = resolved_template,
      params = list(
        reportTitle = reportTitle,
        author = author,
        date = date,
        qcDir = figure_inventory$qc$directory,
        preprocessingDir = figure_inventory$preprocessing$directory,
        postprocessingDir = figure_inventory$postprocessing$directory,
        svaDir = figure_inventory$sva$directory,
        glmDir = figure_inventory$glm$directory,
        glmmDir = figure_inventory$glmm$directory,
        figDir = normalized_fig_dir
      ),
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
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`.
#' @param logs Logical. If `TRUE`, write the same progress messages to a log file.
#' @param logDir Character or `NULL`. Directory for optional log files. When
#'   `NULL`, the directory from `preparedReport$logFile` is reused when available.
#' @param clean Logical. Passed to `rmarkdown::render()` to remove intermediate
#'   files after rendering.
#'
#' @return A list with class `"dnaEPICO_dnamReport_render"`.
#' \describe{
#'   \item{preparedReport}{The prepared input object supplied to
#'   `renderDnamReport()`.}
#'   \item{status}{Render status string such as `"rendered"`, `"skipped"`, or
#'   `"failed"`.}
#'   \item{renderedFile}{Normalized path to the rendered PDF file when rendering
#'   succeeded, otherwise `NULL`.}
#'   \item{errorMessage}{Render error or skip message when available, otherwise
#'   `NULL`.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#' See [dnaEPICO_dnamReport_render-class] for a class-level overview.
#'
#' @description
#' `renderDnamReport()` renders the packaged or custom DNA methylation report
#' template using a prepared report object. The function returns a structured
#' result instead of relying only on file-writing side effects, which makes it
#' easier to inspect rendering status in scripts and tests.
#'
#' @examples
#' ex <- dnaEPICO:::exampleDnamReportStateDnaEpico()
#' render_result <- renderDnamReport(
#'   preparedReport = ex$preparedReport,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' render_result$status
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

  if (is.null(logDir) || !nzchar(logDir)) {
    if (!is.null(preparedReport$logFile) && nzchar(preparedReport$logFile)) {
      logDir <- dirname(preparedReport$logFile)
    } else {
      logDir <- preparedReport$outputDir
    }
  }

  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = logDir,
    log_file = "log_dnamReport.txt"
  )

  dir.create(preparedReport$outputDir, recursive = TRUE, showWarnings = FALSE)
  dir.create(preparedReport$figDir, recursive = TRUE, showWarnings = FALSE)

  rendered_file <- NULL
  error_message <- NULL
  status <- "rendered"

  tiff_available <- requireNamespace("tiff", quietly = TRUE)
  tinytex_available <- requireNamespace("tinytex", quietly = TRUE) &&
    isTRUE(tinytex::is_tinytex())

  if (!tiff_available) {
    status <- "failed"
    error_message <- "Package 'tiff' is required to render report figure sections."
  } else if (!tinytex_available) {
    status <- "skipped"
    error_message <- "TinyTeX/LaTeX is not available; the report was not rendered."
  } else {
    render_attempt <- tryCatch(
      rmarkdown::render(
        input = preparedReport$templatePath,
        output_file = preparedReport$output,
        output_dir = preparedReport$outputDir,
        params = preparedReport$params,
        knit_root_dir = dirname(preparedReport$templatePath),
        quiet = !isTRUE(verbose),
        clean = clean,
        envir = new.env(parent = globalenv())
      ),
      error = function(error) error
    )

    if (inherits(render_attempt, "error")) {
      status <- "failed"
      error_message <- conditionMessage(render_attempt)
    } else {
      rendered_file <- normalizePathDnamReport(render_attempt)
    }
  }

  emitLogMinfiEwasWater(
    c(
      "=======================================================================",
      paste("Report rendering status:        ", status),
      if (is.null(rendered_file)) {
        paste("Report output path:             ", preparedReport$outputFile)
      } else {
        paste("Rendered report path:           ", rendered_file)
      },
      if (is.null(error_message)) {
        "Render message:                 none"
      } else {
        paste("Render message:                 ", error_message)
      },
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  structure(
    list(
      preparedReport = preparedReport,
      status = status,
      renderedFile = rendered_file,
      errorMessage = error_message,
      logFile = log_path
    ),
    class = "dnaEPICO_dnamReport_render"
  )
}

#' Prepare and optionally render a DNA methylation PDF report
#'
#' @param output Character. Name of the output PDF file.
#' @param outputDir Character. Directory where the report would be written.
#' @param qcDir Character. Directory containing ENmix QC figures, typically `.jpg`
#'   files.
#' @param preprocessingDir Character. Directory containing preprocessing QC
#'   figures, typically `.tiff` files.
#' @param postprocessingDir Character. Directory containing postprocessing metric
#'   figures, typically `.tiff` files.
#' @param svaDir Character. Directory containing SVA figures, typically `.tiff`
#'   files.
#' @param glmDir Character. Directory containing GLM figures, including optional
#'   phenotype subdirectories.
#' @param glmmDir Character. Directory containing GLMM figures, including optional
#'   phenotype subdirectories.
#' @param figDir Character. Directory used by `knitr` for temporary report figure
#'   output while rendering.
#' @param reportTitle Character. Title shown on the rendered report.
#' @param author Character. Author name shown on the rendered report.
#' @param date Character. Date string shown on the rendered report.
#' @param render Logical. If `TRUE`, render the PDF report. If `FALSE`, only
#'   prepare and return the report configuration object.
#' @param verbose Logical. If `TRUE`, emit progress messages with `message()`. The
#'   default is `FALSE`, so the function is quiet unless requested.
#' @param logs Logical. If `TRUE`, write the same progress messages to
#'   `file.path(logDir, "log_dnamReport.txt")`.
#' @param logDir Character. Directory for optional log files.
#' @param templatePath Character or `NULL`. Optional path to a custom report
#'   template. When `NULL`, the packaged `dnamReport.Rmd` template is used.
#'
#' @return A list with class `"dnaEPICO_dnamReport"`.
#' \describe{
#'   \item{preparedReport}{Object returned by [prepareDnamReportInputs()].}
#'   \item{renderResult}{Object returned by [renderDnamReport()] when
#'   `render = TRUE`, otherwise `NULL`.}
#'   \item{status}{Final status string such as `"prepared"`, `"rendered"`,
#'   `"skipped"`, or `"failed"`.}
#'   \item{outputFile}{Intended output path when `render = FALSE`, or the
#'   rendered PDF path when rendering succeeded.}
#'   \item{errorMessage}{Final render error or skip message when available,
#'   otherwise `NULL`.}
#'   \item{logFile}{Resolved path to the optional log file, or `NULL` when
#'   logging was disabled.}
#' }
#' See [dnaEPICO_dnamReport-class] for a class-level overview.
#'
#' @description
#' `dnamReport()` is a high-level convenience wrapper for the reporting stage of
#' the `dnaEPICO` workflow. It prepares the report inputs, inventories the
#' available figures from previous steps, and by default renders the packaged PDF
#' report template. Set `render = FALSE` to inspect the prepared report object
#' without writing the PDF.
#'
#' @examples
#' tmp <- tempdir()
#' report_info <- dnamReport(
#'   output = "DNAm_Report.pdf",
#'   outputDir = tmp,
#'   qcDir = tmp,
#'   preprocessingDir = tmp,
#'   postprocessingDir = tmp,
#'   svaDir = tmp,
#'   glmDir = tmp,
#'   glmmDir = tmp,
#'   figDir = file.path(tmp, "report_figures"),
#'   reportTitle = "DNA methylation analysis",
#'   author = "School of Biomedical Sciences",
#'   date = format(Sys.Date(), "%B %d, %Y"),
#'   render = FALSE,
#'   verbose = FALSE,
#'   logs = FALSE
#' )
#' report_info$status
#'
#' @export
dnamReport <- function(
    output = "DNAm_Report.pdf",
    outputDir = "reports",
    qcDir = "figures/preprocessingMinfiEwasWater/enMix",
    preprocessingDir = "figures/preprocessingMinfiEwasWater/qc",
    postprocessingDir = "figures/preprocessingMinfiEwasWater/metrics",
    svaDir = "figures/svaEnmix/sva",
    glmDir = "figures/methylationGLM_T1",
    glmmDir = "figures/methylationGLMM_T1T2",
    figDir = "reports/figures",
    reportTitle = "DNA methylation",
    author = "School of Biomedical Sciences",
    date = format(Sys.Date(), "%B %d, %Y"),
    render = TRUE,
    verbose = FALSE,
    logs = FALSE,
    logDir = outputDir,
    templatePath = NULL
) {
  log_path <- resolveLogPathMinfiEwasWater(
    logs = logs,
    log_dir = logDir,
    log_file = "log_dnamReport.txt"
  )

  emitLogMinfiEwasWater(
    c(
      "==== Starting DNA Methylation Report Step ====",
      paste("Start time:                     ", format(Sys.time())),
      paste("Output directory:               ", outputDir),
      paste("Output file:                    ", output),
      paste("QC figure directory:            ", qcDir),
      paste("Preprocessing directory:        ", preprocessingDir),
      paste("Postprocessing directory:       ", postprocessingDir),
      paste("SVA directory:                  ", svaDir),
      paste("GLM directory:                  ", glmDir),
      paste("GLMM directory:                 ", glmmDir),
      paste("Knitr figure directory:         ", figDir),
      paste("Render report:                  ", isTRUE(render)),
      paste("Verbose messages:               ", isTRUE(verbose)),
      paste("Write logs:                     ", isTRUE(logs)),
      "======================================================================="
    ),
    verbose = verbose,
    log_path = log_path
  )

  withLoggedErrorsMinfiEwasWater(
    expr = {
      prepared_report <- prepareDnamReportInputs(
        output = output,
        outputDir = outputDir,
        qcDir = qcDir,
        preprocessingDir = preprocessingDir,
        postprocessingDir = postprocessingDir,
        svaDir = svaDir,
        glmDir = glmDir,
        glmmDir = glmmDir,
        figDir = figDir,
        reportTitle = reportTitle,
        author = author,
        date = date,
        verbose = verbose,
        logs = logs,
        logDir = logDir,
        templatePath = templatePath
      )

      render_result <- NULL
      status <- "prepared"
      error_message <- NULL
      output_file <- prepared_report$outputFile

      if (isTRUE(render)) {
        render_result <- renderDnamReport(
          preparedReport = prepared_report,
          verbose = verbose,
          logs = logs,
          logDir = logDir
        )
        status <- render_result$status
        error_message <- render_result$errorMessage
        if (!is.null(render_result$renderedFile)) {
          output_file <- render_result$renderedFile
        }
      }

      emitLogMinfiEwasWater(
        c(
          "=======================================================================",
          paste("Finished DNA methylation report step:", format(Sys.time())),
          paste("Final status:                   ", status),
          paste("Output file:                    ", output_file),
          if (is.null(error_message)) {
            "Final message:                  none"
          } else {
            paste("Final message:                  ", error_message)
          },
          "======================================================================="
        ),
        verbose = verbose,
        log_path = log_path
      )

      structure(
        list(
          preparedReport = prepared_report,
          renderResult = render_result,
          status = status,
          outputFile = output_file,
          errorMessage = error_message,
          logFile = log_path
        ),
        class = "dnaEPICO_dnamReport"
      )
    },
    log_path = log_path,
    verbose = verbose,
    context = "dnamReport"
  )
}
