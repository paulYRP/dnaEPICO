#' Logging helpers for preprocessingMinfiEwasWater
#'
#' @param lines Character vector of log lines to emit.
#' @param verbose Logical. If `TRUE`, send log lines to the console with
#'   `message()`.
#' @param log_path Character or `NULL`. Path to the log file. If `NULL`, log
#'   lines are not written to disk.
#'
#' @return Invisibly returns `NULL`.
#' @description
#' Internal helper that writes preprocessing log lines to the console and/or a
#' log file while keeping user-facing functions quiet by default.
#' @keywords internal
#' @noRd
emitLogMinfiEwasWater <- function(lines, verbose = FALSE, log_path = NULL) {
    if (length(lines) == 0L) {
        return(invisible(NULL))
    }

    lines <- as.character(lines)

    if (!is.null(log_path)) {
        con <- file(log_path, open = "at")
        on.exit(close(con), add = TRUE)
        writeLines(lines, con = con, sep = "\n", useBytes = TRUE)
    }

    if (isTRUE(verbose)) {
        for (line in lines) {
            message(line)
        }
    }

    invisible(NULL)
}

#' Format call-stack lines for dnaEPICO workflow error logs
#'
#' @param calls List of calls, typically from `sys.calls()`.
#'
#' @return Character vector of formatted stack-trace lines.
#' @description
#' Internal helper that converts R calls into readable log lines while dropping
#' the noisier `tryCatch()` bookkeeping frames.
#' @keywords internal
#' @noRd
formatCallStackMinfiEwasWater <- function(calls) {
    if (length(calls) == 0L) {
        return("Call stack: unavailable")
    }

    call_text <- vapply(calls, function(call) {
        paste(deparse(call,
            width.cutoff = 500L
        ), collapse = " ")
    }, character(1))

    keep <- !grepl(
        "^(tryCatch|tryCatchList|tryCatchOne|doTryCatch|withCallingHandlers)\\b",
        call_text
    )
    keep <- keep & !grepl(
        "^withLoggedErrorsMinfiEwasWater\\b",
        call_text
    )

    filtered_calls <- call_text[keep]
    if (length(filtered_calls) == 0L) {
        filtered_calls <- call_text
    }

    c("Call stack:", paste0(
        "  ", seq_along(filtered_calls),
        ": ", filtered_calls
    ))
}

#' Wrap a workflow block and log any error before rethrowing
#'
#' @param expr Expression to evaluate.
#' @param log_path Character or `NULL`. Path to the log file.
#' @param verbose Logical. If `TRUE`, mirror the error message with `message()`.
#' @param context Character. Short label describing where the error happened.
#'
#' @return Returns the value of `expr` when successful.
#' @description
#' Internal helper that ensures top-level workflow wrappers record fatal errors
#' in
#' their log files before the error is rethrown to the caller.
#' @keywords internal
#' @noRd
withLoggedErrorsMinfiEwasWater <- function(
    expr, log_path = NULL,
    verbose = FALSE, context = "dnaEPICO workflow"
) {
    tryCatch(expr, error = function(e) {
        failing_call <- conditionCall(e)
        call_lines <- formatCallStackMinfiEwasWater(sys.calls())
        emitLogMinfiEwasWater(
            c(
                "============================================================",
                paste("ERROR in", context, ":"), conditionMessage(e),
                if (is.null(failing_call)) {
                    "Failing call: unavailable"
                } else {
                    paste("Failing call:", paste(deparse(failing_call,
                        width.cutoff = 500L
                    ), collapse = " "))
                }, call_lines,
                    "============================================================"
            ),
            verbose = verbose, log_path = log_path
        )
        stop(e)
    })
}

#' Resolve the log file path for preprocessingMinfiEwasWater helpers
#'
#' @param logs Logical. If `TRUE`, create a log file path.
#' @param log_dir Character or `NULL`. Directory used for log files.
#' @param log_file Character. File name used for the log file.
#'
#' @return Character scalar with the resolved log path, or `NULL` when
#'   `logs = FALSE`.
#' @description
#' Internal helper that creates the target log directory when needed and returns
#' the path used for file-based logs.
#' @keywords internal
#' @noRd
resolveLogPathMinfiEwasWater <- function(
    logs = FALSE, log_dir = NULL,
    log_file = "log_preprocessingMinfiEwasWater.txt"
) {
    if (!isTRUE(logs)) {
        return(NULL)
    }

    if (is.null(log_dir)) {
        log_dir <- getwd()
    }

    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    log_dir <- normalizePath(log_dir, winslash = "/", mustWork = FALSE)

    file.path(log_dir, log_file)
}

#' Resolve the phenotype separator used in preprocessingMinfiEwasWater helpers
#'
#' @param sepType Character or `NULL`. Separator description supplied by the
#'   user.
#'
#' @return Character scalar with the resolved separator, or `NULL` for the
#'   default comma separator used by `utils::read.csv()`.
#' @description
#' Internal helper that standardizes separator handling across modular
#' preprocessing functions.
#' @keywords internal
#' @noRd
resolveSeparatorMinfiEwasWater <- function(sepType = NULL) {
    if (length(sepType) == 0L || is.na(sepType[[1L]])) {
        return(NULL)
    }

    sepType <- as.character(sepType[[1L]])
    sepTypeLabel <- trimws(sepType)

    if (!nzchar(sepType) || identical(
        toupper(sepTypeLabel),
        "NULL"
    )) {
        return(NULL)
    }

    if (identical(sepType, "\\t")) {
        return("\t")
    }

    sepType
}

#' Split character options used in preprocessingMinfiEwasWater helpers
#'
#' @param x Character vector or scalar string to split.
#' @param sep Character scalar used to split `x` when `x` has length 1.
#'
#' @return Character vector with empty values removed and whitespace trimmed.
#' @description
#' Internal helper that standardizes comma- and semicolon-separated option
#' parsing across preprocessing helpers.
#' @keywords internal
#' @noRd
splitOptionMinfiEwasWater <- function(x, sep = ",") {
    if (length(x) == 0L || all(is.na(x))) {
        return(character(0))
    }

    if (length(x) == 1L) {
        pieces <- strsplit(as.character(x), split = sep, fixed = TRUE)[[1]]
    } else {
        pieces <- as.character(x)
    }

    pieces <- trimws(pieces)
    pieces[!is.na(pieces) & nzchar(pieces)]
}

#' Validate a scalar logical option
#'
#' @keywords internal
#' @noRd
validateLogicalScalarDnaEpico <- function(value, name) {
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
        stop(name, " must be either TRUE or FALSE.", call. = FALSE)
    }

    value
}

#' Calculate a finite mean without creating NaN
#'
#' @keywords internal
#' @noRd
meanFiniteOrNADnaEpico <- function(values) {
    finite_values <- values[is.finite(values)]
    if (length(finite_values) == 0L) {
        return(NA_real_)
    }

    mean(finite_values)
}

#' Validate sample identifiers used to align package objects
#'
#' @keywords internal
#' @noRd
validateSampleIdentifiersDnaEpico <- function(ids, label) {
    ids <- as.character(ids)

    if (length(ids) == 0L) {
        stop(label, " must contain at least one sample identifier.",
            call. = FALSE
        )
    }

    missing_ids <- is.na(ids) | !nzchar(trimws(ids))
    if (any(missing_ids)) {
        stop(label, " contains missing or blank sample identifiers.",
            call. = FALSE
        )
    }

    duplicated_ids <- unique(ids[duplicated(ids)])
    if (length(duplicated_ids) > 0L) {
        stop(label, " contains duplicate sample identifiers: ",
            paste(duplicated_ids, collapse = ", "),
            call. = FALSE
        )
    }

    ids
}

#' Match sample identifiers without allowing silent missing rows
#'
#' @keywords internal
#' @noRd
matchSampleIdentifiersDnaEpico <- function(
    query, reference,
    queryLabel = "query sample identifiers",
        referenceLabel = "reference sample identifiers",
    requireSameSet = FALSE
) {
    query <- validateSampleIdentifiersDnaEpico(query, queryLabel)
    reference <- validateSampleIdentifiersDnaEpico(
        reference,
        referenceLabel
    )

    matched <- match(query, reference)
    if (anyNA(matched)) {
        stop(queryLabel, " not found in ", referenceLabel, ": ",
            paste(query[is.na(matched)], collapse = ", "),
            call. = FALSE
        )
    }

    if (isTRUE(requireSameSet) && length(query) != length(reference)) {
        extra_reference <- setdiff(reference, query)
        stop(queryLabel, " and ", referenceLabel,
            " do not contain the same samples",
            if (length(extra_reference) > 0L) {
                paste0(": extra reference samples are ", paste(extra_reference,
                    collapse = ", "
                ))
            } else {
                "."
            },
            call. = FALSE
        )
    }

    matched
}

#' Convert reported or predicted sex values to a common representation
#'
#' @keywords internal
#' @noRd
canonicalizeSexDnaEpico <- function(values) {
    original <- values
    text <- trimws(as.character(values))
    lower <- tolower(text)
    missing <- is.na(values) | !nzchar(text)

    code <- rep(NA_integer_, length(values))
    code[!missing & lower %in% c("f", "female", "0")] <- 0L
    code[!missing & lower %in% c("m", "male", "1")] <- 1L

    label <- rep(NA_character_, length(values))
    label[code == 0L & !is.na(code)] <- "F"
    label[code == 1L & !is.na(code)] <- "M"

    unknown <- unique(text[!missing & is.na(code)])

    list(original = original, code = code, label = label, unknown = unknown)
}

#' Resolve sex labels supplied to methylation normalization
#'
#' @keywords internal
#' @noRd
resolveNormalizationSexDnaEpico <- function(colData, sexColumn) {
    if (is.null(sexColumn) || !(sexColumn %in% colnames(colData))) {
        return(NULL)
    }

    reported_sex <- canonicalizeSexDnaEpico(colData[[sexColumn]])
    sex_labels <- reported_sex$label
    if (anyNA(sex_labels) && "PredSex" %in% colnames(colData)) {
        predicted_sex <- canonicalizeSexDnaEpico(colData$PredSex)$label
        replace_with_predicted <- is.na(sex_labels) & !is.na(predicted_sex)
        sex_labels[replace_with_predicted] <-
            predicted_sex[replace_with_predicted]
    }

    if (anyNA(sex_labels)) {
        stop("The normalization sex information contains missing or unsupported values",
            if (length(reported_sex$unknown) > 0L) {
                paste0(": ", paste(reported_sex$unknown, collapse = ", "))
            } else {
                "."
            },
            call. = FALSE
        )
    }

    sex_labels
}

#' Capture preview lines for logging
#'
#' @param x Object to preview.
#' @param n Integer. Number of rows to show when `x` is a matrix or data frame.
#'
#' @return Character vector containing formatted preview lines.
#' @description
#' Internal helper that standardizes preview output for logs.
#' @keywords internal
#' @noRd
previewLinesMinfiEwasWater <- function(x, n = 6L) {
    if (is.matrix(x) || is.data.frame(x)) {
        return(utils::capture.output(utils::head(x, n = n)))
    }

    utils::capture.output(methods::show(x))
}

#' Draw a stored plot object for preprocessing helpers
#'
#' @param plot_object Plot object to draw.
#'
#' @return Invisibly returns `NULL`.
#' @description
#' Internal helper that draws stored ggplot or grid objects without relying on
#' `print()` in package code.
#' @keywords internal
#' @noRd
drawPlotObjectMinfiEwasWater <- function(plot_object) {
    if (inherits(plot_object, "ggplot")) {
        grid::grid.newpage()
        grid::grid.draw(ggplot2::ggplotGrob(plot_object))
        return(invisible(NULL))
    }

    if (grid::is.grob(plot_object)) {
        grid::grid.newpage()
        grid::grid.draw(plot_object)
        return(invisible(NULL))
    }

    stop("plot_object must inherit from 'ggplot' or be a grid grob.",
        call. = FALSE
    )
}

#' Draw and optionally save a plot for preprocessingMinfiEwasWater helpers
#'
#' @param draw_fun Function with no arguments that draws the plot.
#' @param display Logical. If `TRUE`, draw on the active graphics device.
#' @param file Character or `NULL`. TIFF file path used for saved output.
#' @param width Integer. Plot width in pixels when `file` is supplied.
#' @param height Integer. Plot height in pixels when `file` is supplied.
#' @param res Integer. TIFF resolution in DPI when `file` is supplied.
#'
#' @return Invisibly returns `file` when a TIFF is written, otherwise `NULL`.
#' @description
#' Internal helper that centralizes plot saving for preprocessing helpers.
#' @keywords internal
#' @noRd
runPlotMinfiEwasWater <- function(
    draw_fun, display = FALSE,
    file = NULL, width = 2000L, height = 1000L, res = 150L
) {
    if (is.null(file) && !isTRUE(display)) {
        return(invisible(NULL))
    }

    if (!is.null(file)) {
        dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
        save_plot <- function() {
            grDevices::tiff(
                filename = file, width = width, height = height,
                res = res, type = "cairo"
            )
            on.exit(grDevices::dev.off(), add = TRUE)
            draw_fun()
        }
        save_plot()
    }

    if (isTRUE(display)) {
        draw_fun()
    }

    invisible(file)
}

#' Save a named preprocessingMinfiEwasWater object
#'
#' @param object Object to save.
#' @param object_name Character. Object name to preserve in the `.RData` file.
#' @param file Character. Target file path.
#'
#' @return Invisibly returns `file`.
#' @description
#' Internal helper used by the convenience wrapper to preserve the documented
#' `.RData` object names consumed by other package functions.
#' @keywords internal
#' @noRd
saveNamedObjectMinfiEwasWater <- function(
    object, object_name,
    file
) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)

    save_env <- list2env(stats::setNames(list(object), object_name),
        parent = emptyenv()
    )
    save(list = object_name, file = file, envir = save_env)

    invisible(file)
}
