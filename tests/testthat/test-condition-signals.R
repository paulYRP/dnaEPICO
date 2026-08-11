condition_message_from_error <- function(expression) {
    tryCatch(
        {
            force(expression)
            NA_character_
        },
        error = conditionMessage
    )
}

test_that("condition signals do not construct messages with paste", {
    source_candidates <- c(
        file.path(testthat::test_path(), "..", "..", "R"),
        file.path(
            testthat::test_path(), "..", "..", "00_pkg_src",
            "dnaEPICO", "R"
        )
    )
    source_candidates <- source_candidates[dir.exists(source_candidates)]
    expect_true(length(source_candidates) > 0L)
    source_dir <- normalizePath(source_candidates[[1L]],
        winslash = "/", mustWork = TRUE
    )
    source_files <- list.files(source_dir,
        pattern = "[.]R$", full.names = TRUE
    )
    signal_names <- c("stop", "warning", "message", "signalCondition")

    contains_paste <- function(expression) {
        if (!is.call(expression)) {
            return(FALSE)
        }
        function_name <- if (is.symbol(expression[[1L]])) {
            as.character(expression[[1L]])
        } else {
            ""
        }
        if (function_name %in% c("paste", "paste0")) {
            return(TRUE)
        }
        any(vapply(
            as.list(expression)[-1L], contains_paste, logical(1)
        ))
    }

    inspect_expression <- function(expression, file) {
        if (!is.call(expression)) {
            return(character(0))
        }
        function_name <- if (is.symbol(expression[[1L]])) {
            as.character(expression[[1L]])
        } else {
            ""
        }
        failures <- character(0)
        arguments <- as.list(expression)[-1L]
        if (function_name %in% signal_names &&
            any(vapply(arguments, contains_paste, logical(1)))) {
            failures <- sprintf(
                "%s: %s", basename(file),
                paste(deparse(expression, width.cutoff = 500L), collapse = " ")
            )
        }
        nested <- unlist(lapply(arguments, inspect_expression, file = file),
            use.names = FALSE
        )
        c(failures, nested)
    }

    failures <- unlist(lapply(source_files, function(file) {
        expressions <- parse(file, keep.source = TRUE)
        unlist(lapply(expressions, inspect_expression, file = file),
            use.names = FALSE
        )
    }), use.names = FALSE)
    failures <- as.character(failures)

    expect_identical(failures, character(0))
})

test_that("revised scalar and collapsed condition messages retain their text", {
    expect_identical(
        condition_message_from_error(
            normalizeModelSectionsDnamReport(c("glm", "unsupported"))
        ),
        paste0(
            "Unsupported modelSections value(s): unsupported. ",
            "Expected any of: glm, lme."
        )
    )

    expect_identical(
        condition_message_from_error(
            parsePrsMapMethylationGLM("status:prs,status-prs")
        ),
        paste0(
            "Each prsMap entry must follow the format 'Phenotype:PRS'. ",
            "Invalid entries: status-prs"
        )
    )

    expected_methods <- paste(stats::p.adjust.methods, collapse = ", ")
    expect_identical(
        condition_message_from_error(
            validatePAdjustmentMethodMethylationModels("unsupported")
        ),
        sprintf("padjmethod must be one of: %s", expected_methods)
    )

    expect_identical(
        condition_message_from_error(
            ensurePersonColumnMethylationLME(data.frame(SID = c("P1", "P2")))
        ),
        paste0(
            "Cannot safely derive 'person' from 'SID': every SID must end in ",
            "A or B. Supply an explicit subject identifier column."
        )
    )

    expect_identical(
        condition_message_from_error(
            coerceCorrelationTimeMethylationLME(c("baseline", "followup"))
        ),
        paste0(
            "correlationVar must be numeric or contain numeric text; ",
            "categorical values cannot define AR1/CAR1 spacing."
        )
    )

    expect_identical(
        condition_message_from_error(
            prepareManhattanDataDnaEpico(
                data.frame(IlmnID = "cg1"), pValueColumn = "P.Value"
            )
        ),
        "Manhattan input is missing column(s): chr, pos, P.Value"
    )

    expect_identical(
        condition_message_from_error(
            validateSampleIdentifiersDnaEpico(
                c("sample1", "sample1"), "Sample identifiers"
            )
        ),
        "Sample identifiers contains duplicate sample identifiers: sample1"
    )

    expect_identical(
        condition_message_from_error(
            validatePhenotypeReplacementSvaEnmix(
                original = data.frame(SampleID = "sample1", age = 20),
                updated = data.frame(SampleID = "sample1", PC1 = 1),
                SampleID = "SampleID", pcColumns = "PC1"
            )
        ),
        "The updated phenotype is missing original column(s): age"
    )
})
