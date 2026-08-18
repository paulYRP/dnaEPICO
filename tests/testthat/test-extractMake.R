test_that("extractMake copies the packaged Makefile and returns its path", {
    tmp <- withr::local_tempdir()

    expect_error(
        makefile_path <- extractMake(
            destDir = tmp,
            overwrite = TRUE
        ),
        NA
    )

    expect_true(is.character(makefile_path))
    expect_length(makefile_path, 1)
    expect_true(file.exists(makefile_path))
    expect_match(basename(makefile_path), "^Makefile$")
    makefile <- readLines(makefile_path, warn = FALSE)
    rscript_name <- if (.Platform$OS.type == "windows") {
        "Rscript.exe"
    } else {
        "Rscript"
    }
    rscript_path <- normalizePath(
        file.path(R.home("bin"), rscript_name),
        winslash = "/",
        mustWork = TRUE
    )
    rscript_path <- gsub("$", "$$", rscript_path, fixed = TRUE)
    rscript_path <- gsub("#", "\\#", rscript_path, fixed = TRUE)
    expect_true(any(makefile == paste("RSCRIPT ?=", rscript_path)))
    expect_false(any(grepl("@DNAEPICO_RSCRIPT@", makefile, fixed = TRUE)))
    expect_true(all(c(
        "COVARIATES_GLM = Sex",
        "COVARIATES_LME = Sex",
        "FACTOR_VARS_GLM = Sex,TreatmentGroup",
        "FACTOR_VARS_LME = Sex,TreatmentGroup",
        "SCALE_VARS_GLM = NULL",
        "SCALE_VARS_LME = NULL",
        "R_DIR = NULL"
    ) %in% makefile))
    expect_false(any(grepl(
        "^(COVARIATES|FACTOR_VARS|SCALE_VARS)[[:space:]]*[?+:]?=",
        makefile
    )))
    expect_true(any(makefile == "MODEL ?= model1"))
    expect_true(any(makefile == "MODELS = modelA modelB modelC"))
    expect_true(all(c(
        "VENND_GLM_PHENOTYPES = NULL",
        "VENND_GLM_LABELS = NULL",
        "VENND_GLM_OMNIBUS_PHENOTYPES = NULL",
        "VENND_GLM_OMNIBUS_LABELS = NULL",
        "VENND_LME_PHENOTYPES = NULL",
        "VENND_LME_LABELS = NULL",
        "VENND_LME_OMNIBUS_PHENOTYPES = NULL",
        "VENND_LME_OMNIBUS_LABELS = NULL"
    ) %in% makefile))
    expect_false(any(grepl(
        "^[A-Za-z0-9_ ?:+.-]+=[^#]*#",
        makefile
    )))
    expect_true(any(grepl(
        "GLM_OMNIBUS_TEST = FALSE", makefile,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "LME_OMNIBUS_TEST = FALSE", makefile,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "LME_OMNIBUS_DDF = Satterthwaite", makefile,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "REMOVE_SEX_MISMATCH = FALSE",
        makefile,
        fixed = TRUE
    )))
    expect_true(any(grepl("Makefile.rules.pipeline", makefile, fixed = TRUE)))
    expect_true(any(grepl(
        '$(shell "$(RSCRIPT)" -e', makefile, fixed = TRUE
    )))
    expect_false(any(grepl("ARGUMENT NORMALISATION", makefile, fixed = TRUE)))
    expect_true(any(grepl(
        "DNAPIPE_MK := $(subst $(DNAPIPE_SPACE),\\ ,$(DNAPIPE_MK_RAW))",
        makefile,
        fixed = TRUE
    )))

    ordered_variables <- c(
        "SEED", "SEP_TYPE", "SAMPLE_ID", "N_SAMPLES", "ARRAY_TYPE",
        "ANNOTATION_VERSION", "IDAT_FORCE", "TIFF_WIDTH", "TIFF_HEIGHT",
        "TIFF_RES", "SEX_COLUMN", "REMOVE_SEX_MISMATCH",
        "SENTRIX_ID_COLUMN", "SENTRIX_POSITION_COLUMN", "BASENAME_COLUMN",
        "TIME_VAR", "METHYLATION_SCALE", "PHENO_ORDER", "CPG_PREFIX",
        "CPG_LIMIT", "PRS_MAP", "SUMMARY_PVAL", "N_CORES",
        "SAVE_TXT_SUMMARIES", "RESUME_FROM_SUMMARY", "CHUNK_SIZE",
        "FDR_THRESHOLD", "PADJ_METHOD", "ANNOTATION_PACKAGE",
        "ANNOTATION_COLS", "GENCODE_HUB", "QC_CUTOFF", "DET_PTYPE",
        "DET_PTHRESHOLD", "PVAL_THRESHOLD", "CHR_TO_REMOVE",
        "SNPS_TO_REMOVE", "PROBE_EXCLUSION_FILE",
        "PROBE_EXCLUSION_ID_COLUMN", "USE_EPICV2_MANIFEST",
        "EPICV2_MANIFEST_CH_WGBS_EVIDENCE", "EPICV2_MANIFEST_CH_BLAT",
        "EPICV2_MANIFEST_MISSING_POS", "EPICV2_MANIFEST_MISMATCH_POS",
        "MAF_THRESHOLD", "PLOT_GROUP_VAR", "LC_REF", "CTRL_SVA_PERC_VAR",
        "CTRL_SVA_FLAG", "TIMEPOINTS", "COMBINE_TIMEPOINTS",
        "PHENOTYPES_GLM", "COVARIATES_GLM", "FACTOR_VARS_GLM",
        "SCALE_VARS_GLM", "GLM_LIBS", "INTERACTION_GLM",
        "GLM_OMNIBUS_TEST", "SUMMARY_RESIDUAL_SD",
        "SAVE_SIGNIFICANT_CPGS", "SIGNIFICANT_CPG_PVAL",
        "VENND_GLM_PHENOTYPES", "VENND_GLM_LABELS",
        "VENND_GLM_OMNIBUS_PHENOTYPES", "VENND_GLM_OMNIBUS_LABELS",
        "PHENOTYPES_LME", "COVARIATES_LME", "FACTOR_VARS_LME",
        "SCALE_VARS_LME", "PERSON_VAR", "LME_LIBS",
        "LME_CORRELATION_STRUCTURE", "LME_CORRELATION_VAR",
        "INTERACTION_LME", "LME_OMNIBUS_TEST", "LME_OMNIBUS_DDF",
        "SAVE_SIGNIFICANT_INTERACTIONS", "SIGNIFICANT_INTERACTION_PVAL",
        "VENND_LME_PHENOTYPES", "VENND_LME_LABELS",
        "VENND_LME_OMNIBUS_PHENOTYPES", "VENND_LME_OMNIBUS_LABELS"
    )
    variable_lines <- vapply(ordered_variables, function(variable) {
        matches <- grep(sprintf("^%s[[:space:]]*[?]?=", variable), makefile)
        if (length(matches)) matches[[1L]] else NA_integer_
    }, integer(1))
    expect_false(anyNA(variable_lines))
    expect_identical(unname(variable_lines), sort(unname(variable_lines)))

    rules_file <- system.file(
        "extdata", "make", "Makefile.rules.pipeline",
        package = "dnaEPICO",
        mustWork = TRUE
    )
    rules <- readLines(rules_file, warn = FALSE)
    expect_true(any(grepl(
        "STEP1_PRIMARY = $(RDATA_DIR)/$(MODEL)/$(STEP1)/$(OBJ_DIR)/RGSet.RData",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl("$(STEP1_DERIVED): $(STEP1_PRIMARY)", rules, fixed = TRUE)))
    expect_true(any(grepl("STEP3_OUTPUTS = $(sort", rules, fixed = TRUE)))
    expect_true(any(grepl(
        "STEP3_DERIVED = $(filter-out $(STEP3_PRIMARY),$(STEP3_OUTPUTS))",
        rules, fixed = TRUE
    )))
    expect_true(any(grepl("$(STEP3_DERIVED): $(STEP3_PRIMARY)", rules, fixed = TRUE)))
    expect_true(any(grepl(
        "VENND_GLM_PHENOTYPES_ARG ?= $(call optional_text_arg,$(VENND_GLM_PHENOTYPES))",
        rules, fixed = TRUE
    )))
    expect_true(any(grepl(
        "VENND_LME_OMNIBUS_LABELS_ARG ?= $(call optional_text_arg,$(VENND_LME_OMNIBUS_LABELS))",
        rules, fixed = TRUE
    )))
    expect_false(any(grepl("report-vennD", rules, fixed = TRUE)))
    expect_false(any(grepl("VENND_GLM_MODELS", rules, fixed = TRUE)))
    expect_true(any(grepl("dir.create", rules, fixed = TRUE)))
    expect_true(any(grepl(
        "identical(result[['status']], 'rendered')",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "file.exists(result[['outputFile']])",
        rules,
        fixed = TRUE
    )))
    expect_true(any(rules == "f3: report-f3"))
    expect_true(any(rules == "f4: report-f4"))
    expect_true(any(rules == "f3lme: report-f3lme"))
    expect_true(any(rules == "all: report-all"))
    expect_true(any(rules == "report-f3: $(FIRST3)"))
    expect_true(any(rules == "report-f4: $(FIRST4)"))
    expect_true(any(rules == "report-f3lme: $(F3LME)"))
    expect_true(any(rules == "report-all: $(FIRST4) $(LME_OUTPUTS)"))
    expect_true(any(grepl(
        "modelSections = $(REPORT_MODEL_SECTIONS_ARG)",
        rules,
        fixed = TRUE
    )))
    expect_false(any(grepl(
        "$(STEP6)/$(MODEL)/docs/index.html: $(REPORT_INPUTS)",
        rules,
        fixed = TRUE
    )))
    expect_false(any(grepl("result$$status", rules, fixed = TRUE)))
    expect_true(any(grepl(
        "PHENO_LC = $(DATA_DIR)/$(MODEL)/$(STEP1)/phenoLC.csv",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "SVA_MARKER = $(DATA_DIR)/$(MODEL)/$(STEP2)/.sva_complete",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "$(SVA_MARKER): $(PHENO_LC)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "if (file.exists('$(SVA_MARKER)')) unlink('$(SVA_MARKER)')",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "$(STEP3_PRIMARY): $(SVA_MARKER) $(PHENO_LC)",
        rules,
        fixed = TRUE
    )))
    expect_true(all(c(
        "\t  covariates = $(call optional_text_arg,$(COVARIATES_GLM)), \\",
        "\t  factorVars = $(call optional_text_arg,$(FACTOR_VARS_GLM)), \\",
        "\t  scaleVars = $(call optional_text_arg,$(SCALE_VARS_GLM)), \\",
        "\t  covariates = $(call optional_text_arg,$(COVARIATES_LME)), \\",
        "\t  factorVars = $(call optional_text_arg,$(FACTOR_VARS_LME)), \\",
        "\t  scaleVars = $(call optional_text_arg,$(SCALE_VARS_LME)), \\"
    ) %in% rules))
    expect_true(all(c(
        "SEP_TYPE_ARG ?= $(call optional_text_arg,$(SEP_TYPE))",
        "PRS_MAP_ARG ?= $(call optional_text_arg,$(PRS_MAP))",
        "R_DIR_ARG ?= $(call optional_text_arg,$(R_DIR))",
        "PROBE_EXCLUSION_ID_COLUMN_ARG ?= $(call optional_text_arg,$(PROBE_EXCLUSION_ID_COLUMN))",
        "INTERACTION_GLM_ARG ?= $(call optional_text_arg,$(INTERACTION_GLM))",
        "INTERACTION_LME_ARG ?= $(call optional_text_arg,$(INTERACTION_LME))",
        "LME_CORRELATION_VAR_ARG ?= $(call optional_text_arg,$(LME_CORRELATION_VAR))",
        "COMBINE_TIMEPOINTS_ARG ?= $(call optional_text_arg,$(COMBINE_TIMEPOINTS))"
    ) %in% rules))
    expect_true(any(grepl(
        "combineTimepoints = $(COMBINE_TIMEPOINTS_ARG)",
        rules, fixed = TRUE
    )))
    expect_false(any(grepl(
        "(COVARIATES|FACTOR_VARS|SCALE_VARS_(GLM|LME))_ARG",
        rules
    )))
    expect_true(any(rules == paste0(
        "MODEL_CONFIG_ERROR = The dnaEPICO Makefile configuration has ",
        "changed. Please update your Makefile."
    )))
    expect_true(any(grepl(
        "require_model_config = $(if $(filter undefined,$(origin $(1)))",
        rules,
        fixed = TRUE
    )))
    validation_lines <- rules[grepl(
        '"$(RSCRIPT)" -e "invisible(NULL)"', rules, fixed = TRUE
    )]
    expect_length(validation_lines, 2L)
    expect_false(any(grepl("\t@Rscript -e", rules, fixed = TRUE)))
    expect_true(all(c(
        "COVARIATES_GLM", "FACTOR_VARS_GLM", "SCALE_VARS_GLM",
        "GLM_OMNIBUS_TEST", "COVARIATES_LME", "FACTOR_VARS_LME",
        "SCALE_VARS_LME"
    ) %in% unlist(regmatches(
        validation_lines,
        gregexpr(
            "(COVARIATES|FACTOR_VARS|SCALE_VARS)_(GLM|LME)|GLM_OMNIBUS_TEST",
            validation_lines
        )
    ))))
    expect_true(any(grepl(
        "$(GLM_OUTPUTS) &: $(GLM_PHENO_METHYLATION) | validate-glm-config",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "$(LME_OUTPUTS) &: $(COMBINED_PHENO_METHYLATION) | validate-lme-config",
        rules,
        fixed = TRUE
    )))
    person_argument <- which(grepl(
        "personVar = '$(PERSON_VAR)'", rules, fixed = TRUE
    ))
    expect_length(person_argument, 1L)
    expect_identical(
        rules[[person_argument + 1L]],
        "\t  SampleID = '$(SAMPLE_ID)', \\"
    )
    expect_true(any(grepl(
        "reportAssetsDir = '$(GLM_REPORT_ASSET_DIR)'",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "reportAssetsDir = '$(LME_REPORT_ASSET_DIR)'",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "omnibusTest = $(GLM_OMNIBUS_TEST)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "omnibusTest = $(LME_OMNIBUS_TEST)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "omnibusDdf = '$(LME_OMNIBUS_DDF)'",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "vennDPhenotypes = $(GLM_VENND_PHENOTYPES_RUNTIME_ARG)",
        rules, fixed = TRUE
    )))
    expect_true(any(grepl(
        "vennDOmnibusLabels = $(LME_VENND_OMNIBUS_LABELS_RUNTIME_ARG)",
        rules, fixed = TRUE
    )))
    expect_true(any(grepl(
        "resumeFromSummary = $(RESUME_FROM_SUMMARY)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "RESUME_FROM_SUMMARY = TRUE",
        readLines(makefile_path, warn = FALSE),
        fixed = TRUE
    )))
    expect_false(any(grepl("phenoLC_withSVA", rules, fixed = TRUE)))
    expect_false(any(grepl("outputPhenoFile", rules, fixed = TRUE)))
    expect_false(any(grepl("PRE_SVA_PHENO", rules, fixed = TRUE)))
    expect_false(any(grepl("FINAL_PHENO", rules, fixed = TRUE)))
    expect_true(any(rules == "MODEL_LOG_DIR = $(LOGS_DIR)/$(MODEL)"))
    expect_true(all(c(
        "run-%: FORCE", "runf3-%: FORCE", "runf4-%: FORCE",
        "runf3lme-%: FORCE", "clean-model-%: FORCE"
    ) %in% rules))
    expect_true(any(rules == "FORCE:"))
    expect_true(any(grepl("commandArgs(trailingOnly = TRUE)", rules,
        fixed = TRUE
    )))
    expect_true(any(grepl("statuses <- vapply(targets, unlink", rules,
        fixed = TRUE
    )))
    expect_false(any(grepl("$(abspath $(LOGS_DIR)", rules, fixed = TRUE)))
    expect_false(any(grepl("test -e", rules, fixed = TRUE)))
    expect_false(any(grepl("mkdir -p", rules, fixed = TRUE)))
    expect_false(any(grepl("rm -rf", rules, fixed = TRUE)))
    expect_false(any(grepl("/dev/null", rules, fixed = TRUE)))
    expect_false(any(grepl("for m in", rules, fixed = TRUE)))
    expect_false(any(grepl("$(eval run", rules, fixed = TRUE)))
})

test_that("model-specific Make configuration validation is explicit", {
    make_command <- Sys.which("make")
    skip_if(!nzchar(make_command), "GNU Make is not available")

    rscript_name <- if (.Platform$OS.type == "windows") {
        "Rscript.exe"
    } else {
        "Rscript"
    }
    rscript_path <- normalizePath(
        file.path(R.home("bin"), rscript_name),
        winslash = "/",
        mustWork = TRUE
    )
    rscript_path <- gsub("$", "$$", rscript_path, fixed = TRUE)
    rscript_path <- gsub("#", "\\#", rscript_path, fixed = TRUE)
    rscript_config <- paste("RSCRIPT :=", rscript_path)

    rules_file <- system.file(
        "extdata", "make", "Makefile.rules.pipeline",
        package = "dnaEPICO",
        mustWork = TRUE
    )
    rules <- readLines(rules_file, warn = FALSE)
    common_lines <- rules[grepl(
        "^(MODEL_CONFIG_ERROR|require_model_config)[[:space:]]*=",
        rules
    )]
    glm_start <- match(".PHONY: validate-glm-config", rules)
    glm_end <- match(
        "$(GLM_OUTPUTS) &: $(GLM_PHENO_METHYLATION) | validate-glm-config",
        rules
    ) - 1L
    lme_start <- match(".PHONY: validate-lme-config", rules)
    lme_end <- match(
        "$(LME_OUTPUTS) &: $(COMBINED_PHENO_METHYLATION) | validate-lme-config",
        rules
    ) - 1L
    expect_false(anyNA(c(glm_start, glm_end, lme_start, lme_end)))

    validation_rules <- c(
        common_lines,
        rules[glm_start:glm_end],
        rules[lme_start:lme_end],
        ".PHONY: f3",
        "f3:",
        "\t@\"$(RSCRIPT)\" -e \"invisible(NULL)\""
    )
    run_make <- function(config, target) {
        makefile <- file.path(withr::local_tempdir(), "Makefile")
        writeLines(
            c(rscript_config, config, validation_rules),
            makefile,
            useBytes = TRUE
        )
        output <- suppressWarnings(system2(
            make_command,
            c("--no-print-directory", "-f", makefile, target),
            stdout = TRUE,
            stderr = TRUE
        ))
        list(
            output = output,
            status = if (is.null(attr(output, "status"))) {
                0L
            } else {
                as.integer(attr(output, "status"))
            }
        )
    }

    glm_config <- c(
        "COVARIATES_GLM = NULL",
        "FACTOR_VARS_GLM = NULL",
        "SCALE_VARS_GLM = NULL",
        "GLM_OMNIBUS_TEST = FALSE"
    )
    lme_config <- c(
        "COVARIATES_LME = NULL",
        "FACTOR_VARS_LME = NULL",
        "SCALE_VARS_LME = NULL",
        "COMBINED_PHENO_METHYLATION = phenoBT1T2.RData"
    )

    expect_identical(run_make(glm_config, "validate-glm-config")$status, 0L)
    expect_identical(run_make(lme_config, "validate-lme-config")$status, 0L)

    outdated <- run_make(c(
        "COVARIATES = Sex",
        "FACTOR_VARS = Sex",
        "SCALE_VARS = NULL"
    ), "validate-glm-config")
    expect_gt(outdated$status, 0L)
    expect_match(
        paste(outdated$output, collapse = "\n"),
        "The dnaEPICO Makefile configuration has changed\\. Please update your Makefile\\."
    )

    blank <- run_make(c(
        "COVARIATES_GLM =",
        "FACTOR_VARS_GLM = NULL",
        "SCALE_VARS_GLM = NULL",
        "GLM_OMNIBUS_TEST = FALSE"
    ), "validate-glm-config")
    expect_gt(blank$status, 0L)
    expect_match(
        paste(blank$output, collapse = "\n"),
        "The dnaEPICO Makefile configuration has changed\\. Please update your Makefile\\."
    )

    expect_identical(run_make(character(0), "f3")$status, 0L)
    expect_gt(run_make(glm_config, "validate-lme-config")$status, 0L)
    expect_gt(run_make(lme_config, "validate-glm-config")$status, 0L)

    no_combination <- run_make(c(
        "COVARIATES_LME = NULL",
        "FACTOR_VARS_LME = NULL",
        "SCALE_VARS_LME = NULL",
        "COMBINED_PHENO_METHYLATION ="
    ), "validate-lme-config")
    expect_gt(no_combination$status, 0L)
    expect_match(
        paste(no_combination$output, collapse = "\n"),
        "COMBINE_TIMEPOINTS must contain at least one timepoint"
    )
})

test_that("the exported Makefile runs from a project path containing spaces", {
    make_command <- Sys.which("make")
    skip_if(!nzchar(make_command), "GNU Make is not available")

    path_entries <- strsplit(
        Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE
    )[[1L]]
    contains_rscript <- vapply(path_entries, function(path) {
        if (!nzchar(path)) {
            return(FALSE)
        }
        any(file.exists(file.path(path, c("Rscript", "Rscript.exe"))))
    }, logical(1))
    restricted_path <- paste(
        path_entries[!contains_rscript], collapse = .Platform$path.sep
    )
    withr::local_envvar(c(PATH = restricted_path))
    if (.Platform$OS.type != "windows") {
        expect_false(nzchar(Sys.which("Rscript")))
    }

    tmp <- file.path(withr::local_tempdir(), "dnaEPICO project with spaces")
    dir.create(tmp, recursive = TRUE)
    extractMake(tmp, overwrite = TRUE)

    output <- withr::with_dir(tmp, suppressWarnings(system2(
        make_command,
        c("--no-print-directory", "status", "MODEL=model1"),
        stdout = TRUE,
        stderr = TRUE
    )))
    status <- attr(output, "status")
    if (is.null(status)) {
        status <- 0L
    }

    expect_identical(as.integer(status), 0L)
    expect_true(any(grepl("Pipeline Status", output, fixed = TRUE)))

    pheno_path <- file.path(
        tmp, "data", "preprocessingMinfiEwasWater", "pheno.csv"
    )
    dir.create(dirname(pheno_path), recursive = TRUE)
    file.create(pheno_path)
    grouped_output <- withr::with_dir(tmp, suppressWarnings(system2(
        make_command,
        c("--no-print-directory", "-n", "f3_models", "MODELS=model1"),
        stdout = TRUE,
        stderr = TRUE
    )))
    grouped_status <- attr(grouped_output, "status")
    if (is.null(grouped_status)) {
        grouped_status <- 0L
    }
    expect_identical(as.integer(grouped_status), 0L)
    expect_true(any(grepl(
        ">>> Starting f3 for: model1", grouped_output, fixed = TRUE
    )))
})
