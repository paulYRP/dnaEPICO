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
    expect_true(any(grepl("SCALE_VARS = NULL", makefile, fixed = TRUE)))
    expect_true(any(makefile == "MODEL ?= model1"))
    expect_true(any(makefile == "MODELS = modelA modelB modelC"))
    expect_false(any(grepl(
        "^[A-Za-z0-9_ ?:+.-]+=[^#]*#",
        makefile
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
        "DNAPIPE_MK := $(subst $(DNAPIPE_SPACE),\\ ,$(DNAPIPE_MK_RAW))",
        makefile,
        fixed = TRUE
    )))

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
    expect_true(any(grepl(
        "STEP3_PRIMARY = $(COMBINED_PHENO_METHYLATION)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl("$(STEP3_DERIVED): $(STEP3_PRIMARY)", rules, fixed = TRUE)))
    expect_true(any(grepl(
        "mkdir -p \"$(STEP6)/$(MODEL)\" \"$(MODEL_LOG_DIR)\"",
        rules,
        fixed = TRUE
    )))
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
    expect_true(any(grepl(
        "COVARIATES_ARG = $(call optional_text_arg,$(COVARIATES))",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "covariates = $(COVARIATES_ARG)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "factorVars = $(FACTOR_VARS_ARG)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "scaleVars = $(SCALE_VARS_GLM_ARG)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "scaleVars = $(SCALE_VARS_LME_ARG)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "$(GLM_OUTPUTS) &: $(GLM_PHENO_METHYLATION)",
        rules,
        fixed = TRUE
    )))
    expect_true(any(grepl(
        "$(LME_OUTPUTS) &: $(COMBINED_PHENO_METHYLATION)",
        rules,
        fixed = TRUE
    )))
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
})
