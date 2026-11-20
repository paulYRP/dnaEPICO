# Load resource helpers into an isolated environment for tests.
pipeline_resource_helpers <- function() {
    helpers <- new.env(parent = baseenv())
    sys.source(system.file("extdata", "make", "pipeline-resources.R",
        package = "dnaEPICO", mustWork = TRUE), envir = helpers)
    helpers
}

test_that("model workers share one CPU budget in model order", {
    helpers <- pipeline_resource_helpers()
    plan <- function(cores, models, threads = 1L, available = cores,
        allocated = Inf) {
        helpers$pipelineResourcePlan(cores, threads, models,
            availableCores = available, allocatedCores = allocated)
    }
    large <- plan(64L, "second first third second")
    expect_equal(large$workers, c(second = 22L, first = 21L, third = 21L))
    expect_equal(large$jobs, 3L)
    small <- plan(8L, "second first third")
    expect_equal(small$workers, c(second = 3L, first = 3L, third = 2L))
    queued <- plan(2L, "a b c d e")
    expect_equal(queued$jobs, 2L)
    expect_equal(unname(queued$workers), rep(1L, 5L))
    expect_equal(plan(8L, "a")$workers, c(a = 8L))
    expect_equal(plan(8L, "a b c", threads = 2L)$workers,
        c(a = 2L, b = 1L, c = 1L))
    expect_equal(plan(64L, "a b", available = 16L, allocated = 8L)$cores, 8L)

    for (cores in c(1L, 2L, 3L, 8L, 64L)) {
        for (count in c(1L, 3L, 70L)) {
            models <- paste0("model", seq_len(count), collapse = " ")
            result <- plan(cores, models)
            active <- head(sort(result$workers, decreasing = TRUE),
                result$jobs)
            expect_lte(sum(active), cores)
            expect_true(all(result$workers >= 1L))
            expect_lte(result$jobs, min(count, cores))
        }
    }
    expect_error(plan(1L, "a", threads = 2L), "exceeds")
    expect_error(plan(0L, "a"), "N_CORES")
    expect_error(plan(1.5, "a"), "N_CORES")
    expect_error(plan(8L, "a", threads = "invalid"), "NATIVE_THREADS")
    expect_error(plan(8L, "../outside"), "directory names")
})

test_that("PBS allocation uses CPUs rather than MPI slots", {
    helpers <- pipeline_resource_helpers()
    pbs <- c("Job Id: 123.aqua", "    exec_host = node/0",
        "    Resource_List.ncpus = 64", "    Resource_List.nodect = 1",
        "    Resource_List.select = 1:ncpus=64:mpiprocs=1")
    expect_equal(helpers$pbsPipelineCoresFromQstat(pbs), 64L)
    expect_error(helpers$pbsPipelineCoresFromQstat(
        sub("nodect = 1", "nodect = 2", pbs)), "single-node")
    expect_error(helpers$pbsPipelineCoresFromQstat("exec_host = node/0"),
        "single-node")
    withr::local_envvar(c(PBS_JOBID = "123.aqua", NCPUS = "64",
        PBS_NUM_PPN = "", PBS_NP = "1"))
    expect_equal(helpers$pbsPipelineCores(), 64L)
    withr::local_envvar(c(NCPUS = "", PBS_NUM_PPN = "8"))
    expect_equal(helpers$pbsPipelineCores(), 8L)
    withr::local_envvar(c(NCPUS = "", PBS_NUM_PPN = "",
        PATH = withr::local_tempdir()))
    expect_error(helpers$pbsPipelineCores(), "Cannot determine")
    withr::local_envvar(PBS_JOBID = "")
    expect_equal(helpers$pbsPipelineCores(), Inf)
})

test_that("recursive Make queues models and passes worker and thread limits", {
    make_command <- unname(Sys.which("make"))
    skip_if(!nzchar(make_command), "GNU Make is not available")
    thread_names <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
        "OPENBLAS_DEFAULT_NUM_THREADS", "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
    withr::local_envvar(c(MAKE = make_command,
        stats::setNames(rep("8", length(thread_names)), thread_names)))
    tmp <- file.path(withr::local_tempdir(), "pipeline with spaces")
    dir.create(tmp)
    withr::local_dir(tmp)
    rscript <- normalizePath(file.path(R.home("bin"),
        if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"),
        winslash = "/", mustWork = TRUE)
    rscript <- gsub("$", "$$", rscript, fixed = TRUE)
    rscript <- gsub("#", "\\#", rscript, fixed = TRUE)
    rules <- readLines(system.file("extdata", "make", "Makefile.rules.pipeline",
        package = "dnaEPICO", mustWork = TRUE), warn = FALSE)
    # Run model dispatch with lightweight target fixtures.
    first_step <- grep("^MODEL_LOG_DIR =", rules)[[1L]]
    writeLines(c(paste("RSCRIPT :=", rscript), "N_CORES := 2",
        "MODELS := second first second third", "MODEL := single",
        rules[seq_len(first_step - 1L)],
        ".PHONY: FORCE report-f3 report-f4 report-f3lme report-all",
        "FORCE:", "report-f3 report-f4 report-f3lme report-all:",
        '\t@"$(RSCRIPT)" probe.R "$(MODEL)" "$(DNAEPICO_MODEL_WORKERS)"',
        "probe:",
        '\t@echo PLAN $(DNAEPICO_CPU_BUDGET) $(MODEL_JOBS) $(DNAEPICO_MODEL_SHARES)'),
        "Makefile")
    writeLines(c("args <- commandArgs(TRUE)",
        paste0("threads <- Sys.getenv(",
            paste(deparse(thread_names), collapse = " "), ")"),
        "writeLines(c(args[[2L]], threads), paste0(args[[1L]], '.txt'))"),
        "probe.R")
    for (threads in c(1L, 2L)) {
        args <- if (threads == 1L) character() else
            c("NATIVE_THREADS=2", paste0(thread_names, "=8"))
        output <- suppressWarnings(system2(make_command,
            c("--no-print-directory", "f4_models", args),
            stdout = TRUE, stderr = TRUE))
        status <- attr(output, "status")
        expect_true(is.null(status) || identical(as.integer(status), 0L),
            info = paste(output, collapse = "\n"))
        expect_setequal(list.files(pattern = "[.]txt$"),
            c("second.txt", "first.txt", "third.txt"))
        for (model in c("second", "first", "third")) {
            expect_identical(readLines(paste0(model, ".txt")),
                c("1", rep(as.character(threads), length(thread_names))),
                info = paste(model, "with", threads, "native threads"))
        }
    }
    output <- system2(make_command, c("--no-print-directory", "probe"),
        stdout = TRUE, stderr = TRUE)
    expect_true(any(grepl("PLAN 2 2 second=1 first=1 third=1", output,
        fixed = TRUE)), info = paste(output, collapse = "\n"))
})

test_that("summary chunks honor CPU caps and the serial backend", {
    local_mocked_bindings(makeCluster = function(...) {
        stop("A cluster must not be created under a one-worker limit.")
    }, .package = "parallel")
    withr::local_envvar(c(DNAEPICO_MAX_WORKERS = "1",
        DNAEPICO_PARALLEL_BACKEND = "psock"))
    expect_equal(dnaEPICO:::summaryWorkersMethylationModels(64L, 100L), 1L)
    # Use empty fits to exercise serial GLM and LME summary extraction.
    fits <- stats::setNames(rep(list(NULL), 3L), paste0("cg", seq_len(3L)))
    glm <- dnaEPICO:::runGlmSummaryChunksDnaEpico(fits, "score", NULL,
        FALSE, 64L, 1L, .libPaths(), "glm2")
    lme <- dnaEPICO:::runLmeSummaryChunksDnaEpico(fits, "score", NULL, 64L, 1L)
    expect_s3_class(glm$summary, "data.frame")
    expect_s3_class(lme$summary, "data.frame")
    withr::local_envvar(c(DNAEPICO_MAX_WORKERS = "8",
        DNAEPICO_PARALLEL_BACKEND = "serial"))
    expect_equal(dnaEPICO:::summaryWorkersMethylationModels(8L, 100L), 1L)
    expect_error(dnaEPICO:::runGlmSummaryChunksDnaEpico(fits, "score", NULL,
        FALSE, 8L, 1L, .libPaths(), "glm2"), NA)
    expect_error(dnaEPICO:::runLmeSummaryChunksDnaEpico(
        fits, "score", NULL, 8L, 1L), NA)
})
