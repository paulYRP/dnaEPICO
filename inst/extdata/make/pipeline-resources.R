# CPU core limits and R worker counts for Make workflows.

# Validate a positive integer CPU or thread count.
pipelinePositiveInteger <- function(value, name) {
    number <- suppressWarnings(as.numeric(value))
    if (length(number) != 1L || !is.finite(number) || number < 1 ||
        number > .Machine$integer.max || number != floor(number)) {
        stop(name, " must be a positive integer.", call. = FALSE)
    }
    as.integer(number)
}

# Return the per-node PBS CPU allocation, or Inf outside a PBS job.
pbsPipelineCores <- function() {
    if (!nzchar(Sys.getenv("PBS_JOBID"))) {
        return(Inf)
    }
    # Read the per-node CPU counts provided by PBS.
    per_node <- Sys.getenv(c("NCPUS", "PBS_NUM_PPN"), unset = "")
    per_node <- per_node[nzchar(per_node)]
    if (length(per_node)) {
        return(min(vapply(names(per_node), function(name) {
            pipelinePositiveInteger(per_node[[name]], name)
        }, integer(1))))
    }
    qstat <- Sys.which("qstat")
    if (nzchar(qstat)) {
        output <- suppressWarnings(system2(qstat,
            c("-f", shQuote(Sys.getenv("PBS_JOBID"))),
            stdout = TRUE, stderr = TRUE, timeout = 10
        ))
        status <- attr(output, "status")
        if (is.null(status) || identical(as.integer(status), 0L)) {
            return(pbsPipelineCoresFromQstat(output))
        }
    }
    stop("Cannot determine the PBS CPU allocation. The workflow needs ",
        "NCPUS/PBS_NUM_PPN from PBS or access to qstat -f for this job.",
        call. = FALSE)
}

# Read the CPU count from single-node qstat resource attributes.
pbsPipelineCoresFromQstat <- function(output) {
    resource <- function(name) {
        pattern <- paste0("^[[:space:]]*Resource_List[.]", name,
            "[[:space:]]*=[[:space:]]*")
        value <- grep(pattern, output, value = TRUE)
        if (length(value) != 1L) {
            return(NA_integer_)
        }
        pipelinePositiveInteger(trimws(sub(pattern, "", value)),
            paste0("PBS Resource_List.", name))
    }
    nodes <- resource("nodect")
    cores <- resource("ncpus")
    if (!is.na(nodes) && nodes == 1L && !is.na(cores)) {
        return(cores)
    }
    stop("qstat did not identify a single-node CPU allocation. ",
        "This workflow uses local R workers; a multi-node total is not ",
        "a per-node CPU limit.", call. = FALSE)
}

# Set the maximum R workers for each analysis, in the order of MODELS.
pipelineResourcePlan <- function(nCores, nativeThreads, models,
    availableCores = parallelly::availableCores(),
    allocatedCores = pbsPipelineCores()) {
    requested <- pipelinePositiveInteger(nCores, "N_CORES")
    threads <- pipelinePositiveInteger(nativeThreads, "NATIVE_THREADS")
    budget <- min(requested, availableCores, allocatedCores)
    slots <- floor(budget / threads)
    if (slots < 1L) {
        stop("NATIVE_THREADS exceeds the available CPU core limit.",
            call. = FALSE)
    }
    selected <- unique(strsplit(trimws(models), "[[:space:]]+")[[1L]])
    selected <- selected[nzchar(selected)]
    if (any(!grepl("^[[:alnum:]][[:alnum:]_.-]*$", selected)) ||
        any(selected %in% c("all", ".", ".."))) {
        stop("MODELS must contain analysis directory names, ",
            "separated by spaces.", call. = FALSE)
    }
    jobs <- min(length(selected), slots)
    shares <- if (!length(selected)) {
        integer()
    } else if (length(selected) > slots) {
        rep.int(1L, length(selected))
    } else {
        rep.int(as.integer(slots %/% jobs), jobs) +
            as.integer(seq_len(jobs) <= slots %% jobs)
    }
    list(cores = as.integer(budget), slots = as.integer(slots),
        jobs = max(1L, as.integer(jobs)),
        workers = stats::setNames(shares, selected))
}

# Write the CPU limit, R worker counts and analysis limits for Make.
writePipelineResourcePlan <- function(args) {
    if (length(args) != 3L) {
        stop("Expected N_CORES, NATIVE_THREADS and MODELS.", call. = FALSE)
    }
    plan <- pipelineResourcePlan(args[[1L]], args[[2L]], args[[3L]])
    cat(c(plan$cores, plan$slots, plan$jobs,
        paste(names(plan$workers), plan$workers, sep = "=")), sep = " ")
}
