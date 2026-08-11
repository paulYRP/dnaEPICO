## Rebuild inst/extdata/dnaEPICO.gif from dev/dnaEPICO.mp4.
##
## Run from any directory with:
## Rscript dev/build_gift.R
##
## Set FFMPEG_BIN to the full path of a current FFmpeg executable when
## FFmpeg is not available on PATH.

script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_args)) {
    stop("Run this file with Rscript.", call. = FALSE)
}

script_path <- normalizePath(
    sub("^--file=", "", script_args[[1L]]),
    winslash = "/",
    mustWork = TRUE
)
repository_dir <- normalizePath(
    file.path(dirname(script_path), ".."),
    winslash = "/",
    mustWork = TRUE
)
input_file <- file.path(repository_dir, "dev", "dnaEPICO.mp4")
output_file <- file.path(repository_dir, "inst", "extdata", "dnaEPICO.gif")

if (!file.exists(input_file)) {
    stop(
        sprintf("Input video was not found: %s", input_file),
        call. = FALSE
    )
}

ffmpeg_candidates <- c(
    Sys.getenv("FFMPEG_BIN", unset = ""),
    unname(Sys.which("ffmpeg"))
)

if (.Platform$OS.type == "windows") {
    winget_dir <- file.path(
        Sys.getenv("LOCALAPPDATA", unset = ""),
        "Microsoft",
        "WinGet",
        "Packages"
    )
    if (dir.exists(winget_dir)) {
        ffmpeg_candidates <- c(
            ffmpeg_candidates,
            list.files(
                winget_dir,
                pattern = "^ffmpeg[.]exe$",
                full.names = TRUE,
                recursive = TRUE,
                ignore.case = TRUE
            )
        )
    }
}

ffmpeg_candidates <- unique(ffmpeg_candidates[nzchar(ffmpeg_candidates)])
ffmpeg_candidates <- ffmpeg_candidates[file.exists(ffmpeg_candidates)]

supports_palette <- vapply(ffmpeg_candidates, function(candidate) {
    filters <- suppressWarnings(system2(
        candidate,
        "-filters",
        stdout = TRUE,
        stderr = TRUE
    ))
    all(vapply(
        c("palettegen", "paletteuse"),
        function(filter_name) any(grepl(filter_name, filters, fixed = TRUE)),
        logical(1L)
    ))
}, logical(1L))

if (!any(supports_palette)) {
    stop(
        "A current FFmpeg executable with palettegen and paletteuse is ",
        "required. Install FFmpeg or set FFMPEG_BIN to its full path.",
        call. = FALSE
    )
}

ffmpeg <- ffmpeg_candidates[which(supports_palette)[[1L]]]
message(sprintf("Using FFmpeg: %s", ffmpeg))

## Detect and remove uniform black borders. The most frequently reported crop
## is used so brief transitions do not determine the output dimensions.
crop_log <- suppressWarnings(system2(
    ffmpeg,
    c(
        "-hide_banner", "-loglevel", "info", "-i", shQuote(input_file),
        "-vf", "cropdetect=limit=24:round=16:reset=0",
        "-an", "-f", "null", "-"
    ),
    stdout = TRUE,
    stderr = TRUE
))
crop_matches <- regmatches(
    crop_log,
    gregexpr("crop=[0-9]+:[0-9]+:[0-9]+:[0-9]+", crop_log)
)
crop_values <- sub("^crop=", "", unlist(crop_matches, use.names = FALSE))
crop_values <- crop_values[nzchar(crop_values)]

crop_filter <- ""
if (length(crop_values)) {
    crop_counts <- sort(table(crop_values), decreasing = TRUE)
    detected_crop <- names(crop_counts)[[1L]]
    crop_filter <- sprintf("crop=%s,", detected_crop)
    message(sprintf("Detected crop: %s", detected_crop))
} else {
    message("No uniform border was detected; retaining the complete frame.")
}

frame_filter <- sprintf(
    "%sfps=10,scale=650:-2:flags=lanczos",
    crop_filter
)
filter_graph <- sprintf(
    paste0(
        "[0:v]%s,split[palette_source][gif_source];",
        "[palette_source]palettegen=max_colors=256:stats_mode=diff[palette];",
        "[gif_source][palette]paletteuse=",
        "dither=sierra2_4a:diff_mode=rectangle"
    ),
    frame_filter
)

temporary_output <- tempfile(
    pattern = "dnaEPICO-",
    fileext = ".gif"
)

status <- system2(
    ffmpeg,
    c(
        "-hide_banner", "-loglevel", "warning", "-y",
        "-i", shQuote(input_file),
        "-filter_complex", shQuote(filter_graph),
        "-loop", "0", shQuote(temporary_output)
    )
)

if (!identical(status, 0L) || !file.exists(temporary_output)) {
    stop("FFmpeg did not create the replacement GIF.", call. = FALSE)
}
if (is.na(file.size(temporary_output)) || file.size(temporary_output) == 0) {
    stop("FFmpeg created an empty replacement GIF.", call. = FALSE)
}
if (!file.copy(temporary_output, output_file, overwrite = TRUE)) {
    stop(
        sprintf("Could not replace the GIF at: %s", output_file),
        call. = FALSE
    )
}
unlink(temporary_output)

message(sprintf(
    "Created %s (%.2f MiB).",
    output_file,
    file.size(output_file) / 1024^2
))
