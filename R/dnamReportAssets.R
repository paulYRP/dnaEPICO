## Static dashboard assets used by dnamReport().

reportAssetPathDnaEpico <- function(filename) {
    path <- system.file("extdata", filename, package = "dnaEPICO")
    if (nzchar(path)) {
        return(path)
    }

    source_path <- file.path("inst", "extdata", filename)
    if (file.exists(source_path)) {
        return(normalizePath(source_path, mustWork = TRUE))
    }

    stop(sprintf("The report asset '%s' is unavailable.", filename),
        call. = FALSE
    )
}

siteCssDnamReport <- function() {
    readLines(
        reportAssetPathDnaEpico("dnam-report.css"),
        warn = FALSE,
        encoding = "UTF-8"
    )
}
