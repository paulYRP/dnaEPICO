## Internal numeric conversion helpers ---------------------------------------

.numericTokenPatternDnaEpico <- paste0(
    "^[+-]?(?:",
    "(?:[0-9]+(?:\\.[0-9]*)?)|",
    "(?:\\.[0-9]+)",
    ")(?:[eE][+-]?[0-9]+)?$"
)

coerceNumericDnaEpico <- function(x) {
    if (is.numeric(x)) {
        return(as.numeric(x))
    }

    text <- trimws(as.character(x))
    missing <- is.na(x) | !nzchar(text) |
        toupper(text) %in% c("NA", "N/A", "NULL")
    valid <- !missing & grepl(.numericTokenPatternDnaEpico, text, perl = TRUE)
    output <- rep(NA_real_, length(text))
    output[valid] <- as.numeric(text[valid])
    output
}

coerceIntegerDnaEpico <- function(x) {
    numericValue <- coerceNumericDnaEpico(x)
    valid <- is.finite(numericValue) &
        numericValue == trunc(numericValue) &
        numericValue >= -.Machine$integer.max &
        numericValue <= .Machine$integer.max
    output <- rep(NA_integer_, length(numericValue))
    output[valid] <- as.integer(numericValue[valid])
    output
}
