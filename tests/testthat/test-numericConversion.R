test_that("numeric conversion accepts supported representations", {
    input <- c("1", "-2.5", ".75", "4.", "1e-4", "+3E2", "NA", "", NA)
    expected <- c(1, -2.5, 0.75, 4, 1e-4, 300, NA, NA, NA)

    expect_equal(coerceNumericDnaEpico(input), expected)
    expect_equal(coerceNumericDnaEpico(factor(c("1.5", "2.5"))), c(1.5, 2.5))
})

test_that("numeric conversion rejects invalid tokens without warnings", {
    expect_no_warning(
        converted <- coerceNumericDnaEpico(c("one", "1,000", "2x", "Inf"))
    )
    expect_true(all(is.na(converted)))
})

test_that("integer conversion accepts only finite in-range integers", {
    input <- c("1", "2.0", "-3", "4.5", "1e100", "invalid", NA)
    expected <- c(1L, 2L, -3L, NA_integer_, NA_integer_, NA_integer_, NA_integer_)

    expect_no_warning(converted <- coerceIntegerDnaEpico(input))
    expect_equal(converted, expected)
})
