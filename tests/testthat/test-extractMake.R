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
})
