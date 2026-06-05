test_that("cross-reactive IDs can be read from an unlabeled first column", {
    path <- tempfile(fileext = ".csv")
    writeLines(
        c(
            '"","match47","match48","match49","match50","Total"',
            '"cg00001510",2,0,1,0,3',
            '"cg00002590",5,0,0,0,5',
            '"cg00003969",0,0,2,0,2'
        ),
        con = path
    )

    result <- dnaEPICO:::readCrossReactiveIdsMinfiEwasWater(
        crossReactivePath = path,
        featureNames = c("cg00001510", "cg99999999")
    )

    expect_identical(result$source, "unlabeled first column")
    expect_identical(result$column, "")
    expect_identical(result$ids, c("cg00001510", "cg00002590", "cg00003969"))
    expect_equal(result$overlap, 1)
})

test_that("cross-reactive IDs can be read from TargetID files", {
    path <- tempfile(fileext = ".csv")
    writeLines(
        c(
            "TargetID,47,48,49,50",
            "cg00001510,2,0,1,0",
            "cg00003969,0,0,2,0",
            "cg00004121,0,4,1,0"
        ),
        con = path
    )

    auto_result <- dnaEPICO:::readCrossReactiveIdsMinfiEwasWater(
        crossReactivePath = path,
        featureNames = c("cg00003969", "cg00004121")
    )
    configured_result <- dnaEPICO:::readCrossReactiveIdsMinfiEwasWater(
        crossReactivePath = path,
        crossReactiveIdColumn = "TargetID",
        featureNames = c("cg00003969", "cg00004121")
    )

    expect_identical(auto_result$source, "standard")
    expect_identical(auto_result$column, "TargetID")
    expect_identical(configured_result$column, "TargetID")
    expect_equal(auto_result$overlap, 2)
    expect_equal(configured_result$overlap, 2)
})

test_that("cross-reactive ID column errors are clear", {
    path <- tempfile(fileext = ".csv")
    writeLines(
        c(
            "TargetID,47,48,49,50",
            "cg00001510,2,0,1,0"
        ),
        con = path
    )

    expect_error(
        dnaEPICO:::readCrossReactiveIdsMinfiEwasWater(
            crossReactivePath = path,
            crossReactiveIdColumn = "ProbeID"
        ),
        "crossReactiveIdColumn not found"
    )
})
