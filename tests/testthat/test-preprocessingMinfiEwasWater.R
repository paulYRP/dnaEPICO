test_that("probe-exclusion IDs can be read from an unlabeled first column", {
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

    result <- dnaEPICO:::readProbeExclusionIdsMinfiEwasWater(
        probeExclusionPath = path,
        featureNames = c("cg00001510", "cg99999999")
    )

    expect_identical(result$files$source, "unlabeled first column")
    expect_identical(result$files$column, "<blank>")
    expect_identical(result$ids, c("cg00001510", "cg00002590", "cg00003969"))
    expect_equal(result$overlap, 1)
})

test_that("probe-exclusion IDs can be read from TargetID files", {
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

    auto_result <- dnaEPICO:::readProbeExclusionIdsMinfiEwasWater(
        probeExclusionPath = path,
        featureNames = c("cg00003969", "cg00004121")
    )
    configured_result <- dnaEPICO:::readProbeExclusionIdsMinfiEwasWater(
        probeExclusionPath = path,
        probeExclusionIdColumn = "TargetID",
        featureNames = c("cg00003969", "cg00004121")
    )
    null_string_result <- dnaEPICO:::readProbeExclusionIdsMinfiEwasWater(
        probeExclusionPath = path,
        probeExclusionIdColumn = "NULL",
        featureNames = c("cg00003969", "cg00004121")
    )

    expect_identical(auto_result$files$source, "standard")
    expect_identical(auto_result$files$column, "TargetID")
    expect_identical(configured_result$files$column, "TargetID")
    expect_identical(null_string_result$files$source, "standard")
    expect_identical(null_string_result$files$column, "TargetID")
    expect_equal(auto_result$overlap, 2)
    expect_equal(configured_result$overlap, 2)
    expect_equal(null_string_result$overlap, 2)
})

test_that("multiple probe-exclusion files are combined with per-file columns", {
    path1 <- tempfile(fileext = ".csv")
    path2 <- tempfile(fileext = ".csv")
    writeLines(
        c(
            "TargetID,score",
            "cg00001510,1",
            "cg00003969,1"
        ),
        con = path1
    )
    writeLines(
        c(
            "IlmnID,score",
            "cg00003969,1",
            "cg00004121,1"
        ),
        con = path2
    )

    result <- dnaEPICO:::readProbeExclusionIdsMinfiEwasWater(
        probeExclusionPath = paste(path1, path2, sep = ";"),
        probeExclusionIdColumn = "TargetID;IlmnID",
        featureNames = c("cg00001510", "cg00004121", "cg99999999")
    )

    expect_identical(result$ids, c("cg00001510", "cg00003969", "cg00004121"))
    expect_identical(result$files$column, c("TargetID", "IlmnID"))
    expect_equal(result$nIds, 3)
    expect_equal(result$overlap, 2)
})

test_that("probe-exclusion ID column errors are clear", {
    path <- tempfile(fileext = ".csv")
    writeLines(
        c(
            "TargetID,47,48,49,50",
            "cg00001510,2,0,1,0"
        ),
        con = path
    )

    expect_error(
        dnaEPICO:::readProbeExclusionIdsMinfiEwasWater(
            probeExclusionPath = path,
            probeExclusionIdColumn = "ProbeID"
        ),
        "probeExclusionIdColumn not found"
    )
})

test_that("legacy cross-reactive reader delegates to probe-exclusion reader", {
    path <- tempfile(fileext = ".csv")
    writeLines(
        c(
            "TargetID,score",
            "cg00001510,1"
        ),
        con = path
    )

    result <- dnaEPICO:::readCrossReactiveIdsMinfiEwasWater(
        crossReactivePath = path,
        crossReactiveIdColumn = "TargetID"
    )

    expect_identical(result$ids, "cg00001510")
    expect_s3_class(result, "dnaEPICO_probeExclusion_ids")
})

test_that("EPICv2 manifest flags select expected exclusion IDs", {
    manifest <- data.frame(
        CH_WGBS_evidence = c("Y", "N", "N", "N"),
        CH_BLAT = c("N", "Y", "N", "N"),
        MissingPos = c("N", "N", "Y", "N"),
        MismatchPos = c("N", "N", "N", "Y"),
        row.names = c("cgWGBS", "cgBLAT", "cgMissing", "cgMismatch"),
        stringsAsFactors = FALSE
    )

    default_result <- dnaEPICO:::extractEpicV2ManifestExclusionIdsMinfiEwasWater(
        manifest = manifest,
        featureNames = c("cgWGBS", "cgMismatch")
    )
    mismatch_result <- dnaEPICO:::extractEpicV2ManifestExclusionIdsMinfiEwasWater(
        manifest = manifest,
        epicV2ManifestFlags = c(
            CH_WGBS_evidence = FALSE,
            CH_BLAT = FALSE,
            MissingPos = FALSE,
            MismatchPos = TRUE
        ),
        featureNames = c("cgWGBS", "cgMismatch")
    )

    expect_identical(default_result$ids, c("cgWGBS", "cgBLAT", "cgMissing"))
    expect_equal(default_result$overlap, 1)
    expect_identical(mismatch_result$ids, "cgMismatch")
    expect_equal(mismatch_result$overlap, 1)
})

test_that("EPICv2 manifest extraction can use an IlmnID column", {
    manifest <- data.frame(
        IlmnID = c("cgWGBS", "cgBLAT"),
        CH_WGBS_evidence = c("Y", "N"),
        CH_BLAT = c("N", "Y"),
        MissingPos = c("N", "N"),
        MismatchPos = c("N", "N"),
        stringsAsFactors = FALSE
    )

    result <- dnaEPICO:::extractEpicV2ManifestExclusionIdsMinfiEwasWater(
        manifest = manifest
    )

    expect_identical(result$ids, c("cgWGBS", "cgBLAT"))
})
