makeGencodeHubTestResourceDnaEpico <- function() {
    genes <- GenomicRanges::GRanges(
        seqnames = c("1", "1", "1"),
        ranges = IRanges::IRanges(
            start = c(100L, 500L, 580L),
            end = c(200L, 520L, 600L)
        ),
        strand = c("+", "+", "-")
    )
    GenomeInfoDb::genome(genes) <- "GRCh38"
    S4Vectors::mcols(genes) <- S4Vectors::DataFrame(
        source = rep("HAVANA", 3L), type = rep("gene", 3L),
        score = rep(NA_character_, 3L), phase = rep(NA_character_, 3L),
        gene_id = c("ENSG1", "ENSG2", "ENSG3"),
        gene_type = rep("protein_coding", 3L),
        gene_name = c("DIRECT", "TIE_PLUS", "TIE_MINUS"),
        level = c("1", "2", "2"), tag = rep(NA_character_, 3L),
        hgnc_id = c("HGNC:1", "HGNC:2", "HGNC:3"),
        havana_gene = c("OTTHUMG1", "OTTHUMG2", "OTTHUMG3"),
        artif_dupl = c("", "", "")
    )
    S4Vectors::metadata(genes)$gencode_release <- "50"
    S4Vectors::metadata(genes)$assembly <- "GRCh38.p14"
    S4Vectors::metadata(genes)$source_filename <- "gencode.v50.annotation.gtf.gz"
    S4Vectors::metadata(genes)$source_sha256 <- "test-sha256"
    S4Vectors::metadata(genes)$coordinate_system <- "1-based closed"

    list(
        genes = genes, release = "50", assembly = "GRCh38.p14",
        genomes = "GRCh38", objectMetadata = S4Vectors::metadata(genes),
        annotationHubId = "AH122278", snapshotDate = "2026-08-01",
        hubRecord = data.frame(
            title = "GENCODE release 50 genes", genome = "GRCh38",
            sourceurl = "https://www.gencodegenes.org/",
            stringsAsFactors = FALSE
        )
    )
}

test_that("GENCODEHub annotation keeps direct and nearest tiers aligned", {
    resource <- makeGencodeHubTestResourceDnaEpico()
    results <- data.frame(
        IlmnID = c("cg_direct", "cg_tie", "cg_invalid"),
        chr = c("chr1", "1", NA_character_),
        pos = c(150L, 550L, NA_integer_),
        GencodeV41_Group = c("Body", "TSS1500", ""),
        GencodeV49_RefGene_Name = c("OLD", "OLD", "OLD"),
        check.names = FALSE
    )

    annotated <- dnaEPICO:::appendGencodeHubAnnotationDnaEpico(
        results, resource
    )

    expect_equal(nrow(annotated$data), nrow(results))
    expect_equal(annotated$data$GencodeV41_Group, results$GencodeV41_Group)
    expect_false("GencodeV49_RefGene_Name" %in% names(annotated$data))
    expect_equal(
        annotated$data$GencodeV50_RefGene_Name[[1L]], "DIRECT"
    )
    expect_true(is.na(
        annotated$data$GencodeV50_NonAnnotated_RefGene_Name[[1L]]
    ))
    expect_equal(
        annotated$data$GencodeV50_NonAnnotated_RefGene_Name[[2L]],
        "TIE_PLUS;TIE_MINUS"
    )
    expect_equal(
        annotated$data$GencodeV50_NonAnnotated_RefGene_TSS_Distance[[2L]],
        "50;50"
    )
    expect_equal(
        annotated$data$GencodeV50_Annotation_Status,
        c("annotated", "non_annotated", "unassigned")
    )
    expect_equal(as.integer(annotated$counts), c(1L, 1L, 1L))
})

test_that("GENCODEHub dictionary and provenance use the resource release", {
    resource <- makeGencodeHubTestResourceDnaEpico()
    dictionary <- dnaEPICO:::buildGencodeHubDictionaryDnaEpico("50")
    metadata <- dnaEPICO:::buildGencodeHubMetadataDnaEpico(
        resource,
        c(annotated = 1L, non_annotated = 2L, unassigned = 3L)
    )

    expect_true(all(startsWith(dictionary$Column, "GencodeV50_")))
    expect_equal(metadata$Value[metadata$Key == "gencode.release"], "50")
    expect_equal(
        metadata$Value[metadata$Key == "gencode.annotationhub_id"],
        "AH122278"
    )
})

test_that("GENCODEHub rejects incompatible array builds", {
    expect_error(
        dnaEPICO:::validateGencodeHubArrayCompatibilityDnaEpico(
            "IlluminaHumanMethylation450kanno.ilmn12.hg19"
        ),
        "requires an hg38 array annotation"
    )
    expect_invisible(
        dnaEPICO:::validateGencodeHubArrayCompatibilityDnaEpico(
            "IlluminaHumanMethylationEPICv2anno.20a1.hg38"
        )
    )
})

test_that("the configured AnnotationHub resource is retrievable", {
    skip_if_not(identical(
        Sys.getenv("DNAPICO_RUN_ANNOTATIONHUB_TESTS"), "true"
    ))
    resource <- dnaEPICO:::resolveGencodeHubResourceDnaEpico()
    expect_equal(resource$annotationHubId, "AH122278")
    expect_match(resource$assembly, "^GRCh38")
    expect_true(length(resource$genes) > 0L)
})
