#' @import data.table
NULL

.gencodeHubResourceId <- "AH122278"

gencodeHubFieldSpecificationDnaEpico <- function() {
    data.frame(
        source = c(
            "seqnames", "source", "type", "start", "end", "score",
            "strand", "phase", "gene_id", "gene_type", "gene_name",
            "level", "tag", "hgnc_id", "havana_gene", "artif_dupl"
        ),
        suffix = c(
            "Chr", "Source", "Feature", "Start", "End", "Score",
            "Strand", "Phase", "ID", "Type", "Name", "Level", "Tag",
            "HGNC_ID", "Havana_Gene", "Artificial_Duplicate"
        ),
        description = c(
            "Chromosome in the GENCODE reference assembly.",
            "Annotation source recorded in the GENCODE GTF.",
            "Feature type recorded in the GENCODE GTF.",
            "Gene start coordinate in the GENCODE reference assembly.",
            "Gene end coordinate in the GENCODE reference assembly.",
            "Score recorded for the gene feature in the GENCODE GTF.",
            "Gene transcription strand.",
            "Phase recorded for the gene feature in the GENCODE GTF.",
            "GENCODE gene identifier.",
            "GENCODE gene biotype.",
            "GENCODE gene name or symbol.",
            "GENCODE annotation level.",
            "GENCODE annotation tag.",
            "HGNC gene identifier.",
            "HAVANA gene identifier.",
            "GENCODE artificial duplication attribute."
        ),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}

normalizeGencodeHubChromosomeDnaEpico <- function(x) {
    value <- trimws(as.character(x))
    value[!nzchar(value)] <- NA_character_
    value <- sub("^chr", "", value, ignore.case = TRUE)
    value[toupper(value) == "M"] <- "MT"
    standard <- toupper(value) %in% c("X", "Y", "MT")
    value[standard] <- toupper(value[standard])
    value
}

collapseAlignedGencodeHubDnaEpico <- function(x) {
    value <- as.character(x)
    value[is.na(value) | !nzchar(value)] <- "NA"
    paste(value, collapse = ";")
}

singleGencodeHubMetadataValueDnaEpico <- function(metadata, name) {
    value <- metadata[[name]]
    if (is.null(value) || length(value) != 1L || is.na(value) ||
        !nzchar(trimws(as.character(value)))) {
        stop(
            "GENCODEHub resource metadata must contain one non-empty '",
            name, "' value.",
            call. = FALSE
        )
    }
    trimws(as.character(value))
}

validateGencodeHubGenomeDnaEpico <- function(genes, assembly) {
    if (!startsWith(assembly, "GRCh38")) {
        stop("GENCODEHub annotation requires a GRCh38 resource; received ",
            assembly, ".", call. = FALSE
        )
    }
    resource_genomes <- unique(as.character(GenomeInfoDb::genome(genes)))
    resource_genomes <- resource_genomes[
        !is.na(resource_genomes) & nzchar(resource_genomes)
    ]
    if (length(resource_genomes) &&
        any(!startsWith(resource_genomes, "GRCh38"))) {
        resource_genomes_text <- paste(resource_genomes, collapse = ", ")
        stop(sprintf(
            "GENCODEHub GRanges genome is incompatible with GRCh38: %s.",
            resource_genomes_text
        ), call. = FALSE)
    }
    resource_genomes
}

validateGencodeHubGeneDataDnaEpico <- function(genes) {
    fields <- gencodeHubFieldSpecificationDnaEpico()$source
    gene_data <- as.data.frame(genes, stringsAsFactors = FALSE)
    missing_fields <- setdiff(fields, names(gene_data))
    if (length(missing_fields)) {
        missing_fields_text <- paste(missing_fields, collapse = ", ")
        stop(sprintf(
            "GENCODEHub resource is missing required gene field(s): %s.",
            missing_fields_text
        ), call. = FALSE)
    }
    feature_types <- unique(as.character(gene_data$type))
    feature_types <- feature_types[
        !is.na(feature_types) & nzchar(feature_types)
    ]
    if (length(feature_types) != 1L || !identical(feature_types, "gene")) {
        stop("GENCODEHub resource must contain gene-level features only.",
            call. = FALSE
        )
    }
    invisible(gene_data)
}

validateGencodeHubResourceDnaEpico <- function(genes) {
    if (!methods::is(genes, "GRanges")) {
        stop(
            "The GENCODEHub gene resource must be a GRanges object.",
            call. = FALSE
        )
    }
    if (length(genes) == 0L) {
        stop("The GENCODEHub gene resource is empty.", call. = FALSE)
    }

    object_metadata <- S4Vectors::metadata(genes)
    release <- singleGencodeHubMetadataValueDnaEpico(
        object_metadata, "gencode_release"
    )
    assembly <- singleGencodeHubMetadataValueDnaEpico(
        object_metadata, "assembly"
    )
    resource_genomes <- validateGencodeHubGenomeDnaEpico(genes, assembly)
    validateGencodeHubGeneDataDnaEpico(genes)

    list(
        genes = genes, release = release, assembly = assembly,
        genomes = resource_genomes, objectMetadata = object_metadata
    )
}

annotationHubRecordValueDnaEpico <- function(record, name) {
    if (!is.data.frame(record) || !(name %in% names(record)) ||
        nrow(record) != 1L) {
        return(NA_character_)
    }
    value <- record[[name]][[1L]]
    if (is.null(value) || !length(value) || all(is.na(value))) {
        return(NA_character_)
    }
    paste(as.character(value), collapse = ",")
}

resolveGencodeHubResourceDnaEpico <- function(hub = NULL) {
    if (is.null(hub)) {
        hub <- AnnotationHub::AnnotationHub(ask = FALSE)
    }
    if (!methods::is(hub, "AnnotationHub")) {
        stop("hub must be an AnnotationHub object.", call. = FALSE)
    }
    resource_id <- .gencodeHubResourceId
    if (!(resource_id %in% names(hub))) {
        stop(
            "GENCODEHub resource ", resource_id,
            " is unavailable in AnnotationHub snapshot ",
            as.character(AnnotationHub::snapshotDate(hub)),
            ". Update the AnnotationHub metadata for the installed ",
            "Bioconductor version and try again.",
            call. = FALSE
        )
    }

    hub_record <- as.data.frame(
        S4Vectors::mcols(hub[resource_id]),
        stringsAsFactors = FALSE
    )
    genes <- tryCatch(
        hub[[resource_id]],
        error = function(error) {
            stop(
                "Could not retrieve GENCODEHub resource ", resource_id,
                ": ", conditionMessage(error),
                call. = FALSE
            )
        }
    )
    validated <- validateGencodeHubResourceDnaEpico(genes)

    hub_genome <- annotationHubRecordValueDnaEpico(hub_record, "genome")
    if (!is.na(hub_genome) && !startsWith(hub_genome, "GRCh38")) {
        stop(
            "AnnotationHub metadata for ", resource_id,
            " is incompatible with GRCh38: ", hub_genome, ".",
            call. = FALSE
        )
    }

    validated$annotationHubId <- resource_id
    validated$snapshotDate <- as.character(AnnotationHub::snapshotDate(hub))
    validated$hubRecord <- hub_record
    validated
}

validateGencodeHubArrayCompatibilityDnaEpico <- function(annotationObject) {
    if (is.character(annotationObject) && length(annotationObject) == 1L &&
        !is.na(annotationObject) && nzchar(annotationObject) &&
        !grepl("hg38", annotationObject, ignore.case = TRUE)) {
        stop(
            "GENCODEHub annotation requires an hg38 array annotation; got: ",
            annotationObject, ".",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

prepareGencodeHubCoordinatesDnaEpico <- function(coordinates) {
    coordinates <- as.data.frame(coordinates,
        stringsAsFactors = FALSE, check.names = FALSE
    )
    if (!all(c("chr", "pos") %in% names(coordinates))) {
        stop("coordinates must contain chr and pos columns.", call. = FALSE)
    }
    coordinates$chr <- normalizeGencodeHubChromosomeDnaEpico(coordinates$chr)
    numeric_position <- coerceNumericDnaEpico(coordinates$pos)
    integer_position <- coerceIntegerDnaEpico(numeric_position)
    coordinates$pos <- integer_position
    coordinates$valid <- !is.na(coordinates$chr) &
        !is.na(integer_position) & is.finite(numeric_position) &
        numeric_position == integer_position & integer_position > 0L
    coordinates <- unique(
        coordinates[, c("chr", "pos", "valid"), drop = FALSE]
    )
    coordinates$coordID <- seq_len(nrow(coordinates))
    coordinates <- coordinates[
        , c("coordID", "chr", "pos", "valid"), drop = FALSE
    ]
    if (!any(coordinates$valid)) {
        stop("GENCODEHub annotation requires at least one valid CpG chr/pos ",
            "coordinate.", call. = FALSE
        )
    }
    coordinates
}

prepareGencodeHubRangesDnaEpico <- function(genes, coordinates) {
    specification <- gencodeHubFieldSpecificationDnaEpico()
    gene_data <- as.data.frame(genes, stringsAsFactors = FALSE)
    gene_data$normalized_chr <- normalizeGencodeHubChromosomeDnaEpico(
        gene_data$seqnames
    )
    used_chromosomes <- unique(coordinates$chr[coordinates$valid])
    keep <- gene_data$normalized_chr %in% used_chromosomes &
        !is.na(gene_data$gene_id) & nzchar(as.character(gene_data$gene_id)) &
        !is.na(gene_data$gene_name) &
        nzchar(as.character(gene_data$gene_name)) &
        as.character(gene_data$strand) %in% c("+", "-")
    gene_data <- gene_data[keep, , drop = FALSE]
    gene_data <- gene_data[order(
        gene_data$normalized_chr, gene_data$start, gene_data$gene_id,
        na.last = TRUE
    ), , drop = FALSE]
    rownames(gene_data) <- NULL
    if (!nrow(gene_data)) {
        stop("No GENCODEHub genes match the supplied CpG chromosomes.",
            call. = FALSE
        )
    }
    valid_rows <- which(coordinates$valid)
    cpg_ranges <- GenomicRanges::GRanges(
        seqnames = coordinates$chr[valid_rows],
        ranges = IRanges::IRanges(
            start = coordinates$pos[valid_rows], width = 1L
        )
    )
    S4Vectors::mcols(cpg_ranges)$coordID <- coordinates$coordID[valid_rows]
    gene_ranges <- GenomicRanges::GRanges(
        seqnames = gene_data$normalized_chr,
        ranges = IRanges::IRanges(
            start = gene_data$start, end = gene_data$end
        ), strand = gene_data$strand
    )
    list(
        specification = specification, geneData = gene_data,
        cpgRanges = cpg_ranges, geneRanges = gene_ranges,
        tssRanges = GenomicRanges::resize(
            gene_ranges, width = 1L, fix = "start"
        )
    )
}

directGencodeHubRelationshipsDnaEpico <- function(cpgRanges, geneRanges) {
    hits <- GenomicRanges::findOverlaps(
        cpgRanges, geneRanges, ignore.strand = TRUE
    )
    unique(data.table::data.table(
        coordID = S4Vectors::mcols(cpgRanges)$coordID[
            S4Vectors::queryHits(hits)
        ],
        geneROW = S4Vectors::subjectHits(hits), annotation = "annotated"
    ))
}

seedNearestGencodeHubRangesDnaEpico <- function(
    cpgRanges, tssRanges, nonDirectIds
) {
    query_index <- match(nonDirectIds, S4Vectors::mcols(cpgRanges)$coordID)
    nearest_seed <- GenomicRanges::distanceToNearest(
        cpgRanges[query_index], tssRanges,
        ignore.strand = TRUE, select = "all"
    )
    if (!length(nearest_seed)) {
        return(NULL)
    }
    seed_query <- S4Vectors::queryHits(nearest_seed)
    seed_distance <- as.integer(S4Vectors::mcols(nearest_seed)$distance)
    minimum_distance <- tapply(seed_distance, seed_query, min, na.rm = TRUE)
    retained_query <- as.integer(names(minimum_distance))
    radius <- as.integer(minimum_distance) + 1L
    query_position <- BiocGenerics::start(
        cpgRanges[query_index][retained_query]
    )
    search_ranges <- GenomicRanges::GRanges(
        seqnames = GenomeInfoDb::seqnames(
            cpgRanges[query_index][retained_query]
        ),
        ranges = IRanges::IRanges(
            start = pmax(1L, query_position - radius),
            end = query_position + radius
        )
    )
    list(
        ranges = search_ranges, positions = query_position,
        radius = radius, retainedQuery = retained_query
    )
}

nearestGencodeHubRelationshipsDnaEpico <- function(
    cpgRanges, tssRanges, nonDirectIds
) {
    empty <- data.table::data.table(
        coordID = integer(), geneROW = integer(), annotation = character()
    )
    if (!length(nonDirectIds)) {
        return(empty)
    }
    seed <- seedNearestGencodeHubRangesDnaEpico(
        cpgRanges, tssRanges, nonDirectIds
    )
    if (is.null(seed)) {
        return(empty)
    }
    candidates <- GenomicRanges::findOverlaps(
        seed$ranges, tssRanges, ignore.strand = TRUE
    )
    candidate_query <- S4Vectors::queryHits(candidates)
    candidate_gene <- S4Vectors::subjectHits(candidates)
    candidate_distance <- abs(
        seed$positions[candidate_query] -
            BiocGenerics::start(tssRanges)[candidate_gene]
    )
    keep <- candidate_distance == seed$radius[candidate_query]
    unique(data.table::data.table(
        coordID = nonDirectIds[seed$retainedQuery[candidate_query[keep]]],
        geneROW = candidate_gene[keep], annotation = "non_annotated"
    ))
}

gencodeHubColumnNamesDnaEpico <- function(prefix, specification) {
    tss_suffixes <- c("TSS", "TSS_Distance", "TSS_Direction")
    direct <- paste0(prefix, "_RefGene_", specification$suffix)
    nearest <- paste0(
        prefix, "_NonAnnotated_RefGene_", specification$suffix
    )
    list(
        direct = c(direct, paste0(prefix, "_RefGene_", tss_suffixes)),
        nearest = c(nearest, paste0(
            prefix, "_NonAnnotated_RefGene_", tss_suffixes
        ))
    )
}

buildGencodeHubFieldsDnaEpico <- function(geneData, specification, prefix) {
    fields <- data.table::as.data.table(
        geneData[, specification$source, drop = FALSE]
    )
    data.table::setnames(
        fields, specification$source,
        paste0(prefix, "_RefGene_", specification$suffix)
    )
    data.table::set(fields, j = "geneROW", value = seq_len(nrow(fields)))
    fields
}

addGencodeHubTssColumnsDnaEpico <- function(
    annotation, coordinates, tssRanges, prefix
) {
    tss_names <- paste0(prefix, "_RefGene_", c(
        "TSS", "TSS_Distance", "TSS_Direction"
    ))
    data.table::set(annotation, j = "cpgPOS", value = coordinates$pos[
        match(annotation[["coordID"]], coordinates$coordID)
    ])
    annotation[[tss_names[[1L]]]] <-
        BiocGenerics::start(tssRanges)[annotation$geneROW]
    annotation[[tss_names[[2L]]]] <- abs(
        annotation$cpgPOS - annotation[[tss_names[[1L]]]]
    )
    strand_name <- paste0(prefix, "_RefGene_Strand")
    annotation[[tss_names[[3L]]]] <- ifelse(
        annotation[[tss_names[[2L]]]] == 0L, "at_tss",
        ifelse(
            (annotation[[strand_name]] == "+" &
                annotation$cpgPOS < annotation[[tss_names[[1L]]]]) |
                (annotation[[strand_name]] == "-" &
                    annotation$cpgPOS > annotation[[tss_names[[1L]]]]),
            "upstream", "downstream"
        )
    )
    annotation
}

buildGencodeHubLongAnnotationDnaEpico <- function(
    relationships, ranges, coordinates, prefix
) {
    annotation <- merge(
        relationships,
        buildGencodeHubFieldsDnaEpico(
            ranges$geneData, ranges$specification, prefix
        ),
        by = "geneROW", all.x = TRUE, sort = FALSE
    )
    annotation <- addGencodeHubTssColumnsDnaEpico(
        annotation, coordinates, ranges$tssRanges, prefix
    )
    distance_name <- paste0(prefix, "_RefGene_TSS_Distance")
    data.table::setorderv(annotation, c(
        "coordID", "annotation", distance_name,
        paste0(prefix, "_RefGene_ID")
    ))
    unique(annotation,
        by = c("coordID", "annotation", paste0(prefix, "_RefGene_ID"))
    )
}

collapseGencodeHubRelationshipMapDnaEpico <- function(
    coordinates, relationships, ranges, columnNames, prefix
) {
    coordinate_map <- data.table::as.data.table(coordinates)
    for (column in c(columnNames$direct, columnNames$nearest)) {
        coordinate_map[[column]] <- NA_character_
    }
    if (!nrow(relationships)) {
        return(coordinate_map)
    }
    long_annotation <- buildGencodeHubLongAnnotationDnaEpico(
        relationships, ranges, coordinates, prefix
    )
    direct_map <- long_annotation[
        long_annotation[["annotation"]] == "annotated",
        lapply(.SD, collapseAlignedGencodeHubDnaEpico),
        by = "coordID", .SDcols = columnNames$direct
    ]
    nearest_map <- long_annotation[
        long_annotation[["annotation"]] == "non_annotated",
        lapply(.SD, collapseAlignedGencodeHubDnaEpico),
        by = "coordID", .SDcols = columnNames$direct
    ]
    data.table::setnames(
        nearest_map, columnNames$direct, columnNames$nearest
    )
    coordinate_map <- merge(
        coordinate_map[,
            setdiff(colnames(coordinate_map), columnNames$direct),
            with = FALSE
        ], direct_map, by = "coordID", all.x = TRUE, sort = FALSE
    )
    coordinate_map <- merge(
        coordinate_map[,
            setdiff(colnames(coordinate_map), columnNames$nearest),
            with = FALSE
        ], nearest_map, by = "coordID", all.x = TRUE, sort = FALSE
    )
    data.table::setorderv(coordinate_map, "coordID")
    coordinate_map
}

finalizeGencodeHubCoordinateMapDnaEpico <- function(
    coordinateMap, prefix, release, columnNames
) {
    direct_name <- paste0(prefix, "_RefGene_Name")
    nearest_name <- paste0(prefix, "_NonAnnotated_RefGene_Name")
    status_name <- paste0(prefix, "_Annotation_Status")
    coordinateMap[[status_name]] <- ifelse(
        !is.na(coordinateMap[[direct_name]]), "annotated",
        ifelse(
            !is.na(coordinateMap[[nearest_name]]),
            "non_annotated", "unassigned"
        )
    )
    list(
        data = as.data.frame(coordinateMap, check.names = FALSE),
        release = release, prefix = prefix,
        annotationColumns = c(
            columnNames$direct, columnNames$nearest, status_name
        ),
        statusColumn = status_name
    )
}

annotateCoordinatesGencodeHubDnaEpico <- function(coordinates, resource) {
    if (is.null(resource$genes) || is.null(resource$release)) {
        stop("resource must be returned by a GENCODEHub resource resolver.",
            call. = FALSE
        )
    }
    validated <- validateGencodeHubResourceDnaEpico(resource$genes)
    if (!identical(as.character(resource$release), validated$release)) {
        stop("GENCODEHub resource release metadata are inconsistent.",
            call. = FALSE
        )
    }
    coordinates <- prepareGencodeHubCoordinatesDnaEpico(coordinates)
    ranges <- prepareGencodeHubRangesDnaEpico(resource$genes, coordinates)
    direct <- directGencodeHubRelationshipsDnaEpico(
        ranges$cpgRanges, ranges$geneRanges
    )
    non_direct_ids <- setdiff(
        S4Vectors::mcols(ranges$cpgRanges)$coordID, direct$coordID
    )
    nearest <- nearestGencodeHubRelationshipsDnaEpico(
        ranges$cpgRanges, ranges$tssRanges, non_direct_ids
    )
    relationships <- data.table::rbindlist(
        list(direct, nearest), use.names = TRUE, fill = TRUE
    )
    prefix <- paste0("GencodeV", validated$release)
    column_names <- gencodeHubColumnNamesDnaEpico(
        prefix, ranges$specification
    )
    coordinate_map <- collapseGencodeHubRelationshipMapDnaEpico(
        coordinates, relationships, ranges, column_names, prefix
    )
    finalizeGencodeHubCoordinateMapDnaEpico(
        coordinate_map, prefix, validated$release, column_names
    )
}
countAlignedGencodeHubValuesDnaEpico <- function(x) {
    output <- lengths(strsplit(as.character(x), ";", fixed = TRUE))
    output[is.na(x)] <- 0L
    as.integer(output)
}

prepareGencodeHubResultColumnsDnaEpico <- function(results, resource) {
    results <- as.data.frame(results,
        stringsAsFactors = FALSE, check.names = FALSE
    )
    if (!all(c("chr", "pos") %in% names(results))) {
        stop("GENCODEHub annotation requires chr and pos columns in the ",
            "annotated model results.", call. = FALSE
        )
    }
    original_rows <- nrow(results)
    generated_pattern <- paste0(
        "^GencodeV[0-9]+_(RefGene_|NonAnnotated_RefGene_|",
        "Annotation_Status$)"
    )
    results <- results[
        , !grepl(generated_pattern, names(results)), drop = FALSE
    ]
    chromosome <- normalizeGencodeHubChromosomeDnaEpico(results$chr)
    numeric_position <- coerceNumericDnaEpico(results$pos)
    integer_position <- coerceIntegerDnaEpico(numeric_position)
    valid <- !is.na(chromosome) & !is.na(integer_position) &
        is.finite(numeric_position) & numeric_position == integer_position &
        integer_position > 0L
    coordinate_result <- annotateCoordinatesGencodeHubDnaEpico(
        data.frame(chr = chromosome, pos = integer_position), resource
    )
    coordinate_map <- coordinate_result$data
    coordinate_match <- match(
        paste(chromosome, integer_position, sep = ":"),
        paste(coordinate_map$chr, coordinate_map$pos, sep = ":")
    )
    for (column in coordinate_result$annotationColumns) {
        results[[column]] <- coordinate_map[[column]][coordinate_match]
    }
    results[[coordinate_result$statusColumn]][!valid] <- "unassigned"
    if (nrow(results) != original_rows) {
        stop("GENCODEHub annotation changed the number of result rows.",
            call. = FALSE
        )
    }
    list(data = results, annotation = coordinate_result)
}

validateGencodeHubAnnotationTiersDnaEpico <- function(results, prefix) {
    direct_name <- paste0(prefix, "_RefGene_Name")
    nearest_name <- paste0(prefix, "_NonAnnotated_RefGene_Name")
    if (any(!is.na(results[[direct_name]]) &
        !is.na(results[[nearest_name]]))) {
        stop("Direct and nearest GENCODEHub annotation tiers are not ",
            "mutually exclusive.", call. = FALSE
        )
    }
    invisible(TRUE)
}

validateGencodeHubFieldAlignmentDnaEpico <- function(results, prefix) {
    direct_fields <- grep(
        paste0("^", prefix, "_RefGene_"), names(results), value = TRUE
    )
    nearest_fields <- grep(
        paste0("^", prefix, "_NonAnnotated_RefGene_"),
        names(results), value = TRUE
    )
    if (!nrow(results) || length(direct_fields) <= 1L) {
        return(invisible(TRUE))
    }
    direct_counts <- vapply(
        results[, direct_fields, drop = FALSE],
        countAlignedGencodeHubValuesDnaEpico, integer(nrow(results))
    )
    nearest_counts <- vapply(
        results[, nearest_fields, drop = FALSE],
        countAlignedGencodeHubValuesDnaEpico, integer(nrow(results))
    )
    if (any(rowSums(direct_counts != direct_counts[, 1L]) > 0L) ||
        any(rowSums(nearest_counts != nearest_counts[, 1L]) > 0L)) {
        stop("Collapsed GENCODEHub fields do not preserve gene-to-field ",
            "alignment.", call. = FALSE
        )
    }
    invisible(TRUE)
}

appendGencodeHubAnnotationDnaEpico <- function(results, resource) {
    prepared <- prepareGencodeHubResultColumnsDnaEpico(results, resource)
    results <- prepared$data
    annotation <- prepared$annotation
    validateGencodeHubAnnotationTiersDnaEpico(results, annotation$prefix)
    validateGencodeHubFieldAlignmentDnaEpico(results, annotation$prefix)
    status <- results[[annotation$statusColumn]]
    counts <- table(factor(
        status,
        levels = c("annotated", "non_annotated", "unassigned")
    ))
    list(
        data = results, release = annotation$release,
        prefix = annotation$prefix,
        annotationColumns = annotation$annotationColumns,
        statusColumn = annotation$statusColumn,
        counts = counts
    )
}
buildGencodeHubDictionaryDnaEpico <- function(release) {
    specification <- gencodeHubFieldSpecificationDnaEpico()
    prefix <- paste0("GencodeV", release)
    suffixes <- c(
        specification$suffix, "TSS", "TSS_Distance", "TSS_Direction"
    )
    descriptions <- c(
        specification$description,
        "Strand-specific transcription start site coordinate.",
        paste(
            "Absolute base-pair distance between the CpG and",
            "transcription start site."
        ),
        paste(
            "CpG direction relative to transcription: upstream, downstream,",
            "or at_tss."
        )
    )
    direct_note <- paste(
        "These values describe directly overlapping GENCODE release",
        release, "genes. Multiple values are semicolon separated and aligned",
        "across GENCODE fields."
    )
    nearest_note <- paste(
        "These values describe nearest GENCODE release", release,
        "TSS genes assigned only when no gene span contains the CpG.",
        "Equidistant genes are retained and fields remain aligned."
    )
    rbind(
        data.frame(
            Column = paste0(prefix, "_RefGene_", suffixes),
            Description = paste(descriptions, direct_note), Formula = "",
            stringsAsFactors = FALSE, check.names = FALSE
        ),
        data.frame(
            Column = paste0(prefix, "_NonAnnotated_RefGene_", suffixes),
            Description = paste(descriptions, nearest_note), Formula = "",
            stringsAsFactors = FALSE, check.names = FALSE
        ),
        data.frame(
            Column = paste0(prefix, "_Annotation_Status"),
            Description = paste(
                "GENCODE release", release,
                "assignment route: annotated, non_annotated, or unassigned."
            ),
            Formula = "", stringsAsFactors = FALSE, check.names = FALSE
        )
    )
}

buildGencodeHubMetadataDnaEpico <- function(resource, counts) {
    object_metadata <- resource$objectMetadata
    value_or_na <- function(name) {
        value <- object_metadata[[name]]
        if (is.null(value) || !length(value)) {
            NA_character_
        } else {
            as.character(value[[1L]])
        }
    }
    hub_record <- resource$hubRecord
    data.frame(
        Key = c(
            "gencode.enabled", "gencode.release",
            "gencode.annotationhub_id", "gencode.snapshot_date",
            "gencode.resource_title", "gencode.genome", "gencode.assembly",
            "gencode.source_filename", "gencode.source_url",
            "gencode.source_sha256", "gencode.coordinate_system",
            "gencode.assignment_method", "gencode.direct_coordinates",
            "gencode.nearest_coordinates", "gencode.unassigned_coordinates"
        ),
        Value = c(
            "TRUE", resource$release, resource$annotationHubId,
            resource$snapshotDate,
            annotationHubRecordValueDnaEpico(hub_record, "title"),
            annotationHubRecordValueDnaEpico(hub_record, "genome"),
            resource$assembly, value_or_na("source_filename"),
            annotationHubRecordValueDnaEpico(hub_record, "sourceurl"),
            value_or_na("source_sha256"), value_or_na("coordinate_system"),
            "gene_body_then_nearest_strand_specific_tss",
            as.character(counts[["annotated"]]),
            as.character(counts[["non_annotated"]]),
            as.character(counts[["unassigned"]])
        ),
        stringsAsFactors = FALSE, check.names = FALSE
    )
}
