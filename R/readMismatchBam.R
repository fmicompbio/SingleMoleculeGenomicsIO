#' Read base modifications from "C-to-T" bam file(s)
#'
#' Parse read-to-genome mismatches and return a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object with
#' information on base states. Base mismatches are typically resulting from
#' single molecule genomics experiments, representing the accessibility of
#' individual bases: For bisulfite-sequencing experiments, C-to-T mismatches
#' correspond to unmethylated (inaccessible) bases (\code{readBaseUnmod="T",
#' readBaseMod="C"}, see below), while in deaminase-treatment based experiments,
#' C-to-T mismatches represent modified (accessible) bases
#' (\code{readBaseUnmod="C", readBaseMod="T"}). Secondary and supplementary
#' alignments are ignored.
#'
#' @inheritParams readModBam
#' @inheritParams addSeqContext
#' @param bamfiles Character vector with one or several paths of \code{BAM}
#'     files, containing alignments with specific (e.g., C-to-T) mismatches.
#'     If \code{bamfiles} is a named vector, the names are used
#'     as sample names and prefixes for read names. Otherwise, the prefixes will
#'     be \code{s1}, ..., \code{sN}, where \code{N} is the length of
#'     \code{bamfiles}. All \code{bamfiles} must have an index.
#' @param bamFormat A character scalar giving the format of \code{BAM} files.
#'     Currently supported ar \code{BAM} files that have been created with
#'     \code{"QuasR"} (using \code{\link[QuasR]{qAlign}} in bisulfite mode) or
#'     \code{"Bismark"}.
#' @param level Character scalar specifying the level of returned modification
#'     data. Supported values are:
#'     \describe{
#'         \item{"read"}{: Extracts modification probabilities for individual
#'             reads into an assay called \code{"mod_prob"}, with each column
#'             (sample) consisting of a position-by-read
#'             \code{\link[SparseArray]{NaMatrix}}. This is the
#'             default.}
#'         \item{"quickread"}{: Like "read", but using pileup-based data
#'             rather than parsing individual alignments. This mode does not
#'             support variant labels or read sampling. }
#'         \item{"summary"}{: Counts the total and modified bases for each
#'             position and strand and returns them in assays named
#'             \code{"Nvalid"} and \code{"Nmod"}, respectively, as well
#'             as the \code{Nmod/Nvalid} ratio in the assay \code{"FracMod"}.}
#'     }
#' @param overlapAggregation Character scalar indicating how to aggregate
#'     modification probabilities for positions that are covered by both
#'     mates in a paired-end bam file. Only used if \code{level = "read"}. In
#'     other modes, aggregation corresponding to
#'     \code{overlapAggregation = "maxQscore"} is performed automatically.
#'     Currently supported values are \code{"maxQscore"}.
#' @param sequenceContext A character scalar with an odd number of characters,
#'     specifying the genomic context on the plus strand, for which to report
#'     (mis-)matches. The string may contain IUPAC ambiguity codes, e.g.
#'     \code{sequenceContext="GCH"} would correspond to all GpC dinucleotides,
#'     excluding C's that are in CpG dinucleotides. The call will be made by
#'     comparing the genomic base in the middle of \code{sequenceContext} to
#'     the aligned base in the read and interpreted according to
#'     \code{readBaseUnmod} and \code{readBaseMod}. The genomic sequence context
#'     will be returned in \code{rowData(x)$sequenceContext}.
#' @param readBaseUnmod,readBaseMod Character scalars defining the read bases
#'     that are interpreted as unmodified or modified, respectively, when
#'     aligned to the middle base of \code{sequenceContext}.
#' @param nAlnsToSample A numeric scalar. If non-zero, \code{regions} is ignored
#'     and approximately \code{nAlnsToSample} randomly selected alignments on
#'     \code{seqnamesToSampleFrom} are read from each of the \code{bamfiles}.
#'     In order to make the results reproducible, make sure to set the
#'     \code{RNGseed} argument in the provided \code{BPPARAM} object (see
#'     below). Please note that secondary and supplementary alignments in
#'     \code{bamfiles} contribute to the total number of alignments but will not
#'     be sampled, thus the number of returned alignments may be lower than
#'     \code{nAlnsToSample}. Also, sampled reads that do not overlap
#'     a site with the indicated sequence context will not be returned, which
#'     may further reduce the number of returned alignments. Note that for
#'     paired-end bam files, individual reads are sampled and pairs will not be
#'     complete.
#'
#' @return A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object
#'     with genomic positions in rows and samples in columns. The assays
#'     depend on the value of the \code{level} argument (see above).
#'
#'
#' @examples
#' bamfile <- system.file("extdata", "BisSeq_quasr_single.bam",
#'                        package = "SingleMoleculeGenomicsIO")
#' reffile <- system.file("extdata", "reference.fa.gz",
#'                        package = "SingleMoleculeGenomicsIO")
#' se <- readMismatchBam(bamfiles = bamfile, regions = "chr1:6940000-6955000",
#'                       sequenceReference = reffile, sequenceContext = "NCG",
#'                       readBaseUnmod = "T", readBaseMod = "C",
#'                       verbose = TRUE, BPPARAM = BiocParallel::SerialParam())
#'
#' @author Michael Stadler, Charlotte Soneson
#'
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges colData
#' @importFrom SparseArray NaArray
#' @importFrom GenomicRanges GPos match start end
#' @importFrom IRanges subsetByOverlaps resize
#' @importFrom S4Vectors DataFrame SimpleList make_zero_col_DFrame
#' @importFrom Seqinfo seqnames seqlengths seqlengths<- seqlevelsInUse
#' @importFrom BiocGenerics do.call cbind pos strand sort
#' @importFrom BiocParallel bplapply MulticoreParam bpnworkers bpworkers<-
#'     bpoptions
#' @importFrom Biostrings DNAString DNA_BASES
#'     IUPAC_CODE_MAP vmatchPattern reverseComplement
#' @importFrom BSgenome getSeq
#' @importFrom methods is
#' @importFrom cli cli_abort cli_warn
#' @importFrom stats setNames
#' @importFrom dplyr group_by summarize mutate
#'
#' @export
readMismatchBam <- function(bamfiles,
                            bamFormat = "QuasR",
                            regions = NULL,
                            sequenceContext = "GCH",
                            readBaseUnmod = "T",
                            readBaseMod = "C",
                            level = "read",
                            overlapAggregation = "maxQscore",
                            sampleAnnot = NULL,
                            nAlnsToSample = 0,
                            seqnamesToSampleFrom = "chr19",
                            seqinfo = NULL,
                            sequenceReference = NULL,
                            variantPositions = NULL,
                            trim = FALSE,
                            BPPARAM = MulticoreParam(4L, RNGseed = 42L),
                            verbose = FALSE) {
    # digest arguments
    .assertVector(x = bamfiles, type = "character", rngLen = c(1, Inf))
    i <- !file.exists(bamfiles)
    if (any(i)) {
        cli_abort("not all {.arg bamfiles} exist: {.file {bamfiles[i]}}")
    }
    .assertScalar(x = bamFormat, type = "character", validValues = c("QuasR", "Bismark"))
    if (is.null(names(bamfiles))) {
        names(bamfiles) <- paste0("s", seq_along(bamfiles))
    } else if (anyDuplicated(names(bamfiles)) > 0L) {
        cli_abort("{.code names(bamfiles)} are not unique")
    }
    .assertScalar(x = sequenceContext, type = "character")
    sequenceContext <- toupper(sequenceContext)
    if (identical(nchar(sequenceContext) %% 2, 0)) {
        cli_abort("{.arg sequenceContext} ({sequenceContext}) needs to have an odd number of characters")
    }
    if (!substr(sequenceContext,
                ceiling(nchar(sequenceContext) / 2),
                ceiling(nchar(sequenceContext) / 2)) %in% DNA_BASES) {
        cli_abort("The central base of {.arg sequenceContext} must be A, C, G or T")
    }
    .assertScalar(x = readBaseUnmod, type = "character", validValues = DNA_BASES)
    .assertScalar(x = readBaseMod, type = "character", validValues = DNA_BASES)
    unmodBases <- strsplit(IUPAC_CODE_MAP[readBaseUnmod], "")[[1]]
    modBases <- strsplit(IUPAC_CODE_MAP[readBaseMod], "")[[1]]
    if (length(intersect(unmodBases, modBases)) > 0L) {
        cli_abort("{.arg readBaseUnmod} and {.arg readBaseMod} must not include common bases")
    }
    unmodInteger <- as.integer(sum(c(A = 1L, C = 2L, G = 4L, T = 8L)[unmodBases]))
    unmodIntegerRev <- as.integer(sum(c(A = 8L, C = 4L, G = 2L, T = 1L)[unmodBases]))
    modInteger <- as.integer(sum(c(A = 1L, C = 2L, G = 4L, T = 8L)[modBases]))
    modIntegerRev <- as.integer(sum(c(A = 8L, C = 4L, G = 2L, T = 1L)[modBases]))
    .assertScalar(x = level, type = "character",
                  validValues = c("read", "summary", "quickread"))
    .assertVector(x = sampleAnnot, type = "data.frame", allowNULL = TRUE)
    if (!is.null(sampleAnnot)) {
        if (!("sample" %in% colnames(sampleAnnot))) {
            cli_abort("{.arg sampleAnnot} must have at least a column named 'sample'")
        }
        if (!all(names(bamfiles) %in% sampleAnnot$sample)) {
            cli_abort(paste0(
                "Annotation information missing for some samples: ",
                "{setdiff(names(bamfiles), sampleAnnot$sample)}"))
        }
    }
    if (is.character(regions)) {
        regions <- .regionStringToGRanges(regions = regions,
                                          seqinfo = seqinfo)
    }
    .assertVector(x = regions, type = "GRanges", allowNULL = TRUE)
    .assertScalar(x = nAlnsToSample, type = "numeric", rngIncl = c(0, Inf))
    if (nAlnsToSample > 0 && level %in% c("summary")) {
        cli_abort(paste0("Read sampling is not supported if {.arg level} is set ",
                         "to 'summary'"))
    }
    if (nAlnsToSample > 0) {
        .assertVector(x = seqnamesToSampleFrom, type = "character")
        seqLevelsUsed <- seqnamesToSampleFrom
        if (length(regions) > 0) {
            cli_warn("Ignoring {.arg regions} because {.arg nAlnsToSample} is greater than zero")
        }
        regions <- GRanges()
        if (!is.null(variantPositions)) {
            cli_warn("Ignoring {.arg variantPositions} because {.arg nAlnsToSample} is greater than zero")
        }
        variantPositions <- NULL
    } else {
        if (length(regions) == 0) {
            cli_abort("{.arg regions} must contain at least one genomic range if not in sampling mode")
        }
        seqLevelsUsed <- seqlevelsInUse(regions)
    }
    if (!is.null(seqinfo) &&
        (!is(seqinfo, "Seqinfo") &&
         (!is.numeric(seqinfo) || is.null(names(seqinfo))))) {
        cli_abort(paste0(
            "{.arg seqinfo} must be {.code NULL}, a {.cls Seqinfo} object ",
            "or a named {.cls numeric} vector with genomic sequence lengths."))
    }
    ref <- refargToDNAStringSet(sequenceReference)
    .assertVector(x = variantPositions, type = "GPos", allowNULL = TRUE)
    .assertScalar(x = trim, type = "logical")
    .assertVector(x = BPPARAM, type = "BiocParallelParam")
    .assertScalar(x = verbose, type = "logical")

    # sort and subset variantPositions
    if (length(variantPositions) > 0) {
        variantPositions <- sort(
            sort(subsetByOverlaps(x = variantPositions,
                                  ranges = regions,
                                  ignore.strand = TRUE)),
            ignore.strand = TRUE)
        variantRefNames <- as.character(seqnames(variantPositions))
        # make coordinates zero-based
        variantRefPositions <- pos(variantPositions) - 1L
    } else {
        variantRefNames <- character(0L)
        variantRefPositions <- integer(0L)
    }
    # TODO: lift out as a helper function?

    # determine the number of parallel threads to be used for
    # bam files (preferred) and decompression of bam records (if available)
    # (accept some level of over-subscription)
    ncpuTotal <- bpnworkers(BPPARAM)
    if (is(BPPARAM, "MulticoreParam") || is(BPPARAM, "SnowParam")) {
        ncpuFiles <- min(ncpuTotal, length(bamfiles))
        oversubscriptionRate <- 2.0
        ncpuDecompression <- min(8L, max(1L, as.integer(
            floor(oversubscriptionRate * ncpuTotal / ncpuFiles))))
        bpworkers(BPPARAM) <- ncpuFiles
        on.exit(bpworkers(BPPARAM) <- ncpuTotal)
    } else {
        ncpuDecompression <- 1L
    }
    # TODO: lift out as a helper function?

    # obtain reference sequences
    .message("finding positions with {sequenceContext}")
    seqLevelsUsed <- intersect(seqLevelsUsed, names(ref))
    ref <- ref[seqLevelsUsed]

    # identify positions with `sequenceContext`
    posContext <- vmatchPattern(
        pattern = sequenceContext,
        subject = ref,
        max.mismatch = 0,
        with.indels = FALSE,
        fixed = "subject",
        algorithm = "auto")
    posContextRev <- vmatchPattern(
        pattern = as.character(reverseComplement(DNAString(sequenceContext))),
        subject = ref,
        max.mismatch = 0,
        with.indels = FALSE,
        fixed = "subject",
        algorithm = "auto")

    # convert to list of zero-based indices for each chromosome to use in C++
    posContextList <- lapply(posContext, function(x) {
        start(resize(x = x, width = 1, fix = "center")) - 1L
    })
    posContextRevList <- lapply(posContextRev, function(x) {
        start(resize(x = x, width = 1, fix = "center")) - 1L
    })

    # extract modification states from `bamfiles`
    .message("extracting base mismatches from BAM files", noTimer = TRUE)
    # remark: unsafe to use as.character(GRanges) here, as a length-1 range
    #         would become e.g. "chr1:35000" (no end coordinate), which is
    #         interpreted by htslib as: "read all alignments overlapping chr1:35000-END_OF_chr1"
    regions_str <- paste0(seqnames(regions), ":", start(regions), "-", end(regions))
    resLL <- bplapply(
        setNames(names(bamfiles), names(bamfiles)),
        function(nm,
                 mylevel = level,
                 bamf = bamfiles[nm],
                 mybamFormat = bamFormat,
                 myregions_str = regions_str,
                 myposContextList = posContextList,
                 myposContextRevList = posContextRevList,
                 myunmodInteger = unmodInteger,
                 myunmodIntegerRev = unmodIntegerRev,
                 mymodInteger = modInteger,
                 mymodIntegerRev = modIntegerRev,
                 mynAlnsToSample = nAlnsToSample,
                 myseqnamesToSampleFrom = seqnamesToSampleFrom,
                 myvariantRefNames = variantRefNames,
                 myvariantRefPositions = variantRefPositions,
                 myncpuDecompression = ncpuDecompression,
                 myverbose = if (ncpuTotal > 1) FALSE else verbose) {

            if (mylevel == "read") {
                # extract base mismatches
                resL <- read_mismatchbam_cpp(
                    inname_str = bamf,
                    bam_format = mybamFormat,
                    regions = myregions_str,
                    pos_context_list = myposContextList,
                    pos_context_rev_list = myposContextRevList,
                    unmod_integer = myunmodInteger,
                    unmod_integer_rev = myunmodIntegerRev,
                    mod_integer = mymodInteger,
                    mod_integer_rev = mymodIntegerRev,
                    level = mylevel,
                    n_alns_to_sample = as.integer(mynAlnsToSample),
                    tnames_for_sampling = myseqnamesToSampleFrom,
                    variantRefNames = myvariantRefNames,
                    variantRefPositions = as.integer(myvariantRefPositions),
                    windowSize = 0L,
                    n_threads = as.integer(myncpuDecompression),
                    verbose = myverbose)
            } else if (mylevel == "summary") {
                resL <- pileup_mismatchbam_cpp(
                    inname_str = bamf,
                    bam_format = mybamFormat,
                    regions = myregions_str,
                    pos_context_list = myposContextList,
                    pos_context_rev_list = myposContextRevList,
                    unmod_integer = myunmodInteger,
                    unmod_integer_rev = myunmodIntegerRev,
                    mod_integer = mymodInteger,
                    mod_integer_rev = mymodIntegerRev,
                    level = "summary",
                    n_threads = as.integer(myncpuDecompression),
                    verbose = myverbose
                )
            } else if (mylevel == "quickread") {
                resL <- pileup_mismatchbam_cpp(
                    inname_str = bamf,
                    bam_format = mybamFormat,
                    regions = myregions_str,
                    pos_context_list = myposContextList,
                    pos_context_rev_list = myposContextRevList,
                    unmod_integer = myunmodInteger,
                    unmod_integer_rev = myunmodIntegerRev,
                    mod_integer = mymodInteger,
                    mod_integer_rev = mymodIntegerRev,
                    level = "read",
                    n_threads = as.integer(myncpuDecompression),
                    verbose = myverbose
                )
            }
            resL
        }, BPPARAM = BPPARAM, BPOPTIONS = bpoptions(
            progressbar = (verbose && ncpuTotal > 1)))

    # if level = "read" (results from read_mismatchbam_cpp), resolve overlapping
    # parts of reads
    if (level == "read" && overlapAggregation == "maxQscore") {
        resLL <- lapply(resLL, function(resL) {
            iByReadPos <- split(seq_along(resL$read_id),
                                paste0(resL$read_id, resL$chrom, resL$ref_position,
                                       resL$ref_strand))
            if (any(lengths(iByReadPos) > 1)) {
                resL$read_id <- unlist(unname(lapply(
                    iByReadPos, function(i) unique(resL$read_id[i]))))
                resL$ref_position <- unlist(unname(lapply(
                    iByReadPos, function(i) unique(resL$ref_position[i]))))
                resL$chrom <- unlist(unname(lapply(
                    iByReadPos, function(i) unique(resL$chrom[i]))))
                resL$ref_strand <- unlist(unname(lapply(
                    iByReadPos, function(i) unique(resL$ref_strand[i]))))
                resL$mod_prob <- unlist(unname(lapply(
                    iByReadPos, function(i) resL$mod_prob[i[which.max(resL$qscore[i])]])))
                resL$qscore <- unlist(unname(lapply(
                    iByReadPos, function(i) max(resL$qscore[i]))))
            }
            if (any(duplicated(resL$read_df$read_id))) {
                resL$read_df <- resL$read_df |>
                    mutate(read_id = factor(.data$read_id, levels = unique(.data$read_id))) |>
                    group_by(.data$read_id) |>
                    summarize(qscore = mean(.data$qscore),
                              read_length = sum(.data$read_length),
                              aligned_length = sum(.data$aligned_length),
                              variant_label = ifelse(any(is.na(.data$variant_label)), NA_character_,
                                                     paste(.data$variant_label, collapse = "")),
                              ref_strand = paste(unique(.data$ref_strand, collapse = "/")),
                              .groups = "drop") |>
                    mutate(read_id = as.character(.data$read_id)) |>
                    as.data.frame()
                rownames(resL$read_df) <- seq_len(nrow(resL$read_df))
            }
            resL
        })
    }

    # create GPos objects for each input
    gposL <- bplapply(resLL, function(resL, myseqinfo = seqinfo) {
        GenomicRanges::GPos(seqnames = resL$chrom, pos = resL$ref_position,
                            strand = resL$ref_strand,
                            seqinfo = myseqinfo)
    }, BPPARAM = BPPARAM)

    # create combined GPos, reduce to unique positions
    .message("finding unique genomic positions...")
    gpos <- sort(sort(unique(do.call(c, unname(gposL)))), ignore.strand = TRUE)
    .message("collapsed {sum(lengths(gposL))} positions to {length(gpos)} unique ones")

    # if trim=TRUE, trim GPos to only the indicated region
    if (trim) {
        gpos <- subsetByOverlaps(gpos, regions, ignore.strand = TRUE)
    }

    # add sequence context
    .message("extracting sequence contexts")
    mcols(gpos)$sequenceContext <- extractSeqContext(
        x = as(gpos, "GRanges"),
        sequenceContextWidth = nchar(sequenceContext),
        sequenceReference = ref)

    if (level %in% c("read", "quickread")) {
        # extract unique read names
        readL <- lapply(resLL, function(resL) resL$read_df$read_id)

        # modified probability
        modmat <- make_zero_col_DFrame(nrow = length(gpos))
        readdfL <- SimpleList()
        for (nm in names(bamfiles)) {
            x <- resLL[[nm]]
            if (length(x$read_id) > 0) {
                namat <- NaArray(dim = c(length(gpos), length(readL[[nm]])),
                                 dimnames = list(NULL, paste0(nm, "-", readL[[nm]])),
                                 type = "double")
                i <- match(gposL[[nm]], gpos)
                # if trim=TRUE, not all positions in gposL may be present in gpos
                found <- which(!is.na(i))
                j <- match(x$read_id, readL[[nm]])
                namat[cbind(i[found], j[found])] <- x$mod_prob[found]
                modmat[[nm]] <- namat
                rownames(x$read_df) <- paste0(nm, "-", x$read_df$read_id)
                x$read_df$read_id <- NULL
                x$read_df$aligned_fraction <- x$read_df$aligned_length / x$read_df$read_length
                readdfL[[nm]] <- DataFrame(x$read_df)
            } else {
                modmat[[nm]] <- NaArray(dim = c(length(gpos), 0), type = "double")
                readdfL[[nm]] <- DataFrame(qscore = numeric(0),
                                           read_length = integer(0),
                                           aligned_length = integer(0),
                                           variant_label = character(0),
                                           ref_strand = character(0),
                                           aligned_fraction = numeric(0))
            }
        }
    } else {
        Nmod <- matrix(0, nrow = length(gpos), ncol = length(bamfiles))
        Nvalid <- matrix(0, nrow = length(gpos), ncol = length(bamfiles))
        for (i in seq_along(bamfiles)) {
            nm <- names(bamfiles)[i]
            # if trim=TRUE, not all positions in gposL may be present in gpos
            j <- match(gposL[[nm]], gpos)
            found <- which(!is.na(j))
            Nmod[j[found], i] <- resLL[[nm]]$Nmod[found]
            Nvalid[j[found], i] <- resLL[[nm]]$Nvalid[found]
        }
        FracMod <- Nmod / Nvalid
    }

    # create SummarizedExperiment object
    cdata <- DataFrame(
        row.names = names(bamfiles),
        sample = names(bamfiles)
    )
    if (level %in% c("read", "quickread")) {
        cdata <- cbind(
            cdata,
            DataFrame(n_reads = unlist(lapply(readdfL, nrow), use.names = FALSE),
                      readInfo = readdfL)
        )
    }
    if (!is.null(sampleAnnot) && any(colnames(sampleAnnot) != "sample")) {
        sampleAnnot <- sampleAnnot[match(cdata$sample, sampleAnnot$sample),
                                   colnames(sampleAnnot) != "sample",
                                   drop = FALSE]
        cdata <- cbind(cdata, sampleAnnot)
    }
    if (level %in% c("read", "quickread")) {
        stopifnot(names(modmat) == cdata$sample)
        se <- SummarizedExperiment(
            assays = list(mod_prob = modmat),
            rowRanges = gpos,
            colData = cdata,
            metadata = list(readLevelData = list(assayNames = "mod_prob",
                                                 colDataColumns = "readInfo"),
                            variantPositions = variantPositions,
                            readBaseMod = readBaseMod,
                            readBaseUnmod = readBaseUnmod)
        )
    } else {
        stopifnot(colnames(Nmod) == cdata$sample,
                  colnames(Nvalid) == cdata$sample,
                  colnames(FracMod) == cdata$sample)
        se <- SummarizedExperiment(
            assays = list(Nmod = Nmod,
                          Nvalid = Nvalid,
                          FracMod = FracMod),
            rowRanges = gpos,
            colData = cdata,
            metadata = list(readLevelData = list(assayNames = character(0),
                                                 colDataColumns = character(0)),
                            readBaseMod = readBaseMod,
                            readBaseUnmod = readBaseUnmod)
        )
    }
    if (nrow(se) > 0) {
        rownames(se) <- paste0(
            seqnames(rowRanges(se)), ":", pos(rowRanges(se)), ":",
            strand(rowRanges(se)))
        colnames(se) <- rownames(colData(se))
    }

    # Remove reads with all NA values
    if (level %in% c("read", "quickread")) {
        se <- filterReads(se, readInfoCol = NULL, qcCol = NULL, prune = FALSE)
    }
    se
}
