#' Count pairs of modified bases by distance and modification state for bam files where modifications are indicated by sequence mismatches
#'
#' For all pairs of bases with modification calls in a read (pair), tabulate the
#' number of pairs at a given distance with a given modification state.
#'
#' @inheritParams countStatePairs
#' @inheritParams readMismatchBam
#' @param bamfile Character scalar giving the path to a \code{BAM}
#'     file, containing alignments with specific (e.g., C-to-T) mismatches.
#'     The \code{bamfile} must have an index.
#' @param sequenceContext A character scalar with an odd number of characters,
#'     specifying the genomic context on the plus strand, for which to report
#'     (mis-)matches. The string may contain IUPAC ambiguity codes, e.g.
#'     \code{sequenceContext="GCH"} would correspond to all GpC dinucleotides,
#'     excluding C's that are in CpG dinucleotides. The call will be made by
#'     comparing the genomic base in the middle of \code{sequenceContext} to
#'     the aligned base in the read and interpreted according to
#'     \code{readBaseUnmod} and \code{readBaseMod}.
#'
#' @author Charlotte Soneson, Michael Stadler
#'
#' @return A \code{DataFrame} with \code{windowSize} rows and five columns.
#'
#' @examples
#' bamfile <- system.file("extdata", "BisSeq_quasr_single.bam",
#'                        package = "SingleMoleculeGenomicsIO")
#' reffile <- system.file("extdata", "reference.fa.gz",
#'                        package = "SingleMoleculeGenomicsIO")
#' tbl <- countMismatchStatePairs(bamfile = bamfile, bamFormat = "QuasR",
#'                                regions = "chr1:6940000-6955000",
#'                                sequenceContext = "C",
#'                                readBaseUnmod = "T", readBaseMod = "C",
#'                                windowSize = 200,
#'                                sequenceReference = reffile,
#'                                BPPARAM = BiocParallel::SerialParam(),
#'                                verbose = TRUE)
#' tbl
#'
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel bpnworkers MulticoreParam
#' @importFrom cli cli_abort
#' @importFrom Biostrings IUPAC_CODE_MAP DNA_BASES vmatchPattern
#'     reverseComplement DNAString
#' @importFrom Seqinfo seqlevelsInUse seqnames
#' @importFrom IRanges resize
#' @importFrom GenomicRanges start end
#'
#' @export
countMismatchStatePairs <- function(bamfile,
                                    bamFormat = "QuasR",
                                    regions = ".",
                                    sequenceContext = "GCH",
                                    readBaseUnmod = "T",
                                    readBaseMod = "C",
                                    windowSize = 200,
                                    minMapQ = 0,
                                    minAlignedLength = 0,
                                    seqinfo = NULL,
                                    sequenceReference = NULL,
                                    BPPARAM = MulticoreParam(4L),
                                    verbose = FALSE) {
    # digest arguments
    .assertScalar(x = bamfile, type = "character")
    if (!file.exists(bamfile)) {
        cli_abort("{.arg bamfile} does not exist: {.file {bamfile}}")
    }
    .assertScalar(x = bamFormat, type = "character",
                  validValues = c("QuasR", "Bismark"))
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
    .assertScalar(x = windowSize, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = minMapQ, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minAlignedLength, type = "numeric", rngIncl = c(0, Inf))
    .assertVector(x = BPPARAM, type = "BiocParallelParam")
    .assertScalar(x = verbose, type = "logical")

    ncpuDecompression <- bpnworkers(BPPARAM)

    if (is.character(regions)) {
        regions <- regionStringToGRanges(regions = regions,
                                         seqinfo = seqinfo)
    }
    .assertVector(x = regions, type = "GRanges", allowNULL = TRUE)
    if (length(regions) == 0) {
        cli_abort("{.arg regions} must contain at least one genomic range")
    }
    seqLevelsUsed <- seqlevelsInUse(regions)
    if (!is.null(seqinfo) &&
        (!is(seqinfo, "Seqinfo") &&
         (!is.numeric(seqinfo) || is.null(names(seqinfo))))) {
        cli_abort(paste0(
            "{.arg seqinfo} must be {.code NULL}, a {.cls Seqinfo} object ",
            "or a named {.cls numeric} vector with genomic sequence lengths."))
    }
    ref <- refargToDNAStringSet(sequenceReference)

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
    regions_str <- paste0(seqnames(regions), ":", start(regions), "-", end(regions))

    resL <- read_mismatchbam_cpp(inname_str = bamfile,
                                 bam_format = bamFormat,
                                 regions = regions_str,
                                 pos_context_list = posContextList,
                                 pos_context_rev_list = posContextRevList,
                                 unmod_integer = unmodInteger,
                                 unmod_integer_rev = unmodIntegerRev,
                                 mod_integer = modInteger,
                                 mod_integer_rev = modIntegerRev,
                                 level = "read",
                                 n_alns_to_sample = 0L,
                                 tnames_for_sampling = character(0),
                                 variantRefNames = character(0),
                                 variantRefPositions = integer(0),
                                 windowSize = as.integer(windowSize),
                                 minMapQ = as.integer(minMapQ),
                                 minAlignedLength = as.integer(minAlignedLength),
                                 n_threads = as.integer(ncpuDecompression),
                                 verbose = verbose)

    res <- DataFrame(S = seq.int(windowSize),
                     unmod_unmod = resL$pair_counts[, 1],
                     unmod_mod = resL$pair_counts[, 2],
                     mod_unmod = resL$pair_counts[, 3],
                     mod_mod = resL$pair_counts[, 4])
    return(res)
}
