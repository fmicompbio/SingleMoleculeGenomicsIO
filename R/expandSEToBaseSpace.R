#' Expand the rows of a \code{RangedSummarizedExperiment} to single base resolution.
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#'     object to be expanded to single base resolution.
#' @param regions A \code{\link[GenomicRanges]{GRanges}} object with one or more
#'     regions defining the ranges for expanding \code{se}. Alternatively, the
#'     regions can be specified as character scalars (e.g. "chr1:1200-1300")
#'     that can be coerced into a \code{GRanges} object. If \code{NULL} (the
#'     default), \code{regions} is set to the range of the data in \code{se}.
#' @param seqinfo \code{NULL} or a \code{\link[Seqinfo]{Seqinfo}} object
#'     containing information about the set of genomic sequences (chromosomes).
#'     Alternatively, a named numeric vector with genomic sequence names and
#'     lengths. Used to convert a character \code{region} to a \code{GRanges}
#'     object.
#' @param keepAssays Character vector indicating which (read-level) assays to
#'     expand to base space. Only these assays will be present in the returned
#'     object.
#' @param ignore.strand A logical scalar defining whether to ignore the strand
#'     of the \code{rowRanges(se)}.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment colData assays
#' @importFrom S4Vectors metadata
#' @importFrom IRanges countOverlaps width reduce subsetByOverlaps
#' @importFrom GenomicRanges GPos
#' @importFrom utils modifyList
#' @importFrom cli cli_abort
#' @importFrom BiocGenerics start end nrow ncol colnames unstrand
#' @importFrom SparseArray NaArray nnawhich
#' @importFrom Seqinfo seqnames seqinfo
#'
#' @return A single-base resolution \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#'     corresponding to \code{se}.
#' @author Charlotte Soneson, Michael Stadler
#'
#' @examples
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "SingleMoleculeGenomicsIO")
#' se <- readModBam(bamfiles = modbamfile, regions = "chr1:6940000-6955000",
#'                  modbase = "a", verbose = TRUE,
#'                  BPPARAM = BiocParallel::SerialParam())
#' se_exp <- expandSEToBaseSpace(se)
#' dim(se)
#' dim(se_exp)
#'
#' @export
expandSEToBaseSpace <- function(se,
                                regions = NULL,
                                seqinfo = NULL,
                                keepAssays = getReadLevelAssayNames(se),
                                ignore.strand = TRUE) {
    # check arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    .assertVector(x = keepAssays, type = "character",
                  validValues = getReadLevelAssayNames(se))

    # clean-up metadata
    md <- metadata(se)
    md$readLevelData <- modifyList(x = md$readLevelData,
                                   val = list(assayNames = keepAssays))

    if (nrow(se) > 0) {
        if (is.character(regions)) {
            regions <- regionStringToGRanges(regions = regions,
                                             seqinfo = seqinfo)
        }
        if (is.null(regions)) {
            # infer regions from SE (genomic intervals covered by reads)
            regions <- reduce(do.call(c, lapply(keepAssays, function(assayName) {
                nnaIdx <- nnawhich(as.matrix(assay(se, assayName)), arr.ind = TRUE)
                nnaIdx <- as.data.frame(nnaIdx) |>
                    setNames(c("position", "read")) |>
                    mutate(seqname = as.character(seqnames(se))[.data$position],
                           position = start(se)[.data$position]) |>
                    group_by(read) |>
                    summarize(seqname = seqname[1],
                              minpos = min(position),
                              maxpos = max(position),
                              .groups = "drop")
                reduce(GRanges(seqnames = nnaIdx$seqname,
                               ranges = IRanges(start = nnaIdx$minpos,
                                                end = nnaIdx$maxpos),
                               seqinfo = seqinfo(se)))
            })))
        }
        .assertVector(x = regions, type = "GRanges")
        .assertScalar(x = ignore.strand, type = "logical", validValues = TRUE)

        # note: commenting out for now, in order to allow shrinking of the
        #       range covered by the SE as well
        # if (countOverlaps(query = range(unstrand(rowRanges(se))),
        #                   subject = region, type = "within") == 0L) {
        #     cli_abort("Not all positions in {.arg se} are within {.arg region}")
        # }

        # create base-resolution rowRanges
        if (ignore.strand) {
            rr <- GPos(seqnames = rep(seqnames(regions), width(regions)),
                       pos = unlist(lapply(seq_along(regions), function(i) {
                           seq(start(regions)[i], end(regions)[i])
                       })),
                       strand = "*",
                       seqinfo = seqinfo)
        } else {
            cli_abort("{.arg ignore.strand = FALSE} is not supported yet") # nocov
        }

        # subset SE to only the requested regions
        se <- subsetByOverlaps(se, regions, ignore.strand = TRUE)

        # iterate over assays
        newRowIdx <- match(unstrand(rowRanges(se)), rr)
        assayL <- lapply(assays(se, withDimnames = FALSE)[keepAssays], function(Df) {
            endoapply(Df, function(naArr) {
                naArrBases <- NaArray(
                    dim = c(length(rr), ncol(naArr)),
                    dimnames = list(NULL, colnames(naArr)),
                    type = "double")
                idx <- nnawhich(naArr, arr.ind = TRUE)
                idx[, 1] <- newRowIdx[idx[, 1]]
                naArrBases[idx] <- nnavals(naArr)
                return(naArrBases)
            })
        })
    } else {
        assayL <- assays(se)[keepAssays]
        rr <- GPos(seqinfo = seqinfo)
    }

    # construct new SummarizedExperiment
    suppressWarnings({
        # currently, assigning to assays triggers a deprecation warning
        # (introduced in https://github.com/Bioconductor/IRanges/commit/b4e9e7e8530a822980259c37cef186c652ba8be5)
        # see issue at https://github.com/Bioconductor/SummarizedExperiment/issues/74
        res <- SummarizedExperiment(assays = assayL,
                                    rowRanges = rr,
                                    colData = colData(se),
                                    metadata = md)
    })
    return(res)
}
