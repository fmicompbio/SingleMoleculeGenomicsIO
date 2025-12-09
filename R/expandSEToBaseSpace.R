#' Expand the rows of a \code{RangedSummarizedExperiment} to single base resolution.
#'
#' @param se \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#'     object to be expanded to single base resolution.
#' @param region A \code{\link[GenomicRanges]{GRanges}} object with a single
#'     region defining the range for expanding \code{se}. Alternatively, the
#'     region can be specified as a character scalar (e.g. "chr1:1200-1300")
#'     that can be coerced into a \code{GRanges} object. If \code{NULL} (the
#'     default), \code{region} is set to the range of the data in \code{se}.
#' @param seqinfo \code{NULL} or a \code{\link[Seqinfo]{Seqinfo}} object
#'     containing information about the set of genomic sequences (chromosomes).
#'     Alternatively, a named numeric vector with genomic sequence names and
#'     lengths. Used to convert a character \code{region} to a \code{GRanges}
#'     object.
#' @param ignore.strand A logical scalar defining whether to ignore the strand
#'     of the \code{rowRanges(se)}.
#'
#' @importFrom SummarizedExperiment SummarizedExperiment colData assays
#' @importFrom S4Vectors metadata
#' @importFrom IRanges countOverlaps
#' @importFrom GenomicRanges GPos
#' @importFrom utils modifyList
#' @importFrom cli cli_abort
#' @importFrom BiocGenerics start end ncol colnames unstrand
#' @importFrom SparseArray NaArray nnawhich
#'
#' @return A single-base resolution \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}}
#'     corresponding to \code{se}.
#' @author Charlotte Soneson, Michael Stadler
#'
#' @export
expandSEToBaseSpace <- function(se,
                                region = NULL,
                                seqinfo = NULL,
                                keepAssays = .getReadLevelAssayNames(se),
                                ignore.strand = TRUE) {
    # check arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    if (is.character(region)) {
        region <- .regionStringToGRanges(regions = region,
                                         seqinfo = seqinfo)
    }
    if (is.null(region)) {
        region <- range(rowRanges(se), ignore.strand = TRUE)
    }
    .assertVector(x = region, type = "GRanges")
    if (countOverlaps(query = range(unstrand(rowRanges(se))),
                      subject = region, type = "within") == 0L) {
        cli_abort("Not all positions in {.arg se} are within {.arg region}")
    }
    .assertVector(x = keepAssays, type = "character",
                  validValues = .getReadLevelAssayNames(se))
    .assertScalar(x = ignore.strand, type = "logical", validValues = TRUE)

    # clean-up metadata
    md <- metadata(se)
    md$readLevelData <- modifyList(x = md$readLevelData,
                                   val = list(assayNames = keepAssays))

    # create base-resolution rowRanges
    if (ignore.strand) {
        rr <- GPos(seqnames = seqnames(region)[1],
                   pos = seq(start(region), end(region)),
                   strand = "*",
                   seqinfo = seqinfo)
    } else {
        cli_abort("{.arg ignore.strand = FALSE} is not supported yet") # nocov
    }

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

    # construct new SummarizedExperiment
    res <- SummarizedExperiment(assays = assayL,
                                rowRanges = rr,
                                colData = colData(se),
                                metadata = md)
    rownames(res) <- as.character(rr)
    return(res)
}
