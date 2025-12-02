#' Add read segment annotation
#'
#' Add segemnts for individual reads to the \code{colData} of a
#' \code{SummarizedExperiment} object with read-level data.
#'
#' @param se A \code{\link[SummarizedExperiment]{SummarizedExperiment}} object.
#' @param irlList A list of \code{\link[IRanges]{IRangesList}} objects.
#'     Each list element corresponds to a sample (a columns in \code{se}), with
#'     the ranges in an \code{IRangesList} element giving the segments of an
#'     individual read (coordinates refer to the genome).
#' @param name A character scalar giving the name of the \code{colData(se)}
#'     column, in which the annotated segments should be stored. If \code{name}
#'     already exists in \code{colData(se)}, it will be overwritten.
#'
#' @returns A \code{SummarizedExperiment} object corresponding to \code{se},
#'     with annotated segments added to \code{colData(se)[[name]]}
#'
#' @export
#'
#' @examples
#' library(IRanges)
#' modbamfile <- system.file("extdata", "6mA_1_10reads.bam",
#'                           package = "SingleMoleculeGenomicsIO")
#' se <- readModBam(bamfiles = modbamfile, regions = "chr1:6935400-6936300",
#'                  modbase = "a", verbose = FALSE,
#'                  BPPARAM = BiocParallel::SerialParam())
#'
#' segs <- list(s1 = IRangesList(
#'     "s1-d52a5f6a-a60a-4f85-913e-eada84bfbfb9" = IRanges(
#'         start = 6925834, width = 140),
#'     "s1-6cf74134-e550-4c02-bd2b-91385422ee25" = IRanges(
#'         start = c(6926000, 6926200), width = 140)))
#'
#' se <- annotateReadSegments(se, segs, "nucl")
#' se$nucl
#' lengths(se$nucl$s1)
#'
#' @importFrom SummarizedExperiment colnames colData colData<- assay
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom cli cli_abort
#' @importFrom IRanges IRanges IRangesList
annotateReadSegments <- function(se, irlList, name) {
    # check arguments
    .assertVector(x = se, type = "RangedSummarizedExperiment")
    .assertVector(x = irlList, type = "list")
    if (!identical(colnames(se), names(irlList))) {
        cli_abort("The names of {.arg irlList} must be identical to the column names of {.arg se}")
    }
    for (i in seq_along(irlList)) {
        .assertVector(x = irlList[[i]], type = "IRangesList")
        rla <- .getReadLevelAssayNames(se)
        if (length(rla) == 0) {
            cli_abort("{.arg se} must contain at least one read-level assay")
        }
        if (!all(names(irlList[[i]]) %in% colnames(assay(se, rla[1])[[i]]))) {
            cli_abort("Not all read names in {.arg irlList[[{i}]]} exist in {.arg se} (assay {rla[1]})")
        }
    }
    .assertScalar(x = name, type = "character")

    # make sure irlList contains all reads (fill in missing ones)
    for (i in seq_along(irlList)) {
        missingReads <- setdiff(colnames(assay(se, rla[1])[[i]]),
                                names(irlList[[i]]))
        if (length(missingReads) > 0) {
            missingIR <- setNames(
                IRangesList(rep(list(IRanges()), length(missingReads))),
                missingReads)
            irlList[[i]] <- c(irlList[[i]], missingIR)
        }
        irlList[[i]] <- irlList[[i]][colnames(assay(se, rla[1])[[i]])]
    }

    # add irlList to colData(se)
    cd <- colData(se)
    cd[[name]] <- irlList
    colData(se) <- cd
    metadata(se)$readLevelData$colDataColumns <- union(
        metadata(se)$readLevelData$colDataColumns, name)

    # check validity of se
    # ### TODO

    return(se)
}
