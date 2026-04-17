test_that("pileup_mismatchbam_cpp works", {
    ## helper functions --------------------------------------------------------
    .readExpectedBamHeader <- function(fname) {
        tmp <- Rsamtools::scanBamHeader(fname)[[1]]
        tmp$text <- paste0(names(tmp$text), "\t", lapply(tmp$text, paste, collapse = "\t"))
        return(tmp)
    }
    ## example data ------------------------------------------------------------
    quasr_paired_bamfile <- system.file("extdata", "BisSeq_quasr_paired.bam",
                                        package = "SingleMoleculeGenomicsIO")
    quasr_single_bamfile <- system.file("extdata", "BisSeq_quasr_single.bam",
                                        package = "SingleMoleculeGenomicsIO")
    quasr_single_indel_bamfile <- system.file("extdata", "BisSeq_quasr_single_indels.bam",
                                              package = "SingleMoleculeGenomicsIO")
    quasr_paired_discordant_bamfile <- system.file("extdata", "BisSeq_quasr_paired_discordant.bam",
                                                   package = "SingleMoleculeGenomicsIO")
    bismark_paired_bamfile <- system.file("extdata", "BisSeq_bismark_paired.bam",
                                          package = "SingleMoleculeGenomicsIO")
    true_meth <- data.frame(
        chr = "chr1",
        pos = as.integer(c(6925411, 6925417, 6925426, 6925435, 6925860,
                           6925866, 6925872, 6925875, 6925963, 6925964)),
        strand = c("+", "+", "+", "+", "-", "-", "-", "-", "+", "-"),
        quasr_count_total_paired = as.integer(c(2, 2, 2, 2, 1, 1, 2, 2, 1, 1)),
        quasr_count_meth_paired = as.integer(c(0, 0, 1, 1, 0, 0, 0, 0, 1, 1)),
        quasr_count_total_single = as.integer(c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1)),
        quasr_count_meth_single = as.integer(c(0, 0, 1, 1, 0, 0, 0, 0, 1, 1)),
        bismark_count_total_paired = as.integer(c(2, 2, 2, 2, 1, 1, 2, 2, 1, 1)),
        bismark_count_meth_paired = as.integer(c(0, 0, 1, 1, 0, 0, 0, 0, 1, 1))
    )
    posContextL <- list(chr1 = true_meth$pos[true_meth$strand == "+"] - 1L)
    posContextRevL <- list(chr1 = true_meth$pos[true_meth$strand == "-"] - 1L)
    bisseqIntegers <- c(8L, 1L, 2L, 4L)

    ## invalid arguments -------------------------------------------------------
    # ... non-existing bam file
    expect_error(pileup_mismatchbam_cpp(inname_str = "error", bam_format = "QuasR",
                                        regions = "chr1",
                                        pos_context_list = posContextL,
                                        pos_context_rev_list = posContextRevL,
                                        unmod_integer = bisseqIntegers[1],
                                        unmod_integer_rev = bisseqIntegers[2],
                                        mod_integer = bisseqIntegers[3],
                                        mod_integer_rev = bisseqIntegers[4],
                                        level = "read",
                                        n_threads = 1L, verbose = FALSE),
                 "Could not open input file")

    # ... no bam index
    tmpbam <- tempfile(fileext = ".bam")
    expect_true(file.copy(from = quasr_paired_bamfile, to = tmpbam))
    expect_error(pileup_mismatchbam_cpp(inname_str = tmpbam, bam_format = "QuasR",
                                        regions = "chr1",
                                        pos_context_list = posContextL,
                                        pos_context_rev_list = posContextRevL,
                                        unmod_integer = bisseqIntegers[1],
                                        unmod_integer_rev = bisseqIntegers[2],
                                        mod_integer = bisseqIntegers[3],
                                        mod_integer_rev = bisseqIntegers[4],
                                        level = "read",
                                        n_threads = 1L, verbose = FALSE),
                 "Failed to load the index")
    unlink(tmpbam)

    # ... requesting a region that is not contained in the bam header
    expect_error(pileup_mismatchbam_cpp(inname_str = quasr_paired_bamfile,
                                        bam_format = "QuasR",
                                        regions = "chr2",
                                        pos_context_list = posContextL,
                                        pos_context_rev_list = posContextRevL,
                                        unmod_integer = bisseqIntegers[1],
                                        unmod_integer_rev = bisseqIntegers[2],
                                        mod_integer = bisseqIntegers[3],
                                        mod_integer_rev = bisseqIntegers[4],
                                        level = "read",
                                        n_threads = 1L, verbose = FALSE),
                 "Failed to get bam iterator")

    # ... wrong bam_format
    expect_error(pileup_mismatchbam_cpp(inname_str = quasr_paired_bamfile,
                                        bam_format = "Bismark",
                                        regions = "chr1",
                                        pos_context_list = posContextL,
                                        pos_context_rev_list = posContextRevL,
                                        unmod_integer = bisseqIntegers[1],
                                        unmod_integer_rev = bisseqIntegers[2],
                                        mod_integer = bisseqIntegers[3],
                                        mod_integer_rev = bisseqIntegers[4],
                                        level = "read",
                                        n_threads = 1L, verbose = FALSE),
                 "Invalid Bismark bam format")
    expect_error(pileup_mismatchbam_cpp(inname_str = bismark_paired_bamfile,
                                        bam_format = "QuasR",
                                        regions = "chr1",
                                        pos_context_list = posContextL,
                                        pos_context_rev_list = posContextRevL,
                                        unmod_integer = bisseqIntegers[1],
                                        unmod_integer_rev = bisseqIntegers[2],
                                        mod_integer = bisseqIntegers[3],
                                        mod_integer_rev = bisseqIntegers[4],
                                        level = "read",
                                        n_threads = 1L, verbose = FALSE),
                 "Invalid QuasR bam format")

    # ... region chromosome not existing in bam header
    expect_error(pileup_mismatchbam_cpp(inname_str = quasr_paired_bamfile,
                                        bam_format = "QuasR",
                                        regions = "chr1",
                                        pos_context_list = setNames(posContextL, "chr2"),
                                        pos_context_rev_list = posContextRevL,
                                        unmod_integer = bisseqIntegers[1],
                                        unmod_integer_rev = bisseqIntegers[2],
                                        mod_integer = bisseqIntegers[3],
                                        mod_integer_rev = bisseqIntegers[4],
                                        level = "read",
                                        n_threads = 1L, verbose = FALSE),
                 "Could not find chromosome")
    expect_error(pileup_mismatchbam_cpp(inname_str = quasr_paired_bamfile,
                                        bam_format = "QuasR",
                                        regions = "chr1",
                                        pos_context_list = posContextL,
                                        pos_context_rev_list = setNames(posContextRevL, "chr2"),
                                        unmod_integer = bisseqIntegers[1],
                                        unmod_integer_rev = bisseqIntegers[2],
                                        mod_integer = bisseqIntegers[3],
                                        mod_integer_rev = bisseqIntegers[4],
                                        level = "read",
                                        n_threads = 1L, verbose = FALSE),
                 "Could not find chromosome")

    ## test that pileup_mismatchbam_cpp works by comparing to read_mismatchbam_cpp
    # helper function
    get_expected_result <- function(fname, reg, bformat = "QuasR",
                                    level = "summary", posContextL) {
        tmp0 <- read_mismatchbam_cpp(inname_str = fname, bam_format = bformat,
                                     regions = reg,
                                     n_alns_to_sample = 0,
                                     pos_context_list = posContextL,
                                     pos_context_rev_list = posContextRevL,
                                     unmod_integer = bisseqIntegers[1],
                                     unmod_integer_rev = bisseqIntegers[2],
                                     mod_integer = bisseqIntegers[3],
                                     mod_integer_rev = bisseqIntegers[4],
                                     tnames_for_sampling = "chr1",
                                     variantRefNames = character(0),
                                     variantRefPositions = integer(0),
                                     n_threads = 1, verbose = FALSE)
        tmp <- tmp0[c("chrom", "ref_position", "ref_strand", "mod_prob",
                      "read_id")]
        if (level == "summary") {
            tmp <- tmp |>
                as.data.frame() |>
                dplyr::group_by(chrom, ref_position, ref_strand) |>
                dplyr::summarise(Nmod = sum(mod_prob >= 0.5),
                                 Nvalid = dplyr::n(),
                                 .groups = "drop") |>
                dplyr::arrange(ref_position, ref_strand)
        } else {
            tmp <- tmp |>
                as.data.frame() |>
                dplyr::arrange(ref_position, read_id, dplyr::desc(ref_strand))
        }
        list(lst = c(as.list(tmp), list(bam_header = tmp0$bam_header)),
             df = tmp0$read_df)
    }

    # reading all alignments in a bam file (summary)
    reg <- "."
    res0 <- get_expected_result(quasr_single_bamfile, reg,
                                bformat = "QuasR", level = "summary",
                                posContextL = posContextL)
    # ... compare to pileup_modbam_cpp return value
    suppressMessages({
        res <- pileup_mismatchbam_cpp(inname_str = quasr_single_bamfile,
                                      bam_format = "QuasR",
                                      regions = reg,
                                      pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "summary",
                                      n_threads = 1, verbose = FALSE)
    })
    expect_type(res, "list")
    expect_length(res, 6L)
    expect_named(res, c("chrom", "ref_position", "ref_strand",
                        "Nmod", "Nvalid", "bam_header"))
    expect_identical(res0$lst, res)

    # reading all alignments in a bam file (read)
    reg <- "."
    res0 <- get_expected_result(quasr_single_bamfile, reg,
                                bformat = "QuasR", level = "read",
                                posContextL = posContextL)
    # ... compare to pileup_mismatchbam_cpp return value
    suppressMessages({
        res <- pileup_mismatchbam_cpp(inname_str = quasr_single_bamfile,
                                      bam_format = "QuasR",
                                      regions = reg,
                                      pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read",
                                      n_threads = 1, verbose = FALSE)
    })
    expect_type(res, "list")
    expect_length(res, 7L)
    expect_named(res, c("chrom", "ref_position", "ref_strand",
                        "mod_prob", "read_id", "read_df", "bam_header"))
    expect_identical(res0$df$read_id, res$read_df$read_id)
    expect_identical(names(res0$lst), names(res[c(1:5, 7)]))
    expect_identical(lengths(res0$lst), lengths(res[c(1:5, 7)]))
    res1 <- res[1:5] |>
        as.data.frame() |>
        dplyr::arrange(ref_position, read_id, dplyr::desc(ref_strand)) |>
        as.list()
    expect_identical(res0$lst, c(res1, list(bam_header = res$bam_header)))

    # reading all alignments in a bam file (read, with deletion)
    posContextLtmp <- posContextL
    posContextLtmp$chr1 <- sort(union(posContextLtmp$chr1, 6925373:6925377))
    reg <- "."
    res0 <- get_expected_result(quasr_single_indel_bamfile, reg,
                                bformat = "QuasR", level = "read",
                                posContextL = posContextLtmp)
    # ... compare to pileup_mismatchbam_cpp return value
    suppressMessages({
        res <- pileup_mismatchbam_cpp(inname_str = quasr_single_indel_bamfile,
                                      bam_format = "QuasR",
                                      regions = reg,
                                      pos_context_list = posContextLtmp,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read",
                                      n_threads = 1, verbose = FALSE)
    })
    expect_type(res, "list")
    expect_length(res, 7L)
    expect_named(res, c("chrom", "ref_position", "ref_strand",
                        "mod_prob", "read_id", "read_df", "bam_header"))
    expect_identical(res0$df$read_id, res$read_df$read_id)
    expect_identical(names(res0$lst), names(res[c(1:5, 7)]))
    expect_identical(lengths(res0$lst), lengths(res[c(1:5, 7)]))
    res1 <- res[1:5] |>
        as.data.frame() |>
        dplyr::arrange(ref_position, read_id, dplyr::desc(ref_strand)) |>
        as.list()
    expect_identical(res0$lst, c(res1, list(bam_header = res$bam_header)))

    # reading all alignments in a bam file (read, paired-end)
    reg <- "."
    res0 <- get_expected_result(quasr_paired_bamfile, reg,
                                bformat = "QuasR", level = "read",
                                posContextL = posContextL)
    # ... compare to pileup_mismatchbam_cpp return value
    suppressMessages({
        res <- pileup_mismatchbam_cpp(inname_str = quasr_paired_bamfile,
                                      bam_format = "QuasR",
                                      regions = reg,
                                      pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read",
                                      n_threads = 1, verbose = TRUE)
    })
    expect_type(res, "list")
    expect_length(res, 7L)
    expect_named(res, c("chrom", "ref_position", "ref_strand",
                        "mod_prob", "read_id", "read_df", "bam_header"))
    expect_identical(res0$df$read_id, res$read_df$read_id)
    expect_identical(names(res0$lst), names(res[c(1:5, 7)]))
    expect_true(all(lengths(res0$lst) >= lengths(res[c(1:5, 7)])))
    res1 <- res[1:5] |>
        as.data.frame() |>
        dplyr::arrange(ref_position, read_id, dplyr::desc(ref_strand)) |>
        as.list()
    expect_identical(res0$lst[names(res0$lst) != "bam_header"] |>
                         as.data.frame() |>
                         dplyr::distinct() |>
                         as.list(),
                     res1)

    # reading all alignments in a bam file (read, paired-end,
    # position with discordant call in the two mates)
    reg <- "."
    res0 <- get_expected_result(quasr_paired_discordant_bamfile, reg,
                                bformat = "QuasR", level = "read",
                                posContextL = posContextL)
    # ... compare to pileup_mismatchbam_cpp return value
    suppressMessages({
        res <- pileup_mismatchbam_cpp(inname_str = quasr_paired_discordant_bamfile,
                                      bam_format = "QuasR",
                                      regions = reg,
                                      pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read",
                                      n_threads = 1, verbose = TRUE)
    })
    expect_type(res, "list")
    expect_length(res, 7L)
    expect_named(res, c("chrom", "ref_position", "ref_strand",
                        "mod_prob", "read_id", "read_df", "bam_header"))
    expect_identical(res0$df$read_id, res$read_df$read_id)
    expect_identical(names(res0$lst), names(res[c(1:5, 7)]))
    expect_true(all(lengths(res0$lst) >= lengths(res[c(1:5, 7)])))
    res <- res[c(1,2,3,5)] |> # leave out mod_prob which will be collapsed in the R wrapper
        as.data.frame() |>
        dplyr::arrange(ref_position, read_id, dplyr::desc(ref_strand)) |>
        as.list()
    expect_identical(res0$lst[names(res0$lst) != "bam_header"] |>
                         as.data.frame() |>
                         dplyr::select(-mod_prob) |>
                         dplyr::distinct() |>
                         as.list(),
                     res)

    # reading all alignments in a bam file (read, paired-end, Bismark)
    reg <- "."
    res0 <- get_expected_result(bismark_paired_bamfile, reg,
                                bformat = "Bismark", level = "read",
                                posContextL = posContextL)
    # ... compare to pileup_mismatchbam_cpp return value
    suppressMessages({
        res <- pileup_mismatchbam_cpp(inname_str = bismark_paired_bamfile,
                                      bam_format = "Bismark",
                                      regions = reg,
                                      pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read",
                                      n_threads = 1, verbose = TRUE)
    })
    expect_type(res, "list")
    expect_length(res, 7L)
    expect_named(res, c("chrom", "ref_position", "ref_strand",
                        "mod_prob", "read_id", "read_df", "bam_header"))
    expect_identical(res0$df$read_id, res$read_df$read_id)
    expect_identical(names(res0$lst), names(res[c(1:5, 7)]))
    expect_true(all(lengths(res0$lst) >= lengths(res[c(1:5, 7)])))
    res <- res[1:5] |>
        as.data.frame() |>
        dplyr::arrange(ref_position, read_id, dplyr::desc(ref_strand)) |>
        as.list()
    expect_identical(res0$lst[names(res0$lst) != "bam_header"] |>
                         as.data.frame() |>
                         dplyr::distinct() |>
                         as.list(),
                     res)
})

