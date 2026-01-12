## -------------------------------------------------------------------------- ##
## Checks, read_mismatchbam_cpp
## -------------------------------------------------------------------------- ##
test_that("read_mismatchbam_cpp works", {
    ## example data ------------------------------------------------------------
    quasr_paired_bamfile <- system.file("extdata", "BisSeq_quasr_paired.bam",
                                        package = "SingleMoleculeGenomicsIO")
    quasr_single_bamfile <- system.file("extdata", "BisSeq_quasr_single.bam",
                                        package = "SingleMoleculeGenomicsIO")
    quasr_single_indel_bamfile <- system.file("extdata", "BisSeq_quasr_single_indels.bam",
                                              package = "SingleMoleculeGenomicsIO")
    quasr_paired_indel_bamfile <- system.file("extdata", "BisSeq_quasr_paired_indels.bam",
                                              package = "SingleMoleculeGenomicsIO")
    quasr_paired_bamfile_inconsistentstrand <- system.file(
        "extdata", "BisSeq_quasr_paired_inconsistentstrand.bam",
        package = "SingleMoleculeGenomicsIO")
    bismark_paired_bamfile <- system.file("extdata", "BisSeq_bismark_paired.bam",
                                          package = "SingleMoleculeGenomicsIO")
    bismark_paired_bamfile_inconsistentstrand <- system.file(
        "extdata", "BisSeq_bismark_paired_inconsistentstrand.bam",
        package = "SingleMoleculeGenomicsIO")
    true_meth <- data.frame(
        chr = "chr1",
        pos = as.integer(c(6925411,6925417,6925426,6925435,6925860,6925866,6925872,6925875,6925963,6925964)),
        strand = c("+", "+", "+", "+", "-", "-", "-", "-", "+", "-"),
        quasr_count_total_paired = as.integer(c(2, 2, 2, 2, 1, 1, 2, 2, 1, 1)),
        quasr_count_meth_paired = as.integer(c(0, 0, 1, 1, 0, 0, 0, 0, 1, 1)),
        quasr_count_total_single = as.integer(c(2, 2, 2, 2, 1, 1, 1, 1, 1, 1)),
        quasr_count_meth_single = as.integer(c(0, 0, 1, 1, 0, 0, 0, 0, 1, 1)),
        bismark_count_total_paired = as.integer(c(2, 2, 2, 2, 1, 1, 2, 2, 1, 1)),
        bismark_count_meth_paired = as.integer(c(0, 0, 1, 1, 0, 0, 0, 0, 1, 1))
    )
    true_readlevel_bismark <- data.frame(
        pos = c(6925411, 6925411, 6925417, 6925417, 6925426, 6925426, 6925435,
                6925435, 6925860, 6925866, 6925872, 6925872, 6925875, 6925875,
                6925963, 6925964),
        strand = c("+", "+", "+", "+", "+", "+", "+", "+", "-", "-", "-", "-",
                   "-", "-", "+", "-"),
        read_id = c("NB501735:42:HKFLKBGX3:3:11605:15919:14063_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:1:13103:18758:18374_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:11605:15919:14063_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:1:13103:18758:18374_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:11605:15919:14063_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:1:13103:18758:18374_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:11605:15919:14063_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:1:13103:18758:18374_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:2:12308:22855:7543_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT",
                    "NB501735:53:HNNGGBGX3:2:12308:22855:7543_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:1:23103:22334:7953_1:N:0:GCCAAT",
                    "NB501735:42:HKFLKBGX3:3:22505:15192:19556_1:N:0:GCCAAT"),
        mod_prob = c(0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1)
    )
    posContextL <- list(chr1 = true_meth$pos[true_meth$strand == "+"] - 1L)
    posContextRevL <- list(chr1 = true_meth$pos[true_meth$strand == "-"] - 1L)
    bisseqIntegers <- c(8L, 1L, 2L, 4L)
    ref <- Biostrings::readDNAStringSet(system.file("extdata", "reference.fa.gz",
                                                    package = "SingleMoleculeGenomicsIO"))
    posCpG <- Biostrings::vmatchPattern(pattern = "NCG", subject = ref,
                                        max.mismatch = 0, with.indels = FALSE,
                                        fixed = "subject", algorithm = "auto")
    posCpGRev <- Biostrings::vmatchPattern(pattern = "CGN", subject = ref,
                                           max.mismatch = 0, with.indels = FALSE,
                                           fixed = "subject", algorithm = "auto")
    posCpGL <- lapply(posCpG, function(x) {
        IRanges::start(IRanges::resize(x = x, width = 1, fix = "center")) - 1L
    })
    posCpGRevL <- lapply(posCpGRev, function(x) {
        IRanges::start(IRanges::resize(x = x, width = 1, fix = "center")) - 1L
    })


    ## invalid arguments -------------------------------------------------------
    # ... non-existing bam file
    expect_error(read_mismatchbam_cpp(inname_str = "error", bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Could not open input file")

    # ... no bam index
    tmpbam <- tempfile(fileext = ".bam")
    expect_true(file.copy(from = quasr_paired_bamfile, to = tmpbam))
    expect_error(read_mismatchbam_cpp(inname_str = tmpbam, bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Failed to load the index")
    expect_error(read_mismatchbam_cpp(inname_str = tmpbam, bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE,
                                      windowSize = 40),
                 "Failed to load the index")
    expect_error(read_mismatchbam_cpp(inname_str = tmpbam, bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 30,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Failed to load the index")
    unlink(tmpbam)

    # ... wrong bam_format
    expect_error(read_mismatchbam_cpp(inname_str = quasr_paired_bamfile, bam_format = "Bismark",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Invalid Bismark bam format")
    expect_error(read_mismatchbam_cpp(inname_str = bismark_paired_bamfile, bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Invalid QuasR bam format")
    expect_error(read_mismatchbam_cpp(inname_str = quasr_paired_bamfile_inconsistentstrand,
                                      bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      windowSize = 30,
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Inconsistent strands for mates")
    expect_error(read_mismatchbam_cpp(inname_str = bismark_paired_bamfile_inconsistentstrand,
                                      bam_format = "Bismark",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      windowSize = 30,
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Inconsistent strands for mates")

    # ... context chromosome not existing in bam header
    expect_error(read_mismatchbam_cpp(inname_str = quasr_paired_bamfile, bam_format = "QuasR",
                                      regions = "chr2", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Failed to get bam iterator")
    expect_error(read_mismatchbam_cpp(inname_str = quasr_paired_bamfile, bam_format = "QuasR",
                                      regions = "chr2", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE,
                                      windowSize = 100),
                 "Failed to get bam iterator")

    # ... asking for more alignments than present in bam file
    expect_error(read_mismatchbam_cpp(inname_str = quasr_paired_bamfile, bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 3000,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Cannot sample 3000 alignments from a total of")

    # ... region chromosome not existing in bam header
    expect_error(read_mismatchbam_cpp(inname_str = quasr_paired_bamfile, bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = setNames(posContextL, "chr2"),
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Could not find chromosome")
    expect_error(read_mismatchbam_cpp(inname_str = quasr_paired_bamfile, bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = setNames(posContextRevL, "chr2"),
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Could not find chromosome")

    # ... corrupted bam file
    tmpbam <- tempfile(fileext = ".bam")
    tmpbai <- paste0(tmpbam, ".bai")
    # ... ... copy only part of `quasr_paired_bamfile`
    con_in <- file(quasr_paired_bamfile, "rb")
    data <- readBin(con_in, what = "raw", n = 1e6)
    close(con_in)
    con_out <- file(tmpbam, "wb")
    writeBin(data[seq.int(length(data) - 77)], con_out)
    close(con_out)
    expect_true(file.copy(from = paste0(quasr_paired_bamfile, ".bai"), to = tmpbai))
    expect_error(read_mismatchbam_cpp(inname_str = tmpbam, bam_format = "QuasR",
                                      regions = "chr1", pos_context_list = posContextL,
                                      pos_context_rev_list = posContextRevL,
                                      unmod_integer = bisseqIntegers[1],
                                      unmod_integer_rev = bisseqIntegers[2],
                                      mod_integer = bisseqIntegers[3],
                                      mod_integer_rev = bisseqIntegers[4],
                                      level = "read", n_alns_to_sample = 0,
                                      tnames_for_sampling = "chr1",
                                      variantRefNames = character(0),
                                      variantRefPositions = integer(0),
                                      n_threads = 1L, verbose = FALSE),
                 "Error while reading from")
    unlink(c(tmpbam, tmpbai))

    ## expected results --------------------------------------------------------
    # ... read level
    suppressMessages(expect_message(
        res1 <- read_mismatchbam_cpp(
            inname_str = quasr_paired_bamfile, bam_format = "QuasR",
            regions = "chr1:6925411-6925964", pos_context_list = posContextL,
            pos_context_rev_list = posContextRevL,
            unmod_integer = bisseqIntegers[1], unmod_integer_rev = bisseqIntegers[2],
            mod_integer = bisseqIntegers[3], mod_integer_rev = bisseqIntegers[4],
            level = "read", n_alns_to_sample = 0, tnames_for_sampling = "chr1",
            variantRefNames = character(0), variantRefPositions = integer(0),
            n_threads = 2, verbose = TRUE)
    ))
    res2 <- read_mismatchbam_cpp(
        inname_str = quasr_single_bamfile, bam_format = "QuasR",
        regions = "chr1:6925411-6925964", pos_context_list = posContextL,
        pos_context_rev_list = posContextRevL,
        unmod_integer = bisseqIntegers[1], unmod_integer_rev = bisseqIntegers[2],
        mod_integer = bisseqIntegers[3], mod_integer_rev = bisseqIntegers[4],
        level = "read", n_alns_to_sample = 0, tnames_for_sampling = "chr1",
        variantRefNames = character(0), variantRefPositions = integer(0),
        n_threads = 2, verbose = FALSE)
    res3 <- read_mismatchbam_cpp(
        inname_str = bismark_paired_bamfile, bam_format = "Bismark",
        regions = "chr1:6925411-6925964", pos_context_list = posContextL,
        pos_context_rev_list = posContextRevL,
        unmod_integer = bisseqIntegers[1], unmod_integer_rev = bisseqIntegers[2],
        mod_integer = bisseqIntegers[3], mod_integer_rev = bisseqIntegers[4],
        level = "read", n_alns_to_sample = 0, tnames_for_sampling = "chr1",
        variantRefNames = character(0), variantRefPositions = integer(0),
        n_threads = 2, verbose = FALSE)
    # ... read level (sampling)
    set.seed(1L)
    expect_warning(
        res4a <- read_mismatchbam_cpp(
            inname_str = quasr_single_bamfile, bam_format = "QuasR",
            regions = "chr1:6925411-6925964", pos_context_list = posCpGL,
            pos_context_rev_list = posCpGRevL,
            unmod_integer = bisseqIntegers[1], unmod_integer_rev = bisseqIntegers[2],
            mod_integer = bisseqIntegers[3], mod_integer_rev = bisseqIntegers[4],
            level = "read", n_alns_to_sample = 30, tnames_for_sampling = c("chr1", "error"),
            variantRefNames = character(0), variantRefPositions = integer(0),
            n_threads = 2, verbose = FALSE),
        "Ignoring unknown target name"
    )
    set.seed(1L)
    res4b <- read_mismatchbam_cpp(
        inname_str = quasr_single_bamfile, bam_format = "QuasR",
        regions = "chr1:6925411-6925964", pos_context_list = posCpGL,
        pos_context_rev_list = posCpGRevL,
        unmod_integer = bisseqIntegers[1], unmod_integer_rev = bisseqIntegers[2],
        mod_integer = bisseqIntegers[3], mod_integer_rev = bisseqIntegers[4],
        level = "read", n_alns_to_sample = 30, tnames_for_sampling = "chr1",
        variantRefNames = character(0), variantRefPositions = integer(0),
        n_threads = 2, verbose = FALSE)
    suppressMessages(expect_message(
        res4c <- read_mismatchbam_cpp(
            inname_str = quasr_single_bamfile, bam_format = "QuasR",
            regions = "chr1:6925411-6925964", pos_context_list = posCpGL,
            pos_context_rev_list = posCpGRevL,
            unmod_integer = bisseqIntegers[1], unmod_integer_rev = bisseqIntegers[2],
            mod_integer = bisseqIntegers[3], mod_integer_rev = bisseqIntegers[4],
            level = "read", n_alns_to_sample = 30, tnames_for_sampling = "chr1",
            variantRefNames = character(0), variantRefPositions = integer(0),
            n_threads = 2, verbose = TRUE)
    ))
    res5 <- read_mismatchbam_cpp(
        inname_str = quasr_single_indel_bamfile, bam_format = "QuasR",
        regions = "chr1:6925411-6925964", pos_context_list = posContextL,
        pos_context_rev_list = posContextRevL,
        unmod_integer = bisseqIntegers[1], unmod_integer_rev = bisseqIntegers[2],
        mod_integer = bisseqIntegers[3], mod_integer_rev = bisseqIntegers[4],
        level = "read", n_alns_to_sample = 0, tnames_for_sampling = "chr1",
        variantRefNames = rep("chr1", 3L), variantRefPositions = c(6925369L, 6925370L, 6925372L), # GAT
        n_threads = 2, verbose = FALSE)
    suppressMessages(expect_message(
        res6 <- read_mismatchbam_cpp(
            inname_str = quasr_paired_bamfile, bam_format = "QuasR",
            regions = "chr1:6925411-6925964", pos_context_list = posContextL,
            pos_context_rev_list = posContextRevL, windowSize = 30,
            unmod_integer = bisseqIntegers[1], unmod_integer_rev = bisseqIntegers[2],
            mod_integer = bisseqIntegers[3], mod_integer_rev = bisseqIntegers[4],
            level = "read", n_alns_to_sample = 0, tnames_for_sampling = "chr1",
            variantRefNames = character(0), variantRefPositions = integer(0),
            n_threads = 2, verbose = TRUE)
    ))
    res7 <- read_mismatchbam_cpp(
        inname_str = quasr_paired_indel_bamfile, bam_format = "QuasR",
        regions = "chr1:6925411-6925964", pos_context_list = posContextL,
        pos_context_rev_list = posContextRevL, windowSize = 30,
        unmod_integer = bisseqIntegers[1], unmod_integer_rev = bisseqIntegers[2],
        mod_integer = bisseqIntegers[3], mod_integer_rev = bisseqIntegers[4],
        level = "read", n_alns_to_sample = 0, tnames_for_sampling = "chr1",
        variantRefNames = character(0), variantRefPositions = integer(0),
        n_threads = 2, verbose = FALSE)

    # ... collect all mode 1 and mode 2 results in list
    resL <- list(res1, res2, res3, res4a, res4b, res4c, res5)

    # ... results structure
    invisible(lapply(resL, function(r) expect_type(r, "list")))

    expected_names <- c(
        "read_id", "ref_position", "chrom", "ref_strand", "qscore", "mod_prob", "read_df")
    invisible(lapply(resL, function(r) expect_named(r, expected_names)))

    expected_types <- c(
        "character", "integer", "character", "character", "double", "double", "list")
    for (i in seq_along(expected_names)) {
        invisible(lapply(resL, function(r) {
            expect_type(r[[expected_names[i]]], expected_types[i])
        }))
    }

    invisible(lapply(resL, function(r) expect_s3_class(r[["read_df"]], "data.frame")))

    expected_df_colnames <- c("read_id", "qscore", "read_length",
                              "aligned_length", "variant_label", "ref_strand")
    invisible(lapply(resL, function(r) expect_named(r$read_df, expected_df_colnames)))

    # ... content res1
    expect_identical(sort(unique(res1$ref_position)), true_meth$pos)
    expect_identical(nrow(res1$read_df), 8L)
    expect_true(any(duplicated(res1$read_df$read_id)))
    res1df <- res1[c("ref_position", "chrom", "ref_strand", "mod_prob")] |>
        as.data.frame() |>
        dplyr::group_by(chrom, ref_position, ref_strand) |>
        dplyr::summarize(count_total = dplyr::n(),
                         count_meth = as.integer(sum(mod_prob)),
                         .groups = "drop")
    expect_identical(res1df$ref_position, true_meth$pos)
    expect_identical(res1df$ref_strand, true_meth$strand)
    expect_true(all(res1df$count_total >= true_meth$quasr_count_total_paired))
    expect_identical(res1df$count_meth == 0, true_meth$quasr_count_meth_paired == 0)

    # ... content res2
    expect_identical(sort(unique(res2$ref_position)), true_meth$pos)
    expect_identical(sort(res2$ref_position),
                     rep(as.integer(true_meth$pos), true_meth$quasr_count_total_single))
    expect_true(all(res2$read_df$read_id %in% res1$read_df$read_id))
    expect_identical(nrow(res2$read_df), 4L)
    expect_true(!any(duplicated(res2$read_df$read_id)))
    res2df <- res2[c("ref_position", "chrom", "ref_strand", "mod_prob")] |>
        as.data.frame() |>
        dplyr::group_by(chrom, ref_position, ref_strand) |>
        dplyr::summarize(count_total = dplyr::n(),
                         count_meth = as.integer(sum(mod_prob)),
                         .groups = "drop")
    expect_identical(res2df$ref_position, true_meth$pos)
    expect_identical(res2df$ref_strand, true_meth$strand)
    expect_identical(res2df$count_total, true_meth$quasr_count_total_single)
    expect_identical(res2df$count_meth, true_meth$quasr_count_meth_single)

    # ... content res3
    expect_identical(sort(unique(res3$ref_position)), true_meth$pos)
    expect_identical(nrow(res3$read_df), 8L)
    expect_true(any(duplicated(res3$read_df$read_id)))
    res3df <- res3[c("ref_position", "chrom", "ref_strand", "mod_prob")] |>
        as.data.frame() |>
        dplyr::group_by(chrom, ref_position, ref_strand) |>
        dplyr::summarize(count_total = dplyr::n(),
                         count_meth = as.integer(sum(mod_prob)),
                         .groups = "drop")
    expect_identical(res3df$ref_position, true_meth$pos)
    expect_identical(res3df$ref_strand, true_meth$strand)
    expect_true(all(res3df$count_total >= true_meth$bismark_count_total_paired))
    expect_identical(res3df$count_meth == 0, true_meth$bismark_count_meth_paired == 0)
    res3reads <- res3[c("read_id", "ref_position", "chrom", "ref_strand", "mod_prob")] |>
        as.data.frame() |>
        dplyr::full_join(y = true_readlevel_bismark,
                         by = dplyr::join_by(read_id == read_id,
                                             ref_position == pos,
                                             ref_strand == strand))
    expect_true(!any(is.na(res3reads$ref_position)))
    expect_true(!any(is.na(res3reads$mod_prob.x)))
    expect_true(!any(is.na(res3reads$mod_prob.y)))
    expect_identical(res3reads$mod_prob.x, res3reads$mod_prob.y)

    # ... content of res4a, res4b and res4c
    expect_identical(res4a, res4b)
    expect_length(unique(res4a$read_id), 5L)
    expect_false(identical(res4a, res4c))
    expect_length(unique(res4c$read_id), 9L)
    expect_length(unique(res4a$read_id), nrow(res4a$read_df))
    expect_length(unique(res4c$read_id), nrow(res4c$read_df))

    # ... content of res5
    iByReadId <- split(seq_along(res5$read_id), res5$read_id)
    expect_identical(unname(lengths(iByReadId)), c(4L, 4L, 4L))
    for (j in c(2, 3)) {
        for (nm in setdiff(names(res5), c("read_id", "read_df"))) {
            expect_identical(res5[[nm]][iByReadId[[j]]], res5[[nm]][iByReadId[[1]]])
        }
    }
    expect_identical(nrow(res5$read_df), 3L)
    idx <- match(c("NB501735:42:HKFLKBGX3:3:11605:15919:14063",
                   "NB501735:42:HKFLKBGX3:3:11605:15919:14063_CINS_CSOFT_CLIP",
                   "NB501735:42:HKFLKBGX3:3:11605:15919:14063_CDEL"),
                 res5$read_df$read_id)
    expect_identical(res5$read_df$read_length[idx], c(121L, 126L, 117L))
    expect_identical(res5$read_df$aligned_length[idx], c(121L, 121L, 117L))
    expect_identical(res5$read_df$variant_label[idx], c("GAT", "GAT", "GAT"))

    # ... content of res6
    expect_type(res6, "list")
    expect_named(res6, "pair_counts")
    expect_type(res6$pair_counts, "double")
    expect_identical(dim(res6$pair_counts), c(30L, 4L))

    # ... content of res7
    expect_type(res7, "list")
    expect_named(res7, "pair_counts")
    expect_type(res7$pair_counts, "double")
    expect_identical(dim(res7$pair_counts), c(30L, 4L))
    exp <- matrix(0, nrow = 30L, ncol = 4L)
    exp[1, 1] <- 12   # 3 read pairs, 4 positions in interval, all unmethylated
    exp[7, 1] <- exp[16, 1] <- exp[19, 1] <- exp[25, 1] <- 3
    exp[10, 1] <- 6
    expect_identical(res7$pair_counts, exp)
})
