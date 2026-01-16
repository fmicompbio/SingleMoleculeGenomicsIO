## -------------------------------------------------------------------------- ##
## Checks, concatenate_files
## -------------------------------------------------------------------------- ##
test_that("concatenate_files works", {
    infiles <- tempfile(pattern = paste0("file", seq.int(3)))
    outfile <- tempfile()

    expect_error(concatenate_files(), "missing")
    expect_error(concatenate_files("ERROR"), "missing")
    expect_error(concatenate_files(infiles, "error/error/error"))
    expect_error(concatenate_files(1L, "ERROR"), "string vector")
    expect_error(concatenate_files("in", c("out1", "out2")), "single string")
    expect_error(concatenate_files(infiles, outfile), "Could not open")

    set.seed(1L)
    data <- unlist(lapply(seq_along(infiles), function(i) {
        paste(sample(letters, 10), collapse = "")
    }))

    for (i in seq_along(data)) {
        writeLines(text = data[i], con = infiles[i])
    }

    oL <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
    res <- lapply(oL, function(o) {
        expect_identical(concatenate_files(infiles[o], outfile), outfile)
        expect_identical(data[o], readLines(outfile))
    })

    unlink(c(infiles, outfile))
})

## -------------------------------------------------------------------------- ##
## Checks, concatenate_hts_files
## -------------------------------------------------------------------------- ##
test_that("concatenate_hts_files works", {
    bamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                         "6mA_2_10reads.bam"),
                            package = "SingleMoleculeGenomicsIO")
    samfiles <- tempfile(pattern = paste0("file", seq_along(bamfiles)),
                         fileext = ".sam")
    outsam <- tempfile(fileext = ".sam")
    outbam <- tempfile(fileext = ".bam")
    expect_identical(
        unlist(lapply(seq_along(bamfiles), function(i) {
            Rsamtools::asSam(file = bamfiles[i],
                             destination = sub(".sam$", "", samfiles[i]))
        })),
        samfiles)

    # unknown output extension
    expect_error(
        concatenate_hts_files(bamfiles, sub(".bam$", ".error", outbam), 2L),
        "Unknown .output_file. extension"
    )

    sam_fields <- c("qname", "flag", "rname", "strand", "pos")

    # in: bam, out: bam
    expect_identical(
        concatenate_hts_files(bamfiles, outbam, 2L),
        outbam)
    res <- Rsamtools::scanBam(
        file = outbam, param = Rsamtools::ScanBamParam(what = sam_fields))
    expect_type(res, "list")
    expect_named(res[[1]], sam_fields)
    expect_identical(
        lengths(res[[1]]),
        stats::setNames(rep(20L, length(sam_fields)), sam_fields))
    unlink(outbam)

    # in: bam, out: sam
    expect_identical(
        concatenate_hts_files(bamfiles, outsam, 2L),
        outsam)
    tmp <- readLines(outsam)
    expect_identical(sum(grepl("^@", tmp)), 67L)
    expect_length(tmp, 87L)
    unlink(outsam)

    # in: sam, out: bam
    expect_identical(
        concatenate_hts_files(samfiles, outbam, 2L),
        outbam)
    res <- Rsamtools::scanBam(
        file = outbam, param = Rsamtools::ScanBamParam(what = sam_fields))
    expect_type(res, "list")
    expect_named(res[[1]], sam_fields)
    expect_identical(
        lengths(res[[1]]),
        stats::setNames(rep(20L, length(sam_fields)), sam_fields))
    unlink(outbam)

    # in: sam, out: sam
    expect_identical(
        concatenate_hts_files(samfiles, outsam, 2L),
        outsam)
    tmp <- readLines(outsam)
    expect_identical(sum(grepl("^@", tmp)), 67L)
    expect_length(tmp, 87L)
    unlink(outsam)

    unlink(samfiles)
})

## -------------------------------------------------------------------------- ##
## Checks, getChromosomeNamesFromBam
## -------------------------------------------------------------------------- ##
test_that("getChromosomeNamesFromBam works", {
    bamfiles <- system.file("extdata", c("6mA_1_10reads.bam",
                                         "6mA_2_10reads.bam"),
                            package = "SingleMoleculeGenomicsIO")
    bamfileNoheader <- tempfile(fileext = ".bam")
    res <- filter_modbam_cpp(bamfiles[1], bamfileNoheader, modbase = "a",
                             includeHeader = FALSE)
    expect_identical(res[["retained"]], 10)

    expect_error(getChromosomeNamesFromBam(), "missing")
    expect_error(getChromosomeNamesFromBam(1L), "single string")
    expect_error(getChromosomeNamesFromBam(bamfiles), "single string")
    expect_error(getChromosomeNamesFromBam(bamfileNoheader), "Could not open")

    expected_chrs <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
                       "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
                       "chr16", "chr17", "chr18", "chr19", "chrX", "chrY", "chrM", "GL456210.1",
                       "GL456211.1", "GL456212.1", "GL456219.1", "GL456221.1", "GL456233.2",
                       "GL456239.1", "GL456354.1", "GL456359.1", "GL456360.1", "GL456366.1",
                       "GL456367.1", "GL456368.1", "GL456370.1", "GL456372.1", "GL456378.1",
                       "GL456379.1", "GL456381.1", "GL456382.1", "GL456383.1", "GL456385.1",
                       "GL456387.1", "GL456389.1", "GL456390.1", "GL456392.1", "GL456394.1",
                       "GL456396.1", "JH584295.1", "JH584296.1", "JH584297.1", "JH584298.1",
                       "JH584299.1", "JH584300.1", "JH584301.1", "JH584302.1", "JH584303.1",
                       "JH584304.1", "MU069434.1", "MU069435.1")
    expect_identical(getChromosomeNamesFromBam(bamfiles[1]), expected_chrs)
    expect_identical(getChromosomeNamesFromBam(bamfiles[2]), expected_chrs)

    unlink(bamfileNoheader)
})

## -------------------------------------------------------------------------- ##
## Checks, .assertScalar
## -------------------------------------------------------------------------- ##
test_that(".assertScalar works", {
    expect_error(.assertScalar(1, type = TRUE))
    expect_error(.assertScalar(1, type = 1))
    expect_error(.assertScalar(1, type = c("numeric", "character")))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = TRUE))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = "rng"))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = 1))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = 1:3))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = TRUE))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = "rng"))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = 1))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = 1:3))
    expect_error(.assertScalar(1, type = "numeric", rngIncl = c(0, 2), rngExcl = c(0, 2)))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = 1))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = "rng"))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = NULL))
    expect_error(.assertScalar(1, type = "numeric", allowNULL = c(TRUE, FALSE)))

    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(1, 3)))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(1, 3)))
    expect_true(.assertScalar(1, type = "numeric", rngExcl = c(1, 3), validValues = 1))
    expect_true(.assertScalar(-1, type = "numeric", rngIncl = c(1, 3), validValues = c(-1, 0)))
    expect_error(.assertScalar(-1, type = "numeric", rngIncl = c(1, 3), validValues = 0))
    expect_true(.assertScalar(-1, type = "numeric", validValues = c(-1, 0)))
    expect_error(.assertScalar(-1, type = "numeric", validValues = c(-2, 0)))
    expect_true(.assertScalar(NA_real_, type = "numeric", rngIncl = c(1, 2), validValues = NA_real_))
    expect_error(.assertScalar(NA, type = "numeric", rngIncl = c(1, 2), validValues = NA_real_))
    expect_true(.assertScalar(NA_real_, type = "numeric", rngIncl = c(1, 2), validValues = NA))
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(0, 3), validValues = 3))
    expect_true(.assertScalar(1, rngIncl = c(0, 3), validValues = 3))
    expect_true(.assertScalar(1, type = "numeric", rngIncl = c(0, 1)))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(0, 1)))
    expect_true(.assertScalar(1, type = "numeric", rngExcl = c(0, 1), validValues = 1))
    expect_error(.assertScalar(1, type = "numeric", rngExcl = c(0, 1), validValues = 3:4))
    expect_true(.assertScalar(NULL, type = "numeric", allowNULL = TRUE))
    expect_error(.assertScalar(NULL, type = "numeric", allowNULL = FALSE))
    expect_error(.assertScalar(1, type = "character"))
    expect_error(.assertScalar("x", type = "numeric"))
    expect_error(.assertScalar(FALSE, type = "character"))
    expect_error(.assertScalar(c(1, 2), type = "numeric"))
    test <- "text"
    expect_error(.assertScalar(x = test, type = "numeric"),
                 ".test. must be of class .numeric.")
    expect_error(.assertScalar(x = list(a = 1)$a, type = "logical"),
                 "list(a = 1)$a", fixed = TRUE)
    tmp <- matrix(1:4, ncol = 2)
    expect_error(.assertScalar(x = tmp[, 1], type = "logical"),
                 "tmp[, 1]", fixed = TRUE)
})

## -------------------------------------------------------------------------- ##
## Checks, .assertVector
## -------------------------------------------------------------------------- ##
test_that(".assertVector works", {
    expect_error(.assertVector(1, type = TRUE))
    expect_error(.assertVector(1, type = 1))
    expect_error(.assertVector(1, type = c("numeric", "character")))
    expect_error(.assertVector(1, type = "numeric", rngIncl = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngIncl = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngIncl = 1))
    expect_error(.assertVector(1, type = "numeric", rngIncl = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngExcl = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngExcl = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngExcl = 1))
    expect_error(.assertVector(1, type = "numeric", rngExcl = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngIncl = c(0, 2), rngExcl = c(0, 2)))
    expect_error(.assertVector(1, type = "numeric", allowNULL = 1))
    expect_error(.assertVector(1, type = "numeric", allowNULL = "rng"))
    expect_error(.assertVector(1, type = "numeric", allowNULL = NULL))
    expect_error(.assertVector(1, type = "numeric", allowNULL = c(TRUE, FALSE)))
    expect_error(.assertVector(1, type = "numeric", len = TRUE))
    expect_error(.assertVector(1, type = "numeric", len = "rng"))
    expect_error(.assertVector(1, type = "numeric", len = 1:3))
    expect_error(.assertVector(1, type = "numeric", rngLen = TRUE))
    expect_error(.assertVector(1, type = "numeric", rngLen = "rng"))
    expect_error(.assertVector(1, type = "numeric", rngLen = 1))
    expect_error(.assertVector(1, type = "numeric", rngLen = 1:3))

    expect_true(.assertVector(c(1, 2), type = "numeric", rngIncl = c(1, 3)))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngIncl = c(1, 1.5)))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngExcl = c(1, 3)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngExcl = c(1, 3), validValues = 1))
    expect_error(.assertVector(c(1, 2), type = "numeric", validValues = c(1, 3)))
    expect_true(.assertVector(c(1, 2), type = "numeric", validValues = c(1, 2)))
    expect_error(.assertVector(c(1, 2), type = "numeric", len = 1))
    expect_true(.assertVector(c(1, 2), type = "numeric", len = 2))
    expect_error(.assertVector(c(1, 2), type = "numeric", rngLen = c(3, 5)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngLen = c(2, 5)))
    expect_true(.assertVector(c(1, 2), type = "numeric", rngLen = c(1, 2)))
    expect_error(.assertVector(c("a", "b"), type = "character", validValues = c("A", "B")))
    expect_true(.assertVector(LETTERS[1:2], type = "character", validValues = LETTERS))
    test <- "text"
    expect_error(.assertVector(x = test, type = "numeric"),
                 ".test. must be of class .numeric.")
})

## -------------------------------------------------------------------------- ##
## Checks, .assertPackagesAvailable
## -------------------------------------------------------------------------- ##
test_that(".assertPackagesAvailable works", {
    testfunc <- function(...) .assertPackagesAvailable(...)
    expect_error(testfunc(1L))
    expect_error(testfunc("test", "error"))
    expect_error(testfunc("test", c(TRUE, FALSE)))

    expect_true(testfunc("base"))
    expect_true(testfunc("githubuser/base"))
    expect_true(testfunc(c("base", "methods")))
    expect_error(testfunc(c("error", "error2")), "BiocManager")
    expect_error(testfunc("error1", suggestInstallation = FALSE),
                 "installed[.]")
    rm(testfunc)
})

## -------------------------------------------------------------------------- ##
## Checks, .assertValidModbase
## -------------------------------------------------------------------------- ##
test_that(".assertValidModbase works", {
    good <- c("m", "h", "f", "c", "C", "g", "e", "b", "T", "U", "a", "A", "o", "G", "n", "N")
    bad <- setdiff(c(letters, LETTERS), good)
    expect_true(.assertValidModbase(good))
    for (modbase in good) {
        expect_true(.assertValidModbase(modbase))
    }
    expect_error(.assertValidModbase(bad))
    for (modbase in bad) {
        expect_error(.assertValidModbase(modbase))
    }
})

## -------------------------------------------------------------------------- ##
## Checks, regionStringToGRanges
## -------------------------------------------------------------------------- ##
test_that("regionStringToGRanges works", {
    # supported formats:
    # "REF"
    # "REF:"
    # "REF:START"
    # "REF:-END"
    # "REF:START-END"
    # "."

    expect_error(regionStringToGRanges("chr1", "error"),
                 "or a named .numeric. vector")
    expect_error(regionStringToGRanges(c(".")),
                 ".seqinfo. argument is required")
    expect_error(regionStringToGRanges(c(".", "chr1"), c("chr1" = 100)),
                 "can only be given as a single region")
    expect_error(regionStringToGRanges(c("chr1:1-10:+", "chr1:-", "")),
                 "unrecognized format in 2 regions")

    slens <- c(chr1 = 100, chr2 = 200, chr4 = 400)
    si <- Seqinfo::Seqinfo(seqnames = names(slens), seqlengths = unname(slens))
    reg1 <- c("chr1", "chr1:", "chr1:-",
              "chr2:10", "chr2:10-",
              "chr3:-100",
              "chr4:20-70")

    gr1 <- regionStringToGRanges(regions = reg1, seqinfo = slens)
    gr2 <- regionStringToGRanges(regions = reg1, seqinfo = si)
    gr3 <- regionStringToGRanges(regions = reg1, seqinfo = NULL)

    intmax <- .Machine$integer.max

    expect_identical(gr1, GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr3", "chr4"),
        ranges = IRanges::IRanges(start = c(1, 1, 1, 10, 10, 1, 20),
                                  end = c(100, 100, 100, 200, 200, 100, 70)),
        seqlengths = c(slens, c(chr3 = intmax))[paste0("chr", 1:4)]
    ))
    expect_identical(gr1, gr2)
    expect_identical(gr3, GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1", "chr1", "chr2", "chr2", "chr3", "chr4"),
        ranges = IRanges::IRanges(start = c(1, 1, 1, 10, 10, 1, 20),
                                  end = c(intmax, intmax, intmax, intmax, intmax, 100, 70)),
        seqlengths = stats::setNames(rep(intmax, 4), paste0("chr", 1:4))
    ))

    grall <- regionStringToGRanges(regions = ".", seqinfo = slens)
    expect_identical(grall, GenomicRanges::GRanges(
        seqnames = c("chr1", "chr2", "chr4"),
        ranges = IRanges::IRanges(start = c(1, 1, 1),
                                  end = c(100, 200, 400)),
        seqlengths = slens)
    )

    expect_warning(grtrim <- regionStringToGRanges(regions = c("chr1:10-101",
                                                               "chr3:1-2",
                                                               "chr4:300-500"),
                                                   seqinfo = slens))
    expect_identical(grtrim, GenomicRanges::GRanges(
        seqnames = c("chr1", "chr3", "chr4"),
        ranges = IRanges::IRanges(start = c(10, 1, 300),
                                  end = c(100, 2, 400)),
        seqlengths = c(chr1 = 100, chr3 = intmax, chr4 = 400))
    )
})

## -------------------------------------------------------------------------- ##
## Checks, refargToDNAStringSet
## -------------------------------------------------------------------------- ##
test_that("refargToDNAStringSet works", {
    # temporarily install BSgenome package
    fagnmfile <- system.file("extdata", "reference.fa.gz",
                             package = "SingleMoleculeGenomicsIO")
    bsgnmfile <- system.file("extdata", "BSgenome.Mmusculus.SingleMoleculeGenomicsIO_0.1.0.tar.gz",
                             package = "SingleMoleculeGenomicsIO")
    rlibdir <- tempfile(pattern = "Rlib")
    dir.create(rlibdir)
    install.packages(bsgnmfile, lib = rlibdir, repos = NULL,
                     quiet = TRUE, verbose = FALSE)
    expect_identical(list.files(rlibdir), "BSgenome.Mmusculus.SingleMoleculeGenomicsIO")
    suppressPackageStartupMessages(suppressWarnings(
        library(BSgenome.Mmusculus.SingleMoleculeGenomicsIO, lib.loc = rlibdir, quietly = TRUE)
    ))
    seqs <- Biostrings::DNAStringSet(as.list(BSgenome.Mmusculus.SingleMoleculeGenomicsIO))

    # expected failures
    expect_error(refargToDNAStringSet(TRUE),
                 "must be either a .BSgenome.")
    expect_error(refargToDNAStringSet("non-existent"),
                 "must be either a .BSgenome.")

    # correct results
    expect_identical(refargToDNAStringSet(seqs), seqs)
    expect_identical(refargToDNAStringSet(fagnmfile), seqs)
    expect_identical(refargToDNAStringSet(BSgenome.Mmusculus.SingleMoleculeGenomicsIO), seqs)

    # clean up
    detach("package:BSgenome.Mmusculus.SingleMoleculeGenomicsIO", unload = TRUE,
           character.only = TRUE)
    unlink(rlibdir)
})

