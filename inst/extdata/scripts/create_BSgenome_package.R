# Build BSgenome package for reference genome

library(BSgenomeForge)
library(Biostrings)
library(Seqinfo)
library(rtracklayer)
library(pkgbuild)

# Define variables
genomeFasta <- system.file("extdata", "reference.fa.gz",
                           package = "SingleMoleculeGenomicsIO")
BSgenomeName <- "BSgenome.Mmusculus.SingleMoleculeGenomicsIO"
BSgenomeDir <- tempdir()
genome <- "GRCm38"
release <- "December 2025"
organism <- "Mus musculus"
provider <- "SingleMoleculeGneomicsIO"

# Temporary directory for .2bit file
seqdir <- tempfile(pattern = "BSgenome.seqs", tmpdir = tempdir())
dir.create(seqdir)
seqdir

# Write sequences to .2bit file
chrs <- Biostrings::readDNAStringSet(genomeFasta)
chrs
out_2bit <- file.path(seqdir, "out.2bit")
rtracklayer::export.2bit(chrs, out_2bit)

dir(seqdir)

# Create package seed
pkgseed <- new(
    Class = "BSgenomeDataPkgSeed",
    Package = BSgenomeName,
    Title = paste0("Example sequence from ", genome,
                   ", for use with SingleMoleculeGenomicsIO, ", release),
    Description = paste0("Example sequence from ", genome,
                         ", for use with SingleMoleculeGenomicsIO, ", release),
    Version = "0.1.0",
    Author = "FMI CompBio",
    Maintainer = "FMI CompBio <bioinformatics@fmi.ch>",
    License = "GPL-2",
    organism = organism,
    common_name = "Mouse",
    provider = provider,
    genome = genome,
    release_date = as.character(Sys.Date()),
    organism_biocview = organism,
    BSgenomeObjname = BSgenomeName,
    seqfile_name = "out.2bit",
    # the following fails with an evaluation error (object "chrM" not found)
    #   and deparsing does not work as `chrs` will not be available
    #   --> hard code "chrM", wrapped in quotes for eval (used in both Human and Mouse)
    # circ_seqs = intersect(names(chrs), GenomeInfoDb:::DEFAULT_CIRC_SEQS)
    circ_seqs = "'chrM'"
)

# Build the package
BSgenomeForge::forgeBSgenomeDataPkg(
    x = pkgseed,
    seqs_srcdir = seqdir,
    destdir = BSgenomeDir
)

# Build the source tarball
pkgfname <- pkgbuild::build(path = file.path(BSgenomeDir, BSgenomeName),
                            dest_path = "./inst/extdata", binary = FALSE,
                            vignettes = FALSE, manual = FALSE)
pkgfname

# Session info
date()
sessionInfo()
