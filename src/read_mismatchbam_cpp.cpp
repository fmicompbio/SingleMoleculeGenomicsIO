// #include <cstdlib>
// #include <cstdio>
// #include <cctype>
// #include <vector>
// #include <cstring>
// #include <string>
// #include <stdbool.h>
#include "htslib/sam.h"
#include <Rcpp.h>
#include <cli/progress.h>
#include "utils.h"

//' Read base modifications from "C-to-T" bam file(s) - C++ helper function
//'
//' Parse C-to-T mismatches and return a list of vectors with
//' information on base states. The function implements four distinct reading
//' modes:
//' \enumerate{
//'     \item{Extraction of read-level modification probabilities for
//'         alignments overlapping provided regions. This mode is selected
//'         if \code{n_alns_to_sample = 0} and \code{level = "read"}.}
//'     \item{Extraction of read-level modification probabilities for alignments
//'         randomly sampled from provided chromosomes. This is selected
//'         if \code{n_alns_to_sample > 0} and \code{level = "read"}.}
//'     \item{Counting of pairs of bases by distance and modification state.
//'         This mode is selected if \code{windowSize > 0}.}
//'     \item{Extraction of summary-level modification counts for alignments
//'         overlapping provided regions. This mode is selected if
//'         \code{n_alns_to_sample = 0} and \code{level = "summary"}.}
//' }
//'
//' @param inname_str Character scalar with name of the input bam file.
//' @param regions Character vector specifying the region(s) for which
//'     to extract overlapping reads, in the form \code{"chr:start-end"}
//' @param pos_plus_list,pos_minus_list Named Rcpp::List of positions on each
//'     chromosome to be evaluated regarding mismatches to reads, seperately
//'     for the plus and the minus strand.
//' @param unmod_integer,mod_integer Integers encoding the read bases to be
//'     interpreted as unmodified or modified, respectively. The encoding
//'     scheme corresponds to the one in bam1_seqi from htslib.
//' @param level Character scalar selecting the level of the returned data
//'     (\code{"read"} or \code{"summary"}).
//' @param n_alns_to_sample Integer defining the number of alignments
//'     to randomly sample.
//' @param tnames_for_sampling String vector with target names (chromosomes)
//'     from which to sample \code{n_alns_to_sample} alignments. Ignored if
//'     \code{n_alns_to_sample = 0}.
//' @param variantRefNames Character vector with target names (chromosomes)
//'     of single nucleotide variants.
//' @param variantRefPositions Integer vector with 0-based target positions
//'     of single nucleotide variants. Expected to have identical length and
//'     to be parallel to \code{variantRefNames}.
//' @param windowSize Numeric scalar giving the maximum window size
//'     covering pairs of modified bases to consider in pair-counting mode.
//'     A window size of 1 corresponds to a single base, a size of 2 to
//'     directly adjacent bases, etc.
//' @param minMapQ Numeric scalar giving the minimal mapping quality to include
//'     alignments in pair-counting mode.
//' @param minAlignedLength Numeric scalar giving the minimal alignment length
//'     to include alignments in pair-counting mode.
//' @param n_threads Integer scalar defining the number of threads to
//'     use for decompressing a sam record. Especially useful in sampling mode
//'     (\code{n_alns_to_sample > 0}), where more time is spend reading and
//'     decompressing bam records than processing them.
//' @param verbose Logical scalar. If \code{TRUE}, report on progress.
//'
//' @return For reading modes 1. and 2., a named list with elements \code{"read_id"},
//'     \code{"ref_position"},
//'     \code{"chrom"}, \code{"ref_strand"}, \code{seq_context},
//'     \code{"mod_prob"} and \code{"read_df"}. The meaning of these elements is
//'     similar to the return value of \code{read_modbam_cpp} and described in
//'     https://nanoporetech.github.io/modkit/intro_extract.html,
//'     apart from \code{"mod_prob"}, which is equal to 0 (1) for bases at
//'     C-to-T mismach positions and equal to 1 (0) for C-to-C match positions
//'     for \code{mismatches_are_unmod = TRUE} (\code{mismatches_are_unmod = FALSE}),
//'     and \code{"read_df"}, which is a \code{data.frame} with one row per
//'     read and columns \code{"read_id"} (the read identifier), \code{"qscore"}
//'     (the read quality score recorded in the \code{qs} tag of each bam record),
//'     \code{"read_length"} (the total read length), and \code{"aligned_length"}
//'     (the number of aligned bases), and \code{"ref_position"}, which is
//'     0-based in the output of \code{modkit extract}, but 1-based here.
//'     For reading mode 3., a named list with elements \code{"read_id"},
//'     \code{"ref_position"}, \code{"chrom"}, \code{"ref_strand"},
//'     \code{seq_context}, \code{"Nvalid"} and \code{"Nmod"}.
//'
//' @examples
//' bamfile <- system.file("extdata", "BisSeq_single.bam", package = "SingleMoleculeGenomicsIO")
//' ref <- Biostrings::readDNAStringSet(system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO"))
//' res1 <- read_mismatchbam_cpp(inname_str = bamfile,
//'                         regions = "chr1:6940000-6955000",
//'                         mismatches_are_unmod = TRUE,
//'                         level = "summary",
//'                         n_alns_to_sample = 0,
//'                         tnames_for_sampling = "",
//'                         variantRefNames = "",
//'                         variantRefPositions = 0,
//'                         n_threads = 1,
//'                         verbose = TRUE)
//' str(res1)
//'
//' @author Charlotte Soneson, Michael Stadler
//'
//' @importFrom cli cli_progress_step cli_progress_done cli_alert_info
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List read_mismatchbam_cpp(std::string inname_str,
                                std::vector<std::string> regions,
                                Rcpp::List pos_plus_list,
                                Rcpp::List pos_minus_list,
                                int unmod_integer,
                                int mod_integer,
                                std::string level,
                                int n_alns_to_sample,
                                std::vector<std::string> tnames_for_sampling,
                                std::vector<std::string> variantRefNames,
                                std::vector<int> variantRefPositions,
                                int windowSize = 0,
                                int minMapQ = 0,
                                int minAlignedLength = 0,
                                int n_threads = 2,
                                bool verbose = false) {
    // turn htslib logging off -> handle via Rcpp::warning or Rcpp::stop
    hts_set_log_level(HTS_LOG_OFF);

    // variable declarations
    // ... R functions
    Rcpp::Environment cli = Rcpp::Environment::namespace_env("cli");
    Rcpp::Function cli_alert_info = cli["cli_alert_info"];

    // ... cli progress bar
    Rcpp::RObject bar;

    // ... general variables
    Rcpp::List res;

    // ... general variables
    int c = 0, i = 0, success = 0;
    unsigned long long n_unaligned = 0, n_total = 0;
    bool had_error = false;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;
    unsigned int regcnt = 0, alncnt = 0;
    char **regions_c = NULL, *qseq = NULL;
    int qseq_len = 0;
    int buffer_len = 2000;
    char buffer[2000];
    const char* inname = inname_str.c_str();

    // ... return values for mode 1 or 2
    // ... ... one per modification
    std::vector<std::string> read_id;
    std::vector<int> aligned_length;
    std::vector<char> ref_strand;
    std::vector<std::string> chrom;
    std::vector<int> aligned_read_position;
    std::vector<int> ref_position;
    std::vector<double> mod_prob;

    // ... ... one per aligned read
    std::vector<std::string> df_read_id;
    std::vector<double> df_qscore;
    std::vector<int> df_read_length;
    std::vector<int> df_aligned_length;
    Rcpp::CharacterVector df_variant_label;
    Rcpp::CharacterVector df_ref_strand;

    // ... return value for mode 3
    Rcpp::NumericMatrix pair_counts;

    // initialize bam data storage
    if (!(bamdata = bam_init1())) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to initialize bamdata\n");
        goto end; // # nocov end
    }

    // open input file
    if (verbose) {
        snprintf(buffer, buffer_len, "opening input file {.file %s} using {%d} thread{?s}", inname, n_threads);
        cli_alert_info(buffer);
    }
    if (!(infile = sam_open(inname, "r"))) {
        had_error = true;
        snprintf(buffer, buffer_len, "Could not open input file %s\n", inname);
        goto end;
    }
    if (n_threads > 1) {
        if (hts_set_threads(infile, n_threads)) {
            had_error = true; // # nocov start
            snprintf(buffer, buffer_len, "Error setting htslib threads to %d\n", n_threads);
            goto end; // # nocov end
        }
    }

    // load index file
    if (!(idx = sam_index_load(infile, inname))) {
        had_error = true;
        snprintf(buffer, buffer_len,
                 "Failed to load the index for %s\n", inname);
        goto end;
    }

    // read header
    if (!(in_samhdr = sam_hdr_read(infile))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len,
                 "Failed to read header from file %s\n", inname);
        goto end; // # nocov end
    }

    // ### WAS HERE
    /*
    // start reading according to analysis mode
    if (windowSize > 0) {
        // Mode 3: counting of pairs of bases by distance and modification state
        // ---------------------------------------------------------------------
        pair_counts = Rcpp::NumericMatrix(windowSize, 4);

        // convert regions to C arrays
        regcnt = (unsigned int) regions.size();
        regions_c = (char**) calloc(regcnt, sizeof(char*));
        for (i = 0; i < (int) regcnt; i++) {
            regions_c[i] = (char*) regions[i].c_str();
        }

        // create multi-region iterator
        if (!(iter = sam_itr_regarray(idx, in_samhdr, regions_c, regcnt))) {
            had_error = true;
            snprintf(buffer, buffer_len, "Failed to get bam iterator\n");
            goto end;
        }

        // iterate over regions
        if (verbose) {
            snprintf(buffer, buffer_len,
                     "counting state-pairs for alignments overlapping {%u} region{?s}",
                     regcnt);
            cli_alert_info(buffer);
            bar = cli_progress_bar(NA_REAL,
                                   Rcpp::List::create(Rcpp::_["clear"] = false,
                                                      Rcpp::_["show_after"] = 0.25));
        }

        // read overlapping alignments using iterator
        while ((c = sam_itr_next(infile, iter, bamdata)) >= 0) {
            if (!(bamdata->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) &&
                (bamdata->core.qual >= minMapQ) &&
                (calculate_aligned_bases(bamdata) >= minAlignedLength)) {

                success = count_pairs_bam_record(
                    bamdata,            // bam record
                    alncnt,             // alignment counter
                    qseq,               // buffer for forward read sequence
                    qseq_len,           // allocated length of qseq
                    ms,                 // modification state struct
                    had_error,          // error flag
                    buffer,             // buffer for message
                    buffer_len,         // allocated length of message buffer
                    modbase,            // modified base to analyze
                    threshUnmod,        // mod_prob <  threshUnmod: unmodified
                    threshMod,          // mod_prob >= thresh_unmod: modified
                    pair_counts);       // count matrix for return value

                if (verbose && CLI_SHOULD_TICK) {
                    cli_progress_set(bar, (double)alncnt); // # nocov
                }
                if (alncnt % 100 == 0) { // # nocov start
                    R_CheckUserInterrupt();
                } // # nocov end
                if (success != 0) {
                    goto end;
                }
            }
        }

    } else {
        if (n_alns_to_sample > 0) {
            // Mode 2: random-sampling-based alignment reading
            // ---------------------------------------------------------------------

            // check if tnames_for_sampling exist and count alignments
            uint64_t mapped = 0, unmapped = 0, total_for_sampling = 0;
            std::set<std::string> tnames_for_sampling_set(tnames_for_sampling.begin(), tnames_for_sampling.end());
            std::set<std::string> tnames_existing;
            double rand_val = 0.0;
            regcnt = 0;
            regions_c = (char**) calloc((unsigned int) tnames_for_sampling.size(),
                         sizeof(char*));
            for (i = 0; i < in_samhdr->n_targets; i++) {
                tnames_existing.insert(in_samhdr->target_name[i]);

                // for each target i that is in tnames_for_sampling_set,
                // get the number of mapped and unmapped records
                // and add it to regions_c
                if (tnames_for_sampling_set.find(in_samhdr->target_name[i]) !=
                    tnames_for_sampling_set.end() &&
                    hts_idx_get_stat(idx, i, &mapped, &unmapped) == 0) {
                    total_for_sampling += mapped;
                    regions_c[regcnt] = in_samhdr->target_name[i];
                    regcnt++;
                }
            }
            for (i = 0; i < (int)tnames_for_sampling.size(); i++) {
                if (tnames_existing.find(tnames_for_sampling[i]) == tnames_existing.end()) {
                    Rcpp::warning("Ignoring unknown target name: %s",
                                  tnames_for_sampling[i].c_str());
                }
            }

            // check if we have enough alignments to sample from
            if (total_for_sampling < (uint64_t)n_alns_to_sample) {
                had_error = true;
                snprintf(buffer, buffer_len,
                         "Cannot sample %d alignments from a total of %" PRIu64 "\n",
                         n_alns_to_sample, total_for_sampling);
                goto end;
            }
            double keep_aln_fraction = (double) n_alns_to_sample / total_for_sampling;
            if (verbose) {
                snprintf(buffer, buffer_len, "sampling alignments with probability %g", keep_aln_fraction);
                cli_alert_info(buffer);
            }

            // create multi-region iterator
            if (!(iter = sam_itr_regarray(idx, in_samhdr, regions_c, regcnt))) {
                had_error = true; // # nocov start
                snprintf(buffer, buffer_len, "Failed to get bam iterator\n");
                goto end; // # nocov end
            }

            // iterate over regions
            if (verbose) {
                snprintf(buffer, buffer_len,
                         "reading alignments overlapping {%u} region{?s}",
                         regcnt);
                cli_alert_info(buffer);
                bar = cli_progress_bar(n_alns_to_sample,
                                       Rcpp::List::create(Rcpp::_["clear"] = false,
                                                          Rcpp::_["show_after"] = 0.25));
            }
            // read overlapping alignments using iterator
            while ((c = sam_itr_next(infile, iter, bamdata)) >= 0) {
                rand_val = R::runif(0, 1);
                if (!(bamdata->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) &&
                    rand_val < keep_aln_fraction) {
                    success = process_bam_record(bamdata,          // bam record
                                                 alncnt,           // alignment counter
                                                 qseq,             // buffer for forward read sequence
                                                 qseq_len,         // allocated length of qseq
                                                 ms,               // modification state struct
                                                 had_error,        // error flag
                                                 buffer,           // buffer for message
                                                 buffer_len,       // allocated length of message buffer
                                                 modbase,          // modified base to analyze
                                                 in_samhdr,        // sam file header
                                                 n_unaligned,      // number of unaligned modified bases
                                                 n_total,          // total number of modified bases
                                                 variantRefNames,  // seqnames of SNV sites
                                                 variantRefPositions, // coordinates of SNV sites
                                                 // vectors for return values (per modification)
                                                 read_id,
                                                 call_code,
                                                 canonical_base,
                                                 ref_strand,
                                                 chrom,
                                                 aligned_read_position,
                                                 ref_position,
                                                 mod_prob,
                                                 // vectors for return values (per alignment)
                                                 df_read_id,
                                                 df_qscore,
                                                 df_read_length,
                                                 df_aligned_length,
                                                 df_variant_label,
                                                 df_ref_strand);
                    if (verbose && CLI_SHOULD_TICK) {
                        cli_progress_set(bar, (double)alncnt);
                    }
                    if (alncnt % 100 == 0) { // # nocov start
                        R_CheckUserInterrupt();
                    } // # nocov end
                    if (success != 0) { // # nocov start
                        goto end;
                    } // # nocov end
                }
            }

        } else {
            // Mode 1: region-based alignment reading
            // ---------------------------------------------------------------------
            // convert regions to C arrays
            regcnt = (unsigned int) regions.size();
            regions_c = (char**) calloc(regcnt, sizeof(char*));
            for (i = 0; i < (int) regcnt; i++) {
                regions_c[i] = (char*) regions[i].c_str();
            }

            // create multi-region iterator
            if (!(iter = sam_itr_regarray(idx, in_samhdr, regions_c, regcnt))) {
                had_error = true;
                snprintf(buffer, buffer_len, "Failed to get bam iterator\n");
                goto end;
            }

            // iterate over regions
            if (verbose) {
                snprintf(buffer, buffer_len,
                         "reading alignments overlapping {%u} region{?s}",
                         regcnt);
                cli_alert_info(buffer);
                bar = cli_progress_bar(NA_REAL,
                                       Rcpp::List::create(Rcpp::_["clear"] = false,
                                                          Rcpp::_["show_after"] = 0.25));
            }
            // read overlapping alignments using iterator
            while ((c = sam_itr_next(infile, iter, bamdata)) >= 0) {
                if (!(bamdata->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {
                    success = process_bam_record(bamdata,          // bam record
                                                 alncnt,           // alignment counter
                                                 qseq,             // buffer for forward read sequence
                                                 qseq_len,         // allocated length of qseq
                                                 ms,               // modification state struct
                                                 had_error,        // error flag
                                                 buffer,           // buffer for message
                                                 buffer_len,       // allocated length of message buffer
                                                 modbase,          // modified base to analyze
                                                 in_samhdr,        // sam file header
                                                 n_unaligned,      // number of unaligned modified bases
                                                 n_total,          // total number of modified bases
                                                 variantRefNames,  // seqnames of SNV sites
                                                 variantRefPositions, // coordinates of SNV sites
                                                 // vectors for return values (per modification)
                                                 read_id,
                                                 call_code,
                                                 canonical_base,
                                                 ref_strand,
                                                 chrom,
                                                 aligned_read_position,
                                                 ref_position,
                                                 mod_prob,
                                                 // vectors for return values (per alignment)
                                                 df_read_id,
                                                 df_qscore,
                                                 df_read_length,
                                                 df_aligned_length,
                                                 df_variant_label,
                                                 df_ref_strand);
                    if (verbose && CLI_SHOULD_TICK) {
                        cli_progress_set(bar, (double)alncnt);
                    }
                    if (alncnt % 100 == 0) { // # nocov start
                        R_CheckUserInterrupt();
                    } // # nocov end
                    if (success != 0) {
                        goto end;
                    }
                }
            }
        }
    }
    */

    if (c != -1) {
        had_error = true;
        snprintf(buffer, buffer_len,
                 "Error while reading from %s - aborting\n", inname);
        if (verbose) { // # nocov start
            cli_progress_done(bar);
        } // # nocov end
        goto end;
    }
    if (verbose) {
        cli_progress_done(bar);
        snprintf(buffer, buffer_len,
                 "removed %llu unaligned (e.g. soft-masked) of %llu called bases",
                 n_unaligned, n_total);
        cli_alert_info(buffer);
        snprintf(buffer, buffer_len, "read %u alignments", alncnt);
        cli_alert_info(buffer);
    }

    end:
        //cleanup
        if (qseq) { // # nocov start
            free((void*) qseq);
            qseq = NULL;
        } // # nocov end
        if (regions_c) {
            free((void*) regions_c);
            regions_c = NULL;
        }
        if (in_samhdr) {
            sam_hdr_destroy(in_samhdr);
        }
        if (infile) {
            sam_close(infile);
        }
        if (bamdata) {
            bam_destroy1(bamdata);
        }
        if (ms) {
            hts_base_mod_state_free(ms);
        }
        if (iter) {
            sam_itr_destroy(iter);
        }
        if (idx) {
            hts_idx_destroy(idx);
        }

        if (had_error) {
            // we encountered an error (message in `buffer`) --> stop
            Rcpp::stop(buffer);

        } else {
            Rcpp::List res;

            if (windowSize > 0) {
                // Mode 3
                // create return list
                res = Rcpp::List::create(
                    Rcpp::_["pair_counts"] = pair_counts
                );

            } else {
                // Mode 1 or 2
                // create data.frame for read-level data
                Rcpp::DataFrame df = Rcpp::DataFrame::create(
                    Rcpp::_["read_id"] = df_read_id,
                    Rcpp::_["qscore"] = df_qscore,
                    Rcpp::_["read_length"] = df_read_length,
                    Rcpp::_["aligned_length"] = df_aligned_length,
                    Rcpp::_["variant_label"] = df_variant_label,
                    Rcpp::_["ref_strand"] = df_ref_strand
                );

                // convert 0-based ref_position to 1-based ref_position
                std::for_each(ref_position.begin(),
                              ref_position.end(),
                              [](int &x) { x += 1; });

                // create return list
                res = Rcpp::List::create(
                    Rcpp::_["read_id"] = read_id,
                    Rcpp::_["ref_position"] = ref_position,
                    Rcpp::_["chrom"] = chrom,
                    Rcpp::_["ref_strand"] = ref_strand,
                    Rcpp::_["call_code"] = call_code,
                    Rcpp::_["mod_prob"] = mod_prob,
                    Rcpp::_["read_df"] = df);
            }

            return res;
        }
}
