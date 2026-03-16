// #include <cstdlib>
// #include <cstdio>
// #include <cctype>
// #include <vector>
#include <cstring>
#include <string>
// #include <stdbool.h>
#include <set>
#include <map>
#include "htslib/sam.h"
#include <Rcpp.h>
#include <cli/progress.h>
#include "utils.h"

#define MISMATCHBAM_MODE_READ      1
#define MISMATCHBAM_MODE_STATE     3

// process a pair of bam records (only for pairs mode):
// - increase alignment counter (passed by reference)
// - extract information from record (qscore, modification information, etc.)
// - add information to vectors (passed by reference) for later returning to R
int process_mismatch_bam_record_pair(
        int mode,               // run mode
        bam1_t *bamdata1,       // bam record (mate 1)
        bam1_t *bamdata2,       // second bam record (mate 2)
        std::string bam_format, // format of bam file
        unsigned int &alncnt,   // alignment counter
        bool &had_error,        // error flag
        char *buffer,           // buffer for message
        int &buffer_len,        // allocated length of message buffer
        std::vector<std::set<int>> &pos_context_sets, // which positions to analyse (plus strand)
        std::vector<std::set<int>> &pos_context_rev_sets,// which positions to analyse (minus strand)
        uint8_t unmod_integer,      // what to count as unmodified
        uint8_t unmod_integer_rev,  // what to count as unmodified, opposite strand
        uint8_t mod_integer,        // what to count as modified
        uint8_t mod_integer_rev,    // what to count as modified, opposite strand
        sam_hdr_t *in_samhdr,   // sam file header
        unsigned long long &n_unaligned, // number of unaligned modified bases
        unsigned long long &n_total,     // total number of modified bases
        // vectors for return values (per modification)
        // ... mode == MISMATCHBAM_MODE_STATE
        Rcpp::NumericMatrix &pair_counts) {

    // allocate variable only used inside process_mismatch_bam_record_pair()
    unsigned int i = 0, j = 0, k = 0, mod_pos = 0;
    int ref_pos = 0, read_pos = 0, op = 0, op_len = 0, qscore_pos = 0;
    int unmod_int = 0, mod_int = 0;
    bool useRC = false, useRC2 = false, found = false;
    uint8_t *hitseq = NULL, fwdbase = 0;
    const uint32_t *cigar;
    std::set<int> *pos_set = NULL;
    // ... mode == MISMATCHBAM_MODE_STATE
    std::vector<int> modposref; // reference position of modified bases
    std::vector<int> qscore;    // quality score of modified bases
    std::vector<int> modstate;  // modification states of the bases (0: unmod, 1: mod)
    bam1_t* bams[] = {bamdata1, bamdata2};  // bam records, for iteration later

    // process alignment pair
    alncnt++;

    // determine whether to look at pos_context_sets (ref_strand = "+", e.g. C-T) or
    //     pos_context_rev_sets (ref_strand = "-", e.g. A-G) mismatches
    // for bam_format = "QuasR":
    //     pos_context_sets for plus strand alignments
    //     pos_context_rev_sets for minus strand alignments
    // for bam_format = "Bismark":
    //     pos_context_sets for original top strand alignments (XG="CT")
    //     pos_context_rev_sets for original bottom strand alignments (XG="GA")
    if (bam_format == "QuasR") {
        useRC = bamdata1->core.flag & BAM_FREVERSE;
        useRC2 = bamdata2->core.flag & BAM_FREVERSE;
        if (useRC != useRC2) {
            had_error = true;
            snprintf(buffer, buffer_len, "Inconsistent strands for mates\n");
            return -1;
        }
    } else if (bam_format == "Bismark") {
        useRC = (strcmp(bam_aux2Z(bam_aux_get(bamdata1, "XG")), "GA") == 0) ? true : false;
        useRC2 = (strcmp(bam_aux2Z(bam_aux_get(bamdata2, "XG")), "GA") == 0) ? true : false;
        if (useRC != useRC2) {
            had_error = true;
            snprintf(buffer, buffer_len, "Inconsistent strands for mates\n");
            return -1;
        }
    }
    pos_set = useRC ? &(pos_context_rev_sets[bamdata1->core.tid]) : &(pos_context_sets[bamdata1->core.tid]);
    unmod_int = useRC ? unmod_integer_rev : unmod_integer;
    mod_int = useRC ? mod_integer_rev : mod_integer;

    for (bam1_t* bamdata : bams) {
        // variables
        cigar = bam_get_cigar(bamdata);  // cigar array
        hitseq = bam_get_seq(bamdata);   // query sequence
        ref_pos = bamdata->core.pos;     // reference position (0-based)
        read_pos = 0;                    // read position (0-based)

        // iterate over the CIGAR operations i
        for (i = 0; i < bamdata->core.n_cigar; i++) {
            op = bam_cigar_op(cigar[i]);         // operation type
            op_len = bam_cigar_oplen(cigar[i]);  // operation length

            switch (op) {
            case BAM_CMATCH:  // match or mismatch (M)
            case BAM_CEQUAL:  // match (=)
            case BAM_CDIFF:   // mismatch (X)
                for (j = 0; j < (unsigned int)op_len; j++) {
                    if (pos_set->count(ref_pos) > 0) {
                        // we need to analyze this position
                        // ... check that the read base is either unmod_integer
                        //     or mod_integer (otherwise do nothing)
                        fwdbase = bam_seqi(hitseq, read_pos);
                        if (fwdbase & (unmod_int | mod_int)) {
                            qscore_pos = bam_get_qual(bamdata)[read_pos];
                            mod_pos = fwdbase == unmod_int ? 0 : 1;
                            // check if position has already been seen, and
                            // keep the observation with the highest qscore
                            found = false;
                            for (k = 0; k < modposref.size(); k++) {
                                if (modposref[k] == ref_pos) {
                                    found = true;
                                    if (qscore[k] < qscore_pos) {
                                        qscore[k] = qscore_pos;
                                        modstate[k] = mod_pos;
                                    }
                                    break;
                                }
                            }
                            if (!found) {
                                modposref.push_back(ref_pos);
                                qscore.push_back(qscore_pos);
                                modstate.push_back(mod_pos);
                            }
                        }
                    }
                    ref_pos++;
                    read_pos++;
                }
                break;

            case BAM_CINS:  // insertion (I)
            case BAM_CSOFT_CLIP:  // soft clipping (S)
                // the current read position is within an insertion or soft-clipped region -->
                //     no corresponding reference position
                read_pos += op_len;
                break;

            case BAM_CDEL:       // deletion (D)
            case BAM_CREF_SKIP:  // reference skip (N)
                ref_pos += op_len;
                break;

            case BAM_CHARD_CLIP:  // hard clipping (H) // # nocov start
            case BAM_CPAD:        // padding (P)
                // these do not consume any positions in the read or reference
                break;

            default:
                Rcpp::warning("Unknown CIGAR operation: %d", op);
            } // # nocov end
        }
    }

    // process modposref and modstate to update counter in pair_counts
    int maxdist = pair_counts.nrow() - 1, currdist = 0;
    for (i = 0; i < modposref.size(); i++) {
        for (j = i; j < modposref.size(); j++) {
            // use abs() to protect against situation where modposref
            // is not sorted
            currdist = std::abs(modposref[j] - modposref[i]);
            if (currdist <= maxdist) {
                pair_counts(currdist, 2 * modstate[i] + modstate[j])++;
            }
        }
    }

    return 0;
}

// process a single bam record:
// - increase alignment counter (passed by reference)
// - extract information from record (qscore, modification information, etc.)
// - add information to vectors (passed by reference) for later returning to R
int process_mismatch_bam_record(
        int mode,               // run mode
        bam1_t *bamdata,        // bam record
        std::string bam_format, // format of bam file
        unsigned int &alncnt,   // alignment counter
        bool &had_error,        // error flag
        char *buffer,           // buffer for message
        int &buffer_len,        // allocated length of message buffer
        std::vector<std::set<int>> &pos_context_sets, // which positions to analyse (plus strand)
        std::vector<std::set<int>> &pos_context_rev_sets,// which positions to analyse (minus strand)
        uint8_t unmod_integer,      // what to count as unmodified
        uint8_t unmod_integer_rev,  // what to count as unmodified, opposite strand
        uint8_t mod_integer,        // what to count as modified
        uint8_t mod_integer_rev,    // what to count as modified, opposite strand
        sam_hdr_t *in_samhdr,   // sam file header
        unsigned long long &n_unaligned, // number of unaligned modified bases
        unsigned long long &n_total,     // total number of modified bases
        std::vector<std::string> &variantRefNames, // seqnames of SNV sites
        std::vector<int> &variantRefPositions,     // coordinates of SNV sites
        // vectors for return values (per modification)
        // ... mode == MISMATCHBAM_MODE_READ
        std::vector<std::string> &read_id,
        std::vector<char> &ref_strand,
        std::vector<double> &qscore,
        std::vector<std::string> &chrom,
        std::vector<int> &ref_position,
        std::vector<double> &mod_prob,
        // ... mode == MISMATCHBAM_MODE_STATE
        Rcpp::NumericMatrix &pair_counts,
        // vectors for return values (per alignment)
        std::vector<std::string> &df_read_id,
        std::vector<double> &df_qscore,
        std::vector<int> &df_read_length,
        std::vector<int> &df_aligned_length,
        Rcpp::CharacterVector &df_variant_label,
        Rcpp::CharacterVector &df_ref_strand) {

    // allocate variable only used inside process_mismatch_bam_record()
    unsigned int i = 0, j = 0;
    int ref_pos = 0, read_pos = 0, op = 0, op_len = 0;
    int unmod_int = 0, mod_int = 0;
    int this_read_len = bamdata->core.l_qseq;
    size_t size_before_this_read = read_id.size();
    bool useRC = false;
    uint8_t *hitseq = NULL, fwdbase = 0;
    const uint32_t *cigar;
    std::set<int> *pos_set = NULL;
    // ... mode == MISMATCHBAM_MODE_STATE
    std::vector<int> modposref; // refernce position of modified bases
    std::vector<int> modstate;  // modification states of the bases (0: unmod, 1: mod)

    // process alignment
    alncnt++;

    // determine whether to look at pos_context_sets (ref_strand = "+", e.g. C-T) or
    //     pos_context_rev_sets (ref_strand = "-", e.g. A-G) mismatches
    // for bam_format = "QuasR":
    //     pos_context_sets for plus strand alignments
    //     pos_context_rev_sets for minus strand alignments
    // for bam_format = "Bismark":
    //     pos_context_sets for original top strand alignments (XG="CT")
    //     pos_context_rev_sets for original bottom strand alignments (XG="GA")
    if (bam_format == "QuasR") {
        useRC = bamdata->core.flag & BAM_FREVERSE;
    } else if (bam_format == "Bismark") {
        useRC = (strcmp(bam_aux2Z(bam_aux_get(bamdata, "XG")), "GA") == 0) ? true : false;
    }
    pos_set = useRC ? &(pos_context_rev_sets[bamdata->core.tid]) : &(pos_context_sets[bamdata->core.tid]);
    unmod_int = useRC ? unmod_integer_rev : unmod_integer;
    mod_int = useRC ? mod_integer_rev : mod_integer;

    // variables
    cigar = bam_get_cigar(bamdata);  // cigar array
    hitseq = bam_get_seq(bamdata);   // query sequence
    ref_pos = bamdata->core.pos;     // reference position (0-based)
    read_pos = 0;                    // read position (0-based)

    // iterate over the CIGAR operations i
    for (i = 0; i < bamdata->core.n_cigar; i++) {
        op = bam_cigar_op(cigar[i]);         // operation type
        op_len = bam_cigar_oplen(cigar[i]);  // operation length

        switch (op) {
        case BAM_CMATCH:  // match or mismatch (M)
        case BAM_CEQUAL:  // match (=)
        case BAM_CDIFF:   // mismatch (X)
            for (j = 0; j < (unsigned int)op_len; j++) {
                if (pos_set->count(ref_pos) > 0) {
                    // we need to analyze this position
                    // ... check that the read base is either unmod_integer
                    //     or mod_integer (otherwise do nothing)
                    fwdbase = bam_seqi(hitseq, read_pos);
                    if (fwdbase & (unmod_int | mod_int)) {

                        if (mode == MISMATCHBAM_MODE_READ) {
                            read_id.push_back(bam_get_qname(bamdata));
                            ref_strand.push_back(useRC ? '-' : '+');
                            qscore.push_back(bam_get_qual(bamdata)[read_pos]);
                            chrom.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
                            ref_position.push_back(ref_pos);
                            mod_prob.push_back(fwdbase == unmod_int ? 0.0 : 1.0);

                        } else if (mode == MISMATCHBAM_MODE_STATE) {
                            modposref.push_back(ref_pos);
                            modstate.push_back(fwdbase == unmod_int ? 0.0 : 1.0);

                        }
                    }
                }
                ref_pos++;
                read_pos++;
            }
            break;

        case BAM_CINS:  // insertion (I)
        case BAM_CSOFT_CLIP:  // soft clipping (S)
            // the current read position is within an insertion or soft-clipped region -->
            //     no corresponding reference position
            read_pos += op_len;
            break;

        case BAM_CDEL:       // deletion (D)
        case BAM_CREF_SKIP:  // reference skip (N)
            ref_pos += op_len;
            break;

        case BAM_CHARD_CLIP:  // hard clipping (H) // # nocov start
        case BAM_CPAD:        // padding (P)
            // these do not consume any positions in the read or reference
            break;

        default:
            Rcpp::warning("Unknown CIGAR operation: %d", op);
        } // # nocov end
    }

    // ... extract read-level information if the read had modified bases
    if (mode == MISMATCHBAM_MODE_READ && size_before_this_read < read_id.size()) {
        // ... add to read-level results
        df_read_id.push_back(bam_get_qname(bamdata));
        df_qscore.push_back(extract_qscore(bamdata));
        df_read_length.push_back(this_read_len);
        df_aligned_length.push_back(calculate_aligned_bases(bamdata));
        df_ref_strand.push_back(useRC ? '-' : '+');

        // ... ... variant_label
        if (variantRefNames.size() > 0) {
            df_variant_label.push_back(
                construct_read_label(bamdata, variantRefNames,
                                     variantRefPositions, in_samhdr));
        } else {
            df_variant_label.push_back(NA_STRING);
        }

    } else if (mode == MISMATCHBAM_MODE_STATE) {
        // process modposref and modstate to update counter in pair_counts
        int maxdist = pair_counts.nrow() - 1, currdist = 0;
        for (i = 0; i < modposref.size(); i++) {
            for (j = i; j < modposref.size(); j++) {
                currdist = std::abs(modposref[j] - modposref[i]);
                if (currdist <= maxdist) {
                    pair_counts(currdist, 2 * modstate[i] + modstate[j])++;
                }
            }
        }
    }

    return 0;
}

//' Read base modifications from mismatch bam file(s) - C++ helper function
//'
//' Parse mismatches and return a list of vectors with
//' information on base states. The function implements four distinct reading
//' modes:
//' \enumerate{
//'     \item{Extraction of read-level modification probabilities for
//'         alignments overlapping provided regions. This mode is selected
//'         if \code{n_alns_to_sample = 0} and \code{windowSize = 0}.}
//'     \item{Extraction of read-level modification probabilities for alignments
//'         randomly sampled from provided chromosomes. This is selected
//'         if \code{n_alns_to_sample > 0} and \code{windowSize = 0}.}
//'     \item{Counting of pairs of bases by distance and modification state.
//'         This mode is selected if \code{n_alns_to_sample = 0} and
//'         \code{windowSize > 0}.}
//'     \item{Counting of pairs of bases by distance and modification state
//'         for (pairs of) alignments randomly sampled from provided chromosomes.
//'         This is selected if \code{n_alns_to_sample > 0} and
//'         \code{windowSize > 0}.}
//' }
//'
//' @param inname_str Character scalar with name of the input bam file.
//' @param regions Character vector specifying the region(s) for which
//'     to extract overlapping reads, in the form \code{"chr:start-end"}
//' @param pos_context_list,pos_context_rev_list Named Rcpp::List of positions
//'     on each chromosome to be evaluated regarding mismatches to reads,
//'     seperately for the plus and the minus strand.
//' @param unmod_integer,mod_integer Integers encoding the read bases to be
//'     interpreted as unmodified or modified, respectively. The encoding
//'     scheme corresponds to the one in bam1_seqi from htslib.
//' @param n_alns_to_sample Integer defining the number of alignments
//'     to randomly sample. Note that for paired-end bam files, individual
//'     reads are sampled and pairs will not be complete.
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
//'     \code{"ref_position"}, \code{"chrom"}, \code{"ref_strand"}, \code{"qscore"},
//'     \code{"mod_prob"} and \code{"read_df"}. The meaning of these elements is
//'     similar to the return value of \code{read_modbam_cpp} and described in
//'     https://nanoporetech.github.io/modkit/intro_extract.html,
//'     apart from \code{"mod_prob"}, which is equal to 0 or 1 for bases at
//'     (mis-)match positions controlled by arguments \code{pos_context_list},
//'     \code{unmod_integer}, \code{mod_integer} and their \code{_rev} variants.
//'     \code{"read_df"} is a \code{data.frame} with one row per read and
//'     columns \code{"read_id"} (the read identifier), \code{"qscore"}
//'     (the read quality score recorded in the \code{qs} tag of each bam record),
//'     \code{"read_length"} (the total read length), and \code{"aligned_length"}
//'     (the number of aligned bases), \code{"variant_label"} and
//'     \code{"ref_strand"}. For reading modes 3. and 4., a named list with a
//'     single element called \code{"pair_counts"}, corresponding to a
//'     \code{windowSize}-by-4 matrix with the numbers of pairs of bases at a
//'     given distance (row) and in a given state (columns: 00, 01, 10 and 11).
//'
//' @examples
//' library(Biostrings)
//' bamfile <- system.file("extdata", "BisSeq_quasr_single.bam", package = "SingleMoleculeGenomicsIO")
//' ref <- readDNAStringSet(system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO"))
//' posContext <- vmatchPattern(pattern = "NCG", subject = ref, max.mismatch = 0,
//'                             with.indels = FALSE, fixed = "subject", algorithm = "auto")
//' posContextRev <- vmatchPattern(pattern = "CGN", subject = ref, max.mismatch = 0,
//'                                with.indels = FALSE, fixed = "subject", algorithm = "auto")
//' posContextList <- lapply(posContext, function(x) {
//'     start(resize(x = x, width = 1, fix = "center")) - 1L
//' })
//' posContextRevList <- lapply(posContextRev, function(x) {
//'     start(resize(x = x, width = 1, fix = "center")) - 1L
//' })
//' res1 <- read_mismatchbam_cpp(inname_str = bamfile, bam_format = "QuasR",
//'                              regions = "chr1:6940000-6955000",
//'                              pos_context_list = posContextList,
//'                              pos_context_rev_list = posContextRevList,
//'                              unmod_integer = 8, unmod_integer_rev = 1,
//'                              mod_integer = 2, mod_integer_rev = 4,
//'                              n_alns_to_sample = 0,
//'                              tnames_for_sampling = character(0),
//'                              variantRefNames = character(0),
//'                              variantRefPositions = integer(0),
//'                              n_threads = 1, verbose = TRUE)
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
                                std::string bam_format,
                                std::vector<std::string> regions,
                                Rcpp::List pos_context_list,
                                Rcpp::List pos_context_rev_list,
                                uint8_t unmod_integer,
                                uint8_t unmod_integer_rev,
                                uint8_t mod_integer,
                                uint8_t mod_integer_rev,
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
    int c = 0, success = 0;
    unsigned long long n_unaligned = 0, n_total = 0;
    bool had_error = false;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL, *bamdata2 = NULL, *dup = NULL;
    std::pair<std::map<std::string,bam1_t*>::iterator, bool> inserted;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;
    unsigned int alncnt = 0;
    int buffer_len = 2000;
    char buffer[2000];
    const char* inname = inname_str.c_str();
    std::vector<std::set<int>> pos_context_sets, pos_context_rev_sets;
    Rcpp::List res;
    std::map<std::string,bam1_t*> curr_records; // records waiting for their mate
    std::map<std::string,bam1_t*>::iterator curr_records_it;
    double keep_aln_fraction = 0.0;

    // ... return values for mode 1 or 2
    // ... ... one per modification
    std::vector<std::string> read_id;
    std::vector<int> aligned_length;
    std::vector<char> ref_strand;
    std::vector<double> qscore;
    std::vector<std::string> chrom;
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

    // prepare bam file for reading
    if (verbose) {
        snprintf(buffer, buffer_len, "opening input file {.file %s} using {%d} thread{?s}", inname, n_threads);
        cli_alert_info(buffer);
    }
    success = open_bam_and_read_index_and_header(bamdata, inname, infile, idx,
                                                 in_samhdr, n_threads,
                                                 had_error, buffer_len, buffer);
    if (success != 0) {
        goto end;
    }

    // check if BAM file is conforming to bam_format
    success = check_bam_format(infile, in_samhdr, bamdata, bam_format,
                               had_error, buffer, buffer_len);
    if (success != 0) {
        goto end;
    }

    // store positions to be analyzed for each chromosome in a vector of sets
    // ... for the plus strand
    success = intlist_to_setvector(in_samhdr, pos_context_list, pos_context_sets,
                                   buffer, buffer_len, had_error);
    if (success != 0) {
        goto end;
    }
    // ... and the minus strand
    success = intlist_to_setvector(in_samhdr, pos_context_rev_list, pos_context_rev_sets,
                                   buffer, buffer_len, had_error);
    if (success != 0) {
        goto end;
    }

    // start reading according to analysis mode
    if (windowSize > 0) {
        // Modes 3 or 4 (state-pair counting)
        pair_counts = Rcpp::NumericMatrix(windowSize, 4);

        if (n_alns_to_sample > 0) {
            // Mode 4: random-sampling-based counting of pairs of bases by distance and modification state
            // -------------------------------------------------------------------------------------------
            success = create_multi_region_iterator_for_sampling(
                n_alns_to_sample, tnames_for_sampling,
                keep_aln_fraction, iter, idx, in_samhdr, had_error,
                buffer_len, buffer);

            if (verbose) {
                snprintf(buffer, buffer_len, "sampling alignments with probability %g", keep_aln_fraction);
                cli_alert_info(buffer);
            }

        } else {
            // Mode 3: region-based counting of pairs of bases by distance and modification state
            // ----------------------------------------------------------------------------------
            success = create_multi_region_iterator(regions, iter, idx, in_samhdr,
                                                   had_error, buffer_len, buffer);
        }
        if (success != 0) {
            goto end;
        }

        // iterate over regions
        if (verbose) {
            snprintf(buffer, buffer_len,
                     "counting state-pairs for alignments");
            cli_alert_info(buffer);
            bar = cli_progress_bar(n_alns_to_sample > 0 ? n_alns_to_sample : NA_REAL,
                                   Rcpp::List::create(Rcpp::_["clear"] = false,
                                                      Rcpp::_["show_after"] = 0.25));
        }

        // read overlapping alignments using iterator
        while ((c = sam_itr_next(infile, iter, bamdata)) >= 0) {
            if (!(bamdata->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) &&
                (bamdata->core.qual >= minMapQ) &&
                (calculate_aligned_bases(bamdata) >= minAlignedLength)) {

                if (!(bamdata->core.flag & BAM_FPAIRED) &&
                    ((n_alns_to_sample == 0) || (R::runif(0, 1) < keep_aln_fraction))) {
                    // single-end - process directly
                    success = process_mismatch_bam_record(
                        MISMATCHBAM_MODE_STATE, // run mode
                        bamdata,          // bam record
                        bam_format,       // format of bam file
                        alncnt,           // alignment counter
                        had_error,        // error flag
                        buffer,           // buffer for message
                        buffer_len,       // allocated length of message buffer
                        pos_context_sets, // which positions to analyse (plus strand)
                        pos_context_rev_sets, // which positions to analyse (minus strand)
                        unmod_integer,    // what to count as unmodified
                        unmod_integer_rev, // what to count as unmodified, opposite strand
                        mod_integer,      // what to count as modified
                        mod_integer_rev,  // what to count as modified, opposite strand
                        in_samhdr,        // sam file header
                        n_unaligned,      // number of unaligned positions
                        n_total,          // total number of positions
                        variantRefNames,  // seqnames of SNV sites
                        variantRefPositions, // coordinates of SNV sites
                        // vectors for return values (per modification)
                        read_id,
                        ref_strand,
                        qscore,
                        chrom,
                        ref_position,
                        mod_prob,
                        pair_counts,
                        // vectors for return values (per alignment)
                        df_read_id,
                        df_qscore,
                        df_read_length,
                        df_aligned_length,
                        df_variant_label,
                        df_ref_strand);
                } else {
                    // paired-end and mapped in proper pair
                    // ... check if the mate has already been seen - if so,
                    //     process the pair; if not, store the current mate
                    std::string curr_read_id(bam_get_qname(bamdata));
                    curr_records_it = curr_records.find(curr_read_id);
                    if (curr_records_it == curr_records.end()) {
                        // mate not yet seen - duplicate and store this record
                        // until the mate is seen
                        dup = bam_dup1(bamdata);
                        if (!dup) { // # nocov start
                            had_error = true;
                            snprintf(buffer, buffer_len, "Failed to duplicate bam record for %s\n", curr_read_id.c_str());
                            goto end;
                        } // # nocov end
                        inserted = curr_records.emplace(curr_read_id, dup);
                        if (!inserted.second) { // # nocov start
                            // if insertion failed - destroy record
                            bam_destroy1(dup);
                        } // #nocov end
                    } else {
                        // mate seen - get it from the map and process the pair
                        bamdata2 = curr_records_it->second;

                        if ((n_alns_to_sample == 0) || (R::runif(0, 1) < keep_aln_fraction)) {
                            success = process_mismatch_bam_record_pair(
                                MISMATCHBAM_MODE_STATE, // run mode
                                bamdata2,         // bam record
                                bamdata,
                                bam_format,       // format of bam file
                                alncnt,           // alignment counter
                                had_error,        // error flag
                                buffer,           // buffer for message
                                buffer_len,       // allocated length of message buffer
                                pos_context_sets, // which positions to analyse (plus strand)
                                pos_context_rev_sets, // which positions to analyse (minus strand)
                                unmod_integer,    // what to count as unmodified
                                unmod_integer_rev, // what to count as unmodified, opposite strand
                                mod_integer,      // what to count as modified
                                mod_integer_rev,  // what to count as modified, opposite strand
                                in_samhdr,        // sam file header
                                n_unaligned,      // number of unaligned positions
                                n_total,          // total number of positions
                                // vectors for return values (per modification)
                                pair_counts);
                        }

                        // remove now processed record from map
                        if (bamdata2) {
                            bam_destroy1(bamdata2);
                            bamdata2 = NULL;
                        }
                        curr_records.erase(curr_records_it);
                    }
                }

                if (verbose && CLI_SHOULD_TICK) { // # nocov start
                    cli_progress_set(bar, (double)alncnt);
                } // # nocov end
                if (alncnt % 100 == 0) { // # nocov start
                    Rcpp::checkUserInterrupt();
                } // # nocov end
                if (success != 0) { // # nocov start
                    goto end;
                } // # nocov end
            }
        }
        // process remaining (unpaired) records in curr_records
        while (!curr_records.empty()) {
            curr_records_it = curr_records.begin();
            bamdata2 = curr_records_it->second;

            if ((n_alns_to_sample == 0) || (R::runif(0, 1) < keep_aln_fraction)) {
                success = process_mismatch_bam_record(
                    MISMATCHBAM_MODE_STATE, // run mode
                    bamdata2,         // bam record
                    bam_format,       // format of bam file
                    alncnt,           // alignment counter
                    had_error,        // error flag
                    buffer,           // buffer for message
                    buffer_len,       // allocated length of message buffer
                    pos_context_sets, // which positions to analyse (plus strand)
                    pos_context_rev_sets, // which positions to analyse (minus strand)
                    unmod_integer,    // what to count as unmodified
                    unmod_integer_rev, // what to count as unmodified, opposite strand
                    mod_integer,      // what to count as modified
                    mod_integer_rev,  // what to count as modified, opposite strand
                    in_samhdr,        // sam file header
                    n_unaligned,      // number of unaligned positions
                    n_total,          // total number of positions
                    variantRefNames,  // seqnames of SNV sites
                    variantRefPositions, // coordinates of SNV sites
                    // vectors for return values (per modification)
                    read_id,
                    ref_strand,
                    qscore,
                    chrom,
                    ref_position,
                    mod_prob,
                    pair_counts,
                    // vectors for return values (per alignment)
                    df_read_id,
                    df_qscore,
                    df_read_length,
                    df_aligned_length,
                    df_variant_label,
                    df_ref_strand);
            }

            // remove now processed record from map
            if (bamdata2) {
                bam_destroy1(bamdata2);
                bamdata2 = NULL;
            }
            curr_records.erase(curr_records_it);

            if (verbose && CLI_SHOULD_TICK) { // # nocov start
                cli_progress_set(bar, (double)alncnt);
            } // # nocov end
            if (alncnt % 100 == 0) { // # nocov start
                Rcpp::checkUserInterrupt();
            } // # nocov end
            if (success != 0) { // # nocov start
                goto end;       // currently there are no failure points in process_mismatch_bam_record
            } // # nocov end
        }
    } else {
        // Modes 1 or 2 (reads)
        if (n_alns_to_sample > 0) {
            // Mode 2: random-sampling-based alignment reading
            // ---------------------------------------------------------------------
            success = create_multi_region_iterator_for_sampling(
                n_alns_to_sample, tnames_for_sampling,
                keep_aln_fraction, iter, idx, in_samhdr, had_error,
                buffer_len, buffer);
            if (verbose) {
                snprintf(buffer, buffer_len, "sampling alignments with probability %g", keep_aln_fraction);
                cli_alert_info(buffer);
            }
        } else {
            // Mode 1: region-based alignment reading
            // ---------------------------------------------------------------------
            success = create_multi_region_iterator(regions, iter, idx, in_samhdr,
                                                   had_error, buffer_len, buffer);
        }

        if (success != 0) {
            goto end;
        }

        // iterate over regions
        if (verbose) {
            snprintf(buffer, buffer_len,
                     "reading alignments");
            cli_alert_info(buffer);
            bar = cli_progress_bar(n_alns_to_sample > 0 ? n_alns_to_sample : NA_REAL,
                                   Rcpp::List::create(Rcpp::_["clear"] = false,
                                                      Rcpp::_["show_after"] = 0.25));
        }

        // read overlapping alignments using iterator
        while ((c = sam_itr_next(infile, iter, bamdata)) >= 0) {
            if (!(bamdata->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) &&
                ((n_alns_to_sample == 0) || (R::runif(0, 1) < keep_aln_fraction))) {
                success = process_mismatch_bam_record(
                    MISMATCHBAM_MODE_READ, // run mode
                    bamdata,          // bam record
                    bam_format,       // format of bam file
                    alncnt,           // alignment counter
                    had_error,        // error flag
                    buffer,           // buffer for message
                    buffer_len,       // allocated length of message buffer
                    pos_context_sets, // which positions to analyse (plus strand)
                    pos_context_rev_sets, // which positions to analyse (minus strand)
                    unmod_integer,    // what to count as unmodified
                    unmod_integer_rev, // what to count as unmodified, opposite strand
                    mod_integer,      // what to count as modified
                    mod_integer_rev,  // what to count as modified, opposite strand
                    in_samhdr,        // sam file header
                    n_unaligned,      // number of unaligned positions
                    n_total,          // total number of positions
                    variantRefNames,  // seqnames of SNV sites
                    variantRefPositions, // coordinates of SNV sites
                    // vectors for return values (per modification)
                    read_id,
                    ref_strand,
                    qscore,
                    chrom,
                    ref_position,
                    mod_prob,
                    pair_counts,
                    // vectors for return values (per alignment)
                    df_read_id,
                    df_qscore,
                    df_read_length,
                    df_aligned_length,
                    df_variant_label,
                    df_ref_strand);
                if (verbose && CLI_SHOULD_TICK) { // # nocov start
                    cli_progress_set(bar, (double)alncnt);
                } // # nocov end
                if (alncnt % 100 == 0) { // # nocov start
                    Rcpp::checkUserInterrupt();
                } // # nocov end
                if (success != 0) { // # nocov start
                    goto end;       // currently there are no failure points in process_mismatch_bam_record
                } // # nocov end
            }
        }
    }

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
        snprintf(buffer, buffer_len, "read %u alignments", alncnt);
        cli_alert_info(buffer);
    }

    end:
        //cleanup
        if (in_samhdr) {
            sam_hdr_destroy(in_samhdr);
        }
        if (infile) {
            sam_close(infile);
        }
        if (bamdata) {
            bam_destroy1(bamdata);
        }
        if (iter) {
            sam_itr_destroy(iter);
        }
        if (idx) {
            hts_idx_destroy(idx);
        }
        for (curr_records_it = curr_records.begin();
             curr_records_it != curr_records.end();
             curr_records_it++) { // # nocov start
            if (curr_records_it->second) {
                bam_destroy1(curr_records_it->second);
            }
        } // # nocov end
        if (curr_records.size() > 0) { // # nocov start
            curr_records.clear();
        } // # nocov end

        if (had_error) {
            // we encountered an error (message in `buffer`) --> stop
            Rcpp::stop(buffer);

        } else {
            Rcpp::List res;

            if (windowSize > 0) {
                // Mode 3 or 4
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
                    Rcpp::_["qscore"] = qscore,
                    Rcpp::_["mod_prob"] = mod_prob,
                    Rcpp::_["read_df"] = df);
            }

            return res;
        }
}
