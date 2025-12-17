// #include <cstdlib>
// #include <cstdio>
// #include <cctype>
// #include <vector>
#include <cstring>
#include <string>
// #include <stdbool.h>
#include <set>
#include "htslib/sam.h"
#include <Rcpp.h>
#include <cli/progress.h>
#include "utils.h"

// convert a named Rcpp::List with IntegerVector elements
// to a std::vector<std::set<int>>, where the index in the
// vector corresponds to the target name index as defined in in_samhdr
//
// return 0 if successfull, -1 if a target name was not found in in_samhdr
int intlist_to_setvector(sam_hdr_t *in_samhdr,
                         Rcpp::List &pos_list,
                         std::vector<std::set<int>> &pos_sets,
                         char *buffer,
                         int &buffer_len,
                         bool &had_error) {
    Rcpp::CharacterVector nms;
    Rcpp::IntegerVector vint;
    int i = 0, j = 0, k = 0;

    pos_sets.resize(in_samhdr->n_targets);
    nms = pos_list.names();
    for (i = 0; i < pos_list.size(); i++) {
        j = sam_hdr_name2tid(in_samhdr, ((std::string)nms[i]).c_str());
        if (j >= 0) {
            vint = pos_list[i];
            for (k = 0; k < vint.size(); k++) {
                pos_sets[j].insert((int)vint[k]);
            }
        } else {
            had_error = true;
            snprintf(buffer, buffer_len,
                     "Could not find chromosome %s in bam header\n",
                     ((std::string)nms[i]).c_str());
            return -1;
        }
    }
    return 0;
}

// check if BAM file is conforming to bam_format
// remark: we cannot guarantee in all cases that the bam file is conforming
int check_bam_format(samFile *infile,
                     sam_hdr_t *in_samhdr,
                     bam1_t *bamdata,
                     std::string &bam_format,
                     bool &had_error,
                     char *buffer,
                     int &buffer_len) {
    int ret_r = -1, result = 0;
    while ((ret_r = sam_read1(infile, in_samhdr, bamdata)) >= 0) {
        if (!(bamdata->core.flag & BAM_FUNMAP)) {
            if (bam_format == "Bismark") {
                // XR and XG tags need to exist
                if (bam_aux_get(bamdata, "XR") == NULL || bam_aux_get(bamdata, "XG") == NULL) {
                    had_error = true;
                    snprintf(buffer, buffer_len,
                             "Invalid Bismark bam format (missing XR or XG tags)\n");
                    result = 1;
                }
            } else if (bam_format == "QuasR") {
                // paired alignments need to be on the same strand
                if ((bamdata->core.flag & BAM_FPAIRED) &&
                    (((bamdata->core.flag & BAM_FREVERSE) > 0) != ((bamdata->core.flag & BAM_FMREVERSE) > 0))) {
                    had_error = true;
                    snprintf(buffer, buffer_len,
                             "Invalid QuasR bam format (paired alignments not on same strand)\n");
                    result = 2;
                }
            }
            break;
        }
    }
    return result;
}

// process a single bam record:
// - increase alignment counter (passed by reference)
// - extract information from record (qscore, modification information, etc.)
// - add information to vectors (passed by reference) for later returning to R
int process_mismatch_bam_record(
        bam1_t *bamdata,        // bam record
        std::string bam_format, // format of bam file
        unsigned int &alncnt,   // alignment counter
        bool &had_error,        // error flag
        char *buffer,           // buffer for message
        int &buffer_len,        // allocated length of message buffer
        std::vector<std::set<int>> &pos_plus_sets, // which positions to analyse (plus strand)
        std::vector<std::set<int>> &pos_minus_sets,// which positions to analyse (minus strand)
        int unmod_integer,      // what to count as unmodified
        int mod_integer,        // what to count as modified
        sam_hdr_t *in_samhdr,   // sam file header
        unsigned long long &n_unaligned, // number of unaligned modified bases
        unsigned long long &n_total,     // total number of modified bases
        std::vector<std::string> &variantRefNames, // seqnames of SNV sites
        std::vector<int> &variantRefPositions,     // coordinates of SNV sites
        // vectors for return values (per modification)
        std::vector<std::string> &read_id,
        std::vector<char> &ref_strand,
        std::vector<std::string> &chrom,
        std::vector<int> &ref_position,
        std::vector<double> &mod_prob,
        // vectors for return values (per alignment)
        std::vector<std::string> &df_read_id,
        std::vector<double> &df_qscore,
        std::vector<int> &df_read_length,
        std::vector<int> &df_aligned_length,
        Rcpp::CharacterVector &df_variant_label,
        Rcpp::CharacterVector &df_ref_strand) {

    // allocate variable only used inside process_mismatch_bam_record()
    int i = 0, j = 0, strand = 0, impl = 0, pos = 0, r = 0;
    int useRC = 0, ref_pos = 0, read_pos = 0, op = 0, op_len = 0;
    int this_read_len = bamdata->core.l_qseq;
    size_t size_before_this_read = read_id.size();
    static uint8_t *hitseq = NULL;
    const uint32_t *cigar;
    std::set<int> *pos_set = NULL;

    // process alignment
    alncnt++;

    // determine whether to look at pos_plus_sets (ref_strand = "+", e.g. C-T) or
    //     pos_minus_sets (ref_strand = "-", e.g. A-G) mismatches
    // for bam_format = "QuasR":
    //     pos_plus_sets for plus strand alignments
    //     pos_minus_sets for minus strand alignments
    // for bam_format = "Bismark":
    //     pos_plus_sets for original top strand alignments (XG="CT")
    //     pos_minus_sets for original bottom strand alignments (XG="GA")
    if (bam_format == "QuasR") {
        useRC = BAM_FREVERSE;
    } else if (bam_format == "Bismark") {
        useRC = (strcmp(bam_aux2Z(bam_aux_get(bamdata, "XG")), "GA") == 0) ? !BAM_FUNMAP : BAM_FUNMAP;
    }
    pos_set = useRC ? &(pos_minus_sets[bamdata->core.tid]) : &(pos_plus_sets[bamdata->core.tid]);


    // ### WAS HERE (plan A: copy logic of read_to_reference_pos loop and analyze mismatches?)
    // variables
    cigar = bam_get_cigar(bamdata);  // cigar array
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
            for (j = 0; j < op_len; j++) {
                if (pos_set.count(ref_pos) > 0) {
                    // we need to analyze this position
                    ;
                }
                ref_pos++;
                read_pos++;
            }
            break;

        case BAM_CINS:  // insertion (I)
            if (read_pos + op_len > read_positions[read_positions_index]) {
                // the current read position is within an insertion -->
                //     no corresponding reference position
                while (read_positions_index < read_positions.size() &&
                       read_pos + op_len > read_positions[read_positions_index]) {
                    ref_positions[read_positions_index] = -1;
                    read_positions_index++;
                }
            }
            read_pos += op_len;
            break;

        case BAM_CDEL:       // deletion (D)
        case BAM_CREF_SKIP:  // reference skip (N)
            ref_pos += op_len;
            break;

        case BAM_CSOFT_CLIP:  // soft clipping (S)
            if (read_pos + op_len > read_positions[read_positions_index]) {
                // the current read position is within a soft-clipped region -->
                //     no corresponding reference position
                while (read_positions_index < read_positions.size() &&
                       read_pos + op_len > read_positions[read_positions_index]) {
                    ref_positions[read_positions_index] = -1;
                    read_positions_index++;
                }
            }
            read_pos += op_len;
            break;

        case BAM_CHARD_CLIP:  // hard clipping (H) // # nocov start
        case BAM_CPAD:        // padding (P)
            // these do not consume any positions in the read or reference
            break;

        default:
            Rcpp::warning("Unknown CIGAR operation: %d", op);
        return ref_positions; // # nocov end
        }
    }
/*
    // get aligned sequence of the read
    hitseq = bam_get_seq(bamdata);
    iend = bam_calend(&(bamdata->core), bam1_cigar(bamdata)) - cnt->offset;

    // this part was copied from quantify_methylation (QuasR)
    if ((hit->core.flag & BAM_FPROPER_PAIR) && (hit->core.isize > 0) && (iend > (const uint32_t)(hit->core.mpos) - cnt->offset))
        // left fragment of a paired alignment --> make sure iend does not overlap alignment of right fragment
        iend = (uint32_t)(hit->core.mpos) - cnt->offset;

    if (hit->core.flag & BAM_FREVERSE) {       // alignment on minus strand (reads are reverse complemented, look for G-A mismatches)
        //Rprintf("\nminus strand alignment %d-%d (offset %d), id=%s\n", hit->core.pos+1, bam_calend(&(hit->core), bam1_cigar(hit)), cnt->offset, bam1_qname(hit));
        for(i=(uint32_t)(hit->core.pos)-cnt->offset, j=0; i<iend; i++, j++)
            if(cnt->om[i]) {                   //  target base is 'G'
                //char Twobit2base[] = {'X', 'A', 'C', 'X', 'G', 'X', 'X', 'X', 'T', 'X', 'X', 'X', 'X', 'X', 'X', 'N'};
                //Rprintf("  adding to genomic position %d (read pos %d has %c)\n", i+cnt->offset+1, j+1, Twobit2base[bam1_seqi(hitseq, j)]);
                if(bam1_seqi(hitseq, j)==4) {        //  query base is 'G'
                    cnt->Tm[i]++;
                    cnt->Mm[i]++;
                } else if(bam1_seqi(hitseq, j)==1) { //  query base is 'A'
                    cnt->Tm[i]++;
                }
            }

    } else {                                   // alignment on plus strand (look for C-T mismatches)
        //Rprintf("\nplus strand alignment %d-%d (offset %d), id=%s\n", hit->core.pos+1, bam_calend(&(hit->core), bam1_cigar(hit)), cnt->offset, bam1_qname(hit));
        for(i=(uint32_t)(hit->core.pos)-cnt->offset, j=0; i<iend; i++, j++)
            if(cnt->op[i]) {                    //  target base is 'C'
                //char Twobit2base[] = {'X', 'A', 'C', 'X', 'G', 'X', 'X', 'X', 'T', 'X', 'X', 'X', 'X', 'X', 'X', 'N'};
                //Rprintf("  adding to genomic position %d (read pos %d has %c)\n", i+cnt->offset+1, j+1, Twobit2base[bam1_seqi(hitseq, j)]);
                if(bam1_seqi(hitseq, j)==2) {        //  query base is 'C'
                    cnt->Tp[i]++;
                    cnt->Mp[i]++;
                } else if(bam1_seqi(hitseq, j)==8) { //  query base is 'T'
                    cnt->Tp[i]++;
                }
            }
    }

        // this part is copied from process_bam_record()
        for (i = 0; i < this_read_len; i++) {
            // i is the position in the aligned read (possibly reverse-complemented)
            // pos is the position in the original read (qseq)
            if (bam_is_rev(bamdata)) {
                pos = this_read_len - 1 - i;
            } else{
                pos = i;
            }

            // r: number of found modifications (>=1, 0 or -1 if failed)
            r = bam_mods_at_next_pos(bamdata, ms, mod, sizeof(mod)/sizeof(mod[0]));
            if (r <= -1) {
                had_error = true; // # nocov start
                snprintf(buffer, buffer_len,
                         "Failed to get modifications (read %s)\n",
                         bam_get_qname(bamdata));
                return -2; // # nocov end

            } else if (r > (int)(sizeof(mod) / sizeof(mod[0]))) {
                had_error = true;
                snprintf(buffer, buffer_len,
                         "More modifications than SingleMoleculeGenomicsIO:::read_modbam_cpp can handle (read %s)\n",
                         bam_get_qname(bamdata));
                return -3;

            } else if (!r && impl) {
                // implied base without modification at position i
                if (qseq[pos] == unmodbase) {
                    // base of the right type -> add to results
                    read_id.push_back(bam_get_qname(bamdata));
                    aligned_read_position.push_back(i);
                    forward_read_position.push_back(pos);
                    chrom.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
                    call_code.push_back('-');
                    canonical_base.push_back(canonical);
                    ref_strand.push_back(bam_is_rev(bamdata) ? '-' : '+');
                    mod_prob.push_back(-1.0); // special value of -1.0 indicates inferred unmodified base
                }
            }
            // modifications
            for (j = 0; j < r; j++) {
                if (mod[j].modified_base == modbase) {
                    // found modified base of the right type -> add to results
                    read_id.push_back(bam_get_qname(bamdata));
                    aligned_read_position.push_back(i);
                    forward_read_position.push_back(pos);
                    chrom.push_back(sam_hdr_tid2name(in_samhdr, bamdata->core.tid));
                    call_code.push_back((char) mod[j].modified_base);
                    canonical_base.push_back((char) mod[j].canonical_base);
                    ref_strand.push_back(bam_is_rev(bamdata) == mod[j].strand ? '+' : '-');
                    // `qual` of N corresponds to call probability
                    //     in [N/256, (N+1)/256] -> store midpoint
                    mod_prob.push_back(((double) mod[j].qual + 0.5) / 256.0);
                }
            }
        }
    }

    // ... convert 0-based read positions to 0-based reference coordinates
    //     (a coordinate of -1 means unaligned, e.g. soft-masked)
    std::vector<int> aligned_read_position_converted =
        read_to_reference_pos(bamdata, aligned_read_position);
    ref_position.reserve(ref_position.size() +
        aligned_read_position_converted.size());
    ref_position.insert(ref_position.end(),
                        aligned_read_position_converted.begin(),
                        aligned_read_position_converted.end());
    aligned_read_position.clear();

    // ... remove unaligned (e.g. soft-masked) read-bases
    //     (iterate backwards to avoid messing up indices
    //      when removing elements)
    n_total += ref_position.size();
    for (size_t e = ref_position.size(); e-- > 0;) {
        if (ref_position[e] == -1) {
            n_unaligned++;
            read_id.erase(read_id.begin() + e);
            chrom.erase(chrom.begin() + e);
            forward_read_position.erase(forward_read_position.begin() + e);
            ref_position.erase(ref_position.begin() + e);
            call_code.erase(call_code.begin() + e);
            canonical_base.erase(canonical_base.begin() + e);
            ref_strand.erase(ref_strand.begin() + e);
            mod_prob.erase(mod_prob.begin() + e);
        }
    }

    // ... extract read-level information if the read had modified bases
    if (size_before_this_read < read_id.size()) {
        // ... add to read-level results
        df_read_id.push_back(bam_get_qname(bamdata));
        df_qscore.push_back(extract_qscore(bamdata));
        df_read_length.push_back(this_read_len);
        df_aligned_length.push_back(calculate_aligned_bases(bamdata));
        df_ref_strand.push_back((bamdata->core.flag & BAM_FREVERSE) ? "-" : "+");

        // ... ... variant_label
        if (variantRefNames.size() > 0) {
            df_variant_label.push_back(
                construct_read_label(bamdata, variantRefNames,
                                     variantRefPositions, in_samhdr));
        } else {
            df_variant_label.push_back(NA_STRING);
        }
    }
*/
    return 0;
}

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
//' library(Biostrings)
//' bamfile <- system.file("extdata", "BisSeq_single.bam", package = "SingleMoleculeGenomicsIO")
//' ref <- readDNAStringSet(system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO"))
//' posPlus <- vmatchPattern(pattern = "NCG", subject = ref, max.mismatch = 0,
//'                          with.indels = FALSE, fixed = "subject", algorithm = "auto")
//' posMinus <- vmatchPattern(pattern = "CGN", subject = ref, max.mismatch = 0,
//'                           with.indels = FALSE, fixed = "subject", algorithm = "auto")
//' posPlusList <- lapply(posPlus, function(x) {
//'     start(resize(x = x, width = 1, fix = "center")) - 1L
//' })
//' posMinusList <- lapply(posMinus, function(x) {
//'     start(resize(x = x, width = 1, fix = "center")) - 1L
//' })
//' res1 <- read_mismatchbam_cpp(inname_str = bamfile,
//'                              regions = "chr1:6940000-6955000",
//'                              pos_plus_list = posPlusList,
//'                              pos_minus_list = posMinusList,
//'                              unmod_integer = 2,
//'                              mod_integer = 8,
//'                              level = "summary",
//'                              n_alns_to_sample = 0,
//'                              tnames_for_sampling = "",
//'                              variantRefNames = "",
//'                              variantRefPositions = 0,
//'                              n_threads = 1,
//'                              verbose = TRUE)
//' str(res1)
/*
 bamfile <- system.file("extdata", "BisSeq_single.bam", package = "SingleMoleculeGenomicsIO")
 ref <- readDNAStringSet(system.file("extdata", "reference.fa.gz", package = "SingleMoleculeGenomicsIO"))
 posPlus <- vmatchPattern(pattern = "NCG", subject = ref, max.mismatch = 0,
                         with.indels = FALSE, fixed = "subject", algorithm = "auto")
 posMinus <- vmatchPattern(pattern = "CGN", subject = ref, max.mismatch = 0,
                           with.indels = FALSE, fixed = "subject", algorithm = "auto")
 posPlusList <- lapply(posPlus, function(x) {
     start(resize(x = x, width = 1, fix = "center")) - 1L
 })
 posMinusList <- lapply(posMinus, function(x) {
     start(resize(x = x, width = 1, fix = "center")) - 1L
 })
 res1 <- read_mismatchbam_cpp(inname_str = bamfile,
                               bam_format = "QuasR",
                               regions = "chr1:6940000-6955000",
                               pos_plus_list = posPlusList,
                               pos_minus_list = posMinusList,
                               unmod_integer = 2,
                               mod_integer = 8,
                               level = "summary",
                               n_alns_to_sample = 0,
                               tnames_for_sampling = "",
                               variantRefNames = "",
                               variantRefPositions = 0,
                               n_threads = 1,
                               verbose = TRUE)

 */
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
    int c = 0, i = 0, success = 0;
    unsigned long long n_unaligned = 0, n_total = 0;
    bool had_error = false;
    samFile *infile = NULL;
    sam_hdr_t *in_samhdr = NULL;
    bam1_t *bamdata = NULL;
    hts_idx_t *idx = NULL;
    hts_itr_t *iter = NULL;
    unsigned int regcnt = 0, alncnt = 0;
    char **regions_c = NULL;
    int buffer_len = 2000;
    char buffer[2000];
    const char* inname = inname_str.c_str();
    std::vector<std::set<int>> pos_plus_sets, pos_minus_sets;
    Rcpp::List res;

    // ... return values for mode 1 or 2
    // ... ... one per modification
    std::vector<std::string> read_id;
    std::vector<int> aligned_length;
    std::vector<char> ref_strand;
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

    // check if BAM file is conforming to bam_format
    success = check_bam_format(infile, in_samhdr, bamdata, bam_format,
                               had_error, buffer, buffer_len);
    if (success != 0) {
        goto end;
    }

    // store positions to be analyzed for each chromosome in a vector of sets
    // ... for the plus strand
    success = intlist_to_setvector(in_samhdr, pos_plus_list, pos_plus_sets,
                                   buffer, buffer_len, had_error);
    if (success != 0) {
        goto end;
    }
    // ... and the minus strand
    success = intlist_to_setvector(in_samhdr, pos_minus_list, pos_minus_sets,
                                   buffer, buffer_len, had_error);
    if (success != 0) {
        goto end;
    }

    // start reading according to analysis mode
    if (windowSize > 0) {
    /*
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
    */
    } else {
        if (n_alns_to_sample > 0) {
            /*
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
                    success = process_mismatch_bam_record(bamdata,          // bam record
                                                 alncnt,           // alignment counter
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
                                                 ref_strand,
                                                 chrom,
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
             */

        } else {
            // Mode 1: region-based alignment reading
            // ---------------------------------------------------------------------
            // TODO: lift out parts from below here to reduce redundancy between modes and read_modbam_cpp
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
                    success = process_mismatch_bam_record(
                        bamdata,          // bam record
                        bam_format,       // format of bam file
                        alncnt,           // alignment counter
                        had_error,        // error flag
                        buffer,           // buffer for message
                        buffer_len,       // allocated length of message buffer
                        pos_plus_sets,    // which positions to analyse (plus strand)
                        pos_minus_sets,   // which positions to analyse (minus strand)
                        unmod_integer,    // what to count as unmodified
                        mod_integer,      // what to count as modified
                        in_samhdr,        // sam file header
                        n_unaligned,      // number of unaligned positions
                        n_total,          // total number of positions
                        variantRefNames,  // seqnames of SNV sites
                        variantRefPositions, // coordinates of SNV sites
                        // vectors for return values (per modification)
                        read_id,
                        ref_strand,
                        chrom,
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
                        Rcpp::checkUserInterrupt();
                    } // # nocov end
                    if (success != 0) {
                        goto end;
                    }
                }
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
        snprintf(buffer, buffer_len,
                 "removed %llu unaligned (e.g. soft-masked) of %llu called bases",
                 n_unaligned, n_total);
        cli_alert_info(buffer);
        snprintf(buffer, buffer_len, "read %u alignments", alncnt);
        cli_alert_info(buffer);
    }

    end:
        //cleanup
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
                    Rcpp::_["mod_prob"] = mod_prob,
                    Rcpp::_["read_df"] = df);
            }

            return res;
        }
}
