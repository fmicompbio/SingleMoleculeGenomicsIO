#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <string>
#include <vector>
#include <Rcpp.h>

// for description of the arguments see function definitions in utils.cpp

std::string concatenate_files(std::vector<std::string>, const std::string);
std::string concatenate_hts_files(std::vector<std::string>, const std::string, int);
Rcpp::CharacterVector getChromosomeNamesFromBam(const std::string);
char get_unmodified_base(char);
char complement(char);
int calculate_aligned_bases(bam1_t*);
double extract_qscore(bam1_t*);
int extract_forward_qseq(bam1_t*, char*&, int&);
int extract_mod_probs(bam1_t*,
                      char modbase,
                      char unmodbase,
                      Rcpp::NumericVector* mod_probs,
                      Rcpp::IntegerVector* mod_pos = NULL,  // pass NULL if you don't need positions
                      char* qseq = NULL,
                      hts_base_mod_state* ms = NULL,
                      char* buffer = NULL,
                      int buffer_len = 0);
std::vector<int> read_to_reference_pos(const bam1_t *aln,
                                       const std::vector<int> &read_positions);
std::vector<int> reference_to_read_pos(const bam1_t *aln,
                                       const std::vector<int> &ref_positions);
std::string construct_read_label(const bam1_t *aln,
                                 const std::vector<std::string> &ref_names,
                                 const std::vector<int> &ref_positions,
                                 const sam_hdr_t *hdr);
int open_bam_and_read_index_and_header(bam1_t *&bamdata,
                                       const char *&inname,
                                       samFile *&infile,
                                       hts_idx_t *&idx,
                                       sam_hdr_t *&in_samhdr,
                                       int n_threads,
                                       bool &had_error,
                                       int buffer_len,
                                       char *buffer);
