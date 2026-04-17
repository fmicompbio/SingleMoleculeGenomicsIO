#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <string>
#include <vector>
#include <Rcpp.h>

// for description of the arguments see function definitions in utils.cpp

std::string concatenate_files(std::vector<std::string>, const std::string);
std::string concatenate_hts_files(std::vector<std::string>, const std::string, int);
Rcpp::CharacterVector getChromosomeNamesFromBam(const std::string);
Rcpp::List getTargetsAndTextFromBamHeader(sam_hdr_t *&inbamhdr);
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
int create_multi_region_iterator(std::vector<std::string> &regions,
                                 hts_itr_t *&iter,
                                 hts_idx_t *idx,
                                 sam_hdr_t *in_samhdr,
                                 bool &had_error,
                                 int buffer_len,
                                 char *buffer);
int create_multi_region_iterator_for_sampling(
        int &n_alns_to_sample,
        std::vector<std::string> &tnames_for_sampling,
        double &keep_aln_fraction,
        hts_itr_t *&iter,
        hts_idx_t *idx,
        sam_hdr_t *in_samhdr,
        bool &had_error,
        int buffer_len,
        char *buffer);

int intlist_to_setvector(sam_hdr_t *in_samhdr,
                         Rcpp::List &pos_list,
                         std::vector<std::set<int>> &pos_sets,
                         char *buffer,
                         int &buffer_len,
                         bool &had_error);

typedef struct plpconf {
    const char *inname;
    samFile *infile;
    sam_hdr_t *in_samhdr;
    hts_idx_t *idx;
    hts_itr_t *iter;
} plpconf;

int plpconstructor(void *data, const bam1_t *b, bam_pileup_cd *cd);

int plpdestructor(void *data, const bam1_t *b, bam_pileup_cd *cd);

int readdata(void *data, bam1_t *b);

int check_bam_format(samFile *infile,
                     sam_hdr_t *in_samhdr,
                     bam1_t *bamdata,
                     std::string &bam_format,
                     bool &had_error,
                     char *buffer,
                     int &buffer_len);
