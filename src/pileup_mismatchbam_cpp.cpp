/*
 The pileup_mismatchbam_cpp function is in part based on pileup_mod.c (distributed
 with htslib) and subject to the following copyright and permission notice:

    pileup_mod.c --  showcases the htslib api usage

    Copyright (C) 2023 Genome Research Ltd.

    Author: Vasudeva Sarma <vasudeva.sarma@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE

*/

#include <string>
#include <vector>
#include <map>
#include <unistd.h>
#include <ctype.h>
#include <htslib/sam.h>
#include <Rcpp.h>
#include <cli/progress.h>
#include "utils.h"

//' Read and pile-up base modifications from a mismatch bam file.
//'
//' @param inname_str Character scalar with name of the input bam file.
//' @param regions Character vector specifying the region(s) for which
//'     to extract overlapping reads, in the form \code{"chr:start-end"}.
//'     The strings are interpreted by htslib, which understands:
//'     \describe{
//'         \item{"REF" or "REF:"}{: All reads with RNAME REF}
//'         \item{"REF:START"}{: Reads with RNAME REF overlapping START to end of REF}
//'         \item{"REF:-END"}{: Reads with RNAME REF overlapping start of REF to END}
//'         \item{"REF:START-END"}{: Reads with RNAME REF overlapping START to END}
//'         \item{"."}{: All reads from the start of the file}
//'         \item{"*"}{: Unmapped reads at the end of the file (RNAME '*' in SAM)}
//'     }
//' @param pos_context_list,pos_context_rev_list Named Rcpp::List of positions
//'     on each chromosome to be evaluated regarding mismatches to reads,
//'     seperately for the plus and the minus strand.
//' @param unmod_integer,mod_integer Integers encoding the read bases to be
//'     interpreted as unmodified or modified, respectively. The encoding
//'     scheme corresponds to the one in bam1_seqi from htslib.
//' @param level Character scalar selecting the level of the returned data
//'     (\code{"read"} or \code{"summary"}).
//' @param n_threads Integer scalar defining the number of threads to
//'     use for decompressing a sam record. Especially using in sampling mode
//'     (\code{n_alns_to_sample > 0}), where more time is spend reading and
//'     decompressing bam records than processing them.
//' @param verbose Logical scalar. If \code{TRUE}, report on progress.
//'
//' @return A named list with elements \code{"chrom"} (chromosome name),
//'     \code{"ref_position"} (1-based coordinate on \code{"chrom"}),
//'     \code{"ref_strand"} (the strand from which the original molecule
//'     originated). If \code{level} is \code{"summary"},
//'     the list additionally contains slots \code{"Nmod"} (number of modified
//'     bases) and \code{"Nvalid"} (number of total bases). If \code{level} is
//'     \code{"read"}, it contains slots \code{"mod_prob"} and \code{"read_id"}.
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
//' res <- pileup_mismatchbam_cpp(bamfile, "QuasR", "chr1",
//'                               posContextList, posContextRevList,
//'                               unmod_integer = 8, unmod_integer_rev = 1,
//'                               mod_integer = 2, mod_integer_rev = 4,
//'                               level = "summary", 1, TRUE)
//' str(res)
//'
//' @author Michael Stadler, Charlotte Soneson
//'
//' @importFrom cli cli_progress_step cli_progress_done
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List pileup_mismatchbam_cpp(std::string inname_str,
                                  std::string bam_format,
                                  std::vector<std::string> regions,
                                  Rcpp::List pos_context_list,
                                  Rcpp::List pos_context_rev_list,
                                  uint8_t unmod_integer,
                                  uint8_t unmod_integer_rev,
                                  uint8_t mod_integer,
                                  uint8_t mod_integer_rev,
                                  std::string level,
                                  int n_threads = 2,
                                  bool verbose = false) {
    // turn htslib logging off -> handle via Rcpp::warning or Rcpp::stop
    hts_set_log_level(HTS_LOG_OFF);

    // variable declarations
    bam1_t *bamdata = NULL;
    plpconf conf = {0};
    conf.inname = inname_str.c_str();
    bam_plp_t plpiter = NULL;
    int tid = -1, depth = -1, j = 0, success = 0;
    int refpos = -1;
    const bam_pileup1_t *plp = NULL;
    kstring_t insdata = KS_INITIALIZE; // TODO: need kstring_t here and in pileup_modbam_cpp?
    bool had_error = false;
    int buffer_len = 2000;
    char buffer[2000];
    uint64_t refposcount = 0;
    unsigned int regcnt = 0;
    char **regions_c = NULL;
    bool isRC = false;
    int unmod_int = 0, mod_int = 0;
    uint8_t fwdbase = 0;
    std::vector<std::set<int>> pos_context_sets, pos_context_rev_sets;
    Rcpp::List res;

    // ... return values (one per modification)
    std::vector<int> Nmod;
    std::vector<int> Nvalid;
    std::vector<double> mod_prob;
    std::vector<std::string> read_id;
    int curr_Nmod = 0, curr_Nvalid = 0;
    std::map<std::string,uint8_t[2]> curr_reads; // curr_reads[read_id] = {qscore, state(0:unmod, 1:mod)}
    std::map<std::string,uint8_t[2]>::iterator curr_reads_it;
    std::vector<std::string> chrom;
    std::vector<int> ref_position;
    std::vector<char> ref_strand;

    // ... return values (one per aligned read)
    std::vector<std::string> df_read_id;
    std::vector<double> df_qscore;
    std::vector<int> df_read_length;
    std::vector<int> df_aligned_length;
    Rcpp::CharacterVector df_variant_label;
    Rcpp::CharacterVector df_ref_strand;

    // ... cli progress bar
    Rcpp::RObject bar;

    // prepare bam file for reading
    success = open_bam_and_read_index_and_header(bamdata, conf.inname,
                                                 conf.infile, conf.idx,
                                                 conf.in_samhdr, n_threads,
                                                 had_error, buffer_len, buffer);
    if (success != 0) {
        goto end;
    }

    success = create_multi_region_iterator(regions, regcnt, regions_c,
                                           conf.iter, conf.idx, conf.in_samhdr,
                                           had_error, buffer_len, buffer);
    if (success != 0) {
        goto end;
    }

    // initialize pileup iterator
    if (!(plpiter = bam_plp_init(readdata, &conf))) {
        had_error = true; // # nocov start
        snprintf(buffer, buffer_len, "Failed to initialize pileup data\n");
        goto end; // # nocov end
    }

    // set constructor and destructor callbacks
    bam_plp_constructor(plpiter, plpconstructor);
    bam_plp_destructor(plpiter, plpdestructor);

    if (verbose) {
        bar = cli_progress_bar(
            NA_REAL,
            Rcpp::List::create(
                Rcpp::_["clear"] = false,
                Rcpp::_["show_after"] = 0.25,
                Rcpp::_["format"] = "{cli::pb_spin} {sprintf(\"%.3f\", cli::pb_current / 1e6)} Mio. genomic positions processed ({sprintf(\"%.3f Mio./s\", cli::pb_rate_raw / 1e6)}) [{cli::pb_elapsed}]"));
    }

    // store positions to be analyzed for each chromosome in a vector of sets
    // ... for the plus strand
    success = intlist_to_setvector(conf.in_samhdr, pos_context_list, pos_context_sets,
                                   buffer, buffer_len, had_error);
    if (success != 0) {
        goto end;
    }
    // ... and the minus strand
    success = intlist_to_setvector(conf.in_samhdr, pos_context_rev_list, pos_context_rev_sets,
                                   buffer, buffer_len, had_error);
    if (success != 0) {
        goto end;
    }

    while ((plp = bam_plp_auto(plpiter, &tid, &refpos, &depth))) {
        // should refpos be analyzed?
        if (pos_context_sets[tid].count(refpos) > 0) {
            isRC = false;
            unmod_int = unmod_integer;
            mod_int = mod_integer;

        } else if (pos_context_rev_sets[tid].count(refpos) > 0) {
            isRC = true;
            unmod_int = unmod_integer_rev;
            mod_int = mod_integer_rev;

        } else {
            continue;
        }

        curr_Nmod = 0;
        curr_Nvalid = 0;
        curr_reads.clear();

        // iterate over reads overlapping refpos
        for (j = 0; j < depth; ++j) {
            // is read j on the right strand?
            if ((bam_format == "QuasR" && (isRC != ((plp[j].b->core.flag & BAM_FREVERSE) != 0))) ||
                (bam_format == "Bismark" && (isRC != ((strcmp(bam_aux2Z(bam_aux_get(plp[j].b, "XG")), "GA") == 0) ? true : false)))) {
                continue;
            }

            // if this is the first time the read is seen, add it to the
            // read df vectors
            if (level == "read" && plp[j].is_head &&
                !(plp[j].b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))) {
                df_read_id.push_back(bam_get_qname(plp[j].b));
                df_qscore.push_back(extract_qscore(plp[j].b));
                df_read_length.push_back(plp[j].b->core.l_qseq);
                df_aligned_length.push_back(calculate_aligned_bases(plp[j].b));
                df_variant_label.push_back(NA_STRING);
                df_ref_strand.push_back(isRC ? "-" : "+");
            }

            if (plp[j].is_del || plp[j].is_refskip ||
                (plp[j].b->core.flag & BAM_FSECONDARY) ||
                (plp[j].b->core.flag & BAM_FSUPPLEMENTARY)) {
                continue;
            }

            // ... check that the read base is either unmod_integer
            //     or mod_integer (otherwise do nothing)
            fwdbase = bam_seqi(bam_get_seq(plp[j].b), plp[j].qpos);
            if (!((fwdbase & (unmod_int | mod_int)) != 0)) {
                continue;
            }

            // store read in curr_reads
            std::string curr_read_id(bam_get_qname(plp[j].b));
            curr_reads_it = curr_reads.find(curr_read_id);
            if (curr_reads_it == curr_reads.end()) {
                // add new read
                curr_reads[curr_read_id][0] = (uint8_t)bam_get_qual(plp[j].b)[j];
                curr_reads[curr_read_id][1] = (uint8_t)((fwdbase & mod_int) != 0 ? 1 : 0);
            } else {
                // compare to current record and keep the one with highest qscore
                if ((int)bam_get_qual(plp[j].b)[j] > curr_reads_it->second[0] &&
                    (int)(fwdbase & mod_int ? 1 : 0) != curr_reads_it->second[1]) {
                    curr_reads_it->second[0] = (uint8_t)bam_get_qual(plp[j].b)[j];
                    curr_reads_it->second[1] = (uint8_t)((fwdbase & mod_int) != 0 ? 1 : 0);
                }
            }
        }

        // save data from reads stored for refpos
        if (level == "summary") {
            for (curr_reads_it = curr_reads.begin();
                 curr_reads_it != curr_reads.end();
                 curr_reads_it++) {
                curr_Nvalid++;
                curr_Nmod += curr_reads_it->second[1];
            }
            if (curr_Nvalid > 0) {
                chrom.push_back(sam_hdr_tid2name(conf.in_samhdr, tid));
                ref_position.push_back(refpos + 1);
                ref_strand.push_back(isRC ? '-' : '+');
                Nvalid.push_back(curr_Nvalid);
                Nmod.push_back(curr_Nmod);
            }
        } else if (level == "read") {
            for (curr_reads_it = curr_reads.begin();
                 curr_reads_it != curr_reads.end();
                 curr_reads_it++) {
                chrom.push_back(sam_hdr_tid2name(conf.in_samhdr, tid));
                ref_position.push_back(refpos + 1);
                ref_strand.push_back(isRC ? '-' : '+');
                read_id.push_back(curr_reads_it->first);
                mod_prob.push_back((double) curr_reads_it->second[1]);
            }
        }

        refposcount++;
        if (verbose && CLI_SHOULD_TICK) {
            // # nocov start
            cli_progress_set(bar, (double)refposcount);
            // # nocov end
        }
        if (refposcount % 1000000 == 0) { // # nocov start
            Rcpp::checkUserInterrupt();
        } // # nocov end
    }

end:
    //clean up
    if (conf.in_samhdr) {
        sam_hdr_destroy(conf.in_samhdr);
    }
    if (conf.infile) {
        sam_close(conf.infile);
    }
    if (conf.iter) {
        sam_itr_destroy(conf.iter);
    }
    if (conf.idx) {
        hts_idx_destroy(conf.idx);
    }
    if (bamdata) {
        bam_destroy1(bamdata);
    }
    if (plpiter) {
        bam_plp_destroy(plpiter);
    }
    ks_free(&insdata);

    if (had_error) {
        // we encountered an error (message in `buffer`) --> stop
        // # nocov start
        Rcpp::stop(buffer);
        // # nocov end

    } else {
        // create return list
        if (level == "summary") {
            res = Rcpp::List::create(
                Rcpp::_["chrom"] = chrom,
                Rcpp::_["ref_position"] = ref_position,
                Rcpp::_["ref_strand"] = ref_strand,
                Rcpp::_["Nmod"] = Nmod,
                Rcpp::_["Nvalid"] = Nvalid);
        } else if (level == "read") {
            Rcpp::DataFrame df = Rcpp::DataFrame::create(
                Rcpp::_["read_id"] = df_read_id,
                Rcpp::_["qscore"] = df_qscore,
                Rcpp::_["read_length"] = df_read_length,
                Rcpp::_["aligned_length"] = df_aligned_length,
                Rcpp::_["variant_label"] = df_variant_label,
                Rcpp::_["ref_strand"] = df_ref_strand
            );
            res = Rcpp::List::create(
                Rcpp::_["chrom"] = chrom,
                Rcpp::_["ref_position"] = ref_position,
                Rcpp::_["ref_strand"] = ref_strand,
                Rcpp::_["mod_prob"] = mod_prob,
                Rcpp::_["read_id"] = read_id,
                Rcpp::_["read_df"] = df);
        }

        return res;
    }
}

