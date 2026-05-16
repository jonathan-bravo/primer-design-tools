#pragma once

#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <mutex>
#include <thread>
#include <stdexcept>
#include <cmath>
#include <ctime>
#include <cctype>
#include <iomanip>
#include <numeric>
#include <atomic>
#include <set>
#include <map>
#include <future>
#include <filesystem>

namespace fs = std::filesystem;

extern "C" {
#include "thal.h"
#include "libprimer3.h"
#include "p3_seq_lib.h"
#include "oligotm.h"
}

#define key_t    uint64_t
#define index_t  uint32_t
#define risk_t   double
#define weight_t double

#define ALPHA       0.1
#define INF         1e18
#define KEY_LIMIT   (1 << 26)

#define ALPHABET_SIZE 128
#define KELVIN_ZERO   -273.15

#ifndef PRIMER3_PATH
#define PRIMER3_PATH "./primer3"
#endif

struct Oligo {
    std::string seq;
    int start;
    int length;
    double tm;
    double gc;
};

struct PrimerResult {
    Oligo       left, right;
    std::size_t product_size = 0;
    index_t     pdr_left     = 0;
    index_t     pdr_right    = 0;
    std::size_t segment_id   = 0;
    std::size_t pool         = 0;
};

struct PrimerOutput {
    std::vector<PrimerResult> pairs;
    std::vector<Oligo>        left_oligos;
    std::vector<Oligo>        right_oligos;
    index_t     pdr_left     = 0;
    index_t     pdr_right    = 0;
    std::size_t segment_id  = 0;
};

struct Result {
    std::string label;
    std::string data;
    std::string ref_label;
    std::string ref_seq;
    std::size_t ref_index;
    std::size_t flag;
    double      dg;
};

struct Args {
    // std::string mode;

    std::string input_file;
    std::string output_file;
    std::string ref_file;

    std::size_t len_amp     = 420;
    std::size_t len_amp_min = 252;
    std::size_t len_PDR     = 40;
    std::size_t len_ru      = 80;
    std::size_t len_fu      = 80;
    double      u_max       = 10000;
    double      u_min       = 0.1;

    std::size_t iter        = 1000;
    std::size_t seed        = 42;

    std::size_t kmer_len    = 8;
    std::size_t threshold   = 4;
    std::size_t nthreads    = 8;
    std::size_t block_size  = 1 << 27;
    std::size_t chunk_size  = (std::size_t) 1e8;

    std::size_t p3_opt_size = 20;
    std::size_t p3_min_size = 18;
    std::size_t p3_max_size = 25;
    double      p3_opt_tm   = 60.0;
    double      p3_min_tm   = 57.0;
    double      p3_max_tm   = 63.0;
    double      p3_min_gc   = 40.0;
    double      p3_max_gc   = 60.0;
    std::size_t num_return  = 5;

    double      mv          = 50.0;
    double      dv          = 1.5;
    double      dntp        = 0.6;
    double      dna_conc    = 200.0;
    
    double      dg_thres    = -15000;
    double      temp        = 37.0;

    bool        help        = false;
};

struct SegmentContext {
    std::string               name;
    std::vector<std::string>  labels, sequences;
    std::string               tmpl;
    std::vector<index_t>      pdr_regions;
    std::vector<PrimerOutput> candidate_primers, filtered_primers;
    std::vector<Result>       off_target_hits;
};

struct PipelineContext {
    Args                         args;
    std::vector<SegmentContext>  segments;

    // flat members — stages use these directly, unchanged
    std::string               tmpl;
    std::vector<std::string>  labels, sequences;
    std::vector<index_t>      pdr_regions;
    std::vector<PrimerOutput> candidate_primers, filtered_primers;
    std::vector<Result>       off_target_hits;

    // global
    std::vector<index_t>      dimer_solution;
    std::vector<PrimerResult> solution_primers;

    void load_segment(std::size_t i) {
        current_seg_      = i;
        auto& seg         = segments[i];
        tmpl              = seg.tmpl;
        labels            = seg.labels;
        sequences         = seg.sequences;
        pdr_regions       = seg.pdr_regions;
        candidate_primers = seg.candidate_primers;
        filtered_primers  = seg.filtered_primers;
        off_target_hits   = seg.off_target_hits;
    }

    void save_segment() {
        auto& seg             = segments[current_seg_];
        seg.tmpl              = tmpl;
        seg.labels            = labels;
        seg.sequences         = sequences;
        seg.pdr_regions       = pdr_regions;
        seg.candidate_primers = candidate_primers;
        seg.filtered_primers  = filtered_primers;
        seg.off_target_hits   = off_target_hits;
    }

    std::size_t current_seg_id() const { return current_seg_; }

private:
    std::size_t current_seg_ = 0;
};

std::size_t read_fasta(const std::string& file_name,
                       std::vector<std::string>& labels,
                       std::vector<std::string>& data);

void write_fasta(const std::string& filename,
                 const std::vector<std::string>& headers,
                 const std::vector<std::string>& sequences,
                 int line_width);

void char_freq(const std::string& seq);

std::size_t read_rate(std::string filename,
                      std::vector<risk_t>& input);

std::size_t read_rate_2(std::string filename,
                        std::vector<risk_t>& input,
                        std::string& ref);

std::size_t read_rate_ref(std::string filename,
                          std::vector<risk_t>& input,
                          std::string& ref);

void write_results(const std::string& output_file,
                   const std::vector<Result>& results,
                   const std::vector<std::string>& labels);

void test_primer3();

key_t to_key  (index_t f, index_t r, index_t f_, index_t r_);
void  to_index(key_t k, index_t& f, index_t& r, index_t& f_, index_t& r_);

key_t to_key_2  (index_t r, index_t r_);
void  to_index_2(key_t k, index_t& r, index_t& r_);

bool        is_ws  (char c);
char        comp   (char c);
std::string revcomp(const std::string& s);
bool        is_gc  (char c);

std::string display_PDR(std::vector<index_t> PDR);

std::unordered_map<std::string, std::vector<Result>>
group_by_label(const std::vector<Result>& results);

std::vector<risk_t> random_risk   (std::size_t size);
std::size_t         random_between(std::size_t min, std::size_t max);

std::string msa_consensus(const std::vector<std::string>& msa);

void display_primer_output(const PrimerOutput& out);

void convert(const std::vector<PrimerOutput>& primer_output,
             std::vector<std::string>& labels,
             std::vector<std::string>& sequences);


Oligo extract_oligo(const primer_rec& p,
                    const char* sequence,
                    bool is_right);

PrimerOutput extract_all(const p3retval* retval,
                         const seq_args*  sa);

void print_solution(const std::vector<PrimerResult>& solution);

struct RemovalSets {
    std::set<int> remove_pair_F, remove_pair_R, remove_left, remove_right;
    std::size_t num_return;
    bool all_removed(const PrimerOutput& original) const {
        size_t kept_left  = original.left_oligos.size()  - remove_left.size();
        size_t kept_right = original.right_oligos.size() - remove_right.size();
        return kept_left < num_return || kept_right < num_return;
    }
};

RemovalSets computeRemoval(const std::vector<const Result*>& output_results,
                           double threshold, 
                           std::size_t num_return) ;

std::vector<PrimerOutput> filterByDG_relax(
    const std::vector<Result>& results,
    double dg_threshold,
    const std::vector<PrimerOutput>& original_primers,
    std::size_t num_return);
