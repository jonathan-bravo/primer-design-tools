#pragma once
#include "utility.hpp"
#include "bst.hpp"

class RiskOptimizer {
    public:
        RiskOptimizer(std::vector<std::string> &labels, 
                      std::vector<std::string> &data,
                      index_t len,
                      index_t min,
                      index_t max,
                      index_t fumax,
                      index_t rumax,
                      const Args &args);
        RiskOptimizer(std::vector<risk_t> input, index_t len, index_t min, index_t max);
        ~RiskOptimizer();
        std::vector<index_t> random_search(std::size_t limit);
        std::vector<index_t> search(risk_t l, risk_t r, risk_t eps);
        void validate_PDR(std::vector<index_t>);
        risk_t score(std::vector<index_t> PDR, risk_t alpha);
    private:    
        // Core risk component calculation functions
        double calculate_nucleotide_diversity(const std::vector<char>& bases);
        double calculate_gap_penalty(const std::vector<char>& position_data);
        double calculate_ambiguity_penalty(const std::vector<char>& position_data);
        double calculate_shannon_entropy(const std::vector<char>& bases);
        double calculate_conservation_score(const std::vector<char>& bases);
        
        // Window-based risk component calculation functions
        double calculate_gc_content_deviation(const std::vector<std::string>& sequences, 
                                            size_t position, int window_size);
        double calculate_local_complexity(const std::vector<std::string>& sequences, 
                                        size_t position, int window_size);
        double calculate_indel_risk(const std::vector<std::string>& sequences, 
                                size_t position, int window_size);
        
        // Advanced primer-specific risk component calculation functions
        double calculate_tm_stability_risk(const std::vector<char>& valid_bases, 
                                        const std::vector<std::string>& sequences,
                                        size_t position, int window_size);
        double calculate_secondary_structure_risk(const std::vector<char>& valid_bases,
                                                const std::vector<std::string>& sequences,
                                                size_t position, int window_size);
        
        // Main comprehensive risk computation function
        std::vector<risk_t> compute_risk(const std::vector<std::string>& sequences);
        
        // Risk weights structure
        struct RiskWeights {
            double diversity = 1.0;
            double gaps = 1.0;
            double ambiguity = 1.0;
            double entropy = 0.0;
            double conservation = 0.0;
            double gc_deviation = 0.0;
            double complexity = 0.5;
            double indels = 0.0;
            double tm_stability = 0.0;
            double secondary_structure = 0.0;
        };
        
        std::string tmpl;
        std::unordered_map<key_t, risk_t> memo_umap;
        std::unordered_map<key_t, key_t> prev_umap;
        std::vector<std::pair<index_t, risk_t>> *solution_vector;
        std::vector<index_t> valid_results;
        BST **solutions;
        risk_t *memo;
        key_t *prev;
        index_t size, len, min, max, fumax, rumax;
        risk_t *risk, *prefix_sum;
        std::size_t *gc;
        Args args;
        std::vector<index_t> valid_counts_all_pdrs(std::size_t num_threads);
        std::vector<index_t> random_PDR();
        index_t greedy_random_between(index_t cmin, index_t cmax);
        risk_t opt(index_t f, index_t r, index_t pf, index_t pr, risk_t u, risk_t alpha);
        risk_t cost(index_t p, risk_t u, risk_t alpha);
        std::size_t gc_count(index_t p, index_t l);
        bool gc_valid(index_t p);
        index_t valid_count_in_pdr(index_t pdr_start);
        risk_t top_k_opt(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha);
        risk_t top_k_opt_fast(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha);
        risk_t opt_m(index_t f, index_t r, risk_t u, risk_t alpha);
        risk_t top_k_opt_m(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha);
        bool valid(index_t f, index_t r, index_t f_, index_t r_);
        risk_t top_k_opt_mi(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha);

        std::vector<risk_t> row_prefix_;
        std::vector<key_t>  row_prefix_key_;
        index_t             row_prefix_band_ = 0;

        void init_row_prefix(index_t band) {
            if (band == row_prefix_band_) {
                std::fill(row_prefix_.begin(),     row_prefix_.end(),     INF);
                std::fill(row_prefix_key_.begin(), row_prefix_key_.end(), 0);
                return;
            }
            row_prefix_band_ = band;
            row_prefix_.assign((size + 1) * (band + 1), INF);
            row_prefix_key_.assign((size + 1) * (band + 1), 0);
        }

        risk_t& rp(index_t r, index_t fi) {
            return row_prefix_[r * (row_prefix_band_ + 1) + fi];
        }

        key_t& rpk(index_t r, index_t fi) {
            return row_prefix_key_[r * (row_prefix_band_ + 1) + fi];
        }

        risk_t top_k_opt_mi_fast(risk_t u, std::vector<index_t>& min_PDR, risk_t alpha);
};
