#pragma once
#include "utility.hpp"
#include "thal.hpp"

class KPartiteGraph {
    public:
        KPartiteGraph(const std::vector<PrimerOutput>& input);
        KPartiteGraph(index_t K, index_t N);
        KPartiteGraph(std::string filename);
        ~KPartiteGraph();
        void random_weights();
        index_t part(index_t i) {return i / N;}
        index_t no(index_t i) {return i % N;}
        index_t index(index_t k, index_t n) {return k * N + n;}
        void display();
        weight_t cost(index_t *solution);
        weight_t cost(std::vector<index_t> solution);
        weight_t solve(std::size_t iter);
        std::vector<index_t> solve_fast(std::size_t iter);
        std::vector<index_t> solve_sa(std::size_t iter);
        std::vector<index_t> solve_tabu(std::size_t restarts, std::size_t iterations);
        std::vector<index_t> solve_ga(std::size_t pop_size, std::size_t generations, double mutation_rate);
        index_t get_N() { return N; }
        std::vector<index_t> solve_trivial();
        std::vector<index_t> solve_bottleneck(std::size_t restarts, weight_t& out_threshold);

    private:
        weight_t *vertices;
        weight_t **graph, **opt_graph, **opt_graph_2;
        index_t K, N;
        std::vector<index_t> valid_N;
        weight_t dfs(index_t k, index_t *solution, std::size_t *count);
        void dfs2(index_t k, index_t *solution, weight_t precost, weight_t *max, std::size_t *count);
        weight_t bound(index_t k, index_t *solution);
        void preprocessing();
        weight_t random_search(std::size_t iter);
        std::vector<index_t> random_search_single_fast(std::size_t restarts);
        std::vector<index_t> beam_search_fast(std::size_t beam_width, std::size_t max_iters);
        std::vector<index_t> simulated_annealing(std::size_t restarts);
        std::vector<index_t> tabu_search(std::size_t restarts, std::size_t iterations);
        std::vector<index_t> genetic_algorithm(std::size_t pop_size, std::size_t generations, double mutation_rate);
        std::vector<index_t> bottleneck_search(std::size_t restarts, weight_t& out_threshold);
};
