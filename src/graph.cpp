#include "graph.hpp"

KPartiteGraph::KPartiteGraph(const std::vector<PrimerOutput>& input) {
    K = input.size() * 2;
    N = 0;
    for (const auto& item : input) {
        N = std::max((long unsigned int)N,
            std::max(item.left_oligos.size(),
                    item.right_oligos.size()));
    }

    valid_N.resize(K);
    for (index_t k = 0; k < K; ++k) {
        index_t amp  = k / 2;
        bool is_left = (k % 2 == 0);
        valid_N[k] = is_left ? input[amp].left_oligos.size()
                             : input[amp].right_oligos.size();
    }

    auto get_oligo = [&](index_t idx) -> const Oligo* {
        index_t k = idx / N;
        index_t n = idx % N;
        if (k % 2 == 0) {
            const auto& v = input[k / 2].left_oligos;
            return n < v.size() ? &v[n] : nullptr;
        } else {
            const auto& v = input[k / 2].right_oligos;
            return n < v.size() ? &v[n] : nullptr;
        }
    };

    // collect all edges to compute
    std::vector<std::pair<index_t, index_t>> edges;
    for (index_t i = 0; i < K * N; i++)
        for (index_t j = i + 1; j < K * N; j++)
            if (part(i) != part(j) && get_oligo(i) && get_oligo(j))
                edges.emplace_back(i, j);

    std::size_t total_edges = edges.size();
    const int TW = 14 + 4 + 10 + 6;
    std::cout << "K = " << K << " N = " << N
              << " nodes = " << K * N << "\n"
              << "edges to compute = " << total_edges << "\n\n"
              << "Constructing KPartite Graph...\n"
              << std::string(TW, '-') << "\n"
              << std::right
              << std::setw(14) << "computed"
              << std::setw(4)  << " / "
              << std::setw(10) << "total"
              << std::setw(6)  << "pct"
              << "\n"
              << std::string(TW, '-') << "\n";

    graph = new weight_t*[K * N];
    for (index_t i = 0; i < K * N; i++) {
        graph[i] = new weight_t[K * N];
        for (index_t j = 0; j < K * N; j++)
            graph[i][j] = (part(i) == part(j)) ? 0 : -1e6;
    }

    vertices = new weight_t[K * N];
    for (index_t i = 0; i < K * N; i++)
        vertices[i] = 0;

    // precompute results in parallel
    std::vector<weight_t> results(edges.size());
    std::atomic<std::size_t> count(0);
    std::size_t step = std::max((std::size_t)1, total_edges / 20);

    unsigned int nthreads = std::thread::hardware_concurrency();
    std::vector<std::thread> threads(nthreads);

    for (unsigned int t = 0; t < nthreads; t++) {
        threads[t] = std::thread([&, t]() {
            for (std::size_t e = t; e < edges.size(); e += nthreads) {
                auto [i, j] = edges[e];
                results[e] = Thal::compute_dimer_dg(get_oligo(i)->seq, get_oligo(j)->seq);
                std::size_t c = count.fetch_add(1) + 1;
                if (c % step == 0) {
                    int pct = (int)(100.0 * c / total_edges);
                    // lock for clean cout
                    static std::mutex mtx;
                    std::lock_guard<std::mutex> lock(mtx);
                    std::cout << std::right
                              << std::setw(14) << c
                              << std::setw(4)  << " / "
                              << std::setw(10) << total_edges
                              << std::setw(5)  << pct << "%\n";
                }
            }
        });
    }
    for (auto& th : threads) th.join();

    // write results to graph (single threaded, no races)
    for (std::size_t e = 0; e < edges.size(); e++) {
        auto [i, j] = edges[e];
        graph[i][j] = graph[j][i] = results[e];
    }

    std::cout << std::string(TW, '-') << "\n"
              << "Computed " << count << " edge weights.\n\n";
}

KPartiteGraph::KPartiteGraph(index_t K, index_t N) {
    this->K = K;
    this->N = N;
    graph = new weight_t*[K * N];
    for (index_t i = 0; i < K * N; i ++) 
        graph[i] = new weight_t[K * N];
    vertices = new weight_t[K * N];
}

KPartiteGraph::KPartiteGraph(std::string filename) {
    std::ifstream fin(filename);
    index_t E;
    fin >> K >> N >> E;
    graph = new weight_t*[K * N];
    for (index_t i = 0; i < K * N; i ++) {
        graph[i] = new weight_t[K * N];
        for (index_t j = 0; j < K * N; j ++) 
            graph[i][j] = 0;
    }
    vertices = new weight_t[K * N];
    for (index_t i = 0; i < K * N; i ++) {
        vertices[i] = 0;
    }
    index_t dummy;
    for (index_t i = 0; i < K * N; i ++) 
        fin >> dummy >> vertices[i];
    for (index_t k = 0; k < E; k ++) {
        index_t i, j;
        weight_t w;
        fin >> i >> j >> w;
        graph[i][j] = graph[j][i] = w;
    }
    fin.close();
}

KPartiteGraph::~KPartiteGraph() {
    for (index_t i = 0; i < K * N; i ++) 
        delete [] graph[i];
    delete [] graph;
    delete [] vertices;
}

void KPartiteGraph::random_weights() {
    for (index_t i = 0; i < K * N; i ++) {
        for (index_t j = i; j < K * N; j ++) {
            if (part(i) == part(j)) 
                graph[i][j] = graph[j][i] = 0;
            else 
                graph[i][j] = graph[j][i] = - (weight_t) rand() / RAND_MAX;
        }
    }
    for (index_t i = 0; i < K * N; i ++) {
        vertices[i] = - (weight_t) rand() / RAND_MAX;
    }
}

void KPartiteGraph::display() {
    for (index_t i = 0; i < K * N; i ++) {
        for (index_t j = 0; j < K * N; j ++) {
            std::cout << std::setw(10) << graph[i][j];
        }
        std::cout << std::endl;
    }
}

weight_t KPartiteGraph::cost(index_t *solution) {
    weight_t total = 0;
    for (index_t i = 0; i < K; i ++) {
        total += vertices[index(i, solution[i])];
        for (index_t j = i + 1; j < K; j ++) {
            index_t i_ = index(i, solution[i]), j_ = index(j, solution[j]);
            total += graph[i_][j_];
        }
    }
    return total / K / (K - 1) * 2;
}

weight_t KPartiteGraph::cost(std::vector<index_t> solution) {
    weight_t total = 0;
    for (index_t i = 0; i < K; i ++) {
        total += vertices[index(i, solution[i])];
        for (index_t j = i + 1; j < K; j ++) {
            index_t i_ = index(i, solution[i]), j_ = index(j, solution[j]);
            total += graph[i_][j_];
        }
    }
    return total / K / (K - 1) * 2;
}

weight_t KPartiteGraph::dfs(index_t k, index_t *solution, std::size_t *count) {
    (*count) ++;
    if (k == K) {
        return cost(solution);
    }
    else {
        weight_t max = -INF;
        for (index_t i = 0; i < N; i ++) {
            solution[k] = i;
            weight_t temp = dfs(k + 1, solution, count);
            if (temp > max) max = temp;
        }
        return max;
    }
}

std::vector<index_t> KPartiteGraph::solve_trivial() {
    std::vector<index_t> solution;
    for (index_t i = 0; i < K; i ++) 
        solution.push_back(0);
    return solution;
}

std::vector<index_t> KPartiteGraph::solve_fast(std::size_t iter) {
    return random_search_single_fast(iter);
}

weight_t KPartiteGraph::solve(std::size_t iter) {
    index_t *solution = new index_t[K];
    std::size_t count = 0;
    return dfs(0, solution, &count);
}

std::vector<index_t> KPartiteGraph::solve_sa(std::size_t iter) {
    return simulated_annealing(iter);
}

std::vector<index_t> KPartiteGraph::solve_tabu(std::size_t restart, std::size_t iter) {
    return tabu_search(restart, iter);
}

std::vector<index_t> KPartiteGraph::solve_ga(std::size_t pop_size, std::size_t generations, double mutation_rate) {
    return genetic_algorithm(pop_size, generations, mutation_rate);
}

weight_t KPartiteGraph::random_search(std::size_t iter) {
    index_t *solution = new index_t[K];
    std::size_t i = 0;
    weight_t max = 0;
    while (i < iter) {
        for (index_t k = 0; k < K; k ++) 
            solution[k] = random_between(0, N);
        weight_t current = cost(solution);
        while (true) {
            index_t best_k = -1, best_n;
            weight_t best = current;
            for (index_t j = 0; j < K; j ++) {
                weight_t next = current;
                for (index_t k = 0; k < K; k ++) {
                    if (k == j) continue;
                    next -= graph[index(k, solution[k])][index(j, solution[j])];
                }
                for (index_t n = 0; n < N; n ++) {
                    weight_t delta = 0;
                    for (index_t k = 0; k < K; k ++) {
                        if (k == j) continue;
                        delta += graph[index(k, solution[k])][index(j, n)];
                    }
                    weight_t temp = next + delta;
                    if (temp > best + 1e-6) {
                        best = temp;
                        best_k = j;
                        best_n = n;
                    }
                }
            }
            if (best_k == -1) break;
            current = best;
            solution[best_k] = best_n;
        }
        i += 1;
        if (current > max) {
            max = current;
        }
    }
    delete [] solution;
    return max;
}

std::vector<index_t> KPartiteGraph::random_search_single_fast(std::size_t restarts) {
    // solution[p] stores the within-part index n (0..N-1)
    std::vector<index_t> solution(K), best_solution(K);

    // S[p] = sum_{q!=p} w( index(p,solution[p]), index(q,solution[q]) )
    std::vector<weight_t> S(K, 0);

    // contrib[p][n] = sum_{q!=p} w( index(p,n), index(q,solution[q]) )
    std::vector<std::vector<weight_t>> contrib(K, std::vector<weight_t>(N, 0));

    auto gv = [&](int p) -> index_t { return index(p, solution[p]); };
    auto w  = [&](index_t u, index_t v) -> weight_t { return graph[u][v]; };
    auto v = [&](index_t u) -> weight_t { return vertices[u]; };

    auto init_S_and_contrib = [&]() {
        // S
        for (int p = 0; p < (int)K; ++p) {
            weight_t s = 0; //v(gv(p)) * 2;
            for (int q = 0; q < (int)K; ++q) if (q != p)
                s += w(gv(p), gv(q));
            S[p] = s;
        }
        // contrib
        for (int p = 0; p < (int)K; ++p) {
            for (int n = 0; n < (int)N; ++n) {
                index_t v_ = index(p, n);
                weight_t s = 0; //v(v_) * 2;
                for (int q = 0; q < (int)K; ++q) if (q != p)
                    s += w(v_, gv(q));
                contrib[p][n] = s;
            }
        }
    };

    auto total_from_S = [&]() -> weight_t {
        // each edge counted twice in Σ_p S[p]
        long double sum2 = 0, sum1 = 0;
        for (int p = 0; p < (int)K; ++p) sum2 += (long double)S[p];
        for (int p = 0; p < (int)K; ++p) sum1 += (long double)v(gv(p));
        return (weight_t)(sum2 * 0.5) + sum1;
    };

    weight_t global_best = std::numeric_limits<weight_t>::lowest();

    std::cout << "Random search: " << restarts << " restarts\n"
              << std::string(51, '-') << "\n"
              << std::right
              << std::setw(10) << "restart"
              << std::setw(15) << "local_opt"
              << std::setw(15) << "global_best"
              << std::setw(8)  << "moves"
              << "\n"
              << std::string(51, '-') << "\n";

    for (std::size_t r = 0; r < restarts; ++r) {
        for (int p = 0; p < (int)K; ++p) solution[p] = random_between(0, N - 1);

        init_S_and_contrib();
        weight_t current = total_from_S();
        int moves = 0;

        while (true) {
            int best_p = -1, best_n = -1;
            weight_t best_gain = 0;

            for (int p = 0; p < (int)K; ++p) {
                for (int n = 0; n < (int)N; ++n) {
                    if (n == solution[p]) continue;
                    weight_t gain = contrib[p][n] - S[p] + v(index(p, n)) - v(gv(p));
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_p = p;
                        best_n = n;
                    }
                }
            }

            if (best_p == -1) break;

            const index_t old_v = gv(best_p);
            solution[best_p] = best_n;
            const index_t new_v = gv(best_p);
            current += best_gain;

            moves++;

            for (int q = 0; q < (int)K; ++q) if (q != best_p)
                S[q] += w(new_v, gv(q)) - w(old_v, gv(q));
            S[best_p] = contrib[best_p][best_n];

            for (int q = 0; q < (int)K; ++q) if (q != best_p)
                for (int m = 0; m < (int)N; ++m) {
                    const index_t u = index(q, m);
                    contrib[q][m] += w(u, new_v) - w(u, old_v);
                }

            for (int n = 0; n < (int)N; ++n) {
                const index_t v_ = index(best_p, n);
                weight_t s = 0;
                for (int q = 0; q < (int)K; ++q) if (q != best_p)
                    s += w(v_, gv(q));
                contrib[best_p][n] = s;
            }
        }

        bool improved = current > global_best;
        if (improved) {
            global_best = current;
            best_solution = solution;
        }

        if (improved || r % 1000 == 0) {
            std::cout << std::right << std::fixed
                      << std::setw(10) << r
                      << std::setw(15) << current / K / (K - 1) * 2
                      << std::setw(15) << global_best / K / (K - 1) * 2
                      << std::setw(8)  << moves
                      << (improved ? "  *" : "")
                      // << cost(best_solution) << " "
                      << "\n";
        }
    }

    std::cout << std::string(51, '-') << "\n"
              << "Best solution: " << global_best << "\n";

    return best_solution;
}

std::vector<index_t> KPartiteGraph::simulated_annealing(std::size_t restarts) {
    const std::size_t proposals_per_step = 200;
    const double      alpha              = 0.985;
    const double      Tmin              = 1e-6;

    std::vector<index_t> solution(K, 0);
    std::vector<index_t> best_solution(K, 0);
    std::vector<index_t> global_solution(K, 0);

    auto objective = [&]() -> weight_t {
        weight_t obj = 0.0;
        for (index_t p = 0; p < K; ++p)
            for (index_t q = p + 1; q < K; ++q)
                obj += graph[index(p, solution[p])][index(q, solution[q])];
        return obj;
    };

    auto gain_of = [&](index_t p, index_t n) -> weight_t {
        weight_t delta = 0.0;
        for (index_t q = 0; q < K; ++q) {
            if (q == p) continue;
            delta += graph[index(p, n          )][index(q, solution[q])]
                   - graph[index(p, solution[p])][index(q, solution[q])];
        }
        return delta;
    };

    const int TW = 8 + 16 + 6 + 16 + 14 + 10;
    std::cout << "\nSimulated Annealing  (K=" << K << ", N=" << N
              << ", restarts=" << restarts << ")\n"
              << std::string(TW, '-') << "\n"
              << std::right
              << std::setw(8)  << "Restart"
              << std::setw(16) << "T0"
              << std::setw(6)  << "Steps"
              << std::setw(16) << "Best(run)"
              << std::setw(14) << "Global best"
              << std::setw(10) << "Improved"
              << "\n"
              << std::string(TW, '-') << "\n";

    weight_t global_best = std::numeric_limits<weight_t>::lowest();

    for (std::size_t r = 0; r < restarts; ++r) {
        for (index_t p = 0; p < K; ++p) solution[p] = random_between(0, N - 1);
        weight_t current   = objective();
        weight_t best_this = current;
        best_solution      = solution;

        auto estimate_T0 = [&]() -> double {
            const std::size_t samples = std::min<std::size_t>(1000, (std::size_t)K * N);
            long double s = 0.0L;
            for (std::size_t t = 0; t < samples; ++t) {
                index_t p = random_between(0, K - 1);
                index_t n = solution[p];
                if (N > 1) do { n = random_between(0, N - 1); } while (n == solution[p]);
                s += std::fabs((long double)gain_of(p, n));
            }
            double avg = samples ? (double)(s / samples) : 1.0;
            return avg > 1e-9 ? avg : 1.0;
        };

        double T  = estimate_T0();
        double T0 = T;
        std::size_t steps_done = 0;

        while (T > Tmin) {
            ++steps_done;
            for (std::size_t it = 0; it < proposals_per_step; ++it) {
                index_t p = random_between(0, K - 1);
                index_t n = solution[p];
                if (N > 1) do { n = random_between(0, N - 1); } while (n == solution[p]);

                weight_t gain = gain_of(p, n);
                bool accept = (gain >= 0.0);
                if (!accept)
                    accept = (std::exp((double)gain / T) >
                              (double)random_between(0, 1000000) / 1000000.0);

                if (accept) {
                    solution[p] = n;
                    current    += gain;
                    if (current > best_this) {
                        best_this     = current;
                        best_solution = solution;
                    }
                }
            }
            T *= alpha;
        }

        bool improved = (best_this > global_best);
        if (improved) {
            global_best = best_this;
            global_solution = best_solution;
        }

        if (improved || (r + 1) % 100 == 0 || r == 0) {
            std::cout << std::right << std::fixed << std::setprecision(4)
                      << std::setw(8)  << (r + 1)
                      << std::setw(16) << T0
                      << std::setw(6)  << steps_done
                      << std::setw(16) << best_this / K / (K - 1) * 2
                      << std::setw(14) << global_best / K / (K - 1) * 2
                      << std::setw(10) << (improved ? "✓" : "")
                      << "\n";
        }
    }

    std::cout << std::string(TW, '-') << "\n"
              << "Global best: " << global_best << "\n";

    return global_solution;
}

std::vector<index_t> KPartiteGraph::tabu_search(std::size_t restarts, std::size_t iterations) {
    const index_t tabu_tenure = K * 2;

    auto objective = [&](const std::vector<index_t>& sol) -> weight_t {
        weight_t obj = 0.0;
        for (index_t p = 0; p < K; ++p)
            for (index_t q = p + 1; q < K; ++q)
                obj += graph[index(p, sol[p])][index(q, sol[q])];
        return obj;
    };

    auto gain_of = [&](const std::vector<index_t>& sol, index_t p, index_t n) -> weight_t {
        weight_t delta = 0.0;
        for (index_t q = 0; q < K; ++q) {
            if (q == p) continue;
            delta += graph[index(p, n      )][index(q, sol[q])]
                   - graph[index(p, sol[p])][index(q, sol[q])];
        }
        return delta;
    };

    const int TW = 8 + 12 + 14 + 14 + 14 + 10;
    std::cout << "\nTabu Search  (K=" << K << ", N=" << N
              << ", restarts=" << restarts
              << ", iterations=" << iterations << ")\n"
              << std::string(TW, '-') << "\n"
              << std::right
              << std::setw(8)  << "Restart"
              << std::setw(12) << "Init"
              << std::setw(14) << "Iters"
              << std::setw(14) << "Best(run)"
              << std::setw(14) << "Global best"
              << std::setw(10) << "Improved"
              << "\n"
              << std::string(TW, '-') << "\n";

    std::vector<index_t> global_best_sol(K, 0);
    weight_t             global_best = std::numeric_limits<weight_t>::lowest();

    for (std::size_t r = 0; r < restarts; ++r) {
        std::vector<index_t> solution(K);
        for (index_t p = 0; p < K; ++p) solution[p] = random_between(0, N - 1);

        weight_t             current      = objective(solution);
        weight_t             best_this    = current;
        std::vector<index_t> best_sol_this = solution;

        // tabu_list[p][n] = iteration until (p,n) is no longer tabu
        std::vector<std::vector<std::size_t>> tabu_list(K, std::vector<std::size_t>(N, 0));

        for (std::size_t it = 0; it < iterations; ++it) {
            weight_t best_gain = std::numeric_limits<weight_t>::lowest();
            index_t  best_p   = 0;
            index_t  best_n   = 0;
            bool     found    = false;

            for (index_t p = 0; p < K; ++p) {
                for (index_t n = 0; n < N; ++n) {
                    if (n == solution[p]) continue;

                    bool     is_tabu    = (tabu_list[p][n] > it);
                    weight_t gain       = gain_of(solution, p, n);
                    bool     aspiration = is_tabu && (current + gain > global_best);

                    if (!is_tabu || aspiration) {
                        if (gain > best_gain) {
                            best_gain = gain;
                            best_p    = p;
                            best_n    = n;
                            found     = true;
                        }
                    }
                }
            }

            if (!found) break;  // all moves tabu and no aspiration

            tabu_list[best_p][solution[best_p]] = it + tabu_tenure;
            solution[best_p] = best_n;
            current += best_gain;

            if (current > best_this) {
                best_this     = current;
                best_sol_this = solution;
            }
        }

        bool improved = (best_this > global_best);
        if (improved) {
            global_best     = best_this;
            global_best_sol = best_sol_this;
        }

        if (improved || (r + 1) % 100 == 0 || r == 0) {
            std::cout << std::right << std::fixed << std::setprecision(4)
                      << std::setw(8)  << (r + 1)
                      << std::setw(12) << current
                      << std::setw(14) << iterations
                      << std::setw(14) << best_this / K / (K - 1) * 2
                      << std::setw(14) << global_best / K / (K - 1) * 2
                      << std::setw(10) << (improved ? "✓" : "")
                      << "\n";
        }
    }

    std::cout << std::string(TW, '-') << "\n"
              << "Global best: " << global_best << "\n";

    return global_best_sol;
}

std::vector<index_t> KPartiteGraph::genetic_algorithm(
    std::size_t pop_size, std::size_t generations, double mutation_rate)
{
    using Chromosome = std::vector<index_t>;

    // must take sol as parameter, not capture
    auto objective = [&](const Chromosome& sol) -> weight_t {
        weight_t obj = 0.0;
        for (index_t p = 0; p < K; ++p)
            for (index_t q = p + 1; q < K; ++q)
                obj += graph[index(p, sol[p])][index(q, sol[q])];
        return obj;
    };

    auto gain_of = [&](const Chromosome& sol, index_t p, index_t n) -> weight_t {
        weight_t delta = 0.0;
        for (index_t q = 0; q < K; ++q) {
            if (q == p) continue;
            delta += graph[index(p, n      )][index(q, sol[q])]
                   - graph[index(p, sol[p])][index(q, sol[q])];
        }
        return delta;
    };

    auto fitness = [&](const Chromosome& c) -> weight_t { return objective(c); };

    auto crossover = [&](const Chromosome& a, const Chromosome& b) -> Chromosome {
        Chromosome child(K);
        for (index_t p = 0; p < K; ++p)
            child[p] = (random_between(0, 1) ? a[p] : b[p]);
        return child;
    };

    auto mutate = [&](Chromosome c) -> Chromosome {
        for (index_t p = 0; p < K; ++p)
            if ((double)random_between(0, 1000) / 1000.0 < mutation_rate)
                c[p] = random_between(0, N - 1);
        return c;
    };

    auto tournament = [&](const std::vector<Chromosome>& pop) -> Chromosome {
        index_t a = random_between(0, (index_t)pop.size() - 1);
        index_t b = random_between(0, (index_t)pop.size() - 1);
        return fitness(pop[a]) > fitness(pop[b]) ? pop[a] : pop[b];
    };

    // init population
    std::vector<Chromosome> pop(pop_size, Chromosome(K));
    for (auto& c : pop)
        for (index_t p = 0; p < K; ++p) c[p] = random_between(0, N - 1);

    Chromosome best     = pop[0];
    weight_t   best_obj = fitness(best);

    const int TW = 10 + 14 + 14 + 10;
    std::cout << "\nGenetic Algorithm  (K=" << K << ", N=" << N
              << ", pop=" << pop_size << ", gen=" << generations << ")\n"
              << std::string(TW, '-') << "\n"
              << std::setw(10) << "Gen"
              << std::setw(14) << "Best(gen)"
              << std::setw(14) << "Global best"
              << std::setw(10) << "Improved"
              << "\n" << std::string(TW, '-') << "\n";

    for (std::size_t g = 0; g < generations; ++g) {
        std::vector<Chromosome> next_pop;
        next_pop.reserve(pop_size);

        // elitism: keep best
        next_pop.push_back(best);

        while (next_pop.size() < pop_size) {
            auto child = crossover(tournament(pop), tournament(pop));
            child      = mutate(child);
            next_pop.push_back(child);
        }
        pop = std::move(next_pop);

        // find best in generation
        weight_t gen_best = std::numeric_limits<weight_t>::lowest();
        for (auto& c : pop) {
            weight_t f = fitness(c);
            if (f > gen_best) gen_best = f;
            if (f > best_obj) { best_obj = f; best = c; }
        }

        bool improved = (gen_best > best_obj);
        if (improved || (g + 1) % 500 == 0 || g == 0) {
            std::cout << std::fixed << std::setprecision(4)
                      << std::setw(10) << (g + 1)
                      << std::setw(14) << gen_best / K / (K - 1) * 2
                      << std::setw(14) << best_obj / K / (K - 1) * 2
                      << std::setw(10) << (improved ? "✓" : "")
                      << "\n";
        }
    }

    std::cout << std::string(TW, '-') << "\n"
              << "Global best: " << best_obj << "\n";
    return best;
}

std::vector<index_t> KPartiteGraph::solve_bottleneck(std::size_t restarts, weight_t& out_threshold) {
    return bottleneck_search(restarts, out_threshold);
}

std::vector<index_t> KPartiteGraph::bottleneck_search(std::size_t restarts, weight_t& out_threshold) {

    // Guard: if any partition is empty the problem is infeasible
    for (index_t p = 0; p < K; ++p) {
        if (valid_N[p] == 0) {
            std::cerr << "bottleneck_search: partition " << p << " has no valid nodes\n";
            out_threshold = std::numeric_limits<weight_t>::lowest();
            return std::vector<index_t>(K, 0);
        }
    }

    auto full_bottleneck = [&](const std::vector<index_t>& sol, index_t& bn_p) -> weight_t {
        weight_t min_w = std::numeric_limits<weight_t>::max();
        for (index_t p = 0; p < K; ++p) {
            // Guard: sol[p] must be within valid range
            assert(sol[p] < valid_N[p]);
            for (index_t q = p + 1; q < K; ++q) {
                assert(sol[q] < valid_N[q]);
                weight_t w = graph[index(p, sol[p])][index(q, sol[q])];
                if (w < min_w) { min_w = w; bn_p = p; }
            }
        }
        return min_w;
    };

    auto evaluate_bn = [&](std::vector<index_t>& sol, index_t bn_p, index_t n) -> weight_t {
        index_t  old_n = sol[bn_p];
        sol[bn_p]      = n;
        index_t  dummy = 0;
        weight_t val   = full_bottleneck(sol, dummy);
        sol[bn_p]      = old_n;
        return val;
    };

    std::vector<index_t> best_solution(K, 0);
    out_threshold = std::numeric_limits<weight_t>::lowest();

    const int TW = 10 + 15 + 15 + 8;
    std::cout << "\nBottleneck Search  (K=" << K << ", N=" << N
              << ", restarts=" << restarts << ")\n"
              << std::string(TW, '-') << "\n"
              << std::right
              << std::setw(10) << "restart"
              << std::setw(15) << "local_opt"
              << std::setw(15) << "global_best"
              << std::setw(8)  << "moves"
              << "\n"
              << std::string(TW, '-') << "\n";

    for (std::size_t r = 0; r < restarts; ++r) {
        std::vector<index_t> solution(K);
        for (index_t p = 0; p < K; ++p)
            solution[p] = random_between(0, valid_N[p] - 1);

        index_t  bn_p    = 0;
        weight_t current = full_bottleneck(solution, bn_p);
        int moves = 0;

        while (true) {
            int      best_n   = -1;
            weight_t best_val = current;

            for (index_t n = 0; n < valid_N[bn_p]; ++n) {
                if (n == solution[bn_p]) continue;
                weight_t candidate = evaluate_bn(solution, bn_p, n);
                if (candidate > best_val) {
                    best_val = candidate;
                    best_n   = (int)n;
                }
            }

            if (best_n == -1) break;

            solution[bn_p] = best_n;
            current = full_bottleneck(solution, bn_p);
            ++moves;
        }

        bool improved = current > out_threshold;
        if (improved) {
            out_threshold = current;
            best_solution = solution;
        }

        if (improved || r % 1000 == 0) {
            std::cout << std::right << std::fixed << std::setprecision(4)
                      << std::setw(10) << r
                      << std::setw(15) << current
                      << std::setw(15) << out_threshold
                      << std::setw(8)  << moves
                      << (improved ? "  *" : "")
                      << "\n";
        }
    }

    std::cout << std::string(TW, '-') << "\n"
              << "Best bottleneck: " << out_threshold << "\n";

    return best_solution;
}