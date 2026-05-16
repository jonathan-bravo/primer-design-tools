#include "risk_optimizer.hpp"

RiskOptimizer::RiskOptimizer(std::vector<std::string> &labels, 
                             std::vector<std::string> &data,
                             index_t len,
                             index_t min,
                             index_t max,
                             index_t fumax,
                             index_t rumax,
                             const Args &args) {
    std::vector<risk_t> input = compute_risk(data);
    this->tmpl = msa_consensus(data);
    size = input.size();
    risk = new risk_t[size];
    prefix_sum = new risk_t[size + 1];
    prefix_sum[0] = 0;
    gc = new std::size_t[size + 1];
    gc[0] = 0;
    for (index_t i = 0; i < input.size(); i ++) {
        risk[i] = input[i];
        prefix_sum[i + 1] = prefix_sum[i] + risk[i];
        gc[i + 1] = gc[i] + (is_gc(tmpl[i]) ? 1 : 0);
    }
    this->len = len;
    this->min = min;
    this->max = max;
    this->fumax = fumax;
    this->rumax = rumax;
    memo = new risk_t[KEY_LIMIT];
    prev = new key_t[KEY_LIMIT];
    solution_vector = new std::vector<std::pair<index_t, risk_t>>[size];
    solutions = new BST*[size];
    this->args = args;
    valid_results = valid_counts_all_pdrs(4);
}

RiskOptimizer::RiskOptimizer(std::vector<risk_t> input, 
                             index_t len, 
                             index_t min, 
                             index_t max) {
    size = input.size();
    risk = new risk_t[size];
    prefix_sum = new risk_t[size + 1];
    prefix_sum[0] = 0;
    for (index_t i = 0; i < input.size(); i ++) {
        risk[i] = input[i];
        prefix_sum[i + 1] = prefix_sum[i] + risk[i];
    }
    this->len = len;
    this->min = min;
    this->max = max;
    memo = new risk_t[KEY_LIMIT];
    prev = new key_t[KEY_LIMIT];
    solution_vector = new std::vector<std::pair<index_t, risk_t>>[size];
    solutions = new BST*[size];
}

risk_t RiskOptimizer::cost(index_t p, risk_t u, risk_t alpha) {
    risk_t sum = prefix_sum[p + len] - prefix_sum[p];
    risk_t sum2 = sum * sum - u;
    if (sum2 < 0) return u * alpha;
    return sum2 + u * alpha;
}

std::size_t RiskOptimizer::gc_count(index_t p, index_t l) {
    std::size_t count = gc[p + l] - gc[p];
    return count;
}

index_t RiskOptimizer::valid_count_in_pdr(index_t pdr_start) {
    index_t valid_count = 0;
    for (std::size_t pl = args.p3_min_size; pl <= args.p3_max_size; ++pl) {
        for (index_t p = pdr_start; p + pl < pdr_start + len; ++p) {

            double r = (double) gc_count(p, pl) / pl * 100;
            if (r < args.p3_min_gc || r > args.p3_max_gc) continue;

            std::string primer_seq = tmpl.substr(p, pl);

            tm_ret tm_result = oligotm(
                primer_seq.c_str(),
                args.dna_conc,
                args.mv,
                args.dv,
                args.dntp,
                0.0,
                0.6,
                0.0,
                santalucia_auto,
                santalucia,
                -10.0
            );

            if (tm_result.Tm < args.p3_min_tm || tm_result.Tm > args.p3_max_tm) continue;

            ++ valid_count;
        }
    }
    return valid_count;
}

bool RiskOptimizer::gc_valid(index_t p) {
    index_t count = valid_results[p];
    return count >= args.num_return;
}

std::vector<index_t> RiskOptimizer::valid_counts_all_pdrs(std::size_t num_threads)
{
    if (num_threads == 0) num_threads = std::thread::hardware_concurrency();

    std::cout << "Validating PDRs...\n";

    if (tmpl.size() < len) {
        std::cerr << "  ERROR: template shorter than PDR length\n";
        return {};
    }

    std::size_t num_pdrs = tmpl.size() - len;
    std::size_t chunk    = (num_pdrs + num_threads - 1) / num_threads;
    std::vector<index_t> results(num_pdrs);
    std::vector<std::future<void>> futures;

    for (std::size_t t = 0; t < num_threads; ++t) {
        std::size_t start = t * chunk;
        std::size_t end   = std::min(start + chunk, num_pdrs);
        if (start >= end) break;

        futures.push_back(std::async(std::launch::async,
            [this, &results, start, end]() {
                for (std::size_t i = start; i < end; ++i)
                    results[i] = valid_count_in_pdr(i);
            }
        ));
    }
    for (auto& f : futures) f.get();

    index_t     min_c    = *std::min_element(results.begin(), results.end());
    index_t     max_c    = *std::max_element(results.begin(), results.end());
    double      avg_c    = (double) std::accumulate(results.begin(), results.end(), 0ULL) / num_pdrs;
    std::size_t zero_cnt = std::count(results.begin(), results.end(), 0);

    std::cout << "Total PDRs             : " << num_pdrs << "\n"
              << "Valid count min/avg/max : " << min_c << " / " << avg_c << " / " << max_c << "\n"
              << "PDRs with 0 valid      : " << zero_cnt
              << "(" << (100.0 * zero_cnt / num_pdrs) << " %)\n\n";

    return results;
}

RiskOptimizer::~RiskOptimizer() {
    delete [] risk;
    delete [] prefix_sum;
    if (max == min) {
        delete [] memo;
        delete [] prev;
    }
    delete [] solutions;
    delete [] solution_vector;
}

std::vector<index_t> RiskOptimizer::random_PDR() {
    std::vector<index_t> PDR;
    while (true) {
        index_t k = PDR.size();
        index_t cmin, cmax;
        if (k % 2 == 0) {
            if (k == 0) {
                cmin = 0;
                cmax = 3 * len - 1;
            }
            else if (k == 2) {
                if (max <= PDR[k - 1] + 3 * len)
                    cmin = std::max(PDR[k - 2] + len, PDR[k - 1] + 3 * len - max);
                else 
                    cmin = PDR[k - 2] + len;
                cmax = PDR[k - 1] - 1;
            }
            else {
                if (max <= PDR[k - 1] + 3 * len)
                    cmin = std::max(PDR[k - 3] + len, PDR[k - 1] + 3 * len - max);
                else 
                    cmin = PDR[k - 3] + len;
                cmax = PDR[k - 1] - 1;
            }
        }
        else {
            if (k == 1) {
                cmin = PDR[k - 1] + min - len;
                cmax = PDR[k - 1] + max - 1;
            }
            else {
                cmin = std::max(PDR[k - 1] - len + min, PDR[k - 2] + 2 * len);
                cmax = PDR[k - 1] + max - 1;
            }
            if (cmax >= size) {
                PDR.pop_back();
                break;
            }
        }
        PDR.push_back(greedy_random_between(cmin, cmax - len + 2));
    }
    return PDR;
}

std::vector<index_t> RiskOptimizer::random_search(std::size_t limit) {
    risk_t min = -1;
    std::vector<index_t> min_PDR;
    std::cout << std::setw(8) << "i" << std::setw(12) << "loss" << std::endl;
    for (std::size_t i = 0; i < limit; i ++) {
        std::vector<index_t> PDR = random_PDR();
        risk_t temp = score(PDR, ALPHA);
        if (min < 0 || temp < min) {
            min = temp;
            min_PDR = PDR;
            std::cout << std::setw(8) << i << std::setw(12) << min << std::endl;
        }
    }
    return min_PDR;
}

index_t RiskOptimizer::greedy_random_between(index_t cmin, index_t cmax) {
    std::vector<std::pair<index_t, risk_t>> all_sum;
    for (index_t i = cmin; i < cmax; i ++) {
        risk_t sum = 0;
        for (index_t j = i; j < i + len; j ++) 
            sum += risk[j];
        all_sum.push_back(std::make_pair(i, sum * sum));
    }
    std::sort(all_sum.begin(), all_sum.end(), [](const std::pair<index_t, risk_t> &a, const std::pair<index_t, risk_t> &b) {
        return a.second < b.second;
    });
    index_t bound = (index_t) ((cmax - cmin) * 0.3);
    if (bound == 0) bound = 1;
    return all_sum[random_between(0, bound)].first;
}

void RiskOptimizer::validate_PDR(std::vector<index_t> PDR) {
    for (index_t i = 0; i < PDR.size(); i ++) {
        if (! gc_valid(PDR[i])) {
            std::cout << i << " " << PDR[i] << "\n";
        }
        assert(gc_valid(PDR[i]));
        if (i % 2 == 0) {
            if (i == 0) {
                continue;
            }
            else if (i == 2) {
                assert(PDR[i] - PDR[i - 2] >= len);
            }
            else {
                assert(PDR[i] - PDR[i - 3] >= len);
                
            }
            index_t l = PDR[i + 1] - PDR[i];
            assert(l + len <= max);
            assert(l + len >= min);
        }
        else {
            if (i == PDR.size() - 1) {
                assert(PDR[i] - PDR[i - 2] >= len);
            }
            else {
                assert(PDR[i] - PDR[i + 1] >= len);
            }
        }

    }
    std::cout << "PDRs: " << PDR.size() << std::endl;
    std::cout << "uncovered: " << PDR[0] << " (front), " << size - PDR[PDR.size() - 1] - len << " (rear)" << std::endl;
    assert(PDR[0] <= len * 2);
    // assert(size - PDR[PDR.size() - 1] <= 3 * len);
    assert(size - PDR[PDR.size() - 1] >= len);
}

risk_t RiskOptimizer::score(std::vector<index_t> PDR, risk_t alpha) {
    std::vector<risk_t> all_sum;
    for (index_t i = 0; i < PDR.size(); i ++) {
        risk_t sum = 0;
        for (index_t j = PDR[i]; j < PDR[i] + len; j ++) 
            sum += risk[j];
        all_sum.push_back(sum * sum);
    }
    std::sort(all_sum.begin(), all_sum.end());
    risk_t total = 0;
    for (index_t i = 0; i < std::ceil(alpha * PDR.size()); i ++) 
        total += all_sum[PDR.size() - i - 1];
    return total;
}

risk_t RiskOptimizer::opt(index_t f, index_t r, index_t f_, index_t r_, risk_t u, risk_t alpha) {
    if (r < r_ + len || r_ < f + len || f < f_ + len) return INF;
    key_t k = to_key(f, r, f_, r_);
    if (memo_umap.find(k) == memo_umap.end()) {
        risk_t min_risk;
        key_t min_key;
        if (f_ <= len * 2) {
            min_risk = cost(f, u, alpha) + cost(r, u, alpha) + cost(f_, u, alpha) + cost(r_, u, alpha);
            min_key = 0;
        }
        else {
            min_risk = INF;
            for (index_t r__ = f_ + len; f >= r__ + len; r__ ++) {
                for (index_t f__ = (r__ >= max ? r__ - max : 0); r__  >= f__ + min; f__ ++) {
                    risk_t temp = opt(f_, r_, f__, r__, u, alpha);
                    if (temp < min_risk) {
                        min_risk = temp;
                        min_key = to_key(f_, r_, f__, r__);
                    }
                }
            }
            min_risk += cost(f, u, alpha) + cost(r, u, alpha);
        }
        memo_umap[k] = min_risk;
        prev_umap[k] = min_key;
    }
    return memo_umap[k];
}

risk_t RiskOptimizer::top_k_opt(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha) {
    min_PDR.clear();
    memo_umap.clear();
    prev_umap.clear();
    risk_t min_risk = INF;
    key_t min_key = 0;
    for (index_t r = size - len * 3; r <= size - len; r ++) {
        for (index_t f = (r >= max ? r - max : 0); r >= f + min; f ++) {
            for (index_t r_ = f + len; r >= r_ + len; r_ ++) {
                for (index_t f_ = (r_ >= max ? r_ - max : 0); r_ >= f_ + min; f_ ++) {
                    risk_t temp = opt(f, r, f_, r_, u, alpha);
                    if (temp < min_risk) {
                        min_risk = temp;
                        min_key = to_key(f, r, f_, r_);
                    }
                }
            }
        }
    }
    while (min_key != 0) {
        index_t a, b, c, d;
        to_index(min_key, a, b, c, d);
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev_umap[min_key];
        if (min_key == 0) {
            min_PDR.push_back(d);
            min_PDR.push_back(c);
        }
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    return min_risk;
}

risk_t RiskOptimizer::top_k_opt_fast(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha) {
    min_PDR.clear();
    for (index_t i = 0; i < size; i ++) 
        solution_vector[i].clear();
    risk_t all_min_risk = INF;
    key_t all_min_key = 0;
    for (index_t r = max; r <= size - len; r ++) {
        index_t f = r - max;
        for (index_t r_ = f + len; r_ <= r - len; r_ ++) {
            if (r_ < max) continue;
            index_t f_ = r_ - max;
            if (r < r_ + len || r_ < f + len || f < f_ + len) continue;
            if (! gc_valid(r) || ! gc_valid(r_) || ! gc_valid(f) || ! gc_valid(f_)) continue;
            key_t k = to_key_2(r, r_);
            risk_t min_risk;
            key_t min_key;
            if (f_ <= fumax) {
                min_risk = cost(f, u, alpha) + cost(r, u, alpha) + cost(f_, u, alpha) + cost(r_, u, alpha);
                min_key = 0;
            }
            else {
                min_risk = INF;
                index_t r__min = f_ + len, r__max = f - len;
                if (r__min < max) r__min = max;
                if (r__min < f_ + len) r__min = f_ + len;
                if (r__max > r_ - len) r__max = r_ - len;
                if (r__min <= r__max) {
                    std::pair<index_t, risk_t> temp = solutions[r_]->range_query(r__min, r__max);
                    min_risk = temp.second;
                    min_key = to_key_2(r_, temp.first);
                }
                /*
                for (auto elem : solution_vector[r_]) {
                    index_t r__ = elem.first;
                    if (r__min <= r__ && r__ <= r__max) {
                        risk_t temp_risk = elem.second;
                        key_t temp_key = to_key_2(r_, r__);
                        if (temp_risk < min_risk) {
                            min_risk = temp_risk;
                            min_key = temp_key;
                        }
                    }
                }
                */
                min_risk += cost(f, u, alpha) + cost(r, u, alpha);
            }
            solution_vector[r].push_back(std::make_pair(r_, min_risk));
            prev[k] = min_key;
            if (r >= size - rumax) {
                if (min_risk < all_min_risk) {
                    all_min_risk = min_risk;
                    all_min_key = k;
                }
            }
        }
        solutions[r] = new BST(solution_vector[r]);
    }
    risk_t min_risk = all_min_risk;
    key_t min_key = all_min_key;
    while (min_key != 0) {
        index_t a, b, c, d;
        to_index_2(min_key, b, d);
        a = b - max;
        c = d - max;
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev[min_key];
        if (min_key == 0) {
            min_PDR.push_back(d);
            min_PDR.push_back(c);
        }
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    index_t i = min_PDR.size();
    for (index_t r = max; r <= size - len; r ++) 
        delete solutions[r];
    return min_risk;
}

risk_t RiskOptimizer::opt_m(index_t f, index_t r, risk_t u, risk_t alpha) {
    key_t k = to_key_2(r, f);
    if (memo_umap.find(k) == memo_umap.end()) {
        risk_t min_risk = INF;
        key_t min_key;
        if (f <= len * 2) {
            min_risk = min_key = 0;
        }
        else {
            index_t r_min = std::max(max, f + len), r_max = (f + r + len) / 2 - len;
            for (index_t r_ = r_min; r_ <= r_max; r_ ++) {
                index_t f_min = (r_ >= max ? r_ - max : 0);
                index_t f_max = std::min((r_ >= min ? r_ - min : 0), (2 * f >= len + r_? 2 * f - len - r_ : 0));
                for (index_t f_ = f_min; f_ <= f_max; f_ ++) {
                    key_t k_ = to_key_2(r_, f_);
                    risk_t temp = opt_m(f_, r_, u, alpha);
                    if (temp < min_risk) {
                        min_risk = temp;
                        min_key = k_;
                    }
                }
            }
        }
        min_risk += cost(f, u, alpha) + cost(r, u, alpha);
        memo_umap[k] = min_risk;
        prev_umap[k] = min_key;
    }
    return memo_umap[k];
}

risk_t RiskOptimizer::top_k_opt_m(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha) {
    min_PDR.clear();
    memo_umap.clear();
    prev_umap.clear();
    risk_t min_risk = INF;
    key_t min_key = 0;
    for (index_t r = size - len * 3; r <= size - len; r ++) {
        for (index_t f = (r >= max ? r - max : 0); r >= f + min; f ++) {
            risk_t temp = opt_m(f, r, u, alpha);
            if (temp < min_risk) {
                min_risk = temp;
                min_key = to_key_2(r, f);
            }
        }
    }
    while (min_key != 0) {
        index_t a, b;
        to_index_2(min_key, b, a);
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev_umap[min_key];
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    return min_risk;
}

risk_t RiskOptimizer::top_k_opt_mi(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha) {
    min_PDR.clear();
    risk_t all_min_risk = INF;
    key_t all_min_key = 0;
    for (index_t r = max; r <= size - len; r ++) {
        for (index_t f = (r >= max ? r - max : 0); r >= f + min; f ++) {
            risk_t min_risk = INF;
            key_t min_key, k = to_key_2(r, f);
            if (f <= len * 2) {
                min_risk = min_key = 0;
            }
            else {
                index_t r_min = std::max(max, f + len), r_max = (f + r + len) / 2 - len;
                for (index_t r_ = r_min; r_ <= r_max; r_ ++) {
                    index_t f_min = (r_ >= max ? r_ - max : 0);
                    index_t f_max = std::min((r_ >= min ? r_ - min : 0), (2 * f >= len + r_? 2 * f - len - r_ : 0));
                    for (index_t f_ = f_min; f_ <= f_max; f_ ++) {
                        key_t k_ = to_key_2(r_, f_);
                        risk_t temp = memo[k_];
                        if (temp < min_risk) {
                            min_risk = temp;
                            min_key = k_;
                        }
                    }
                }
            }
            min_risk += cost(f, u, alpha) + cost(r, u, alpha);
            memo[k] = min_risk;
            prev[k] = min_key;
            if (r >= size - len * 3) {
                if (min_risk < all_min_risk) {
                    all_min_risk = min_risk;
                    all_min_key = k;
                }
            }
        }
    }
    risk_t min_risk = all_min_risk;
    key_t min_key = all_min_key;
    while (min_key != 0) {
        index_t a, b;
        to_index_2(min_key, b, a);
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev[min_key];
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    return min_risk;
}

risk_t RiskOptimizer::top_k_opt_mi_fast(risk_t u,
                                    std::vector<index_t>& min_PDR,
                                    risk_t alpha) {
    min_PDR.clear();

    const index_t band = max - min + 1;
    init_row_prefix(band);  // reuse allocation across calls

    risk_t all_min_risk = INF;
    key_t  all_min_key  = 0;

    auto f_base_fn = [&](index_t r_) -> index_t {
        return (r_ >= (index_t)max ? r_ - max : 0);
    };

    for (index_t r = (index_t)max; r <= size - (index_t)len; r++) {
        index_t fb_r = f_base_fn(r);

        for (index_t f = fb_r; r >= f + (index_t)min; f++) {
            if (!gc_valid(r) || !gc_valid(f)) continue;

            risk_t min_risk = INF;
            key_t  min_key  = 0;
            key_t  k        = to_key_2(r, f);

            if (f <= (index_t)(len * 2)) {
                min_risk = min_key = 0;
            } else {
                index_t r_min = std::max((index_t)max, f + (index_t)len);
                index_t r_max = (f + r + (index_t)len) / 2 - (index_t)len;

                for (index_t r_ = r_min; r_ <= r_max; r_++) {
                    if (!gc_valid(r_)) continue;
                    index_t fb_r_ = f_base_fn(r_);

                    // exact fi_max for this r_
                    index_t f_max_ = std::min(
                        (r_ >= (index_t)min ? r_ - (index_t)min : (index_t)0),
                        (2*f >= (index_t)(len + r_) ?
                         2*f - (index_t)len - r_ : (index_t)0)
                    );

                    if (f_max_ < fb_r_) continue;
                    index_t fi_max_ = std::min(f_max_ - fb_r_, band - 1);

                    // O(1) prefix min lookup for row r_
                    if (rp(r_, fi_max_) < min_risk) {
                        min_risk = rp(r_,  fi_max_);
                        min_key  = rpk(r_, fi_max_);
                    }
                }
            }

            min_risk += cost(f, u, alpha) + cost(r, u, alpha);
            memo[k]   = min_risk;
            prev[k]   = min_key;

            // update row prefix for row r at fi
            index_t fi = f - fb_r;
            if (min_risk < rp(r, fi)) {
                rp(r,  fi) = min_risk;
                rpk(r, fi) = k;
            }

            if (r >= (index_t)(size - (index_t)len * 3)) {
                if (min_risk < all_min_risk) {
                    all_min_risk = min_risk;
                    all_min_key  = k;
                }
            }
        }

        // build prefix min along fi for completed row r
        index_t fb_r_ = f_base_fn(r);
        for (index_t fi = 1; fi < band; fi++) {
            if (rp(r, fi-1) < rp(r, fi)) {
                rp(r,  fi) = rp(r,  fi-1);
                rpk(r, fi) = rpk(r, fi-1);
            }
        }
    }

    // traceback
    key_t min_key = all_min_key;
    while (min_key != 0) {
        index_t a, b;
        to_index_2(min_key, b, a);
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev[min_key];
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    return all_min_risk;
}

std::vector<index_t> RiskOptimizer::search(risk_t lower, 
                                           risk_t upper, 
                                           risk_t eps) {
    const int TW = 8 + 12 + 12 + 12; // = 44
    std::cout << "Search for optimal PDR...\n"
              << std::string(TW, '-') << "\n"
              << std::right
              << std::setw(8)  << "iter"
              << std::setw(12) << "u'"
              << std::setw(12) << "loss'"
              << std::setw(12) << "loss"
              << "\n"
              << std::string(TW, '-') << "\n";

    max -= len; min -= len;
    std::vector<index_t> min_PDR, final_PDR, temp_PDR;
    risk_t alpha = ALPHA, l = lower, r = upper;

    for (index_t i = 1; r - l > eps; i++) {
        risk_t u = (l + r) / 2, temp;
        if (min == max)
            temp = top_k_opt_fast(u, min_PDR, alpha);
        else
            temp = top_k_opt_mi_fast(u, min_PDR, alpha);

        if (final_PDR.size() == 0 || score(min_PDR, alpha) < score(final_PDR, alpha))
            final_PDR = min_PDR;

        risk_t grad = alpha * min_PDR.size();
        for (index_t j = 0; j < min_PDR.size(); j++)
            if (u < cost(min_PDR[j], 0, 0)) grad--;
        (grad < 0 ? l : r) = u;

        std::cout << std::right
                  << std::setw(8)  << i
                  << std::setw(12) << u
                  << std::setw(12) << temp
                  << std::setw(12) << score(min_PDR, alpha)
                  << "\n";
    }

    std::cout << std::string(TW, '-') << "\n";
    max += len; min += len;
    std::cout << "Search completed!\n\n";
    return final_PDR;
}

