#include "pipeline.hpp"

void PDRStage::run(PipelineContext& ctx) {
    srand(ctx.args.seed);

    RiskOptimizer ro(ctx.labels,
                     ctx.sequences, 
                     ctx.args.len_PDR, 
                     ctx.args.len_amp, 
                     ctx.args.len_amp,
                     ctx.args.len_fu,
                     ctx.args.len_ru,
                     ctx.args);

    auto PDR = ro.search(0, ctx.args.u_max, ctx.args.u_min);

    ro.validate_PDR(PDR);
    
    risk_t score = ro.score(PDR, ALPHA);
    std::cout << "loss: " << score << std::endl;

    std::ofstream fout(ctx.args.output_file + "." + short_name());
    fout << display_PDR(PDR) << std::endl;
    fout.close();

    ctx.pdr_regions = PDR;
}

void PrimerSelStage::run(PipelineContext& ctx) {
    std::cout << "Starting primer selection for "
              << ctx.pdr_regions.size() / 2 << " PDRs\n";

    // --- 1. Global settings ---
    p3_global_settings *pa = p3_create_global_settings();
    p3_set_gs_primer_opt_size(pa, ctx.args.p3_opt_size);
    p3_set_gs_primer_min_size(pa, ctx.args.p3_min_size);
    p3_set_gs_primer_max_size(pa, ctx.args.p3_max_size);
    p3_set_gs_primer_opt_tm(pa, ctx.args.p3_opt_tm);
    p3_set_gs_primer_min_tm(pa, ctx.args.p3_min_tm);
    p3_set_gs_primer_max_tm(pa, ctx.args.p3_max_tm);
    p3_set_gs_primer_min_gc(pa, ctx.args.p3_min_gc);
    p3_set_gs_primer_max_gc(pa, ctx.args.p3_max_gc);

    pa->p_args.salt_conc     = ctx.args.mv;
    pa->p_args.divalent_conc = ctx.args.dv;
    pa->p_args.dntp_conc     = ctx.args.dntp;
    pa->p_args.dna_conc      = ctx.args.dna_conc;

    pa->pr_min[0] = ctx.args.len_amp_min;
    pa->pr_max[0] = ctx.args.len_amp + ctx.args.len_PDR;
    pa->num_intervals = 1;

    pa->pick_left_primer    = 1;
    pa->pick_right_primer   = 1;
    pa->pick_internal_oligo = 0;
    pa->num_return          = ctx.args.num_return;

    std::cout << "Product size: [" << pa->pr_min[0] 
              << ", " << pa->pr_max[0] << "] bp"
              << "  |  Template length: " << ctx.tmpl.size() << " bp\n";

    // --- 2. Sequence arguments ---
    seq_args *sa = create_seq_arg();
    sa->sequence      = strdup(ctx.tmpl.c_str());
    sa->sequence_name = strdup("my_template");

    sa->ok_regions.any_left  = 0;
    sa->ok_regions.any_right = 0;
    sa->ok_regions.any_pair  = 0;
    sa->ok_regions.count     = 1;
    sa->ok_regions.left_pairs[0][1]  = ctx.args.len_PDR;
    sa->ok_regions.right_pairs[0][1] = ctx.args.len_PDR;

    // --- 3. Run primer design ---
    const int TW = 10 + 20 + 10 + 10 + 10; // = 60
    std::cout << std::string(TW, '-') << "\n"
              << std::right
              << std::setw(10) << "Index"
              << std::setw(20) << "PDRs"
              << std::setw(10) << "Left"
              << std::setw(10) << "Right"
              << std::setw(10) << "Pairs"
              << "\n"
              << std::string(TW, '-') << "\n";

    for (std::size_t i = 0; i < ctx.pdr_regions.size(); i += 2) {
        int region_idx = i / 2;
        sa->ok_regions.left_pairs[0][0]  = ctx.pdr_regions[i];
        sa->ok_regions.right_pairs[0][0] = ctx.pdr_regions[i + 1];

        p3retval *retval = choose_primers(pa, sa);

        if (retval == nullptr) {
            std::cerr << "ERROR: choose_primers() returned NULL for region "
                      << region_idx << "\n";
            return;
        }
        if (retval->glob_err.data)
            std::cerr << "Global error: "   
                      << retval->glob_err.data << "\n";
        if (retval->per_sequence_err.data)
            std::cerr << "Sequence error: " 
                      << retval->per_sequence_err.data << "\n";

        PrimerOutput out  = extract_all(retval, sa);
        out.pdr_left      = ctx.pdr_regions[i];
        out.pdr_right     = ctx.pdr_regions[i + 1];
        out.segment_id = ctx.current_seg_id();
        ctx.candidate_primers.push_back(out);

        for (const Oligo& o : out.left_oligos) {
            std::string tmpl_sub = ctx.tmpl.substr(o.start, o.length);
            if (tmpl_sub != o.seq)
                std::cout << "  LEFT  MISMATCH start=" << o.start
                          << " len=" << o.length
                          << " seq=  " << o.seq
                          << " tmpl= " << tmpl_sub
                          << "\n";
        }
        for (const Oligo& o : out.right_oligos) {
            std::string tmpl_sub = ctx.tmpl.substr(o.start, o.length);
            if (tmpl_sub != o.seq)
                std::cout << "  RIGHT MISMATCH start=" << o.start
                          << " len=" << o.length
                          << " seq=  " << o.seq
                          << " tmpl= " << tmpl_sub
                          << "\n";
        }

        int n_pairs = retval->best_pairs.num_pairs;
        std::string range = "[" + std::to_string(ctx.pdr_regions[i]) +
                            ", " + std::to_string(ctx.pdr_regions[i+1]) + "]";
        std::cout << std::right
                  << std::setw(10) << region_idx
                  << std::setw(20) << range
                  << std::setw(10) << retval->fwd.num_elem
                  << std::setw(10) << retval->rev.num_elem
                  << std::setw(10) << n_pairs
                  << "\n";

        destroy_p3retval(retval);
    }

    std::cout << std::string(TW, '-') << "\n"
              << "Total PDRs processed: " 
              << ctx.pdr_regions.size() / 2 << "\n";

    destroy_seq_args(sa);
    p3_destroy_global_settings(pa);
}

void OffTargetStage::run(PipelineContext& ctx) {
    if (ctx.args.ref_file == "") {
        std::cout << "No ref file provided, skip off target search\n";
        for (auto& seg : ctx.segments)
            seg.filtered_primers = seg.candidate_primers;
        return;
    }

    Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config",
               ctx.args.mv, ctx.args.dv, ctx.args.dntp,
               ctx.args.dna_conc, ctx.args.temp);

    // flatten all candidate_primers from all segments — track offsets
    std::vector<PrimerOutput> all_candidates;
    std::vector<std::size_t>  seg_offsets;  // seg_offsets[s] = start index in all_candidates

    for (auto& seg : ctx.segments) {
        seg_offsets.push_back(all_candidates.size());
        for (auto& po : seg.candidate_primers)
            all_candidates.push_back(po);
    }

    std::cout << "Off-target search: " << all_candidates.size()
              << " PDRs across " << ctx.segments.size() << " segments\n";

    // single automaton + single genome scan
    std::vector<std::string> labels, sequences;
    convert(all_candidates, labels, sequences);

    Automaton *ac = new Automaton(labels, sequences, ctx.args.kmer_len);
    auto results = ac->search(ctx.args.ref_file,
                              ctx.args.kmer_len - 1,
                              ctx.args.threshold,
                              ctx.args.dg_thres,
                              ctx.args.chunk_size,
                              ctx.args.block_size,
                              ctx.args.nthreads);
    delete ac;

    write_results(ctx.args.output_file + "." + short_name(), results, labels);

    // single filterByDG_relax call on the flat list
    auto all_filtered = filterByDG_relax(results, ctx.args.dg_thres,
                                         all_candidates, ctx.args.num_return);

    // split back to segments using offsets
    for (std::size_t s = 0; s < ctx.segments.size(); ++s) {
        std::size_t start = seg_offsets[s];
        std::size_t end   = (s + 1 < ctx.segments.size())
                          ? seg_offsets[s + 1]
                          : all_filtered.size();

        ctx.segments[s].filtered_primers.clear();
        for (std::size_t j = start; j < end; ++j) {
            PrimerOutput po    = all_filtered[j];
            // restore metadata lost by filterByDG_relax
            po.pdr_left        = all_candidates[j].pdr_left;
            po.pdr_right       = all_candidates[j].pdr_right;
            po.segment_id      = all_candidates[j].segment_id;
            ctx.segments[s].filtered_primers.push_back(po);
        }
    }
}

/*
void DimerStage::run(PipelineContext& ctx) {
    srand(ctx.args.seed);

    Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config", 
        ctx.args.mv, 
        ctx.args.dv, 
        ctx.args.dntp, 
        ctx.args.dna_conc, 
        ctx.args.temp);

    KPartiteGraph g(ctx.filtered_primers);

    using Clock = std::chrono::steady_clock;
    using Ms    = std::chrono::milliseconds;

    struct Result {
        std::string          name;
        std::vector<index_t> solution;
        long                 ms;
        weight_t             cost;
    };

    auto timed_run = [&](const std::string& name, auto fn) -> Result {
        auto t0  = Clock::now();
        auto sol = fn();
        auto ms  = std::chrono::duration_cast<Ms>(Clock::now() - t0).count();
        return { name, sol, ms, g.cost(sol) };
    };

    weight_t bottleneck_threshold = 0;

    std::vector<Result> results = {
        timed_run("Bottleneck",    [&]{ return g.solve_bottleneck(ctx.args.iter * 10,
                                                               bottleneck_threshold);        }),
        // timed_run("Random Search", [&]{ return g.solve_fast(ctx.args.iter);    }),
        // timed_run("SA",            [&]{ return g.solve_sa(ctx.args.iter);             }),
        // timed_run("Tabu Search",   [&]{ return g.solve_tabu(20, ctx.args.iter);  }),
        // timed_run("Genetic",       [&]{ return g.solve_ga(1000, ctx.args.iter * 10, 0.05); }),
    };

    // summary table
    const int TW = 16 + 14 + 10 + 10;
    std::cout << "\n" << std::string(TW, '-') << "\n"
              << std::left  << std::setw(16) << "Algorithm"
              << std::right << std::setw(14) << "Cost"
                            << std::setw(10) << "Time(ms)"
                            << std::setw(10) << "Winner"
              << "\n" << std::string(TW, '-') << "\n";

    auto& winner = *std::max_element(results.begin(), results.end(),
        [](const Result& a, const Result& b){ return a.cost < b.cost; });

    for (auto& r : results) {
        std::cout << std::left  << std::setw(16) << r.name
                  << std::right << std::fixed << std::setprecision(4)
                                << std::setw(14) << r.cost
                                << std::setw(10) << r.ms
                                << std::setw(10) << (&r == &winner ? "✓" : "")
                  << "\n";
    }

    std::cout << std::string(TW, '-') << "\n"
              << " Winner: " << winner.name 
              << "  cost=" << winner.cost
              << "  (" << winner.ms << " ms)\n";

    ctx.dimer_solution = results[0].solution;

    for (std::size_t amp = 0; amp < ctx.filtered_primers.size(); amp++) {
        index_t left_n  = ctx.dimer_solution[amp * 2];
        index_t right_n = ctx.dimer_solution[amp * 2 + 1];

        const PrimerOutput& po = ctx.filtered_primers[amp];
        PrimerResult result;
        result.left        = po.left_oligos[left_n];
        result.right       = po.right_oligos[right_n];
        result.product_size = result.right.start - result.left.start + result.right.length;
        result.pdr_left      = po.pdr_left;
        result.pdr_right     = po.pdr_right;
        result.segment_id = po.segment_id;

        ctx.solution_primers.push_back(result);
    }
}
*/

void DimerStage::run(PipelineContext& ctx) {
    srand(ctx.args.seed);

    Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config", 
        ctx.args.mv, 
        ctx.args.dv, 
        ctx.args.dntp, 
        ctx.args.dna_conc, 
        ctx.args.temp);

    using Clock = std::chrono::steady_clock;
    using Ms    = std::chrono::milliseconds;

    // Split filtered_primers into two pools by alternating per segment
    // seg_order[segment_id] = running count so we can alternate
    std::unordered_map<std::size_t, int> seg_order;
    std::vector<PrimerOutput> pool1_primers, pool2_primers;
    // track original index so we can reconstruct later
    std::vector<std::size_t>  pool1_idx,    pool2_idx;

    for (std::size_t i = 0; i < ctx.filtered_primers.size(); i++) {
        const PrimerOutput& po = ctx.filtered_primers[i];
        int& cnt = seg_order[po.segment_id];
        if (cnt % 2 == 0) {
            pool1_primers.push_back(po);
            pool1_idx.push_back(i);
        } else {
            pool2_primers.push_back(po);
            pool2_idx.push_back(i);
        }
        cnt++;
    }

    struct Result {
        std::string          name;
        std::vector<index_t> solution;
        long                 ms;
        weight_t             cost;
    };

    auto solve_pool = [&](const std::vector<PrimerOutput>& primers,
                      int pool_num) -> std::vector<index_t>
    {
        KPartiteGraph g(primers);

        weight_t bn_threshold   = 0;
        weight_t tabu_threshold = 0;
        weight_t sa_threshold   = 0;

        // --- Bottleneck Search ---
        auto t0     = Clock::now();
        auto bn_sol = g.solve_bottleneck(ctx.args.iter * 100, bn_threshold);
        auto bn_ms  = std::chrono::duration_cast<Ms>(Clock::now() - t0).count();

        // --- Tabu Search ---
        auto t1       = Clock::now();
        auto tabu_sol = g.solve_bottleneck_tabu(20, g.get_K() / 3, ctx.args.iter / 2, tabu_threshold);
        auto tabu_ms  = std::chrono::duration_cast<Ms>(Clock::now() - t1).count();

        // --- SA Search ---
        // init_temp scaled to score range, cooling chosen to freeze after max_iter steps
        double      init_temp = 1000.0;
        std::size_t sa_iter   = ctx.args.iter * 5;
        double      cooling   = std::exp(std::log(1e-6 / init_temp) / (double)sa_iter);
        auto t2     = Clock::now();
        auto sa_sol = g.solve_bottleneck_sa(4000, sa_iter, init_temp, cooling, sa_threshold);
        auto sa_ms  = std::chrono::duration_cast<Ms>(Clock::now() - t2).count();

        // --- Pick winner ---
        weight_t best = std::max({ bn_threshold, tabu_threshold, sa_threshold });

        // --- Summary table ---
        const int TW = 16 + 14 + 10 + 10;
        std::cout << "\nPool " << pool_num << " Summary\n"
                << std::string(TW, '-') << "\n"
                << std::left  << std::setw(16) << "Method"
                << std::right << std::setw(14) << "Bottleneck"
                                << std::setw(10) << "Time(ms)"
                                << std::setw(10) << "Winner"
                << "\n"
                << std::string(TW, '-') << "\n"
                << std::left  << std::setw(16) << "Bottleneck"
                << std::right << std::fixed << std::setprecision(4)
                                << std::setw(14) << bn_threshold
                                << std::setw(10) << bn_ms
                                << std::setw(10) << (bn_threshold == best ? "✓" : "")
                << "\n"
                << std::left  << std::setw(16) << "Tabu"
                << std::right << std::fixed << std::setprecision(4)
                                << std::setw(14) << tabu_threshold
                                << std::setw(10) << tabu_ms
                                << std::setw(10) << (tabu_threshold == best ? "✓" : "")
                << "\n"
                << std::left  << std::setw(16) << "SA"
                << std::right << std::fixed << std::setprecision(4)
                                << std::setw(14) << sa_threshold
                                << std::setw(10) << sa_ms
                                << std::setw(10) << (sa_threshold == best ? "✓" : "")
                << "\n"
                << std::string(TW, '-') << "\n";

        // Return the best solution
        if (sa_threshold == best) {
            std::cout << " Winner: SA  (cost=" << best << ")\n";
            return sa_sol;
        } else if (tabu_threshold == best) {
            std::cout << " Winner: Tabu  (cost=" << best << ")\n";
            return tabu_sol;
        } else {
            std::cout << " Winner: Bottleneck  (cost=" << best << ")\n";
            return bn_sol;
        }
    };

    auto sol1 = solve_pool(pool1_primers, 1);
    auto sol2 = solve_pool(pool2_primers, 2);

    // Reconstruct PrimerResult for pool 1
    for (std::size_t amp = 0; amp < pool1_primers.size(); amp++) {
        index_t left_n  = sol1[amp * 2];
        index_t right_n = sol1[amp * 2 + 1];

        const PrimerOutput& po = pool1_primers[amp];
        PrimerResult result;
        result.left         = po.left_oligos[left_n];
        result.right        = po.right_oligos[right_n];
        result.product_size = result.right.start - result.left.start + result.right.length;
        result.pdr_left     = po.pdr_left;
        result.pdr_right    = po.pdr_right;
        result.segment_id   = po.segment_id;
        result.pool         = 1;

        ctx.solution_primers.push_back(result);
    }

    // Reconstruct PrimerResult for pool 2
    for (std::size_t amp = 0; amp < pool2_primers.size(); amp++) {
        index_t left_n  = sol2[amp * 2];
        index_t right_n = sol2[amp * 2 + 1];

        const PrimerOutput& po = pool2_primers[amp];
        PrimerResult result;
        result.left         = po.left_oligos[left_n];
        result.right        = po.right_oligos[right_n];
        result.product_size = result.right.start - result.left.start + result.right.length;
        result.pdr_left     = po.pdr_left;
        result.pdr_right    = po.pdr_right;
        result.segment_id   = po.segment_id;
        result.pool         = 2;

        ctx.solution_primers.push_back(result);
    }

    // Sort solution_primers back to original amplicon order
    std::vector<std::size_t> all_idx;
    all_idx.insert(all_idx.end(), pool1_idx.begin(), pool1_idx.end());
    all_idx.insert(all_idx.end(), pool2_idx.begin(), pool2_idx.end());

    std::vector<PrimerResult> sorted(ctx.solution_primers.size());
    for (std::size_t i = 0; i < all_idx.size(); i++)
        sorted[all_idx[i]] = ctx.solution_primers[i];
    ctx.solution_primers = std::move(sorted);
}

void Pipeline::run(PipelineContext& ctx) {
    const int W = 80;
    std::cout << "\n"
              << std::string(W, '=') << "\n"
              << std::string((W - 26) / 2, ' ') 
              << "PIPELINE EXECUTION STARTED\n"
              << std::string(W, '=') << "\n";

    for (size_t i = 0; i < stages_.size(); ++i) {

        if (stages_[i]->short_name() == "dim") {
            // flatten all filtered_primers from all segments
            ctx.filtered_primers.clear();
            for (auto& seg : ctx.segments)
                ctx.filtered_primers.insert(ctx.filtered_primers.end(),
                    seg.filtered_primers.begin(),
                    seg.filtered_primers.end());

            std::string running   = "[ " + std::to_string(i+1) + "/" +
                                    std::to_string(stages_.size()) +
                                    " Running: " + stages_[i]->name() + " (global) ]";
            std::string completed = "[ Completed: " + stages_[i]->name() + " ]";
            int rpad = (W - running.size())   / 2;
            int cpad = (W - completed.size()) / 2;

            std::cout << "\n\n"
                      << std::string(rpad, '=') << running
                      << std::string(W - rpad - running.size(), '=') << "\n";

            stages_[i]->run(ctx);

            std::cout << std::string(cpad, '=') << completed
                      << std::string(W - cpad - completed.size(), '=') << "\n\n";

        } else if (stages_[i]->short_name() == "offt") {
            // global — operates directly on ctx.segments, no flatten needed
            std::string running   = "[ " + std::to_string(i+1) + "/" +
                                    std::to_string(stages_.size()) +
                                    " Running: " + stages_[i]->name() + " (global) ]";
            std::string completed = "[ Completed: " + stages_[i]->name() + " ]";
            int rpad = (W - running.size())   / 2;
            int cpad = (W - completed.size()) / 2;

            std::cout << "\n\n"
                      << std::string(rpad, '=') << running
                      << std::string(W - rpad - running.size(), '=') << "\n";

            stages_[i]->run(ctx);

            std::cout << std::string(cpad, '=') << completed
                      << std::string(W - cpad - completed.size(), '=') << "\n\n";

        } else {
            for (std::size_t s = 0; s < ctx.segments.size(); ++s) {
                ctx.load_segment(s);

                std::string seg_tag   = ctx.segments[s].name +
                                        " [" + std::to_string(s+1) +
                                        "/" + std::to_string(ctx.segments.size()) + "]";
                std::string running   = "[ " + std::to_string(i+1) + "/" +
                                        std::to_string(stages_.size()) +
                                        " Running: " + stages_[i]->name() +
                                        " | " + seg_tag + " ]";
                std::string completed = "[ Completed: " + stages_[i]->name() +
                                        " | " + seg_tag + " ]";
                int rpad = (W - running.size())   / 2;
                int cpad = (W - completed.size()) / 2;

                std::cout << "\n\n"
                          << std::string(rpad, '=') << running
                          << std::string(W - rpad - running.size(), '=') << "\n";

                stages_[i]->run(ctx);
                ctx.save_segment();

                std::cout << std::string(cpad, '=') << completed
                          << std::string(W - cpad - completed.size(), '=') << "\n\n";
            }
        }
    }
}