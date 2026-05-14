#include "automaton.hpp"

AutomatonNode::AutomatonNode() {
    for (std::size_t i = 0; i < ALPHABET_SIZE; i ++) 
        children[i] = NULL;
    parent = NULL;
    prev = NULL;
}

AutomatonNode::~AutomatonNode() {
    for (auto node : children) {
        if (node != NULL) 
            delete node;
    }
}

Automaton::Automaton(const std::vector<std::string> &labels, 
                     const std::vector<std::string> &seqs, 
                     std::size_t kmer_len) : k(kmer_len) {
    auto parse_label = [](const std::string& label, 
                        std::size_t& pdr_id, 
                        bool& is_left) {
        pdr_id  = 0;
        is_left = true;
        if (label.size() < 3 || label.substr(0,3) != "OUT") return;
        std::size_t p = 3;
        while (p < label.size() && std::isdigit(label[p])) p++;
        pdr_id = std::stoull(label.substr(3, p - 3));
        // find _L or _R or _P..._F/_R
        std::size_t l_pos = label.find("_L");
        std::size_t r_pos = label.find("_R");
        std::size_t f_pos = label.find("_F");
        if (l_pos != std::string::npos || f_pos != std::string::npos)
            is_left = true;
        else if (r_pos != std::string::npos)
            is_left = false;
    };
    std::vector<std::string> data;
    for (std::size_t i = 0; i < labels.size(); i ++) {
        std::size_t pdr_id; bool is_left;
        parse_label(labels[i], pdr_id, is_left);
        const std::string &forward = seqs[i], &reverse = revcomp(seqs[i]);
        for (std::size_t j = 0; j + k <= forward.size(); j ++) {
            data.push_back(forward.substr(j, k));
            index2label.push_back(labels[i]);
            index2data.push_back(forward);
            index2pos.push_back(j);
            index2strand.push_back(i * 2);
            index2pdr.push_back(pdr_id);
            index2lr.push_back(is_left);
        }
        for (std::size_t j = 0; j + k <= reverse.size(); j ++) {
            data.push_back(reverse.substr(j, k));
            index2label.push_back(labels[i]);
            index2data.push_back(reverse);
            index2pos.push_back(j);
            index2strand.push_back(i * 2 + 1);
            index2pdr.push_back(pdr_id);
            index2lr.push_back(is_left);
        }
        label2index[labels[i]] = i;
    }

    root = new AutomatonNode();
    for (std::size_t i = 0; i < data.size(); i ++) {
        AutomatonNode *current = root;
        for (std::size_t j = 0; j < data[i].size(); j ++) {
            if (current->children[data[i][j]] == NULL) {
                AutomatonNode *next = new AutomatonNode();
                current->children[data[i][j]] = next;
                next->parent = current;
                next->ch = data[i][j];
            }
            current = current->children[data[i][j]];
        }
        current->index.push_back(i);
    }

    std::vector<AutomatonNode *> queue;
    std::size_t head = 0;
    for (std::size_t i = 0; i < ALPHABET_SIZE; i ++) {
        if (root->children[i] != NULL) {
            queue.push_back(root->children[i]);
            root->children[i]->prev = root;
        }
    }

    while (head < queue.size()) {
        AutomatonNode *current = queue[head];
        for (std::size_t i = 0; i < ALPHABET_SIZE; i ++) {
            if (current->children[i] != NULL) {
                AutomatonNode *temp = current->prev;
                while (temp != root && temp->children[i] == NULL)
                    temp = temp->prev;
                if (temp->children[i] != NULL) 
                    current->children[i]->prev = temp->children[i];
                else 
                    current->children[i]->prev = root;
                queue.push_back(current->children[i]);
            }
        }
        head ++;
    }
}

Automaton::~Automaton() {
    delete root;
}

std::vector<Result> Automaton::search(const std::string &path, 
                                      std::size_t overlap, 
                                      std::size_t threshold,
                                      double dg_thres,
                                      std::size_t chunk_size, 
                                      std::size_t io_block_size, 
                                      std::size_t nthreads) {
    return process_fasta_chunks(path, 
                                overlap,
                                threshold,
                                dg_thres,
                                chunk_size,
                                io_block_size,
                                nthreads);
}

std::vector<Result> Automaton::search_paired(const std::string &path, 
                                      std::size_t overlap, 
                                      std::size_t threshold,
                                      double dg_thres,
                                      std::size_t amp_min,
                                      std::size_t amp_max,
                                      std::size_t chunk_size, 
                                      std::size_t io_block_size, 
                                      std::size_t nthreads) {
    return process_fasta_chunks_paired(path, 
                                overlap,
                                threshold,
                                dg_thres,
                                amp_min,
                                amp_max,
                                chunk_size,
                                io_block_size,
                                nthreads);
}

std::size_t Automaton::search_chunk(const std::string& contig,
                                    const char* data,
                                    std::size_t len,
                                    std::size_t index,
                                    std::size_t threshold,
                                    double dg_thres,
                                    std::vector<Result> &local_results) {
    /*
    std::cout << contig << "\t"
              << index + len - this->k + 1
              << "\n";
    */
    
    std::size_t count = 0;
    AutomatonNode *current = root;
    std::unordered_set<CandidateKey, CandidateKeyHash> active;
    std::deque<ActiveEntry> q;

    for (std::size_t i = 0, j = index; i < len; i ++, j ++) {
        while (!q.empty() && q.front().expire_at < i) {
            active.erase(q.front().key);
            q.pop_front();
        }

        while (current != root && current->children[data[i]] == NULL) 
            current = current->prev;
        if (current->children[data[i]] != NULL) {
            current = current->children[data[i]];
        }

        for (std::size_t k = 0; k < current->index.size(); k ++) {
            std::size_t seq_id = current->index[k], pos = index2pos[seq_id];
            std::string seq = index2data[seq_id], seq_ = seq;
            if (this->k + pos > i + 1) continue;

            std::size_t start = i + 1 - this->k - pos;
            if (start + seq.size() > len) continue;

            CandidateKey key{index2strand[seq_id], start};
            if (active.find(key) != active.end()) continue;
            active.insert(key);
            q.push_back(ActiveEntry{key, start + seq.size() - 1});

            std::size_t mismatch[2] = {0, 0};
            std::size_t i_, j_;
            for (i_ = 0, j_ = start; i_ < seq.size() - 1; i_ ++, j_ ++) {
                if (data[j_] != seq[i_]) {
                    ++ mismatch[0];
                    ++ mismatch[1];
                    seq_[i_] = data[j_];
                }
                else if (data[j_ + 1] != seq[i_ + 1]) {
                    ++ mismatch[1];
                }
            }
            if (data[j_] != seq[i_]) {
                ++ mismatch[0];
                seq_[i_] = data[j_];
            }
            if (mismatch[threshold % 2] > threshold / 2) continue;

            double dg = Thal::compute_dimer_dg(seq, seq_);
            if (dg < dg_thres) {
                local_results.emplace_back(
                    Result{
                        index2label[seq_id],
                        seq,
                        contig,
                        seq_,
                        i + 1 - this->k - pos,
                        index2strand[seq_id] % 2,
                        dg
                    }
                );
            }
            count ++;
        }
    }
    return count;
}

std::size_t Automaton::search_chunk_paired(const std::string& contig,
                                            const char* data,
                                            std::size_t len,
                                            std::size_t index,
                                            std::size_t threshold,
                                            double dg_thres,
                                            std::size_t amp_min,
                                            std::size_t amp_max,
                                            std::vector<Result>& local_results) {
    std::size_t count = 0;
    AutomatonNode* current = root;
    std::unordered_set<CandidateKey, CandidateKeyHash> active;
    std::deque<ActiveEntry> q;

    // pending left hits per PDR — deque sorted by pos
    struct PendingHit {
        std::size_t pos;
        std::size_t seq_id;
        std::string seq;
        std::string seq_;
    };
    std::unordered_map<std::size_t, std::deque<PendingHit>> left_pending;

    for (std::size_t i = 0, j = index; i < len; i++, j++) {
        // evict active deduplication entries
        while (!q.empty() && q.front().expire_at < i) {
            active.erase(q.front().key);
            q.pop_front();
        }

        // evict left_pending entries outside amp_max window
        for (auto& [pdr_id, dq] : left_pending)
            while (!dq.empty() && dq.front().pos + amp_max < i)
                dq.pop_front();

        // advance automaton
        while (current != root && current->children[data[i]] == nullptr)
            current = current->prev;
        if (current->children[data[i]] != nullptr)
            current = current->children[data[i]];

        for (std::size_t k = 0; k < current->index.size(); k++) {
            std::size_t seq_id = current->index[k];
            std::size_t pos    = index2pos[seq_id];
            std::string seq    = index2data[seq_id];
            std::string seq_   = seq;

            if (this->k + pos > i + 1) continue;
            std::size_t start = i + 1 - this->k - pos;
            if (start + seq.size() > len) continue;

            CandidateKey key{index2strand[seq_id], start};
            if (active.find(key) != active.end()) continue;
            active.insert(key);
            q.push_back(ActiveEntry{key, start + seq.size() - 1});

            // Hamming check — unchanged
            std::size_t mismatch[2] = {0, 0};
            std::size_t i_, j_;
            for (i_ = 0, j_ = start; i_ < seq.size() - 1; i_++, j_++) {
                if (data[j_] != seq[i_]) {
                    ++mismatch[0]; ++mismatch[1];
                    seq_[i_] = data[j_];
                } else if (data[j_ + 1] != seq[i_ + 1]) {
                    ++mismatch[1];
                }
            }
            if (data[j_] != seq[i_]) { ++mismatch[0]; seq_[i_] = data[j_]; }
            if (mismatch[threshold % 2] > threshold / 2) continue;

            count++;

            std::size_t pdr_id = index2pdr[seq_id];
            bool        is_left = index2lr[seq_id];

            // left primer on forward strand — store as pending
            if (is_left && index2strand[seq_id] % 2 == 0) {
                left_pending[pdr_id].push_back({start, seq_id, seq, seq_});
            }
            // right primer on reverse strand — check pairing
            else if (!is_left && index2strand[seq_id] % 2 == 1) {
                auto it = left_pending.find(pdr_id);
                if (it == left_pending.end()) continue;

                for (auto& lhit : it->second) {
                    std::size_t dist = start > lhit.pos ? 
                                       start - lhit.pos : 0;
                    if (dist < amp_min || dist > amp_max) continue;

                    // PAIRED — compute ΔG only now
                    double dg_l = Thal::compute_dimer_dg(lhit.seq, lhit.seq_);
                    double dg_r = Thal::compute_dimer_dg(seq, seq_);

                    if (dg_l < dg_thres)
                        local_results.emplace_back(Result{
                            index2label[lhit.seq_id],
                            lhit.seq, contig, lhit.seq_,
                            lhit.pos, 0, dg_l
                        });

                    if (dg_r < dg_thres)
                        local_results.emplace_back(Result{
                            index2label[seq_id],
                            seq, contig, seq_,
                            start, 1, dg_r
                        });
                }
            }
            // other strand combinations — skip entirely (no pairing possible)
        }
    }
    return count;
}

std::vector<Result> Automaton::process_fasta_chunks_paired(const std::string& path,
                                                    std::size_t overlap,
                                                    std::size_t threshold,
                                                    double dg_thres,
                                                    std::size_t amp_min,
                                                    std::size_t amp_max,
                                                    std::size_t chunk_size,
                                                    std::size_t io_block_size,
                                                    std::size_t nthreads) {
    
    TaskQueue q(std::max<std::size_t>(2, nthreads * 2));

    std::vector<std::thread> workers;
    workers.reserve(nthreads);

    std::vector<std::vector<Result>> buckets(nthreads);
    std::vector<std::size_t> counts(nthreads);

    // Thread-safe progress tracking
    std::atomic<std::size_t> processed_chunks{0};
    std::atomic<std::size_t> total_candidates{0};
    std::atomic<std::size_t> total_chars_processed{0};  // Add character counter

    for (std::size_t t = 0; t < nthreads; ++t) {
        workers.emplace_back([this, 
                              &q, 
                              &local_results = buckets[t], 
                              &local_counts = counts[t],
                              &threshold, 
                              &dg_thres,
                              &amp_min,
                              &amp_max,
                              &processed_chunks,
                              &total_candidates,
                              &total_chars_processed]() {  // Add to capture
            Task x;
            while (q.pop(x)) {
                std::size_t chunk_candidates = this->search_chunk_paired(
                    x.contig, 
                    x.seq.data(), 
                    x.seq.size(), 
                    x.start, 
                    threshold,
                    dg_thres,
                    amp_min,
                    amp_max,
                    local_results
                );
                local_counts += chunk_candidates;
                processed_chunks.fetch_add(1);
                total_candidates.fetch_add(chunk_candidates);
                total_chars_processed.fetch_add(x.seq.size());  // Track characters processed
            }
        });
    }

    std::thread prod([&, this]() {
        std::ifstream in(path, std::ios::binary);
        if (!in) {
            std::cerr << "Error: Cannot open " << path << std::endl;
            q.close();
            return;
        }

        // Progress tracking setup
        std::size_t contigs_processed = 0;
        std::size_t chunks_queued = 0;
        std::size_t chars_queued = 0;  // Track characters queued
        auto start_time = std::chrono::steady_clock::now();
        
        std::cout << "Searching " << path << " (threads: " << nthreads << ")" << std::endl;

        std::vector<char> io(io_block_size);
        std::string contig, seq_buf, header_acc;
        seq_buf.reserve(chunk_size + overlap);

        bool in_header = false;
        bool have_contig = false;
        std::size_t last_index = 0;

        auto flush_buf_as_last = [&]() {
            if (have_contig && !seq_buf.empty()) {
                q.push(Task{contig, seq_buf, last_index, true});
                chunks_queued++;
                chars_queued += seq_buf.size();  // Track characters
                seq_buf.clear();
            }
        };

        // Progress reporting with characters per second
        auto maybe_print_progress = [&]() {
            static bool header_printed = false;
            static auto last_print = start_time;
            
            auto now = std::chrono::steady_clock::now();
            auto elapsed_since_last = std::chrono::duration_cast<std::chrono::seconds>(now - last_print).count();
            
            if (elapsed_since_last >= 10) {
                if (!header_printed) {
                    std::cout << std::string(58, '-') << std::endl
                              << std::right 
                              << std::setw(10) << "Contigs" 
                              << std::setw(12) << "Chunks"
                              << std::setw(14) << "Candidates"  
                              << std::setw(10) << "Time(s)"
                              << std::setw(12) << "Rate(ch/s)" << std::endl;
                    std::cout << std::string(58, '-') << std::endl;
                    header_printed = true;
                }
                
                auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
                
                // Calculate characters per second rate
                double char_rate = 0.0;
                if (total_elapsed > 0) {
                    char_rate = static_cast<double>(total_chars_processed.load()) / total_elapsed;
                }
                
                // Format rate with appropriate units (K, M, G)
                std::string rate_str;
                if (char_rate >= 1e9) {
                    rate_str = std::to_string(char_rate / 1e9).substr(0, 4) + "G";
                } else if (char_rate >= 1e6) {
                    rate_str = std::to_string(char_rate / 1e6).substr(0, 4) + "M";
                } else if (char_rate >= 1e3) {
                    rate_str = std::to_string(char_rate / 1e3).substr(0, 4) + "K";
                } else {
                    rate_str = std::to_string(static_cast<int>(char_rate));
                }
                
                std::cout << std::right 
                          << std::setw(10) << contigs_processed
                          << std::setw(12) << chunks_queued  
                          << std::setw(14) << total_candidates.load()
                          << std::setw(10) << total_elapsed
                          << std::setw(12) << rate_str << std::endl;
                
                last_print = now;
            }
        };

        while (in) {
            in.read(io.data(), io.size());
            std::streamsize n = in.gcount();
            if (n <= 0) break;

            for (std::streamsize i = 0; i < n; i++) {
                char c = io[i];

                if (in_header) {
                    if (c == '\n' || c == '\r') {
                        size_t pos = header_acc.find_first_of(" \t");
                        contig = (pos != std::string::npos) ? 
                            header_acc.substr(0, pos) : header_acc;

                        if (have_contig) contigs_processed++;
                        
                        header_acc.clear();
                        in_header = false;
                        have_contig = true;
                        seq_buf.clear();
                        last_index = 0;
                    } else {
                        header_acc.push_back(c);
                    }
                    continue;
                }

                if (c == '>') {
                    flush_buf_as_last();
                    header_acc.clear();
                    in_header = true;
                    continue;
                }

                if (!have_contig || is_ws(c)) continue;

                c = std::toupper(static_cast<unsigned char>(c));
                seq_buf.push_back(c);

                if (seq_buf.size() >= chunk_size + overlap) {
                    std::string chunk(seq_buf.data(), chunk_size + overlap);
                    std::size_t chunk_chars = chunk.size();
                    
                    q.push(Task{contig, std::move(chunk), last_index, false});
                    chunks_queued++;
                    chars_queued += chunk_chars;

                    if (overlap > 0) {
                        seq_buf.erase(0, chunk_size);
                        last_index += chunk_size;
                    } else {
                        seq_buf.clear();
                        last_index += chunk_size;
                    }
                    
                    maybe_print_progress();
                }
            }
        }

        flush_buf_as_last();
        if (have_contig) contigs_processed++;
        
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - start_time).count();
        
        // Final rate calculation
        double final_char_rate = 0.0;
        if (elapsed > 0) {
            final_char_rate = static_cast<double>(total_chars_processed.load()) / elapsed;
        }
        
        std::cout << std::string(58, '-') << std::endl;
        std::cout << "Completed: " << contigs_processed << " contigs, " 
                 << chunks_queued << " chunks, " 
                 << (total_chars_processed.load() / 1e6) << "M chars (" 
                 << elapsed << "s, " 
                 << (final_char_rate / 1e6) << "M chars/s)" << std::endl;
        
        q.close();
    });

    prod.join();
    for (auto& th : workers) th.join();

    std::size_t total = 0;
    for (auto& count : counts) 
        total += count;
    
    std::cout << "Found " << total << " candidates in reference" << std::endl;

    std::vector<Result> results;
    for (auto& b : buckets) {
        results.insert(results.end(),
                    std::make_move_iterator(b.begin()),
                    std::make_move_iterator(b.end()));
    }
    return results;
}


std::vector<Result> Automaton::process_fasta_chunks(const std::string& path,
                                                    std::size_t overlap,
                                                    std::size_t threshold,
                                                    double dg_thres,
                                                    std::size_t chunk_size,
                                                    std::size_t io_block_size,
                                                    std::size_t nthreads) {
    
    TaskQueue q(std::max<std::size_t>(2, nthreads * 2));

    std::vector<std::thread> workers;
    workers.reserve(nthreads);

    std::vector<std::vector<Result>> buckets(nthreads);
    std::vector<std::size_t> counts(nthreads);

    // Thread-safe progress tracking
    std::atomic<std::size_t> processed_chunks{0};
    std::atomic<std::size_t> total_candidates{0};
    std::atomic<std::size_t> total_chars_processed{0};  // Add character counter

    for (std::size_t t = 0; t < nthreads; ++t) {
        workers.emplace_back([this, 
                              &q, 
                              &local_results = buckets[t], 
                              &local_counts = counts[t],
                              &threshold, 
                              &dg_thres,
                              &processed_chunks,
                              &total_candidates,
                              &total_chars_processed]() {  // Add to capture
            Task x;
            while (q.pop(x)) {
                std::size_t chunk_candidates = this->search_chunk(
                    x.contig, 
                    x.seq.data(), 
                    x.seq.size(), 
                    x.start, 
                    threshold,
                    dg_thres,
                    local_results
                );
                local_counts += chunk_candidates;
                processed_chunks.fetch_add(1);
                total_candidates.fetch_add(chunk_candidates);
                total_chars_processed.fetch_add(x.seq.size());  // Track characters processed
            }
        });
    }

    std::thread prod([&, this]() {
        std::ifstream in(path, std::ios::binary);
        if (!in) {
            std::cerr << "Error: Cannot open " << path << std::endl;
            q.close();
            return;
        }

        // Progress tracking setup
        std::size_t contigs_processed = 0;
        std::size_t chunks_queued = 0;
        std::size_t chars_queued = 0;  // Track characters queued
        auto start_time = std::chrono::steady_clock::now();
        
        std::cout << "Searching " << path << " (threads: " << nthreads << ")" << std::endl;

        std::vector<char> io(io_block_size);
        std::string contig, seq_buf, header_acc;
        seq_buf.reserve(chunk_size + overlap);

        bool in_header = false;
        bool have_contig = false;
        std::size_t last_index = 0;

        auto flush_buf_as_last = [&]() {
            if (have_contig && !seq_buf.empty()) {
                q.push(Task{contig, seq_buf, last_index, true});
                chunks_queued++;
                chars_queued += seq_buf.size();  // Track characters
                seq_buf.clear();
            }
        };

        // Progress reporting with characters per second
        auto maybe_print_progress = [&]() {
            static bool header_printed = false;
            static auto last_print = start_time;
            
            auto now = std::chrono::steady_clock::now();
            auto elapsed_since_last = std::chrono::duration_cast<std::chrono::seconds>(now - last_print).count();
            
            if (elapsed_since_last >= 10) {
                if (!header_printed) {
                    std::cout << std::string(58, '-') << std::endl
                              << std::right 
                              << std::setw(10) << "Contigs" 
                              << std::setw(12) << "Chunks"
                              << std::setw(14) << "Candidates"  
                              << std::setw(10) << "Time(s)"
                              << std::setw(12) << "Rate(ch/s)" << std::endl;
                    std::cout << std::string(58, '-') << std::endl;
                    header_printed = true;
                }
                
                auto total_elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start_time).count();
                
                // Calculate characters per second rate
                double char_rate = 0.0;
                if (total_elapsed > 0) {
                    char_rate = static_cast<double>(total_chars_processed.load()) / total_elapsed;
                }
                
                // Format rate with appropriate units (K, M, G)
                std::string rate_str;
                if (char_rate >= 1e9) {
                    rate_str = std::to_string(char_rate / 1e9).substr(0, 4) + "G";
                } else if (char_rate >= 1e6) {
                    rate_str = std::to_string(char_rate / 1e6).substr(0, 4) + "M";
                } else if (char_rate >= 1e3) {
                    rate_str = std::to_string(char_rate / 1e3).substr(0, 4) + "K";
                } else {
                    rate_str = std::to_string(static_cast<int>(char_rate));
                }
                
                std::cout << std::right 
                          << std::setw(10) << contigs_processed
                          << std::setw(12) << chunks_queued  
                          << std::setw(14) << total_candidates.load()
                          << std::setw(10) << total_elapsed
                          << std::setw(12) << rate_str << std::endl;
                
                last_print = now;
            }
        };

        while (in) {
            in.read(io.data(), io.size());
            std::streamsize n = in.gcount();
            if (n <= 0) break;

            for (std::streamsize i = 0; i < n; i++) {
                char c = io[i];

                if (in_header) {
                    if (c == '\n' || c == '\r') {
                        size_t pos = header_acc.find_first_of(" \t");
                        contig = (pos != std::string::npos) ? 
                            header_acc.substr(0, pos) : header_acc;

                        if (have_contig) contigs_processed++;
                        
                        header_acc.clear();
                        in_header = false;
                        have_contig = true;
                        seq_buf.clear();
                        last_index = 0;
                    } else {
                        header_acc.push_back(c);
                    }
                    continue;
                }

                if (c == '>') {
                    flush_buf_as_last();
                    header_acc.clear();
                    in_header = true;
                    continue;
                }

                if (!have_contig || is_ws(c)) continue;

                c = std::toupper(static_cast<unsigned char>(c));
                seq_buf.push_back(c);

                if (seq_buf.size() >= chunk_size + overlap) {
                    std::string chunk(seq_buf.data(), chunk_size + overlap);
                    std::size_t chunk_chars = chunk.size();
                    
                    q.push(Task{contig, std::move(chunk), last_index, false});
                    chunks_queued++;
                    chars_queued += chunk_chars;

                    if (overlap > 0) {
                        seq_buf.erase(0, chunk_size);
                        last_index += chunk_size;
                    } else {
                        seq_buf.clear();
                        last_index += chunk_size;
                    }
                    
                    maybe_print_progress();
                }
            }
        }

        flush_buf_as_last();
        if (have_contig) contigs_processed++;
        
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            std::chrono::steady_clock::now() - start_time).count();
        
        // Final rate calculation
        double final_char_rate = 0.0;
        if (elapsed > 0) {
            final_char_rate = static_cast<double>(total_chars_processed.load()) / elapsed;
        }
        
        std::cout << std::string(58, '-') << std::endl;
        std::cout << "Completed: " << contigs_processed << " contigs, " 
                 << chunks_queued << " chunks, " 
                 << (total_chars_processed.load() / 1e6) << "M chars (" 
                 << elapsed << "s, " 
                 << (final_char_rate / 1e6) << "M chars/s)" << std::endl;
        
        q.close();
    });

    prod.join();
    for (auto& th : workers) th.join();

    std::size_t total = 0;
    for (auto& count : counts) 
        total += count;
    
    std::cout << "Found " << total << " candidates in reference" << std::endl;

    std::vector<Result> results;
    for (auto& b : buckets) {
        results.insert(results.end(),
                    std::make_move_iterator(b.begin()),
                    std::make_move_iterator(b.end()));
    }
    return results;
}
