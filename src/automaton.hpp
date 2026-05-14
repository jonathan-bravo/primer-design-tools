#pragma once
#include "utility.hpp"
#include "task.hpp"
#include "thal.hpp"

struct CandidateKey {
    std::size_t seq_id;
    std::size_t start;

    bool operator==(const CandidateKey& other) const {
        return seq_id == other.seq_id && start == other.start;
    }
};

struct CandidateKeyHash {
    std::size_t operator()(const CandidateKey& x) const {
        return std::hash<std::size_t>{}(x.seq_id) ^
               (std::hash<std::size_t>{}(x.start) << 1);
    }
};

struct ActiveEntry {
    CandidateKey key;
    std::size_t expire_at;
};

class AutomatonNode {
public:
    AutomatonNode();
    ~AutomatonNode();

private:
    friend class Automaton;

    AutomatonNode *children[ALPHABET_SIZE], *prev, *parent;
    std::vector<std::size_t> index;
    std::size_t ch;
};

class Automaton {
public:
    Automaton(const std::vector<std::string> &labels, 
              const std::vector<std::string> &seqs, 
              std::size_t k);
    ~Automaton();
    std::vector<Result> search(const std::string &path, 
                               std::size_t overlap, 
                               std::size_t threshold,
                               double dg_thres,
                               std::size_t chunk_size, 
                               std::size_t io_block_size, 
                               std::size_t nthreads);
    std::vector<Result> search_paired(const std::string& path,
                                        std::size_t overlap,
                                        std::size_t threshold,
                                        double dg_thres,
                                        std::size_t amp_min,
                                        std::size_t amp_max,
                                        std::size_t chunk_size,
                                        std::size_t io_block_size,
                                        std::size_t nthreads);

private:
    std::size_t search_chunk(const std::string& contig,
                             const char* data,
                             std::size_t len,
                             std::size_t index,
                             std::size_t threshold,
                             double dg_thres,
                             std::vector<Result> &local_results);
    std::vector<Result> process_fasta_chunks(const std::string &path, 
                                             std::size_t overlap, 
                                             std::size_t threshold,
                                             double dg_thres,
                                             std::size_t chunk_size, 
                                             std::size_t io_block_size,
                                             std::size_t nthreads);
    std::size_t search_chunk_paired(const std::string& contig,
                                                const char* data,
                                                std::size_t len,
                                                std::size_t index,
                                                std::size_t threshold,
                                                double dg_thres,
                                                std::size_t amp_min,
                                                std::size_t amp_max,
                                                std::vector<Result>& local_results);
    std::vector<Result> process_fasta_chunks_paired(const std::string &path, 
                                             std::size_t overlap, 
                                             std::size_t threshold,
                                             double dg_thres,
                                             std::size_t amp_min,
                                             std::size_t amp_max,
                                             std::size_t chunk_size, 
                                             std::size_t io_block_size,
                                             std::size_t nthreads);
private:
    AutomatonNode *root;
    std::vector<std::string> index2label, index2data;
    std::vector<std::size_t> index2pos, index2strand;
    std::unordered_map<std::string, std::size_t> label2index;
    std::vector<std::size_t> index2pdr;
    std::vector<bool>        index2lr;
    std::size_t k;  
};
