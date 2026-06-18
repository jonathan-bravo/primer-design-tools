#include "utility.hpp"
#include "automaton.hpp"
#include "risk_optimizer.hpp"
#include "graph.hpp"
#include "pipeline.hpp"

static void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " -i <input> -o <output> [options]\n"
        << "\n"
        /*
        << "Modes (required):\n"
        << "  pdr               Run PDR optimization\n"
        << "  dimer             Run dimer optimization\n"
        << "  off               Run off-target search\n"
        << "\n"
        */
        << "I/O (required):\n"
        << "  -i <file>         Input file\n"
        << "  -o <file>         Output file\n"
        << "  -r <file>         Reference file (for off target search)\n"
        << "\n"
        << "PDR / risk optimization parameters:\n"
        << "  -Ln <int>         Amplicon max length (default 420)\n"
        << "  -Lx <int>         Amplicon min length (default 252)\n"
        << "  -Lp <int>         PDR length (default 40)\n"
        << "  -Ux <double>      Max PDR score (default 10000)\n"
        << "  -Un <double>      Min PDR score (default 0.1)\n"
        << "\n"
        << "Primer3 parameters:\n"
        << "  --p3-opt-size <int>     Primer optimal size (default 20)\n"
        << "  --p3-min-size <int>     Primer min size (default 18)\n"
        << "  --p3-max-size <int>     Primer max size (default 25)\n"
        << "  --p3-opt-tm <double>    Primer optimal Tm (default 60.0)\n"
        << "  --p3-min-tm <double>    Primer min Tm (default 57.0)\n"
        << "  --p3-max-tm <double>    Primer max Tm (default 63.0)\n"
        << "  --p3-min-gc <double>    Primer min GC% (default 40.0)\n"
        << "  --p3-max-gc <double>    Primer max GC% (default 60.0)\n"
        << "  --p3-num-return <int>   Number of candidates (default 10)\n"
        << "\n"
        << "Off-target search parameters:\n"
        << "  -k <int>          K-mer length (default 8)\n"
        << "  -H <int>          Hamming threshold (default 2)\n"
        << "  -G <double>       Delta G threshold (default 1e8)\n"
        << "  -t <int>          Number of threads (default 8)\n"
        << "  -C <int>          Chunk size (default 100000000)\n"
        << "  -B <int>          Block size (default 134217728)\n"
        << "\n"
        << "Dimer parameters:\n"
        << "  -I <int>          Iterations for solver (default 1000)\n"
        << "  -S <int>          Random seed (default 42)\n"
        << "\n"
        << "Thermodynamic parameters:\n"
        << "  --mv <double>     Monovalent conc (default 50.0)\n"
        << "  --dv <double>     Divalent conc (default 1.5)\n"
        << "  --dntp <double>   dNTP conc (default 0.6)\n"
        << "  --dna <double>    DNA conc (default 200.0)\n"
        << "  --temp <double>   Temperature in C (default 37.0)\n"
        << "Other:\n"
        << "  -h, --help        Show this help message\n"
        << "\n";
}

static bool is_flag(const char* s, 
                    const char* shortf, 
                    const char* longf = nullptr) {
    if (!s) return false;
    if (shortf && std::strcmp(s, shortf) == 0) return true;
    if (longf && std::strcmp(s, longf) == 0) return true;
    return false;
}

static const char* require_value(int& i, 
                                 int argc, 
                                 char** argv, 
                                 const char* opt) {
    if (i + 1 >= argc) {
        throw std::runtime_error(
            std::string("Missing value for option: ") + opt
        );
    }
    return argv[++i];
}

static std::size_t require_integer(int& i, 
                                   int argc, 
                                   char** argv, 
                                   const char* opt) {
    const char* v = require_value(i, argc, argv, opt);
    try {
        std::size_t pos = 0;
        std::size_t tmp = std::stoull(v, &pos);
        if (pos != std::strlen(v)) {
            throw std::runtime_error("Seed contains non-numeric characters");
        }
        return tmp;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("Invalid number for ") + opt + ": " + 
            v + " (" + e.what() + ")"
        );
    }
}

static double require_numeric(int& i, 
                              int argc, 
                              char** argv, 
                              const char* opt) {
    const char* v = require_value(i, argc, argv, opt);
    try {
        std::size_t pos = 0;
        double tmp = std::stod(v, &pos);
        if (pos != std::strlen(v)) {
            throw std::runtime_error("contains non-numeric characters");
        }
        return tmp;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("Invalid number for ") + opt + ": " + 
            v + " (" + e.what() + ")"
        );
    }
}

static Args parse_args_file(const std::string& path, Args a) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open args file: " + path);

    std::string line;
    while (std::getline(f, line)) {
        // strip comments and blank lines
        auto pos = line.find('#');
        if (pos != std::string::npos) line = line.substr(0, pos);
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

        std::istringstream ss(line);
        std::string key, eq, val;
        if (!(ss >> key >> eq >> val) || eq != "=")
            throw std::runtime_error("Malformed args file line: " + line);

        if      (key == "input_file")   a.input_file   = val;
        else if (key == "output_file")  a.output_file  = val;
        else if (key == "ref_file")     a.ref_file     = val;
        else if (key == "seed")         a.seed         = std::stoull(val);
        else if (key == "len_amp")      a.len_amp      = std::stoull(val);
        else if (key == "len_amp_min")  a.len_amp_min  = std::stoull(val);
        else if (key == "len_PDR")      a.len_PDR      = std::stoull(val);
        else if (key == "u_max")        a.u_max        = std::stod(val);
        else if (key == "u_min")        a.u_min        = std::stod(val);
        else if (key == "iter")         a.iter         = std::stoull(val);
        else if (key == "kmer_len")     a.kmer_len     = std::stoull(val);
        else if (key == "threshold")    a.threshold    = std::stoull(val) * 2;
        else if (key == "nthreads")     a.nthreads     = std::stoull(val);
        else if (key == "block_size")   a.block_size   = std::stoull(val);
        else if (key == "chunk_size")   a.chunk_size   = std::stoull(val);
        else if (key == "dg_thres")     a.dg_thres     = std::stod(val);
        else if (key == "mv")           a.mv           = std::stod(val);
        else if (key == "dv")           a.dv           = std::stod(val);
        else if (key == "dntp")         a.dntp         = std::stod(val);
        else if (key == "dna_conc")     a.dna_conc     = std::stod(val);
        else if (key == "temp")         a.temp         = std::stod(val);
        else if (key == "p3_opt_size")  a.p3_opt_size  = std::stoull(val);
        else if (key == "p3_min_size")  a.p3_min_size  = std::stoull(val);
        else if (key == "p3_max_size")  a.p3_max_size  = std::stoull(val);
        else if (key == "p3_opt_tm")    a.p3_opt_tm    = std::stod(val);
        else if (key == "p3_min_tm")    a.p3_min_tm    = std::stod(val);
        else if (key == "p3_max_tm")    a.p3_max_tm    = std::stod(val);
        else if (key == "p3_min_gc")    a.p3_min_gc    = std::stod(val);
        else if (key == "p3_max_gc")    a.p3_max_gc    = std::stod(val);
        else if (key == "num_return")   a.num_return   = std::stoull(val);
        else throw std::runtime_error("Unknown key in args file: " + key);
    }
    return a;
}

static void save_args(const Args& a, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open args save file: " + path);

    auto now = std::chrono::system_clock::now();
    auto t   = std::chrono::system_clock::to_time_t(now);
    f << "# Saved args -- " << std::ctime(&t);

    f << "\n# I/O\n";
    f << "input_file  = " << a.input_file  << "\n";
    f << "output_file = " << a.output_file << "\n";
    if (!a.ref_file.empty())
    f << "ref_file    = " << a.ref_file    << "\n";

    f << "\n# PDR\n";
    f << "len_amp     = " << a.len_amp     << "\n";
    f << "len_amp_min = " << a.len_amp_min << "\n";
    f << "len_PDR     = " << a.len_PDR     << "\n";
    f << "u_max       = " << a.u_max       << "\n";
    f << "u_min       = " << a.u_min       << "\n";
    f << "seed        = " << a.seed        << "\n";

    f << "\n# Dimer\n";
    f << "iter        = " << a.iter        << "\n";

    f << "\n# Off-target search\n";
    f << "kmer_len    = " << a.kmer_len    << "\n";
    f << "threshold   = " << a.threshold/2 << "\n";
    f << "nthreads    = " << a.nthreads    << "\n";
    f << "chunk_size  = " << a.chunk_size  << "\n";
    f << "block_size  = " << a.block_size  << "\n";

    f << "\n# Thermodynamics\n";
    f << "dg_thres    = " << a.dg_thres    << "\n";
    f << "mv          = " << a.mv          << "\n";
    f << "dv          = " << a.dv          << "\n";
    f << "dntp        = " << a.dntp        << "\n";
    f << "dna_conc    = " << a.dna_conc    << "\n";
    f << "temp        = " << a.temp        << "\n";

    f << "\n# Primer3\n";
    f << "p3_opt_size = " << a.p3_opt_size << "\n";
    f << "p3_min_size = " << a.p3_min_size << "\n";
    f << "p3_max_size = " << a.p3_max_size << "\n";
    f << "p3_opt_tm   = " << a.p3_opt_tm   << "\n";
    f << "p3_min_tm   = " << a.p3_min_tm   << "\n";
    f << "p3_max_tm   = " << a.p3_max_tm   << "\n";
    f << "p3_min_gc   = " << a.p3_min_gc   << "\n";
    f << "p3_max_gc   = " << a.p3_max_gc   << "\n";
    f << "num_return  = " << a.num_return  << "\n";
}

static Args parse_args(int argc, char** argv) {
    Args a;
 
    // Print command line start
    std::cout << "\nCommand executed:" << std::endl;
    std::cout << argv[0] << " " << argv[1];
 
    for (int i = 1; i < argc; ++i) {
        const char* cur = argv[i];
 
        // Print current argument with formatting
        std::cout << " \\" << std::endl << "  " << cur;
 
        if (is_flag(cur, "-h", "--help")) {
            a.help = true;
            std::cout << std::endl << std::endl;
            return a;
        } else if (is_flag(cur, "--args")) {
            const char* path = require_value(i, argc, argv, "--args");
            a = parse_args_file(path, a);
            std::cout << " " << path << "  # Args file";
        } else if (is_flag(cur, "-i")) {
            a.input_file = require_value(i, argc, argv, "-i");
            std::cout << " " << a.input_file << "  # Input file";
        } else if (is_flag(cur, "-o")) {
            a.output_file = require_value(i, argc, argv, "-o");
            std::cout << " " << a.output_file << "  # Output file";
        } else if (is_flag(cur, "-r")) {
            a.ref_file = require_value(i, argc, argv, "-r");
            std::cout << " " << a.ref_file << "  # Reference file";
        } else if (is_flag(cur, "-k")) {
            a.kmer_len = require_integer(i, argc, argv, "-k");
            std::cout << " " << a.kmer_len << "  # Kmer length";
        } else if (is_flag(cur, "-H")) {
            a.threshold = require_integer(i, argc, argv, "-H") * 2;
            std::cout << " " << (a.threshold/2) << "  # Hamming threshold";
        } else if (is_flag(cur, "-H2")) {
            a.threshold = require_integer(i, argc, argv, "-H2") * 2 + 1;
            std::cout << " " << ((a.threshold-1)/2) << "  # Hamming threshold (odd)";
        } else if (is_flag(cur, "-t")) {
            a.nthreads = require_integer(i, argc, argv, "-t");
            std::cout << " " << a.nthreads << "  # Threads";
        } else if (is_flag(cur, "-C")) {
            a.chunk_size = require_integer(i, argc, argv, "-C");
            std::cout << " " << a.chunk_size << "  # Chunk size";
        } else if (is_flag(cur, "-B")) {
            a.block_size = require_integer(i, argc, argv, "-B");
            std::cout << " " << a.block_size << "  # Block size";
        } else if (is_flag(cur, "-G")) {
            a.dg_thres = require_numeric(i, argc, argv, "-G");
            std::cout << " " << a.dg_thres << "  # dG threshold";
        } else if (is_flag(cur, "--mv")) {
            a.mv = require_numeric(i, argc, argv, "--mv");
            std::cout << " " << a.mv << "  # Monovalent concentration";
        } else if (is_flag(cur, "--dv")) {
            a.dv = require_numeric(i, argc, argv, "--dv");
            std::cout << " " << a.dv << "  # Divalent concentration";
        } else if (is_flag(cur, "--dntp")) {
            a.dntp = require_numeric(i, argc, argv, "--dntp");
            std::cout << " " << a.dntp << "  # dNTP concentration";
        } else if (is_flag(cur, "--dna")) {
            a.dna_conc = require_numeric(i, argc, argv, "--dna");
            std::cout << " " << a.dna_conc << "  # DNA concentration";
        } else if (is_flag(cur, "--temp")) {
            a.temp = require_numeric(i, argc, argv, "--temp");
            std::cout << " " << a.temp << "  # Temperature";
        } else if (is_flag(cur, "-S")) {
            a.seed = require_integer(i, argc, argv, "-S");
            std::cout << " " << a.seed << "  # Random seed";
        } else if (is_flag(cur, "-Ln")) {
            a.len_amp = require_integer(i, argc, argv, "-Ln");
            std::cout << " " << a.len_amp << "  # Amplicon length max";
        } else if (is_flag(cur, "-Lx")) {
            a.len_amp_min = require_integer(i, argc, argv, "-Lx");
            std::cout << " " << a.len_amp_min << "  # Amplicon length min";
        } else if (is_flag(cur, "-Lp")) {
            a.len_PDR = require_integer(i, argc, argv, "-Lp");
            std::cout << " " << a.len_PDR << "  # PDR length";
        } else if (is_flag(cur, "-Ux")) {
            a.u_max = require_numeric(i, argc, argv, "-Ux");
            std::cout << " " << a.u_max << "  # Max PDR score";
        } else if (is_flag(cur, "-Un")) {
            a.u_min = require_numeric(i, argc, argv, "-Un");
            std::cout << " " << a.u_min << "  # Min PDR score";
        } else if (is_flag(cur, "-I")) {
            a.iter = require_integer(i, argc, argv, "-I");
            std::cout << " " << a.iter << "  # Iterations";
        }  else if (is_flag(cur, "--p3-opt-size")) {
            a.p3_opt_size = require_integer(i, argc, argv, "--p3-opt-size");
            std::cout << " " << a.p3_opt_size << "  # Primer optimal size";
        } else if (is_flag(cur, "--p3-min-size")) {
            a.p3_min_size = require_integer(i, argc, argv, "--p3-min-size");
            std::cout << " " << a.p3_min_size << "  # Primer min size";
        } else if (is_flag(cur, "--p3-max-size")) {
            a.p3_max_size = require_integer(i, argc, argv, "--p3-max-size");
            std::cout << " " << a.p3_max_size << "  # Primer max size";
        } else if (is_flag(cur, "--p3-opt-tm")) {
            a.p3_opt_tm = require_numeric(i, argc, argv, "--p3-opt-tm");
            std::cout << " " << a.p3_opt_tm << "  # Primer optimal Tm";
        } else if (is_flag(cur, "--p3-min-tm")) {
            a.p3_min_tm = require_numeric(i, argc, argv, "--p3-min-tm");
            std::cout << " " << a.p3_min_tm << "  # Primer min Tm";
        } else if (is_flag(cur, "--p3-max-tm")) {
            a.p3_max_tm = require_numeric(i, argc, argv, "--p3-max-tm");
            std::cout << " " << a.p3_max_tm << "  # Primer max Tm";
        } else if (is_flag(cur, "--p3-min-gc")) {
            a.p3_min_gc = require_numeric(i, argc, argv, "--p3-min-gc");
            std::cout << " " << a.p3_min_gc << "  # Primer min GC%";
        } else if (is_flag(cur, "--p3-max-gc")) {
            a.p3_max_gc = require_numeric(i, argc, argv, "--p3-max-gc");
            std::cout << " " << a.p3_max_gc << "  # Primer max GC%";
        } else if (is_flag(cur, "--p3-num-return")) {
            a.num_return = require_integer(i, argc, argv, "--p3-num-return");
            std::cout << " " << a.num_return << "  # Primers to return";
        } else {
            throw std::runtime_error(std::string("Unknown option: ") + cur);
        }
    }
 
    // End command line display
    std::cout << std::endl << std::endl;

    if (a.input_file.empty() || a.output_file.empty()) {
        throw std::runtime_error("Missing required -i / -o");
    }
    /*
    if (a.mode == "off" && a.ref_file.empty()) {
        throw std::runtime_error("Missing required -r");
    }
    */
    std::cout << "\nProgram Parameters\n";
    std::cout << "  Input file      : " << a.input_file  << "\n";
    std::cout << "  Output file     : " << a.output_file << "\n";

    if (!a.ref_file.empty())
        std::cout << "  Reference file  : " << a.ref_file    << "\n";

    std::cout << "\nPDR Parameters\n";
    std::cout << "  Amp length max  : " << a.len_amp     << "\n";
    std::cout << "  Amp length min  : " << a.len_amp_min << "\n";
    std::cout << "  PDR length      : " << a.len_PDR     << "\n";
    std::cout << "  Max PDR score   : " << a.u_max       << "\n";
    std::cout << "  Min PDR score   : " << a.u_min       << "\n";
    std::cout << "  Random seed     : " << a.seed        << "\n";

    std::cout << "\nPrimer3 Parameters\n";
    std::cout << "  Opt size        : " << a.p3_opt_size << "\n";
    std::cout << "  Min size        : " << a.p3_min_size << "\n";
    std::cout << "  Max size        : " << a.p3_max_size << "\n";
    std::cout << "  Opt Tm          : " << a.p3_opt_tm   << "\n";
    std::cout << "  Min Tm          : " << a.p3_min_tm   << "\n";
    std::cout << "  Max Tm          : " << a.p3_max_tm   << "\n";
    std::cout << "  Min GC%         : " << a.p3_min_gc   << "\n";
    std::cout << "  Max GC%         : " << a.p3_max_gc   << "\n";
    std::cout << "  Num return      : " << a.num_return  << "\n";

    std::cout << "\nOff-target Search Parameters\n";
    std::cout << "  Kmer length     : " << a.kmer_len    << "\n";
    std::cout << "  Hamming thres   : " << a.threshold/2 << "\n";
    std::cout << "  Threads         : " << a.nthreads    << "\n";
    std::cout << "  Chunk size      : " << a.chunk_size  << "\n";
    std::cout << "  Block size      : " << a.block_size  << "\n";

    std::cout << "\nDimer Parameters\n";
    std::cout << "  Iterations      : " << a.iter        << "\n";
    std::cout << "  Random seed     : " << a.seed        << "\n";

    std::cout << "\nThermo parameters\n";
    std::cout << "  dG thres        : " << a.dg_thres    << "\n";
    std::cout << "  Monovalent (mv) : " << a.mv          << "\n";
    std::cout << "  Divalent (dv)   : " << a.dv          << "\n";
    std::cout << "  dNTP            : " << a.dntp        << "\n";
    std::cout << "  DNA conc        : " << a.dna_conc    << "\n";
    std::cout << "  Temp (C)        : " << a.temp        << "\n";

    std::cout << std::endl;

    std::string save_path = a.output_file + ".args";
    save_args(a, save_path);
    std::cout << "Args saved to: " << save_path << "\n\n";
    return a;
}

int main(int argc, char** argv) {
    try {
        Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config", 50.0, 1.5, 0.6, 250.0, 60.0);
        Thal::run_tests();          // default Ta=60, Te=72
        Thal::run_tests(58.0, 72.0); // custom temperatures
        PipelineContext ctx;

        ctx.args = parse_args(argc, argv);
        
        if (ctx.args.help) {
            print_usage(argv[0]);
            return 0;
        }

        std::vector<std::string> fasta_files;

        if (fs::is_directory(ctx.args.input_file)) {
            for (auto& entry : fs::directory_iterator(ctx.args.input_file)) {
                auto ext = entry.path().extension().string();
                if (ext == ".fasta" || ext == ".fa") {
                    fasta_files.push_back(entry.path().string());
                }
            }
            std::sort(fasta_files.begin(), fasta_files.end());
        } else {
            // backward-compatible: single FASTA file
            fasta_files.push_back(ctx.args.input_file);
        }

        ctx.segments.resize(fasta_files.size());

        for (std::size_t i = 0; i < fasta_files.size(); ++i) {
            auto& seg = ctx.segments[i];
            seg.name = fs::path(fasta_files[i]).stem().string();

            std::size_t n = read_fasta(fasta_files[i], seg.labels, seg.sequences);
            if (n == 0) {
                throw std::runtime_error("No sequences loaded from: " + fasta_files[i]);
            }

            seg.tmpl = msa_consensus(seg.sequences);
            char_freq(seg.tmpl);
            write_fasta(ctx.args.output_file + "." + seg.name + ".tmpl.fasta",
                        {"template"}, {seg.tmpl}, 80);

            std::cout << "  [" << (i+1) << "/" << fasta_files.size() << "] "
                    << seg.name << " — " << n << " sequences, "
                    << seg.tmpl.size() << " bp\n";
        }

        Pipeline p;
        p.add(std::make_unique<PDRStage>());
        p.add(std::make_unique<PrimerSelStage>());
        p.add(std::make_unique<OffTargetStage>());
        p.add(std::make_unique<DimerStage>());
        p.run(ctx);

        print_solution(ctx.solution_primers);

        save_solution(ctx.segments, ctx.solution_primers, ctx.args.output_file + ".bed");

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n";
        print_usage(argv[0]);
        return 2;
    }
}
