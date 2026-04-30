
# DPRO - Primer Design and Off-Target Analysis Tool

DPRO is a comprehensive C++ tool for multiplex PCR primer design optimization and tiling that implements a four-stage pipeline: PDR (Primer Design Region) optimization to identify high-quality, conserved regions with sufficient primer candidates; primer selection using the established Primer3 algorithm to generate high-quality primer pairs with proper thermodynamic constraints; off-target screening against reference databases using k-mer based searching to filter primers with potential off-target binding; and dimer minimization using multiple optimization algorithms to select primer sets with minimal primer-dimer formation potential, ensuring both specificity and compatibility for successful multiplex PCR amplification and tiling coverage.

## Installation

### Prerequisites
- **C++17 compatible compiler** (GCC 7+, Clang 5+)
- **Primer3 library** (https://github.com/primer3-org/primer3)

### Build
```bash
git clone <repo_url>
cd primer-design-tools
make PRIMER3=/path/to/primer3
```
The executable will be generated as `bin/dpro`.

## Usage

```bash
./bin/dpro -i <input> -o <output> [options]
```

### Core Arguments
```bash
-i <file>    Input sequence file (required)
-o <file>    Output directory (required)
-r <file>    Reference genome file (optional)
-h, --help   Show help message
```

## Parameters

### Search Parameters
```bash
-k <int>       Kmer length for indexing (default: 15)
-H <int>       Hamming distance threshold (default: 3)  
-H2 <int>      Alternative Hamming threshold (default: 3)
-t <int>       Number of threads (default: 8)
-C <int>       Processing chunk size (default: 100000)
-B <int>       I/O block size (default: 8192)
```

### Optimization Parameters
```bash
-Ln <int>     Maximum amplicon length (default: 420)
-Lx <int>     Minimum amplicon length (default: 252)  
-Lp <int>     Primer design region length (default: 40)
-Ux <double>  Maximum risk threshold (default: 10000.0)
-Un <double>  Minimum risk threshold (default: 0.1)
-I <int>      Optimization iterations (default: 1000)
-S <int>      Random seed (default: 42)
```

### Thermodynamic Parameters  
```bash
-G <double>     dG threshold for filtering (default: -20000.0)
--mv <double>   Monovalent concentration mM (default: 50.0)
--dv <double>   Divalent concentration mM (default: 1.5) 
--dntp <double> dNTP concentration mM (default: 0.6)
--dna <double>  DNA concentration nM (default: 50.0)
--temp <double> Temperature °C (default: 60.0)
```

## Examples

### Basic Pipeline Analysis
```bash
# Complete pipeline with reference genome screening
./bin/dpro \
  -i sequences.fasta \
  -o results/ \
  -r reference_genome.fna

# Pipeline without off-target screening
./bin/dpro \
  -i sequences.fasta \
  -o results/

# Custom parameters for stringent filtering
./bin/dpro \
  -i sequences.fasta \
  -o results/ \
  -r genome.fna \
  -G -25000.0 \
  -t 16 \
  -H 2
```

### Advanced Configuration
```bash
# High-throughput analysis
./bin/dpro all \
  -i large_dataset.fasta \
  -o analysis_results/ \
  -r combined_genomes.fna \
  -t 32 \
  -C 200000 \
  -B 16384

# Custom thermodynamic conditions
./bin/dpro all \
  -i primers.fasta \
  -o filtered_results/ \
  -r host_genome.fna \
  --temp 65.0 \
  --mv 100.0 \
  --dna 25.0
```

## Input Format

- **Input**: FASTA format with target sequences
- **Reference**: FASTA format reference genomes (optional)

## Output Structure
### Pipeline Stages
1. **PDR (Primer Design Region) Optimization** - Identifies high-quality, conserved regions with sufficient primer candidates and achieves optimal tiling coverage through constrained optimization
2. **Primer Selection** - Uses Primer3 algorithm to generate high-quality primer pairs with proper thermodynamic constraints within optimized regions
3. **Off-target Screening** - Reference genome analysis using k-mer based searching with thermodynamic filtering to remove primers with potential off-target binding
4. **Dimer Minimization** - Employs multiple optimization algorithms (Random Search, Simulated Annealing, Tabu Search, Genetic Algorithm) to select primer sets with minimal primer-dimer formation potential

### Output Structure

The pipeline generates comprehensive analysis reports:

### PDR Optimization Results
```
PDRs: 20
uncovered: 3 (front), 308 (rear)
loss: 109.869

Search for optimal PDR...
--------------------------------------------
    iter          u'       loss'        loss
--------------------------------------------
       1        5000        9000     257.973
       2        2500        4500     257.973
       ......   ......      ......    ......
      17     52.5665     109.869     109.869
--------------------------------------------
Search completed!
```

### Primer Selection Summary
```
Starting primer selection for 10 PDRs
Product size: [252, 460] bp  |  Template length: 3063 bp
------------------------------------------------------------
     Index                PDRs      Left     Right     Pairs
------------------------------------------------------------
         0            [3, 383]        46        73         5
         1          [252, 632]        80        19         5
         2          [585, 965]        47        87         5
         3         [800, 1180]        50        64         5
         4        [1122, 1502]        49        20         5
    ......              ......    ......    ......    ......
------------------------------------------------------------
```

### Off-target Screening Results
```
Searching ./scripts/ref/livestock_combined.fna (threads: 8)
----------------------------------------------------------
   Contigs      Chunks    Candidates   Time(s)  Rate(ch/s)
----------------------------------------------------------
         6          13           269        11       14.1M
        29          41          2929        35       44.7M
    ......       ......       ......    ......      ......
      2358        2387         10943       116       47.0M
----------------------------------------------------------
Found 16756 candidates in reference
Total primers to remove: 5057/5057 (dG <= -15000)

Summary
--------------------------------------------------------------------------------
PDR pair  Original (P/L/R)   Removed (P/L/R)     Final (P/L/R)         Threshold
--------------------------------------------------------------------------------
       0           5/46/73           5/16/27           0/30/46    -15000.000000
       1           5/80/19           5/27/14            0/53/5    -15000.000000
  ......            ......            ......            ......           ......
       5            5/6/48             0/0/0            5/6/48    -19000.000000*
--------------------------------------------------------------------------------
   TOTAL               930               411               519
--------------------------------------------------------------------------------
```

### Dimer Minimization Results
```
--------------------------------------------------
Algorithm                 Cost  Time(ms)    Winner
--------------------------------------------------
Random Search       -3223.0680      2485       ✓
SA                  -3223.0680     65825
Tabu Search         -3223.0680      5422
Genetic             -3223.0680     32909
--------------------------------------------------
 Winner: Random Search  cost=-3223.0680  (2485 ms)
```

### Final Solution Primers
```
Solution Primers
-------------------------------------------------------------------------------------------------------------------
   idx             PDR                    left seq     Tm    GC%                   right seq     Tm    GC%   size
-------------------------------------------------------------------------------------------------------------------
     0        [3, 383]        TTAAAATAGTGTCGCCGACG   58.2   45.0      GGATGTTCAGCCTCTAAAGGTT   59.9   45.5    419
     1      [252, 632]     CGAGTGTTAGATATCACACTGAG   58.0   43.5   TATAATAGTGAACGCTGAAAGGAGG   60.1   40.0    411
     2      [585, 965]      GTATTGAAATGTTTCACGCAGC   58.9   40.9      GATACGGAGGAGGAACCTAAGG   60.9   54.5    403
     ......      ......              ......       ......   ......              ......       ......   ......    ......
-------------------------------------------------------------------------------------------------------------------
```
## License

This project is licensed under the MIT License.

