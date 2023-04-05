# Genome assembly statistics

A python script to calculate basic contiguity statistics of a genome assembly from a file in fasta format. Gaps are stretches of >=3 ambiguous sites in a row (i.e., NNN). Ambiguous sites (N) are excluded from contig-length and genome-size calculations.

## Requirements 
- Biopython

## Usage
python3 assembly_basic_stats.py [assembly_fasta]
