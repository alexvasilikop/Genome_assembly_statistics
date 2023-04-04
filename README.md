# Genome assembly statistics

A python script to calculate basic contiguity statistics of a genome assembly from a file in fasta format. Any number of Ns in a row (>=1) are counted as gaps. Ambiguous sites are excluded from contig length and genome size calculations

## Requirements 
- Biopython

## Usage
python3 assembly_basic_stats.py [assembly_fasta]
