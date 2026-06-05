# HUNT - Homology-based Utility for Nucleotide Target extraction
## Description
HUNT is a lightweight command-line tool designed for homolog identification from bacterial genomes. It extracts homologous nucleotide regions from query sequences (e.g., draft genomes) by searching against one or more reference gene sequences (nucleotide or protein) using homology search tools. In current version, only homolog scan using blastn is implemented, and the tool is still under active development. The workflow handles multiple references, deduplicates overlapping hits, removes fully contained intervals, and outputs clean FASTA sequences along with filtered BLAST tables — all without external Python dependencies beyond the standard library.

## Dependencies
* Python 3.10+ (tested on 3.10.12) – no extra Python packages required
* BLAST+ (makeblastdb, blastn) – must be installed and available in $PATH
* Tested on Ubuntu 22.04 LTS.

## Installation
```bash
git clone https://github.com/yourusername/hunt.git
cd hunt
chmod +x run_blastn.sh
```

## Usage
./run_blastn.sh -q /path/to/query/sequences -r /path/to/reference/genes -b blast-out -e extracted-homolog

| Option | Long form         | Description                                   | Default             |
|--------|-------------------|-----------------------------------------------|---------------------|
| `-q`   | `--query`         | Directory containing query FASTA files       | required            |
| `-r`   | `--ref`           | Directory containing reference FASTA files   | required            |
| `-b`   | `--blast-out`     | Output directory for BLAST results           | `oud_blastn`        |
| `-e`   | `--extract-out`   | Output directory for extracted FASTA/TSV     | `extracted_homologs`|
| `-h`   | `--help`          | Show help message                            | -                   |

## Example
Test files are provided in the test_files/ directory. You can run the tool with the following command:
```
bash main.sh -q test_files/query -r test_files/reference -b oud_blastn -e extracted_homologs
```

**Note**: this is optimized for complete circular chromosomes or near‑complete genome assemblies when aligning against full‑length reference gene sequences (or similarly sized regions). When applied to draft genome assemblies (e.g., contigs or scaffolds), the extracted homologous regions may be shorter or fragmented due to assembly gaps or incomplete gene coverage.
