#!/bin/bash
set -e

# ==================== Description ====================
# version: v0.10
# update: 20260605
# implementation: 1) input query and reference are bot nucleic acid sequences, homolog scanning employs blastn

# ==================== Show help ====================
show_help() {
    cat << EOF
Usage: $0 -q <query directory> -r <reference directory> [-b <BLAST output directory>] [-e <extraction output directory>] [-h]

Required arguments:
  -q, --query     Directory containing query FASTA files (.fasta/.fna/.fa)
  -r, --ref       Directory containing reference FASTA files (.fasta/.fna/.fa)

Optional arguments:
  -b, --blast-out   BLAST output directory (default: oud_blastn)
  -e, --extract-out Output directory for extracted FASTA and filtered table (default: extracted_homologs)
  -h, --help        Show this help message

Examples:
  $0 -q ./query -r ./reference
  $0 -q ./query -r ./reference -b my_blast -e my_extract
EOF
    exit 0
}

# ==================== Parse arguments ====================
# Default values
blast_out="oud_blastn"
extract_out="extracted_homologs"
qued=""
refd=""

# Use getopt to support long options
OPTS=$(getopt -o q:r:b:e:h --long query:,ref:,blast-out:,extract-out:,help -n "$0" -- "$@")
if [ $? != 0 ]; then
    echo "Failed to parse arguments. Use -h for help."
    exit 1
fi
eval set -- "$OPTS"

while true; do
    case "$1" in
        -q|--query)
            qued="$2"
            shift 2
            ;;
        -r|--ref)
            refd="$2"
            shift 2
            ;;
        -b|--blast-out)
            blast_out="$2"
            shift 2
            ;;
        -e|--extract-out)
            extract_out="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Internal error!"
            exit 1
            ;;
    esac
done

# Check required arguments
if [ -z "$qued" ] || [ -z "$refd" ]; then
    echo "Error: Both query directory (-q) and reference directory (-r) must be specified."
    echo "Use -h for help."
    exit 1
fi

# Create output directories
mkdir -p "$blast_out" "$extract_out"

# ==================== Collect files ====================
# Support common FASTA extensions: .fasta / .fna / .fa
shopt -s nullglob
ref_files=("$refd"/*.fasta "$refd"/*.fna "$refd"/*.fa)
query_files=("$qued"/*.fasta "$qued"/*.fna "$qued"/*.fa)
shopt -u nullglob

if [ ${#ref_files[@]} -eq 0 ]; then
    echo "Error: No reference gene files found in $refd (supported: .fasta/.fna/.fa)"
    exit 1
fi
if [ ${#query_files[@]} -eq 0 ]; then
    echo "Error: No query FASTA files found in $qued (supported: .fasta/.fna/.fa)"
    exit 1
fi

# BLAST output column headers (consistent with -outfmt order)
header="qaccver\tsaccver\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle\tqacc\tqseqid"

# ==================== Processing function ====================
process_ref() {
    local ref_file=$1
    local ref_name=$2
    local ref_db_prefix="$blast_out/$ref_name"

    if [ ! -f "$ref_db_prefix.nhr" ]; then
        echo "[makeblastdb] Building database for reference $ref_name ..."
        makeblastdb -in "$ref_file" -out "$ref_db_prefix" -dbtype nucl
    fi

    for qfile in "${query_files[@]}"; do
        # Remove path and extension (supports .fasta/.fna/.fa)
        qname=$(basename "$qfile")
        qname="${qname%.fasta}"
        qname="${qname%.fna}"
        qname="${qname%.fa}"
        blast_tsv="$blast_out/${qname}_vs_${ref_name}.tsv"

        if [ -f "$blast_tsv" ]; then
            echo "[skip] BLAST result $blast_tsv already exists. Skipping alignment."
        else
            echo "[blastn] Aligning $qname against reference $ref_name ..."
            (echo -e "$header" && \
             blastn -query "$qfile" -db "$ref_db_prefix" \
                -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle qacc qseqid' \
                -num_threads "$(nproc)" \
                -word_size 21) > "$blast_tsv"
        fi

        prefix="${extract_out}/${qname}_${ref_name}"
        extract_fasta="${prefix}.fasta"
        if [ -f "$extract_fasta" ]; then
            echo "[skip] Extraction result $extract_fasta already exists. Skipping extraction."
        else
            echo "[extract] Extracting regions matching $ref_name from $qname ..."
            python3 extract_homolog_regions.py -b "$blast_tsv" -q "$qfile" -o "$prefix"
        fi
    done
}

# ==================== Main loop ====================
for ref_file in "${ref_files[@]}"; do
    # Remove path and extension (supports .fasta/.fna/.fa)
    ref_basename=$(basename "$ref_file")
    ref_basename="${ref_basename%.fasta}"
    ref_basename="${ref_basename%.fna}"
    ref_basename="${ref_basename%.fa}"
    echo "========== Processing reference gene: $ref_basename =========="
    process_ref "$ref_file" "$ref_basename"
done

echo "All done! Extracted sequences are saved in $extract_out/"
