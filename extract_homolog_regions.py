#!/usr/bin/env python3
"""
extract_homolog_regions.py

Extract homologous regions from BLAST results (pure standard library, no Biopython dependency):
- Filter alignments with length < 50 and coverage (scov = length/slen) < 0.05
- Group by qaccver
- Deduplicate identical coordinates (keep best scov * pident)
- Remove intervals fully contained within another (keep outer interval)
- Keep all partially overlapping intervals (do not discard)
- Output extracted FASTA and filtered BLAST table

Usage:
    python extract_homolog_regions.py -b blast.tsv -q query.fasta -o output_prefix
"""

import sys
import argparse
from collections import defaultdict

# ------------------------------ FASTA pharse ---------------------------------
def read_fasta_dict(fasta_file):
    """
    Read FASTA file and return a dictionary: {seq_id: seq_string}
    Assumes sequence ID is the first non-whitespace string after '>' (until space or newline)
    """
    seq_dict = {}
    current_id = None
    seq_lines = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Save previous sequence
                if current_id is not None:
                    seq_dict[current_id] = ''.join(seq_lines)
                # Parse new ID: take first word after '>' (split by space or tab)
                header = line[1:].strip()
                current_id = header.split()[0]   # Keep only first word as ID
                seq_lines = []
            else:
                seq_lines.append(line)
        # Last sequence
        if current_id is not None:
            seq_dict[current_id] = ''.join(seq_lines)
    return seq_dict

def write_fasta(records, output_file):
    """
    Write list of (id, seq) to FASTA file
    """
    with open(output_file, 'w') as f:
        for seq_id, seq in records:
            f.write(f">{seq_id}\n")
            # Write 80 bases per line (optional)
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

# ------------------------------ BLAST processing  ---------------------------------
def parse_blast_with_headers(blast_file):
    """Read BLAST result file with header line (tab-separated)"""
    with open(blast_file) as f:
        header_line = f.readline().strip()
        if not header_line:
            raise ValueError("BLAST BLAST file is empty")
        columns = header_line.split('\t')
        col_idx = {col: i for i, col in enumerate(columns)}
        required = ['qaccver', 'saccver', 'pident', 'length', 'qstart', 'qend', 'slen']
        missing = [c for c in required if c not in col_idx]
        if missing:
            raise ValueError(f"BLAST file is in lack of essential columns: {missing}")
        # stitle is optional; set index to -1 if not present
        if 'stitle' not in col_idx:
            col_idx['stitle'] = -1
        rows = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < len(required):
                continue
            rows.append((line, parts))
    return col_idx, rows

def filter_by_length_and_scov(rows, col_idx, min_len=50, min_scov=0.05):
    """Filter: length >= min_len and scov = length/slen >= min_scov"""
    filtered = []
    for line, parts in rows:
        length = float(parts[col_idx['length']])
        slen = float(parts[col_idx['slen']])
        if length < min_len:
            continue
        scov = length / slen
        if scov < min_scov:
            continue
        filtered.append((line, parts, scov))
    return filtered

def group_by_qaccver(filtered_rows, col_idx):
    """Group by qaccver; each item contains (line, parts, scov, qaccver, qstart, qend, score)"""
    groups = defaultdict(list)
    for line, parts, scov in filtered_rows:
        qaccver = parts[col_idx['qaccver']]
        qstart = int(parts[col_idx['qstart']])
        qend = int(parts[col_idx['qend']])
        if qstart > qend:
            qstart, qend = qend, qstart
        pident = float(parts[col_idx['pident']])
        score = scov * pident
        groups[qaccver].append((line, parts, scov, qaccver, qstart, qend, score))
    return groups

def deduplicate_by_coords(intervals):
    """Deduplicate identical coordinates, keep the one with highest score"""
    best = {}
    for itv in intervals:
        key = (itv[4], itv[5])  # (qstart, qend)
        if key not in best or itv[6] > best[key][6]:
            best[key] = itv
    return list(best.values())

def remove_contained(intervals):
    """
    Remove intervals that are fully contained within another interval.
    Note: Only fully contained intervals are removed; partially overlapping ones are kept.
    Returns the remaining list.
    """
    if len(intervals) <= 1:
        return intervals
    intervals_sorted = sorted(intervals, key=lambda x: (x[5] - x[4]), reverse=True)
    keep = []
    for itv in intervals_sorted:
        contained = False
        for kept in keep:
            if kept[4] <= itv[4] and kept[5] >= itv[5]:
                contained = True
                break
        if not contained:
            keep.append(itv)
    return keep

def process_one_qaccver(intervals):
    """Deduplicate -> remove contained -> return kept intervals"""
    if not intervals:
        return []
    step1 = deduplicate_by_coords(intervals)
    step2 = remove_contained(step1)
    return step2

def extract_sequences(intervals, query_fasta_dict, col_idx):
    """
    Extract subsequences from query FASTA dict based on kept intervals.
    Returns list of [(seq_id, seq_string), ...]
    """
    records = []
    for itv in intervals:
        qaccver = itv[3]
        qstart = itv[4]
        qend = itv[5]
        saccver = itv[1][col_idx['saccver']]
        if qaccver not in query_fasta_dict:
            print(f"Warning: {qaccver} not found in query FASTA, skipping", file=sys.stderr)
            continue
        full_seq = query_fasta_dict[qaccver]
        # 边界保护
        if qstart < 1:
            qstart = 1
        if qend > len(full_seq):
            qend = len(full_seq)
        sub_seq = full_seq[qstart-1:qend]   # 1-based to 0-based Python slice
        seq_id = f"{qaccver}_{saccver}_{qstart}_{qend}"
        records.append((seq_id, sub_seq))
    return records

def output_filtered_table(intervals, original_headers, output_tsv):
    """Output filtered BLAST table (preserve original header)"""
    with open(output_tsv, 'w') as f:
        f.write(original_headers + '\n')
        for itv in intervals:
            f.write(itv[0] + '\n')
    print(f"Filtered table -> {output_tsv} ( {len(intervals)} )")

# ------------------------------ main function ---------------------------------
def main():
    parser = argparse.ArgumentParser(description="Extract homologous regions from BLAST results (pure standard library, no Biopython)")
    parser.add_argument("-b", "--blast", required=True, help="BLAST result file (tab-separated with header line)")
    parser.add_argument("-q", "--query", required=True, help="Query sequence FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output prefix (will generate .fasta and .filtered.tsv)")
    parser.add_argument("--min_len", type=int, default=50, help="Minimum alignment length, default 50")
    parser.add_argument("--min_scov", type=float, default=0.05, help="Minimum coverage (length/slen), default 0.05")
    args = parser.parse_args()

    # 1. Read BLAST file
    print("Reading BLAST file ...")
    col_idx, rows = parse_blast_with_headers(args.blast)
    original_headers = open(args.blast).readline().strip()

    # 2. Filter by length and coverage
    print(f"Filtering  length<{args.min_len} and scov<{args.min_scov} ...")
    filtered = filter_by_length_and_scov(rows, col_idx, args.min_len, args.min_scov)
    if not filtered:
        print("Warning: No rows passed filtering, output empty")
        open(args.output + ".fasta", 'w').close()
        with open(args.output + ".filtered.tsv", 'w') as f:
            f.write(original_headers + '\n')
        return
    print(f"  {len(filtered)} rows remaining after filtering")

    # 3. Group by qaccver
    groups = group_by_qaccver(filtered, col_idx)

    # 4. Process each group
    all_selected = []
    for qaccver, intervals in groups.items():
        selected = process_one_qaccver(intervals)
        all_selected.extend(selected)
        print(f"  {qaccver}: original {len(intervals)} -> dedup+contained removal {len(selected)} intervals")

    print("Reading query FASTA...")
    query_dict = read_fasta_dict(args.query)

    print("Extracting subsequences...")
    records = extract_sequences(all_selected, query_dict, col_idx)
    
    # Write FASTA
    fasta_out = args.output + ".fasta"
    write_fasta(records, fasta_out)
    print(f"Extracted {len(records)} sequences -> {fasta_out}")

    # Write filtered table
    output_filtered_table(all_selected, original_headers, args.output + ".filtered.tsv")

    print("Done")

if __name__ == "__main__":
    main()
