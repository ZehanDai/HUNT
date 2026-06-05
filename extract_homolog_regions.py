#!/usr/bin/env python3
"""
extract_homolog_regions.py

基于 BLAST 结果提取同源区段（纯标准库，无 Biopython 依赖）：
- 过滤 length<50 且 scov<0.05
- 按 qaccver 分组
- 相同坐标去重（保留 scov*pident 最大）
- 移除被完全包含的区间（保留外层区间）
- 部分重叠的区间全部保留（不剔除）
- 输出提取的 FASTA 和过滤后的 BLAST 表格

用法：
    python extract_homolog_regions.py -b blast.tsv -q query.fasta -o output_prefix
"""

import sys
import argparse
from collections import defaultdict

# ------------------------------ FASTA 解析（纯手写）---------------------------------
def read_fasta_dict(fasta_file):
    """
    读取 FASTA 文件，返回字典：{seq_id: seq_string}
    假设序列 ID 为 '>' 后第一个非空白连续字符串（直到遇到空格或换行）
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
                # 保存上一条序列
                if current_id is not None:
                    seq_dict[current_id] = ''.join(seq_lines)
                # 解析新 ID：取 '>' 后的第一个单词（空格、制表符分隔）
                header = line[1:].strip()
                current_id = header.split()[0]   # 只保留第一个单词作为 ID
                seq_lines = []
            else:
                seq_lines.append(line)
        # 最后一条序列
        if current_id is not None:
            seq_dict[current_id] = ''.join(seq_lines)
    return seq_dict

def write_fasta(records, output_file):
    """
    将 [(id, seq), ...] 写入 FASTA 文件
    """
    with open(output_file, 'w') as f:
        for seq_id, seq in records:
            f.write(f">{seq_id}\n")
            # 每行 80 个碱基（可选）
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + "\n")

# ------------------------------ BLAST 处理函数 ---------------------------------
def parse_blast_with_headers(blast_file):
    """读取带列名行的 BLAST 结果文件"""
    with open(blast_file) as f:
        header_line = f.readline().strip()
        if not header_line:
            raise ValueError("BLAST 文件为空")
        columns = header_line.split('\t')
        col_idx = {col: i for i, col in enumerate(columns)}
        required = ['qaccver', 'saccver', 'pident', 'length', 'qstart', 'qend', 'slen']
        missing = [c for c in required if c not in col_idx]
        if missing:
            raise ValueError(f"BLAST 文件缺少必要列: {missing}")
        # stitle 不是必须的，如果没有则设索引为 -1
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
    """过滤：length >= min_len 且 scov = length/slen >= min_scov"""
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
    """按 qaccver 分组，每个元素包含 (line, parts, scov, qaccver, qstart, qend, score)"""
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
    """相同坐标去重，保留 score 最大的"""
    best = {}
    for itv in intervals:
        key = (itv[4], itv[5])  # (qstart, qend)
        if key not in best or itv[6] > best[key][6]:
            best[key] = itv
    return list(best.values())

def remove_contained(intervals):
    """
    移除被其他区间完全包含的区间。
    注意：不处理部分重叠，被包含的才移除。
    返回剩余的列表。
    """
    if len(intervals) <= 1:
        return intervals
    # 按区间长度降序排序（长的在前）
    intervals_sorted = sorted(intervals, key=lambda x: (x[5] - x[4]), reverse=True)
    keep = []
    for itv in intervals_sorted:
        contained = False
        for kept in keep:
            # 如果 kept 完全包含 itv
            if kept[4] <= itv[4] and kept[5] >= itv[5]:
                contained = True
                break
        if not contained:
            keep.append(itv)
    return keep

def process_one_qaccver(intervals):
    """去重 -> 去包含 -> 返回保留区间"""
    if not intervals:
        return []
    step1 = deduplicate_by_coords(intervals)
    step2 = remove_contained(step1)
    return step2

def extract_sequences(intervals, query_fasta_dict, col_idx):
    """
    根据保留的区间，从 query FASTA 字典中提取子序列。
    返回列表 [(seq_id, seq_string), ...]
    """
    records = []
    for itv in intervals:
        qaccver = itv[3]
        qstart = itv[4]
        qend = itv[5]
        saccver = itv[1][col_idx['saccver']]
        if qaccver not in query_fasta_dict:
            print(f"警告: {qaccver} 在查询 FASTA 中未找到，跳过", file=sys.stderr)
            continue
        full_seq = query_fasta_dict[qaccver]
        # 边界保护
        if qstart < 1:
            qstart = 1
        if qend > len(full_seq):
            qend = len(full_seq)
        sub_seq = full_seq[qstart-1:qend]   # 1-based 转 Python 0-based
        seq_id = f"{qaccver}_{saccver}_{qstart}_{qend}"
        records.append((seq_id, sub_seq))
    return records

def output_filtered_table(intervals, original_headers, output_tsv):
    """输出过滤后的 BLAST 表格（保留原表头）"""
    with open(output_tsv, 'w') as f:
        f.write(original_headers + '\n')
        for itv in intervals:
            f.write(itv[0] + '\n')
    print(f"过滤后表格 -> {output_tsv} (共 {len(intervals)} 行)")

# ------------------------------ 主函数 ---------------------------------
def main():
    parser = argparse.ArgumentParser(description="基于 BLAST 结果提取同源区段（纯标准库，无 Biopython）")
    parser.add_argument("-b", "--blast", required=True, help="BLAST 结果文件（带列名行，制表符分隔）")
    parser.add_argument("-q", "--query", required=True, help="查询序列 FASTA 文件")
    parser.add_argument("-o", "--output", required=True, help="输出文件前缀（将生成 .fasta 和 .filtered.tsv）")
    parser.add_argument("--min_len", type=int, default=50, help="最小比对长度，默认 50")
    parser.add_argument("--min_scov", type=float, default=0.05, help="最小覆盖度 (length/slen)，默认 0.05")
    args = parser.parse_args()

    # 1. 读取 BLAST 文件
    print("读取 BLAST 文件...")
    col_idx, rows = parse_blast_with_headers(args.blast)
    original_headers = open(args.blast).readline().strip()

    # 2. 过滤长度和覆盖度
    print(f"过滤 length<{args.min_len} 且 scov<{args.min_scov} ...")
    filtered = filter_by_length_and_scov(rows, col_idx, args.min_len, args.min_scov)
    if not filtered:
        print("警告：没有通过过滤的行，输出为空")
        open(args.output + ".fasta", 'w').close()
        with open(args.output + ".filtered.tsv", 'w') as f:
            f.write(original_headers + '\n')
        return
    print(f"  过滤后剩余 {len(filtered)} 行")

    # 3. 按 qaccver 分组
    groups = group_by_qaccver(filtered, col_idx)

    # 4. 对每个组进行处理
    all_selected = []
    for qaccver, intervals in groups.items():
        selected = process_one_qaccver(intervals)
        all_selected.extend(selected)
        print(f"  {qaccver}: 原始 {len(intervals)} -> 去重+去包含后 {len(selected)} 个区间")

    # 5. 读取查询 FASTA
    print("读取查询 FASTA...")
    query_dict = read_fasta_dict(args.query)

    # 6. 提取序列
    print("提取子序列...")
    records = extract_sequences(all_selected, query_dict, col_idx)

    # 7. 输出 FASTA
    fasta_out = args.output + ".fasta"
    write_fasta(records, fasta_out)
    print(f"提取 {len(records)} 条序列 -> {fasta_out}")

    # 8. 输出过滤后的 BLAST 表格
    output_filtered_table(all_selected, original_headers, args.output + ".filtered.tsv")

    print("完成。")

if __name__ == "__main__":
    main()