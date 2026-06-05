#!/bin/bash
set -e

# ==================== 显示帮助 ====================
show_help() {
    cat << EOF
用法: $0 -q <查询序列目录> -r <参考基因目录> [-b <BLAST输出目录>] [-e <提取输出目录>] [-h]

必需参数:
  -q, --query     查询序列所在目录（包含 .fasta/.fna/.fa 文件）
  -r, --ref       参考基因所在目录（包含 .fasta/.fna/.fa 文件）

可选参数:
  -b, --blast-out BLAST 结果输出目录（默认: oud_blastn）
  -e, --extract-out 提取的 FASTA 和过滤表输出目录（默认: extracted_homologs）
  -h, --help      显示此帮助信息

示例:
  $0 -q ./query -r ./reference
  $0 -q ./query -r ./reference -b my_blast -e my_extract
EOF
    exit 0
}

# ==================== 解析参数 ====================
# 默认值
blast_out="oud_blastn"
extract_out="extracted_homologs"
qued=""
refd=""

# 使用 getopt 支持长选项
OPTS=$(getopt -o q:r:b:e:h --long query:,ref:,blast-out:,extract-out:,help -n "$0" -- "$@")
if [ $? != 0 ]; then
    echo "参数解析失败，请使用 -h 查看帮助"
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
            echo "内部错误！"
            exit 1
            ;;
    esac
done

# 检查必需参数
if [ -z "$qued" ] || [ -z "$refd" ]; then
    echo "错误：必须指定查询序列目录 (-q) 和参考基因目录 (-r)"
    echo "使用 -h 查看帮助"
    exit 1
fi

# 创建输出目录
mkdir -p "$blast_out" "$extract_out"

# ==================== 收集文件 ====================
# 支持 .fasta / .fna / .fa 等常见后缀
shopt -s nullglob
ref_files=("$refd"/*.fasta "$refd"/*.fna "$refd"/*.fa)
query_files=("$qued"/*.fasta "$qued"/*.fna "$qued"/*.fa)
shopt -u nullglob

if [ ${#ref_files[@]} -eq 0 ]; then
    echo "错误：在 $refd 中未找到参考基因文件（*.fasta/fna/fa）"
    exit 1
fi
if [ ${#query_files[@]} -eq 0 ]; then
    echo "错误：在 $qued 中未找到查询fasta文件（*.fasta/fna/fa）"
    exit 1
fi

# BLAST 输出列名（与 -outfmt 顺序一致）
header="qaccver\tsaccver\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tstitle\tqacc\tqseqid"

# ==================== 处理函数 ====================
process_ref() {
    local ref_file=$1
    local ref_name=$2
    local ref_db_prefix="$blast_out/$ref_name"

    if [ ! -f "$ref_db_prefix.nhr" ]; then
        echo "[makeblastdb] 为参考基因 $ref_name 建立数据库 ..."
        makeblastdb -in "$ref_file" -out "$ref_db_prefix" -dbtype nucl
    fi

    for qfile in "${query_files[@]}"; do
        # 去除路径和扩展名（支持 .fasta/.fna/.fa）
        qname=$(basename "$qfile")
        qname="${qname%.fasta}"
        qname="${qname%.fna}"
        qname="${qname%.fa}"
        blast_tsv="$blast_out/${qname}_vs_${ref_name}.tsv"

        if [ -f "$blast_tsv" ]; then
            echo "[skip] BLAST结果 $blast_tsv 已存在，跳过比对"
        else
            echo "[blastn] 将 $qname 比对到参考基因 $ref_name ..."
            (echo -e "$header" && \
             blastn -query "$qfile" -db "$ref_db_prefix" \
                -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle qacc qseqid' \
                -num_threads "$(nproc)" \
                -word_size 21) > "$blast_tsv"
        fi

        prefix="${extract_out}/${qname}_${ref_name}"
        extract_fasta="${prefix}.fasta"
        if [ -f "$extract_fasta" ]; then
            echo "[skip] 提取结果 $extract_fasta 已存在，跳过提取"
        else
            echo "[extract] 从 $qname 中提取匹配 $ref_name 的区间 ..."
            python3 extract_homolog_regions.py -b "$blast_tsv" -q "$qfile" -o "$prefix"
        fi
    done
}

# ==================== 主循环 ====================
for ref_file in "${ref_files[@]}"; do
    # 去除路径和扩展名（支持 .fasta/.fna/.fa）
    ref_basename=$(basename "$ref_file")
    ref_basename="${ref_basename%.fasta}"
    ref_basename="${ref_basename%.fna}"
    ref_basename="${ref_basename%.fa}"
    echo "========== 处理参考基因：$ref_basename =========="
    process_ref "$ref_file" "$ref_basename"
done

echo "全部完成！提取的序列保存在 $extract_out/ 目录中。"