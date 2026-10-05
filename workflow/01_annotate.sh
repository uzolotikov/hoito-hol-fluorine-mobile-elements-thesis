#!/usr/bin/env bash
# Аннотация одной сборки: контиги >= 1000 п. н., гены Prodigal, профили Pfam с порогами GA,
# рибопереключатель RF01734, дерево белков CLC с референсами.
# Запуск: workflow/01_annotate.sh <имя сборки> <FASTA с контигами>
source "$(dirname "$0")/common.sh"
NAME=$1; FA=$2; D="$OUT/$NAME"
mkdir -p "$D"

seqkit seq -m 500 "$FA" | seqkit stats -a -T > "$D/stats500.tsv"
seqkit seq -m 1000 "$FA" > "$D/contigs1k.fna"
seqkit fx2tab -n -i -l "$D/contigs1k.fna" | cut -f1,2 > "$D/contig_lengths.tsv"
log "$NAME: контигов >= 1000 п. н.: $(wc -l < "$D/contig_lengths.tsv")"

# Prodigal обрабатывает каждый контиг независимо, поэтому разбиение на части на результат не влияет.
rm -rf "$D/split" "$D/prod"; mkdir "$D/prod"
seqkit split2 -p "$THREADS" -O "$D/split" "$D/contigs1k.fna" > /dev/null 2>&1
ls "$D/split" | xargs -P "$THREADS" -I{} prodigal -p meta -q -i "$D/split/{}" -a "$D/prod/{}.faa" -o /dev/null
cat "$D"/prod/*.faa > "$D/proteins.faa"
rm -rf "$D/split" "$D/prod"
log "$NAME: белков $(grep -c '^>' "$D/proteins.faa")"

for m in crcb clc mge_146; do
  hmmsearch --cut_ga --cpu "$THREADS" --tblout "$D/${m%_146}.tbl" -o /dev/null "$DB/$m.hmm" "$D/proteins.faa"
done
# -Z фиксирует размер базы для E-значений; отбор находок идёт по порогам GA и от -Z не зависит.
cmsearch --cut_ga -Z 2958934 --cpu "$THREADS" --tblout "$D/ribo.tbl" -o /dev/null "$DB/RF01734.cm" "$D/contigs1k.fna"
log "$NAME: CrcB $(grep -vc '^#' "$D/crcb.tbl"), CLC $(grep -vc '^#' "$D/clc.tbl"), RF01734 $(grep -vc '^#' "$D/ribo.tbl")"

# Дерево CLC: клада CLC^F отделяется от хлоридных каналов по референсам из resources/clc_refs.faa.
grep -v '^#' "$D/clc.tbl" | awk '{print $1}' | sort -u > "$D/clc_ids.txt"
# Белки сортируются по имени: порядок входа влияет на дерево, а без сортировки он зависел бы от THREADS.
seqkit grep -f "$D/clc_ids.txt" "$D/proteins.faa" | seqkit seq -i | seqkit sort -n -N | sed 's/\*$//' > "$D/clc_hits.faa"
cat resources/clc_refs.faa "$D/clc_hits.faa" > "$D/clc_all.faa"
mafft --auto --thread "$THREADS" "$D/clc_all.faa" > "$D/clc_all.aln" 2> /dev/null
FastTree -lg -quiet "$D/clc_all.aln" > "$D/clc_all.nwk" 2> /dev/null
log "$NAME: аннотация готова"
