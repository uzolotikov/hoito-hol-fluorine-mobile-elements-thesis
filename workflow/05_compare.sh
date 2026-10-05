#!/usr/bin/env bash
# Сверка локусов между сборками: участок каждого локуса с флангами 300 п. н. ищется blastn в
# контигах каждой другой сборки. Итоговые таблицы пишутся в $OUT/results.
source "$(dirname "$0")/common.sh"
NAMES=(); for a in $ASSEMBLIES; do NAMES+=("${a%%:*}"); done
C="$OUT/compare"; mkdir -p "$C"
F="6 qseqid sseqid pident length qstart qend sstart send qlen slen evalue bitscore"
for B in "${NAMES[@]}"; do
  [ -e "$OUT/$B/blastdb.ndb" ] || makeblastdb -in "$OUT/$B/contigs1k.fna" -dbtype nucl -out "$OUT/$B/blastdb" > /dev/null
done
for A in "${NAMES[@]}"; do
  for B in "${NAMES[@]}"; do
    [ "$A" = "$B" ] && continue
    blastn -query "$OUT/$A/regions.fna" -db "$OUT/$B/blastdb" -outfmt "$F" -evalue 1e-50 -max_target_seqs 5 \
      -num_threads "$THREADS" > "$C/${A}_vs_${B}.tsv"
  done
done
python3 workflow/05_compare.py "$OUT" "${NAMES[@]}"
