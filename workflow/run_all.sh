#!/usr/bin/env bash
# Весь анализ: базы, затем каждая сборка из config.sh, затем сверка сборок и итоговые таблицы.
source "$(dirname "$0")/common.sh"
bash workflow/00_databases.sh
for a in $ASSEMBLIES; do
  NAME=${a%%:*}; FA=${a#*:}
  bash workflow/01_annotate.sh "$NAME" "$FA"
  python3 workflow/02_loci.py "$OUT/$NAME" | tee "$OUT/$NAME/loci.log"
  bash workflow/03_genomad.sh "$NAME"
  python3 workflow/04_plasmids.py "$OUT/$NAME" | tee "$OUT/$NAME/plasmids.log"
done
bash workflow/05_compare.sh | tee "$OUT/compare.log"
