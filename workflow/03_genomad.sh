#!/usr/bin/env bash
# geNomad на контигах с находками и контрольных контигах; маркерные гены плазмид на контигах с локусами.
# Запуск: workflow/03_genomad.sh <имя сборки>
source "$(dirname "$0")/common.sh"
NAME=$1; D="$OUT/$NAME"

seqkit grep -f "$D/genomad_ids.txt" "$D/contigs1k.fna" > "$D/genomad_input.fna"
rm -rf "$D/genomad_out"
genomad end-to-end --threads "$THREADS" --cleanup "$D/genomad_input.fna" "$D/genomad_out" "$DB/genomad_db" > "$D/genomad.log" 2>&1
log "$NAME: geNomad готов"

# Белки контигов с локусами и 29 профилей маркеров плазмид: репликация, мобилизация, конъюгация, разделение.
sed 's/[.[\*^$]/\\&/g; s/^/^/; s/$/_[0-9]+ /' "$D/locus_contigs.txt" > "$D/locus_contigs.regex"
seqkit grep -n -r -f "$D/locus_contigs.regex" "$D/proteins.faa" > "$D/locus_proteins.faa"
hmmsearch --cut_ga --cpu "$THREADS" --tblout "$D/plasmid_markers.tbl" -o /dev/null "$DB/plasmid_markers_29.hmm" "$D/locus_proteins.faa"
log "$NAME: маркеры плазмид: $(grep -vc '^#' "$D/plasmid_markers.tbl") находок"
