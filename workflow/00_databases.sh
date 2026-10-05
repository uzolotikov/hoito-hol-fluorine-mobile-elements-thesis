#!/usr/bin/env bash
# Базы для конвейера.
# Профили Pfam и модель Rfam лежат в resources/models в том виде, в каком использованы в работе:
#   crcb, clc, plasmid_markers_29 — Pfam 33.1; mge_146 — 140 профилей Pfam 33.1 и 6 из более позднего релиза
#   (номера и версии в resources/mge_profiles.tsv); RF01734 — Rfam 15.1. Пороги GA берутся из этих файлов.
# Скачивается только база geNomad 1.9; если в config.sh задан GENOMAD_DB, используется она.
source "$(dirname "$0")/common.sh"
mkdir -p "$DB"
for f in crcb.hmm clc.hmm mge_146.hmm plasmid_markers_29.hmm RF01734.cm; do
  gunzip -c "resources/models/$f.gz" > "$DB/$f"
done
( cd "$DB" && md5sum -c --quiet "$ROOT/resources/models.md5" ) || { echo "ОШИБКА: контрольные суммы моделей не совпали" >&2; exit 1; }

if [ -n "${GENOMAD_DB:-}" ]; then
  ln -sfn "$(cd "$GENOMAD_DB" && pwd)" "$DB/genomad_db"
elif [ ! -d "$DB/genomad_db" ]; then
  log "скачиваю базу geNomad 1.9 (Zenodo, запись 14886553, архив около 850 МБ)"
  wget -q -O - "https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz?download=1" | tar -xz -C "$DB"
fi
log "базы готовы в $DB"
