#!/usr/bin/env bash
# Как получены три сборки. Команды восстановлены по журналам SPAdes и Flye, заголовкам файлов выравнивания
# и скрипту полировки. Скрипт приведён для справки и в run_all.sh не входит: основной анализ начинается
# с готовых контигов (см. data/README.md).
#
# Входные прочтения:
#   short_R1.fastq, short_R2.fastq — Illumina, образцы KI0305_0_45, KI0309_0_22, KI0309_0_45, объединены
#     и очищены от последовательностей человека картированием bwa на GRCh38; 38,2 млн пар.
#   long.fastq.gz — Oxford Nanopore, объединённые библиотеки; 10,05 млн прочтений, 14,5 млрд н.
set -euo pipefail
T=${THREADS:-16}

# 1. Только короткие прочтения: metaSPAdes 4.2.0
spades.py --meta -t "$T" -m 120 -1 short_R1.fastq -2 short_R2.fastq -k 21,33,55,77,99 -o spades_short

# 2. Гибридная сборка: metaSPAdes 4.2.0, те же короткие прочтения и все длинные прочтения
spades.py --meta -t "$T" -m 120 -1 short_R1.fastq -2 short_R2.fastq --nanopore long.fastq.gz \
  -k 21,33,55,77,99 -o spades_hybrid

# 3. Только длинные прочтения: metaFlye 2.9.6, три раунда Racon и три раунда Polypolish
# 3.1. Фильтрация длинных прочтений: только по длине, от 1000 н., без порога качества.
# Чем фильтровали, не записано; эта команда даёт тот же набор: 5 331 134 прочтения, 11 814 861 378 н.
seqkit seq -m 1000 long.fastq.gz | gzip > long_filtered.fastq.gz

# 3.2. Сборка
flye --nano-raw long_filtered.fastq.gz --out-dir flye --meta --threads "$T"

# 3.3. Три раунда Racon 1.5.0 по длинным прочтениям, выравнивание minimap2 2.28.
# Команда racon не сохранилась; ниже запуск с параметрами по умолчанию.
draft=flye/assembly.fasta
for i in 1 2 3; do
  minimap2 -ax map-ont -t "$T" "$draft" long_filtered.fastq.gz > aln$i.sam
  racon -t "$T" long_filtered.fastq.gz aln$i.sam "$draft" > polished_r$i.fasta
  draft=polished_r$i.fasta
done

# 3.4. Три раунда Polypolish 0.6.1 по коротким прочтениям, выравнивание bwa 0.7.19 со всеми выравниваниями (-a)
for i in 1 2 3; do
  bwa index "$draft"
  bwa mem -t "$T" -a "$draft" short_R1.fastq > alignments_1.sam
  bwa mem -t "$T" -a "$draft" short_R2.fastq > alignments_2.sam
  polypolish filter --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
  polypolish polish "$draft" filtered_1.sam filtered_2.sam > polished_p$i.fasta
  draft=polished_p$i.fasta
done
# polished_p3.fasta — сборка metaFlye, использованная в анализе
