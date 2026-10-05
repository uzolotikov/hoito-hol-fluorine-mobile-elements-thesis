# Входные данные

Каталог пуст. Скачайте сюда три сборки с Zenodo; в `config.example.sh` уже указаны их имена. Распаковывать не нужно.

| Файл | Что это | Где взять |
|---|---|---|
| `HoitoGol_water_metaSPAdes_short-reads_contigs.fasta.gz` | metaSPAdes 4.2.0, только прочтения Illumina | _[ЗАПОЛНИТЬ: номер в ENA или DOI на Zenodo]_ |
| `HoitoGol_water_metaSPAdes_hybrid_Illumina-Nanopore_contigs.fasta.gz` | metaSPAdes 4.2.0, прочтения Illumina и Nanopore | _[ЗАПОЛНИТЬ]_ |
| `HoitoGol_water_metaFlye_Nanopore_Racon3-Polypolish3_contigs.fasta.gz` | metaFlye 2.9.6, прочтения Nanopore от 1000 н., по три раунда полировки Racon и Polypolish | _[ЗАПОЛНИТЬ]_ |

Исходные прочтения: _[ЗАПОЛНИТЬ: номер BioProject и номера прочтений]_.

- Illumina, парные: образцы KI0305_0_45, KI0309_0_22, KI0309_0_45.
- Oxford Nanopore: библиотеки тех же проб, две из них после амплификации phi29.

Сборки можно получить из прочтений командами из `workflow/assembly.sh`. Команды восстановлены по журналам
программ; не сохранились только параметры Racon, там указаны параметры по умолчанию. Поэтому для точного
повторения лучше брать депонированные сборки: при пересборке могут измениться и контиги, и их имена. От имён
контигов зависит только столбец таксона в итоговой таблице, он берётся из `resources/locus_taxonomy.tsv`.

Таксоны локусов взяты из MAG, полученных биннингом объединённого набора контигов гибридной сборки и сборки metaFlye. Биннинг в этот конвейер не входит. MAG: _[ЗАПОЛНИТЬ: номера или DOI, если депонированы]_.

Контрольные суммы MD5 распакованных файлов:

```
ea7948cfec00c17787e0e06b6fda3a0e  HoitoGol_water_metaSPAdes_short-reads_contigs.fasta
990e9938020dc0bc6443597e702f642b  HoitoGol_water_metaSPAdes_hybrid_Illumina-Nanopore_contigs.fasta
265b3dd43e6794cf7ee84f962e492870  HoitoGol_water_metaFlye_Nanopore_Racon3-Polypolish3_contigs.fasta
```
