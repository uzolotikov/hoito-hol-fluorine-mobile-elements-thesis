# Настройки конвейера. Скопируйте файл в config.sh и укажите свои пути:
#   cp config.example.sh config.sh
# Относительные пути считаются от корня репозитория.

THREADS=16

# Каталог баз. Его заполняет workflow/00_databases.sh.
DB=db

# Каталог результатов. Для каждой сборки создаётся подкаталог с её именем.
OUT=out

# Сборки в формате «имя:путь к FASTA», через пробел.
# Имена short, hybrid и flye используются в итоговых таблицах; порядок задаёт порядок столбцов.
ASSEMBLIES="short:data/HoitoGol_water_metaSPAdes_short-reads_contigs.fasta.gz hybrid:data/HoitoGol_water_metaSPAdes_hybrid_Illumina-Nanopore_contigs.fasta.gz flye:data/HoitoGol_water_metaFlye_Nanopore_Racon3-Polypolish3_contigs.fasta.gz"

# Необязательно: уже скачанная база geNomad 1.9 вместо скачивания в 00_databases.sh.
# GENOMAD_DB=/path/to/genomad_db
