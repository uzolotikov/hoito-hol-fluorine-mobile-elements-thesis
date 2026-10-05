# Общие функции: читается остальными скриптами конвейера.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
if [ ! -f config.sh ]; then
  echo "Нет config.sh. Скопируйте config.example.sh в config.sh и укажите пути." >&2
  exit 1
fi
# shellcheck source=/dev/null
source config.sh
log() { echo "[$(date '+%F %T')] $*" >&2; }
