#!/usr/bin/env bash
#SBATCH --job-name=step_salsa
#SBATCH --partition=batch
#SBATCH --constraint=intel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=256G
#SBATCH --time=336:00:00
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail
MODULES="yahs/1.1 samtools/1.16.1 bedtools/2.30.0"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ORIG_SCRIPT="${PIPELINE_ROOT}/scripts/step_salsa.sh"

CONFIG=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      [[ $# -ge 2 ]] || { echo "[ERROR] --config requires a YAML path" >&2; exit 1; }
      CONFIG="$2"; shift 2
      ;;
    *) echo "[ERROR] Unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -n "${CONFIG}" && -s "${CONFIG}" ]] || { echo "[ERROR] Missing --config YAML" >&2; exit 1; }
[[ -s "${ORIG_SCRIPT}" ]] || { echo "[ERROR] Missing source script: ${ORIG_SCRIPT}" >&2; exit 1; }

if [[ -n "${MODULES}" ]]; then
  source /etc/profile.d/modules.sh 2>/dev/null || true
  if declare -F module >/dev/null 2>&1 || command -v module >/dev/null 2>&1; then
    for m in ${MODULES}; do module load "$m"; done
  else
    echo "[WARN] Environment modules are unavailable; expecting tools on PATH."
  fi
fi

python3 -c "import yaml" >/dev/null 2>&1 || { echo "[ERROR] PyYAML missing" >&2; exit 1; }

eval "$(python3 - "$CONFIG" <<'PYEOF'
import sys, yaml, shlex
cfg = yaml.safe_load(open(sys.argv[1])) or {}
params = cfg.get('params', {})
for k,v in params.items():
    if isinstance(v, list):
        vv='|'.join(str(x) for x in v)
    elif v is None:
        vv=''
    else:
        vv=str(v)
    print(f"{k}={shlex.quote(vv)}")
print('log_file=' + shlex.quote(str(cfg.get('log_file','./logs/step_salsa.log'))))
PYEOF
)"

mkdir -p "$(dirname "$log_file")"
exec > >(tee -a "$log_file") 2>&1

# Export all keys loaded from YAML params
for key in $(python3 - "$CONFIG" <<'PYEOF'
import sys, yaml
cfg = yaml.safe_load(open(sys.argv[1])) or {}
for k in (cfg.get('params', {}) or {}).keys():
    print(k)
PYEOF
); do
  export "$key"
done

bash "${ORIG_SCRIPT}"
