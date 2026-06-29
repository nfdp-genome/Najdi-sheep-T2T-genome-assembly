#!/bin/bash
set -euo pipefail

ensure_module_command() {
  if command -v module >/dev/null 2>&1; then
    return 0
  fi
  if [[ -f /etc/profile.d/modules.sh ]]; then
    # shellcheck disable=SC1091
    source /etc/profile.d/modules.sh
  elif [[ -f /usr/share/Modules/init/bash ]]; then
    # shellcheck disable=SC1091
    source /usr/share/Modules/init/bash
  fi
  command -v module >/dev/null 2>&1
}

load_required_modules_or_fail() {
  local module_list="${1:-}"
  local script_name="${2:-script}"
  [[ -n "${module_list}" ]] || return 0

  if ! ensure_module_command; then
    echo "[ERROR] ${script_name}: environment module command is unavailable." >&2
    return 1
  fi

  local module_name
  for module_name in ${module_list}; do
    if ! module load "${module_name}"; then
      echo "[ERROR] ${script_name}: failed to load module '${module_name}'." >&2
      return 1
    fi
  done
}
