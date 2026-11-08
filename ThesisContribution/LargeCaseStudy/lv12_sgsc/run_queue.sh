#!/usr/bin/env bash
# Usage:  ./run_queue.sh queue_tonight.txt 4
# Each queue line:  <CONFIG> [ENV=VALUE ...]   (lines starting with # are skipped)
set -u
QUEUE=${1:?queue file}; N=${2:-4}; JL_PROJ=${JL_PROJ:-.}
mkdir -p logs
grep -v '^[[:space:]]*#' "$QUEUE" | grep -v '^[[:space:]]*$' | \
xargs -P "$N" -L 1 bash -c '
  cfg=$0; shift
  tag=$(echo "$cfg ${*:-base}" | tr " =" "__")
  echo "$(date +%H:%M) start $tag"
  env LV12_NINIT=5 LV12_NNEW=5 LV12_CONFIGS="$cfg" "$@" \
      julia --project="'"$JL_PROJ"'" lv12_opf.jl > "logs/$tag.log" 2>&1
  echo "$(date +%H:%M) done  $tag (exit $?)"
'