#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "usage: run_production_batch.sh RUN_DIR PARALLELISM [--resume]" >&2
  exit 2
fi

RUN_DIR="$1"
[[ "$RUN_DIR" = /* ]] || RUN_DIR="$PWD/$RUN_DIR"
PARALLELISM="$2"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROTOCOL_DIR="$SCRIPT_DIR"
REPO_DIR="$(cd "$PROTOCOL_DIR/../.." && pwd)"

RESUME="${3:-}"
if [[ -n "$RESUME" && "$RESUME" != "--resume" ]]; then exit 2; fi
if [[ -e "$RUN_DIR" && "$RESUME" != "--resume" ]]; then
 echo "Output exists; use --resume explicitly" >&2; exit 3
fi
if [[ "$RESUME" == "--resume" ]]; then
 cmp "$PROTOCOL_DIR/CODE_SHA256SUMS.txt" "$RUN_DIR/frozen_SHA256SUMS.txt"
 cmp "$PROTOCOL_DIR/repeated_simulation_spec.yml" "$RUN_DIR/frozen_spec.yml"
fi
if ! [[ "$PARALLELISM" =~ ^[1-9][0-9]*$ ]]; then
  echo "parallelism must be a positive integer" >&2
  exit 4
fi

cd "$REPO_DIR"
(cd "$PROTOCOL_DIR" && shasum -a 256 -c "CODE_SHA256SUMS.txt")
Rscript "$SCRIPT_DIR/R/preflight_production.R"
Rscript "$SCRIPT_DIR/R/test_rank_auc_orientation.R"

mkdir -p "$RUN_DIR/logs"
cp "$PROTOCOL_DIR/repeated_simulation_spec.yml" "$RUN_DIR/frozen_spec.yml"
cp "$PROTOCOL_DIR/CODE_SHA256SUMS.txt" "$RUN_DIR/frozen_SHA256SUMS.txt"
cp "$PROTOCOL_DIR/source_snapshot/DGP_generation.R" "$RUN_DIR/DGP_generation_professor_source.R"
cp "$PROTOCOL_DIR/README.md" "$RUN_DIR/PROTOCOL_README.md"
shasum -a 256 \
  "$PROTOCOL_DIR/source_snapshot/DGP_generation.R" \
  "$PROTOCOL_DIR/repeated_simulation_spec.yml" \
  "$PROTOCOL_DIR/R/generate_professor_dgps.R" \
  "$PROTOCOL_DIR/R/validate_professor_generators.R" \
  "$SCRIPT_DIR/R/production_common.R" \
  "$SCRIPT_DIR/R/test_rank_auc_orientation.R" \
  "$SCRIPT_DIR/R/preflight_production.R" \
  "$SCRIPT_DIR/R/run_one_replication.R" \
  "$SCRIPT_DIR/R/validate_stage.R" \
  "$SCRIPT_DIR/R/summarize_production.R" \
  "$SCRIPT_DIR/run_production_batch.sh" > "$RUN_DIR/production_code_SHA256.txt"
[[ -f "$RUN_DIR/STARTED_AT.txt" ]] || date -u +%Y-%m-%dT%H:%M:%SZ > "$RUN_DIR/STARTED_AT.txt"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

run_stage() {
  local dgp="$1"
  local reps="$2"
  mkdir -p "$RUN_DIR/dgp${dgp}" "$RUN_DIR/logs/dgp${dgp}"
  seq 1 "$reps" | xargs -P "$PARALLELISM" -I {} bash -c '
    rep="$1"
    dgp="$2"
    run_dir="$3"
    script="$4"
    log=$(printf "%s/logs/dgp%s/rep_%04d.log" "$run_dir" "$dgp" "$rep")
    completed=$(printf "%s/dgp%s/rep_%04d/COMPLETE" "$run_dir" "$dgp" "$rep")
    if [[ ! -f "$completed" ]]; then
      partial="${completed%/COMPLETE}.partial"
      if [[ -e "$partial" ]]; then
        archive="$run_dir/incomplete/dgp$dgp"
        mkdir -p "$archive"
        mv "$partial" "$archive/$(basename "$partial").$(date +%s).$$"
      fi
      Rscript "$script" --dgp="$dgp" --replication="$rep" --run_dir="$run_dir" >"$log" 2>&1
    fi
  ' _ {} "$dgp" "$RUN_DIR" "$SCRIPT_DIR/R/run_one_replication.R"
  Rscript "$SCRIPT_DIR/R/validate_stage.R" --dgp="$dgp" --run_dir="$RUN_DIR" | tee "$RUN_DIR/logs/dgp${dgp}_validation.log"
}

# Stages remain sequential; a failed validation prevents the next DGP.
run_stage 2 100
run_stage 3 100
run_stage 4 100

Rscript "$SCRIPT_DIR/R/summarize_production.R" --run_dir="$RUN_DIR" | tee "$RUN_DIR/logs/summary.log"
date -u +%Y-%m-%dT%H:%M:%SZ > "$RUN_DIR/COMPLETED_AT.txt"
touch "$RUN_DIR/PRODUCTION_COMPLETE"
