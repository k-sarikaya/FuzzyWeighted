#!/usr/bin/env bash
# Reproduction driver (Linux/macOS).
#   ./reproduce.sh          unit tests + fast artifacts + figures + manuscript
#   ./reproduce.sh --full   also re-solve every grid cell and the scaling family (hours)
set -euo pipefail
cd "$(dirname "$0")"
export PYTHONPATH="$PWD/src"
run() { echo ">>> $*"; python3 "experiments/$1" "${@:2}"; }

python3 -m unittest discover -s tests -v
if [[ "${1:-}" == "--full" ]]; then
  run verify_pharma_oracle.py
  run verify_fyffe_dp.py
  run run_campaign_v2.py
  run generate_benchmark_tables_v6.py
  run parametric_stability.py --cpp
  run run_scaling_campaign.py
  run run_scaling_campaign_v4.py --verify-reduction
  run ablation_v5.py
  run run_baseline_matched_v6.py
  run run_uncertainty_baseline_v6.py
  run run_power_rts_unit_specific_v7.py   # unreduced search in 30 cells + decomposition: about 1 h on two cores
  run rts_rtarget_sensitivity_v9.py
  run rts_unit_specific_matched_v9.py
  run rts_true_frontier_v9.py
fi
for s in verify_fyffe_original.py verify_unit_specific_bruteforce_v8.py verify_fyffe_provenance_v6.py run_measure_family.py run_audits_v6.py \
         run_power_rts_campaign_v7.py run_delta_sensitivity_v7.py generate_figures_v6.py \
         plot_parametric_partition_v7.py policy_menu_v8.py crisp_vs_fuzzy_v8.py \
         plot_true_frontier_v9.py replication_table_v10.py baselines_v11.py nesting_rts_example_v11.py supp_trace_v11.py; do run "$s"; done
(cd manuscript && for d in manuscript supplementary; do latexmk -pdf -interaction=nonstopmode "$d.tex"; done)
echo "Done."
