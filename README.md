# Dual-fuzzy weighted k-out-of-n:G systems

Code, data and generated results for the manuscript *Reliability Optimization of
Weighted k-out-of-n:G Systems under Dual Fuzziness* (prepared for *Fuzzy Optimization and
Decision Making*). Every number, table and figure in `manuscript/` is produced by a script in
`experiments/` from the model in `src/`. The repository contains only what the article and its
supplement use.

Author: Kadir Sarıkaya, Department of Industrial Engineering, Tokat
Gaziosmanpaşa University.

## Scope

Component capacities and the success threshold are triangular fuzzy numbers,
evaluated through an α-cut level and a Hurwicz optimism index γ. For each
(α, γ) cell an exact branch-and-bound solver returns the design of minimum
feasible cost and, among those, maximum reliability. The benchmarks comprise:
1. **IEEE RTS-79 Power Generation Fleet (`power_rts_adapted`):** An engineering case study
   parameterised from the IEEE RTS-79 benchmark with $n=10$ generation bays, MW dispatch capacities, peak load forecast threshold $\tilde{k} = (162, 180, 198)$~MW,
   target availability $R_{\text{target}}=0.95$, availability probabilities $p_i \in [0.90, 0.99]$
   calculated from published IEEE RTS-79 MTTF/MTTR outage statistics, and an assumed ambient thermal derating (RTS-79 reports no derating).
2. **Pharmaceutical-derived and Fyffe-derived instances:** catalogs derived from two published allocation problems; they are
   not the original optimization problems. Every transformation is documented in `docs/` and `results/generated/fyffe_provenance_v6.json`.

## Layout

| Path | Content |
|---|---|
| `src/fuzzy_wkon/` | model, reliability kernels, exact solver, C++ driver, uncertainty-theory oracle |
| `cpp/` | C++ exact backend (built automatically on first use by `cpp_backend.py`) |
| `tests/` | unit, property and regression tests, including checks that every number quoted in the article and the supplement matches the generated results (`pytest tests`) |
| `experiments/` | one script per artifact; see the table below |
| `results/generated/` | JSON/CSV artifacts written by the scripts |
| `figures/` | every figure of the article and the supplement (PDF) |
| `manuscript/` | `manuscript.tex` (article, Springer Nature template) and `supplementary.tex` with their PDFs; the `table_*.tex` files are written by the scripts below |
| `data/` | source tables and reconciliations (source PDFs are kept locally, not redistributed) |
| `docs/` | data provenance notes |

## Script → artifact map

Script and artifact names carry the manuscript version in which they were introduced (for example
`_v6`, `_v11`); they are the names listed in Supplementary Material S2, and every one of them is used by
the current article or supplement.

| Script | Artifact | Used for |
|---|---|---|
| `run_power_rts_campaign_v7.py` | `benchmark_power_rts.json`, `manuscript/table_power_rts.tex` | IEEE RTS-79 30-cell exact optimization campaign |
| `run_delta_sensitivity_v7.py` | `delta_sensitivity_v7.json`, `figures/fig_delta_sensitivity.*` | continuous spread sensitivity ($\delta \in [0.02, 0.20]$) and seasonal asymmetric derating |
| `run_power_rts_unit_specific_v7.py` | `rts_unit_specific_v7.json`, `manuscript/table_rts_unit_specific.tex` | unit-specific spreads on RTS-79: all 30 cells by the unreduced search (600 s limit), thermal/hydro decomposition, hydro and thermal spread sweeps; about 1 h on two cores (`--full` only) |
| `verify_unit_specific_bruteforce_v8.py` | `unit_specific_bruteforce_v8.json` | independent exhaustive check (exact rationals) of four unit-specific RTS-79 optima and of the $171M fleet outside the crisp library |
| `policy_menu_v8.py` | `policy_menu_v8.json`, `manuscript/table_policy_menu.tex` | design menu of the RTS-79 parametric solution (cost, R_N, R_Pi, area share), budget inversion and gamma intervals from Corollary 2 (Sections 5.1 and 6.1) |
| `crisp_vs_fuzzy_v8.py` | `crisp_vs_fuzzy_v8.json`, `manuscript/table_crisp_vs_fuzzy.tex` | crisp versus dual-fuzzy reliability of the seven RTS-79 designs; capacities only, demand only and both fuzzified, in the necessity and possibility cells (Section 5.2) |
| `verify_pharma_oracle.py` | `pharma_oracle_audit.json` | exact solver vs exhaustive oracle, 30 cells |
| `verify_fyffe_dp.py` | `fyffe_exact_audit.json` | exact Fyffe campaign, 30 cells |
| `run_campaign_v2.py` | `campaign_v2_default.json` | 60-cell grid campaign |
| `generate_benchmark_tables_v6.py` | `benchmark_tables_v6.json`, `manuscript/table_*.tex` | Tables 4–5 |
| `run_measure_family.py` | `measure_family.json` | possibility / necessity / credibility; R_N, R_Π columns |
| `parametric_stability.py` | `parametric_stability*.json` | complete parametric solution (22 and 31 crisp problems) |
| `run_scaling_campaign.py`, `run_scaling_campaign_v4.py` | `scaling_campaign*.json` | scaling family, n ≤ 20 |
| `ablation_v5.py` | `ablation_v5.json` | ablation table; `ablation_v6.json` adds two re-measured cells (see its `correction_note`) |
| `run_audits_v6.py` | `audits_v6.json` | every verification claim (closed forms, correspondence sweep, factorisation, credibility, node accounting) |
| `run_baseline_matched_v6.py` | `baseline_matched_v6.json` | budget-matched genetic-algorithm comparison |
| `run_uncertainty_baseline_v6.py` | `uncertainty_baseline_v6.json` | uncertainty-theory (Liu) oracle substitution |
| `verify_fyffe_original.py`, `verify_fyffe_provenance_v6.py` | `fyffe_original_reproduction.json`, `fyffe_provenance_v6.json` | reproduction of Fyffe, Hines and Lee (1968) |
| `generate_figures_v6.py` | `figures/fig2`–`fig5`, `figS_fyffe_sensitivity` | policy-effect and cost figures of the article and the supplement |
| `plot_parametric_partition_v7.py` | `figures/fig_parametric_partition.*` | complete parametric solution of RTS-79 over the policy square (checks Corollary 2 boundaries and all 30 grid cells) |
| `rts_true_frontier_v9.py` | `rts_true_frontier_v9.json`, `manuscript/table_adequate_frontier.tex` | exhaustive enumeration of the 59,049 RTS-79 fleets: cost--necessity frontier with and without crisp adequacy (Section 5.4) |
| `plot_true_frontier_v9.py` | `figures/fig_true_pareto_frontier.*` | figure of the cost--necessity frontier (Section 5.4) |
| `rts_rtarget_sensitivity_v9.py` | `rts_rtarget_sensitivity_v9.json`, `manuscript/table_rts_rtarget_sensitivity.tex` | optimal fleet for reliability targets 0.90--0.99 at T = 148, 180 and 220 MW (Section 5.5, Supplementary S6) |
| `rts_unit_specific_matched_v9.py` | `rts_unit_specific_matched_v9.json` | common, unit-specific and matched common spreads on RTS-79, exact rational brute force (Supplementary S6) |
| `replication_table_v10.py` | `replication_v10.json`, `manuscript/table_replication.tex` | the three instances side by side and the area measures of what the 30-cell grid misses (Section 5.1) |
| `baselines_v11.py` | `baselines_v11.json`, `manuscript/table_baselines*.tex` | model baselines: crisp, one-sided fuzzy and probabilistic counterpart, evaluated under the declared rule (Sections 4.2 and 5.3, Supplementary S7) |
| `nesting_rts_example_v11.py` | `nesting_rts_example_v11.json`, `manuscript/table_nesting_rts_supp.tex` | readings of the nesting table mapped to cells, thresholds and designs on the RTS-79 fleet (no optimization) |
| `supp_trace_v11.py` | `supp_trace_v11.json` | optimization trace of Supplementary Material S1: Algorithm 1 on the three-subsystem example (independent transcription, exact rationals) |

## Reproduce

```text
pip install -r requirements.txt
.\reproduce.ps1            # Windows: tests, fast artifacts, figures, manuscript
./reproduce.sh             # Linux/macOS
.\reproduce.ps1 -Full      # also re-solve every grid cell and the scaling family (hours)
```

A C++17 compiler (`g++` or `clang++`) must be on `PATH` for the exact backend.
Individual scripts run from the repository root with `PYTHONPATH=src`.
