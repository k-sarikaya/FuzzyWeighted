# Reproduction driver (Windows PowerShell).
#   .\reproduce.ps1            unit tests + fast artifacts + figures + manuscript
#   .\reproduce.ps1 -Full      also re-solve every grid cell and the scaling family (hours)
param([switch]$Full, [switch]$SkipLatex)

$ErrorActionPreference = "Stop"
$root = (Resolve-Path $PSScriptRoot).Path
$python = (Get-Command python -ErrorAction SilentlyContinue).Source
if (-not $python) { throw "Python interpreter not found on PATH." }
$env:PYTHONPATH = Join-Path $root "src"
Set-Location $root

function Run($script, [string[]]$scriptArgs = @()) {
    Write-Output "`n>>> $script $scriptArgs"
    & $python (Join-Path "experiments" $script) @scriptArgs
    if ($LASTEXITCODE -ne 0) { throw "$script failed" }
}

Write-Output "[1] unit tests"
& $python -m unittest discover -s tests -v
if ($LASTEXITCODE -ne 0) { throw "Tests failed" }

if ($Full) {
    Write-Output "[2] exact solves (long)"
    Run "verify_pharma_oracle.py"
    Run "verify_fyffe_dp.py"
    Run "run_campaign_v2.py"
    Run "generate_benchmark_tables_v6.py"
    Run "parametric_stability.py" @("--cpp")
    Run "run_scaling_campaign.py"
    Run "run_scaling_campaign_v4.py" @("--verify-reduction")
    Run "ablation_v5.py"
    Run "run_baseline_matched_v6.py"
    Run "run_uncertainty_baseline_v6.py"
    Run "run_power_rts_unit_specific_v7.py"   # about 1 h on two cores
    Run "rts_rtarget_sensitivity_v9.py"
    Run "rts_unit_specific_matched_v9.py"
    Run "rts_true_frontier_v9.py"
}

Write-Output "[3] fast artifacts"
Run "verify_fyffe_original.py"
Run "verify_unit_specific_bruteforce_v8.py"
Run "verify_fyffe_provenance_v6.py"
Run "run_measure_family.py"
Run "run_audits_v6.py"
Run "run_power_rts_campaign_v7.py"
Run "run_delta_sensitivity_v7.py"
Run "generate_figures_v6.py"
Run "plot_parametric_partition_v7.py"
Run "policy_menu_v8.py"
Run "crisp_vs_fuzzy_v8.py"
Run "plot_true_frontier_v9.py"
Run "replication_table_v10.py"
Run "baselines_v11.py"
Run "nesting_rts_example_v11.py"
Run "supp_trace_v11.py"

if (-not $SkipLatex) {
    Write-Output "[4] manuscript"
    Set-Location (Join-Path $root "manuscript")
    foreach ($doc in "manuscript", "supplementary") {
        latexmk -pdf -interaction=nonstopmode "$doc.tex"
        if ($LASTEXITCODE -ne 0) { throw "LaTeX failed on $doc" }
    }
}
Write-Output "`nDone."
