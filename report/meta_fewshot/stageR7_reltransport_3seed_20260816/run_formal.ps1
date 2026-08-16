# Stage 3 (R7) formal development runner

# Three seeds, full 1200-step budget, double-cold meta_val. A0 is the frozen
# Stage R3/R4 incumbent checkpoints (1200 steps, similarity_only); no A0
# training is rerun. meta_test stays sealed. Run from D:\MetaSieve.
$ErrorActionPreference = "Stop"
$env:PYTHONIOENCODING = "utf-8"

$root = "D:\MetaSieve\report\meta_fewshot\stageR7_reltransport_3seed_20260816"
$split = "D:\MetaSieve\dataset\processed\meta_fewshot\bindingdb_ki_double_cold_v1"
$python = "D:\anaconda\envs\drug\python.exe"
$seeds = @(20260815, 20260816, 20260817)
$a0root = "D:\MetaSieve\report\meta_fewshot\stageR3R4_level_shape_20260815"

foreach ($seed in $seeds) {
    $out = "$root\A1_ordinary_seed$seed"
    if (-not (Test-Path "$out\RESULT.json")) {
        & $python -m scripts.train_reltransport --ordinary --steps 1200 `
            --episodes-per-step 3 --val-interval 200 --seed $seed `
            --split-directory $split --output $out
    }
    $out = "$root\A2_full_seed$seed"
    if (-not (Test-Path "$out\RESULT.json")) {
        & $python -m scripts.train_reltransport --steps 1200 `
            --episodes-per-step 3 --val-interval 200 --seed $seed `
            --split-directory $split --output $out
    }
    $out = "$root\A3_nogate_seed$seed"
    if (-not (Test-Path "$out\RESULT.json")) {
        & $python -m scripts.train_reltransport --no-gate --steps 1200 `
            --episodes-per-step 3 --val-interval 200 --seed $seed `
            --split-directory $split --output $out
    }
}

# Paired comparison on the identical bank: candidate arms vs the frozen
# incumbent checkpoints, with the corrected donor contract.
$armArgs = @("--split-directory", $split)
foreach ($seed in $seeds) {
    $armArgs += @("--arm", "A0=$a0root\A0_incumbent_seed$seed\checkpoint.pt")
}
foreach ($seed in $seeds) {
    $armArgs += @("--arm", "A1=$root\A1_ordinary_seed$seed\checkpoint.pt")
}
foreach ($seed in $seeds) {
    $armArgs += @("--arm", "A2=$root\A2_full_seed$seed\checkpoint.pt")
}
foreach ($seed in $seeds) {
    $armArgs += @("--arm", "A3=$root\A3_nogate_seed$seed\checkpoint.pt")
}
$armArgs += @("--reference", "A0", "--output", "$root\COMPARE_R7_meta_val.json")
& $python -m scripts.stageR6_compare_arms @armArgs
