# Stage 2 (R6) screening runner

# Reproducible launch of the four-arm short training. Run from D:\MetaSieve
# with the drug env:  D:\anaconda\envs\drug\python.exe
$ErrorActionPreference = "Stop"
$env:PYTHONIOENCODING = "utf-8"

$root = "D:\MetaSieve\report\meta_fewshot\stageR6_reltransport_screening_20260816"
$split = "D:\MetaSieve\dataset\processed\meta_fewshot\bindingdb_ki_double_cold_v1"
$python = "D:\anaconda\envs\drug\python.exe"
$seeds = @(20260815, 20260816, 20260817)

foreach ($seed in $seeds) {
    # A0: incumbent similarity_only grammar, R3/R4 recipe at short budget
    $out = "$root\A0_incumbent_seed$seed"
    if (-not (Test-Path "$out\RESULT.json")) {
        & $python -m scripts.train_qpsmp --arch similarity_only --steps 300 `
            --episodes-per-step 3 --learning-rate 6e-4 --lr-schedule cosine `
            --val-interval 50 --seed $seed --split-directory $split --output $out
    }
    # A1: new architecture, ordinary training
    $out = "$root\A1_ordinary_seed$seed"
    if (-not (Test-Path "$out\RESULT.json")) {
        & $python -m scripts.train_reltransport --ordinary --steps 300 `
            --episodes-per-step 3 --val-interval 100 --seed $seed --split-directory $split --output $out
    }
    # A2: new architecture, full training method
    $out = "$root\A2_full_seed$seed"
    if (-not (Test-Path "$out\RESULT.json")) {
        & $python -m scripts.train_reltransport --steps 300 `
            --episodes-per-step 3 --val-interval 100 --seed $seed --split-directory $split --output $out
    }
    # A3: A2 without the delta gate
    $out = "$root\A3_nogate_seed$seed"
    if (-not (Test-Path "$out\RESULT.json")) {
        & $python -m scripts.train_reltransport --no-gate --steps 300 `
            --episodes-per-step 3 --val-interval 100 --seed $seed --split-directory $split --output $out
    }
    # A4: A2 without counterfactual contrasts
    $out = "$root\A4_nocounterfactual_seed$seed"
    if (-not (Test-Path "$out\RESULT.json")) {
        & $python -m scripts.train_reltransport --no-counterfactual --steps 300 `
            --episodes-per-step 3 --val-interval 100 --seed $seed --split-directory $split --output $out
    }
}

# Paired comparison on the identical bank (after all arms finish)
& $python -m scripts.stageR6_compare_arms `
    --split-directory $split `
    --arm A0=$root\A0_incumbent_seed20260815\checkpoint.pt `
    --arm A0=$root\A0_incumbent_seed20260816\checkpoint.pt `
    --arm A0=$root\A0_incumbent_seed20260817\checkpoint.pt `
    --arm A1=$root\A1_ordinary_seed20260815\checkpoint.pt `
    --arm A1=$root\A1_ordinary_seed20260816\checkpoint.pt `
    --arm A1=$root\A1_ordinary_seed20260817\checkpoint.pt `
    --arm A2=$root\A2_full_seed20260815\checkpoint.pt `
    --arm A2=$root\A2_full_seed20260816\checkpoint.pt `
    --arm A2=$root\A2_full_seed20260817\checkpoint.pt `
    --arm A3=$root\A3_nogate_seed20260815\checkpoint.pt `
    --arm A3=$root\A3_nogate_seed20260816\checkpoint.pt `
    --arm A3=$root\A3_nogate_seed20260817\checkpoint.pt `
    --arm A4=$root\A4_nocounterfactual_seed20260815\checkpoint.pt `
    --arm A4=$root\A4_nocounterfactual_seed20260816\checkpoint.pt `
    --arm A4=$root\A4_nocounterfactual_seed20260817\checkpoint.pt `
    --reference A0 --output $root\COMPARE_R6_meta_val.json
