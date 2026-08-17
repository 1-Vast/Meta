import json, sys
for name in sys.argv[1:]:
    d = json.load(open(f"tools/research/stageF_pairwise/{name}_meta_val.rows.summary.json"))
    print(name + ":")
    for k in ("0", "1", "2", "3", "5"):
        c = d["aggregates"][k]["correct"]
        print("  k=" + k, {f: round(c[f], 4) for f in
              ("mse_pk", "level_squared", "centered_mse_pk", "spearman",
               "pearson", "ci", "cliff_sign")})
