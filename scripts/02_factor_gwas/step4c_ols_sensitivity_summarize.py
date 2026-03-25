from pathlib import Path

import pandas as pd


BASE = Path(r"D:\文章\GS\GWAS\ols_sensitivity_mini_compare")


def main() -> None:
    merged = pd.read_csv(BASE / "sumstats_false_vs_true_merged.tsv.gz", sep="\t")
    false_df = pd.read_csv(BASE / "ols_false" / "ols_false_sumstats.tsv.gz", sep="\t")
    true_df = pd.read_csv(BASE / "ols_true" / "ols_true_sumstats.tsv.gz", sep="\t")

    cols = set(merged.columns)
    metrics = [{"metric": "n_rows_merged", "value": len(merged)}]

    for cfalse, ctrue, label in [
        ("MAF_false", "MAF_true", "MAF"),
    ]:
        if cfalse in cols and ctrue in cols:
            x = merged[cfalse]
            y = merged[ctrue]
            metrics.append({"metric": f"cor_{label}", "value": x.corr(y)})
            metrics.append({"metric": f"mean_abs_delta_{label}", "value": (x - y).abs().mean()})
            metrics.append({"metric": f"median_abs_delta_{label}", "value": (x - y).abs().median()})

    trait_names = sorted(
        {
            c.replace("beta.", "").replace("_false", "")
            for c in cols
            if c.startswith("beta.") and c.endswith("_false")
        }
    )
    trait_rows = []
    for trait in trait_names:
        bf = merged[f"beta.{trait}_false"]
        bt = merged[f"beta.{trait}_true"]
        sf = merged[f"se.{trait}_false"]
        st = merged[f"se.{trait}_true"]
        trait_rows.append(
            {
                "trait": trait,
                "cor_beta": bf.corr(bt),
                "mean_abs_delta_beta": (bf - bt).abs().mean(),
                "median_abs_delta_beta": (bf - bt).abs().median(),
                "mean_abs_ratio_beta": (bf.abs() / bt.abs().replace(0, pd.NA)).dropna().mean(),
                "cor_se": sf.corr(st),
                "mean_abs_delta_se": (sf - st).abs().mean(),
                "median_abs_delta_se": (sf - st).abs().median(),
                "mean_ratio_se_false_over_true": (sf / st.replace(0, pd.NA)).dropna().mean(),
            }
        )

    false_n = len(false_df)
    true_n = len(true_df)
    metrics.extend(
        [
            {"metric": "n_rows_ols_false", "value": false_n},
            {"metric": "n_rows_ols_true", "value": true_n},
            {"metric": "overlap_ratio_false", "value": len(merged) / false_n},
            {"metric": "overlap_ratio_true", "value": len(merged) / true_n},
        ]
    )

    out = pd.DataFrame(metrics)
    out.to_csv(BASE / "sumstats_false_vs_true_summary.tsv", sep="\t", index=False)
    pd.DataFrame(trait_rows).to_csv(BASE / "sumstats_false_vs_true_by_trait.tsv", sep="\t", index=False)
    print(out.to_string(index=False))
    print(pd.DataFrame(trait_rows).to_string(index=False))


if __name__ == "__main__":
    main()
