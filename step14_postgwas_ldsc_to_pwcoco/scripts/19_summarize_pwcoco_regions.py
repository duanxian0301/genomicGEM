from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE = Path(r"D:\文章\GS\postgwas\08_pwcoco")

CONFIGS = [
    ("R01", "F1_AD", "sensitivity_p1e-5_p2_1e-4", BASE / "runs" / "R01" / "F1_AD_sensitivity" / "F1_AD_R01_p1e5_p2_1e4.coloc"),
    ("R03", "F1_PD", "sensitivity_p1=5e-5", BASE / "runs" / "R03" / "F1_PD_sensitivity2" / "F1_PD_R03_p5e5_force.coloc"),
    ("R06", "F1_PD", "sensitivity_p2=1e-5", BASE / "runs" / "R06" / "F1_PD_sensitivity" / "F1_PD_R06_p2e5.coloc"),
    ("R12", "F1_PD", "sensitivity_p1=1e-5", BASE / "runs" / "R12" / "F1_PD_sensitivity" / "F1_PD_R12_p1e5_force.coloc"),
    ("R15", "F1_PD", "sensitivity_p1=2e-5_p2=1e-5", BASE / "runs" / "R15" / "F1_PD_sensitivity2" / "F1_PD_R15_p2e5_p2_1e5.coloc"),
    ("R17", "F1_AD", "sensitivity_p1=1e-5", BASE / "runs" / "R17" / "F1_AD_sensitivity" / "F1_AD_R17_p1e5.coloc"),
    ("R17", "F2_AD", "primary", BASE / "runs" / "R17" / "F2_AD" / "F2_AD_R17.coloc"),
    ("R17", "F1_PD", "primary", BASE / "runs" / "R17" / "F1_PD" / "F1_PD_R17.coloc"),
    ("R17", "F2_PD", "primary", BASE / "runs" / "R17" / "F2_PD" / "F2_PD_R17.coloc"),
    ("R18", "F1_PD", "primary", BASE / "runs" / "R18" / "F1_PD" / "F1_PD_R18.coloc"),
    ("R18", "F2_PD", "primary", BASE / "runs" / "R18" / "F2_PD" / "F2_PD_R18.coloc"),
]


def classify(row: pd.Series) -> str:
    if row["SNP1"] == "unconditioned" and row["SNP2"] == "unconditioned":
        return "unconditioned"
    if row["SNP1"] != "unconditioned" and row["SNP2"] != "unconditioned":
        return "conditioned_both"
    if row["SNP1"] != "unconditioned":
        return "conditioned_trait"
    return "conditioned_disease"


def main() -> None:
    frames = []
    for region_id, pair, run_label, path in CONFIGS:
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "region_id", region_id)
        df.insert(1, "pair", pair)
        df.insert(2, "run_label", run_label)
        df["result_type"] = df.apply(classify, axis=1)
        frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(BASE / "pwcoco_region_summary.tsv", sep="\t", index=False)

    unconditioned = out[out["result_type"] == "unconditioned"].copy()
    conditioned = out[out["result_type"] != "unconditioned"].copy()

    best_uncond = (
        unconditioned.sort_values(["region_id", "pair", "H4"], ascending=[True, True, False])
        .groupby(["region_id", "pair"], as_index=False)
        .head(1)
        .reset_index(drop=True)
    )
    best_cond = (
        conditioned.sort_values(["region_id", "pair", "H4"], ascending=[True, True, False])
        .groupby(["region_id", "pair"], as_index=False)
        .head(1)
        .reset_index(drop=True)
    )

    merged = best_uncond.merge(
        best_cond[
            [
                "region_id",
                "pair",
                "run_label",
                "SNP1",
                "SNP2",
                "H3",
                "H4",
                "result_type",
            ]
        ],
        on=["region_id", "pair"],
        how="left",
        suffixes=("_uncond", "_best_cond"),
    )
    merged.to_csv(BASE / "pwcoco_region_best_summary.tsv", sep="\t", index=False)
    print(BASE / "pwcoco_region_summary.tsv")
    print(BASE / "pwcoco_region_best_summary.tsv")


if __name__ == "__main__":
    main()
