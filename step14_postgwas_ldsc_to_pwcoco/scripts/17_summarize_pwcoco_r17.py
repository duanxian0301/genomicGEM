from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE = Path(r"D:\文章\GS\postgwas\08_pwcoco")
RUNS = BASE / "runs" / "R17"

CONFIGS = [
    ("F1_AD", RUNS / "F1_AD" / "F1_AD_R17_final.coloc"),
    ("F2_AD", RUNS / "F2_AD" / "F2_AD_R17.coloc"),
    ("F1_PD", RUNS / "F1_PD" / "F1_PD_R17.coloc"),
    ("F2_PD", RUNS / "F2_PD" / "F2_PD_R17.coloc"),
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
    for pair, path in CONFIGS:
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "pair", pair)
        df["result_type"] = df.apply(classify, axis=1)
        frames.append(df)

    out = pd.concat(frames, ignore_index=True)
    out.to_csv(BASE / "pwcoco_R17_summary.tsv", sep="\t", index=False)

    best = (
        out.sort_values(["pair", "H4"], ascending=[True, False])
        .groupby("pair", as_index=False)
        .head(1)
        .reset_index(drop=True)
    )
    best.to_csv(BASE / "pwcoco_R17_best_h4.tsv", sep="\t", index=False)
    print(BASE / "pwcoco_R17_summary.tsv")
    print(BASE / "pwcoco_R17_best_h4.tsv")


if __name__ == "__main__":
    main()
