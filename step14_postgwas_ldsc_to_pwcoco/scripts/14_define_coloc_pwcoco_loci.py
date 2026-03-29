from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE_RESULTS = Path(r"D:\文章\GS\postgwas\03_pleiofdr_all_pairs\results")
OUT_DIR = Path(r"D:\文章\GS\postgwas\06_coloc_pwcoco_candidate_loci")
PAIRS = ["F1_AD", "F2_AD", "F1_PD", "F2_PD"]
WINDOW_BP = 500_000


def load_loci() -> pd.DataFrame:
    frames = []
    for pair in PAIRS:
        trait, disease = pair.split("_")
        df = pd.read_csv(BASE_RESULTS / pair / f"{pair}_conjfdr_0.05_loci.csv")
        df["pair"] = pair
        df["trait"] = trait
        df["disease"] = disease
        frames.append(df.loc[:, ["pair", "trait", "disease", "locusnum", "snpid", "chrnum", "chrpos", "min_conjfdr"]])
    return pd.concat(frames, ignore_index=True)


def build_sentinel_table(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby(["snpid", "chrnum", "chrpos"], as_index=False)
        .agg(
            pair_count=("pair", "nunique"),
            pairs=("pair", lambda x: ";".join(sorted(set(x)))),
            traits=("trait", lambda x: ";".join(sorted(set(x)))),
            diseases=("disease", lambda x: ";".join(sorted(set(x)))),
            best_conjfdr=("min_conjfdr", "min"),
        )
        .sort_values(["chrnum", "chrpos", "best_conjfdr"], ignore_index=True)
    )
    grouped.insert(0, "sentinel_id", [f"S{i:02d}" for i in range(1, len(grouped) + 1)])
    grouped["region_start_500kb"] = (grouped["chrpos"] - WINDOW_BP).clip(lower=1).astype(int)
    grouped["region_end_500kb"] = (grouped["chrpos"] + WINDOW_BP).astype(int)
    return grouped


def build_region_table(sentinels: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for chrnum, sub in sentinels.groupby("chrnum", sort=True):
        sub = sub.sort_values("chrpos").reset_index(drop=True)
        current = None
        for _, row in sub.iterrows():
            start = int(row["region_start_500kb"])
            end = int(row["region_end_500kb"])
            if current is None:
                current = {
                    "chrnum": int(chrnum),
                    "region_start": start,
                    "region_end": end,
                    "sentinel_ids": [row["sentinel_id"]],
                    "sentinel_snps": [row["snpid"]],
                    "pairs": set(row["pairs"].split(";")),
                    "traits": set(row["traits"].split(";")),
                    "diseases": set(row["diseases"].split(";")),
                    "best_conjfdr": float(row["best_conjfdr"]),
                }
            elif start <= current["region_end"]:
                current["region_end"] = max(current["region_end"], end)
                current["sentinel_ids"].append(row["sentinel_id"])
                current["sentinel_snps"].append(row["snpid"])
                current["pairs"].update(row["pairs"].split(";"))
                current["traits"].update(row["traits"].split(";"))
                current["diseases"].update(row["diseases"].split(";"))
                current["best_conjfdr"] = min(current["best_conjfdr"], float(row["best_conjfdr"]))
            else:
                rows.append(current)
                current = {
                    "chrnum": int(chrnum),
                    "region_start": start,
                    "region_end": end,
                    "sentinel_ids": [row["sentinel_id"]],
                    "sentinel_snps": [row["snpid"]],
                    "pairs": set(row["pairs"].split(";")),
                    "traits": set(row["traits"].split(";")),
                    "diseases": set(row["diseases"].split(";")),
                    "best_conjfdr": float(row["best_conjfdr"]),
                }
        if current is not None:
            rows.append(current)

    out = pd.DataFrame(rows)
    out.insert(0, "region_id", [f"R{i:02d}" for i in range(1, len(out) + 1)])
    out["n_sentinels"] = out["sentinel_ids"].apply(len)
    out["sentinel_ids"] = out["sentinel_ids"].apply(lambda x: ";".join(x))
    out["sentinel_snps"] = out["sentinel_snps"].apply(lambda x: ";".join(x))
    out["pairs"] = out["pairs"].apply(lambda x: ";".join(sorted(x)))
    out["traits"] = out["traits"].apply(lambda x: ";".join(sorted(x)))
    out["diseases"] = out["diseases"].apply(lambda x: ";".join(sorted(x)))
    out["region_width_bp"] = out["region_end"] - out["region_start"] + 1
    out["priority"] = out["n_sentinels"].apply(lambda n: "PWCoCo_priority" if n >= 2 else "coloc_first")
    out["notes"] = out["n_sentinels"].apply(
        lambda n: "Multiple nearby sentinel signals; prioritize PWCoCo after coloc." if n >= 2 else "Single sentinel region; classical coloc is a reasonable first pass."
    )
    return out.loc[:, ["region_id", "chrnum", "region_start", "region_end", "region_width_bp", "n_sentinels", "sentinel_ids", "sentinel_snps", "pairs", "traits", "diseases", "best_conjfdr", "priority", "notes"]]


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    raw = load_loci()
    raw.to_csv(OUT_DIR / "conjfdr_factor_loci_raw_23rows.tsv", sep="\t", index=False)

    sentinel = build_sentinel_table(raw)
    sentinel.to_csv(OUT_DIR / "conjfdr_factor_loci_sentinel_21.tsv", sep="\t", index=False)

    region = build_region_table(sentinel)
    region.to_csv(OUT_DIR / "conjfdr_factor_loci_regions_500kb.tsv", sep="\t", index=False)

    pair_rows = []
    for _, row in region.iterrows():
        for pair in str(row["pairs"]).split(";"):
            trait, disease = pair.split("_")
            pair_rows.append(
                {
                    "region_id": row["region_id"],
                    "pair": pair,
                    "trait": trait,
                    "disease": disease,
                    "chrnum": row["chrnum"],
                    "region_start": row["region_start"],
                    "region_end": row["region_end"],
                    "n_sentinels": row["n_sentinels"],
                    "sentinel_snps": row["sentinel_snps"],
                    "priority": row["priority"],
                }
            )
    pair_df = pd.DataFrame(pair_rows).sort_values(["disease", "trait", "chrnum", "region_start"]).reset_index(drop=True)
    pair_df.to_csv(OUT_DIR / "conjfdr_factor_loci_region_pair_tasks.tsv", sep="\t", index=False)

    summary_lines = [
        "# Candidate loci for coloc and PWCoCo",
        "",
        f"- Raw conjFDR locus rows across F1/F2 x AD/PD: {len(raw)}",
        f"- Deduplicated sentinel loci: {len(sentinel)}",
        f"- Region-level merged loci using +/-500 kb windows: {len(region)}",
        "",
        "Interpretation:",
        "- Use the sentinel table for reporting unique sentinel hits.",
        "- Use the region table for extraction of locus-level summary statistics for coloc/PWCoCo.",
        "- Use the region-pair task table to run coloc/PWCoCo per trait-disease pair.",
        "- Regions marked `PWCoCo_priority` contain multiple nearby sentinel signals and are the best first candidates for conditional colocalization.",
    ]
    (OUT_DIR / "README.md").write_text("\n".join(summary_lines), encoding="utf-8")
    print(f"Wrote locus files to: {OUT_DIR}")


if __name__ == "__main__":
    main()
