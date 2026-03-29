from __future__ import annotations

from pathlib import Path

import pandas as pd


BASE_OUT = Path(r"D:\文章\GS\postgwas\03b_pleiofdr_alps_family_vs_ndd")
RESULTS_DIR = BASE_OUT / "results"

ALPS_GROUP = {
    "aALPS": "original",
    "Left_ALPS": "original",
    "mALPS": "original",
    "pALPS": "original",
    "Right_ALPS": "original",
    "Mean_ALPS": "summary",
    "tALPS": "summary",
}


def count_rows(path: Path) -> int:
    if not path.exists():
        return 0
    try:
        df = pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return 0
    return int(len(df))


def main() -> None:
    rows = []
    for pair_dir in sorted(RESULTS_DIR.iterdir()):
        if not pair_dir.is_dir():
            continue
        trait, disease = pair_dir.name.rsplit("_", 1)
        prefix = f"{trait}_{disease}"
        conj_all = pair_dir / f"{prefix}_conjfdr_0.05_all.csv"
        conj_loci = pair_dir / f"{prefix}_conjfdr_0.05_loci.csv"
        z_all = pair_dir / f"{prefix}_zscore_conjfdr_0.05_all.csv"
        z_loci = pair_dir / f"{prefix}_zscore_conjfdr_0.05_loci.csv"
        rows.append(
            {
                "pair": pair_dir.name,
                "trait": trait,
                "trait_group": ALPS_GROUP.get(trait, "unknown"),
                "disease": disease,
                "conjfdr_all_n": count_rows(conj_all),
                "conjfdr_loci_n": count_rows(conj_loci),
                "zscore_all_n": count_rows(z_all),
                "zscore_loci_n": count_rows(z_loci),
                "all_file": str(conj_all),
                "loci_file": str(conj_loci),
            }
        )

    df = pd.DataFrame(rows).sort_values(["disease", "trait_group", "trait"])
    df.to_csv(BASE_OUT / "pleiofdr_alps_family_summary.tsv", sep="\t", index=False)

    notes = [
        "# pleioFDR ALPS-family vs NDD summary",
        "",
        "This table summarizes conjFDR<0.05 results for original ALPS traits and summary ALPS traits against AD, PD, and LBD.",
        "Counts are based on actual rows in the exported CSV files.",
    ]
    (BASE_OUT / "pleiofdr_alps_family_notes.md").write_text("\n".join(notes), encoding="utf-8")
    print(f"Wrote summary to {BASE_OUT}")


if __name__ == "__main__":
    main()
