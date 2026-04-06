from __future__ import annotations

from pathlib import Path

import pandas as pd


def main() -> None:
    out = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\batch_status\brainmeta_batch_launch.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for trait in ["AD", "PD", "LBD"]:
        for chrom in range(1, 7):
            rows.append(
                {
                    "trait": trait,
                    "chromosome": chrom,
                    "mode": "background",
                    "launch_batch": "initial_ndd_chr1_6",
                    "status_note": "launched_from_powershell_background",
                }
            )

    pd.DataFrame(rows).to_csv(out, sep="\t", index=False)
    print(out)


if __name__ == "__main__":
    main()
