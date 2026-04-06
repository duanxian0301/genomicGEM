from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


BASE_OUT = Path(r"D:\文章\GS\postgwas\03_pleiofdr_all_pairs")
TXT_DIR = BASE_OUT / "inputs_txt"

TRAITS = {
    "F1": Path(r"D:\文章\GS\GWAS\step11_factor_gwas_native_results\merged\standard_txt\ALPS_F1_factorGWAS_native_standard.txt"),
    "F2": Path(r"D:\文章\GS\GWAS\step11_factor_gwas_native_results\merged\standard_txt\ALPS_F2_factorGWAS_native_standard.txt"),
    "AD": Path(r"D:\文章\4NDD\NDDGWAS\AD.txt"),
    "PD": Path(r"D:\文章\4NDD\NDDGWAS\PD.txt"),
    "LBD": Path(r"D:\文章\4NDD\NDDGWAS\LBD.txt"),
}


def prepare_trait(name: str, path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{name}: missing columns {missing}")

    out = df.loc[:, required].copy()
    out["Z"] = out["BETA"] / out["SE"]
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out[
        out["SNP"].notna()
        & out["CHR"].notna()
        & out["BP"].notna()
        & out["A1"].notna()
        & out["A2"].notna()
        & out["SE"].notna()
        & (out["SE"] > 0)
        & out["P"].notna()
        & (out["P"] > 0)
        & (out["P"] <= 1)
        & out["Z"].notna()
        & out["N"].notna()
    ].copy()

    out = out.rename(columns={"P": "PVAL"})
    out = out.loc[:, ["SNP", "CHR", "BP", "A1", "A2", "Z", "PVAL", "N"]]
    out["CHR"] = out["CHR"].astype(int)
    out["BP"] = out["BP"].astype(int)
    out = out.sort_values(["CHR", "BP", "SNP"]).drop_duplicates(subset=["SNP"], keep="first")
    return out


def main() -> None:
    TXT_DIR.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for trait, path in TRAITS.items():
        df = prepare_trait(trait, path)
        out_path = TXT_DIR / f"{trait}_fdr.txt"
        df.to_csv(out_path, sep="\t", index=False)
        manifest_rows.append(
            {
                "trait": trait,
                "source_file": str(path),
                "output_txt": str(out_path),
                "n_rows": len(df),
                "min_chr": int(df["CHR"].min()),
                "max_chr": int(df["CHR"].max()),
                "min_bp": int(df["BP"].min()),
                "max_bp": int(df["BP"].max()),
            }
        )

    pd.DataFrame(manifest_rows).to_csv(BASE_OUT / "input_manifest.tsv", sep="\t", index=False)
    print(f"Wrote inputs to {TXT_DIR}")


if __name__ == "__main__":
    main()
