from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm


BASE_OUT = Path(r"D:\文章\GS\postgwas\03b_pleiofdr_alps_family_vs_ndd")
TXT_DIR = BASE_OUT / "inputs_txt"
REF_PATH = Path(r"D:\pleioFDR\pleiofdr-master\9545380.ref")

TRAITS_BETA = {
    "aALPS": Path(r"D:\文章\GS\GWAS\aALPS.txt"),
    "Left_ALPS": Path(r"D:\文章\GS\GWAS\Left_ALPS.txt"),
    "mALPS": Path(r"D:\文章\GS\GWAS\mALPS.txt"),
    "pALPS": Path(r"D:\文章\GS\GWAS\pALPS.txt"),
    "Right_ALPS": Path(r"D:\文章\GS\GWAS\Right_ALPS.txt"),
}

TRAITS_Z_ONLY = {
    "Mean_ALPS": Path(r"D:\文章\GS\GWAS\step10_ldsc_validation\sumstats_inputs\Mean_ALPS.sumstats.gz"),
    "tALPS": Path(r"D:\文章\GS\GWAS\step10_ldsc_validation\sumstats_inputs\tALPS.sumstats.gz"),
}

NDD_TRAITS = {
    "AD": Path(r"D:\文章\4NDD\NDDGWAS\AD.txt"),
    "PD": Path(r"D:\文章\4NDD\NDDGWAS\PD.txt"),
    "LBD": Path(r"D:\文章\4NDD\NDDGWAS\LBD.txt"),
}


def load_ref() -> pd.DataFrame:
    ref = pd.read_csv(REF_PATH, sep=r"\s+", usecols=["SNP", "CHR", "BP"])
    ref = ref.drop_duplicates(subset=["SNP"], keep="first")
    return ref


def finalize_frame(df: pd.DataFrame) -> pd.DataFrame:
    out = df.replace([np.inf, -np.inf], np.nan).copy()
    out = out[
        out["SNP"].notna()
        & out["CHR"].notna()
        & out["BP"].notna()
        & out["A1"].notna()
        & out["A2"].notna()
        & out["Z"].notna()
        & out["PVAL"].notna()
        & (out["PVAL"] > 0)
        & (out["PVAL"] <= 1)
        & out["N"].notna()
    ].copy()
    out["CHR"] = out["CHR"].astype(int)
    out["BP"] = out["BP"].astype(int)
    out = out.loc[:, ["SNP", "CHR", "BP", "A1", "A2", "Z", "PVAL", "N"]]
    out = out.sort_values(["CHR", "BP", "SNP"]).drop_duplicates(subset=["SNP"], keep="first")
    return out


def prepare_beta_trait(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    out = df.loc[:, ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]].copy()
    out["Z"] = out["BETA"] / out["SE"]
    out = out.rename(columns={"P": "PVAL"})
    return finalize_frame(out)


def prepare_z_only_trait(path: Path, ref: pd.DataFrame) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", compression="infer")
    out = df.loc[:, ["SNP", "A1", "A2", "Z", "N"]].copy()
    out["PVAL"] = 2 * norm.sf(np.abs(out["Z"].astype(float)))
    out["PVAL"] = np.maximum(out["PVAL"], np.finfo(float).tiny)
    out = out.merge(ref, on="SNP", how="inner")
    return finalize_frame(out)


def prepare_ndd_trait(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    out = df.loc[:, ["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"]].copy()
    out["Z"] = out["BETA"] / out["SE"]
    out = out.rename(columns={"P": "PVAL"})
    return finalize_frame(out)


def main() -> None:
    TXT_DIR.mkdir(parents=True, exist_ok=True)
    ref = load_ref()

    manifest_rows = []

    for trait, path in TRAITS_BETA.items():
        df = prepare_beta_trait(path)
        out_path = TXT_DIR / f"{trait}_fdr.txt"
        df.to_csv(out_path, sep="\t", index=False)
        manifest_rows.append({"trait": trait, "group": "original", "source_file": str(path), "output_txt": str(out_path), "n_rows": len(df)})

    for trait, path in TRAITS_Z_ONLY.items():
        df = prepare_z_only_trait(path, ref)
        out_path = TXT_DIR / f"{trait}_fdr.txt"
        df.to_csv(out_path, sep="\t", index=False)
        manifest_rows.append({"trait": trait, "group": "summary", "source_file": str(path), "output_txt": str(out_path), "n_rows": len(df)})

    for trait, path in NDD_TRAITS.items():
        df = prepare_ndd_trait(path)
        out_path = TXT_DIR / f"{trait}_fdr.txt"
        df.to_csv(out_path, sep="\t", index=False)
        manifest_rows.append({"trait": trait, "group": "ndd", "source_file": str(path), "output_txt": str(out_path), "n_rows": len(df)})

    pd.DataFrame(manifest_rows).to_csv(BASE_OUT / "input_manifest.tsv", sep="\t", index=False)
    print(f"Wrote inputs to {TXT_DIR}")


if __name__ == "__main__":
    main()
