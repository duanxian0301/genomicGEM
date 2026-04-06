from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


TRAIT_INPUTS = {
    "F1": Path(r"D:\文章\GS\GWAS\step11_factor_gwas_native_results\merged\standard_txt\ALPS_F1_factorGWAS_native_standard.txt"),
    "F2": Path(r"D:\文章\GS\GWAS\step11_factor_gwas_native_results\merged\standard_txt\ALPS_F2_factorGWAS_native_standard.txt"),
    "AD": Path(r"D:\文章\4NDD\NDDGWAS\AD.txt"),
    "PD": Path(r"D:\文章\4NDD\NDDGWAS\PD.txt"),
    "LBD": Path(r"D:\文章\4NDD\NDDGWAS\LBD.txt"),
}

TRAIT_FREQ_MODE = {
    "F1": "maf",
    "F2": "maf",
    "AD": "effect",
    "PD": "effect",
    "LBD": "effect",
}

SMR_COLUMNS = {
    "SNP": "SNP",
    "A1": "A1",
    "A2": "A2",
    "BETA": "b",
    "SE": "se",
    "P": "p",
    "N": "n",
}


def prepare_trait(trait: str, output_dir: Path) -> Path:
    input_path = TRAIT_INPUTS[trait]
    if not input_path.exists():
        raise FileNotFoundError(f"Input not found for {trait}: {input_path}")

    df = pd.read_csv(input_path, sep="\t", dtype={"CHR": str})
    missing = [col for col in SMR_COLUMNS if col not in df.columns]
    if missing:
        raise ValueError(f"{trait} is missing required columns: {missing}")
    if "FREQ" not in df.columns:
        raise ValueError(f"{trait} is missing required column: FREQ")

    out = df[list(SMR_COLUMNS)].rename(columns=SMR_COLUMNS).copy()
    freq = pd.to_numeric(df["FREQ"], errors="coerce")
    # Factor GWAS standard exports store MAF rather than effect-allele frequency, which is not
    # directly suitable for SMR frequency QC. We keep them on a separate branch and do not use
    # the same preparation logic as NDD traits.
    if TRAIT_FREQ_MODE[trait] == "effect":
        out["freq"] = freq
    elif TRAIT_FREQ_MODE[trait] == "maf":
        out["freq"] = freq
    else:
        raise ValueError(f"Unknown frequency mode for {trait}: {TRAIT_FREQ_MODE[trait]}")
    out = out.dropna(subset=["SNP", "A1", "A2", "freq", "b", "se", "p", "n"])
    out = out.loc[out["SNP"].astype(str).str.startswith("rs")].copy()
    out = out.loc[out["A1"].astype(str).str.len().between(1, 50)]
    out = out.loc[out["A2"].astype(str).str.len().between(1, 50)]
    out = out.drop_duplicates(subset=["SNP"], keep="first")
    out = out[["SNP", "A1", "A2", "freq", "b", "se", "p", "n"]]

    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{trait}.smr.txt"
    out.to_csv(output_path, sep="\t", index=False, float_format="%.10g")
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare GWAS summary statistics for SMR.")
    parser.add_argument(
        "--traits",
        nargs="+",
        default=["F1", "F2", "AD", "PD", "LBD"],
        choices=sorted(TRAIT_INPUTS),
        help="Traits to export in SMR format.",
    )
    parser.add_argument(
        "--output-dir",
        default=r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\input\smr_sumstats",
        help="Output directory for prepared SMR summary files.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    for trait in args.traits:
        output_path = prepare_trait(trait, output_dir)
        print(f"{trait}\t{output_path}")


if __name__ == "__main__":
    main()
