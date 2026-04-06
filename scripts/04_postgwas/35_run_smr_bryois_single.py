from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


SMR_EXE = Path(r"D:\SMR\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe")
LD_BFILE = Path(r"D:\SMR\g1000\g1000_eur")
SUMSTAT_DIR = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\input\smr_sumstats")
BRYOIS_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\refs\bryois2022_celltype_eqtl\prepared")
OUT_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\single_trait")

TRAIT_FREQ_MODE = {
    "F1": "limited_no_freq_qc",
    "F2": "limited_no_freq_qc",
    "AD": "standard",
    "PD": "standard",
    "LBD": "standard",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run SMR for a single trait against one Bryois cell-type chromosome reference.")
    parser.add_argument("--trait", required=True, choices=sorted(TRAIT_FREQ_MODE))
    parser.add_argument("--celltype", required=True)
    parser.add_argument("--chrom", required=True, type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    gwas_file = SUMSTAT_DIR / f"{args.trait}.smr.txt"
    ref_prefix = BRYOIS_ROOT / args.celltype / f"chr{args.chrom}" / f"{args.celltype}_chr{args.chrom}_filtered"
    out_dir = OUT_ROOT / args.trait / "bryois2022_celltype" / args.celltype
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / f"{args.trait}_{args.celltype}_chr{args.chrom}"
    log_path = Path(str(out_prefix) + ".run.log")

    cmd = [
        str(SMR_EXE),
        "--bfile",
        str(LD_BFILE),
        "--gwas-summary",
        str(gwas_file),
        "--beqtl-summary",
        str(ref_prefix),
        "--out",
        str(out_prefix),
    ]
    if TRAIT_FREQ_MODE[args.trait] == "limited_no_freq_qc":
        cmd.append("--disable-freq-ck")

    with log_path.open("w", encoding="utf-8") as log_f:
        proc = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, text=True)
    print(f"exit_code\t{proc.returncode}")
    print(f"out_prefix\t{out_prefix}")


if __name__ == "__main__":
    main()
