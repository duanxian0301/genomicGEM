from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


SMR_EXE = Path(r"D:\SMR\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe")
LD_BFILE = Path(r"D:\SMR\g1000\g1000_eur")
SUMSTAT_DIR = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\input\smr_sumstats")
GTEX_BRAIN_DIR = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\refs\smr_gtex_v8_brain_hg19\brain_lite\eQTL_besd_lite")
BRAINMETA_DIR = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\refs\smr_brainmeta_v2_hg19_complete\BrainMeta_cis_eqtl_summary")
OUT_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd")

GTEX_BRAIN_TISSUES = [
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra",
]

TRAIT_FREQ_MODE = {
    "F1": "limited_no_freq_qc",
    "F2": "limited_no_freq_qc",
    "AD": "standard",
    "PD": "standard",
    "LBD": "standard",
}


def run_command(cmd: list[str], log_path: Path) -> int:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as log_f:
        proc = subprocess.run(cmd, stdout=log_f, stderr=subprocess.STDOUT, text=True)
    return proc.returncode


def gtex_prefix(tissue: str) -> Path:
    return GTEX_BRAIN_DIR / f"{tissue}.lite"


def brainmeta_prefix(chrom: int) -> Path:
    return BRAINMETA_DIR / f"BrainMeta_cis_eQTL_chr{chrom}"


def build_cmd(trait: str, gwas_file: Path, beqtl_prefix: Path, out_prefix: Path) -> list[str]:
    cmd = [
        str(SMR_EXE),
        "--bfile",
        str(LD_BFILE),
        "--gwas-summary",
        str(gwas_file),
        "--beqtl-summary",
        str(beqtl_prefix),
        "--out",
        str(out_prefix),
    ]
    if TRAIT_FREQ_MODE[trait] == "limited_no_freq_qc":
        cmd.append("--disable-freq-ck")
    return cmd


def main() -> None:
    parser = argparse.ArgumentParser(description="Run SMR for a single trait and reference panel.")
    parser.add_argument("--trait", required=True, choices=sorted(TRAIT_FREQ_MODE))
    parser.add_argument("--panel", required=True, choices=["gtex", "brainmeta"])
    parser.add_argument("--tissue", help="GTEx tissue name")
    parser.add_argument("--chrom", type=int, help="BrainMeta chromosome")
    args = parser.parse_args()

    gwas_file = SUMSTAT_DIR / f"{args.trait}.smr.txt"
    if not gwas_file.exists():
        raise FileNotFoundError(f"GWAS sumstats not found: {gwas_file}")

    if args.panel == "gtex":
        if args.tissue not in GTEX_BRAIN_TISSUES:
            raise ValueError(f"Invalid or missing --tissue for GTEx panel: {args.tissue}")
        ref_prefix = gtex_prefix(args.tissue)
        out_dir = OUT_ROOT / "single_trait" / args.trait / "gtex_v8_brain"
        out_prefix = out_dir / f"{args.trait}_{args.tissue}"
    else:
        if args.chrom is None or not (1 <= args.chrom <= 22):
            raise ValueError("BrainMeta runs require --chrom 1..22")
        ref_prefix = brainmeta_prefix(args.chrom)
        out_dir = OUT_ROOT / "single_trait" / args.trait / "brainmeta_v2_cortex"
        out_prefix = out_dir / f"{args.trait}_BrainMeta_chr{args.chrom}"

    if not (
        Path(str(ref_prefix) + ".besd").exists()
        and Path(str(ref_prefix) + ".epi").exists()
        and Path(str(ref_prefix) + ".esi").exists()
    ):
        raise FileNotFoundError(f"Missing BESD triplet for reference prefix: {ref_prefix}")

    cmd = build_cmd(args.trait, gwas_file, ref_prefix, out_prefix)
    log_path = Path(str(out_prefix) + ".run.log")
    code = run_command(cmd, log_path)
    print(f"exit_code\t{code}")
    print(f"log\t{log_path}")
    print(f"out_prefix\t{out_prefix}")


if __name__ == "__main__":
    main()
