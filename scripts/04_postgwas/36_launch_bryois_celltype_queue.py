from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import psutil
import requests


PYTHON = sys.executable
PREP_REF = Path(r"D:\codex\GenomicSEM\scripts\04_postgwas\34_prepare_bryois_reference.py")
RUN_SMR = Path(r"D:\codex\GenomicSEM\scripts\04_postgwas\35_run_smr_bryois_single.py")
BRYOIS_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\refs\bryois2022_celltype_eqtl")
OUT_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\single_trait")
TRAITS = ["AD", "PD", "LBD", "F1", "F2"]
CHROMS = list(range(1, 23))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Launch Bryois cell-type SMR queue with bounded concurrency.")
    parser.add_argument("--celltype", required=True)
    parser.add_argument("--max-total-smr", type=int, default=8)
    return parser.parse_args()


def count_active_smr() -> int:
    total = 0
    for proc in psutil.process_iter(["name"]):
        try:
            name = (proc.info.get("name") or "").lower()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
        if "smr" in name:
            total += 1
    return total


def ensure_download(celltype: str, chrom: int) -> Path:
    download_dir = BRYOIS_ROOT / "downloads" / celltype
    download_dir.mkdir(parents=True, exist_ok=True)
    path = download_dir / f"{celltype}.{chrom}.gz"
    if not path.exists():
        url = f"https://zenodo.org/records/7276971/files/{celltype}.{chrom}.gz?download=1"
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with path.open("wb") as fout:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fout.write(chunk)
    snp_pos = BRYOIS_ROOT / "downloads" / "snp_pos.txt.gz"
    if not snp_pos.exists():
        url = "https://zenodo.org/records/7276971/files/snp_pos.txt.gz?download=1"
        with requests.get(url, stream=True, timeout=120) as r:
            r.raise_for_status()
            with snp_pos.open("wb") as fout:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:
                        fout.write(chunk)
    return path


def ensure_reference(celltype: str, chrom: int) -> None:
    prefix = BRYOIS_ROOT / "prepared" / celltype / f"chr{chrom}" / f"{celltype}_chr{chrom}_filtered.besd"
    if prefix.exists():
        return
    subprocess.run([PYTHON, str(PREP_REF), "--celltype", celltype, "--chrom", str(chrom)], check=True)


def output_exists(celltype: str, trait: str, chrom: int) -> bool:
    out = OUT_ROOT / trait / "bryois2022_celltype" / celltype / f"{trait}_{celltype}_chr{chrom}.smr"
    return out.exists() and out.stat().st_size > 0


def launch_one(celltype: str, trait: str, chrom: int) -> None:
    subprocess.Popen(
        [PYTHON, str(RUN_SMR), "--trait", trait, "--celltype", celltype, "--chrom", str(chrom)],
        creationflags=subprocess.CREATE_NEW_PROCESS_GROUP,
    )


def main() -> None:
    args = parse_args()
    launched = 0
    for chrom in CHROMS:
        ensure_download(args.celltype, chrom)
        ensure_reference(args.celltype, chrom)
        for trait in TRAITS:
            if output_exists(args.celltype, trait, chrom):
                continue
            while count_active_smr() >= args.max_total_smr:
                time.sleep(20)
            launch_one(args.celltype, trait, chrom)
            launched += 1
            time.sleep(2)
    print(f"celltype\t{args.celltype}")
    print(f"launched_total\t{launched}")


if __name__ == "__main__":
    main()
