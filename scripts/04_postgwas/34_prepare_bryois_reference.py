from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
import csv


PYTHON = sys.executable
PREP_FASTQTL = Path(r"D:\codex\GenomicSEM\scripts\04_postgwas\32_prepare_bryois_celltype_fastqtl.py")
SMR_EXE = Path(r"D:\SMR\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64\smr-1.3.1-win.exe")
BRYOIS_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\refs\bryois2022_celltype_eqtl")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare Bryois cell-type eQTL BESD reference for one cell type and chromosome.")
    parser.add_argument("--celltype", required=True)
    parser.add_argument("--chrom", required=True, type=int)
    return parser.parse_args()


def run(cmd: list[str]) -> None:
    subprocess.run(cmd, check=True)


def has_invalid_bp(path: Path, kind: str) -> bool:
    if not path.exists():
        return True
    bp_index = 3
    with path.open("r", encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if len(row) <= bp_index:
                return True
            try:
                bp = int(float(row[bp_index]))
            except (TypeError, ValueError):
                return True
            if bp <= 0:
                return True
    return False


def delete_prefix(prefix: Path) -> None:
    for suffix in [".besd", ".epi", ".esi", ".bak.epi", ".bak.esi"]:
        path = Path(str(prefix) + suffix)
        if path.exists():
            path.unlink()


def main() -> None:
    args = parse_args()
    out_dir = BRYOIS_ROOT / "prepared" / args.celltype / f"chr{args.chrom}"
    out_dir.mkdir(parents=True, exist_ok=True)

    filtered_fastqtl = out_dir / f"{args.celltype}.{args.chrom}.filtered.fastqtl.txt"
    besd_prefix = out_dir / f"{args.celltype}_chr{args.chrom}_filtered"
    update_epi = out_dir / f"{args.celltype}_chr{args.chrom}_filtered_update.epi"
    update_esi = out_dir / f"{args.celltype}_chr{args.chrom}_filtered_update.esi"

    run(
        [
            PYTHON,
            str(PREP_FASTQTL),
            "--celltype",
            args.celltype,
            "--chrom",
            str(args.chrom),
        ]
    )

    if has_invalid_bp(update_epi, "epi") or has_invalid_bp(update_esi, "esi"):
        raise RuntimeError(f"Filtered update files still contain NA coordinates: {update_epi} / {update_esi}")

    if (
        not Path(str(besd_prefix) + ".besd").exists()
        or has_invalid_bp(Path(str(besd_prefix) + ".epi"), "epi")
        or has_invalid_bp(Path(str(besd_prefix) + ".esi"), "esi")
    ):
        delete_prefix(besd_prefix)
        run(
            [
                str(SMR_EXE),
                "--eqtl-summary",
                str(filtered_fastqtl),
                "--fastqtl-nominal-format",
                "--make-besd",
                "--out",
                str(besd_prefix),
            ]
        )

    run(
        [
            str(SMR_EXE),
            "--beqtl-summary",
            str(besd_prefix),
            "--update-epi",
            str(update_epi),
            "--update-esi",
            str(update_esi),
        ]
    )

    print(f"prepared_prefix\t{besd_prefix}")


if __name__ == "__main__":
    main()
