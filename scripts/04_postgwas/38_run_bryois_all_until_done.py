from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

import psutil


PYTHON = sys.executable
PREP_REF = Path(r"D:\codex\GenomicSEM\scripts\04_postgwas\34_prepare_bryois_reference.py")
RUN_SMR = Path(r"D:\codex\GenomicSEM\scripts\04_postgwas\35_run_smr_bryois_single.py")
BRYOIS_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\refs\bryois2022_celltype_eqtl")
OUT_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\single_trait")

TRAITS = ["AD", "PD", "LBD", "F1", "F2"]
CELLTYPES = [
    "Astrocytes",
    "Microglia",
    "Endothelial.cells",
    "Excitatory.neurons",
    "Inhibitory.neurons",
    "OPCs...COPs",
    "Oligodendrocytes",
    "Pericytes",
]
CHROMS = list(range(1, 23))
MAX_TOTAL_SMR = 8
POLL_SECONDS = 20


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


def active_task_keys() -> set[tuple[str, str, int]]:
    keys: set[tuple[str, str, int]] = set()
    for proc in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            name = (proc.info.get("name") or "").lower()
            cmd = proc.info.get("cmdline") or []
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
        if not cmd:
            continue
        joined = " ".join(cmd)
        if "35_run_smr_bryois_single.py" not in joined:
            continue
        trait = None
        cell = None
        chrom = None
        for i, token in enumerate(cmd):
            if token == "--trait" and i + 1 < len(cmd):
                trait = cmd[i + 1]
            elif token == "--celltype" and i + 1 < len(cmd):
                cell = cmd[i + 1]
            elif token == "--chrom" and i + 1 < len(cmd):
                try:
                    chrom = int(cmd[i + 1])
                except ValueError:
                    chrom = None
        if trait and cell and chrom:
            keys.add((trait, cell, chrom))
    return keys


def output_exists(trait: str, cell: str, chrom: int) -> bool:
    out = OUT_ROOT / trait / "bryois2022_celltype" / cell / f"{trait}_{cell}_chr{chrom}.smr"
    return out.exists() and out.stat().st_size > 0


def ensure_reference(cell: str, chrom: int) -> None:
    subprocess.run([PYTHON, str(PREP_REF), "--celltype", cell, "--chrom", str(chrom)], check=True)


def launch_task(trait: str, cell: str, chrom: int) -> None:
    subprocess.Popen(
        [PYTHON, str(RUN_SMR), "--trait", trait, "--celltype", cell, "--chrom", str(chrom)],
        creationflags=subprocess.CREATE_NEW_PROCESS_GROUP,
    )


def missing_tasks() -> list[tuple[str, str, int]]:
    todo: list[tuple[str, str, int]] = []
    for cell in CELLTYPES:
        for chrom in CHROMS:
            for trait in TRAITS:
                if not output_exists(trait, cell, chrom):
                    todo.append((trait, cell, chrom))
    return todo


def main() -> None:
    rounds = 0
    while True:
        rounds += 1
        todo = missing_tasks()
        active = active_task_keys()
        ensured_refs: set[tuple[str, int]] = set()
        if not todo and count_active_smr() == 0:
            print("status\tcomplete")
            print(f"rounds\t{rounds}")
            break

        launched = 0
        for trait, cell, chrom in todo:
            if output_exists(trait, cell, chrom):
                continue
            if (trait, cell, chrom) in active:
                continue
            ref_key = (cell, chrom)
            if ref_key not in ensured_refs:
                ensure_reference(cell, chrom)
                ensured_refs.add(ref_key)
            while count_active_smr() >= MAX_TOTAL_SMR:
                time.sleep(POLL_SECONDS)
                active = active_task_keys()
                if (trait, cell, chrom) in active or output_exists(trait, cell, chrom):
                    break
            if output_exists(trait, cell, chrom) or (trait, cell, chrom) in active_task_keys():
                continue
            launch_task(trait, cell, chrom)
            launched += 1
            time.sleep(2)

        print(f"round\t{rounds}\tmissing\t{len(todo)}\tlaunched\t{launched}\tactive_smr\t{count_active_smr()}")
        time.sleep(POLL_SECONDS)


if __name__ == "__main__":
    main()
