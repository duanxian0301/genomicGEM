from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

import psutil


SCRIPT = Path(r"D:\codex\GenomicSEM\scripts\04_postgwas\26_run_smr_single_trait.py")
OUT_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\single_trait")

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

TRAITS = ["AD", "PD", "LBD", "F1", "F2"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Launch and maintain a modest-concurrency GTEx SMR queue.")
    parser.add_argument("--max-total-smr", type=int, default=3, help="Maximum simultaneous smr.exe processes system-wide.")
    parser.add_argument("--poll-seconds", type=int, default=30, help="Polling interval while waiting for slots.")
    parser.add_argument("--max-launch", type=int, default=6, help="Maximum new tasks to launch during this run.")
    parser.add_argument(
        "--traits",
        nargs="*",
        default=TRAITS,
        choices=TRAITS,
        help="Subset of traits to queue.",
    )
    parser.add_argument(
        "--tissues",
        nargs="*",
        default=GTEX_BRAIN_TISSUES,
        choices=GTEX_BRAIN_TISSUES,
        help="Subset of GTEx tissues to queue.",
    )
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


def output_exists(trait: str, tissue: str) -> bool:
    out = OUT_ROOT / trait / "gtex_v8_brain" / f"{trait}_{tissue}.smr"
    return out.exists() and out.stat().st_size > 0


def build_queue(traits: list[str], tissues: list[str]) -> list[tuple[str, str]]:
    queue: list[tuple[str, str]] = []
    for trait in traits:
        for tissue in tissues:
            if not output_exists(trait, tissue):
                queue.append((trait, tissue))
    return queue


def launch_one(trait: str, tissue: str) -> None:
    cmd = [
        sys.executable,
        str(SCRIPT),
        "--trait",
        trait,
        "--panel",
        "gtex",
        "--tissue",
        tissue,
    ]
    subprocess.Popen(cmd, creationflags=subprocess.CREATE_NEW_PROCESS_GROUP)


def main() -> None:
    args = parse_args()
    queue = build_queue(list(args.traits), list(args.tissues))
    print(f"pending_tasks\t{len(queue)}")
    launched = 0

    while queue and launched < args.max_launch:
        active = count_active_smr()
        if active >= args.max_total_smr:
            time.sleep(args.poll_seconds)
            continue
        trait, tissue = queue.pop(0)
        launch_one(trait, tissue)
        launched += 1
        print(f"launched\t{trait}\t{tissue}")
        time.sleep(3)

    print(f"launched_total\t{launched}")
    print(f"remaining_pending\t{len(queue)}")


if __name__ == "__main__":
    main()
