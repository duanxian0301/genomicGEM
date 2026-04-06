from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

import pandas as pd


RUNNER = Path(r"D:\codex\GenomicSEM\scripts\04_postgwas\26_run_smr_single_trait.py")
OUT_ROOT = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\single_trait")
STATUS_DIR = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\batch_status")
STATUS_FILE = STATUS_DIR / "brainmeta_resume_status.tsv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Resume missing BrainMeta SMR jobs with throttled parallelism.")
    parser.add_argument("--traits", nargs="+", default=["AD", "PD", "LBD", "F1", "F2"])
    parser.add_argument("--chromosomes", nargs="+", type=int, default=list(range(1, 23)))
    parser.add_argument("--max-parallel", type=int, default=4)
    parser.add_argument("--poll-seconds", type=int, default=15)
    return parser.parse_args()


def missing_jobs(traits: list[str], chromosomes: list[int]) -> list[dict]:
    jobs = []
    for trait in traits:
        for chrom in chromosomes:
            out_prefix = OUT_ROOT / trait / "brainmeta_v2_cortex" / f"{trait}_BrainMeta_chr{chrom}"
            if out_prefix.with_suffix(".smr").exists():
                continue
            jobs.append({"trait": trait, "chromosome": chrom, "out_prefix": out_prefix})
    return jobs


def launch_job(job: dict) -> subprocess.Popen:
    cmd = [
        "python",
        str(RUNNER),
        "--trait",
        job["trait"],
        "--panel",
        "brainmeta",
        "--chrom",
        str(job["chromosome"]),
    ]
    return subprocess.Popen(cmd)


def main() -> None:
    args = parse_args()
    STATUS_DIR.mkdir(parents=True, exist_ok=True)

    queue = missing_jobs(args.traits, args.chromosomes)
    print(f"missing_jobs\t{len(queue)}")
    if not queue:
        pd.DataFrame().to_csv(STATUS_FILE, sep="\t", index=False)
        print(f"saved_status\t{STATUS_FILE}")
        return

    running: dict[int, tuple[subprocess.Popen, dict]] = {}
    completed: list[dict] = []

    while queue or running:
        while queue and len(running) < args.max_parallel:
            job = queue.pop(0)
            proc = launch_job(job)
            meta = {
                "trait": job["trait"],
                "chromosome": job["chromosome"],
                "pid": proc.pid,
                "launch_time": pd.Timestamp.now().isoformat(timespec="seconds"),
                "out_prefix": str(job["out_prefix"]),
                "status": "running",
            }
            running[proc.pid] = (proc, meta)
            print(f"launched\t{job['trait']}\tchr{job['chromosome']}\tpid={proc.pid}")

        time.sleep(args.poll_seconds)

        for pid in list(running):
            proc, meta = running[pid]
            code = proc.poll()
            if code is None:
                continue
            out_prefix = Path(meta["out_prefix"])
            smr_file = out_prefix.with_suffix(".smr")
            log_file = out_prefix.with_suffix(".run.log")
            meta["status"] = "completed" if smr_file.exists() else "failed_or_incomplete"
            meta["return_code"] = code
            meta["finish_time"] = pd.Timestamp.now().isoformat(timespec="seconds")
            meta["smr_exists"] = smr_file.exists()
            meta["log_exists"] = log_file.exists()
            completed.append(meta)
            print(f"finished\t{meta['trait']}\tchr{meta['chromosome']}\t{meta['status']}\trc={code}")
            del running[pid]

    pd.DataFrame(completed).sort_values(["trait", "chromosome"]).to_csv(STATUS_FILE, sep="\t", index=False)
    print(f"saved_status\t{STATUS_FILE}")


if __name__ == "__main__":
    main()
