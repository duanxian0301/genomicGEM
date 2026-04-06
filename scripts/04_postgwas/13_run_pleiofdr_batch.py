from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


PFDIR = Path(r"D:\pleioFDR\pleiofdr-master")


def pair_done(cfg_path: Path) -> bool:
    pair = cfg_path.name.replace(".config.txt", "")
    result_dir = cfg_path.parent.parent / "results" / pair
    out_csv = result_dir / f"{pair}_conjfdr_0.05_all.csv"
    loci_csv = result_dir / f"{pair}_conjfdr_0.05_loci.csv"
    manhattan = result_dir / f"{pair}_conjfdr_0.05_manhattan.png"
    return out_csv.exists() and loci_csv.exists() and manhattan.exists()


def run_config(cfg_path: Path) -> None:
    cfg_forward = cfg_path.as_posix()
    command = [
        "matlab",
        "-batch",
        f"cd('{PFDIR.as_posix()}'); clear config; config='{cfg_forward}'; run('runme.m');",
    ]
    subprocess.run(command, cwd=PFDIR, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("configs", nargs="+")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    for raw_cfg in args.configs:
        cfg_path = Path(raw_cfg)
        if not cfg_path.exists():
            raise FileNotFoundError(cfg_path)
        if not args.force and pair_done(cfg_path):
            print(f"Skipping completed pair: {cfg_path.name}")
            continue
        print(f"Running pair: {cfg_path.name}")
        run_config(cfg_path)
        print(f"Completed pair: {cfg_path.name}")


if __name__ == "__main__":
    main()
