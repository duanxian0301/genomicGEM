import csv
import gzip
import math
from pathlib import Path


OUT_DIR = Path(r"D:\codex\GenomicSEM\postgwas_mixer_f1f2_ndd\inputs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# AD sample definition is inferred from the Bellenguez et al. meta-analysis
# underlying GCST90027158: 39,106 clinical cases + 46,828 proxy cases + 401,577 controls.
TRAITS = [
    {
        "trait": "F1",
        "source": Path(r"D:\文章\GS\GWAS\step11_factor_gwas_native_results\merged\standard_txt\ALPS_F1_factorGWAS_native_standard.txt"),
        "n": 30998.0,
    },
    {
        "trait": "F2",
        "source": Path(r"D:\文章\GS\GWAS\step11_factor_gwas_native_results\merged\standard_txt\ALPS_F2_factorGWAS_native_standard.txt"),
        "n": 31629.0,
    },
    {
        "trait": "AD",
        "source": Path(r"D:\文章\4NDD\NDDGWAS\AD.txt"),
        "n": 4.0 / (1.0 / (39106.0 + 46828.0) + 1.0 / 401577.0),
    },
    {
        "trait": "PD",
        "source": Path(r"D:\文章\4NDD\NDDGWAS\PD.txt"),
        "n": 4.0 / (1.0 / (63555.0 + 17700.0) + 1.0 / 1746386.0),
    },
    {
        "trait": "LBD",
        "source": Path(r"D:\文章\4NDD\NDDGWAS\LBD.txt"),
        "n": 4.0 / (1.0 / 2591.0 + 1.0 / 4027.0),
    },
]


def parse_chr(value: str):
    if value is None:
        return None
    text = value.strip().upper().replace("CHR", "")
    if text in {"X", "Y", "MT", "M"}:
        return None
    try:
        chrom = int(float(text))
    except ValueError:
        return None
    if 1 <= chrom <= 22:
        return chrom
    return None


def to_float(value: str):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def transform_trait(trait_cfg):
    source = trait_cfg["source"]
    target = OUT_DIR / f"{trait_cfg['trait']}.sumstats.gz"
    kept = 0
    skipped = 0

    with source.open("r", encoding="utf-8", errors="ignore", newline="") as fin, gzip.open(target, "wt", encoding="utf-8", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t", lineterminator="\n")
        writer.writerow(["SNP", "CHR", "BP", "A1", "A2", "N", "Z"])

        for row in reader:
            snp = (row.get("SNP") or "").strip()
            a1 = (row.get("A1") or "").strip().upper()
            a2 = (row.get("A2") or "").strip().upper()
            chrom = parse_chr(row.get("CHR"))
            bp = to_float(row.get("BP"))
            beta = to_float(row.get("BETA"))
            se = to_float(row.get("SE"))

            if not snp or not a1 or not a2 or chrom is None or bp is None or beta is None or se is None:
                skipped += 1
                continue
            if se <= 0 or not math.isfinite(beta) or not math.isfinite(se):
                skipped += 1
                continue

            z = beta / se
            if not math.isfinite(z):
                skipped += 1
                continue

            writer.writerow([
                snp,
                chrom,
                int(bp),
                a1,
                a2,
                f"{trait_cfg['n']:.6f}",
                f"{z:.10f}",
            ])
            kept += 1

    return {
        "trait": trait_cfg["trait"],
        "source": str(source),
        "output": str(target),
        "N_used": trait_cfg["n"],
        "rows_kept": kept,
        "rows_skipped": skipped,
    }


def main():
    manifest_path = OUT_DIR.parent / "input_manifest.tsv"
    rows = [transform_trait(cfg) for cfg in TRAITS]
    with manifest_path.open("w", encoding="utf-8", newline="") as fout:
        writer = csv.DictWriter(
            fout,
            fieldnames=["trait", "source", "output", "N_used", "rows_kept", "rows_skipped"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote manifest: {manifest_path}")
    for row in rows:
        print(row)


if __name__ == "__main__":
    main()
