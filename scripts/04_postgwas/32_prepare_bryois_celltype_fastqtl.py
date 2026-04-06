from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import pandas as pd


DEFAULT_BRYOIS_DIR = Path(r"D:\文章\GS\postgwas\06_smr_f1f2_ndd\refs\bryois2022_celltype_eqtl")
DEFAULT_VALID_GENES = Path(r"D:\sLDSC\data\valid_genes.tsv.gz")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare Bryois 2022 cell-type fastQTL summary stats for SMR BESD conversion."
    )
    parser.add_argument("--celltype", required=True, help="Cell type label, e.g. Microglia")
    parser.add_argument("--chrom", required=True, type=int, help="Chromosome number, e.g. 6")
    parser.add_argument(
        "--bryois-dir",
        default=str(DEFAULT_BRYOIS_DIR),
        help="Root directory for Bryois downloads and prepared outputs.",
    )
    parser.add_argument(
        "--valid-genes",
        default=str(DEFAULT_VALID_GENES),
        help="Gene-position reference used to fill probe coordinates.",
    )
    parser.add_argument(
        "--snp-pos",
        default="",
        help="Optional override for snp_pos.txt(.gz). Defaults to file inside Bryois root.",
    )
    parser.add_argument(
        "--download-dir",
        default="downloads",
        help="Subdirectory under Bryois root for downloaded fastQTL files.",
    )
    parser.add_argument(
        "--out-dir",
        default="prepared",
        help="Subdirectory under Bryois root for prepared outputs.",
    )
    return parser.parse_args()


def read_pairs(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes: set[str] = set()
    snps: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split()
            if len(parts) < 5:
                continue
            genes.add(parts[0])
            snps.add(parts[1])
    return pd.DataFrame({"probe": sorted(genes)}), pd.DataFrame({"snp": sorted(snps)})


def write_filtered_fastqtl(
    in_path: Path, out_path: Path, valid: pd.DataFrame, snp_pos: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    valid = valid.copy()
    valid["start"] = pd.to_numeric(valid["start"], errors="coerce")
    valid_set = set(valid.loc[valid["start"].notna() & (valid["start"] > 0), "ensgid"])
    snp_pos = snp_pos.copy()
    if "pos_hg19" in snp_pos.columns:
        snp_pos["bp_candidate"] = pd.to_numeric(snp_pos["pos_hg19"], errors="coerce")
    elif "SNP_id_hg19" in snp_pos.columns:
        snp_pos["bp_candidate"] = pd.to_numeric(
            snp_pos["SNP_id_hg19"].astype(str).str.extract(r"^chr[^:]+:(?P<bp>\d+)$")["bp"],
            errors="coerce",
        )
    else:
        snp_pos["bp_candidate"] = pd.NA
    snp_set = set(snp_pos.loc[snp_pos["bp_candidate"].notna() & (snp_pos["bp_candidate"] > 0), "SNP"])
    probes: set[str] = set()
    snps: set[str] = set()
    kept = 0
    with gzip.open(in_path, "rt", encoding="utf-8") as fin, out_path.open("w", encoding="utf-8") as fout:
        for line in fin:
            parts = line.rstrip("\n").split()
            if len(parts) < 5:
                continue
            probe, snp = parts[0], parts[1]
            ensg = probe.rsplit("_", 1)[-1]
            if ensg not in valid_set or snp not in snp_set:
                continue
            fout.write(line)
            kept += 1
            probes.add(probe)
            snps.add(snp)
    return pd.DataFrame({"probe": sorted(probes)}), pd.DataFrame({"snp": sorted(snps)}), kept


def make_update_epi(probe_df: pd.DataFrame, valid: pd.DataFrame) -> pd.DataFrame:
    valid = valid.copy()
    valid["start"] = pd.to_numeric(valid["start"], errors="coerce")
    valid["chr_num"] = valid["chr"].astype(str).str.replace("^chr", "", regex=True)
    probe_df = probe_df.copy()
    probe_df["ensgid"] = probe_df["probe"].str.extract(r"(ENSG\d+)$")
    probe_df["gene_symbol"] = probe_df["probe"].str.replace(r"_ENSG\d+$", "", regex=True)
    merged = probe_df.merge(valid, on="ensgid", how="left")
    out = pd.DataFrame(
        {
            "chr": merged["chr_num"],
            "probe": merged["probe"],
            "genetic_distance": 0,
            "bp": merged["start"],
            "gene": merged["gene_symbol"],
            "strand": "NA",
        }
    )
    return out.loc[out["bp"].notna() & (pd.to_numeric(out["bp"], errors="coerce") > 0)].copy()


def make_update_esi(snp_df: pd.DataFrame, snp_pos: pd.DataFrame) -> pd.DataFrame:
    snp_pos = snp_pos.copy()
    if "chr" in snp_pos.columns:
        snp_pos["chr_num"] = snp_pos["chr"].astype(str).str.replace("^chr", "", regex=True)
    elif "SNP_id_hg19" in snp_pos.columns:
        chr_bp = snp_pos["SNP_id_hg19"].astype(str).str.extract(r"^chr(?P<chr>[^:]+):(?P<bp>\d+)$")
        snp_pos["chr_num"] = chr_bp["chr"]
        snp_pos["pos_hg19"] = pd.to_numeric(chr_bp["bp"], errors="coerce")
    else:
        raise KeyError("snp_pos must contain either chr/pos_hg19 or SNP_id_hg19 columns.")

    if "pos_hg19" not in snp_pos.columns:
        raise KeyError("snp_pos must contain hg19 position information.")

    merged = snp_df.merge(snp_pos, left_on="snp", right_on="SNP", how="left")
    out = pd.DataFrame(
        {
            "chr": merged["chr_num"],
            "snp": merged["snp"],
            "genetic_distance": 0,
            "bp": merged["pos_hg19"],
            "a1": merged["effect_allele"],
            "a2": merged["other_allele"],
            "freq": "NA",
        }
    )
    return out.loc[out["bp"].notna() & (pd.to_numeric(out["bp"], errors="coerce") > 0)].copy()


def resolve_snp_pos(root: Path, override: str) -> Path:
    if override:
        return Path(override)
    txt = root / "downloads" / "snp_pos.txt"
    gz = root / "downloads" / "snp_pos.txt.gz"
    legacy_txt = root / "test_download" / "snp_pos.txt"
    legacy_gz = root / "test_download" / "snp_pos.txt.gz"
    for path in (txt, gz, legacy_txt, legacy_gz):
        if path.exists():
            return path
    raise FileNotFoundError("Could not locate snp_pos.txt(.gz) under Bryois root.")


def main() -> None:
    args = parse_args()
    bryois_root = Path(args.bryois_dir)
    valid_genes = pd.read_csv(args.valid_genes, sep="\t", compression="gzip")
    snp_pos_path = resolve_snp_pos(bryois_root, args.snp_pos)
    snp_pos = pd.read_csv(snp_pos_path, sep="\t", compression="gzip" if snp_pos_path.suffix == ".gz" else None)

    download_dir = bryois_root / args.download_dir / args.celltype
    in_path = download_dir / f"{args.celltype}.{args.chrom}.gz"
    if not in_path.exists():
        raise FileNotFoundError(f"Input fastQTL file not found: {in_path}")

    out_dir = bryois_root / args.out_dir / args.celltype / f"chr{args.chrom}"
    out_dir.mkdir(parents=True, exist_ok=True)

    raw_probe_df, raw_snp_df = read_pairs(in_path)
    raw_epi = make_update_epi(raw_probe_df, valid_genes)
    raw_esi = make_update_esi(raw_snp_df, snp_pos)
    raw_epi_path = out_dir / f"{args.celltype}_chr{args.chrom}_update.epi"
    raw_esi_path = out_dir / f"{args.celltype}_chr{args.chrom}_update.esi"
    raw_epi.to_csv(raw_epi_path, sep="\t", index=False, header=False, na_rep="NA")
    raw_esi.to_csv(raw_esi_path, sep="\t", index=False, header=False, na_rep="NA")

    filtered_fastqtl_path = out_dir / f"{args.celltype}.{args.chrom}.filtered.fastqtl.txt"
    filtered_probe_df, filtered_snp_df, kept = write_filtered_fastqtl(
        in_path, filtered_fastqtl_path, valid_genes, snp_pos
    )
    filtered_epi = make_update_epi(filtered_probe_df, valid_genes)
    filtered_esi = make_update_esi(filtered_snp_df, snp_pos)
    filtered_epi_path = out_dir / f"{args.celltype}_chr{args.chrom}_filtered_update.epi"
    filtered_esi_path = out_dir / f"{args.celltype}_chr{args.chrom}_filtered_update.esi"
    filtered_epi.to_csv(filtered_epi_path, sep="\t", index=False, header=False, na_rep="NA")
    filtered_esi.to_csv(filtered_esi_path, sep="\t", index=False, header=False, na_rep="NA")

    print(f"input_fastqtl\t{in_path}")
    print(f"snp_pos\t{snp_pos_path}")
    print(f"out_dir\t{out_dir}")
    print(f"raw_probes\t{len(raw_probe_df)}")
    print(f"raw_snps\t{len(raw_snp_df)}")
    print(f"raw_missing_probe_bp\t{int(raw_epi['bp'].isna().sum())}")
    print(f"raw_missing_snp_bp\t{int(raw_esi['bp'].isna().sum())}")
    print(f"filtered_fastqtl\t{filtered_fastqtl_path}")
    print(f"filtered_rows\t{kept}")
    print(f"filtered_probes\t{len(filtered_probe_df)}")
    print(f"filtered_snps\t{len(filtered_snp_df)}")
    print(f"filtered_missing_probe_bp\t{int(filtered_epi['bp'].isna().sum())}")
    print(f"filtered_missing_snp_bp\t{int(filtered_esi['bp'].isna().sum())}")
    print(f"filtered_epi\t{filtered_epi_path}")
    print(f"filtered_esi\t{filtered_esi_path}")


if __name__ == "__main__":
    main()
