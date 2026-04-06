from __future__ import annotations

from pathlib import Path

import pandas as pd
from openpyxl import load_workbook


POSTGWAS = Path("D:/文章/GS/postgwas")
SMR_CURATED = POSTGWAS / "06_smr_f1f2_ndd" / "summary_curated"
SMR_SUMMARY = POSTGWAS / "06_smr_f1f2_ndd" / "summary"
CTWAS_SUMMARY = POSTGWAS / "05_mctwas_f1f2_ndd" / "summary" / "ctwas_trait_level_gene_results.tsv"
COLOC_SUMMARY = POSTGWAS / "07_coloc_factor_ndd" / "coloc_factor_ndd_summary.tsv"
PWCOCO_SUMMARY = POSTGWAS / "08_pwcoco" / "pwcoco_region_best_summary.tsv"
WORKBOOK = POSTGWAS / "Supplementary_LDSC_MiXeR_pleioFDR_tables.xlsx"
OUT_DIR = POSTGWAS / "09_candidate_gene_integration"

PAIR_TO_TRAITS = {
    "F1_AD": ("F1", "AD"),
    "F2_AD": ("F2", "AD"),
    "F1_PD": ("F1", "PD"),
    "F2_PD": ("F2", "PD"),
    "F1_LBD": ("F1", "LBD"),
    "F2_LBD": ("F2", "LBD"),
}


def autosize_worksheet(ws) -> None:
    for column_cells in ws.columns:
        values = [str(cell.value) if cell.value is not None else "" for cell in column_cells]
        max_len = max((len(v) for v in values), default=0)
        ws.column_dimensions[column_cells[0].column_letter].width = min(max(max_len + 2, 10), 40)


def replace_or_add_sheet(xlsx: Path, sheet_name: str, df: pd.DataFrame) -> None:
    with pd.ExcelWriter(xlsx, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        safe = df.copy()
        for col in safe.columns:
            if pd.api.types.is_bool_dtype(safe[col]):
                safe[col] = safe[col].astype(str)
        safe.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    wb = load_workbook(xlsx)
    ws = wb[sheet_name[:31]]
    ws.freeze_panes = "A2"
    autosize_worksheet(ws)
    wb.save(xlsx)


def load_smr_support(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["trait", "Gene"])
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["trait", "Gene"])
    agg = (
        df.sort_values(["trait", "Gene", "p_SMR"], kind="stable")
        .groupby(["trait", "Gene"], as_index=False)
        .first()[["trait", "Gene", "context", "topSNP", "p_SMR", "p_HEIDI"]]
        .rename(
            columns={
                "context": f"{label}_best_context",
                "topSNP": f"{label}_topSNP",
                "p_SMR": f"{label}_best_p_SMR",
                "p_HEIDI": f"{label}_best_p_HEIDI",
            }
        )
    )
    agg[f"{label}_supported"] = True
    return agg


def load_ctwas() -> pd.DataFrame:
    if not CTWAS_SUMMARY.exists():
        return pd.DataFrame(columns=["trait", "Gene"])
    df = pd.read_csv(CTWAS_SUMMARY, sep="\t")
    if df.empty:
        return pd.DataFrame(columns=["trait", "Gene"])
    df["Gene"] = df["gene_symbol"]
    df["ctwas_primary_supported"] = pd.to_numeric(df["max_pip"], errors="coerce") >= 0.5
    df["ctwas_secondary_supported"] = (
        pd.to_numeric(df["min_FDR"], errors="coerce") < 0.05
    ) & (~df["ctwas_primary_supported"])
    return df[
        [
            "trait",
            "Gene",
            "best_tissue",
            "best_feature",
            "best_p",
            "min_FDR",
            "max_pip",
            "priority_label",
            "ctwas_primary_supported",
            "ctwas_secondary_supported",
        ]
    ]


def load_pair_support() -> pd.DataFrame:
    coloc = pd.read_csv(COLOC_SUMMARY, sep="\t")
    pw = pd.read_csv(PWCOCO_SUMMARY, sep="\t")

    coloc_ok = coloc.loc[coloc["status"].eq("ok"), ["pair", "coloc_interpretation"]].copy()
    coloc_ok["coloc_supported"] = coloc_ok["coloc_interpretation"].isin(["Strong_H4", "Moderate_H4"])

    pw = pw[["pair", "H4_uncond", "H4_best_cond"]].copy()
    pw["pwcoco_supported"] = (
        pd.to_numeric(pw["H4_uncond"], errors="coerce").fillna(0) >= 0.5
    ) | (
        pd.to_numeric(pw["H4_best_cond"], errors="coerce").fillna(0) >= 0.5
    )

    rows = []
    for pair, (trait1, trait2) in PAIR_TO_TRAITS.items():
        rows.append({"pair": pair, "trait": trait1})
        rows.append({"pair": pair, "trait": trait2})
    support = pd.DataFrame(rows)
    support = support.merge(coloc_ok, on="pair", how="left").merge(pw, on="pair", how="left")
    out = support.groupby("trait", as_index=False).agg(
        coloc_pairs=("pair", lambda s: ";".join(sorted(set(x for x in s.dropna())))),
        any_coloc_supported=("coloc_supported", "max"),
        any_pwcoco_supported=("pwcoco_supported", "max"),
    )
    return out


def assign_tier(row: pd.Series) -> str:
    has_ctwas = bool(row.get("ctwas_primary_supported", False) or row.get("ctwas_secondary_supported", False))
    has_smr_bulk = bool(row.get("smr_bulk_bonf_supported", False) or row.get("smr_bulk_fdr_supported", False))
    has_smr_cell = bool(row.get("smr_cell_bonf_supported", False) or row.get("smr_cell_fdr_supported", False))
    has_snp = bool(row.get("any_coloc_supported", False) or row.get("any_pwcoco_supported", False))
    evidence_count = sum([has_ctwas, has_smr_bulk or has_smr_cell, has_snp])

    if bool(row.get("ctwas_primary_supported", False)) and (has_smr_bulk or has_smr_cell) and has_snp:
        return "A_high_convergent"
    if evidence_count >= 2:
        return "B_multi_source"
    if evidence_count == 1:
        return "C_single_source"
    return "D_other"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    bulk_bonf = load_smr_support(SMR_CURATED / "SMR_bulk_bonferroni.tsv", "smr_bulk_bonf")
    bulk_fdr = load_smr_support(SMR_CURATED / "SMR_bulk_fdr.tsv", "smr_bulk_fdr")
    cell_bonf = load_smr_support(SMR_CURATED / "SMR_celltype_bonferroni.tsv", "smr_cell_bonf")
    cell_fdr = load_smr_support(SMR_CURATED / "SMR_celltype_fdr.tsv", "smr_cell_fdr")
    ctwas = load_ctwas()
    pair_support = load_pair_support()

    frames = [df for df in [bulk_bonf, bulk_fdr, cell_bonf, cell_fdr, ctwas] if not df.empty]
    master = pd.concat(frames, ignore_index=True, sort=False) if frames else pd.DataFrame(columns=["trait", "Gene"])
    master = master.groupby(["trait", "Gene"], as_index=False).first()
    master = master.merge(pair_support, on="trait", how="left")

    bool_cols = [c for c in master.columns if c.endswith("_supported") or c.startswith("any_")]
    for col in bool_cols:
        master[col] = master[col].fillna(False).astype(bool)

    master["candidate_tier"] = master.apply(assign_tier, axis=1)
    master["evidence_ctwas"] = master["ctwas_primary_supported"] | master["ctwas_secondary_supported"]
    master["evidence_smr_bulk"] = master["smr_bulk_bonf_supported"] | master["smr_bulk_fdr_supported"]
    master["evidence_smr_celltype"] = master["smr_cell_bonf_supported"] | master["smr_cell_fdr_supported"]
    master["evidence_snp_locus"] = master["any_coloc_supported"] | master["any_pwcoco_supported"]

    master = master.sort_values(
        ["candidate_tier", "trait", "max_pip", "smr_bulk_bonf_best_p_SMR", "smr_cell_bonf_best_p_SMR"],
        ascending=[True, True, False, True, True],
        kind="stable",
    )

    shortlist = master.loc[master["candidate_tier"].isin(["A_high_convergent", "B_multi_source"])].copy()

    tier_summary = (
        master.groupby(["trait", "candidate_tier"], as_index=False)
        .size()
        .rename(columns={"size": "n_genes"})
        .sort_values(["trait", "candidate_tier"], kind="stable")
    )

    master_path = OUT_DIR / "candidate_gene_master.tsv"
    shortlist_path = OUT_DIR / "candidate_gene_shortlist.tsv"
    tier_path = OUT_DIR / "candidate_gene_tier_summary.tsv"
    master.to_csv(master_path, sep="\t", index=False)
    shortlist.to_csv(shortlist_path, sep="\t", index=False)
    tier_summary.to_csv(tier_path, sep="\t", index=False)

    replace_or_add_sheet(WORKBOOK, "Candidate_gene_master", master)
    replace_or_add_sheet(WORKBOOK, "Candidate_gene_shortlist", shortlist)
    replace_or_add_sheet(WORKBOOK, "Candidate_gene_tiers", tier_summary)

    print(f"master_rows\t{len(master)}")
    print(f"shortlist_rows\t{len(shortlist)}")
    print(f"tier_rows\t{len(tier_summary)}")
    print(f"master_path\t{master_path}")


if __name__ == "__main__":
    main()
