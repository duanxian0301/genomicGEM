from __future__ import annotations

from pathlib import Path

import pandas as pd
from openpyxl import load_workbook


BASE = Path(r"D:\文章\GS\postgwas")
SMR_BASE = BASE / "06_smr_f1f2_ndd"
SMR_SINGLE = SMR_BASE / "single_trait"
SMR_SUMMARY = SMR_BASE / "summary"
CTWAS_SUMMARY = BASE / "05_mctwas_f1f2_ndd" / "summary"
COLOC_SUMMARY = BASE / "07_coloc_factor_ndd" / "coloc_factor_ndd_summary.tsv"
PWCOCO_BEST = BASE / "08_pwcoco" / "pwcoco_region_best_summary.tsv"
SUPP_XLSX = BASE / "Supplementary_LDSC_MiXeR_pleioFDR_tables.xlsx"


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


def collect_smr_results() -> pd.DataFrame:
    frames = []
    for path in sorted(SMR_SINGLE.rglob("*.smr")):
        if path.stat().st_size == 0:
            continue
        trait = path.parent.parent.name
        reference = path.parent.name
        stem = path.stem
        chrom = int(stem.rsplit("chr", 1)[1]) if "chr" in stem else pd.NA
        try:
            df = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        df["trait"] = trait
        df["reference"] = reference
        df["chrom"] = chrom
        df["source_file"] = str(path)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = pd.concat(frames, ignore_index=True)
    for col in ["p_SMR", "p_HEIDI", "b_SMR", "se_SMR", "p_GWAS", "p_eQTL"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out["heidi_pass"] = out["p_HEIDI"].isna() | (out["p_HEIDI"] >= 0.01)
    out["smr_sig_5e8"] = out["p_SMR"] < 5e-8
    out["smr_sig_1e6"] = out["p_SMR"] < 1e-6
    out["smr_sig_1e4"] = out["p_SMR"] < 1e-4
    return out


def summarize_trait_level(smr: pd.DataFrame) -> pd.DataFrame:
    if smr.empty:
        return pd.DataFrame()
    rows = []
    for trait, sub in smr.groupby("trait"):
        best = sub.sort_values("p_SMR").iloc[0]
        rows.append(
            {
                "trait": trait,
                "n_results": len(sub),
                "n_chr_completed": sub["chrom"].nunique(),
                "n_sig_5e8": int(sub["smr_sig_5e8"].sum()),
                "n_sig_5e8_heidi": int((sub["smr_sig_5e8"] & sub["heidi_pass"]).sum()),
                "n_sig_1e6": int(sub["smr_sig_1e6"].sum()),
                "n_sig_1e6_heidi": int((sub["smr_sig_1e6"] & sub["heidi_pass"]).sum()),
                "n_sig_1e4": int(sub["smr_sig_1e4"].sum()),
                "n_sig_1e4_heidi": int((sub["smr_sig_1e4"] & sub["heidi_pass"]).sum()),
                "best_gene": best["Gene"],
                "best_probe": best["probeID"],
                "best_chr": best["chrom"],
                "best_topSNP": best["topSNP"],
                "best_p_SMR": best["p_SMR"],
                "best_p_HEIDI": best["p_HEIDI"],
                "best_source_file": best["source_file"],
            }
        )
    return pd.DataFrame(rows).sort_values("trait")


def load_ctwas_genes() -> pd.DataFrame:
    path = CTWAS_SUMMARY / "ctwas_trait_level_gene_results.tsv"
    df = pd.read_csv(path, sep="\t")
    return df


def load_convergence() -> tuple[pd.DataFrame, pd.DataFrame]:
    coloc = pd.read_csv(COLOC_SUMMARY, sep="\t")
    pw = pd.read_csv(PWCOCO_BEST, sep="\t")
    return coloc, pw


def load_conjfdr_loci() -> pd.DataFrame:
    frames = []
    base = BASE / "03_pleiofdr_all_pairs" / "results"
    for pair, (trait, disease) in PAIR_TO_TRAITS.items():
        loci = base / pair / f"{pair}_conjfdr_0.05_loci.csv"
        if loci.exists():
            df = pd.read_csv(loci)
            df["pair"] = pair
            df["trait"] = trait
            df["disease"] = disease
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def build_candidate_table(
    smr: pd.DataFrame, ctwas: pd.DataFrame, coloc: pd.DataFrame, pw: pd.DataFrame
) -> pd.DataFrame:
    if smr.empty:
        return pd.DataFrame()
    smr_best = (
        smr.sort_values(["trait", "p_SMR"])
        .groupby(["trait", "Gene"], as_index=False)
        .first()[["trait", "Gene", "chrom", "topSNP", "p_SMR", "p_HEIDI", "heidi_pass", "source_file"]]
    )
    ct = ctwas.rename(columns={"gene_symbol": "Gene"})[
        ["trait", "Gene", "priority_label", "best_tissue", "best_p", "max_pip"]
    ].copy()
    ct["ctwas_supported"] = True

    coloc_pairs = coloc.loc[coloc["status"] == "ok", ["pair", "coloc_interpretation"]].copy()
    coloc_pairs["coloc_supported"] = coloc_pairs["coloc_interpretation"].isin(["Strong_H4", "Moderate_H4"])
    pw_pairs = pw[["pair", "H4_uncond", "H4_best_cond", "result_type_best_cond"]].copy()
    pw_pairs["pwcoco_supported"] = (pw_pairs["H4_uncond"] >= 0.5) | (pw_pairs["H4_best_cond"] >= 0.5)

    pair_support_rows = []
    for pair, (trait, disease) in PAIR_TO_TRAITS.items():
        for member in (trait, disease):
            pair_support_rows.append({"pair": pair, "trait": member})
    pair_support = pd.DataFrame(pair_support_rows)
    pair_support = pair_support.merge(coloc_pairs, on="pair", how="left").merge(pw_pairs, on="pair", how="left")

    merged = smr_best.merge(ct, on=["trait", "Gene"], how="left")
    merged = merged.merge(pair_support.groupby("trait", as_index=False).agg(
        any_coloc_supported=("coloc_supported", "max"),
        any_pwcoco_supported=("pwcoco_supported", "max"),
        coloc_pairs=("pair", lambda s: ";".join(sorted(set(x for x in s.dropna())))),
    ), on="trait", how="left")

    def assign_tier(row: pd.Series) -> str:
        strong_smr = bool(row["p_SMR"] < 1e-4 and (pd.isna(row["p_HEIDI"]) or row["p_HEIDI"] >= 0.01))
        ctwas_ok = bool(row.get("ctwas_supported", False))
        coloc_ok = bool(row.get("any_coloc_supported", False) or row.get("any_pwcoco_supported", False))
        if strong_smr and ctwas_ok and coloc_ok:
            return "A_high_convergent"
        if strong_smr or ctwas_ok or coloc_ok:
            return "B_moderate_support"
        return "C_exploratory"

    merged["candidate_tier"] = merged.apply(assign_tier, axis=1)
    return merged.sort_values(["candidate_tier", "trait", "p_SMR"])


def replace_or_add_sheet(xlsx: Path, sheet_name: str, df: pd.DataFrame) -> None:
    with pd.ExcelWriter(xlsx, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

    wb = load_workbook(xlsx)
    ws = wb[sheet_name[:31]]
    autosize_worksheet(ws)
    ws.freeze_panes = "A2"
    wb.save(xlsx)


def main() -> None:
    SMR_SUMMARY.mkdir(parents=True, exist_ok=True)
    smr = collect_smr_results()
    trait_summary = summarize_trait_level(smr)
    ctwas = load_ctwas_genes()
    coloc, pw = load_convergence()
    conjfdr = load_conjfdr_loci()
    candidates = build_candidate_table(smr, ctwas, coloc, pw)

    smr.to_csv(SMR_SUMMARY / "smr_all_results.tsv", sep="\t", index=False)
    trait_summary.to_csv(SMR_SUMMARY / "smr_trait_summary.tsv", sep="\t", index=False)
    candidates.to_csv(SMR_SUMMARY / "smr_candidate_convergence.tsv", sep="\t", index=False)
    conjfdr.to_csv(SMR_SUMMARY / "conjfdr_pair_loci_compiled.tsv", sep="\t", index=False)

    replace_or_add_sheet(SUPP_XLSX, "SMR_trait_summary", trait_summary)
    replace_or_add_sheet(SUPP_XLSX, "SMR_all_results", smr)
    replace_or_add_sheet(SUPP_XLSX, "SMR_candidates", candidates)
    replace_or_add_sheet(SUPP_XLSX, "SMR_conjFDR_loci", conjfdr)

    print(f"Wrote SMR summaries to: {SMR_SUMMARY}")
    print(f"Updated workbook: {SUPP_XLSX}")


if __name__ == "__main__":
    main()
