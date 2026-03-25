from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(r"D:\文章\GS\GWAS\step11_factor_gwas_native_results\merged")
OUTDIR = ROOT / "q_snp_analysis"
OUTDIR.mkdir(parents=True, exist_ok=True)

F1 = ROOT / "ALPS_F1_factorGWAS_native_merged.tsv.gz"
F2 = ROOT / "ALPS_F2_factorGWAS_native_merged.tsv.gz"


def lead_loci_count(df: pd.DataFrame, p_col: str, window_bp: int = 500_000) -> int:
    sig = df.loc[df[p_col] < 5e-8, ["CHR", "BP", p_col]].copy()
    if sig.empty:
        return 0
    sig = sig.sort_values(["CHR", p_col, "BP"]).reset_index(drop=True)
    lead_n = 0
    for chr_, sub in sig.groupby("CHR", sort=True):
        used = np.zeros(len(sub), dtype=bool)
        sub = sub.reset_index(drop=True)
        for i in range(len(sub)):
            if used[i]:
                continue
            lead_n += 1
            pos = sub.loc[i, "BP"]
            used |= (sub["BP"].between(pos - window_bp, pos + window_bp)).to_numpy()
    return lead_n


def exclude_near_q_hits(factor_df: pd.DataFrame, q_df: pd.DataFrame, factor_p_col: str, q_p_col: str, window_bp: int = 500_000) -> int:
    f = factor_df.loc[factor_df[factor_p_col] < 5e-8, ["SNP", "CHR", "BP", factor_p_col]].copy()
    q = q_df.loc[q_df[q_p_col] < 5e-8, ["CHR", "BP", q_p_col]].copy()
    if f.empty:
        return 0
    if q.empty:
        return len(f)
    keep = []
    for _, row in f.iterrows():
        subq = q.loc[q["CHR"] == row["CHR"]]
        if subq.empty:
            keep.append(True)
            continue
        near = ((subq["BP"] - row["BP"]).abs() <= window_bp).any()
        keep.append(not near)
    return int(np.sum(keep))


def main() -> None:
    f1 = pd.read_csv(F1, sep="\t", compression="gzip")
    f2 = pd.read_csv(F2, sep="\t", compression="gzip")

    # Q_SNP columns are present in both factor outputs, but they are not necessarily identical
    # because the userGWAS run used sub = c("F1~SNP", "F2~SNP"), yielding factor-specific result tables.
    q_cols = ["SNP", "CHR", "BP", "A1", "A2", "MAF", "chisq", "chisq_df", "chisq_pval", "Q_SNP", "Q_SNP_df", "Q_SNP_pval"]
    q1 = f1[q_cols].copy()
    q2 = f2[q_cols].copy()

    q_compare = pd.DataFrame({
        "metric": [
            "same_n_rows",
            "same_snp_order",
            "max_abs_delta_Q_SNP",
            "max_abs_delta_Q_SNP_pval",
        ],
        "value": [
            len(q1) == len(q2),
            bool((q1["SNP"].values == q2["SNP"].values).all()),
            float((q1["Q_SNP"] - q2["Q_SNP"]).abs().max()),
            float((q1["Q_SNP_pval"] - q2["Q_SNP_pval"]).abs().max()),
        ],
    })
    q_compare.to_csv(OUTDIR / "q_snp_consistency_check.tsv", sep="\t", index=False)

    q1_table = q1.rename(columns={
        "MAF": "FREQ",
    })
    q2_table = q2.rename(columns={
        "MAF": "FREQ",
    })
    q1_table.to_csv(OUTDIR / "ALPS_F1_Q_SNP_total_table.tsv.gz", sep="\t", index=False, compression="gzip")
    q2_table.to_csv(OUTDIR / "ALPS_F2_Q_SNP_total_table.tsv.gz", sep="\t", index=False, compression="gzip")

    summary = pd.DataFrame([
        {"metric": "n_total_snps", "value": len(q1_table)},
        {"metric": "F1_qsnp_sig_p_lt_5e-8", "value": int((q1_table["Q_SNP_pval"] < 5e-8).sum())},
        {"metric": "F2_qsnp_sig_p_lt_5e-8", "value": int((q2_table["Q_SNP_pval"] < 5e-8).sum())},
        {"metric": "F1_qsnp_suggestive_p_lt_1e-5", "value": int((q1_table["Q_SNP_pval"] < 1e-5).sum())},
        {"metric": "F2_qsnp_suggestive_p_lt_1e-5", "value": int((q2_table["Q_SNP_pval"] < 1e-5).sum())},
        {"metric": "F1_qsnp_lead_loci_500kb_proxy", "value": lead_loci_count(q1_table, "Q_SNP_pval")},
        {"metric": "F2_qsnp_lead_loci_500kb_proxy", "value": lead_loci_count(q2_table, "Q_SNP_pval")},
        {"metric": "F1_sig_hits_p_lt_5e-8", "value": int((f1["Pval_Estimate"] < 5e-8).sum())},
        {"metric": "F2_sig_hits_p_lt_5e-8", "value": int((f2["Pval_Estimate"] < 5e-8).sum())},
        {"metric": "F1_lead_loci_500kb_proxy", "value": lead_loci_count(f1, "Pval_Estimate")},
        {"metric": "F2_lead_loci_500kb_proxy", "value": lead_loci_count(f2, "Pval_Estimate")},
    ])
    summary.to_csv(OUTDIR / "q_snp_summary.tsv", sep="\t", index=False)

    overlap = pd.DataFrame([
        {
            "comparison": "F1_hits_vs_F1_QSNP_exact_snp_overlap",
            "n_overlap_sig_snps": int(((f1["Pval_Estimate"] < 5e-8) & (q1_table["Q_SNP_pval"] < 5e-8)).sum()),
        },
        {
            "comparison": "F2_hits_vs_F2_QSNP_exact_snp_overlap",
            "n_overlap_sig_snps": int(((f2["Pval_Estimate"] < 5e-8) & (q2_table["Q_SNP_pval"] < 5e-8)).sum()),
        },
        {
            "comparison": "F1_unique_sig_hits_excluding_500kb_from_F1_QSNP_sig",
            "n_overlap_sig_snps": exclude_near_q_hits(f1, q1_table, "Pval_Estimate", "Q_SNP_pval"),
        },
        {
            "comparison": "F2_unique_sig_hits_excluding_500kb_from_F2_QSNP_sig",
            "n_overlap_sig_snps": exclude_near_q_hits(f2, q2_table, "Pval_Estimate", "Q_SNP_pval"),
        },
    ])
    overlap.to_csv(OUTDIR / "factor_qsnp_overlap_summary.tsv", sep="\t", index=False)

    # Top tables
    q1_table.sort_values("Q_SNP_pval").head(1000).to_csv(
        OUTDIR / "ALPS_F1_Q_SNP_top1000.tsv", sep="\t", index=False
    )
    q2_table.sort_values("Q_SNP_pval").head(1000).to_csv(
        OUTDIR / "ALPS_F2_Q_SNP_top1000.tsv", sep="\t", index=False
    )

    print(f"Created: {OUTDIR / 'ALPS_F1_Q_SNP_total_table.tsv.gz'}")
    print(f"Created: {OUTDIR / 'ALPS_F2_Q_SNP_total_table.tsv.gz'}")
    print(f"Created: {OUTDIR / 'q_snp_summary.tsv'}")
    print(f"Created: {OUTDIR / 'factor_qsnp_overlap_summary.tsv'}")
    print(q_compare.to_string(index=False))
    print(summary.to_string(index=False))
    print(overlap.to_string(index=False))


if __name__ == "__main__":
    main()
