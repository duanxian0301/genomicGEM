from pathlib import Path

import pandas as pd


SUPP_DIR = Path(r"D:\文章\GS\GWAS\publication_supplement")
BASE_XLSX = SUPP_DIR / "ALPS_GenomicSEM_Supplementary_Tables_updated_v3_usermodel.xlsx"
OUT_XLSX = SUPP_DIR / "ALPS_GenomicSEM_Supplementary_Tables_updated_v4_qsnp.xlsx"

QDIR = Path(r"D:\文章\GS\GWAS\step11_factor_gwas_native_results\merged\q_snp_analysis")
CLUMP_DIR = QDIR / "clump"


def main() -> None:
    xls = pd.ExcelFile(BASE_XLSX)
    q_summary = pd.read_csv(QDIR / "q_snp_summary.tsv", sep="\t")
    q_overlap_proxy = pd.read_csv(QDIR / "factor_qsnp_overlap_summary.tsv", sep="\t")
    q_consistency = pd.read_csv(QDIR / "q_snp_consistency_check.tsv", sep="\t")
    clump_summary = pd.read_csv(CLUMP_DIR / "clumped_lead_loci_summary.tsv", sep="\t")
    clump_overlap = pd.read_csv(CLUMP_DIR / "clumped_factor_qsnp_overlap.tsv", sep="\t")
    f1_q_leads = pd.read_csv(CLUMP_DIR / "F1_Q_lead_loci.tsv", sep="\t")
    f2_q_leads = pd.read_csv(CLUMP_DIR / "F2_Q_lead_loci.tsv", sep="\t")

    with pd.ExcelWriter(OUT_XLSX, engine="openpyxl") as writer:
        for sheet in xls.sheet_names:
            df = pd.read_excel(BASE_XLSX, sheet_name=sheet)
            df.to_excel(writer, sheet_name=sheet, index=False)

        q_summary.to_excel(writer, sheet_name="S19_QSNP_summary", index=False)
        clump_summary.to_excel(writer, sheet_name="S20_QSNP_clumped_leads", index=False)
        clump_overlap.to_excel(writer, sheet_name="S21_QSNP_factor_overlap", index=False)
        q_overlap_proxy.to_excel(writer, sheet_name="S22_QSNP_overlap_proxy", index=False)
        q_consistency.to_excel(writer, sheet_name="S23_QSNP_consistency", index=False)
        f1_q_leads.to_excel(writer, sheet_name="S24_F1_QSNP_leads", index=False)
        f2_q_leads.to_excel(writer, sheet_name="S25_F2_QSNP_leads", index=False)

    print(f"Created workbook: {OUT_XLSX}")


if __name__ == "__main__":
    main()
