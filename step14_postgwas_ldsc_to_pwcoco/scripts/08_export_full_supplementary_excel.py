from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


BASE_POSTGWAS = Path(r"D:\文章\GS\postgwas")
LDSC_MAIN = BASE_POSTGWAS / "01_ldsc_f1f2_ndd"
LDSC_EXT = BASE_POSTGWAS / "01b_ldsc_alps_family_vs_ndd"
MIXER_SUMMARY = BASE_POSTGWAS / "02_mixer_f1f2_ndd"
MIXER_RUN = Path(r"D:\codex\GenomicSEM\postgwas_mixer_f1f2_ndd\runs\rep1")
PLEIO_BASE = BASE_POSTGWAS / "03_pleiofdr_all_pairs"
PLEIO_RESULTS = PLEIO_BASE / "results"
PLEIO_ALPS_BASE = BASE_POSTGWAS / "03b_pleiofdr_alps_family_vs_ndd"
PLEIO_ALPS_RESULTS = PLEIO_ALPS_BASE / "results"
COLOC_BASE = BASE_POSTGWAS / "07_coloc_factor_ndd"
PWCOCO_BASE = BASE_POSTGWAS / "08_pwcoco"
OUTPUT_XLSX = BASE_POSTGWAS / "Supplementary_LDSC_MiXeR_pleioFDR_tables.xlsx"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def markdown_to_df(path: Path, section: str) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8").splitlines()
    return pd.DataFrame({"section": section, "line_no": range(1, len(lines) + 1), "text": lines})


def parse_mixer_univariate(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8"))
    trait = path.name.replace(".fit.rep1.json", "")
    ci = obj["ci"]
    params = obj["params"]
    return {
        "trait": trait,
        "analysis_file": str(path),
        "n": obj["options"].get("trait1_nval"),
        "num_snp": obj["options"].get("num_snp"),
        "num_tag": obj["options"].get("num_tag"),
        "sum_weights": obj["options"].get("sum_weights"),
        "pi": ci["pi"]["point_estimate"],
        "nc": ci["nc"]["point_estimate"],
        "nc_p9": ci["nc@p9"]["point_estimate"],
        "sig2_beta": ci["sig2_beta"]["point_estimate"],
        "sig2_zero": ci["sig2_zero"]["point_estimate"],
        "h2": ci["h2"]["point_estimate"],
        "param_pi": params["pi"],
        "param_sig2_beta": params["sig2_beta"],
        "param_sig2_zero": params["sig2_zero"],
    }


def parse_mixer_bivariate(path: Path) -> dict:
    obj = json.loads(path.read_text(encoding="utf-8"))
    pair = path.name.replace(".test.rep1.json", "")
    ci = obj["ci"]
    params = obj["params"]
    trait1 = Path(obj["options"]["trait1_file"]).stem
    trait2 = Path(obj["options"]["trait2_file"]).stem
    return {
        "pair": pair.replace("_vs_", "-"),
        "analysis_file": str(path),
        "trait1": trait1,
        "trait2": trait2,
        "n_trait1": obj["options"].get("trait1_nval"),
        "n_trait2": obj["options"].get("trait2_nval"),
        "num_snp": obj["options"].get("num_snp"),
        "num_tag": obj["options"].get("num_tag"),
        "sum_weights": obj["options"].get("sum_weights"),
        "rg": ci["rg"]["point_estimate"],
        "dice": ci["dice"]["point_estimate"],
        "pi1": ci["pi1"]["point_estimate"],
        "pi2": ci["pi2"]["point_estimate"],
        "pi12": ci["pi12"]["point_estimate"],
        "pi1u": ci["pi1u"]["point_estimate"],
        "pi2u": ci["pi2u"]["point_estimate"],
        "pi12_over_pi1u": ci["pi12_over_pi1u"]["point_estimate"],
        "pi12_over_pi2u": ci["pi12_over_pi2u"]["point_estimate"],
        "pi1_over_totalpi": ci["pi1_over_totalpi"]["point_estimate"],
        "pi2_over_totalpi": ci["pi2_over_totalpi"]["point_estimate"],
        "pi12_over_totalpi": ci["pi12_over_totalpi"]["point_estimate"],
        "nc1": ci["nc1"]["point_estimate"],
        "nc2": ci["nc2"]["point_estimate"],
        "nc12": ci["nc12"]["point_estimate"],
        "totalnc": ci["totalnc"]["point_estimate"],
        "rho_beta": ci["rho_beta"]["point_estimate"],
        "rho_zero": ci["rho_zero"]["point_estimate"],
        "sig2_zero_t1": ci["sig2_zero_T1"]["point_estimate"],
        "sig2_zero_t2": ci["sig2_zero_T2"]["point_estimate"],
        "sig2_beta_t1": ci["sig2_beta_T1"]["point_estimate"],
        "sig2_beta_t2": ci["sig2_beta_T2"]["point_estimate"],
        "h2_t1": ci["h2_T1"]["point_estimate"],
        "h2_t2": ci["h2_T2"]["point_estimate"],
        "param_pi1": params["pi"][0],
        "param_pi2": params["pi"][1],
        "param_pi12": params["pi"][2],
        "param_rho_beta": params["rho_beta"],
        "param_rho_zero": params["rho_zero"],
    }


def autosize_worksheet(ws) -> None:
    for column_cells in ws.columns:
        values = [str(cell.value) if cell.value is not None else "" for cell in column_cells]
        max_len = max((len(v) for v in values), default=0)
        ws.column_dimensions[column_cells[0].column_letter].width = min(max(max_len + 2, 10), 40)


def pair_files(base_results: Path) -> list[tuple[str, Path, Path]]:
    out = []
    for pair_dir in sorted(base_results.iterdir()):
        if not pair_dir.is_dir():
            continue
        all_file = next((p for p in pair_dir.glob("*conjfdr_0.05_all.csv") if "zscore" not in p.name), None)
        loci_file = next((p for p in pair_dir.glob("*conjfdr_0.05_loci.csv") if "zscore" not in p.name), None)
        if all_file and loci_file:
            out.append((pair_dir.name, all_file, loci_file))
    return out


def combine_markdown_notes() -> pd.DataFrame:
    dfs = [
        markdown_to_df(LDSC_MAIN / "LDSC_summary.md", "LDSC_main"),
        markdown_to_df(LDSC_EXT / "LDSC_extended", "LDSC_extended")
        if (LDSC_EXT / "LDSC_extended").exists()
        else markdown_to_df(LDSC_EXT / "LDSC_comparison_summary.md", "LDSC_extended"),
        markdown_to_df(MIXER_SUMMARY / "MiXeR_summary.md", "MiXeR"),
        markdown_to_df(PLEIO_BASE / "pleiofdr_notes.md", "pleioFDR_factor"),
        markdown_to_df(PLEIO_ALPS_BASE / "pleiofdr_alps_family_notes.md", "pleioFDR_ALPS_family"),
        markdown_to_df(COLOC_BASE / "README.md", "coloc"),
        markdown_to_df(PWCOCO_BASE / "pwcoco_region_best_summary.tsv", "PWCoCo_best_summary"),
    ]
    return pd.concat(dfs, ignore_index=True)


def combine_ldsc_extended_rankings() -> pd.DataFrame:
    frames = []
    for disease in ["AD", "PD", "LBD"]:
        df = read_tsv(LDSC_EXT / f"{disease}_ranked_rg.tsv")
        if "disease" not in df.columns:
            df.insert(0, "disease", disease)
        frames.append(df)
    return pd.concat(frames, ignore_index=True)


def combine_pleio_results(base_results: Path, family_label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_frames = []
    loci_frames = []
    for pair_name, all_file, loci_file in pair_files(base_results):
        all_df = read_csv(all_file)
        all_df.insert(0, "pair", pair_name)
        all_df.insert(1, "analysis_family", family_label)
        all_frames.append(all_df)

        loci_df = read_csv(loci_file)
        loci_df.insert(0, "pair", pair_name)
        loci_df.insert(1, "analysis_family", family_label)
        loci_frames.append(loci_df)
    return pd.concat(all_frames, ignore_index=True), pd.concat(loci_frames, ignore_index=True)


def combine_coloc_snp_results() -> pd.DataFrame:
    per_task_dir = COLOC_BASE / "per_task"
    frames = []
    for path in sorted(per_task_dir.glob("*_snp_results.tsv.gz")):
        frames.append(read_tsv(path))
    return pd.concat(frames, ignore_index=True)


def main() -> None:
    sheets: list[tuple[str, pd.DataFrame]] = []

    readme = pd.DataFrame(
        {
            "sheet_name": [
                "README",
                "Notes",
                "LDSC_main_pairs",
                "LDSC_main_h2",
                "LDSC_ext_compare",
                "LDSC_ext_h2",
                "LDSC_ext_rank",
                "MiXeR_bivariate",
                "MiXeR_univariate",
                "pleio_factor_summary",
                "pleio_factor_loci",
                "pleio_alps_summary",
                "pleio_alps_compare",
                "pleio_alps_rank",
                "pleio_alps_loci",
                "coloc_summary",
                "pwcoco_best",
                "pwcoco_all",
            ],
            "description": [
                "Workbook overview and sheet guide",
                "Integrated method/result notes extracted from markdown summaries",
                "Primary LDSC pairwise results for F1/F2 vs AD/PD/LBD",
                "Primary LDSC heritability diagnostics",
                "Extended LDSC comparison across F1/F2, original ALPS, Mean_ALPS, and tALPS vs NDDs",
                "Extended LDSC heritability diagnostics",
                "Extended LDSC disease-specific rankings combined into one table",
                "MiXeR bivariate results for F1/F2 vs AD/PD/LBD",
                "MiXeR univariate trait architecture estimates",
                "Factor-level pleioFDR summary table",
                "Factor-level pleioFDR loci results combined across all pairs",
                "ALPS-family pleioFDR summary table",
                "Factor vs ALPS-family pleioFDR comparison table",
                "ALPS-family disease-specific rankings combined into one table",
                "ALPS-family pleioFDR loci results combined across all pairs",
                "coloc summary across all candidate regions",
                "PWCoCo best-summary table across tested regions",
                "PWCoCo all tested models across tested regions",
            ],
        }
    )
    sheets.append(("README", readme))

    sheets.append(("Notes", combine_markdown_notes()))

    sheets.append(("LDSC_main_pairs", read_tsv(LDSC_MAIN / "f1f2_vs_ndd_requested_pairs.tsv")))
    sheets.append(("LDSC_main_h2", read_tsv(LDSC_MAIN / "f1f2_vs_ndd_ldsc_h2.tsv")))
    sheets.append(("LDSC_ext_compare", read_tsv(LDSC_EXT / "alps_family_vs_ndd_comparison_table.tsv")))
    sheets.append(("LDSC_ext_h2", read_tsv(LDSC_EXT / "alps_family_vs_ndd_ldsc_h2.tsv")))
    sheets.append(("LDSC_ext_rank", combine_ldsc_extended_rankings()))

    univariate_files = sorted([p for p in MIXER_RUN.glob("*.fit.rep1.json") if "_vs_" not in p.name])
    bivariate_files = sorted(MIXER_RUN.glob("*_vs_*.test.rep1.json"))
    sheets.append(("MiXeR_bivariate", pd.DataFrame(parse_mixer_bivariate(p) for p in bivariate_files).sort_values(["trait2", "trait1"])))
    sheets.append(("MiXeR_univariate", pd.DataFrame(parse_mixer_univariate(p) for p in univariate_files).sort_values("trait")))

    pleio_factor_snps, pleio_factor_loci = combine_pleio_results(PLEIO_RESULTS, "factor")
    sheets.append(("pleio_factor_summary", read_tsv(PLEIO_BASE / "pleiofdr_pair_summary.tsv")))
    sheets.append(("pleio_factor_loci", pleio_factor_loci))

    pleio_alps_snps, pleio_alps_loci = combine_pleio_results(PLEIO_ALPS_RESULTS, "alps_family")
    alps_rank_frames = []
    for disease in ["AD", "PD", "LBD"]:
        df = read_tsv(PLEIO_ALPS_BASE / f"{disease}_factor_vs_alps_family_ranked.tsv")
        if "disease" not in df.columns:
            df.insert(0, "disease", disease)
        alps_rank_frames.append(df)
    sheets.append(("pleio_alps_summary", read_tsv(PLEIO_ALPS_BASE / "pleiofdr_alps_family_summary.tsv")))
    sheets.append(("pleio_alps_compare", read_tsv(PLEIO_ALPS_BASE / "pleiofdr_factor_vs_alps_family_comparison.tsv")))
    sheets.append(("pleio_alps_rank", pd.concat(alps_rank_frames, ignore_index=True)))
    sheets.append(("pleio_alps_loci", pleio_alps_loci))

    sheets.append(("coloc_summary", read_tsv(COLOC_BASE / "coloc_factor_ndd_summary.tsv")))
    sheets.append(("pwcoco_best", read_tsv(PWCOCO_BASE / "pwcoco_region_best_summary.tsv")))
    sheets.append(("pwcoco_all", read_tsv(PWCOCO_BASE / "pwcoco_region_summary.tsv")))

    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
        for sheet_name, df in sheets:
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)
        for ws in writer.book.worksheets:
            autosize_worksheet(ws)
            ws.freeze_panes = "A2"

    print(f"Wrote workbook to: {OUTPUT_XLSX}")
    print(f"Total sheets: {len(sheets)}")


if __name__ == "__main__":
    main()
