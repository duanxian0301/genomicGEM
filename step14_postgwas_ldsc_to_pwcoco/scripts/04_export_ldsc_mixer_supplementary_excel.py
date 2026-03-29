from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


BASE_POSTGWAS = Path(r"D:\文章\GS\postgwas")
LDSC_MAIN = BASE_POSTGWAS / "01_ldsc_f1f2_ndd"
LDSC_EXT = BASE_POSTGWAS / "01b_ldsc_alps_family_vs_ndd"
MIXER_SUMMARY = BASE_POSTGWAS / "02_mixer_f1f2_ndd"
MIXER_RUN = Path(r"D:\codex\GenomicSEM\postgwas_mixer_f1f2_ndd\runs\rep1")
OUTPUT_XLSX = BASE_POSTGWAS / "Supplementary_LDSC_MiXeR_tables.xlsx"


def read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def markdown_to_df(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8").splitlines()
    return pd.DataFrame({"line_no": range(1, len(lines) + 1), "text": lines})


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


def main() -> None:
    sheets: list[tuple[str, pd.DataFrame]] = []

    sheets.append(("README", pd.DataFrame(
        {
            "item": [
                "Workbook",
                "Purpose",
                "Output file",
                "Generated from",
                "Contains",
            ],
            "value": [
                "Supplementary_LDSC_MiXeR_tables.xlsx",
                "Integrated supplementary tables for LDSC and MiXeR results",
                str(OUTPUT_XLSX),
                str(Path(__file__)),
                "LDSC main, LDSC extended, MiXeR summary, MiXeR univariate, MiXeR bivariate, source manifests, text summaries",
            ],
        }
    )))

    # LDSC main
    sheets.append(("LDSC_main_pairs", read_tsv(LDSC_MAIN / "f1f2_vs_ndd_requested_pairs.tsv")))
    sheets.append(("LDSC_main_h2", read_tsv(LDSC_MAIN / "f1f2_vs_ndd_ldsc_h2.tsv")))
    sheets.append(("LDSC_main_rg_matrix", read_tsv(LDSC_MAIN / "f1f2_vs_ndd_ldsc_rg_matrix.tsv")))
    sheets.append(("LDSC_main_manifest", read_tsv(LDSC_MAIN / "input_manifest.tsv")))
    sheets.append(("LDSC_main_summary", markdown_to_df(LDSC_MAIN / "LDSC_summary.md")))

    # LDSC extended
    sheets.append(("LDSC_ext_compare", read_tsv(LDSC_EXT / "alps_family_vs_ndd_comparison_table.tsv")))
    sheets.append(("LDSC_ext_pairs", read_tsv(LDSC_EXT / "alps_family_vs_ndd_requested_pairs.tsv")))
    sheets.append(("LDSC_ext_h2", read_tsv(LDSC_EXT / "alps_family_vs_ndd_ldsc_h2.tsv")))
    sheets.append(("LDSC_ext_rg_matrix", read_tsv(LDSC_EXT / "alps_family_vs_ndd_ldsc_rg_matrix.tsv")))
    sheets.append(("LDSC_ext_AD_rank", read_tsv(LDSC_EXT / "AD_ranked_rg.tsv")))
    sheets.append(("LDSC_ext_PD_rank", read_tsv(LDSC_EXT / "PD_ranked_rg.tsv")))
    sheets.append(("LDSC_ext_LBD_rank", read_tsv(LDSC_EXT / "LBD_ranked_rg.tsv")))
    sheets.append(("LDSC_ext_manifest", read_tsv(LDSC_EXT / "input_manifest.tsv")))
    sheets.append(("LDSC_ext_summary", markdown_to_df(LDSC_EXT / "LDSC_comparison_summary.md")))

    # MiXeR summary
    sheets.append(("MiXeR_summary", read_tsv(MIXER_SUMMARY / "mixer_f1f2_ndd_summary.tsv")))
    sheets.append(("MiXeR_summary_text", markdown_to_df(MIXER_SUMMARY / "MiXeR_summary.md")))

    univariate_files = sorted(MIXER_RUN.glob("*.fit.rep1.json"))
    univariate_files = [p for p in univariate_files if "_vs_" not in p.name]
    bivariate_files = sorted(MIXER_RUN.glob("*_vs_*.test.rep1.json"))

    univariate_df = pd.DataFrame(parse_mixer_univariate(path) for path in univariate_files).sort_values("trait")
    bivariate_df = pd.DataFrame(parse_mixer_bivariate(path) for path in bivariate_files).sort_values(["trait2", "trait1"])

    sheets.append(("MiXeR_univariate", univariate_df))
    sheets.append(("MiXeR_bivariate", bivariate_df))
    sheets.append((
        "MiXeR_files",
        pd.DataFrame(
            {
                "file": [str(p) for p in sorted(MIXER_RUN.glob("*.json"))],
                "size_bytes": [p.stat().st_size for p in sorted(MIXER_RUN.glob("*.json"))],
            }
        ),
    ))

    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
        for sheet_name, df in sheets:
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

        for ws in writer.book.worksheets:
            autosize_worksheet(ws)
            ws.freeze_panes = "A2"

    print(f"Wrote workbook to: {OUTPUT_XLSX}")


if __name__ == "__main__":
    main()
