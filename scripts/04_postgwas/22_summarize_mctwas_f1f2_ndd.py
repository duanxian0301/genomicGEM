from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


PAIR_DEFS = [
    ("F1", "AD"),
    ("F2", "AD"),
    ("F1", "PD"),
    ("F2", "PD"),
    ("F1", "LBD"),
    ("F2", "LBD"),
]


def load_trait_table(trait_dir: Path, trait: str) -> pd.DataFrame:
    path = trait_dir / trait / "results" / f"{trait}_finemap_annotated.tsv"
    if not path.exists():
        raise FileNotFoundError(f"Missing result file: {path}")
    df = pd.read_csv(path, sep="\t")
    df["trait"] = trait
    df["source_file"] = str(path)
    if "feature" in df.columns:
        df["gene"] = df["feature"].astype(str).str.replace(r"_expression$", "", regex=True)
    else:
        df["gene"] = pd.NA
    if "pip" not in df.columns:
        df["pip"] = pd.NA
    return df


def choose_p_col(df: pd.DataFrame) -> str | None:
    for col in ("p", "pvalue", "p_value", "P", "pval"):
        if col in df.columns:
            return col
    return None


def summarize_pairs(all_traits: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[dict] = []
    for left, right in PAIR_DEFS:
        left_df = all_traits[left].copy()
        right_df = all_traits[right].copy()
        merged = left_df.merge(
            right_df,
            on="gene",
            how="outer",
            suffixes=(f"_{left}", f"_{right}"),
        )
        for _, rec in merged.iterrows():
            rows.append(
                {
                    "pair": f"{left}_vs_{right}",
                    "gene / feature": rec.get("gene"),
                    "tissue / model": "; ".join(
                        sorted(
                            {
                                str(x)
                                for x in [rec.get(f"tissue_{left}"), rec.get(f"tissue_{right}")]
                                if pd.notna(x) and str(x) != ""
                            }
                        )
                    )
                    or pd.NA,
                    "z / effect": "; ".join(
                        [
                            f"{left}:{rec.get(f'z_{left}')}" if pd.notna(rec.get(f"z_{left}")) else "",
                            f"{right}:{rec.get(f'z_{right}')}" if pd.notna(rec.get(f"z_{right}")) else "",
                        ]
                    ).strip("; "),
                    "p": "; ".join(
                        [
                            f"{left}:{rec.get(f'p_{left}')}" if pd.notna(rec.get(f"p_{left}")) else "",
                            f"{right}:{rec.get(f'p_{right}')}" if pd.notna(rec.get(f"p_{right}")) else "",
                        ]
                    ).strip("; "),
                    "FDR": "; ".join(
                        [
                            f"{left}:{rec.get(f'FDR_{left}')}" if pd.notna(rec.get(f"FDR_{left}")) else "",
                            f"{right}:{rec.get(f'FDR_{right}')}" if pd.notna(rec.get(f"FDR_{right}")) else "",
                        ]
                    ).strip("; "),
                    "whether significant": "; ".join(
                        [
                            f"{left}:{rec.get(f'whether significant_{left}')}"
                            if pd.notna(rec.get(f"whether significant_{left}"))
                            else "",
                            f"{right}:{rec.get(f'whether significant_{right}')}"
                            if pd.notna(rec.get(f"whether significant_{right}"))
                            else "",
                        ]
                    ).strip("; "),
                    "source files": "; ".join(
                        [
                            str(rec.get(f"source_file_{left}")) if pd.notna(rec.get(f"source_file_{left}")) else "",
                            str(rec.get(f"source_file_{right}")) if pd.notna(rec.get(f"source_file_{right}")) else "",
                        ]
                    ).strip("; "),
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()

    root = Path(args.output_root)
    single_trait_dir = root / "single_trait"
    summary_dir = root / "summary"
    summary_dir.mkdir(parents=True, exist_ok=True)

    all_traits: dict[str, pd.DataFrame] = {}
    single_trait_rows: list[pd.DataFrame] = []

    for trait in ["F1", "F2", "AD", "PD", "LBD"]:
        df = load_trait_table(single_trait_dir, trait)
        p_col = choose_p_col(df)
        if p_col is not None and "p" not in df.columns:
            df["p"] = df[p_col]
        if "pip" in df.columns and "whether significant" not in df.columns:
            df["whether significant"] = df["pip"].fillna(0).ge(0.5)
        if "FDR" not in df.columns:
            df["FDR"] = pd.NA
        all_traits[trait] = df
        single_trait_rows.append(df)

    single_trait_df = pd.concat(single_trait_rows, ignore_index=True)
    single_trait_out = summary_dir / "mctwas_single_trait_results.tsv"
    single_trait_df.to_csv(single_trait_out, sep="\t", index=False)

    pair_df = summarize_pairs(all_traits)
    pair_out = summary_dir / "mctwas_f1f2_ndd_summary.tsv"
    pair_df.to_csv(pair_out, sep="\t", index=False)

    lines = [
        "# M-cTWAS summary",
        "",
        "Analysis scheme: official single-trait cTWAS with 9 Broad GTEx brain tissues in hg19/GRCh37, followed by pairwise comparison across F1/F2 and AD/PD/LBD.",
        "",
        "Reference weights: Broad FUSION GTEx brain hg19 panel (9 tissues).",
        "",
        "Output tables:",
        f"- {single_trait_out}",
        f"- {pair_out}",
    ]
    (summary_dir / "mctwas_f1f2_ndd_summary.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
