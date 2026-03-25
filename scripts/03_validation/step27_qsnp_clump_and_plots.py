from pathlib import Path
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(r"D:\文章\GS\GWAS\step11_factor_gwas_native_results\merged")
QDIR = ROOT / "q_snp_analysis"
PLOT_DIR = Path(r"D:\文章\GS\GWAS\publication_supplement\figures")
PLOT_DIR.mkdir(parents=True, exist_ok=True)

CLUMP_DIR = QDIR / "clump"
CLUMP_DIR.mkdir(parents=True, exist_ok=True)

REF_BFILE_WSL = "/mnt/d/GNOVA/GNOVA-master/EUR1000g/EUR.QC"

FILES = {
    "F1": ROOT / "ALPS_F1_factorGWAS_native_merged.tsv.gz",
    "F2": ROOT / "ALPS_F2_factorGWAS_native_merged.tsv.gz",
    "F1_Q": QDIR / "ALPS_F1_Q_SNP_total_table.tsv.gz",
    "F2_Q": QDIR / "ALPS_F2_Q_SNP_total_table.tsv.gz",
}


def to_wsl(p: Path) -> str:
    return "/mnt/" + p.drive[0].lower() + p.as_posix()[2:]


def prep_assoc(infile: Path, label: str, pcol: str) -> Path:
    df = pd.read_csv(infile, sep="\t", compression="gzip")
    assoc = df[["SNP", "Pval_Estimate" if pcol == "Pval_Estimate" else "Q_SNP_pval"]].copy()
    assoc.columns = ["SNP", "P"]
    assoc = assoc.dropna().drop_duplicates(subset=["SNP"])
    out = CLUMP_DIR / f"{label}.assoc"
    assoc.to_csv(out, sep="\t", index=False)
    return out


def run_clump(label: str, assoc_file: Path) -> Path:
    out_prefix = CLUMP_DIR / label
    cmd = [
        "wsl",
        "bash",
        "-lc",
        f"/usr/local/bin/plink --bfile {REF_BFILE_WSL} --clump {to_wsl(assoc_file)} "
        f"--clump-p1 5e-8 --clump-p2 1e-5 --clump-r2 0.1 --clump-kb 500 --out {to_wsl(out_prefix)}"
    ]
    subprocess.run(cmd, check=True)
    return out_prefix.with_suffix(".clumped")


def parse_clumped(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["analysis", "CHR", "SNP", "BP", "P"])
    df = pd.read_csv(path, delim_whitespace=True, comment="#")
    if df.empty:
        return pd.DataFrame(columns=["analysis", "CHR", "SNP", "BP", "P"])
    keep = [c for c in ["CHR", "SNP", "BP", "P"] if c in df.columns]
    df = df[keep].copy()
    df.insert(0, "analysis", label)
    return df


def manhattan_plot(df: pd.DataFrame, p_col: str, title: str, outfile: Path) -> None:
    plot_df = df[["CHR", "BP", p_col]].dropna().copy()
    plot_df = plot_df[(plot_df[p_col] > 0) & (plot_df["CHR"].between(1, 22))]
    plot_df["logp"] = -np.log10(plot_df[p_col])
    plot_df = plot_df.sort_values(["CHR", "BP"])
    offsets = {}
    current = 0
    ticks = []
    labels = []
    xs = []
    for chrom, sub in plot_df.groupby("CHR", sort=True):
        offsets[chrom] = current
        vals = sub["BP"] + current
        xs.extend(vals.tolist())
        ticks.append(vals.median())
        labels.append(str(chrom))
        current = vals.max() + 1_000_000
    plot_df["x"] = xs

    fig, ax = plt.subplots(figsize=(12, 4.8), dpi=300)
    colors = ["#3C6E71", "#D9A441"]
    for i, (chrom, sub) in enumerate(plot_df.groupby("CHR", sort=True)):
        ax.scatter(sub["x"], sub["logp"], s=4, color=colors[i % 2], alpha=0.75, linewidths=0)
    ax.axhline(-np.log10(5e-8), color="#B22222", linestyle="--", linewidth=1)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("-log10(P)")
    ax.set_xlabel("Chromosome")
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


def qq_plot(df: pd.DataFrame, p_col: str, title: str, outfile: Path) -> None:
    p = df[p_col].dropna()
    p = p[(p > 0) & (p <= 1)].sort_values().to_numpy()
    n = len(p)
    exp = -np.log10((np.arange(1, n + 1) - 0.5) / n)
    obs = -np.log10(p)
    m = max(exp.max(), obs.max())

    fig, ax = plt.subplots(figsize=(4.8, 4.8), dpi=300)
    ax.scatter(exp, obs, s=4, color="#2F4858", alpha=0.65, linewidths=0)
    ax.plot([0, m], [0, m], color="#B22222", linestyle="--", linewidth=1)
    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(outfile, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    assoc_files = {
        "F1": prep_assoc(FILES["F1"], "F1", "Pval_Estimate"),
        "F2": prep_assoc(FILES["F2"], "F2", "Pval_Estimate"),
        "F1_Q": prep_assoc(FILES["F1_Q"], "F1_Q", "Q_SNP_pval"),
        "F2_Q": prep_assoc(FILES["F2_Q"], "F2_Q", "Q_SNP_pval"),
    }

    clumped = {}
    for label, assoc in assoc_files.items():
        clumped[label] = parse_clumped(run_clump(label, assoc), label)
        clumped[label].to_csv(CLUMP_DIR / f"{label}_lead_loci.tsv", sep="\t", index=False)

    summary = pd.DataFrame([
        {"analysis": label, "n_lead_loci_clumped": len(df)}
        for label, df in clumped.items()
    ])

    f1_leads = set(clumped["F1"]["SNP"])
    f2_leads = set(clumped["F2"]["SNP"])
    f1q_leads = set(clumped["F1_Q"]["SNP"])
    f2q_leads = set(clumped["F2_Q"]["SNP"])
    overlap = pd.DataFrame([
        {"comparison": "F1_leads_vs_F1Q_leads_exact", "n_overlap": len(f1_leads & f1q_leads)},
        {"comparison": "F2_leads_vs_F2Q_leads_exact", "n_overlap": len(f2_leads & f2q_leads)},
        {"comparison": "F1_unique_leads_excluding_F1Q_leads", "n_overlap": len(f1_leads - f1q_leads)},
        {"comparison": "F2_unique_leads_excluding_F2Q_leads", "n_overlap": len(f2_leads - f2q_leads)},
    ])
    summary.to_csv(CLUMP_DIR / "clumped_lead_loci_summary.tsv", sep="\t", index=False)
    overlap.to_csv(CLUMP_DIR / "clumped_factor_qsnp_overlap.tsv", sep="\t", index=False)

    f1q = pd.read_csv(FILES["F1_Q"], sep="\t", compression="gzip")
    f2q = pd.read_csv(FILES["F2_Q"], sep="\t", compression="gzip")
    manhattan_plot(f1q, "Q_SNP_pval", "F1 Q_SNP Manhattan Plot", PLOT_DIR / "Fig_S8_F1_QSNP_manhattan.png")
    manhattan_plot(f2q, "Q_SNP_pval", "F2 Q_SNP Manhattan Plot", PLOT_DIR / "Fig_S9_F2_QSNP_manhattan.png")
    qq_plot(f1q, "Q_SNP_pval", "F1 Q_SNP QQ Plot", PLOT_DIR / "Fig_S10_F1_QSNP_QQ.png")
    qq_plot(f2q, "Q_SNP_pval", "F2 Q_SNP_pval QQ Plot", PLOT_DIR / "Fig_S11_F2_QSNP_QQ.png")

    print(summary.to_string(index=False))
    print(overlap.to_string(index=False))


if __name__ == "__main__":
    main()
