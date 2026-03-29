from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image, ImageChops


BASE = Path(r"D:\文章\GS\postgwas")
FIG_BASE = BASE / "04_integrated_figures"
MAIN_DIR = FIG_BASE / "Main_Figures"
SUPP_DIR = FIG_BASE / "Supplementary_Figures"
SUPP_LDSC = SUPP_DIR / "LDSC"
SUPP_MIXER = SUPP_DIR / "MiXeR"
SUPP_PLEIO = SUPP_DIR / "pleioFDR"
SUPP_PLEIO_ALPS = SUPP_PLEIO / "ALPS_family_vs_NDD"
SUPP_COLOC = SUPP_DIR / "coloc"
SUPP_PWCOCO = SUPP_DIR / "PWCoCo"


def ensure_dirs() -> None:
    for p in [MAIN_DIR, SUPP_DIR, SUPP_LDSC, SUPP_MIXER, SUPP_PLEIO, SUPP_PLEIO_ALPS, SUPP_COLOC, SUPP_PWCOCO]:
        p.mkdir(parents=True, exist_ok=True)


def copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        shutil.copy2(src, dst)


def make_ldsc_mixer_main() -> None:
    ldsc = pd.read_csv(BASE / "01_ldsc_f1f2_ndd" / "f1f2_vs_ndd_requested_pairs.tsv", sep="\t")
    mixer = pd.read_csv(BASE / "02_mixer_f1f2_ndd" / "mixer_f1f2_ndd_summary.tsv", sep="\t")
    mixer["pair"] = mixer["pair"].str.replace("_", "-", regex=False)

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    colors = ["#0b6e4f" if "AD" in p else "#c97b26" if "PD" in p else "#607d8b" for p in ldsc["trait2"]]
    labels = ldsc["trait1"] + "-" + ldsc["trait2"]

    axes[0].bar(range(len(ldsc)), ldsc["rg"], yerr=ldsc["rg_se"], color=colors, capsize=3)
    axes[0].axhline(0, color="black", linewidth=0.8)
    axes[0].set_xticks(range(len(ldsc)))
    axes[0].set_xticklabels(labels, rotation=45, ha="right")
    axes[0].set_title("LDSC rg")
    axes[0].set_ylabel("Genetic correlation")

    axes[1].bar(range(len(mixer)), mixer["rg"], color=colors, alpha=0.9)
    axes[1].axhline(0, color="black", linewidth=0.8)
    axes[1].set_xticks(range(len(mixer)))
    axes[1].set_xticklabels(mixer["pair"], rotation=45, ha="right")
    axes[1].set_title("MiXeR rg")
    axes[1].set_ylabel("Polygenic overlap rg")

    axes[2].bar(range(len(mixer)), mixer["pi12_over_pi2u"], color=colors, alpha=0.9)
    axes[2].set_xticks(range(len(mixer)))
    axes[2].set_xticklabels(mixer["pair"], rotation=45, ha="right")
    axes[2].set_ylim(0, 1)
    axes[2].set_title("MiXeR shared fraction in disease")
    axes[2].set_ylabel("pi12 / pi2u")

    fig.suptitle("Global Genetic Overlap of F1/F2 with AD, PD, and LBD", fontsize=14)
    fig.tight_layout()
    fig.savefig(MAIN_DIR / "Figure1_LDSC_MiXeR_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_ldsc_supp() -> None:
    ext = pd.read_csv(BASE / "01b_ldsc_alps_family_vs_ndd" / "alps_family_vs_ndd_comparison_table.tsv", sep="\t")
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharey=True)
    diseases = ["AD", "PD", "LBD"]
    palette = {"factor": "#0b6e4f", "original": "#c97b26", "summary": "#607d8b"}

    for ax, disease in zip(axes, diseases):
        sub = ext[ext["trait2"] == disease].copy().sort_values("rg")
        ax.barh(sub["trait1"], sub["rg"], color=[palette.get(p, "#999999") for p in sub["panel"]])
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title(disease)
        ax.set_xlabel("rg")
    axes[0].set_ylabel("Trait")
    fig.suptitle("Extended LDSC comparison across ALPS-derived traits", fontsize=14)
    fig.tight_layout()
    fig.savefig(SUPP_LDSC / "Supplementary_LDSC_extended_comparison.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_mixer_supp() -> None:
    mixer = pd.read_csv(BASE / "02_mixer_f1f2_ndd" / "mixer_f1f2_ndd_summary.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = ["#0b6e4f" if "AD" in p else "#c97b26" if "PD" in p else "#607d8b" for p in mixer["pair"]]

    axes[0].scatter(mixer["rg"], mixer["dice"], s=120, c=colors)
    for _, row in mixer.iterrows():
        axes[0].text(row["rg"], row["dice"], row["pair"], fontsize=8, ha="left", va="bottom")
    axes[0].set_xlabel("rg")
    axes[0].set_ylabel("dice")
    axes[0].set_title("MiXeR rg vs dice")

    axes[1].bar(range(len(mixer)), mixer["nc12"], color=colors)
    axes[1].set_xticks(range(len(mixer)))
    axes[1].set_xticklabels(mixer["pair"], rotation=45, ha="right")
    axes[1].set_ylabel("nc12")
    axes[1].set_title("Estimated shared causal component")

    fig.tight_layout()
    fig.savefig(SUPP_MIXER / "Supplementary_MiXeR_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def load_png(path: Path) -> Image.Image:
    return Image.open(path).convert("RGB")


def crop_whitespace(img: Image.Image, pad: int = 16) -> Image.Image:
    bg = Image.new(img.mode, img.size, "white")
    diff = ImageChops.difference(img, bg)
    bbox = diff.getbbox()
    if bbox is None:
        return img
    left = max(bbox[0] - pad, 0)
    top = max(bbox[1] - pad, 0)
    right = min(bbox[2] + pad, img.size[0])
    bottom = min(bbox[3] + pad, img.size[1])
    return img.crop((left, top, right, bottom))


def make_pleio_main() -> None:
    plt.rcParams.update({"font.family": "DejaVu Serif", "figure.facecolor": "white"})
    fig = plt.figure(figsize=(14, 10), facecolor="white")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.65, 1.0], hspace=0.16, wspace=0.08)

    panel_specs = [
        ("A", "F2-AD Manhattan", BASE / "03_pleiofdr_all_pairs" / "results" / "F2_AD" / "F2_AD_conjfdr_0.05_manhattan.png"),
        ("B", "F2-AD Conditional QQ", BASE / "03_pleiofdr_all_pairs" / "results" / "F2_AD" / "F2_vs_AD_qq.png"),
        ("C", "F1-PD Manhattan", BASE / "03_pleiofdr_all_pairs" / "results" / "F1_PD" / "F1_PD_conjfdr_0.05_manhattan.png"),
        ("D", "F1-PD Conditional QQ", BASE / "03_pleiofdr_all_pairs" / "results" / "F1_PD" / "F1_vs_PD_qq.png"),
    ]

    for idx, (label, title, path) in enumerate(panel_specs):
        ax = fig.add_subplot(gs[idx // 2, idx % 2])
        ax.imshow(crop_whitespace(load_png(path)))
        ax.set_axis_off()
        ax.text(0.0, 1.03, label, transform=ax.transAxes, fontsize=16, fontweight="bold", ha="left", va="bottom")
        ax.text(0.07, 1.03, title, transform=ax.transAxes, fontsize=12, ha="left", va="bottom")

    fig.savefig(MAIN_DIR / "Figure2_pleioFDR_main.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def archive_pleio_figures() -> None:
    factor_base = BASE / "03_pleiofdr_all_pairs" / "results"
    for pair_dir in sorted(factor_base.iterdir()):
        if not pair_dir.is_dir():
            continue
        dst = SUPP_PLEIO / pair_dir.name
        dst.mkdir(parents=True, exist_ok=True)
        for ext in ("*.png", "*.svg", "*.csv", "*.log"):
            for src in pair_dir.glob(ext):
                copy_if_exists(src, dst / src.name)

    alps_base = BASE / "03b_pleiofdr_alps_family_vs_ndd" / "results"
    if alps_base.exists():
        for pair_dir in sorted(alps_base.iterdir()):
            if not pair_dir.is_dir():
                continue
            dst = SUPP_PLEIO_ALPS / pair_dir.name
            dst.mkdir(parents=True, exist_ok=True)
            for ext in ("*.png", "*.svg", "*.csv", "*.log"):
                for src in pair_dir.glob(ext):
                    copy_if_exists(src, dst / src.name)


def make_pleio_supp_summary() -> None:
    df = pd.read_csv(BASE / "03_pleiofdr_all_pairs" / "pleiofdr_pair_summary.tsv", sep="\t")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    colors = ["#0b6e4f" if "AD" in p else "#c97b26" if "PD" in p else "#607d8b" for p in df["pair"]]

    axes[0].bar(range(len(df)), df["conjfdr_all_n"], color=colors)
    axes[0].set_xticks(range(len(df)))
    axes[0].set_xticklabels(df["pair"], rotation=45, ha="right")
    axes[0].set_ylabel("Significant SNPs")
    axes[0].set_title("conjFDR<0.05 SNP count")

    axes[1].bar(range(len(df)), df["conjfdr_loci_n"], color=colors)
    axes[1].set_xticks(range(len(df)))
    axes[1].set_xticklabels(df["pair"], rotation=45, ha="right")
    axes[1].set_ylabel("Independent loci")
    axes[1].set_title("conjFDR<0.05 loci count")

    fig.tight_layout()
    fig.savefig(SUPP_PLEIO / "Supplementary_pleioFDR_counts.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_coloc_figures() -> None:
    plt.rcParams.update({"font.family": "DejaVu Serif", "figure.facecolor": "white"})
    df = pd.read_csv(BASE / "07_coloc_factor_ndd" / "coloc_factor_ndd_summary.tsv", sep="\t")
    df["pair_label"] = df["pair"].str.replace("_", "-", regex=False)
    df["coord_label"] = (
        "chr"
        + df["chrnum"].astype(str)
        + ":"
        + (df["region_start"] / 1e6).round(1).map(lambda x: f"{x:.1f}")
        + "-"
        + (df["region_end"] / 1e6).round(1).map(lambda x: f"{x:.1f}")
        + " Mb"
    )
    df["task_label"] = df["pair_label"] + " | " + df["coord_label"]
    df["priority_group"] = df["priority"].map(
        {"PWCoCo_priority": "Complex region", "coloc_first": "Standard region"}
    ).fillna("Standard region")
    palette = {"AD": "#0b6e4f", "PD": "#c97b26", "LBD": "#607d8b"}

    top = (
        df.sort_values(["PP.H4.abf", "PP.H3.abf"], ascending=[False, False])
        .head(8)
        .sort_values("PP.H4.abf", ascending=True)
        .copy()
    )

    fig = plt.figure(figsize=(15.5, 9), facecolor="white")
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.55], wspace=0.34)

    ax1 = fig.add_subplot(gs[0, 0])
    y = range(len(top))
    ax1.hlines(
        y,
        xmin=0,
        xmax=top["PP.H4.abf"],
        color=[palette.get(x, "#666666") for x in top["disease"]],
        linewidth=3,
        alpha=0.85,
    )
    ax1.scatter(
        top["PP.H4.abf"],
        y,
        s=90,
        color=[palette.get(x, "#666666") for x in top["disease"]],
        edgecolor="white",
        linewidth=0.9,
        zorder=3,
    )
    ax1.axvline(0.5, color="#9e9e9e", linestyle="--", linewidth=1.0)
    ax1.axvline(0.8, color="#444444", linestyle=":", linewidth=1.2)
    ax1.set_yticks(list(y))
    ax1.set_yticklabels(top["task_label"], fontsize=9)
    ax1.set_xlim(0, 1.02)
    ax1.set_xlabel("PP.H4.abf")
    ax1.set_title("Colocalization support across top regions", fontsize=12)
    ax1.grid(axis="x", color="#e5e5e5", linewidth=0.8)

    ax2 = fig.add_subplot(gs[0, 1])
    stack = top[["PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf"]].copy()
    stack.index = top["pair_label"] + "\n" + top["coord_label"]
    colors = ["#d9d9d9", "#bdbdbd", "#9ecae1", "#f4a261", "#2a9d8f"]
    left = [0.0] * len(stack)
    for idx, col in enumerate(stack.columns):
        ax2.barh(
            stack.index,
            stack[col],
            left=left,
            color=colors[idx],
            edgecolor="white",
            linewidth=0.6,
            label=col.replace(".abf", ""),
        )
        left = [l + v for l, v in zip(left, stack[col])]
    ax2.set_xlim(0, 1)
    ax2.set_xlabel("Posterior probability")
    ax2.set_title("Posterior decomposition in selected regions", fontsize=12)
    ax2.legend(frameon=False, loc="lower right", ncol=1, fontsize=9)
    ax2.tick_params(axis="y", labelsize=9, pad=2)

    fig.text(0.02, 0.98, "A", fontsize=16, fontweight="bold", va="top")
    fig.text(0.51, 0.98, "B", fontsize=16, fontweight="bold", va="top")
    fig.subplots_adjust(left=0.14, right=0.98)
    fig.savefig(MAIN_DIR / "Figure3_coloc_main.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    heat = df.copy()
    heat = heat.sort_values(["disease", "trait", "chrnum", "region_start"]).reset_index(drop=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, max(6, 0.38 * len(heat))), sharey=True)
    matrices = [("PP.H4.abf", "PP.H4"), ("PP.H3.abf", "PP.H3")]
    cmap_colors = {"PP.H4.abf": "#2a9d8f", "PP.H3.abf": "#f4a261"}

    for ax, (col, title) in zip(axes, matrices):
        ax.barh(
            heat["task_label"],
            heat[col],
            color=[cmap_colors[col]] * len(heat),
            edgecolor="white",
            linewidth=0.4,
        )
        ax.axvline(0.5, color="#9e9e9e", linestyle="--", linewidth=1.0)
        ax.axvline(0.8, color="#444444", linestyle=":", linewidth=1.0)
        ax.set_xlim(0, 1.0)
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("Posterior probability")
        ax.grid(axis="x", color="#efefef", linewidth=0.8)

    axes[0].set_ylabel("Pair and region")
    fig.suptitle("Region-level coloc summary across all factor-NDD candidate loci", fontsize=14)
    fig.tight_layout()
    fig.savefig(SUPP_COLOC / "Supplementary_coloc_overview.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    strong = df.loc[
        (df["PP.H4.abf"] >= 0.5) | (df["PP.H3.abf"] >= 0.5),
        [
            "pair_label",
            "region_id",
            "chrnum",
            "region_start",
            "region_end",
            "priority_group",
            "PP.H3.abf",
            "PP.H4.abf",
            "coloc_interpretation",
            "top_snp_h4",
        ],
    ].copy()
    strong.to_csv(SUPP_COLOC / "coloc_high_priority_regions.tsv", sep="\t", index=False)


def make_pwcoco_figures() -> None:
    df = pd.read_csv(BASE / "08_pwcoco" / "pwcoco_region_best_summary.tsv", sep="\t")
    df["pair_label"] = df["pair"].str.replace("_", "-", regex=False)
    df["task_label"] = df["pair_label"] + " | " + df["region_id"]
    df = df.sort_values(["region_id", "pair"]).reset_index(drop=True)

    fig, axes = plt.subplots(1, 2, figsize=(13, max(4.5, 0.75 * len(df))), sharey=True)
    y = range(len(df))
    colors = ["#0b6e4f" if "AD" in p else "#c97b26" for p in df["pair_label"]]

    axes[0].barh(y, df["H4_uncond"], color=colors, edgecolor="white", linewidth=0.6)
    axes[0].axvline(0.5, color="#9e9e9e", linestyle="--", linewidth=1.0)
    axes[0].axvline(0.8, color="#444444", linestyle=":", linewidth=1.0)
    axes[0].set_yticks(list(y))
    axes[0].set_yticklabels(df["task_label"])
    axes[0].set_xlim(0, 1)
    axes[0].set_xlabel("Unconditioned H4")
    axes[0].set_title("PWCoCo initial colocalization")
    axes[0].grid(axis="x", color="#efefef", linewidth=0.8)

    axes[1].barh(y, df["H4_best_cond"], color="#607d8b", edgecolor="white", linewidth=0.6, label="Best conditioned H4")
    axes[1].barh(y, df["H3_best_cond"], color="#f4a261", edgecolor="white", linewidth=0.6, alpha=0.75, label="Best conditioned H3")
    axes[1].axvline(0.5, color="#9e9e9e", linestyle="--", linewidth=1.0)
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("Conditioned posterior")
    axes[1].set_title("PWCoCo best conditioned model")
    axes[1].legend(frameon=False, loc="lower right")
    axes[1].grid(axis="x", color="#efefef", linewidth=0.8)

    fig.suptitle("PWCoCo summary across tested candidate regions", fontsize=14)
    fig.tight_layout()
    fig.savefig(SUPP_PWCOCO / "Supplementary_PWCoCo_summary.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ensure_dirs()
    make_ldsc_mixer_main()
    make_ldsc_supp()
    make_mixer_supp()
    make_pleio_main()
    archive_pleio_figures()
    make_pleio_supp_summary()
    make_coloc_figures()
    make_pwcoco_figures()
    print(f"Wrote figures to: {FIG_BASE}")


if __name__ == "__main__":
    main()
