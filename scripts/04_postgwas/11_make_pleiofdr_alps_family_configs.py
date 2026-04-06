from __future__ import annotations

from pathlib import Path


BASE_OUT = Path(r"D:\文章\GS\postgwas\03b_pleiofdr_alps_family_vs_ndd")
MAT_DIR = BASE_OUT / "inputs_mat"
CONFIG_DIR = BASE_OUT / "configs"
OUTPUT_DIR = BASE_OUT / "results"
PFDIR = Path(r"D:\pleioFDR\pleiofdr-master")
REFMAT = PFDIR / "ref9545380_1kgPhase3eur_LDr2p1.mat"

ALPS_TRAITS = ["aALPS", "Left_ALPS", "mALPS", "pALPS", "Right_ALPS", "Mean_ALPS", "tALPS"]
NDD_TRAITS = ["AD", "PD", "LBD"]


def config_text(trait1: str, trait2: str) -> str:
    return f"""# Auto-generated config for GS ALPS-family pleioFDR analysis
reffile={REFMAT}
traitfolder=
traitfile1={MAT_DIR / f"{trait1}_fdr.mat"}
traitname1={trait1}
traitfiles={{'{MAT_DIR / f"{trait2}_fdr.mat"}'}}
traitnames={{'{trait2}'}}
outputdir={OUTPUT_DIR / f"{trait1}_{trait2}"}
stattype=conjfdr
fdrthresh=0.05
randprune=true
randprune_n=20
exclude_chr_pos=[6 25119106 33854733]
manh_fontsize_genenames=12
manh_yspace=0.75
manh_ymargin=0.25
manh_colorlist=[1 0 0; 1 0.5 0 ; 0 0.75 0.75; 0 0.5 0; 0.75 0 0.75; 0 0 1; 0 1 0; 0 1 1]
refinfo={PFDIR / '9545380.ref'}
reset_pruneidx=true
randprune_repeats=default
pthresh=1
perform_gc=true
use_standard_gc=false
randprune_gc=true
exclude_from_discovery=false
mafthresh = 0.005
exclude_ambiguous_snps = true
onscreen = false
"""


def main() -> None:
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for trait1 in ALPS_TRAITS:
        for trait2 in NDD_TRAITS:
            cfg_path = CONFIG_DIR / f"{trait1}_{trait2}.config.txt"
            cfg_path.write_text(config_text(trait1, trait2), encoding="utf-8")
    print(f"Wrote configs to {CONFIG_DIR}")


if __name__ == "__main__":
    main()
