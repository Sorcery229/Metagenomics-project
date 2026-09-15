#!/usr/bin/env python3
"""Reproduce the strain- and species-level overlap Venn diagrams from database metadata.

Every region count is computed as a set operation on committed metadata files; nothing
is hard-coded. Two figures are produced:

  figures/venn_strain_overlap.png    Bacteria (strain level): RefSeq vs BV-BRC
  figures/venn_species_overlap.png   Viruses + Fungi (species level)

Set definitions (key -> source file, all under data/):

  Bacteria (strain)  representative (taxid, strain) keys
     RefSeq   data/venn/bacteria_refseq_strain_keys.csv.gz
              (unique taxid + strain of RefSeq representative genomes, from
               NCBI assembly_summary.txt; strain from infraspecific_name)
     BV-BRC   data/venn/bacteria_bvbrc_strain_keys.csv.gz
              (unique taxon_id + strain from the BV-BRC genome_metadata table)

  Viruses (species)  NCBI virus taxonomy id
     RefSeq        data/assembly_summaries/assembly_summary_viral.txt  column 'taxid'
     VirusHostDB   data/venn/virushostdb.tsv.gz                        column 'virus tax id'
                   (full VirusHostDB list; source https://www.genome.jp/ftp/db/virushostdb/virushostdb.tsv)

  Fungi (species)    NCBI (species) taxonomy id
     Ensembl   data/venn/species_EnsemblFungi.txt                            column 'taxonomy_id'
     RefSeq    data/assembly_summaries/assembly_summary_fungi_refseq.txt     column 'taxid'
     FungiDB   data/venn/fungidb_organisms.csv                              column 'species_ncbi_tax_id'
               (VEuPathDB/FungiDB AllOrganisms report; source https://fungidb.org)

Counts drift slightly from the manuscript screenshot because VirusHostDB and FungiDB
are living databases that have grown since the paper; the key set intersections
(RefSeq/BV-BRC strain, Ensembl/RefSeq fungal) reproduce the manuscript exactly.
"""
import os, re, gzip
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib_venn import venn2, venn2_circles, venn3, venn3_circles

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = lambda *p: os.path.join(ROOT, *p)
FIGS = D("figures"); os.makedirs(FIGS, exist_ok=True)


def strain_keys(path):
    df = pd.read_csv(path, dtype=str)
    return set(zip(df["taxid"].astype(str).str.strip(), df["strain"].astype(str).str.strip()))


def parse_ensembl(path):
    # header row is one column short of the data namespace on load; parse by hand
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    hdr = lines[0].lstrip("#").split("\t")
    recs = [dict(zip(hdr, l.split("\t"))) for l in lines[1:] if l.strip()]
    df = pd.DataFrame(recs)
    return set(df["taxonomy_id"].dropna().astype(str).str.strip()) - {""}


def col_set(path, col, sep="\t", skiprows=0, gz=False):
    op = gzip.open(path, "rt") if gz else open(path, encoding="utf-8", errors="replace")
    df = pd.read_csv(op, sep=sep, dtype=str, skiprows=skiprows)
    return set(df[col].dropna().astype(str).str.strip()) - {""}


def compute_regions():
    # Bacteria (strain)
    rs = strain_keys(D("data/venn/bacteria_refseq_strain_keys.csv.gz"))
    bv = strain_keys(D("data/venn/bacteria_bvbrc_strain_keys.csv.gz"))
    BAC = dict(refseq_only=len(rs - bv), shared=len(rs & bv), bvbrc_only=len(bv - rs))

    # Viruses (species)
    rv = col_set(D("data/assembly_summaries/assembly_summary_viral.txt"), "taxid", skiprows=1)
    vh = col_set(D("data/venn/virushostdb.tsv.gz"), "virus tax id", gz=True)
    VIR = dict(vhdb_only=len(vh - rv), shared=len(vh & rv), refseq_only=len(rv - vh))

    # Fungi (species)
    E = parse_ensembl(D("data/venn/species_EnsemblFungi.txt"))
    R = col_set(D("data/assembly_summaries/assembly_summary_fungi_refseq.txt"), "taxid", skiprows=1)
    F = col_set(D("data/venn/fungidb_organisms.csv"), "species_ncbi_tax_id", sep=",")
    tri = len(E & R & F)
    FUN = dict(E_only=len(E - R - F), R_only=len(R - E - F), F_only=len(F - E - R),
               ER=len(E & R) - tri, EF=len(E & F) - tri, RF=len(R & F) - tri, triple=tri)
    return BAC, VIR, FUN


def setnum(v, mp, fs):
    for lab, txt in mp.items():
        o = v.get_label_by_id(lab)
        if o:
            o.set_text(f"{txt}"); o.set_fontweight("bold"); o.set_fontsize(fs)


def fig_strain(BAC):
    fig, ax = plt.subplots(1, 1, figsize=(4.6, 4.4))
    s = (BAC["refseq_only"], BAC["bvbrc_only"], BAC["shared"])
    v = venn2(subsets=s, set_labels=("", ""), ax=ax, set_colors=("#7B9FD4", "#8FCB9B"), alpha=0.55)
    venn2_circles(subsets=s, ax=ax, lw=0.8, color="grey")
    setnum(v, {"10": s[0], "11": s[2], "01": s[1]}, 13)
    v.get_label_by_id("10").set_position((-0.62, 0.0))
    v.get_label_by_id("11").set_position((-0.28, 0.0))
    ax.text(-0.02, 1.04, "a", transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")
    ax.text(0.30, 1.02, "Bacteria", transform=ax.transAxes, fontsize=15)
    ax.text(-0.55, -0.60, "RefSeq", fontsize=11); ax.text(0.34, -0.85, "BV-BRC", fontsize=11)
    ax.text(0.5, -0.26, "Strain overlap", transform=ax.transAxes, ha="center", fontsize=13)
    fig.subplots_adjust(bottom=0.14, top=0.86)
    fig.savefig(D("figures/venn_strain_overlap.png"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def fig_species(VIR, FUN):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.8))
    ax = axes[0]
    sv = (VIR["vhdb_only"], VIR["refseq_only"], VIR["shared"])
    v = venn2(subsets=sv, set_labels=("", ""), ax=ax, set_colors=("#D98C8C", "#B79FD4"), alpha=0.55)
    venn2_circles(subsets=sv, ax=ax, lw=0.8, color="grey")
    setnum(v, {"10": sv[0], "11": sv[2], "01": sv[1]}, 12)
    ax.text(-0.02, 1.05, "a", transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")
    ax.text(0.5, 1.04, "Viruses", transform=ax.transAxes, ha="center", fontsize=15)
    ax.text(0.14, 0.12, "VirusHostDB", transform=ax.transAxes, fontsize=11)
    ax.text(0.72, 0.12, "RefSeq", transform=ax.transAxes, fontsize=11)
    ax = axes[1]
    s3 = (FUN["E_only"], FUN["R_only"], FUN["ER"], FUN["F_only"], FUN["EF"], FUN["RF"], FUN["triple"])
    v3 = venn3(subsets=s3, set_labels=("", "", ""), ax=ax, set_colors=("#E8D96A", "#7B9FD4", "#B79FD4"), alpha=0.55)
    venn3_circles(subsets=s3, ax=ax, lw=0.8, color="grey")
    setnum(v3, {"100": s3[0], "010": s3[1], "110": s3[2], "001": s3[3],
                "101": s3[4], "011": s3[5], "111": s3[6]}, 11)
    ax.text(-0.02, 1.05, "b", transform=ax.transAxes, fontsize=18, fontweight="bold", va="top")
    ax.text(0.5, 1.04, "Fungi", transform=ax.transAxes, ha="center", fontsize=15)
    ax.text(0.10, 0.90, "Ensembl", transform=ax.transAxes, fontsize=11)
    ax.text(0.66, 0.90, "RefSeq", transform=ax.transAxes, fontsize=11)
    ax.text(0.42, 0.02, "FungiDB", transform=ax.transAxes, fontsize=11)
    fig.text(0.5, 0.02, "Species overlap", ha="center", fontsize=13)
    fig.subplots_adjust(bottom=0.12, top=0.88, wspace=0.15)
    fig.savefig(D("figures/venn_species_overlap.png"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    BAC, VIR, FUN = compute_regions()
    rows = [
        ("Bacteria", "strain", "RefSeq only", BAC["refseq_only"], 58793),
        ("Bacteria", "strain", "RefSeq ∩ BV-BRC", BAC["shared"], 295055),
        ("Bacteria", "strain", "BV-BRC only", BAC["bvbrc_only"], 617380),
        ("Viruses", "species", "VirusHostDB only", VIR["vhdb_only"], 20688),
        ("Viruses", "species", "VirusHostDB ∩ RefSeq", VIR["shared"], 14294),
        ("Viruses", "species", "RefSeq only", VIR["refseq_only"], 207),
        ("Fungi", "species", "Ensembl only", FUN["E_only"], 785),
        ("Fungi", "species", "RefSeq only", FUN["R_only"], 290),
        ("Fungi", "species", "FungiDB only", FUN["F_only"], 156),
        ("Fungi", "species", "Ensembl ∩ RefSeq", FUN["ER"], 318),
        ("Fungi", "species", "Ensembl ∩ FungiDB", FUN["EF"], 47),
        ("Fungi", "species", "RefSeq ∩ FungiDB", FUN["RF"], 16),
        ("Fungi", "species", "Ensembl ∩ RefSeq ∩ FungiDB", FUN["triple"], 35),
    ]
    ct = pd.DataFrame(rows, columns=["domain", "level", "region",
                                     "reproduced_from_metadata", "screenshot_manuscript"])
    ct.to_csv(D("figures/venn_overlap_counts.csv"), index=False)
    print(ct.to_string(index=False))
    fig_strain(BAC)
    fig_species(VIR, FUN)
    print("\nWrote figures to", FIGS)


if __name__ == "__main__":
    main()
