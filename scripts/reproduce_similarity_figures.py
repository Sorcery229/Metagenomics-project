#!/usr/bin/env python3
"""Reproduce the cross-database per-base similarity distribution figures.

Reads the committed per-domain similarity tables under data/similarity/ and
regenerates, deterministically:

  * figures/fig_similarity_strain.png    (Bacteria, strain-level)
  * figures/fig_similarity_species.png   (Viruses + Fungi, species-level)
  * figures/similarity_band_distribution.csv   (the numbers behind both figures)

Similarity is recomputed from the primitive columns and cross-checked:

    similarity == matched_bases / subject_length

  matched_bases = n_of_1 + n_of_5 for bacteria (per-base identical + N-corrected
  positions from the BLAST alignment); = aligned matches for viruses and fungi.

Band definitions (identical across all three domains):

    100%      : similarity >= 0.999999
    95-<100%  : 0.95 <= similarity < 0.999999
    <95%      : similarity < 0.95

Run from the repository root:

    python scripts/reproduce_similarity_figures.py
"""
import os
import pandas as pd
import matplotlib.pyplot as plt

THR_100 = 0.999999
THR_95 = 0.95
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data", "similarity")
FIGS = os.path.join(HERE, "figures")
os.makedirs(FIGS, exist_ok=True)

BANDS = ["100%", "95-<100%", "<95%"]
CBAND = {"100%": "#2166ac", "95-<100%": "#67a9cf", "<95%": "#ef8a62"}
META_GREY = "#8a8f98"


def load(name):
    df = pd.read_csv(os.path.join(DATA, name))
    recomputed = df["matched_bases"] / df["subject_length"]
    assert (recomputed - df["similarity"]).abs().max() < 1e-9, \
        f"{name}: similarity != matched_bases / subject_length"
    return df["similarity"].astype(float)


def band_counts(s):
    n = len(s)
    c100 = int((s >= THR_100).sum())
    c95 = int(((s >= THR_95) & (s < THR_100)).sum())
    clt = int((s < THR_95).sum())
    return n, c100, c95, clt, float(s.mean()), float(s.std())


def plot_panel(ax, counts, n, mean, std, title, sublabel, letter=None):
    pct = [100 * c / n for c in counts]
    xs = range(3)
    ax.bar(xs, pct, color=[CBAND[b] for b in BANDS], width=0.68,
           edgecolor="white", linewidth=0.8)
    ax.set_xticks(list(xs))
    ax.set_xticklabels(BANDS)
    ax.set_ylim(0, 116)
    ax.set_ylabel("genome pairs (%)")
    ax.set_title(title, loc="center")
    for i, (p, c) in enumerate(zip(pct, counts)):
        if p >= 15:
            ax.text(i, p + 2.5, f"{p:.1f}%", ha="center", va="bottom",
                    fontsize=9, fontweight="bold")
            ax.text(i, p - 5, f"n={c:,}", ha="center", va="center",
                    fontsize=7.5, color="white")
        else:
            ax.text(i, p + 1.0, f"n={c:,}", ha="center", va="bottom",
                    fontsize=7, color="#555555")
            ax.text(i, p + 8.5, f"{p:.1f}%", ha="center", va="bottom",
                    fontsize=9, fontweight="bold")
    ax.text(0.5, -0.30, sublabel, transform=ax.transAxes, ha="center",
            va="top", fontsize=8.5, color=META_GREY)
    ax.text(0.98, 0.99, f"N = {n:,}\nmean {mean:.3f} +/- {std:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
            color="#333333")
    for sp in ["top", "right"]:
        ax.spines[sp].set_visible(False)
    if letter:
        ax.text(-0.18, 1.08, letter, transform=ax.transAxes, fontsize=13,
                fontweight="bold", va="top", ha="right")


def main():
    plt.rcParams.update({"font.family": "sans-serif", "font.size": 9,
                         "axes.titlesize": 10, "svg.fonttype": "none"})
    domains = {
        "Bacteria": (load("bacteria_similarity.csv"), "strain"),
        "Viruses": (load("viruses_similarity.csv"), "species"),
        "Fungi": (load("fungi_similarity.csv"), "species"),
    }
    rows, D = [], {}
    for dom, (s, lvl) in domains.items():
        n, c100, c95, clt, mu, sd = band_counts(s)
        D[dom] = ([c100, c95, clt], n, mu, sd)
        rows.append(dict(domain=dom, match_level=lvl, n_pairs=n,
                         n_100=c100, n_95to100=c95, n_below95=clt,
                         pct_100=round(100 * c100 / n, 1),
                         pct_95to100=round(100 * c95 / n, 1),
                         pct_below95=round(100 * clt / n, 1),
                         mean_similarity=round(mu, 4),
                         std_similarity=round(sd, 4)))
    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(FIGS, "similarity_band_distribution.csv"),
                 index=False)
    print(table.to_string(index=False))

    # species-level figure (Viruses, Fungi)
    figA, axes = plt.subplots(1, 2, figsize=(7.0, 3.4))
    for ax, dom, lt in [(axes[0], "Viruses", "a"), (axes[1], "Fungi", "b")]:
        cnts, n, mu, sd = D[dom]
        plot_panel(ax, cnts, n, mu, sd, dom, "species-level overlap", letter=lt)
    figA.suptitle("Cross-database per-base similarity - species-level comparisons",
                  y=1.02, fontsize=10)
    figA.tight_layout()
    figA.savefig(os.path.join(FIGS, "fig_similarity_species.png"),
                 dpi=300, bbox_inches="tight")

    # strain-level figure (Bacteria)
    figB, ax = plt.subplots(1, 1, figsize=(3.8, 3.6))
    cnts, n, mu, sd = D["Bacteria"]
    plot_panel(ax, cnts, n, mu, sd, "Bacteria", "strain-level overlap", letter="a")
    figB.suptitle("Cross-database per-base similarity - strain-level comparison",
                  y=1.02, fontsize=10)
    figB.tight_layout()
    figB.savefig(os.path.join(FIGS, "fig_similarity_strain.png"),
                 dpi=300, bbox_inches="tight")
    print("\nWrote figures to", FIGS)


if __name__ == "__main__":
    main()
