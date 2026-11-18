"""Plot the real NRSA survey results produced by analyze.R."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path("figures/results")
INK = "#1e3b4b"
TEAL = "#077a73"
RUST = "#ba6f42"
NUTRIENTS = {"nitrate": ("Nitrate", "mg N/L"),
             "phosphorus": ("Total phosphorus", "mg/L"),
             "nitrogen": ("Total nitrogen", "mg N/L")}


def read_table(name):
    with (Path("runs") / f"{name}.csv").open(encoding="utf-8", newline="") as source:
        return list(csv.DictReader(source))


def style(fig, axes, title, note, *, bottom=0.19):
    fig.set_facecolor("#fcfdfb")
    fig.subplots_adjust(left=0.09, right=0.98, bottom=bottom, top=0.76, wspace=0.33)
    fig.text(0.046, 0.92, title, fontsize=18, weight="bold", color=INK)
    fig.text(0.046, 0.040, note, fontsize=9, color="#4e6874")
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", alpha=0.20)
        ax.set_axisbelow(True)


def region_comparison(estimates):
    fig, axes = plt.subplots(1, 3, figsize=(12, 5.4))
    style(fig, axes, "Survey weights change mean nutrient concentrations",
          "EPA NRSA 2018–2019 · matched observations with half-limit substitution · weighted 95% stratified bootstrap intervals")
    regions = sorted({r["zone"] for r in estimates} - {"National"})
    for ax, (code, (title, unit)) in zip(axes, NUTRIENTS.items()):
        rows = {r["zone"]: r for r in estimates if r["analyte"] == code}
        rows = [rows[region] for region in regions]
        pos = np.arange(len(rows))
        value = np.array([float(r["weighted"]) for r in rows])
        lo = np.array([float(r["ci_low"]) for r in rows])
        hi = np.array([float(r["ci_high"]) for r in rows])
        ax.errorbar(value, pos, xerr=[value - lo, hi - value], fmt="o", ms=5, color=TEAL,
                    capsize=2, label="Weighted")
        ax.scatter([float(r["unweighted"]) for r in rows], pos + 0.15, color=RUST, s=22, marker="D",
                   label="Unweighted")
        ax.set(yticks=pos, yticklabels=regions, xlabel=unit, title=title)
        ax.invert_yaxis()
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=2, frameon=False, loc="upper right", bbox_to_anchor=(0.97, 0.87))
    fig.savefig(OUT / "ecoregion-weighting.png", dpi=175)
    plt.close(fig)


def censoring_comparison(estimates):
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    style(fig, axes, "Censoring sensitivity and retained survey weight",
          "EPA NRSA 2018–2019 · each point represents one nutrient and one ecoregion · valid values keep the same weight")
    for code, (label, _), marker in zip(NUTRIENTS, NUTRIENTS.values(), ("o", "s", "^")):
        rows = [r for r in estimates if r["analyte"] == code and r["zone"] != "National"]
        censored = [100 * float(r["censored"]) / float(r["stations"]) for r in rows]
        span = [100 * (float(r["upper_sub_mean"]) - float(r["lower_sub_mean"])) / float(r["weighted"]) for r in rows]
        axes[0].scatter(censored, span, label=label, marker=marker, s=45, alpha=0.8)
        axes[1].scatter([float(r["kish_n"]) for r in rows],
                        [100 * float(r["weight_coverage"]) for r in rows], label=label, marker=marker, s=45, alpha=0.8)
    axes[0].set(xlabel="Censored observations (%)", ylabel="Zero-to-limit spread (% of estimated mean)")
    axes[1].set(xlabel="Kish effective sample size", ylabel="Retained population weight (%)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, ncol=3, frameon=False, loc="upper right", bbox_to_anchor=(0.97, 0.86))
    fig.savefig(OUT / "censoring-and-coverage.png", dpi=175)
    plt.close(fig)


def pca_summary():
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    style(fig, axes, "Chemical variation in the weighted complete-case sample",
          "EPA NRSA 2018–2019 · weighted PCA of log(1 + concentration) · loading signs depend on the axis convention", bottom=0.22)
    explained = read_table("pca-explained")
    labels = [f"PC{r['component']}" for r in explained]
    perc = [100 * float(r["explained"]) for r in explained]
    bars = axes[0].bar(labels, perc, color=["#067f6e", "#69a68e", "#8ca9a3", "#afc0b9"])
    axes[0].bar_label(bars, fmt="%.1f%%", padding=3, color=INK)
    axes[0].set(ylabel="Variance explained (%)", ylim=(0, max(perc) * 1.30))
    loadings = read_table("pca-loadings")
    pos = np.arange(len(loadings))
    for i in range(2):
        axes[1].barh(pos + (i - .5) * 0.30, [float(r[f"X{i + 1}"]) for r in loadings],
                     height=.27, label=f"PC{i + 1}", color=(TEAL, RUST)[i])
    axes[1].axvline(0, color=INK, linewidth=.8)
    axes[1].set(yticks=pos, yticklabels=["Nitrate", "Total P", "Total N", "Dissolved C"], xlabel="Eigenvector coefficient")
    axes[1].invert_yaxis()
    axes[1].legend(frameon=False, ncol=2, loc="lower right", bbox_to_anchor=(1, 1.01))
    fig.savefig(OUT / "weighted-nutrient-pca.png", dpi=175)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    data = read_table("nutrient-estimates")
    region_comparison(data)
    censoring_comparison(data)
    pca_summary()
