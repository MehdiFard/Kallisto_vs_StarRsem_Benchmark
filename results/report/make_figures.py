#!/usr/bin/env python3
"""Figures comparing STAR+RSEM and kallisto (gene and isoform level, within and across samples).

Run from the project root after the pipeline has finished:   python results/report/make_figures.py
Needs python >= 3.10 with pandas, numpy, scipy and matplotlib. Reads results/summary, results/expression,
results/phenotypes and results/eqtl; writes results/report/figures/Fig*.png and results/report/figures/key_numbers.json.
"""
import gzip
import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

OUT = "results/report/figures"
os.makedirs(OUT, exist_ok=True)

# ----------------------------------------------------------------------------------------- style
RSEM, KAL = "#2a78d6", "#eb6834"               # tool identity (categorical slots 1, 2)
GENE, ISO = "#1baf7a", "#4a3aa7"               # level identity when both share a panel (slots 3, 7)
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
SHARED = "#a8a69e"                             # "both tools" in overlap bars (neutral)
DENSITY = LinearSegmentedColormap.from_list("density", ["#e6e3f6", "#9085e9", "#4a3aa7", "#211763"])
TOOL = {"rsem": ("STAR+RSEM", RSEM), "kallisto": ("kallisto", KAL)}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
    "legend.fontsize": 6.5, "legend.frameon": False,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.6, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "xtick.color": INK2, "ytick.color": INK2, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "text.color": INK,
    "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight", "pdf.fonttype": 42, "svg.fonttype": "none",
})


def panel(ax, letter, title=None):
    ax.text(-0.16, 1.08, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="left")
    if title:
        ax.set_title(title, pad=6)


def save(fig, name):
    fig.savefig(f"{OUT}/{name}.png", facecolor="white")
    plt.close(fig)
    print("wrote", name)


def hexdens(ax, x, y, gridsize=90, extent=None):
    hb = ax.hexbin(x, y, gridsize=gridsize, bins="log", mincnt=1, cmap=DENSITY, linewidths=0, extent=extent,
                   rasterized=True)
    return hb


def identity(ax, lo, hi):
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.7, zorder=0)


def binned_trend(ax, x, y, bins, color=INK, label=None):
    """Median and IQR of y within x bins (needs >= 30 points per bin)."""
    idx = np.digitize(x, bins)
    cx, med, q1, q3 = [], [], [], []
    for i in range(1, len(bins)):
        v = y[idx == i]
        if len(v) >= 30:
            cx.append((bins[i - 1] + bins[i]) / 2)
            med.append(np.median(v)); q1.append(np.quantile(v, 0.25)); q3.append(np.quantile(v, 0.75))
    ax.fill_between(cx, q1, q3, color=color, alpha=0.15, lw=0)
    ax.plot(cx, med, color=color, lw=1.6, label=label)


def strip_box(ax, data, positions, colors, width=0.5, s=3, alpha=0.35, seed=1):
    rng = np.random.default_rng(seed)
    for d, p, c in zip(data, positions, colors):
        d = np.asarray(d)[~np.isnan(d)]
        ax.scatter(p + rng.uniform(-width / 2.6, width / 2.6, len(d)), d, s=s, color=c, alpha=alpha, lw=0,
                   rasterized=True, zorder=2)
        bp = ax.boxplot(d, positions=[p], widths=width, showfliers=False, patch_artist=True, zorder=3,
                        medianprops=dict(color=INK, lw=1.2), whiskerprops=dict(color=INK2, lw=0.7),
                        capprops=dict(color=INK2, lw=0.7), boxprops=dict(edgecolor=INK2, lw=0.7))
        bp["boxes"][0].set_facecolor("none")


def stacked(ax, i, parts, tot):
    """Horizontal 100% bar; counts inside wide segments, above narrow ones."""
    left = 0
    for n, c in parts:
        w = n / tot * 100
        ax.barh(i, w, left=left, color=c, height=0.5, edgecolor="white", lw=1.2)
        if w >= 20:
            ax.text(left + w / 2, i, f"{n:,}", ha="center", va="center", fontsize=6.3,
                    color=INK if c == SHARED else "white")
        else:
            ax.text(left + w / 2, i - 0.3, f"{n:,}", ha="center", va="bottom", fontsize=6.3, color=c)
        left += w


def colorbar(fig, hb, ax):
    cb = fig.colorbar(hb, ax=ax, fraction=0.05, pad=0.03)
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.5, length=1.5)
    return cb


def sp(x, y):
    m = ~(np.isnan(x) | np.isnan(y))
    return stats.spearmanr(x[m], y[m])[0]


# ----------------------------------------------------------------------------------------- data
S = "results/summary"
qc = pd.read_csv(f"{S}/mapping_qc.tsv", sep="\t", index_col="sample_id")
comp = pd.read_csv(f"{S}/compute_benchmark.tsv", sep="\t", index_col="sample_id")
ps = {lv: pd.read_csv(f"{S}/quant_correlation_{lv}_per_sample.tsv", sep="\t", index_col="sample_id")
      for lv in ("gene", "isoform")}
pf = {lv: pd.read_csv(f"{S}/quant_correlation_{lv}_per_feature.tsv.gz", sep="\t", index_col=0)
      for lv in ("gene", "isoform")}
tested = {(t, lv): set(pd.read_csv(f"results/phenotypes/{t}/{lv}.groups.tsv", sep="\t", header=None)[0])
          for t in TOOL for lv in ("gene", "isoform")}
cis = {(t, lv): pd.read_csv(f"results/eqtl/{t}/{lv}/cis.tsv.gz", sep="\t", index_col=0)
       for t in TOOL for lv in ("gene", "isoform")}
FDR = 0.05
K = {}  # key numbers, written to key_numbers.json

# ========================================================================= Fig 1: mapping QC
fig, axs = plt.subplots(1, 4, figsize=(7.4, 2.1), gridspec_kw=dict(wspace=0.62, width_ratios=[1, 1.3, 1, 1]))
ax = axs[0]
ax.hist(qc["star_input_reads"] / 1e6, bins=40, color=INK2, lw=0)
ax.set_xlabel("Read pairs per sample (M)"); ax.set_ylabel("Samples")
panel(ax, "a", "Depth")
ax = axs[1]
cols = [("star_unique_pct", "STAR\nunique", RSEM), ("star_multi_pct", "STAR\nmulti", RSEM),
        ("kallisto_pseudoaligned_pct", "kal.\naligned", KAL), ("kallisto_unique_pct", "kal.\n1 tx", KAL)]
strip_box(ax, [qc[c].values for c, _, _ in cols], range(4), [c for _, _, c in cols], s=2.5)
ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l, _ in cols]); ax.set_ylabel("% of reads"); ax.set_ylim(0, 100)
panel(ax, "b", "Mapping rates")
ax = axs[2]
ax.scatter(qc["star_unique_pct"], qc["kallisto_pseudoaligned_pct"], s=5, color=INK2, alpha=0.6, lw=0)
identity(ax, 40, 100); ax.set_xlim(45, 97); ax.set_ylim(40, 97)
ax.set_xlabel("STAR uniquely mapped (%)"); ax.set_ylabel("kallisto pseudoaligned (%)")
ax.text(0.04, 0.96, f"ρ = {sp(qc['star_unique_pct'].values, qc['kallisto_pseudoaligned_pct'].values):.2f}",
        transform=ax.transAxes, va="top", color=INK2)
for sid in qc["kallisto_pseudoaligned_pct"].nsmallest(2).index:
    xs = qc.loc[sid, "star_unique_pct"]
    ax.annotate(sid, (xs, qc.loc[sid, "kallisto_pseudoaligned_pct"]), fontsize=5.8, color=INK2,
                xytext=(-4 if xs > 70 else 4, -2), textcoords="offset points", ha="right" if xs > 70 else "left")
panel(ax, "c", "STAR vs kallisto")
ax = axs[3]
r, k = ps["gene"]["rsem_total_counts"] / 1e6, ps["gene"]["kallisto_total_counts"] / 1e6
ax.scatter(r, k, s=5, color=INK2, alpha=0.6, lw=0)
hi = max(r.max(), k.max()) * 1.05; identity(ax, 0, hi); ax.set_xlim(0, hi); ax.set_ylim(0, hi)
ax.set_xlabel("RSEM assigned fragments (M)"); ax.set_ylabel("kallisto assigned fragments (M)")
ratio = (k / r).median()
ax.text(0.96, 0.04, f"median ratio\nkallisto / RSEM\n= {ratio:.2f}", transform=ax.transAxes, ha="right",
        va="bottom", color=INK2, fontsize=6.3)
ax.set_xlabel("RSEM fragments (M)"); ax.set_ylabel("kallisto fragments (M)")
panel(ax, "d", "Assigned fragments")
save(fig, "Fig1_mapping_qc")
K["mapping"] = {c: float(qc[c].median()) for c, _, _ in cols} | {
    "reads_median_M": float(qc["star_input_reads"].median() / 1e6), "assigned_ratio_kal_over_rsem": float(ratio)}

# ========================================================================= Fig 2: within-sample agreement
rep = (ps["gene"]["pearson_log2tpm"] - ps["gene"]["pearson_log2tpm"].median()).abs().idxmin()
fig = plt.figure(figsize=(7.4, 4.6))
gs = fig.add_gridspec(2, 3, wspace=0.62, hspace=0.62)
for j, lv in enumerate(("gene", "isoform")):
    ax = fig.add_subplot(gs[0, j])
    idcol = "transcript_id" if lv == "isoform" else "gene_id"
    x = pd.read_csv(f"results/expression/rsem/{lv}.tpm.tsv.gz", sep="\t", index_col=0, usecols=[idcol, rep]).iloc[:, 0]
    y = pd.read_csv(f"results/expression/kallisto/{lv}.tpm.tsv.gz", sep="\t", index_col=0,
                    usecols=[idcol, rep]).iloc[:, 0].reindex(x.index)
    lx, ly = np.log2(x + 1), np.log2(y + 1)
    keep = (x > 0) | (y > 0)
    hb = hexdens(ax, lx[keep], ly[keep], gridsize=80, extent=(0, 16, 0, 16))
    identity(ax, 0, 16); ax.set_xlim(0, 16); ax.set_ylim(0, 16); ax.set_aspect("equal")
    ax.set_xlabel("STAR+RSEM  log2(TPM + 1)"); ax.set_ylabel("kallisto  log2(TPM + 1)")
    r_ = ps[lv].loc[rep, "pearson_log2tpm"]; rho_ = ps[lv].loc[rep, "spearman_tpm"]
    ax.text(0.04, 0.96, f"r = {r_:.3f}\nρ = {rho_:.3f}\nn = {keep.sum():,}", transform=ax.transAxes, va="top",
            fontsize=6.8, color=INK2)
    panel(ax, "ab"[j], f"{lv.capitalize()}, sample {rep}")
    colorbar(fig, hb, ax)
ax = fig.add_subplot(gs[0, 2])
data = [ps["gene"]["pearson_log2tpm"], ps["isoform"]["pearson_log2tpm"], ps["gene"]["spearman_tpm"],
        ps["isoform"]["spearman_tpm"]]
strip_box(ax, [d.values for d in data], [0, 1, 2.5, 3.5], [GENE, ISO, GENE, ISO], width=0.7, s=2.5)
ax.set_xticks([0.5, 3]); ax.set_xticklabels(["Pearson r", "Spearman ρ"])
ax.set_ylabel("Correlation per sample"); ax.set_ylim(0.55, 1.0)
ax.legend(handles=[mpl.patches.Patch(color=GENE, label="gene"), mpl.patches.Patch(color=ISO, label="isoform")],
          loc="lower left", handlelength=1, handleheight=0.8)
panel(ax, "c", "All 450 samples")
# drivers of within-sample agreement
drivers = [("star_multi_pct", "STAR multi-mapped reads (%)", qc), ("star_input_reads", "Read pairs (M)", qc)]
for j, (col, lab, src) in enumerate(drivers):
    ax = fig.add_subplot(gs[1, j])
    xv = src[col].reindex(ps["gene"].index) / (1e6 if col == "star_input_reads" else 1)
    for lv, c in (("gene", GENE), ("isoform", ISO)):
        yv = ps[lv]["pearson_log2tpm"]
        ax.scatter(xv, yv, s=4, color=c, alpha=0.5, lw=0, label=f"{lv} (ρ = {sp(xv.values, yv.values):.2f})")
    ax.set_xlabel(lab); ax.set_ylabel("Pearson r, log2(TPM+1)")
    if col == "star_multi_pct":
        ax.set_xscale("log"); ax.set_xticks([5, 10, 20, 40]); ax.set_xticklabels(["5", "10", "20", "40"])
    ax.legend(loc="lower left", handletextpad=0.2, markerscale=2.5)
    panel(ax, "de"[j], ["vs multi-mapping", "vs sequencing depth"][j])
ax = fig.add_subplot(gs[1, 2])
for lv, c, mk in (("gene", GENE, "o"), ("isoform", ISO, "o")):
    d = ps[lv]
    ax.scatter(d["rsem_detected"] / 1e3, d["kallisto_detected"] / 1e3, s=4, color=c, alpha=0.5, lw=0, label=lv)
lo, hi = 10, 125; identity(ax, lo, hi); ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
for a in (ax.xaxis, ax.yaxis):
    a.set_major_locator(mpl.ticker.FixedLocator([10, 20, 50, 100])); a.set_major_formatter(mpl.ticker.ScalarFormatter())
    a.set_minor_formatter(mpl.ticker.NullFormatter())
ax.set_xlabel("RSEM detected (×1000)"); ax.set_ylabel("kallisto detected (×1000)")
ax.legend(loc="upper left", markerscale=2.5, handletextpad=0.2)
panel(ax, "f", "Detected (TPM ≥ 0.1)")
save(fig, "Fig2_within_sample")
K["within"] = {lv: {"pearson_median": float(ps[lv]["pearson_log2tpm"].median()),
                    "pearson_min": float(ps[lv]["pearson_log2tpm"].min()),
                    "spearman_median": float(ps[lv]["spearman_tpm"].median()),
                    "rsem_detected_median": float(ps[lv]["rsem_detected"].median()),
                    "kallisto_detected_median": float(ps[lv]["kallisto_detected"].median()),
                    "rho_vs_multimap": float(sp(qc["star_multi_pct"].reindex(ps[lv].index).values,
                                                ps[lv]["pearson_log2tpm"].values))}
               for lv in ("gene", "isoform")} | {"representative_sample": rep}

# ========================================================================= Fig 3: across samples (per feature)
expr = {lv: pf[lv].index.isin(tested[("rsem", lv)] | tested[("kallisto", lv)]) for lv in pf}
for lv in pf:
    pf[lv]["mean_log2tpm"] = np.log2((pf[lv]["rsem_median_tpm"] + pf[lv]["kallisto_median_tpm"]) / 2 + 1)
fig = plt.figure(figsize=(7.4, 4.8))
gs = fig.add_gridspec(2, 3, wspace=0.62, hspace=0.62, width_ratios=[1, 1, 1.15])
ax = fig.add_subplot(gs[0, 0])
for lv, c in (("gene", GENE), ("isoform", ISO)):
    v = np.sort(pf[lv].loc[expr[lv], "spearman_tpm"].dropna().values)
    ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=c, lw=1.6, label=f"{lv} (n = {len(v):,})")
    K.setdefault("across", {})[lv] = {"n_expressed": int(len(v)), "median": float(np.median(v)),
                                      "frac_gt_0.9": float((v > 0.9).mean()), "frac_lt_0.5": float((v < 0.5).mean())}
ax.axvline(0.9, color=GRID, lw=0.8, zorder=0)
ax.set_xlim(-0.2, 1); ax.set_ylim(0, 1); ax.set_xlabel("Spearman ρ across 450 samples")
ax.set_ylabel("Cumulative fraction"); ax.legend(loc="upper left", handlelength=1.2)
panel(ax, "a", "Expressed features")
for j, lv in enumerate(("gene", "isoform")):
    ax = fig.add_subplot(gs[0, j + 1])
    d = pf[lv][expr[lv]].dropna(subset=["spearman_tpm"])
    hb = hexdens(ax, d["mean_log2tpm"], d["spearman_tpm"], gridsize=60, extent=(0, 12, -0.4, 1))
    binned_trend(ax, d["mean_log2tpm"].values, d["spearman_tpm"].values, np.arange(0, 12.5, 0.5))
    ax.set_xlim(0, 12); ax.set_ylim(-0.4, 1.02)
    ax.set_xlabel("Median log2(TPM + 1)"); ax.set_ylabel("Spearman ρ across samples")
    panel(ax, "bc"[j], f"{lv.capitalize()} vs expression")
    colorbar(fig, hb, ax)
# isoform complexity
ax = fig.add_subplot(gs[1, 0])
edges = [0, 1, 2, 5, 10, 20, 1e9]; labels = ["1", "2", "3–5", "6–10", "11–20", ">20"]
for k_, (lv, c, col) in enumerate((("gene", GENE, "n_isoforms"), ("isoform", ISO, "n_isoforms_of_gene"))):
    d = pf[lv][expr[lv]].dropna(subset=["spearman_tpm"])
    groups = pd.cut(d[col], edges, labels=labels)
    vals = [d.loc[groups == l, "spearman_tpm"].values for l in labels]
    pos = np.arange(len(labels)) + (k_ - 0.5) * 0.38
    bp = ax.boxplot(vals, positions=pos, widths=0.32, showfliers=False, patch_artist=True,
                    medianprops=dict(color="white", lw=1.2), whiskerprops=dict(color=c, lw=0.7),
                    capprops=dict(color=c, lw=0.7), boxprops=dict(edgecolor=c, lw=0))
    for b in bp["boxes"]:
        b.set_facecolor(c)
ax.set_xticks(range(len(labels))); ax.set_xticklabels(labels, fontsize=5.8, rotation=35, ha="right", rotation_mode="anchor"); ax.set_xlabel("Annotated isoforms per gene")
ax.set_ylabel("Spearman ρ across samples"); ax.set_ylim(-0.4, 1.02)
ax.legend(handles=[mpl.patches.Patch(color=GENE, label="gene"), mpl.patches.Patch(color=ISO, label="isoform")],
          loc="upper center", bbox_to_anchor=(0.5, -0.38), ncol=2, handlelength=1, handleheight=0.8)
panel(ax, "d", "Isoforms per gene")
ax = fig.add_subplot(gs[1, 1])
d = pf["isoform"][expr["isoform"]].dropna(subset=["spearman_tpm"])
lbins = [0, 500, 1000, 2000, 3000, 5000, 1e9]; llab = ["<0.5", "0.5–1", "1–2", "2–3", "3–5", ">5"]
g = pd.cut(d["length"], lbins, labels=llab)
bp = ax.boxplot([d.loc[g == l, "spearman_tpm"].values for l in llab], widths=0.55, showfliers=False,
                patch_artist=True, medianprops=dict(color="white", lw=1.2), whiskerprops=dict(color=ISO, lw=0.7),
                capprops=dict(color=ISO, lw=0.7), boxprops=dict(edgecolor=ISO, lw=0))
for b in bp["boxes"]:
    b.set_facecolor(ISO)
ax.set_xticklabels(llab, fontsize=5.8, rotation=35, ha="right", rotation_mode="anchor"); ax.set_xlabel("Transcript length (kb)"); ax.set_ylabel("Spearman ρ across samples")
ax.set_ylim(-0.4, 1.02); panel(ax, "e", "Isoform length")
ax = fig.add_subplot(gs[1, 2])
top = d["transcript_biotype"].value_counts().index[:6]
vals = [d.loc[d["transcript_biotype"] == b, "spearman_tpm"].values for b in top]
order = np.argsort([np.median(v) for v in vals])[::-1]
bp = ax.boxplot([vals[i] for i in order], vert=False, widths=0.55, showfliers=False, patch_artist=True,
                medianprops=dict(color="white", lw=1.2), whiskerprops=dict(color=ISO, lw=0.7),
                capprops=dict(color=ISO, lw=0.7), boxprops=dict(edgecolor=ISO, lw=0))
for b in bp["boxes"]:
    b.set_facecolor(ISO)
SHORT = {"protein_coding": "protein coding", "retained_intron": "retained intron", "lncRNA": "lncRNA",
         "protein_coding_CDS_not_defined": "CDS not defined", "nonsense_mediated_decay": "NMD",
         "processed_pseudogene": "proc. pseudogene", "processed_transcript": "processed tx"}
ax.set_yticklabels([f"{SHORT.get(top[i], top[i])}\n(n = {len(vals[i]):,})" for i in order], fontsize=5.8)
ax.invert_yaxis(); ax.set_xlabel("Spearman ρ across samples"); ax.set_xlim(-0.4, 1.02)
panel(ax, "f", "Isoform biotype")
save(fig, "Fig3_across_samples")

# ========================================================================= Fig 4: detection and systematic bias
fig = plt.figure(figsize=(7.4, 2.4))
gs = fig.add_gridspec(1, 3, wspace=0.55, width_ratios=[1.1, 1, 1])
ax = fig.add_subplot(gs[0, 0])
for i, lv in enumerate(("gene", "isoform")):
    a, b = tested[("rsem", lv)], tested[("kallisto", lv)]
    stacked(ax, i, [(len(a - b), RSEM), (len(a & b), SHARED), (len(b - a), KAL)], len(a | b))
    K.setdefault("expressed_overlap", {})[lv] = {"rsem_only": len(a - b), "both": len(a & b), "kallisto_only": len(b - a)}
ax.set_yticks([0, 1]); ax.set_yticklabels(["gene", "isoform"]); ax.set_xlim(0, 100)
ax.set_xlabel("% of expressed features (union)"); ax.set_ylim(1.55, -0.65)
ax.legend(handles=[mpl.patches.Patch(color=RSEM, label="RSEM only"), mpl.patches.Patch(color=SHARED, label="both"),
                   mpl.patches.Patch(color=KAL, label="kallisto only")],
          loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=3, handlelength=1, columnspacing=0.8)
panel(ax, "a", "Expressed features")
for j, lv in enumerate(("gene", "isoform")):
    ax = fig.add_subplot(gs[0, j + 1])
    d = pf[lv]
    keep = (d["rsem_median_tpm"] > 0) | (d["kallisto_median_tpm"] > 0)
    d = d[keep]
    A = d["mean_log2tpm"].values
    M = np.log2((d["kallisto_median_tpm"] + 1) / (d["rsem_median_tpm"] + 1)).values
    hb = hexdens(ax, A, M, gridsize=60, extent=(0, 12, -4, 4))
    binned_trend(ax, A, M, np.arange(0, 12.5, 0.5))
    ax.axhline(0, color=MUTED, lw=0.7, zorder=0)
    ax.set_xlim(0, 12); ax.set_ylim(-4, 4)
    ax.set_xlabel("Median log2(TPM + 1)"); ax.set_ylabel("log2 ratio, kallisto / RSEM")
    ax.text(0.97, 0.96, "higher in kallisto", transform=ax.transAxes, ha="right", va="top", fontsize=6.3, color=KAL)
    ax.text(0.97, 0.04, "higher in RSEM", transform=ax.transAxes, ha="right", va="bottom", fontsize=6.3, color=RSEM)
    panel(ax, "bc"[j], f"{lv.capitalize()}: bias vs expression")
    K.setdefault("bias", {})[lv] = {"frac_abs_log2fc_gt1": float((np.abs(M) > 1).mean()),
                                    "frac_abs_log2fc_gt1_expressed": float((np.abs(M[A > 1]) > 1).mean())}
save(fig, "Fig4_detection_bias")

# ========================================================================= Fig 5: eQTL
def egenes(t, lv):
    d = cis[(t, lv)]
    d = d[d["qval"] <= FDR]
    return set(d["group_id"]) if "group_id" in d else set(d.index)


def lead(t, lv):
    d = cis[(t, lv)].copy()
    d["gene"] = d["group_id"] if "group_id" in d else d.index
    return d.set_index("gene")


fig = plt.figure(figsize=(7.4, 4.6))
gs = fig.add_gridspec(2, 3, wspace=0.55, hspace=0.75)
ax = fig.add_subplot(gs[0, 0])
w = 0.36
for k_, t in enumerate(TOOL):
    vals = [len(egenes(t, lv)) for lv in ("gene", "isoform")]
    bars = ax.bar(np.arange(2) + (k_ - 0.5) * (w + 0.03), vals, width=w, color=TOOL[t][1], label=TOOL[t][0])
    for b_, v in zip(bars, vals):
        ax.text(b_.get_x() + w / 2, v + 150, f"{v / 1000:.1f}k", ha="center", va="bottom", fontsize=6, color=INK2)
ax.set_xticks([0, 1]); ax.set_xticklabels(["gene", "isoform\n(per gene)"]); ax.set_ylabel("eGenes (FDR 5%)")
ax.set_ylim(0, 18000); ax.legend(loc="upper left", ncol=2, handlelength=1, columnspacing=0.8)
panel(ax, "a", "cis-eGenes discovered")
ax = fig.add_subplot(gs[0, 1])
for i, lv in enumerate(("gene", "isoform")):
    a, b = egenes("rsem", lv), egenes("kallisto", lv)
    stacked(ax, i, [(len(a - b), RSEM), (len(a & b), SHARED), (len(b - a), KAL)], len(a | b))
    K.setdefault("egene_overlap", {})[lv] = {"rsem_only": len(a - b), "both": len(a & b), "kallisto_only": len(b - a),
                                             "jaccard": len(a & b) / len(a | b)}
ax.set_yticks([0, 1]); ax.set_yticklabels(["gene", "isoform"]); ax.set_xlim(0, 100); ax.set_ylim(1.55, -0.65)
ax.set_xlabel("% of eGenes (union)")
ax.legend(handles=[mpl.patches.Patch(color=RSEM, label="RSEM only"), mpl.patches.Patch(color=SHARED, label="both"),
                   mpl.patches.Patch(color=KAL, label="kallisto only")],
          loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=3, handlelength=1, columnspacing=0.8)
panel(ax, "b", "eGene overlap")
# tool-specific eGenes: q-value in the other tool
ax = fig.add_subplot(gs[0, 2])
bins = np.linspace(0, 0.5, 26)
for t, o in (("rsem", "kallisto"), ("kallisto", "rsem")):
    only = egenes(t, "gene") - egenes(o, "gene")
    q_other = lead(o, "gene").reindex(list(only))["qval"]
    untested = q_other.isna().sum()
    ax.hist(q_other.dropna(), bins=bins, histtype="step", lw=1.4, color=TOOL[t][1],
            label=f"{'RSEM' if t == 'rsem' else 'kallisto'}-only (n = {len(only):,};\n{untested:,} not tested in other)")
    K.setdefault("tool_specific_egenes", {})[t] = {"n": len(only), "untested_in_other": int(untested),
                                                   "other_q_lt_0.2": int((q_other < 0.2).sum())}
ax.axvline(FDR, color=MUTED, lw=0.7)
ax.set_xlim(0, 0.5); ax.set_xlabel("q-value in the other tool"); ax.set_ylabel("eGenes (gene level)")
ax.legend(loc="upper right", fontsize=5.8, handlelength=1.2)
panel(ax, "c", "Tool-specific eGenes")
for j, lv in enumerate(("gene", "isoform")):
    ax = fig.add_subplot(gs[1, j])
    a, b = lead("rsem", lv), lead("kallisto", lv)
    both = sorted(egenes("rsem", lv) & egenes("kallisto", lv))
    x = -np.log10(a.loc[both, "pval_nominal"].clip(lower=1e-300)); y = -np.log10(b.loc[both, "pval_nominal"].clip(lower=1e-300))
    hi = np.quantile(np.r_[x, y], 0.995)
    hb = hexdens(ax, x.clip(upper=hi), y.clip(upper=hi), gridsize=55, extent=(0, hi, 0, hi))
    identity(ax, 0, hi); ax.set_xlim(0, hi); ax.set_ylim(0, hi); ax.set_aspect("equal")
    same = (a.loc[both, "variant_id"] == b.loc[both, "variant_id"]).mean()
    ax.text(0.04, 0.96, f"ρ = {sp(x.values, y.values):.2f}\nsame lead SNP: {same:.0%}\nn = {len(both):,}",
            transform=ax.transAxes, va="top", fontsize=6.8, color=INK2)
    ax.set_xlabel("RSEM  −log10 p, lead SNP"); ax.set_ylabel("kallisto  −log10 p, lead SNP")
    panel(ax, "de"[j], f"{lv.capitalize()}: shared eGenes")
    K.setdefault("lead_concordance", {})[lv] = {"spearman": float(sp(x.values, y.values)), "same_lead_snp": float(same)}
ax = fig.add_subplot(gs[1, 2])
for k_, t in enumerate(TOOL):
    ind = pd.read_csv(f"results/eqtl/{t}/gene/cis_independent.tsv.gz", sep="\t")
    n = ind.groupby("phenotype_id").size().clip(upper=5).value_counts().sort_index()
    n = n.reindex(range(1, 6), fill_value=0)
    ax.bar(np.arange(1, 6) + (k_ - 0.5) * 0.4, n.values / n.sum() * 100, width=0.38, color=TOOL[t][1], label=TOOL[t][0])
ax.set_xticks(range(1, 6)); ax.set_xticklabels(["1", "2", "3", "4", "≥5"])
ax.set_xlabel("Independent signals per eGene"); ax.set_ylabel("% of eGenes (gene level)")
ax.legend(loc="upper right", handlelength=1)
panel(ax, "f", "Allelic heterogeneity")
save(fig, "Fig5_eqtl")
K["eqtl"] = {f"{t}_{lv}": {"tested": int(len(cis[(t, lv)])), "egenes": len(egenes(t, lv))}
             for t in TOOL for lv in ("gene", "isoform")}

# ========================================================================= Fig 6: compute
fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.3), gridspec_kw=dict(wspace=0.55))
steps = [("star", "STAR", RSEM), ("rsem", "RSEM", RSEM), ("star_rsem", "STAR+\nRSEM", RSEM), ("kallisto", "kallisto", KAL)]
for ax, (unit, lab, key, scale) in zip(axs[:2], [("min", "Wall time per sample (min)", "wall_s", 60),
                                                  ("h", "CPU time per sample (h)", "cpu_s", 3600)]):
    data = [comp[f"{s}_{key}"].values / scale for s, _, _ in steps]
    strip_box(ax, data, range(4), [c for _, _, c in steps], s=2.5)
    ax.set_yscale("log"); ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l, _ in steps]); ax.set_ylabel(lab)
    lo_, hi_ = ax.get_ylim(); ax.set_ylim(lo_, hi_ * 4)
    med = [np.nanmedian(d) for d in data]
    ax.text(0.98, 0.98, f"STAR+RSEM: {med[2] / med[3]:.0f}× kallisto", transform=ax.transAxes, ha="right", va="top",
            fontsize=6.3, color=INK2)
    K.setdefault("compute", {})[key] = {s: float(m) for (s, _, _), m in zip(steps, med)}
panel(axs[0], "a", "Wall time"); panel(axs[1], "b", "CPU time")
ax = axs[2]
mem = [comp[f"{s}_max_rss_mb"].median() / 1024 for s in ("star", "rsem", "kallisto")]
bars = ax.bar(range(3), mem, color=[RSEM, RSEM, KAL], width=0.6)
for b_, v in zip(bars, mem):
    ax.text(b_.get_x() + 0.3, v + 0.5, f"{v:.1f}", ha="center", va="bottom", fontsize=6.5, color=INK2)
ax.set_xticks(range(3)); ax.set_xticklabels(["STAR", "RSEM", "kallisto"]); ax.set_ylabel("Peak memory (GB, median)")
panel(ax, "c", "Memory")
save(fig, "Fig6_compute")
K["compute"]["max_rss_gb"] = dict(zip(("star", "rsem", "kallisto"), map(float, mem)))
K["compute"]["total_cpu_h"] = {"star_rsem": float(comp["star_rsem_cpu_s"].sum() / 3600),
                               "kallisto": float(comp["kallisto_cpu_s"].sum() / 3600)}

# ========================================================================= Fig 0: one-page summary (README)
fig = plt.figure(figsize=(7.4, 4.7))
gs = fig.add_gridspec(2, 3, wspace=0.62, hspace=0.75)
# a. within-sample agreement, all samples
ax = fig.add_subplot(gs[0, 0])
data = [ps["gene"]["pearson_log2tpm"], ps["isoform"]["pearson_log2tpm"], ps["gene"]["spearman_tpm"],
        ps["isoform"]["spearman_tpm"]]
strip_box(ax, [d.values for d in data], [0, 1, 2.5, 3.5], [GENE, ISO, GENE, ISO], width=0.7, s=2.5)
ax.set_xticks([0.5, 3]); ax.set_xticklabels(["Pearson r\nlog2 TPM", "Spearman ρ\nTPM"])
ax.set_ylabel("RSEM vs kallisto, per sample"); ax.set_ylim(0.55, 1.0)
ax.legend(handles=[mpl.patches.Patch(color=GENE, label="gene"), mpl.patches.Patch(color=ISO, label="isoform")],
          loc="lower left", handlelength=1, handleheight=0.8)
panel(ax, "a", "Within samples (n = 450)")
# b. across-sample agreement per feature
ax = fig.add_subplot(gs[0, 1])
for lv, c in (("gene", GENE), ("isoform", ISO)):
    v = np.sort(pf[lv].loc[expr[lv], "spearman_tpm"].dropna().values)
    ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=c, lw=1.6, label=f"{lv} (median {np.median(v):.2f})")
ax.axvline(0.9, color=GRID, lw=0.8, zorder=0)
ax.set_xlim(-0.2, 1); ax.set_ylim(0, 1); ax.set_xlabel("Spearman ρ across samples")
ax.set_ylabel("Cumulative fraction of features"); ax.legend(loc="upper left", handlelength=1.2)
panel(ax, "b", "Across samples")
# c. isoform agreement vs expression
ax = fig.add_subplot(gs[0, 2])
d = pf["isoform"][expr["isoform"]].dropna(subset=["spearman_tpm"])
hb = hexdens(ax, d["mean_log2tpm"], d["spearman_tpm"], gridsize=50, extent=(0, 12, -0.4, 1))
binned_trend(ax, d["mean_log2tpm"].values, d["spearman_tpm"].values, np.arange(0, 12.5, 0.5))
ax.set_xlim(0, 12); ax.set_ylim(-0.4, 1.02)
ax.set_xlabel("Median log2(TPM + 1)"); ax.set_ylabel("Spearman ρ across samples")
panel(ax, "c", "Isoforms vs expression"); colorbar(fig, hb, ax)
# d. eGene overlap
ax = fig.add_subplot(gs[1, 0])
for i, lv in enumerate(("gene", "isoform")):
    a, b = egenes("rsem", lv), egenes("kallisto", lv)
    stacked(ax, i, [(len(a - b), RSEM), (len(a & b), SHARED), (len(b - a), KAL)], len(a | b))
ax.set_yticks([0, 1]); ax.set_yticklabels(["gene", "isoform"]); ax.set_xlim(0, 100); ax.set_ylim(1.55, -0.65)
ax.set_xlabel("% of cis-eGenes (union, FDR 5%)")
ax.legend(handles=[mpl.patches.Patch(color=RSEM, label="RSEM only"), mpl.patches.Patch(color=SHARED, label="both"),
                   mpl.patches.Patch(color=KAL, label="kallisto only")],
          loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=3, handlelength=1, columnspacing=0.6)
panel(ax, "d", "eGene overlap")
# e. lead-SNP strength for shared eGenes (gene level)
ax = fig.add_subplot(gs[1, 1])
a, b = lead("rsem", "gene"), lead("kallisto", "gene")
both = sorted(egenes("rsem", "gene") & egenes("kallisto", "gene"))
x = -np.log10(a.loc[both, "pval_nominal"].clip(lower=1e-300)); y = -np.log10(b.loc[both, "pval_nominal"].clip(lower=1e-300))
hi = np.quantile(np.r_[x, y], 0.995)
hexdens(ax, x.clip(upper=hi), y.clip(upper=hi), gridsize=45, extent=(0, hi, 0, hi))
identity(ax, 0, hi); ax.set_xlim(0, hi); ax.set_ylim(0, hi); ax.set_aspect("equal")
ax.text(0.04, 0.96, f"ρ = {sp(x.values, y.values):.2f}\nn = {len(both):,}", transform=ax.transAxes, va="top",
        fontsize=6.8, color=INK2)
ax.set_xlabel("RSEM  −log10 p"); ax.set_ylabel("kallisto  −log10 p")
panel(ax, "e", "Shared eGenes, lead SNP")
# f. wall time per sample
ax = fig.add_subplot(gs[1, 2])
steps_ = [("star", "STAR", RSEM), ("rsem", "RSEM", RSEM), ("kallisto", "kallisto", KAL)]
strip_box(ax, [comp[f"{s_}_wall_s"].values / 60 for s_, _, _ in steps_], range(3), [c for _, _, c in steps_], s=2.5)
ax.set_yscale("log"); lo_, hi_ = ax.get_ylim(); ax.set_ylim(lo_, hi_ * 4)
ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l, _ in steps_]); ax.set_ylabel("Wall time per sample (min)")
ax.text(0.98, 0.98, f"STAR+RSEM: {K['compute']['wall_s']['star_rsem'] / K['compute']['wall_s']['kallisto']:.0f}× kallisto",
        transform=ax.transAxes, ha="right", va="top", fontsize=6.3, color=INK2)
panel(ax, "f", "Compute")
save(fig, "Fig0_summary")

with open(f"{OUT}/key_numbers.json", "w") as f:
    json.dump(K, f, indent=1)
print(json.dumps(K, indent=1))
