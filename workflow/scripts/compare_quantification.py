"""RSEM vs kallisto agreement at one level (gene or isoform).

per_sample:  correlation between tools within each sample (across features)
per_feature: correlation between tools for each feature (across samples) + expression and annotation
"""
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

level = snakemake.wildcards.level
inp = snakemake.input


def load(path):
    return pd.read_csv(path, sep="\t", index_col=0)


r_cnt, r_tpm = load(inp.rsem_counts), load(inp.rsem_tpm)
k_cnt, k_tpm = load(inp.kallisto_counts), load(inp.kallisto_tpm)

# Both tools must report the same features and samples (same transcript set by construction)
for name, df in [("rsem tpm", r_tpm), ("kallisto counts", k_cnt), ("kallisto tpm", k_tpm)]:
    if not (df.index.sort_values().equals(r_cnt.index.sort_values()) and df.columns.equals(r_cnt.columns)):
        sys.exit(f"{name}: features/samples differ from rsem counts")
k_cnt, k_tpm, r_tpm = (df.loc[r_cnt.index] for df in (k_cnt, k_tpm, r_tpm))

DETECT_TPM = 0.1


def rowwise_pearson(a, b):
    a = a - a.mean(axis=1, keepdims=True)
    b = b - b.mean(axis=1, keepdims=True)
    denom = np.sqrt((a**2).sum(axis=1) * (b**2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return (a * b).sum(axis=1) / denom


# ---- per sample (columns)
log_r, log_k = np.log2(r_tpm + 1), np.log2(k_tpm + 1)
rows = []
for s in r_cnt.columns:
    rows.append({
        "sample_id": s,
        "pearson_log2tpm": stats.pearsonr(log_r[s], log_k[s])[0],
        "spearman_tpm": stats.spearmanr(r_tpm[s], k_tpm[s])[0],
        "pearson_log2count": stats.pearsonr(np.log2(r_cnt[s] + 1), np.log2(k_cnt[s] + 1))[0],
        "rsem_detected": int((r_tpm[s] >= DETECT_TPM).sum()),
        "kallisto_detected": int((k_tpm[s] >= DETECT_TPM).sum()),
        "both_detected": int(((r_tpm[s] >= DETECT_TPM) & (k_tpm[s] >= DETECT_TPM)).sum()),
        "rsem_total_counts": r_cnt[s].sum(),
        "kallisto_total_counts": k_cnt[s].sum(),
    })
pd.DataFrame(rows).to_csv(snakemake.output.per_sample, sep="\t", index=False, float_format="%.5g")

# ---- per feature (rows)
rank_r = stats.rankdata(r_tpm.values, axis=1)
rank_k = stats.rankdata(k_tpm.values, axis=1)
feat = pd.DataFrame(index=r_cnt.index)
feat["spearman_tpm"] = rowwise_pearson(rank_r, rank_k)
feat["pearson_log2tpm"] = rowwise_pearson(log_r.values, log_k.values)
feat["rsem_median_tpm"] = r_tpm.median(axis=1)
feat["kallisto_median_tpm"] = k_tpm.median(axis=1)
feat["rsem_median_count"] = r_cnt.median(axis=1)
feat["kallisto_median_count"] = k_cnt.median(axis=1)
feat["rsem_frac_detected"] = (r_tpm >= DETECT_TPM).mean(axis=1)
feat["kallisto_frac_detected"] = (k_tpm >= DETECT_TPM).mean(axis=1)

if level == "isoform":
    tx = pd.read_csv(inp.transcripts, sep="\t", index_col="transcript_id", dtype={"chr": str})
    feat = feat.join(tx[["gene_id", "transcript_biotype", "length", "n_isoforms_of_gene"]])
else:
    genes = pd.read_csv(inp.genes, sep="\t", index_col="gene_id", dtype={"chr": str})
    tx = pd.read_csv(inp.transcripts, sep="\t", dtype={"chr": str})
    feat = feat.join(genes[["gene_name", "gene_biotype"]])
    feat["n_isoforms"] = tx.groupby("gene_id").size().reindex(feat.index)
feat.index.name = f"{level}_id"
feat.to_csv(snakemake.output.per_feature, sep="\t", float_format="%.5g")
print(f"{level}: {len(feat)} features, {r_cnt.shape[1]} samples; "
      f"median per-feature spearman = {feat['spearman_tpm'].median():.3f}")
