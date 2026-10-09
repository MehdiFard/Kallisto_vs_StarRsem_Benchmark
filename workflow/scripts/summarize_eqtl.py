"""cis-eQTL discovery per tool and level, and eGene overlap between RSEM and kallisto.

For isoforms (mapped grouped by gene) an eGene is a gene with a significant isoform group.
"""
import sys

import pandas as pd

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

tools, levels, fdr = snakemake.params.tools, snakemake.params.levels, snakemake.params.fdr
cis_files = iter(snakemake.input.cis)
indep_files = iter(snakemake.input.independent)
signif_files = iter(snakemake.input.signif)


def gene_ids(df):
    return df["group_id"] if "group_id" in df.columns else df.index.to_series()


rows, egenes = [], {}
for tool in tools:  # same order as the expand() in the rule
    for level in levels:
        cis = pd.read_csv(next(cis_files), sep="\t", index_col=0)
        try:
            indep = pd.read_csv(next(indep_files), sep="\t")
        except pd.errors.EmptyDataError:
            indep = pd.DataFrame()
        signif = pd.read_csv(next(signif_files), sep="\t")
        sig = cis[cis["qval"] <= fdr]
        egenes[(tool, level)] = set(gene_ids(sig))
        n_signals = indep.groupby(indep["group_id"] if "group_id" in indep else indep["phenotype_id"]).size() \
            if len(indep) else pd.Series(dtype=int)
        rows.append({
            "tool": tool,
            "level": level,
            "n_tested": len(cis),
            "n_egenes": len(egenes[(tool, level)]),
            "pct_egenes": 100 * len(sig) / len(cis),
            "n_independent_signals": len(indep),
            "n_egenes_multiple_signals": int((n_signals > 1).sum()),
            "n_signif_pairs": len(signif),
            "n_signif_variants": signif["variant_id"].nunique() if len(signif) else 0,
        })
pd.DataFrame(rows).to_csv(snakemake.output.summary, sep="\t", index=False, float_format="%.4g")

rows = []
for level in levels:
    a, b = egenes[("rsem", level)], egenes[("kallisto", level)]
    union = a | b
    rows.append({"level": level, "rsem_only": len(a - b), "kallisto_only": len(b - a), "both": len(a & b),
                 "jaccard": len(a & b) / len(union) if union else float("nan")})
pd.DataFrame(rows).to_csv(snakemake.output.overlap, sep="\t", index=False, float_format="%.4g")
print(pd.DataFrame(rows).to_string(index=False))
