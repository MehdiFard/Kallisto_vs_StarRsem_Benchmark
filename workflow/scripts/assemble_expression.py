"""Combine per-sample RSEM or kallisto output into isoform and gene matrices (features x samples).

Gene values are the sum over the gene's isoforms for both tools, so gene-level differences come
from quantification only (RSEM's own genes.results is the same sum, kept in results/rsem/).
"""
import sys

import pandas as pd

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

tool = snakemake.wildcards.tool
samples = snakemake.params.samples
id_col, count_col, tpm_col = {
    "rsem": ("transcript_id", "expected_count", "TPM"),
    "kallisto": ("target_id", "est_counts", "tpm"),
}[tool]

counts, tpm, index = {}, {}, None
for sample, path in zip(samples, snakemake.input.quant):
    d = pd.read_csv(path, sep="\t", usecols=[id_col, count_col, tpm_col], index_col=id_col)
    if index is None:
        index = d.index
    elif not d.index.equals(index):
        sys.exit(f"{path}: transcript order differs from {snakemake.input.quant[0]}")
    counts[sample] = d[count_col]
    tpm[sample] = d[tpm_col]

counts = pd.DataFrame(counts)
tpm = pd.DataFrame(tpm)
counts.index.name = tpm.index.name = "transcript_id"

tx = pd.read_csv(snakemake.input.transcripts, sep="\t", usecols=["transcript_id", "gene_id"],
                 index_col="transcript_id")
unknown = counts.index.difference(tx.index)
if len(unknown):
    sys.exit(f"{len(unknown)} quantified transcripts are not in the annotation, e.g. {list(unknown[:5])}")
print(f"{tool}: {counts.shape[0]} transcripts x {counts.shape[1]} samples "
      f"({len(tx.index.difference(counts.index))} annotated transcripts not quantified)")

gene_of = tx.loc[counts.index, "gene_id"]
gene_counts = counts.groupby(gene_of).sum()
gene_tpm = tpm.groupby(gene_of).sum()
gene_counts.index.name = gene_tpm.index.name = "gene_id"
print(f"{tool}: {gene_counts.shape[0]} genes")

counts.to_csv(snakemake.output.isoform_counts, sep="\t", float_format="%.6g")
tpm.to_csv(snakemake.output.isoform_tpm, sep="\t", float_format="%.6g")
gene_counts.to_csv(snakemake.output.gene_counts, sep="\t", float_format="%.6g")
gene_tpm.to_csv(snakemake.output.gene_tpm, sep="\t", float_format="%.6g")
