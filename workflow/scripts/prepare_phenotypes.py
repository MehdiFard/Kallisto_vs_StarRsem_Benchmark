"""Normalized phenotype BED for tensorQTL (GTEx expression pipeline approach).

1. keep features on the selected chromosomes with TPM >= t and counts >= c in >= f of samples
2. TMM-normalized CPM (library sizes from all features), 3. rank-based inverse normal transform per feature
Isoforms are placed at their gene's TSS and written next to each other, with an isoform -> gene
groups file, so tensorQTL can map them grouped by gene.
"""
import sys

import numpy as np
import pandas as pd

sys.stdout = sys.stderr = open(snakemake.log[0], "w")
sys.path.insert(0, snakemake.scriptdir)
import rnaseqnorm  # noqa: E402

level = snakemake.wildcards.level
p = snakemake.params
chroms = [str(c) for c in p.chromosomes]

counts = pd.read_csv(snakemake.input.counts, sep="\t", index_col=0)
tpm = pd.read_csv(snakemake.input.tpm, sep="\t", index_col=0)
assert counts.index.equals(tpm.index) and counts.columns.equals(tpm.columns)
n_samples = counts.shape[1]

genes = pd.read_csv(snakemake.input.genes, sep="\t", index_col="gene_id", dtype={"chr": str})
if level == "gene":
    gene_of = pd.Series(counts.index, index=counts.index)
else:
    tx = pd.read_csv(snakemake.input.transcripts, sep="\t", index_col="transcript_id", dtype={"chr": str})
    gene_of = tx.loc[counts.index, "gene_id"]
chrom = genes.loc[gene_of, "chr"].set_axis(counts.index)
tss = genes.loc[gene_of, "tss"].set_axis(counts.index)

min_samples = p.sample_frac_threshold * n_samples
expressed = ((tpm >= p.tpm_threshold).sum(axis=1) >= min_samples) & \
            ((counts >= p.count_threshold).sum(axis=1) >= min_samples)
on_chrom = chrom.isin(chroms)
keep = expressed & on_chrom
print(f"{level}: {len(counts)} features, {expressed.sum()} expressed, "
      f"{keep.sum()} expressed on chromosomes {chroms[0]}..{chroms[-1]}")

cpm = rnaseqnorm.edgeR_cpm(counts, normalized_lib_sizes=True).loc[keep]
cpm = cpm[cpm.nunique(axis=1) > 1]  # constant rows give a constant INT
norm = rnaseqnorm.inverse_normal_transform(cpm)
print(f"{level}: {len(norm)} phenotypes written")

bed = pd.DataFrame({
    "#chr": chrom[norm.index].to_numpy(),
    "start": tss[norm.index].to_numpy() - 1,
    "end": tss[norm.index].to_numpy(),
    "phenotype_id": norm.index.to_numpy(),
    "gene_id": gene_of[norm.index].to_numpy(),
})
bed["chr_order"] = bed["#chr"].map({c: i for i, c in enumerate(chroms)})
# gene_id as a tie-breaker keeps isoforms of one gene adjacent (required for grouped mapping)
bed = bed.sort_values(["chr_order", "start", "gene_id", "phenotype_id"])

bed[["phenotype_id", "gene_id"]].to_csv(snakemake.output.groups, sep="\t", index=False, header=False)
out = pd.concat([bed[["#chr", "start", "end", "phenotype_id"]].reset_index(drop=True),
                 norm.loc[bed["phenotype_id"]].reset_index(drop=True)], axis=1)
out.to_csv(snakemake.output.bed, sep="\t", index=False, float_format="%.6g")
