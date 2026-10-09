"""Gene and transcript tables from an Ensembl GTF.

genes.tsv:       gene_id, gene_name, gene_biotype, chr, start, end, strand, tss (1-based)
transcripts.tsv: transcript_id, gene_id, transcript_name, transcript_biotype, chr, start, end, strand,
                 tss (1-based), length (sum of exon lengths), n_isoforms_of_gene
"""
import re
import sys
from collections import defaultdict

import pandas as pd

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

ATTR = re.compile(r'(\S+) "([^"]*)"')

genes, transcripts = [], []
exon_length = defaultdict(int)
with open(snakemake.input.gtf) as f:
    for line in f:
        if line.startswith("#"):
            continue
        chrom, _, feature, start, end, _, strand, _, attributes = line.rstrip("\n").split("\t")
        if feature not in ("gene", "transcript", "exon"):
            continue
        a = dict(ATTR.findall(attributes))
        start, end = int(start), int(end)
        if feature == "gene":
            genes.append((a["gene_id"], a.get("gene_name", a["gene_id"]), a.get("gene_biotype", ""),
                          chrom, start, end, strand))
        elif feature == "transcript":
            transcripts.append((a["transcript_id"], a["gene_id"], a.get("transcript_name", a["transcript_id"]),
                                a.get("transcript_biotype", ""), chrom, start, end, strand))
        else:
            exon_length[a["transcript_id"]] += end - start + 1

genes = pd.DataFrame(genes, columns=["gene_id", "gene_name", "gene_biotype", "chr", "start", "end", "strand"])
genes["tss"] = genes["start"].where(genes["strand"] == "+", genes["end"])

tx = pd.DataFrame(transcripts, columns=["transcript_id", "gene_id", "transcript_name", "transcript_biotype",
                                        "chr", "start", "end", "strand"])
tx["tss"] = tx["start"].where(tx["strand"] == "+", tx["end"])
tx["length"] = tx["transcript_id"].map(exon_length).astype(int)
tx["n_isoforms_of_gene"] = tx.groupby("gene_id")["transcript_id"].transform("size")

assert genes["gene_id"].is_unique and tx["transcript_id"].is_unique
assert tx["gene_id"].isin(genes["gene_id"]).all(), "transcripts with a gene_id missing from gene records"
print(f"{len(genes)} genes, {len(tx)} transcripts")

genes.to_csv(snakemake.output.genes, sep="\t", index=False)
tx.to_csv(snakemake.output.transcripts, sep="\t", index=False)
