#!/usr/bin/env python3
"""tensorQTL wrapper: cis (permutations + q-values), cis_independent, cis_nominal (+ significant pairs), trans.

Example:
  python run_tensorqtl.py cis --genotypes results/genotype/genotypes \
      --phenotypes results/phenotypes/rsem/gene.bed.gz --covariates results/covariates/rsem/gene.covariates.tsv \
      --output cis.tsv.gz
"""
import argparse
import os
import sys

import pandas as pd
import torch
import tensorqtl
from tensorqtl import cis, genotypeio, post, trans

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("mode", choices=["cis", "cis_independent", "cis_nominal", "trans"])
parser.add_argument("--genotypes", required=True, help="PLINK2 pgen/pvar/psam prefix")
parser.add_argument("--phenotypes", required=True, help="phenotype BED(.gz)")
parser.add_argument("--covariates", required=True, help="covariates x samples TSV")
parser.add_argument("--groups", help="phenotype_id -> group (gene) TSV without header")
parser.add_argument("--output", required=True)
parser.add_argument("--cis-results", help="output of the cis mode (cis_independent, cis_nominal)")
parser.add_argument("--nominal-dir", help="cis_nominal: directory for per-chromosome parquet files")
parser.add_argument("--prefix", default="tensorqtl", help="cis_nominal: file prefix")
parser.add_argument("--window", type=int, default=1_000_000)
parser.add_argument("--maf-threshold", type=float, default=0.01)
parser.add_argument("--nperm", type=int, default=10_000)
parser.add_argument("--seed", type=int, default=None)
parser.add_argument("--fdr", type=float, default=0.05)
parser.add_argument("--qvalue-lambda", type=float, default=None)
parser.add_argument("--pval-threshold", type=float, default=1e-5, help="trans")
parser.add_argument("--trans-cis-window", type=int, default=5_000_000, help="trans: drop pairs within this window")
parser.add_argument("--batch-size", type=int, default=20_000, help="trans")
parser.add_argument("--allow-cpu", action="store_true", help="run on CPU if no GPU is visible (default: fail)")
args = parser.parse_args()

if not torch.cuda.is_available() and not args.allow_cpu:
    sys.exit("no GPU visible to torch; run on the GPU partition (or pass --allow-cpu)")
print(f"tensorQTL {tensorqtl.__version__}, device: {'GPU ' + torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'}",
      flush=True)

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(args.phenotypes)
covariates_df = pd.read_csv(args.covariates, sep="\t", index_col=0).T
assert covariates_df.index.equals(phenotype_df.columns), "covariate and phenotype samples differ"
group_s = None
if args.groups:
    group_s = pd.read_csv(args.groups, sep="\t", header=None, index_col=0).squeeze("columns")
    group_s = group_s.loc[phenotype_df.index]

genotype_df, variant_df = genotypeio.load_genotypes(args.genotypes, select_samples=phenotype_df.columns.tolist())
genotype_df = genotype_df[phenotype_df.columns]  # same sample order as phenotypes

# only phenotypes on genotyped chromosomes
on_geno_chrom = phenotype_pos_df["chr"].isin(variant_df["chrom"].astype(str).unique())
if not on_geno_chrom.all():
    print(f"dropping {(~on_geno_chrom).sum()} phenotypes on chromosomes without genotypes")
    phenotype_df, phenotype_pos_df = phenotype_df[on_geno_chrom], phenotype_pos_df[on_geno_chrom]
    if group_s is not None:
        group_s = group_s.loc[phenotype_df.index]

common = dict(covariates_df=covariates_df, maf_threshold=args.maf_threshold, window=args.window)

if args.mode == "cis":
    res = cis.map_cis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, group_s=group_s,
                      nperm=args.nperm, seed=args.seed, **common)
    post.calculate_qvalues(res, fdr=args.fdr, qvalue_lambda=args.qvalue_lambda)
    res.to_csv(args.output, sep="\t", float_format="%.6g")

elif args.mode == "cis_independent":
    cis_df = pd.read_csv(args.cis_results, sep="\t", index_col=0)
    if (cis_df["qval"] <= args.fdr).any():
        res = cis.map_independent(genotype_df, variant_df, cis_df, phenotype_df, phenotype_pos_df,
                                  group_s=group_s, fdr=args.fdr, nperm=args.nperm, seed=args.seed, **common)
    else:
        print(f"no phenotype with qval <= {args.fdr}; writing an empty table")
        res = pd.DataFrame()
    res.to_csv(args.output, sep="\t", index=False, float_format="%.6g")

elif args.mode == "cis_nominal":
    cis_df = pd.read_csv(args.cis_results, sep="\t", index_col=0)
    os.makedirs(args.nominal_dir, exist_ok=True)
    cis.map_nominal(genotype_df, variant_df, phenotype_df, phenotype_pos_df, args.prefix,
                    group_s=group_s, output_dir=args.nominal_dir, **common)
    signif = post.get_significant_pairs(cis_df, os.path.join(args.nominal_dir, args.prefix),
                                        group_s=group_s, fdr=args.fdr)
    signif.to_csv(args.output, sep="\t", index=False, float_format="%.6g")

elif args.mode == "trans":
    pairs = trans.map_trans(genotype_df, phenotype_df, covariates_df=covariates_df, return_sparse=True,
                            pval_threshold=args.pval_threshold, maf_threshold=args.maf_threshold,
                            batch_size=args.batch_size)
    pairs = trans.filter_cis(pairs, phenotype_pos_df, variant_df, window=args.trans_cis_window)
    pairs.to_csv(args.output, sep="\t", index=False, float_format="%.6g")

print("done", flush=True)
