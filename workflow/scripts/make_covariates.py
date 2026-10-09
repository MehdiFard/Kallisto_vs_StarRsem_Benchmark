"""tensorQTL covariates: top genotype PCs (plink2 --pca) + top expression PCs of the phenotype BED.

Output: one row per covariate, one column per sample (first column "ID").
"""
import sys

import numpy as np
import pandas as pd

sys.stdout = sys.stderr = open(snakemake.log[0], "w")

n_geno = snakemake.params.n_genotype_pcs
n_expr = snakemake.params.n_expression_pcs

bed = pd.read_csv(snakemake.input.bed, sep="\t", index_col=3)
expr = bed.iloc[:, 3:]  # features x samples
samples = expr.columns

# expression PCs: samples are observations, features (already INT, ~unit variance) are variables
x = expr.T.values
x = x - x.mean(axis=0)
u, s, _ = np.linalg.svd(x, full_matrices=False)
expr_pcs = pd.DataFrame(u[:, :n_expr] * s[:n_expr], index=samples,
                        columns=[f"expr_PC{i + 1}" for i in range(n_expr)])
explained = s**2 / np.sum(s**2)
print("variance explained by expression PCs:", np.round(explained[:n_expr], 4).tolist())

eig = pd.read_csv(snakemake.input.eigenvec, sep="\t")
id_col = "IID" if "IID" in eig.columns else "#IID"
eig = eig.set_index(id_col)
missing = samples.difference(eig.index)
if len(missing):
    sys.exit(f"{len(missing)} samples have no genotype PCs, e.g. {list(missing[:5])}")
geno_pcs = eig.loc[samples, [f"PC{i + 1}" for i in range(n_geno)]]
geno_pcs.columns = [f"geno_PC{i + 1}" for i in range(n_geno)]

covariates = pd.concat([geno_pcs, expr_pcs], axis=1).T
covariates.index.name = "ID"
covariates.to_csv(snakemake.output[0], sep="\t", float_format="%.6g")
print(f"{covariates.shape[0]} covariates x {covariates.shape[1]} samples")
