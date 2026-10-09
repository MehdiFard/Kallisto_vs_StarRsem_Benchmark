# results/

Pipeline outputs. The folders are in git (each holds an empty `.gitkeep`); their contents are not, except the
figure script and PNG figures in `report/`.

| Folder | Content |
|---|---|
| `star/<sample>/` | STAR transcriptome BAM (input to RSEM), `Log.final.out`, `ReadsPerGene.out.tab`, `SJ.out.tab` |
| `rsem/<sample>/` | `isoforms.results.gz`, `genes.results.gz`, `<sample>.stat/` |
| `kallisto/<sample>/` | `abundance.tsv`, `abundance.h5`, `run_info.json` |
| `expression/{rsem,kallisto}/` | `{gene,isoform}.{counts,tpm}.tsv.gz`: features × samples (gene = sum of isoforms) |
| `genotype/` | `genotypes.{pgen,pvar,psam}` (QC'd SNPs of the RNA-seq samples), `pca/` (genotype PCs) |
| `phenotypes/{tool}/` | `{level}.bed.gz` normalized tensorQTL input, `{level}.groups.tsv` isoform → gene |
| `covariates/{tool}/` | `{level}.covariates.tsv`: genotype PCs + expression PCs |
| `eqtl/{tool}/{level}/` | `cis.tsv.gz` (permutation pass, q-values), `cis_independent.tsv.gz`, `cis_signif_pairs.tsv.gz`, `nominal/` |
| `summary/` | `mapping_qc.tsv`, `compute_benchmark.tsv`, `quant_correlation_{gene,isoform}_per_{sample,feature}`, `eqtl_summary.tsv`, `eqtl_egene_overlap.tsv` |
| `report/` | Not made by the pipeline: `make_figures.py` reads the tables above and writes `figures/Fig*.png` and `figures/key_numbers.json` |
