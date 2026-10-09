# data/

Pipeline inputs. Only this file and `metadata/` are in git; download the rest (or symlink existing copies) to the
paths below, or change the paths in `config/config.yaml`.

| Path | Content | Source |
|---|---|---|
| `fastq/{run}_{1,2}.fastq.gz` | Geuvadis mRNA-seq, paired-end 2×75 bp, one pair per sequencing run (621 runs for the 450 samples in `config/samples.tsv`) | ENA study PRJEB3366 (E-GEUV-1) |
| `genotype_1kgp30x/chr{1..22}.vcf.gz` (+ `.tbi`) | 1000 Genomes 30x on GRCh38 (NYGC), all 3,202 samples, per chromosome (`20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr{N}.recalibrated_variants.vcf.gz`, renamed or symlinked) | ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/ |
| `reference/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa` | GRCh38 primary assembly, soft-masked | Ensembl release 110 |
| `reference/Homo_sapiens.GRCh38.110.chr.gtf` | Ensembl 110 annotation (chromosomes only) | Ensembl release 110 |
| `reference/trimmomatic_adapters_overrepseq.fasta` | adapter FASTA for Trimmomatic; only for `trimming: true` | e.g. `adapters/TruSeq3-PE-2.fa` from Trimmomatic, saved under this name (or set `trimmomatic.adapters` in `config/config.yaml`) |
| `metadata/geuvadis_fastq_ena.tsv` | ENA download URL and md5 of every Geuvadis FASTQ, with sample and population | IGSR data portal (Geuvadis data collection) |
| `metadata/igsr_samples.tsv` | population, superpopulation and sex of the samples (sample annotation; not read by the pipeline) | IGSR data portal, Geuvadis ∩ 1kGP 30x |

Sample sheet: `config/samples.tsv` (sample → comma-separated run accessions). Samples are the Geuvadis individuals
with 1kGP 30x genotypes (EUR + YRI).
