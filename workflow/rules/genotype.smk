## Genotypes for the 450 RNA-seq samples from the complete 1kGP 30x per-chromosome VCFs.
## All QC filters are per-variant, so QC per chromosome and then merging = QC on the merged set.


rule genotype_sample_list:
    output:
        "results/genotype/samples.txt",
    run:
        with open(output[0], "w") as f:
            f.write("\n".join(SAMPLES) + "\n")


rule plink_chrom:
    """One chromosome, straight from the 3202-sample VCF: keep our samples + PASS biallelic SNPs, QC, write pgen.
    plink2 reads the VCF multithreaded and parses only GT, so no intermediate VCF is written. MAF / geno / HWE are
    computed on the kept samples. Priority: starts as soon as the run starts, alongside quantification."""
    input:
        vcf=config["genotype_vcf"],
        tbi=config["genotype_vcf"] + ".tbi",
        samples="results/genotype/samples.txt",
    output:
        temp(multiext("results/genotype/plink/chr{chrom}", ".pgen", ".pvar", ".psam")),
    params:
        plink_mem=lambda w, resources: int(resources.mem_mb * 0.8),  # MB; rest = plink2 overhead
        prefix="results/genotype/plink/chr{chrom}",
        filters=config["genotype"]["vcf_filters"],
        qc=config["genotype"]["plink_qc"],
    priority: 50
    threads: 16
    resources:
        mem_mb=64000,
        runtime=1440,
    log:
        "logs/genotype/plink_chr{chrom}.log",
    conda:
        "../envs/tools.yaml"
    shell:
        # --output-chr MT: Ensembl-style names (1..22) to match the GTF/phenotype BEDs
        "plink2 --vcf {input.vcf} --vcf-half-call missing --keep {input.samples} {params.filters} "
        "--output-chr MT --set-all-var-ids '@:#_$r/$a' --rm-dup force-first {params.qc} "
        "--threads {threads} --memory {params.plink_mem} "
        "--make-pgen --out {params.prefix} > {log} 2>&1"


rule merge_genotypes:
    priority: 50
    input:
        expand("results/genotype/plink/chr{chrom}.{ext}", chrom=CHROMS, ext=["pgen", "pvar", "psam"]),
    output:
        multiext(GENO_PREFIX, ".pgen", ".pvar", ".psam"),
    params:
        plink_mem=lambda w, resources: int(resources.mem_mb * 0.8),  # MB; rest = plink2 overhead
        prefixes=expand("results/genotype/plink/chr{chrom}", chrom=CHROMS),
        merge_list="results/genotype/plink/merge_list.txt",
        prefix=GENO_PREFIX,
    threads: 8
    resources:
        mem_mb=64000,
        runtime=240,
    log:
        "logs/genotype/merge_genotypes.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        printf '%s\\n' {params.prefixes} > {params.merge_list}
        plink2 --pmerge-list {params.merge_list} pfile --output-chr MT \
            --threads {threads} --memory {params.plink_mem} \
            --make-pgen --out {params.prefix} > {log} 2>&1
        """


rule genotype_pca:
    priority: 50
    input:
        multiext(GENO_PREFIX, ".pgen", ".pvar", ".psam"),
    output:
        eigenvec="results/genotype/pca/genotype_pcs.eigenvec",
        eigenval="results/genotype/pca/genotype_pcs.eigenval",
    params:
        plink_mem=lambda w, resources: int(resources.mem_mb * 0.8),  # MB; rest = plink2 overhead
        geno=GENO_PREFIX,
        prune=config["genotype"]["ld_prune"],
        n_pcs=config["genotype"]["n_pcs"],
        pruned="results/genotype/pca/ld_pruned",
        out="results/genotype/pca/genotype_pcs",
    threads: 8
    resources:
        mem_mb=32000,
        runtime=120,
    log:
        "logs/genotype/genotype_pca.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        plink2 --pfile {params.geno} --indep-pairwise {params.prune} \
            --threads {threads} --memory {params.plink_mem} --out {params.pruned} > {log} 2>&1
        plink2 --pfile {params.geno} --extract {params.pruned}.prune.in --pca {params.n_pcs} \
            --threads {threads} --memory {params.plink_mem} --out {params.out} >> {log} 2>&1
        """
