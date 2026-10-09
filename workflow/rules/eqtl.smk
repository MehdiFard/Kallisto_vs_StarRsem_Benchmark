## Phenotypes -> covariates -> tensorQTL -> eQTL summary.
## Both tools go through identical code, so differences come from quantification only.
## tensorQTL runs on the GPU partition (profiles/slurm/config.yaml) from the conda env envs/tensorqtl.yaml.
## Isoforms are mapped grouped by gene,
## so an isoform "eGene" is a gene with at least one significant isoform.
import os

GENO_FILES = multiext(GENO_PREFIX, ".pgen", ".pvar", ".psam")
TQTL = config["tensorqtl"]
# relative to the working directory, so moving the project does not change the rule params
TQTL_SCRIPT = os.path.relpath(os.path.join(workflow.basedir, "scripts", "run_tensorqtl.py"))
# Optional, not set in config.yaml: tensorqtl.container (+ apptainer, bind, container_setup) runs tensorQTL inside
# an apptainer image with --nv instead of the conda env; "$@" passes the script and its arguments through to python.
if TQTL.get("container"):
    _bind = f" --bind {TQTL['bind']}" if TQTL.get("bind") else ""
    _setup = f"{TQTL['container_setup']} && " if TQTL.get("container_setup") else ""
    TQTL_GPU_RUN = (f"{TQTL.get('apptainer', 'apptainer')} exec --nv{_bind} {TQTL['container']} "
                    f"bash -c '{_setup}python \"$@\"' _")
    TQTL_ENV = None
else:
    TQTL_GPU_RUN = "python"
    TQTL_ENV = "../envs/tensorqtl.yaml"


# ------------------------------------------------------------------------------------------------ phenotypes + covariates


rule prepare_phenotypes:
    """Filter -> TMM CPM -> inverse normal transform -> BED sorted by TSS (+ isoform->gene groups)."""
    input:
        counts="results/expression/{tool}/{level}.counts.tsv.gz",
        tpm="results/expression/{tool}/{level}.tpm.tsv.gz",
        genes="resources/annotation/genes.tsv",
        transcripts="resources/annotation/transcripts.tsv",
    output:
        bed=temp("results/phenotypes/{tool}/{level}.bed"),
        groups="results/phenotypes/{tool}/{level}.groups.tsv",
    params:
        chromosomes=CHROMS,
        tpm_threshold=config["expression"]["tpm_threshold"],
        count_threshold=config["expression"]["count_threshold"],
        sample_frac_threshold=config["expression"]["sample_frac_threshold"],
    resources:
        mem_mb=32000,
        runtime=120,
    log:
        "logs/expression/prepare_{tool}_{level}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/prepare_phenotypes.py"


rule index_phenotype_bed:
    input:
        "results/phenotypes/{tool}/{level}.bed",
    output:
        bed="results/phenotypes/{tool}/{level}.bed.gz",
        tbi="results/phenotypes/{tool}/{level}.bed.gz.tbi",
    conda:
        "../envs/tools.yaml"
    shell:
        "bgzip -c {input} > {output.bed} && tabix -p bed {output.bed}"


rule covariates:
    """Genotype PCs + expression PCs (computed separately for each tool/level)."""
    input:
        bed="results/phenotypes/{tool}/{level}.bed.gz",
        eigenvec="results/genotype/pca/genotype_pcs.eigenvec",
    output:
        "results/covariates/{tool}/{level}.covariates.tsv",
    params:
        n_genotype_pcs=config["covariates"]["n_genotype_pcs"],
        n_expression_pcs=config["covariates"]["n_expression_pcs"],
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        "logs/expression/covariates_{tool}_{level}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/make_covariates.py"

# ------------------------------------------------------------------------------------------------ tensorQTL


def tensorqtl_common(wildcards, input):
    args = (f"--genotypes {GENO_PREFIX} --phenotypes {input.bed} --covariates {input.covariates} "
            f"--window {TQTL['window']} --maf-threshold {TQTL['maf_threshold']}")
    if wildcards.level == "isoform":
        args += f" --groups {input.groups}"
    return args


rule tensorqtl_cis:
    """Permutation pass: top variant per phenotype (or per gene for isoforms) + q-values."""
    input:
        GENO_FILES,
        bed="results/phenotypes/{tool}/{level}.bed.gz",
        groups="results/phenotypes/{tool}/{level}.groups.tsv",
        covariates="results/covariates/{tool}/{level}.covariates.tsv",
    output:
        "results/eqtl/{tool}/{level}/cis.tsv.gz",
    params:
        common=tensorqtl_common,
        gpu_run=TQTL_GPU_RUN,
        script=TQTL_SCRIPT,
        extra=(f"--nperm {TQTL['nperm']} --seed {TQTL['seed']} --fdr {TQTL['fdr']} "
               f"--qvalue-lambda {TQTL['qvalue_lambda']}"),
    threads: 8
    resources:
        mem_mb=96000,
        runtime=1440,
    log:
        "logs/eqtl/cis_{tool}_{level}.log",
    conda:
        TQTL_ENV
    shell:
        "{params.gpu_run} {params.script} cis {params.common} --output {output} {params.extra} > {log} 2>&1"


rule tensorqtl_cis_independent:
    """Conditionally independent cis-eQTLs (forward-backward stepwise regression)."""
    input:
        GENO_FILES,
        bed="results/phenotypes/{tool}/{level}.bed.gz",
        groups="results/phenotypes/{tool}/{level}.groups.tsv",
        covariates="results/covariates/{tool}/{level}.covariates.tsv",
        cis="results/eqtl/{tool}/{level}/cis.tsv.gz",
    output:
        "results/eqtl/{tool}/{level}/cis_independent.tsv.gz",
    params:
        common=tensorqtl_common,
        gpu_run=TQTL_GPU_RUN,
        script=TQTL_SCRIPT,
        extra=f"--nperm {TQTL['nperm']} --seed {TQTL['seed']} --fdr {TQTL['fdr']}",
    threads: 8
    resources:
        mem_mb=96000,
        runtime=1440,
    log:
        "logs/eqtl/cis_independent_{tool}_{level}.log",
    conda:
        TQTL_ENV
    shell:
        "{params.gpu_run} {params.script} cis_independent {params.common} --cis-results {input.cis} "
        "--output {output} {params.extra} > {log} 2>&1"


rule tensorqtl_cis_nominal:
    """All cis variant-phenotype pairs (parquet per chromosome) + the significant subset."""
    input:
        GENO_FILES,
        bed="results/phenotypes/{tool}/{level}.bed.gz",
        groups="results/phenotypes/{tool}/{level}.groups.tsv",
        covariates="results/covariates/{tool}/{level}.covariates.tsv",
        cis="results/eqtl/{tool}/{level}/cis.tsv.gz",
    output:
        nominal=directory("results/eqtl/{tool}/{level}/nominal"),
        signif="results/eqtl/{tool}/{level}/cis_signif_pairs.tsv.gz",
    params:
        common=tensorqtl_common,
        gpu_run=TQTL_GPU_RUN,
        script=TQTL_SCRIPT,
        extra=f"--fdr {TQTL['fdr']}",
    threads: 8
    resources:
        mem_mb=96000,
        runtime=1440,
    log:
        "logs/eqtl/cis_nominal_{tool}_{level}.log",
    conda:
        TQTL_ENV
    shell:
        "{params.gpu_run} {params.script} cis_nominal {params.common} --cis-results {input.cis} "
        "--nominal-dir {output.nominal} --prefix {wildcards.tool}_{wildcards.level} "
        "--output {output.signif} {params.extra} > {log} 2>&1"


rule tensorqtl_trans:
    """Optional (target `trans`): trans associations with p < threshold, cis (+/-5 Mb) removed."""
    input:
        GENO_FILES,
        bed="results/phenotypes/{tool}/{level}.bed.gz",
        groups="results/phenotypes/{tool}/{level}.groups.tsv",
        covariates="results/covariates/{tool}/{level}.covariates.tsv",
    output:
        "results/eqtl/{tool}/{level}/trans.tsv.gz",
    params:
        common=tensorqtl_common,
        gpu_run=TQTL_GPU_RUN,
        script=TQTL_SCRIPT,
        extra=(f"--pval-threshold {TQTL['trans_pval_threshold']} "
               f"--trans-cis-window {TQTL['trans_cis_window']} --batch-size {TQTL['trans_batch_size']}"),
    threads: 8
    resources:
        mem_mb=256000,
        runtime=1440,
    log:
        "logs/eqtl/trans_{tool}_{level}.log",
    conda:
        TQTL_ENV
    shell:
        "{params.gpu_run} {params.script} trans {params.common} --output {output} {params.extra} > {log} 2>&1"

# ------------------------------------------------------------------------------------------------ summary


rule eqtl_summary:
    input:
        cis=expand("results/eqtl/{tool}/{level}/cis.tsv.gz", tool=TOOLS, level=LEVELS),
        independent=expand("results/eqtl/{tool}/{level}/cis_independent.tsv.gz", tool=TOOLS, level=LEVELS),
        signif=expand("results/eqtl/{tool}/{level}/cis_signif_pairs.tsv.gz", tool=TOOLS, level=LEVELS),
    output:
        summary="results/summary/eqtl_summary.tsv",
        overlap="results/summary/eqtl_egene_overlap.tsv",
    params:
        tools=TOOLS,
        levels=LEVELS,
        fdr=config["tensorqtl"]["fdr"],
    resources:
        mem_mb=16000,
        runtime=30,
    log:
        "logs/summary/eqtl_summary.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_eqtl.py"
