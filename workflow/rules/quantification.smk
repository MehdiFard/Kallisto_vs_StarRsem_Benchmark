## References -> (trimming -> STAR) -> RSEM, and kallisto -> expression matrices -> quantification benchmark tables.
## Both tools quantify the same transcript set: kallisto is indexed on the transcript FASTA written by
## rsem-prepare-reference (252,835 transcripts, Ensembl 110).
## STAR and kallisto read the raw FASTQs, or Trimmomatic output when trimming: true (config.yaml).
## Samples with two sequencing runs are aligned/quantified jointly; each run gets its own read group.


def sample_reads(sample, mate):
    """FASTQs given to STAR and kallisto: trimmed here (trimming: true) or raw."""
    pattern = "results/trimmed/{run}_R{mate}.fq.gz" if config["trimming"] else config["fastq"][f"r{mate}"]
    return [pattern.format(run=run, mate=mate) for run in RUNS[sample]]


def star_output(sample, suffix):
    return f"results/star/{sample}/{sample}.{suffix}"


# ------------------------------------------------------------------------------------------------ references


rule star_index:
    input:
        fasta=config["reference"]["genome_fasta"],
        gtf=config["reference"]["gtf"],
    output:
        directory(f"resources/star/{REF}"),
    params:
        overhang=config["reference"]["read_length"] - 1,
        tmp=f"resources/star/{REF}_tmp",
    threads: 16
    resources:
        mem_mb=64000,
        runtime=360,
    log:
        "logs/reference/star_index.log",
    conda:
        "../envs/tools.yaml"
    shell:
        """
        rm -rf {params.tmp} && mkdir -p {output}
        STAR --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {output} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.overhang} \
            --outTmpDir {params.tmp} \
            --outFileNamePrefix {output}/ > {log} 2>&1
        rm -rf {params.tmp}
        """


rule rsem_reference:
    input:
        fasta=config["reference"]["genome_fasta"],
        gtf=config["reference"]["gtf"],
    output:
        multiext(RSEM_REF, ".grp", ".ti", ".seq", ".chrlist", ".idx.fa", ".n2g.idx.fa", ".transcripts.fa"),
    params:
        prefix=RSEM_REF,
    threads: 8
    resources:
        mem_mb=32000,
        runtime=240,
    log:
        "logs/reference/rsem_reference.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "rsem-prepare-reference --gtf {input.gtf} --num-threads {threads} "
        "{input.fasta} {params.prefix} > {log} 2>&1"


rule kallisto_index:
    input:
        f"{RSEM_REF}.transcripts.fa",
    output:
        f"resources/kallisto/{REF}.idx",
    threads: 16
    resources:
        mem_mb=64000,
        runtime=240,
    log:
        "logs/reference/kallisto_index.log",
    conda:
        "../envs/tools.yaml"
    shell:
        "kallisto index --threads {threads} --index {output} {input} > {log} 2>&1"


rule annotation_tables:
    """Gene (with TSS) and transcript tables from the GTF, used to build phenotype BEDs."""
    input:
        gtf=config["reference"]["gtf"],
    output:
        genes="resources/annotation/genes.tsv",
        transcripts="resources/annotation/transcripts.tsv",
    resources:
        mem_mb=16000,
        runtime=60,
    log:
        "logs/reference/annotation_tables.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/gtf_to_tables.py"


# ------------------------------------------------------------------------------------------------ reads (trimming: true only)


rule trimmomatic:
    """Adapter / quality trimming, one job per sequencing run (trimming: true)."""
    input:
        r1=lambda w: config["fastq"]["r1"].format(run=w.run),
        r2=lambda w: config["fastq"]["r2"].format(run=w.run),
        adapters=config["trimmomatic"]["adapters"],
    output:
        r1=keep_or_temp("results/trimmed/{run}_R1.fq.gz", "keep_trimmed_fastq"),
        r2=keep_or_temp("results/trimmed/{run}_R2.fq.gz", "keep_trimmed_fastq"),
        u1=temp("results/trimmed/{run}_R1.unpaired.fq.gz"),
        u2=temp("results/trimmed/{run}_R2.unpaired.fq.gz"),
    params:
        steps=lambda w, input: config["trimmomatic"]["steps"].format(adapters=input.adapters),
    threads: 8
    resources:
        mem_mb=16000,
        runtime=480,
    log:
        "logs/trimmomatic/{run}.log",
    benchmark:
        "benchmarks/trimmomatic/{run}.tsv"
    conda:
        "../envs/tools.yaml"
    shell:
        "trimmomatic PE -threads {threads} -phred33 {input.r1} {input.r2} "
        "{output.r1} {output.u1} {output.r2} {output.u2} {params.steps} > {log} 2>&1"


# ------------------------------------------------------------------------------------------------ quantification


def star_read_groups(wildcards):
    # STAR syntax: one RG per comma-separated input file, entries separated by " , "
    return " , ".join(f"ID:{run} SM:{wildcards.sample} LB:{wildcards.sample} PL:ILLUMINA PU:{run}"
                      for run in RUNS[wildcards.sample])


STAR_GENOME_BAM = (
    {"genome_bam": "results/star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"}
    if config["keep_star_genome_bam"] else {}
)


rule star_align:
    """STAR 2-pass; writes the transcriptome BAM for RSEM (no genome BAM unless keep_star_genome_bam)."""
    input:
        r1=lambda w: sample_reads(w.sample, 1),
        r2=lambda w: sample_reads(w.sample, 2),
        index=f"resources/star/{REF}",
    output:
        tx_bam=keep_or_temp("results/star/{sample}/{sample}.Aligned.toTranscriptome.out.bam",
                            "keep_star_transcriptome_bam"),
        log_final="results/star/{sample}/{sample}.Log.final.out",
        gene_counts="results/star/{sample}/{sample}.ReadsPerGene.out.tab",
        sj="results/star/{sample}/{sample}.SJ.out.tab",
        **STAR_GENOME_BAM,
    params:
        r1=lambda w, input: ",".join(input.r1),
        r2=lambda w, input: ",".join(input.r2),
        rg=star_read_groups,
        prefix="results/star/{sample}/{sample}.",
        tmp="results/star/{sample}/_STARtmp",
        sam_type=("BAM SortedByCoordinate --outSAMunmapped Within --limitBAMsortRAM 30000000000"
                  if config["keep_star_genome_bam"] else "None"),
        extra=config["star"]["extra"],
    threads: 16
    resources:
        mem_mb=48000,
        runtime=480,
    log:
        "logs/star/{sample}.log",
    benchmark:
        "benchmarks/star/{sample}.tsv"
    conda:
        "../envs/tools.yaml"
    shell:
        """
        rm -rf {params.tmp}
        STAR --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {input.index} \
            --readFilesIn {params.r1} {params.r2} \
            --readFilesCommand zcat \
            --twopassMode Basic \
            --quantMode TranscriptomeSAM GeneCounts \
            --outSAMtype {params.sam_type} \
            --outSAMattrRGline {params.rg} \
            --outTmpDir {params.tmp} \
            --outFileNamePrefix {params.prefix} \
            {params.extra} > {log} 2>&1
        rm -rf {params.prefix}_STARgenome {params.prefix}_STARpass1 {params.tmp}
        """


rule rsem_quant:
    input:
        bam=lambda w: star_output(w.sample, "Aligned.toTranscriptome.out.bam"),
        ref=multiext(RSEM_REF, ".grp", ".ti", ".seq", ".idx.fa"),
    output:
        isoforms="results/rsem/{sample}/{sample}.isoforms.results.gz",
        genes="results/rsem/{sample}/{sample}.genes.results.gz",
        cnt="results/rsem/{sample}/{sample}.stat/{sample}.cnt",
    params:
        ref=RSEM_REF,
        prefix="results/rsem/{sample}/{sample}",
        strandedness=config["rsem"]["strandedness"],
        seed=config["rsem"]["seed"],
        extra=config["rsem"]["extra"],
    threads: 16
    resources:
        mem_mb=32000,
        runtime=720,
    log:
        "logs/rsem/{sample}.log",
    benchmark:
        "benchmarks/rsem/{sample}.tsv"
    conda:
        "../envs/tools.yaml"
    shell:
        """
        rsem-calculate-expression \
            --paired-end \
            --alignments \
            --no-bam-output \
            --num-threads {threads} \
            --strandedness {params.strandedness} \
            --seed {params.seed} \
            --temporary-folder {params.prefix}.temp \
            {params.extra} \
            {input.bam} {params.ref} {params.prefix} > {log} 2>&1
        gzip -f {params.prefix}.isoforms.results {params.prefix}.genes.results
        """


rule kallisto_quant:
    input:
        r1=lambda w: sample_reads(w.sample, 1),
        r2=lambda w: sample_reads(w.sample, 2),
        index=f"resources/kallisto/{REF}.idx",
    output:
        abundance="results/kallisto/{sample}/abundance.tsv",
        h5="results/kallisto/{sample}/abundance.h5",
        run_info="results/kallisto/{sample}/run_info.json",
    params:
        # kallisto takes mate pairs adjacently: run1_R1 run1_R2 run2_R1 run2_R2
        reads=lambda w, input: " ".join(f"{a} {b}" for a, b in zip(input.r1, input.r2)),
        outdir="results/kallisto/{sample}",
        seed=config["kallisto"]["seed"],
        extra=config["kallisto"]["extra"],
    threads: 8
    resources:
        mem_mb=16000,
        runtime=240,
    log:
        "logs/kallisto/{sample}.log",
    benchmark:
        "benchmarks/kallisto/{sample}.tsv"
    conda:
        "../envs/tools.yaml"
    shell:
        "kallisto quant --index {input.index} --output-dir {params.outdir} "
        "--threads {threads} --seed {params.seed} {params.extra} "
        "{params.reads} > {log} 2>&1"


# ------------------------------------------------------------------------------------------------ matrices + tables


def quant_files(wildcards):
    if wildcards.tool == "rsem":
        return [f"results/rsem/{s}/{s}.isoforms.results.gz" for s in SAMPLES]
    return [f"results/kallisto/{s}/abundance.tsv" for s in SAMPLES]


rule assemble_expression:
    """Isoform matrices from the tool output; gene = sum of its isoforms (same rule for both tools)."""
    input:
        quant=quant_files,
        transcripts="resources/annotation/transcripts.tsv",
    output:
        isoform_counts="results/expression/{tool}/isoform.counts.tsv.gz",
        isoform_tpm="results/expression/{tool}/isoform.tpm.tsv.gz",
        gene_counts="results/expression/{tool}/gene.counts.tsv.gz",
        gene_tpm="results/expression/{tool}/gene.tpm.tsv.gz",
    params:
        samples=SAMPLES,
    resources:
        mem_mb=32000,
        runtime=120,
    log:
        "logs/expression/assemble_{tool}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/assemble_expression.py"


rule mapping_qc:
    """STAR / RSEM / kallisto mapping stats and compute cost per sample (tables only, no plots)."""
    input:
        star=[star_output(s, "Log.final.out") for s in SAMPLES],
        rsem=expand("results/rsem/{s}/{s}.stat/{s}.cnt", s=SAMPLES),
        kallisto=expand("results/kallisto/{s}/run_info.json", s=SAMPLES),
    output:
        qc="results/summary/mapping_qc.tsv",
        compute="results/summary/compute_benchmark.tsv",
    params:
        samples=SAMPLES,
        runs=RUNS,
    resources:
        mem_mb=8000,
        runtime=30,
    log:
        "logs/summary/mapping_qc.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/summarize_mapping.py"


rule compare_quantification:
    """RSEM vs kallisto agreement per sample and per feature."""
    input:
        rsem_counts="results/expression/rsem/{level}.counts.tsv.gz",
        rsem_tpm="results/expression/rsem/{level}.tpm.tsv.gz",
        kallisto_counts="results/expression/kallisto/{level}.counts.tsv.gz",
        kallisto_tpm="results/expression/kallisto/{level}.tpm.tsv.gz",
        transcripts="resources/annotation/transcripts.tsv",
        genes="resources/annotation/genes.tsv",
    output:
        per_sample="results/summary/quant_correlation_{level}_per_sample.tsv",
        per_feature="results/summary/quant_correlation_{level}_per_feature.tsv.gz",
    resources:
        mem_mb=32000,
        runtime=120,
    log:
        "logs/summary/compare_quantification_{level}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/compare_quantification.py"
