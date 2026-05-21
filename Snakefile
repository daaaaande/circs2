"""
circrna_pipeline  –  v2.0
Snakemake circRNA detection & quantification pipeline

Detection tools
  - CIRCexplorer2  (STAR chimeric alignment)
  - CIRI2          (BWA-MEM alignment)

Only circRNAs detected by BOTH tools (2-of-2 consensus) are retained,
then normalised to RPM (back-splice junction reads / mapped reads × 10⁶).
"""

import pandas as pd
from pathlib import Path

# ── configuration ────────────────────────────────────────────────────────────
configfile: "config.yaml"

samples_df  = pd.read_csv(config["samplesheet"])
SAMPLES     = samples_df["sample"].tolist()
OUTPUT_DIR  = config["output_dir"]

# ── helper: final outputs ─────────────────────────────────────────────────────
def get_all_outputs():
    outputs = []
   # outputs += expand(
   #     "{output_dir}/qc/{sample}_R1_fastqc.html",
   #     output_dir=OUTPUT_DIR, sample=SAMPLES
   # )
    outputs += expand(
        "{output_dir}/trimmed/{sample}_R1_trimmed.fq.gz",
        output_dir=OUTPUT_DIR, sample=SAMPLES
    )
    outputs += expand(
        "{output_dir}/star/{sample}/Aligned.sortedByCoord.out.bam",
        output_dir=OUTPUT_DIR, sample=SAMPLES
    )
    outputs += expand(
        "{output_dir}/bwa/{sample}/{sample}.sorted.bam",
        output_dir=OUTPUT_DIR, sample=SAMPLES
    )
    outputs += expand(
        "{output_dir}/circexplorer2/{sample}/{sample}_circexplorer2.txt",
        output_dir=OUTPUT_DIR, sample=SAMPLES
    )
    # "{output_dir}/ciri_full/{sample}/{sample}.ciri_full.tsv",
    outputs += expand(
        "{output_dir}/ciri_full/{sample}/{sample}.ciri_full.tsv",
        output_dir=OUTPUT_DIR, sample=SAMPLES
    )
    outputs += expand(
        "{output_dir}/consensus/{sample}/{sample}_consensus.tsv",
        output_dir=OUTPUT_DIR, sample=SAMPLES
    )
    outputs += [f"{OUTPUT_DIR}/final/all_samples_normalized.tsv"]
    outputs += [f"{OUTPUT_DIR}/multiqc_report.html"]
    return outputs


rule all:
    input:
        get_all_outputs()


# ── QC ────────────────────────────────────────────────────────────────────────
rule fastqc:
    input:
        r1 = lambda wc: samples_df.loc[samples_df["sample"] == wc.sample, "fastq_r1"].values[0],
        r2 = lambda wc: samples_df.loc[samples_df["sample"] == wc.sample, "fastq_r2"].values[0]
    output:
        html_r1 = "{output_dir}/qc/{sample}_R1_fastqc.html",
        html_r2 = "{output_dir}/qc/{sample}_R2_fastqc.html",
        zip_r1  = "{output_dir}/qc/{sample}_R1_fastqc.zip",
        zip_r2  = "{output_dir}/qc/{sample}_R2_fastqc.zip"
    params:
        outdir = "{output_dir}/qc"
    conda:
        "envs/qc.yaml"
    log:
        "{output_dir}/logs/fastqc/{sample}.log"
    resources:
        threads  = 2,
        mem_gb   = lambda wildcards, attempt: 4 + (attempt * 2),
        time_hrs = lambda wildcards, attempt: attempt * 1
    message:
        "FastQC on {wildcards.sample} ..."
    shell:
        """
        fastqc -t {resources.threads} -o {params.outdir} \
            {input.r1} {input.r2} >{log} 2>&1
        # rename to canonical names expected downstream        mv {params.outdir}/$(basename {input.r1} .fastq.gz)_fastqc.html {output.html_r1}
        mv {params.outdir}/$(basename {input.r1} .fastq.gz)_fastqc.zip  {output.zip_r1}
        mv {params.outdir}/$(basename {input.r2} .fastq.gz)_fastqc.html {output.html_r2}
        mv {params.outdir}/$(basename {input.r2} .fastq.gz)_fastqc.zip  {output.zip_r2}
        """


# ── adapter trimming ──────────────────────────────────────────────────────────
rule trim_galore:
    input:
        r1 = lambda wc: samples_df.loc[samples_df["sample"] == wc.sample, "fastq_r1"].values[0],
        r2 = lambda wc: samples_df.loc[samples_df["sample"] == wc.sample, "fastq_r2"].values[0]
    output:
        r1      = "{output_dir}/trimmed/{sample}_R1_trimmed.fq.gz",
        r2      = "{output_dir}/trimmed/{sample}_R2_trimmed.fq.gz"
    params:
        outdir = "{output_dir}/trimmed",
        extra  = config.get("trim_galore_extra"),
        long_r1="{output_dir}/trimmed/{sample}_R1_001_val_1.fq.gz",
        long_r2="{output_dir}/trimmed/{sample}_R2_001_val_2.fq.gz"
    conda:
        "envs/qc.yaml"
    log:
        "{output_dir}/logs/trim_galore/{sample}.log"
    resources:
        threads  = lambda wildcards, attempt: attempt * 8,
        mem_gb   = lambda wildcards, attempt: 12 + (attempt * 4),
        time_hrs = lambda wildcards, attempt: attempt * 2
    message:
        "Trimming adapters for {wildcards.sample} ..."
    shell:
        """
        trim_galore {params.extra} \
            --cores {resources.threads} \
            -o {params.outdir} \
            --paired \
            {input.r1} {input.r2} >{log} --fastqc 2>&1
            # renaming the files for the run
        mv {params.long_r1} {output.r1} >{log}  2>&1
        mv {params.long_r2} {output.r2} >{log}  2>&1
        """


# ═══════════════════════════════════════════════════════════════════════════════
#  CIRCexplorer2 branch  (STAR chimeric alignment → CIRCexplorer2 parse/annotate)
# ═══════════════════════════════════════════════════════════════════════════════

rule star_align_chimeric:
    """
    STAR alignment with chimeric reads enabled for CIRCexplorer2.
    --chimSegmentMin 10 is the recommended minimum for back-splice detection.
    """
    input:
        r1    = "{output_dir}/trimmed/{sample}_R1_trimmed.fq.gz",
        r2    = "{output_dir}/trimmed/{sample}_R2_trimmed.fq.gz",
    output:
        bam       = "{output_dir}/star/{sample}/Aligned.sortedByCoord.out.bam",
        bai       = "{output_dir}/star/{sample}/Aligned.sortedByCoord.out.bam.bai",
        chimeric  = "{output_dir}/star/{sample}/Chimeric.out.junction",
        log_final = "{output_dir}/star/{sample}/Log.final.out"
    params:
        prefix      = "{output_dir}/star/{sample}/",
        gtf         = config["gtf"],
        extra       = config.get("star_extra", ""),
        index = config["star_index"]

    conda:
        "envs/star.yaml"
    log:
        "{output_dir}/logs/star/{sample}.log"
    resources:
        threads  = lambda wildcards, attempt: attempt * 8,
        mem_gb   = lambda wildcards, attempt: 32 + (attempt * 16),
        time_hrs = lambda wildcards, attempt: attempt * 4
    message:
        "STAR chimeric alignment for {wildcards.sample} ..."
    shell:
        """
        STAR \
            --runThreadN {resources.threads} \
            --genomeDir {params.index} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD \
            --outFileNamePrefix {params.prefix} \
            --sjdbGTFfile {params.gtf} \
            --chimSegmentMin 10 \
            --chimScoreMin 1 \
            --chimJunctionOverhangMin 10 \
            --chimOutType Junctions SeparateSAMold \
            --alignSJDBoverhangMin 1 \
            --alignSJoverhangMin 5 \
            --outFilterMismatchNmax 2 \
            --outFilterMultimapNmax 20 \
            --runRNGseed 1234 \
            {params.extra} >{log} 2>&1
        samtools index {output.bam} >>{log} 2>&1
        """


rule circexplorer2_parse:
    """Convert STAR chimeric junctions to the CIRCexplorer2 BED format."""
    input:
        chimeric = "{output_dir}/star/{sample}/Chimeric.out.junction"
    output:
        bed = "{output_dir}/circexplorer2/{sample}/back_spliced_junction.bed"
    conda:
        "envs/circexplorer2.yaml"
    log:
        "{output_dir}/logs/circexplorer2/{sample}_parse.log"
    resources:
        threads  = 1,
        mem_gb   = lambda wildcards, attempt: 4 + attempt,
        time_hrs = lambda wildcards, attempt: attempt * 1
    message:
        "CIRCexplorer2 parse for {wildcards.sample} ..."
    shell:
        """
        CIRCexplorer2 parse \
            -t STAR \
            -b {output.bed} \
            {input.chimeric} >{log} 2>&1
        """


rule circexplorer2_annotate:
    """Annotate back-splice junctions with gene model information."""
    input:
        bed       = "{output_dir}/circexplorer2/{sample}/back_spliced_junction.bed",
        ref       = config["gene_pred"],        # UCSC genePred / refFlat format
        genome    = config["genome_fasta"]
    output:
        circ_txt  = "{output_dir}/circexplorer2/{sample}/{sample}_circexplorer2.txt"
    conda:
        "envs/circexplorer2.yaml"
    log:
        "{output_dir}/logs/circexplorer2/{sample}_annotate.log"
    resources:
        threads  = 1,
        mem_gb   = lambda wildcards, attempt: 8 + (attempt * 4),
        time_hrs = lambda wildcards, attempt: attempt * 2
    message:
        "CIRCexplorer2 annotate for {wildcards.sample} ..."
    shell:
        """
        CIRCexplorer2 annotate \
            -r {input.ref} \
            -g {input.genome} \
            -b {input.bed} \
            -o {output.circ_txt} >{log} 2>&1
        """


# ═══════════════════════════════════════════════════════════════════════════════
#  CIRI2 branch  (BWA-MEM → CIRI2)
# ═══════════════════════════════════════════════════════════════════════════════

rule bwa_align:
    """
    BWA-MEM alignment required by CIRI2.
    CIRI2 reads the SAM directly (it needs the SA tag for chimeric detection),
    so we keep the unsorted SAM and also produce a sorted BAM for QC / reuse.
    """
    input:
        r1    = "{output_dir}/trimmed/{sample}_R1_trimmed.fq.gz",
        r2    = "{output_dir}/trimmed/{sample}_R2_trimmed.fq.gz",
    output:
        sam = temp("{output_dir}/bwa/{sample}/{sample}.sam"),
        bam = "{output_dir}/bwa/{sample}/{sample}.sorted.bam",
        bai = "{output_dir}/bwa/{sample}/{sample}.sorted.bam.bai"
    params:
        rg = r"@RG\tID:{sample}\tSM:{sample}\tPL:ILLUMINA",
        index = config["bwa_index"]         # prefix (genome.fa → genome.fa.bwt etc.)
    conda:
        "envs/bwa.yaml"
    log:
        "{output_dir}/logs/bwa/{sample}.log"
    resources:
        threads  = lambda wildcards, attempt: attempt * 8,
        mem_gb   = lambda wildcards, attempt: 26 + (attempt * 8),
        time_hrs = lambda wildcards, attempt: attempt * 3
    message:
        "BWA-MEM alignment for {wildcards.sample} ..."
    shell:
        """
        bwa mem \
            -t {resources.threads} \
            -T 19 \
            -R '{params.rg}' \
            {params.index} {input.r1} {input.r2} \
            > {output.sam} 2>{log}

        samtools sort -@ {resources.threads} -o {output.bam} {output.sam} >>{log} 2>&1
        samtools index {output.bam} >>{log} 2>&1
        """

rule ciri_full:
    input:
        sam    = "{output_dir}/bwa/{sample}/{sample}.sam",
        genome = config["genome_fasta"],
        gtf    = config["gtf"]
    output:
        tsv    = "{output_dir}/ciri_full/{sample}/{sample}.ciri_full.tsv"
    params:
        outdir    = "{output_dir}/ciri_full/{sample}",
        min_reads = config.get("ciri_full_min_reads", 2)
    conda:
        "envs/ciri_full.yaml"
    log:
        "{output_dir}/logs/ciri_full/{sample}.log"
    resources:
        threads  = lambda wildcards, attempt: attempt * 4,
        mem_gb   = lambda wildcards, attempt: 16 + (attempt * 8),
        time_hrs = lambda wildcards, attempt: attempt * 3
    message:
        "CIRI-full circRNA detection for {wildcards.sample} ..."
    shell:
        """
        CIRI_Full.py pipeline \
            -1 {input.sam} \
            -r {input.genome} \
            -a {input.gtf} \
            -o {params.outdir}/{wildcards.sample} \
            -t {resources.threads} \
            --min_junction_reads {params.min_reads} \
            >{log} 2>&1

        # CIRI-full writes <prefix>.ciri; symlink to expected name
        ln -sf {params.outdir}/{wildcards.sample}.ciri {output.tsv}
        """

# ═══════════════════════════════════════════════════════════════════════════════
#  Consensus & normalisation
# ═══════════════════════════════════════════════════════════════════════════════

rule consensus_filter:
    """
    Keep only circRNAs detected by BOTH CIRCexplorer2 AND CIRI2 (2-of-2 vote).
    Outputs a unified TSV with BSJ coordinates + read counts from each tool.
    Normalisation: RPM = (BSJ_reads / total_uniquely_mapped) × 1e6
    """
    input:
        cx2        = "{output_dir}/circexplorer2/{sample}/{sample}_circexplorer2.txt",
        ciri2    = "{output_dir}/ciri_full/{sample}/{sample}.ciri_full.tsv",
        star_log   = "{output_dir}/star/{sample}/Log.final.out"
    output:
        consensus  = "{output_dir}/consensus/{sample}/{sample}_consensus.tsv"
    conda:
        "envs/python.yaml"
    log:
        "{output_dir}/logs/consensus/{sample}.log"
    resources:
        threads  = 1,
        mem_gb   = lambda wildcards, attempt: 4 + attempt,
        time_hrs = lambda wildcards, attempt: attempt * 1
    message:
        "2-of-2 consensus filter & RPM normalisation for {wildcards.sample} ..."
    script:
        "scripts/consensus_filter.py"


rule merge_samples:
    """
    Merge all per-sample consensus TSVs into one wide matrix
    (rows = circRNA, columns = samples, values = RPM).
    """
    input:
        expand(
            "{output_dir}/consensus/{sample}/{sample}_consensus.tsv",
            output_dir=OUTPUT_DIR,
            sample=SAMPLES
        )
    output:
        matrix = "{output_dir}/final/all_samples_normalized.tsv"
    conda:
        "envs/python.yaml"
    log:
        "{output_dir}/logs/merge_samples.log"
    resources:
        threads  = 1,
        mem_gb   = 8,
        time_hrs = 1
    message:
        "Merging all samples into final normalised matrix ..."
    script:
        "scripts/merge_samples.py"


# ── MultiQC ───────────────────────────────────────────────────────────────────
rule multiqc:
    input:
     #   expand("{output_dir}/qc/{sample}_R1_fastqc.zip",   output_dir=OUTPUT_DIR, sample=SAMPLES),
      #  expand("{output_dir}/qc/{sample}_R2_fastqc.zip",   output_dir=OUTPUT_DIR, sample=SAMPLES),
        expand("{output_dir}/star/{sample}/Log.final.out", output_dir=OUTPUT_DIR, sample=SAMPLES)
    output:
        report = "{output_dir}/multiqc_report.html"
    params:
        outdir   = OUTPUT_DIR,
        indir    = OUTPUT_DIR
    conda:
        "envs/qc.yaml"
    log:
        "{output_dir}/logs/multiqc.log"
    resources:
        threads  = 1,
        mem_gb   = 4,
        time_hrs = 1
    message:
        "MultiQC report ..."
    shell:
        """
        multiqc {params.indir} \
            --outdir {params.outdir} \
            --filename multiqc_report.html \
            --force >{log} 2>&1
        """
