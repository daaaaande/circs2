
import pandas as pd
from pathlib import Path

# ── configuration ────────────────────────────────────────────────────────────
configfile: "config.yaml"

samples_df  = pd.read_csv(config["samplesheet"])
SAMPLES     = samples_df["sample"].tolist()
OUTPUT_DIR  = config["output_dir"]
# adding tools to later enable a 1-4/4 to 2/2 vote with any of the 4 tools, maybe dropping it to 3 max if one tool fails to be usable 
# for now not including it since the extra complexity does not seem to be worth the extra output quality.

rule circminer_detect:
    """
    CircMiner: Sensitive and fast detection of circular RNAs using splice-aware pseudo-alignment.
    """
    input:
        r1  = "{output_dir}/trimmed/{sample}_R1_trimmed.fq.gz",
        r2  = "{output_dir}/trimmed/{sample}_R2_trimmed.fq.gz",
    output:
        circ = "{output_dir}/circminer/{sample}/circrna_results.txt" # Adjust to match exact tool output name
    params:
        outdir = "{output_dir}/circminer/{sample}",
        ref_fa = config["ref_fasta"],
        gtf    = config["gtf"],
        extra  = config.get("circminer_extra", "")
    conda:
        "envs/circminer.yaml"
    log:
        "{output_dir}/logs/circminer/{sample}.log"
    resources:
        threads  = lambda wildcards, attempt: attempt * 8,
        mem_gb   = lambda wildcards, attempt: 16 + (attempt * 8),
        time_hrs = lambda wildcards, attempt: attempt * 2
    message:
        "Running CircMiner detection for {wildcards.sample} ..."
    shell:
        """
        circminer \
            -t {resources.threads} \
            -r {params.ref_fa} \
            -g {params.gtf} \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            {params.extra} >{log} 2>&1
        """


rule psirc_detect:
    """
    psirc: Detection of back-splicing junctions and full-length circular/linear isoforms.
    Requires a pre-built psirc custom transcriptome FASTA index from the GTF.
    """
    input:
        r1  = "{output_dir}/trimmed/{sample}_R1_trimmed.fq.gz",
        r2  = "{output_dir}/trimmed/{sample}_R2_trimmed.fq.gz",
        transcriptome_fa = config["psirc_transcriptome_fa"] # Built beforehand from GTF via psirc utilities
    output:
        out_dir = directory("{output_dir}/psirc/{sample}")
    params:
        kallisto_bin = config.get("psirc_kallisto", "kallisto"),
        extra        = config.get("psirc_extra", "")
    conda:
        "envs/psirc.yaml"
    log:
        "{output_dir}/logs/psirc/{sample}.log"
    resources:
        threads  = lambda wildcards, attempt: attempt * 8,
        mem_gb   = lambda wildcards, attempt: 24 + (attempt * 8),
        time_hrs = lambda wildcards, attempt: attempt * 3
    message:
        "Running psirc pipeline for {wildcards.sample} ..."
    shell:
        """
        mkdir -p {output.out_dir}
        perl psirc_v1.0.pl \
            -f \
            -o {output.out_dir} \
            {input.transcriptome_fa} \
            {params.kallisto_bin} \
            {input.r1} \
            {input.r2} \
            {params.extra} >{log} 2>&1
        """