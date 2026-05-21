# circrna_pipeline v2.0

A Snakemake pipeline for circular RNA detection, consensus filtering, and RPM normalisation.

## Detection tools

| Tool | Aligner | Approach |
|------|---------|----------|
| **CIRCexplorer2** | STAR (chimeric) | Annotation-guided, annotation-first |
| **CIRI2** | BWA-MEM | Annotation-independent, SA-tag based |

Only circRNAs detected by **both** tools (2-of-2 consensus) are kept — this keeps precision high while the orthogonal aligners improve sensitivity vs. either tool alone.  
See: Vromman et al., *Nat Methods* 2023 (large-scale benchmark, 16 tools).

### Why not DCC?
DCC was dropped because it requires a full additional STAR run per mate (3× STAR jobs total per sample) and relies on junction files from all three runs — high compute cost for marginal sensitivity gain over the CIRCexplorer2 + CIRI2 combination.

### Why not find_circ?
The original find_circ is Python 2 / Bowtie2-based and unmaintained. The CIRI2 BWA-MEM leg covers the same "aligner-agnostic BSJ split read" niche with an actively maintained tool.

## Normalisation

```
RPM = (BSJ_reads / uniquely_mapped_reads) × 1_000_000
```

Both CIRCexplorer2 and CIRI2 BSJ counts are normalised independently. The final matrix uses the **mean RPM** of the two tools per circRNA per sample. A raw BSJ count matrix (`all_samples_raw_BSJ.tsv`) is also produced for DESeq2 / edgeR input.

## Quick start

### 1. Install Snakemake

```bash
conda create -n snakemake -c conda-forge -c bioconda snakemake>=8 pandas
conda activate snakemake
```

### 2. Edit config.yaml

Fill in all paths (genome FASTA, GTF, STAR index, BWA index, genePred annotation).

### 3. Create samplesheet.csv

```csv
sample,fastq_r1,fastq_r2
ctrl1,/data/ctrl1_R1.fastq.gz,/data/ctrl1_R2.fastq.gz
treat1,/data/treat1_R1.fastq.gz,/data/treat1_R2.fastq.gz
```

### 4. Prepare reference indices

```bash
# STAR index (once per genome/annotation combination)
STAR --runMode genomeGenerate \
     --genomeDir /path/to/star_index \
     --genomeFastaFiles genome.fa \
     --sjdbGTFfile annotation.gtf \
     --sjdbOverhang 149 \
     --runThreadN 16

# BWA index
bwa index genome.fa

# genePred for CIRCexplorer2 annotate
gtfToGenePred -genePredExt annotation.gtf /dev/stdout | \
    awk 'BEGIN{OFS="\t"}{print $12,$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' \
    > annotation.genePred.txt
```

### 5. Run

```bash
# Dry run
snakemake -n --use-conda

# Local (e.g. workstation)
snakemake --use-conda --cores 32

# SLURM cluster
snakemake --use-conda \
          --cluster "sbatch --time {resources.time_hrs}:00:00 --mem {resources.mem_gb}G --cpus-per-task {resources.threads}" \
          --jobs 50 \
          --latency-wait 60
```

## Output structure

```
results/
├── qc/                          FastQC reports
├── trimmed/                     Adapter-trimmed reads
├── star/<sample>/               STAR BAM + Chimeric.out.junction
├── bwa/<sample>/                BWA-MEM BAM + raw SAM (temp)
├── circexplorer2/<sample>/      CIRCexplorer2 annotated output
├── ciri2/<sample>/              CIRI2 output
├── consensus/<sample>/          Per-sample 2-of-2 consensus + RPM
├── final/
│   ├── all_samples_normalized.tsv   mean-RPM wide matrix
│   └── all_samples_raw_BSJ.tsv      mean-BSJ count wide matrix
└── multiqc_report.html
```

## Dependencies (conda environments in envs/)

| env file | tools |
|----------|-------|
| `qc.yaml` | fastqc, trim-galore, multiqc |
| `star.yaml` | STAR 2.7.11, samtools 1.20 |
| `circexplorer2.yaml` | circexplorer2 2.3 |
| `bwa.yaml` | bwa 0.7.18, samtools 1.20 |
| `ciri2.yaml` | ciri2 2.0, perl 5.32 |
| `python.yaml` | python 3.11, pandas 2.2, numpy 1.26 |
