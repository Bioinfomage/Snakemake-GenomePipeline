# Snakemake-GenomePipeline
This repository contains a Snakemake pipeline for aligning reads to a reference genome, Snakemake-GenomePipeline is a scalable and reproducible workflow for genome analysis using the Snakemake workflow management system. This pipeline is designed to automate the key steps of processing sequencing data, including read alignment, BAM file sorting, indexing, and variant calling, from raw FASTQ files to high-quality VCF outputs.

## Introduction to Snakemake

Snakemake is a workflow management system that enables reproducible and scalable data analyses. It uses a Python-based language to define rules for executing tasks and their dependencies.
 
# Features
1.Automated Workflow: Streamlines the entire genome analysis process with minimal manual intervention.
2.Scalability: Handles multiple samples efficiently using Snakemake's dependency management.
3.Reproducibility: Ensures consistent results by tracking dependencies and intermediate outputs.
4.Customizable: Easily extendable to include additional steps, such as quality control or annotation.
5.Visualization: Generate directed acyclic graphs (DAG) to visualize the workflow.

**Prerequisites**

**Install the following tools:**
--BWA
--Samtools
--BCFtools
**Install Python and Snakemake**
--pip install snakemake

### Basic Concepts in Snakemake

1. **Rules**: Each rule represents a step in the pipeline.
2. **Wildcards**: Placeholders for variable parts of file names.
3. **Inputs and Outputs**: Define dependencies and results of a rule.
4. **Shell Commands**: The actual commands that perform the task.
5. **Target Rule**: The final outputs the pipeline aims to generate.

## Workflow Overview

### Input Data

- Reference genome: `data/genome.fa`
- Read files: `data/samples/{sample}.fastq` for samples A, B, and C.
  
## Steps in the Pipeline

1. **Read Alignment**: Align reads to the reference genome using BWA.
2. **Convert SAM to BAM**: Convert uncompressed SAM files to BAM files using `samtools view`.
3. **Sort BAM Files**: Sort BAM files using `samtools sort`.
4. **Index BAM Files**: Create BAM index files using `samtools index`.
5. **Variant Calling**: Generate a VCF file using `bcftools`.
6. **Plot Quality Metrics**: Produces visualizations of variant quality.

### Output Data
- Aligned BAM files (output of bwa_map) : 'mapped_reads/A.bam, B.bam, C.bam'
- Sorted BAM files (output of samtools_sort): `sorted_reads/A.bam, B.bam, C.bam`
- BAM index files (output of samtools_index): `sorted_reads/A.bam.bai, B.bam.bai, C.bam.bai`
- Variant call file (output of call_variants): `calls/all.vcf`
- Quality metrics plot (output of plot_quals): 'plots/quals.pdf'

## Commands to Execute the Workflow

**Option1**
### Dry-run (Check Workflow)
snakemake --cores 1 --dry-run
### Execute Workflow
snakemake --cores 1
### Generate Workflow Graph
snakemake --dag | dot -Tpdf > workflow_dag.pdf
**Option2**
### Run on Multiple Cores
snakemake --cores 4
**Option3**
### Clean Intermediate Files
snakemake --cores 1 --delete-temp-output

**Notes**
--Snakemake automatically creates output directories.
--Modify the SAMPLES list in the Snakefile to include your sample names.
--Custom scripts can be added to extend the workflow.

