# Snakefile: A Snakemake script for aligning reads, converting, sorting, and indexing BAM files.

# Define sample names
SAMPLES = ["A", "B", "C"]

# Rule to align reads using BWA and convert SAM to BAM
rule bwa_map:
    input:
        ref="data/genome.fa",
        query="data/samples/{sample}.fastq"
    output:
        bam="mapped_reads/{sample}.bam"
    shell:
        """
        bwa mem {input.ref} {input.query} | samtools view -Sb - > {output.bam}
        """

# Rule to sort BAM files
rule samtools_sort:
    input:
        bam="mapped_reads/{sample}.bam"
    output:
        sorted_bam="sorted_reads/{sample}.bam"
    shell:
        """
        samtools sort -T sorted_reads/{wildcards.sample} -O bam {input.bam} > {output.sorted_bam}
        """

# Rule to index BAM files
rule samtools_index:
    input:
        sorted_bam="sorted_reads/{sample}.bam"
    output:
        bai="sorted_reads/{sample}.bam.bai"
    shell:
        """
        samtools index {input.sorted_bam}
        """

# Rule to aggregate BAM files into VCF
rule call_variants:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        vcf="calls/all.vcf"
    shell:
        """
        bcftools mpileup -f {input.fa} {input.bam} | bcftools call -mv - > {output.vcf}
        """

# Example of adding custom scripts
rule plot_quals:
    input:
        vcf="calls/all.vcf"
    output:
        plot="plots/quals.pdf"
    script:
        "scripts/plot_quals.py"

# Target rule for final output
rule all:
    input:
        expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES),
        "calls/all.vcf",
        "plots/quals.pdf"
