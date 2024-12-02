SAMPLE = ['A','B','C']
rule target:
    input:
        'plots/quals.pdf'
rule bwa_map:
    input:
        Ref ='data/genome.fa',
        fastq ='data/samples/{sample}.fastq',
    output:
        'mapped_reads/{sample}.bam'
    shell:
        'bwa mem {input}| samtools view -Sb -> {output}'  
rule samtools_sort:
    input:
        'mapped_reads/{sample}.bam'
    output:
        'sorted_reads/{sample}.bam'
    shell:
        'samtools sort -T sorted_reads/{wildcards.sample} -O bam {input} > {output}'
rule samtools_index:
    input:
        'sorted_reads/{sample}.bam'
    output:
        'sorted_reads/{sample}.bam.bai'
    shell:
        'samtools index {input}'
rule bcftools_call:
    input:
        fa = 'data/genome.fa',
        bam = expand('sorted_reads/{sample}.bam', sample = SAMPLE),
        bai = expand('sorted_reads/{sample}.bam.bai', sample = SAMPLE),
    output:
        'calls/all.vcf'
    shell:
        'bcftools mpileup -f {input.fa} {input.bam} | '
        'bcftools call -mv -> {output}'   
rule plot_quals:
    input:
        'calls/all.vcf'
    output:
        'plots/quals.pdf'
    script:
        'Scripts/plot_quals.py'
