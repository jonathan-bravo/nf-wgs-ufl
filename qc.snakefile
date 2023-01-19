fastqc_suffixes = ('1P_fastqc.html', '2P_fastqc.html', '1P_fastqc.zip', '2P_fastqc.zip')
qc_suffixes = ['trimmomatic.stats.log', 'md_metrics.txt', 'wgs_metrics.txt']
qc_suffixes.extend(fastqc_suffixes)


rule run_fastqc: #fastqc
    input:
        expand(outdir + "{sample}/{sample}.{suffix}.fastq.gz", suffix = ('1P', '2P'))
    output:
        expand(outdir + "{sample}/QC/{sample}.{suffix}", suffix = fastqc_suffixes)
    params:
        fastqc_out = outdir + "{sample}/Fastqc/"
    conda:
        "envs/qc.yaml"
    threads:
        2
    resources:
        mem_mb = 2048
    shell:
        "mkdir -p {params.fastqc_out}; "
        "fastqc -t {threads} -o {params.fastqc_out} -f fastq {input}"


rule collect_wgs_metrics: #picard
    input:
        expand(ref_genome + suffix, suffix = ref_suffixes),
        outdir + "{sample}/Alignment/{sample}.sorted.md.bam.bai",
        bam = outdir + "{sample}/Alignment/{sample}.sorted.md.bam",
        reference = ref_genome,
    output:
        outdir + "{sample}/QC/{sample}.wgs_metrics.txt"
    conda:
        "envs/qc.yaml"
    threads:
        4
    resources:
        mem_mb = 63488
    shell:
        "picard CollectWgsMetrics -I {input.bam} -O {output} -R {reference}"


rule run_multiqc:
    input:
        expand(outdir + "{sample}/QC/{sample}.{suffix}", suffix = qc_suffixes)
    output:
        expand(outdir + "{sample}/MultiQC/{sample}{suffix}", suffix = multiqc_suffixes)
    conda:
        "envs/qc.yaml"
    threads:
        2
    resources:
        mem_mb = 2048
    shell:
        "multiqc -n {wildcards.sample} ."