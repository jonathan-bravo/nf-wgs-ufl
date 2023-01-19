rule run_fastqc: #fastqc
    input:
        outdir + "{sample}/Trimmomatic/Paired/{sample}.1P.fastq.gz",
        outdir + "{sample}/Trimmomatic/Paired/{sample}.2P.fastq.gz"
    output:
        outdir + "{sample}/Fastqc/{sample}.1P_fastqc.html",
        outdir + "{sample}/Fastqc/{sample}.2P_fastqc.html",
        outdir + "{sample}/Fastqc/{sample}.1P_fastqc.zip",
        outdir + "{sample}/Fastqc/{sample}.2P_fastqc.zip"
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
        ref_genome + ".fai",
        ref_genome + ".gzi",
        reference = ref_genome,
        bam = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai"
    output:
        outdir + "{sample}/WgsMetrics/{sample}.wgs_metrics.txt"
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
        outdir + "{sample}/Trimmomatic/{sample}.trimmomatic.stats.log",
        outdir + "{sample}/MarkDuplcatesAlign/{sample}.md_metrics.txt",
        outdir + "{sample}/Fastqc/{sample}.1P_fastqc.html",
        outdir + "{sample}/Fastqc/{sample}.2P_fastqc.html",
        outdir + "{sample}/Fastqc/{sample}.1P_fastqc.zip",
        outdir + "{sample}/Fastqc/{sample}.2P_fastqc.zip",
        outdir + "{sample}/WgsMetrics/{sample}.wgs_metrics.txt"
    output:
        outdir + "{sample}/MultiQC/{sample}.html",
        outdir + "{sample}/MultiQC/{sample}_data",
    conda:
        "envs/qc.yaml"
    threads:
        2
    resources:
        mem_mb = 2048
    shell:
        "multiqc -n {wildcards.sample} ."