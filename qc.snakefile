rule run_fastqc: #fastqc
    input:
        OUTDIR + "{sample}/Trimmomatic/Paired/{sample}.1P.fastq.gz",
        OUTDIR + "{sample}/Trimmomatic/Paired/{sample}.2P.fastq.gz"
    output:
        OUTDIR + "{sample}/Fastqc/{sample}.1P_fastqc.html",
        OUTDIR + "{sample}/Fastqc/{sample}.2P_fastqc.html",
        OUTDIR + "{sample}/Fastqc/{sample}.1P_fastqc.zip",
        OUTDIR + "{sample}/Fastqc/{sample}.2P_fastqc.zip"
    params:
        fastqc_out = OUTDIR + "{sample}/Fastqc/"
    conda:
        config["FASTQC"]["ENV"]
    threads:
        config["FASTQC"]["THREADS"]
    shell:
        "mkdir -p {params.fastqc_out}; "
        "fastqc -t {threads} -o {params.fastqc_out} -f fastq {input}"


rule collect_wgs_metrics: #picard
    input:
        REF_GENOME + ".fai",
        REF_GENOME + ".gzi",
        reference = REF_GENOME,
        bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai"
    output:
        OUTDIR + "{sample}/WgsMetrics/{sample}.wgs_metrics.txt"
    conda:
        config["PICARD"]["ENV"]
    shell:
        "picard CollectWgsMetrics -I {input.bam} -O {output} -R {reference}"


rule run_multiqc:
    input:
        OUTDIR + "{sample}/Trimmomatic/{sample}.trimmomatic.stats.log",
        OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.md_metrics.txt",
        OUTDIR + "{sample}/Fastqc/{sample}.1P_fastqc.html",
        OUTDIR + "{sample}/Fastqc/{sample}.2P_fastqc.html",
        OUTDIR + "{sample}/Fastqc/{sample}.1P_fastqc.zip",
        OUTDIR + "{sample}/Fastqc/{sample}.2P_fastqc.zip",
        OUTDIR + "{sample}/WgsMetrics/{sample}.wgs_metrics.txt"
    output:
        OUTDIR + "{sample}/MultiQC/{sample}.html",
        OUTDIR + "{sample}/MultiQC/{sample}_data",
    conda:
        config["MULTIQC"]["ENV"]
    shell:
        "multiqc -n {wildcards.sample} ."