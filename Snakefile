configfile: "config.json"

OUTDIR = config["WORKFLOW"]["OUTPUT"]
REF_GENOME = config["BWA"]["REF"]
SAMPLES, = glob_wildcards(config["WORKFLOW"]["READS_SOURCE"] + "/{sample}_R1.fastq.gz")

# input:

rule all:
    input:
        ""


rule trim_reads: #trimmomatic
    input:
        f_read = config["WORKFLOW"]["READS_SOURCE"] + "{sample}_R1.fastq.gz",
        r_read = config["WORKFLOW"]["READS_SOURCE"] + "{sample}_R2.fastq.gz"
    output:
        p1 = OUTDIR + "{sample}/Trimmomatic/Paired/{sample}.1P.fastq.gz",
        p2 = OUTDIR + "{sample}Trimmomatic/Paired/{sample}.2P.fastq.gz",
        u1 = OUTDIR + "{sample}Trimmomatic/Unpaired/{sample}.1U.fastq.gz",
        u2 = OUTDIR + "{sample}Trimmomatic/Unpaired/{sample}.2U.fastq.gz",
        trim_log = OUTDIR + "{sample}Trimmomatic/{sample}.trimmomatic.stats.log"
    params:
        illumina_clip = "ILLUMINACLIP:" + config["TRIMMOMATIC"]["ADAPTERS"] + ":2:30:10:2:keepBothReads",
        sliding_window = "SLIDINGWINDOW:" + config["TRIMMOMATIC"]["SLIDING_WINDOW"],
        crop = "CROP:" + config["TRIMMOMATIC"]["CROP"],
        minlen = "MINLEN:" + config["TRIMMOMATIC"]["MINLEN"]
    conda:
        config["TRIMMOMATIC"]["ENV"]
    threads:
        config["TRIMMOMATIC"]["THREADS"]
    shell:
        "trimmomatic PE -threads {threads} "
        "{input.f_read} {input.r_read} "
        "{output.p1} {output.u1} {output.p2} {output.u2} "
        "{params.illumina_clip} {params.sliding_window} "
        "{params.crop} {params.minlen} 2> {output.trim_log}"


rule run_fastqc: #fastqc
    input:
        OUTDIR + "{sample}/Trimmomatic/Paired/{sample}.1P.fastq.gz",
        OUTDIR + "{sample}Trimmomatic/Paired/{sample}.2P.fastq.gz"
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


rule align_reads: #bwa
    input:
        REF_GENOME + ".amb",
        REF_GENOME + ".ann",
        REF_GENOME + ".bwt",
        REF_GENOME + ".pac",
        REF_GENOME + ".sa",
        reference = REF_GENOME,
        p1 = OUTDIR + "{sample}/Trimmomatic/Paired/{sample}.1P.fastq.gz",
        p2 = OUTDIR + "{sample}Trimmomatic/Paired/{sample}.2P.fastq.gz",
    output:
        temp(OUTDIR + "{sample}.sam")
    conda:
        config["BWA"]["ENV"]
    threads:
        config["BWA"]["THREADS"]
    shell:
        "bwa mem -t {threads} {input.reference} "
        "{input.p1} {input.p2} > {output}"


rule sam_to_bam: #samtools
    input:
        OUTDIR + "{sample}.sam"
    output:
        bam = temp(OUTDIR + "{sample}/AlignReadsToHost/{sample}.sorted.bam"),
        bai = temp(OUTDIR + "{sample}/AlignReadsToHost/{sample}.sorted.bam.bai")
    conda:
        config["SAMTOOLS"]["ENV"]
    threads:
        config["SAMTOOLS"]["THREADS"]
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output.bam}; "
        "samtools index -@ {threads} -o {output.bai} {output.bam}"


rule mark_duplicates: #picard
    input:
        bam = OUTDIR + "{sample}/AlignReadsToHost/{sample}.sorted.bam",
        bai = OUTDIR + "{sample}/AlignReadsToHost/{sample}.sorted.bam.bai"
    output:
        md_bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        md_metrics = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.md_metrics.txt"
    params:
        tagging = "--TAGGING_POLICY " + config["PICARD"]["TAGGING"],
    conda:
        config["PICARD"]["ENV"]
    shell:
        "picard MarkDuplicates {params.tagging} "
        "-I {input.bam} -O {output.md_bam} -M {output.md_metrics}"


rule index_dupaware_bam: #samtools
    input:
        OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam"
    output:
        OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.bam.bai"
    conda:
        config["SAMTOOLS"]["ENV"]
    threads:
        config["SAMTOOLS"]["THREADS"]
    shell:
        "samtools index -@ {threads} -o {output} {input}"


rule collect_wgs_metrics: #picard
    input:
        REF_GENOME + ".fai",
        REF_GENOME + ".gzi",
        reference = REF_GENOME,
        bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.bam.bai"
    output:
        OUTDIR + "{sample}/WgsMetrics/{sample}.wgs_metrics.txt"
    conda:
        config["PICARD"]["ENV"]
    shell:
        "picard CollectWgsMetrics -I {input.bam} -O {output} -R {reference}"


rule call_snvs: #strelka2
    input:
        REF_GENOME + ".fai",
        REF_GENOME + ".gzi",
        reference = REF_GENOME,
        bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.bam.bai"
    output:
        OUTDIR + "{sample}/CallSNVs/variants.vcf.gz"
    params:
        strelka_out = OUTDIR + "{sample}/CallSNVs/"
    conda:
        config["STRELKA"]["ENV"]
    threads:
        config["STRELKA"]["THREADS"]
    shell:
        "mkdir -p {params.strelka_out}; "
        "configureStrelkaGermlineWorkflow.py --bam {input.bam} "
        "--referenceFasta {input.reference} --runDir {sample}_strelka2; "
        "{sample}_strelka2/runWorkflow.py -m local -j {threads}; "
        "mv {sample}_strelka2/results/variants/variants.vcf.gz "
        "{params.strelka_out}variants.vcf.gz; "
        "rm -rf {sample}_strelka2/"


rule call_indels: #manta
    input:
        REF_GENOME + ".fai",
        REF_GENOME + ".gzi",
        reference = REF_GENOME,
        bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.bam.bai"
    output:
        vcf = OUTDIR + "{sample}/CallInDels/diploidSV.vcf.gz",
        tbi = OUTDIR + "{sample}/CallInDels/diploidSV.vcf.gz.tbi"
    params:
        manta_out = OUTDIR + "{sample}/CallInDels/"
    conda:
        config["MANTA"]["ENV"]
    threads:
        config["MANTA"]["THREADS"]
    shell:
        "mkdir -p {params.manta_out}; "
        "configManta.py --bam {input.bam} "
        "--referenceFasta {input.reference} --runDir {sample}_manta; "
        "{sample}_manta/runWorkflow.py -j {threads}; "
        "mv {sample}_manta/results/variants/diploidSV.vcf.gz {output.vcf}; "
        "mv {sample}_manta/results/variants/diploidSV.vcf.gz.tbi {output.tbi}; "
        "rm -rf {sample}_manta"


rule call_cnvs: #cn.mops
    input:
        OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.bam.bai",
        config["CNMOPS"]["CONTROL"],
        header = config["CNMOPS"]["HEADER"]
    output:
        temp(OUTDIR + "{sample}/CallCNVs/{sample}.cnv.vcf")
    params:
        cnmops_out = "{sample}/CallCNVs/"
    conda:
        config["CNMOPS"]["ENV"]
    shell:
        "mkdir -p {params.cnmops_out}; "
        "bin/callCNV.R {wildcards.sample}; "
        "bin/csvToVCF.sh {wildcards.sample} {input.header}; "
        "mv {wildcards.sample}.cnv.vcf > {output}"


rule annotate_cnvs:
    input:
        bed = config["CNMOPS"]["BED"],
        vcf = OUTDIR + "{sample}/CallCNVs/{sample}.cnv.vcf"
    output:
        OUTDIR + "{sample}/CallCNVs/{sample}.cnv.ann.vcf"
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/annotate_cnv.py -v {input.vcf} -b {input.bed}"
        "mv {wildcards.sample}.cnv.ann.vcf {output}"


rule call_expansions: #expansion_hunter
    input:
        REF_GENOME + ".fai",
        REF_GENOME + ".gzi",
        reference = REF_GENOME,
        bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.bam.bai",
        variant_catalog = config["EXPANSIONHUNTER"]["CATALOG"]
    output:
        temp(OUTDIR + "{sample}/CallExpansions/{sample}.eh.vcf")
    params:
        eh_out = OUTDIR + "{sample}/CallExpansions/"
    conda:
        config["EXPANSIONHUNTER"]["ENV"]
    shell:
        "mkdir -p {params.eh_out}; "
        "ExpansionHunter  --reads {input.bam}  "
        "--reference {input.reference} "
        "--variant-catalog {input.variant_catalog} "
        "--output-prefix {wildcards.sample}.eh; "
        "mv {wildcards.sample}.eh.vcf > {output}"


rule annotate_expansions:
    input:
        vcf = OUTDIR + "{sample}/CallExpansions/{sample}.eh.vcf",
        variant_catalog = config["EXPANSIONHUNTER"]["CATALOG"]
    output:
        OUTDIR + "{sample}/CallExpansions/{sample}.eh.ann.vcf"
    conda:
        config["WORKFLOW"]["ENV"]
    shell:
        "bin/annotate_eh.py -v {input.vcf} -c {input.variant_catalog}; "
        "mv {wildcards.sample}.eh.ann.vcf {output}"


#rule gba_check: #gauchian


#rule cyp2d6_check: #cyrius


#rule merge_vcf: #bcftools


#rule vcf_to_parquet: