# from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider()
# f_read = S3.remote("nf-wgs-ufl/snakemake/test-143-10x_L001_R1.fastq.gz"),
# for running on a single ec2 instance
# look at step functions
configfile: "config.json"

OUTDIR = config["WORKFLOW"]["OUTPUT"]
REF_GENOME = config["BWA"]["REF_GENOME"]
SAMPLES, = glob_wildcards(config["WORKFLOW"]["READS_SOURCE"] + "/{sample}_R1.fastq.gz")

all_input = [
    expand(OUTDIR + "{sample}/ParquetVCF/{sample}.SNV.parquet", sample = SAMPLES),
    expand(OUTDIR + "{sample}/ParquetVCF/{sample}.CNV.parquet", sample = SAMPLES),
    expand(OUTDIR + "{sample}/ParquetVCF/{sample}.EXP.parquet", sample = SAMPLES),
    expand(OUTDIR + "{sample}/ParquetVCF/{sample}.INDEL.parquet", sample = SAMPLES)
]

if config["WORKFLOW"]["QC"].upper() == "TRUE":
    include: "qc.snakefile"
    qc_results = [
        expand(OUTDIR + "{sample}/MultiQC/{sample}.html", sample = SAMPLES),
        expand(OUTDIR + "{sample}/MultiQC/{sample}_data", sample = SAMPLES)
    ]
    all_input.append(qc_results)


rule all:
    input:
        all_input


rule trim_reads: #trimmomatic
    input:
        f_read = config["WORKFLOW"]["READS_SOURCE"] + "{sample}_R1.fastq.gz",
        r_read = config["WORKFLOW"]["READS_SOURCE"] + "{sample}_R2.fastq.gz"
    output:
        p1 = OUTDIR + "{sample}/Trimmomatic/Paired/{sample}.1P.fastq.gz",
        p2 = OUTDIR + "{sample}/Trimmomatic/Paired/{sample}.2P.fastq.gz",
        u1 = OUTDIR + "{sample}/Trimmomatic/Unpaired/{sample}.1U.fastq.gz",
        u2 = OUTDIR + "{sample}/Trimmomatic/Unpaired/{sample}.2U.fastq.gz",
        trim_log = OUTDIR + "{sample}/Trimmomatic/{sample}.trimmomatic.stats.log"
    params:
        illumina_clip = "ILLUMINACLIP:" + config["TRIMMOMATIC"]["ADAPTERS"] + ":2:30:10:2:keepBothReads",
        sliding_window = "SLIDINGWINDOW:" + config["TRIMMOMATIC"]["SLIDING_WINDOW"],
        crop = "CROP:" + config["TRIMMOMATIC"]["CROP"],
        minlen = "MINLEN:" + config["TRIMMOMATIC"]["MINLEN"]
    conda:
        "envs/read_trimming.yaml"
    threads:
        config["TRIMMOMATIC"]["THREADS"]
    shell:
        "trimmomatic PE -threads {threads} "
        "{input.f_read} {input.r_read} "
        "{output.p1} {output.u1} {output.p2} {output.u2} "
        "{params.illumina_clip} {params.sliding_window} "
        "{params.crop} {params.minlen} 2> {output.trim_log}"


rule align_reads: #bwa
    input:
        REF_GENOME + ".amb",
        REF_GENOME + ".ann",
        REF_GENOME + ".bwt",
        REF_GENOME + ".pac",
        REF_GENOME + ".sa",
        reference = REF_GENOME,
        p1 = OUTDIR + "{sample}/Trimmomatic/Paired/{sample}.1P.fastq.gz",
        p2 = OUTDIR + "{sample}/Trimmomatic/Paired/{sample}.2P.fastq.gz",
    output:
        temp(OUTDIR + "{sample}.sam")
    conda:
        "envs/alignment.yaml"
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
        "envs/alignment.yaml"
    threads:
        config["SAMTOOLS"]["THREADS"]
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output.bam}; "
        "samtools index -@ {threads} {output.bam} {output.bai}"


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
        "envs/qc.yaml"
    shell:
        "picard MarkDuplicates {params.tagging} "
        "-I {input.bam} -O {output.md_bam} -M {output.md_metrics}"


rule index_dupaware_bam: #samtools
    input:
        OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam"
    output:
        OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai"
    conda:
        "envs/alignment.yaml"
    threads:
        config["SAMTOOLS"]["THREADS"]
    shell:
        "samtools index -@ {threads} {input} {output}"


rule call_snvs: #strelka2
    input:
        REF_GENOME + ".fai",
        REF_GENOME + ".gzi",
        reference = REF_GENOME,
        bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai"
    output:
        OUTDIR + "{sample}/CallSNVs/variants.vcf.gz"
    params:
        strelka_out = OUTDIR + "{sample}/CallSNVs/",
        temp_out = "{sample}_strelka2"
    conda:
        "envs/variant_calling.yaml"
    threads:
        config["STRELKA"]["THREADS"]
    shell:
        "mkdir -p {params.strelka_out}; "
        "configureStrelkaGermlineWorkflow.py --bam {input.bam} "
        "--referenceFasta {input.reference} --runDir {params.temp_out}; "
        "{params.temp_out}/runWorkflow.py -m local -j {threads}; "
        "mv {params.temp_out}/results/variants/variants.vcf.gz "
        "{params.strelka_out}variants.vcf.gz; "
        "rm -rf {params.temp_out}/"


rule call_indels: #manta
    input:
        REF_GENOME + ".fai",
        REF_GENOME + ".gzi",
        reference = REF_GENOME,
        bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai"
    output:
        vcf = OUTDIR + "{sample}/CallInDels/diploidSV.vcf.gz",
        tbi = OUTDIR + "{sample}/CallInDels/diploidSV.vcf.gz.tbi"
    params:
        manta_out = OUTDIR + "{sample}/CallInDels/",
        temp_out = "{sample}_manta"
    conda:
        "envs/variant_calling.yaml"
    threads:
        config["MANTA"]["THREADS"]
    shell:
        "mkdir -p {params.manta_out}; "
        "configManta.py --bam {input.bam} "
        "--referenceFasta {input.reference} --runDir {params.temp_out}; "
        "{params.temp_out}/runWorkflow.py -j {threads}; "
        "mv {params.temp_out}/results/variants/diploidSV.vcf.gz {output.vcf}; "
        "mv {params.temp_out}/results/variants/diploidSV.vcf.gz.tbi {output.tbi}; "
        "rm -rf {params.temp_out}"


rule call_cnvs: #cn.mops
    input:
        bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai",
        r_control = config["CNMOPS"]["CONTROL"],
        header = config["CNMOPS"]["HEADER"]
    output:
        csv = temp(OUTDIR + "{sample}/CallCNVs/{sample}.cnv.csv"),
        head = temp(OUTDIR + "{sample}/CallCNVs/{sample}.cnv.head"),
        tmp = temp(OUTDIR + "{sample}/CallCNVs/{sample}.cnv.tmp"),
        vcf = temp(OUTDIR + "{sample}/CallCNVs/{sample}.cnv.vcf")
    conda:
        "envs/variant_calling.yaml"
    shell:
        "bin/callCNV.R {input.r_control} {input.bam} {output.csv}; "
        "bin/csvToVCF.sh {input.header} {output.head} {output.tmp} "
        "{output.csv} {output.vcf}; "


rule annotate_cnvs:
    input:
        bed = config["CNMOPS"]["BED"],
        vcf = OUTDIR + "{sample}/CallCNVs/{sample}.cnv.vcf"
    output:
        ann_vcf = OUTDIR + "{sample}/CallCNVs/{sample}.cnv.ann.vcf",
        gz_vcf = temp(OUTDIR + "{sample}/CallCNVs/{sample}.cnv.vcf.gz"),
        tbi = temp(OUTDIR + "{sample}/CallCNVs/{sample}.cnv.vcf.gz.tbi")
    conda:
        "envs/default.yaml"
    shell:
        "bin/annotate_cnv.py "
        "-v {input.vcf} "
        "-b {input.bed} "
        "-o {output.ann_vcf}"


rule call_expansions: #expansion_hunter
    input:
        REF_GENOME + ".fai",
        REF_GENOME + ".gzi",
        reference = REF_GENOME,
        bam = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = OUTDIR + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai",
        variant_catalog = config["EXPANSIONHUNTER"]["CATALOG"]
    output:
        vcf = temp(OUTDIR + "{sample}/CallExpansions/{sample}.eh.vcf"),
        json = temp(OUTDIR + "{sample}/CallExpansions/{sample}.eh.json"),
        bam = temp(OUTDIR + "{sample}/CallExpansions/{sample}.eh_realigned.bam")
    params:
        eh_out = OUTDIR + "{sample}/CallExpansions/{sample}.eh"
    conda:
        "envs/variant_calling.yaml"
    shell:
        "ExpansionHunter "
        "--reads {input.bam} "
        "--reference {input.reference} "
        "--variant-catalog {input.variant_catalog} "
        "--output-prefix {params.eh_out}"


rule annotate_expansions:
    input:
        vcf = OUTDIR + "{sample}/CallExpansions/{sample}.eh.vcf",
        variant_catalog = config["EXPANSIONHUNTER"]["CATALOG"]
    output:
        ann_vcf = OUTDIR + "{sample}/CallExpansions/{sample}.eh.ann.vcf",
        gz_vcf = temp(OUTDIR + "{sample}/CallExpansions/{sample}.eh.vcf.gz"),
        tbi = temp(OUTDIR + "{sample}/CallExpansions/{sample}.eh.vcf.gz.tbi")
    conda:
        "envs/default.yaml"
    shell:
        "bin/annotate_eh.py "
        "-v {input.vcf} "
        "-c {input.variant_catalog} "
        "-o {output.ann_vcf}"


if config["WORKFLOW"]["QC"].upper() == "TRUE":
    include: "qc.snakefile"


rule zip_vcf:
    input:
        cnv_vcf = OUTDIR + "{sample}/CallCNVs/{sample}.cnv.ann.vcf",
        exp_vcf = OUTDIR + "{sample}/CallExpansions/{sample}.eh.ann.vcf"
    output:
        OUTDIR + "{sample}/CallCNVs/{sample}.cnv.ann.vcf",
        OUTDIR + "{sample}/CallExpansions/{sample}.eh.ann.vcf"
    conda:
        "envs/variant_calling.yaml"
    shell:
        "tabix {input.cnv_vcf} "
        "tabix {input.exp_vcf}"


rule index_vcf:
    input:
        snv_vcf = OUTDIR + "{sample}/CallSNVs/variants.vcf.gz",
        cnv_vcf = OUTDIR + "{sample}/CallCNVs/{sample}.cnv.ann.vcf.gz",
        exp_vcf = OUTDIR + "{sample}/CallExpansions/{sample}.eh.ann.vcf.gz"
    output:
        OUTDIR + "{sample}/CallSNVs/variants.vcf.gz.tbi",
        OUTDIR + "{sample}/CallCNVs/{sample}.cnv.ann.vcf.gz.tbi",
        OUTDIR + "{sample}/CallExpansions/{sample}.eh.ann.vcf.gz.tbi"
    conda:
        "envs/variant_calling.yaml"
    shell:
        "bcftools index --tbi {input.snv_vcf} "
        "bcftools index --tbi {input.cnv_vcf} "
        "bcftools index --tbi {input.exp_vcf}"


rule vcf_to_parquet:
    input:
        snv_vcf = OUTDIR + "{sample}/CallSNVs/variants.vcf.gz",
        snv_tbi = OUTDIR + "{sample}/CallSNVs/variants.vcf.gz.tbi",
        indel_vcf = OUTDIR + "{sample}/CallInDels/diploidSV.vcf.gz",
        indel_tbi = OUTDIR + "{sample}/CallInDels/diploidSV.vcf.gz.tbi",
        cnv_vcf = OUTDIR + "{sample}/CallCNVs/{sample}.cnv.ann.vcf",
        cnv_tbi = OUTDIR + "{sample}/CallCNVs/{sample}.cnv.ann.vcf.gz.tbi",
        exp_vcf = OUTDIR + "{sample}/CallExpansions/{sample}.eh.ann.vcf",
        exp_tbi = OUTDIR + "{sample}/CallExpansions/{sample}.eh.ann.vcf.gz.tbi"
    output:
        OUTDIR + "{sample}/ParquetVCF/{sample}.SNV.parquet",
        OUTDIR + "{sample}/ParquetVCF/{sample}.CNV.parquet",
        OUTDIR + "{sample}/ParquetVCF/{sample}.EXP.parquet",
        OUTDIR + "{sample}/ParquetVCF/{sample}.INDEL.parquet"
    params:
        parquet_out = OUTDIR + "{sample}/ParquetVCF"
    conda:
        "envs/default.yaml"
    shell:
        "bin/vcf_to_parquet.py "
        "-s {wildcards.sample} "
        "-v {input.snv_vcf} {input.indel_vcf} {input.cnv_vcf} {input.exp_vcf} "
        "mv {sample}.SNV.parquet {params.parquet_out}/{sample}.SNV.parquet "
        "mv {sample}.CNV.parquet {params.parquet_out}/{sample}.CNV.parquet "
        "mv {sample}.EXP.parquet {params.parquet_out}/{sample}.EXP.parquet "
        "mv {sample}.INDEL.parquet {params.parquet_out}/{sample}.INDEL.parquet"