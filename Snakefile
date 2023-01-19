# from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
# S3 = S3RemoteProvider()
# f_read = S3.remote("nf-wgs-ufl/snakemake/test-143-10x_L001_R1.fastq.gz"),
# for running on a single ec2 instance
# look at step functions
configfile: "config.yaml"

reads = config["reads"]
outdir = config["output"]
ref_genome = config["ref_genome"]
samples, = glob_wildcards(config["reads"] + "/{sample}_R1.fastq.gz")

all_input = [
    expand(outdir + "{sample}/ParquetVCF/{sample}.SNV.parquet", sample = samples),
    expand(outdir + "{sample}/ParquetVCF/{sample}.CNV.parquet", sample = samples),
    expand(outdir + "{sample}/ParquetVCF/{sample}.EXP.parquet", sample = samples),
    expand(outdir + "{sample}/ParquetVCF/{sample}.INDEL.parquet", sample = samples)
]

if config["WORKFLOW"]["QC"].upper() == "TRUE":
    include: "qc.snakefile"
    qc_results = [
        expand(outdir + "{sample}/MultiQC/{sample}.html", sample = samples),
        expand(outdir + "{sample}/MultiQC/{sample}_data", sample = samples)
    ]
    all_input.append(qc_results)


rule all:
    input:
        all_input


rule trim_reads: #trimmomatic
    input:
        f_read = reads + "{sample}_R1.fastq.gz",
        r_read = reads + "{sample}_R2.fastq.gz"
    output:
        p1 = outdir + "{sample}/Trimmomatic/Paired/{sample}.1P.fastq.gz",
        p2 = outdir + "{sample}/Trimmomatic/Paired/{sample}.2P.fastq.gz",
        u1 = outdir + "{sample}/Trimmomatic/Unpaired/{sample}.1U.fastq.gz",
        u2 = outdir + "{sample}/Trimmomatic/Unpaired/{sample}.2U.fastq.gz",
        trim_log = outdir + "{sample}/Trimmomatic/{sample}.trimmomatic.stats.log"
    params:
        illumina_clip = "ILLUMINACLIP:" + config["adapters"] + ":2:30:10:2:keepBothReads",
        sliding_window = "SLIDINGWINDOW:" + config["sliding_window"],
        crop = "CROP:" + config["crop"],
        minlen = "MINLEN:" + config["MINLEN"]
    conda:
        "envs/read_trimming.yaml"
    threads:
        8
    resources:
        mem_mb = 14336
    shell:
        "trimmomatic PE -threads {threads} "
        "{input.f_read} {input.r_read} "
        "{output.p1} {output.u1} {output.p2} {output.u2} "
        "{params.illumina_clip} {params.sliding_window} "
        "{params.crop} {params.minlen} 2> {output.trim_log}"


rule align_reads: #bwa
    input:
        ref_genome + ".amb",
        ref_genome + ".ann",
        ref_genome + ".bwt",
        ref_genome + ".pac",
        ref_genome + ".sa",
        reference = ref_genome,
        p1 = outdir + "{sample}/Trimmomatic/Paired/{sample}.1P.fastq.gz",
        p2 = outdir + "{sample}/Trimmomatic/Paired/{sample}.2P.fastq.gz",
    output:
        temp(outdir + "{sample}.sam")
    conda:
        "envs/alignment.yaml"
    threads:
        32
    resources:
        mem_mb = 63488
    shell:
        "bwa mem -t {threads} {input.reference} "
        "{input.p1} {input.p2} > {output}"


rule sam_to_bam: #samtools
    input:
        outdir + "{sample}.sam"
    output:
        bam = temp(outdir + "{sample}/AlignReadsToHost/{sample}.sorted.bam"),
        bai = temp(outdir + "{sample}/AlignReadsToHost/{sample}.sorted.bam.bai")
    conda:
        "envs/alignment.yaml"
    threads:
        8
    resources:
        mem_mb = 63488
    shell:
        "samtools view -Sb {input} | "
        "samtools sort -@ {threads} -o {output.bam}; "
        "samtools index -@ {threads} {output.bam} {output.bai}"


rule mark_duplicates: #picard
    input:
        bam = outdir + "{sample}/AlignReadsToHost/{sample}.sorted.bam",
        bai = outdir + "{sample}/AlignReadsToHost/{sample}.sorted.bam.bai"
    output:
        md_bam = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        md_metrics = outdir + "{sample}/MarkDuplcatesAlign/{sample}.md_metrics.txt"
    params:
        tagging = "--TAGGING_POLICY " + config["tagging"],
    conda:
        "envs/qc.yaml"
    threads:
        4
    resources:
        mem_mb = 63488
    shell:
        "picard MarkDuplicates {params.tagging} "
        "-I {input.bam} -O {output.md_bam} -M {output.md_metrics}"


rule index_dupaware_bam: #samtools
    input:
        outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam"
    output:
        outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai"
    conda:
        "envs/alignment.yaml"
    threads:
        8
    resources:
        mem_mb = 63488
    shell:
        "samtools index -@ {threads} {input} {output}"


rule call_snvs: #strelka2
    input:
        ref_genome + ".fai",
        ref_genome + ".gzi",
        reference = ref_genome,
        bam = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai"
    output:
        outdir + "{sample}/CallSNVs/variants.vcf.gz"
    params:
        strelka_out = outdir + "{sample}/CallSNVs/",
        temp_out = "{sample}_strelka2"
    conda:
        "envs/variant_calling.yaml"
    threads:
        8
    resources:
        mem_mb = 14336
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
        ref_genome + ".fai",
        ref_genome + ".gzi",
        reference = ref_genome,
        bam = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai"
    output:
        vcf = outdir + "{sample}/CallInDels/diploidSV.vcf.gz",
        tbi = outdir + "{sample}/CallInDels/diploidSV.vcf.gz.tbi"
    params:
        manta_out = outdir + "{sample}/CallInDels/",
        temp_out = "{sample}_manta"
    conda:
        "envs/variant_calling.yaml"
    threads:
        8
    resources:
        mem_mb = 14336
    shell:
        "mkdir -p {params.manta_out}; "
        "configManta.py --bam {input.bam} "
        "--referenceFasta {input.reference} --runDir {params.temp_out}; "
        "{params.temp_out}/runWorkflow.py -j {threads}; "
        "mv {params.temp_out}/results/variants/diploidSV.vcf.gz {output.vcf}; "
        "mv {params.temp_out}/results/variants/diploidSV.vcf.gz.tbi {output.tbi}; "
        "rm -rf {params.temp_out}"


rule call_cnvs: #cn.mops MODIFY THIS RULE AFTER MODDING SCRIPT
    input:
        bam = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai",
        r_control = config["controls"],
        # header = config["CNMOPS"]["HEADER"]
    output:
        csv = temp(outdir + "{sample}/CallCNVs/{sample}.cnv.csv"),
        head = temp(outdir + "{sample}/CallCNVs/{sample}.cnv.head"),
        tmp = temp(outdir + "{sample}/CallCNVs/{sample}.cnv.tmp"),
        vcf = temp(outdir + "{sample}/CallCNVs/{sample}.cnv.vcf")
    conda:
        "envs/variant_calling.yaml"
    threads:
        8
    resources:
        mem_mb = 14336
    shell:
        "bin/callCNV.R {input.r_control} {input.bam} {output.csv}; "
        "bin/csvToVCF.sh {input.header} {output.head} {output.tmp} "
        "{output.csv} {output.vcf}; "


rule annotate_cnvs: # MODIFY THIS RULE AFTER MODDING SCRIPT
    input:
        bed = config["bed"],
        vcf = outdir + "{sample}/CallCNVs/{sample}.cnv.vcf"
    output:
        ann_vcf = outdir + "{sample}/CallCNVs/{sample}.cnv.ann.vcf",
        gz_vcf = temp(outdir + "{sample}/CallCNVs/{sample}.cnv.vcf.gz"),
        tbi = temp(outdir + "{sample}/CallCNVs/{sample}.cnv.vcf.gz.tbi")
    conda:
        "envs/default.yaml"
    threads:
        4
    resources:
        mem_mb = 8192
    shell:
        "bin/annotate_cnv.py "
        "-v {input.vcf} "
        "-b {input.bed} "
        "-o {output.ann_vcf}"


rule call_expansions: #expansion_hunter
    input:
        ref_genome + ".fai",
        ref_genome + ".gzi",
        reference = ref_genome,
        bam = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam",
        bai = outdir + "{sample}/MarkDuplcatesAlign/{sample}.sorted.md.bam.bai",
        variant_catalog = config["catalog"]
    output:
        vcf = temp(outdir + "{sample}/CallExpansions/{sample}.eh.vcf"),
        json = temp(outdir + "{sample}/CallExpansions/{sample}.eh.json"),
        bam = temp(outdir + "{sample}/CallExpansions/{sample}.eh_realigned.bam")
    params:
        eh_out = outdir + "{sample}/CallExpansions/{sample}.eh"
    conda:
        "envs/variant_calling.yaml"
    threads:
        4
    resources:
        mem_mb = 8192
    shell:
        "ExpansionHunter "
        "--reads {input.bam} "
        "--reference {input.reference} "
        "--variant-catalog {input.variant_catalog} "
        "--output-prefix {params.eh_out}"


rule annotate_expansions:
    input:
        vcf = outdir + "{sample}/CallExpansions/{sample}.eh.vcf",
        variant_catalog = config["catalog"]
    output:
        ann_vcf = outdir + "{sample}/CallExpansions/{sample}.eh.ann.vcf",
        gz_vcf = temp(outdir + "{sample}/CallExpansions/{sample}.eh.vcf.gz"),
        tbi = temp(outdir + "{sample}/CallExpansions/{sample}.eh.vcf.gz.tbi")
    conda:
        "envs/default.yaml"
    threads:
        4
    resources:
        mem_mb = 8192
    shell:
        "bin/annotate_eh.py "
        "-v {input.vcf} "
        "-c {input.variant_catalog} "
        "-o {output.ann_vcf}"


if config["qc"].upper() == "TRUE":
    include: "qc.snakefile"


rule zip_vcf:
    input:
        cnv_vcf = outdir + "{sample}/CallCNVs/{sample}.cnv.ann.vcf",
        exp_vcf = outdir + "{sample}/CallExpansions/{sample}.eh.ann.vcf"
    output:
        outdir + "{sample}/CallCNVs/{sample}.cnv.ann.vcf",
        outdir + "{sample}/CallExpansions/{sample}.eh.ann.vcf"
    conda:
        "envs/variant_calling.yaml"
    threads:
        2
    resources:
        mem_mb = 2048
    shell:
        "tabix {input.cnv_vcf} "
        "tabix {input.exp_vcf}"


rule index_vcf:
    input:
        snv_vcf = outdir + "{sample}/CallSNVs/variants.vcf.gz",
        cnv_vcf = outdir + "{sample}/CallCNVs/{sample}.cnv.ann.vcf.gz",
        exp_vcf = outdir + "{sample}/CallExpansions/{sample}.eh.ann.vcf.gz"
    output:
        outdir + "{sample}/CallSNVs/variants.vcf.gz.tbi",
        outdir + "{sample}/CallCNVs/{sample}.cnv.ann.vcf.gz.tbi",
        outdir + "{sample}/CallExpansions/{sample}.eh.ann.vcf.gz.tbi"
    conda:
        "envs/variant_calling.yaml"
    threads:
        2
    resources:
        mem_mb = 2048
    shell:
        "bcftools index --tbi {input.snv_vcf} "
        "bcftools index --tbi {input.cnv_vcf} "
        "bcftools index --tbi {input.exp_vcf}"


rule vcf_to_parquet:
    input:
        snv_vcf = outdir + "{sample}/CallSNVs/variants.vcf.gz",
        snv_tbi = outdir + "{sample}/CallSNVs/variants.vcf.gz.tbi",
        indel_vcf = outdir + "{sample}/CallInDels/diploidSV.vcf.gz",
        indel_tbi = outdir + "{sample}/CallInDels/diploidSV.vcf.gz.tbi",
        cnv_vcf = outdir + "{sample}/CallCNVs/{sample}.cnv.ann.vcf",
        cnv_tbi = outdir + "{sample}/CallCNVs/{sample}.cnv.ann.vcf.gz.tbi",
        exp_vcf = outdir + "{sample}/CallExpansions/{sample}.eh.ann.vcf",
        exp_tbi = outdir + "{sample}/CallExpansions/{sample}.eh.ann.vcf.gz.tbi"
    output:
        outdir + "{sample}/ParquetVCF/{sample}.SNV.parquet",
        outdir + "{sample}/ParquetVCF/{sample}.CNV.parquet",
        outdir + "{sample}/ParquetVCF/{sample}.EXP.parquet",
        outdir + "{sample}/ParquetVCF/{sample}.INDEL.parquet"
    params:
        parquet_out = outdir + "{sample}/ParquetVCF"
    conda:
        "envs/default.yaml"
    threads:
        4
    resources:
        mem_mb = 8192
    shell:
        "bin/vcf_to_parquet.py "
        "-s {wildcards.sample} "
        "-v {input.snv_vcf} {input.indel_vcf} {input.cnv_vcf} {input.exp_vcf} "
        "mv {sample}.SNV.parquet {params.parquet_out}/{sample}.SNV.parquet "
        "mv {sample}.CNV.parquet {params.parquet_out}/{sample}.CNV.parquet "
        "mv {sample}.EXP.parquet {params.parquet_out}/{sample}.EXP.parquet "
        "mv {sample}.INDEL.parquet {params.parquet_out}/{sample}.INDEL.parquet"