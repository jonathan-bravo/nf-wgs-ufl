configfile: "config.yaml"

reads = config["reads"]
outdir = config["output"]
ref_genome = config["ref_genome"]
samples, = glob_wildcards(reads + "/{sample}_R1.fastq.gz")

# suffixes for files in the pipeline
align_suffixes = ('.amb', 'ann', 'bwt', '.pac', '.sa')
ref_suffixes = ('.fai', '.gzi')
variant_suffixes = ('snv.vcf.gz', 'indel.vcf.gz', 'cnv.vcf.gz', 'eh.vcf.gz')
parquet_suffixes = ('SNV', 'CNV', 'EXP', 'INDEL')
multiqc_suffixes = ('.html', '_data')


all_input = [
    expand(
        outdir + "{sample}/ParquetVCF/{sample}.{suffix}.parquet",
        sample = samples,
        suffix = parquet_suffixes
    )
]

if config["qc"].upper() == "TRUE":
    include: "qc.snakefile"
    qc_results = [
        expand(
            outdir + "{sample}/MultiQC/{sample}{suffix}",
            sample = samples,
            suffix = multiqc_suffixes
        )
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
        p1 = temp(outdir + "{sample}/{sample}.1P.fastq.gz"),
        p2 = temp(outdir + "{sample}/{sample}.2P.fastq.gz"),
        u1 = temp(outdir + "{sample}/{sample}.1U.fastq.gz"),
        u2 = temp(outdir + "{sample}/{sample}.2U.fastq.gz"),
        trim_log = outdir + "{sample}/QC/{sample}.trimmomatic.stats.log"
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
        expand(ref_genome + suffix, suffix = align_suffixes)
        reference = ref_genome,
        p1 = outdir + "{sample}/{sample}.1P.fastq.gz",
        p2 = outdir + "{sample}/{sample}.2P.fastq.gz",
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
        bam = temp(outdir + "{sample}/{sample}.sorted.bam"),
        bai = temp(outdir + "{sample}/{sample}.sorted.bam.bai")
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
        outdir + "{sample}/{sample}.sorted.bam.bai",
        bam = outdir + "{sample}/{sample}.sorted.bam",
    output:
        md_bam = outdir + "{sample}/Alignment/{sample}.sorted.md.bam",
        md_metrics = outdir + "{sample}/QC/{sample}.md_metrics.txt"
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
        outdir + "{sample}/Alignment/{sample}.sorted.md.bam"
    output:
        outdir + "{sample}/Alignment/{sample}.sorted.md.bam.bai"
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
        expand(ref_genome + suffix, suffix = ref_suffixes),
        outdir + "{sample}/Alignment/{sample}.sorted.md.bam.bai",
        bam = outdir + "{sample}/Alignment/{sample}.sorted.md.bam",
        reference = ref_genome
    output:
        outdir + "{sample}/Variants/{sample}.snv.vcf.gz"
    params:
        strelka_out = outdir + "{sample}_strelka2"
    conda:
        "envs/variant_calling.yaml"
    threads:
        8
    resources:
        mem_mb = 14336
    shell:
        "mkdir -p {params.strelka_out}; "
        "configureStrelkaGermlineWorkflow.py "
        "--bam {input.bam} "
        "--referenceFasta {input.reference} "
        "--runDir {params.strelka_out}; "
        "{params.strelka_out}/runWorkflow.py -m local -j {threads}; "
        "mv {params.strelka_out}/results/variants/variants.vcf.gz {output}; "
        "rm -rf {params.strelka_out}/"


rule call_indels: #manta
    input:
        expand(ref_genome + suffix, suffix = ref_suffixes),
        outdir + "{sample}/Alignment/{sample}.sorted.md.bam.bai",
        bam = outdir + "{sample}/Alignment/{sample}.sorted.md.bam",
        reference = ref_genome
    output:
        vcf = outdir + "{sample}/Variants/{sample}.indel.vcf.gz",
        tbi = outdir + "{sample}/Variants/{sample}.indel.vcf.gz.tbi"
    params:
        manta_out = outdir + "{sample}_manta"
    conda:
        "envs/variant_calling.yaml"
    threads:
        8
    resources:
        mem_mb = 14336
    shell:
        "mkdir -p {params.manta_out}; "
        "configManta.py --bam {input.bam} "
        "--referenceFasta {input.reference} "
        "--runDir {params.temp_out}; "
        "{params.manta_out}/runWorkflow.py -j {threads}; "
        "mv {params.manta_out}/results/variants/diploidSV.vcf.gz {output.vcf}; "
        "mv {params.manta_out}/results/variants/diploidSV.vcf.gz.tbi {output.tbi}; "
        "rm -rf {params.manta_out}"


rule call_cnvs: #cn.mops
    input:
        outdir + "{sample}/Alignment/{sample}.sorted.md.bam.bai",
        bam = outdir + "{sample}/Alignment/{sample}.sorted.md.bam",
        r_control = config["controls"]
    output:
        csv = temp(outdir + "{sample}/Variants/{sample}.cnv.csv"),
    conda:
        "envs/variant_calling.yaml"
    threads:
        8
    resources:
        mem_mb = 14336
    shell:
        "bin/callCNV.R {input.r_control} {input.bam} {output.csv}"


rule annotate_cnvs:
    input:
        bed = config["bed"],
        csv = outdir + "{sample}/Variants/{sample}.cnv.csv"
    output:
        ann_vcf = temp(outdir + "{sample}/Variants/{sample}.cnv.ann.vcf")
    params:
        vcf = outdir + "{sample}/Variants/{sample}.cnv.vcf",
        gz = outdir + "{sample}/Variants/{sample}.cnv.vcf.gz",
        tbi = outdir + "{sample}/Variants/{sample}.cnv.vcf.gz.tbi"
    conda:
        "envs/default.yaml"
    threads:
        4
    resources:
        mem_mb = 8192
    shell:
        "bin/annotate_cnv.py "
        "-s {wildecards.sample} "
        "-c {input.csv} "
        "-b {input.bed} "
        "-o {output.ann_vcf}; "
        "rm {params.vcf} {params.gz} {params.tbi}"


rule call_expansions: #expansion_hunter
    input:
        expand(ref_genome + suffix, suffix = ref_suffixes),
        outdir + "{sample}/Alignment/{sample}.sorted.md.bam.bai",
        bam = outdir + "{sample}/Alignment/{sample}.sorted.md.bam",
        reference = ref_genome,
        variant_catalog = config["catalog"]
    output:
        vcf = temp(outdir + "{sample}/Variants/{sample}.eh.vcf"),
    params:
        eh_out = outdir + "{sample}/Variants/{sample}.eh"
        eh_json = outdir + "{sample}/Variants/{sample}.eh.json",
        eh_bam = outdir + "{sample}/Variants/{sample}.eh_realigned.bam"
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
        "--output-prefix {params.eh_out}; "
        "rm {params.eh_json} {params.eh_bam}"


rule annotate_expansions:
    input:
        vcf = outdir + "{sample}/Variants/{sample}.eh.vcf",
        variant_catalog = config["catalog"]
    output:
        ann_vcf = temp(outdir + "{sample}/Variants/{sample}.eh.ann.vcf"),
    params:
        gz = outdir + "{sample}/Variants/{sample}.eh.vcf.gz",
        tbi = outdir + "{sample}/Variants/{sample}.eh.vcf.gz.tbi"
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
        "-o {output.ann_vcf}; "
        "rm {params.gz} {params.tbi}"


rule zip_vcf:
    input:
        cnv_vcf = outdir + "{sample}/Variants/{sample}.cnv.ann.vcf",
        exp_vcf = outdir + "{sample}/Variants/{sample}.eh.ann.vcf"
    output:
        outdir + "{sample}/Variants/{sample}.cnv.ann.vcf.gz",
        outdir + "{sample}/Variants/{sample}.eh.ann.vcf.gz"
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
        snv_vcf = outdir + "{sample}/Variants/{sample}.snv.vcf.gz",
        cnv_vcf = outdir + "{sample}/Variants/{sample}.cnv.ann.vcf.gz",
        exp_vcf = outdir + "{sample}/Variants/{sample}.eh.ann.vcf.gz"
    output:
        outdir + "{sample}/Variants/{sample}.snv.vcf.gz.tbi",
        outdir + "{sample}/Variants/{sample}.cnv.ann.vcf.gz.tbi",
        outdir + "{sample}/Variants/{sample}.eh.ann.vcf.gz.tbi"
    conda:
        "envs/variant_calling.yaml"
    threads:
        2
    resources:
        mem_mb = 2048
    shell:
        "bcftools index --tbi {input.snv_vcf}; "
        "bcftools index --tbi {input.cnv_vcf}; "
        "bcftools index --tbi {input.exp_vcf}"


rule vcf_to_parquet:
    input:
        expand(outdir + "{sample}/Variants/{sample}.{suffix}.tbi", suffix = variant_suffixes),
        vcfs = expand(outdir + "{sample}/Variants/{sample}.{suffix}", suffix = variant_suffixes)
    output:
        expand(outdir + "{sample}/Parquet/{sample}.{suffix}.parquet", suffix = parquet_suffixes)
    params:
        parquet_out = outdir + "{sample}/Parquet/{sample}"
    conda:
        "envs/default.yaml"
    threads:
        4
    resources:
        mem_mb = 8192
    shell:
        "bin/vcf_to_parquet.py "
        "-s {wildcards.sample} "
        "-v {input.vcfs}; "
        "mv {wildcards.sample}.SNV.parquet {params.parquet_out}.SNV.parquet; "
        "mv {wildcards.sample}.CNV.parquet {params.parquet_out}.CNV.parquet; "
        "mv {wildcards.sample}.EXP.parquet {params.parquet_out}.EXP.parquet; "
        "mv {wildcards.sample}.INDEL.parquet {params.parquet_out}.INDEL.parquet"