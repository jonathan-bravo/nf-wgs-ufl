#!/usr/bin/env nextflow

params.single_lane = ""

params.bucket      = ""
params.run_id      = ""
params.ref_dir     = "${params.bucket}/Pipeline/Reference"

params.reference   = "${params.ref_dir}/hg19/hg19.fa"
params.bwa_amb     = "${params.reference}.amb"
params.bwa_ann     = "${params.reference}.ann"
params.bwa_bwt     = "${params.reference}.bwt"
params.bwa_pac     = "${params.reference}.pac"
params.bwa_sa      = "${params.reference}.sa"
params.ref_fai     = "${params.reference}.fai"

params.cnv_control = "${params.ref_dir}/cnv/sliding_windows_sim_control.RData"
params.cnv_vcf_header = "${params.ref_dir}/cnv/cnv_vcf_header.tsv"

params.dbnsfp      = "${params.ref_dir}/dbnsfp/dbNSFP4.1a.txt.gz"
params.dbnsfp_tbi  = "${params.dbnsfp}.tbi"
params.dbnsfp_dt   = "${params.dbnsfp}.data_types"

params.outdir      = "${params.bucket}/Pipeline_Output"


// Will be using dudeML and DeepVariant in the future

if (params.single_lane == "YES"){

    params.match = ""

    params.reads_path = "${params.bucket}/Fastqs/${params.run_id}*${params.match}"

    log.info """\

         W G S - U F L    P I P E L I N E
         ================================
         reference     : ${params.reference}
         reads         : ${params.reads_path}
         cnv control   : ${params.cnv_control}
         dbnsfp        : ${params.dbnsfp}
         outdir        : ${params.outdir}

         """
         .stripIndent()
    
    reads_ch1 = Channel.fromFilePairs(params.reads_path)
    reads_ch2 = Channel.fromFilePairs(params.reads_path)


    process fastqc_single {

        tag "${sample_id}"
        publishDir "${params.outdir}/${params.run_id}/${sample_id}/", mode: 'copy'
        label 'small_process'

        input:
        tuple sample_id, path(reads) from reads_ch1

        output:
        path "fastqc_${sample_id}_logs" into fastqc_ch

        script:
        """
        mkdir fastqc_${sample_id}_logs
        fastqc \
        -t ${task.cpus} \
        -o fastqc_${sample_id}_logs \
        -f fastq \
        ${reads[0]} \
        ${reads[1]}
        """
    }

    process trimReads_single {

        tag "${sample_id}"
        publishDir "${params.outdir}/${params.run_id}/${sample_id}/Trimmomatic", mode: 'copy'
        label 'small_process'

        input:
        tuple sample_id, path(reads) from reads_ch2

        output:
        tuple sample_id, file("${sample_id}_forward-paired.fastq.gz"), file("${sample_id}_reverse-paired.fastq.gz"), file("${sample_id}_forward-unpaired.fastq.gz"), file("${sample_id}_reverse-unpaired.fastq.gz") into trimmed_ch
        tuple sample_id, file("${sample_id}_trim_out.log") into trim_log_ch

        script:
        """
        TrimmomaticPE -threads ${task.cpus} \
        ${reads[0]} \
        ${reads[1]} \
        ${sample_id}_forward-paired.fastq.gz \
        ${sample_id}_forward-unpaired.fastq.gz \
        ${sample_id}_reverse-paired.fastq.gz \
        ${sample_id}_reverse-unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 MINLEN:36 2> ${sample_id}_trim_out.log
        """
    }
}
else {

    if (params.exome == "YES"){
        params.reads1 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_{L001,L002}_R1_001.fastq.gz"
        params.reads2 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_{L001,L002}_R2_001.fastq.gz"
    }
    else {
        params.reads1 = "${params.bucket}/Fastqs/${params.run_id}*_{L001,L002}_R1_001.fastq.gz"
        params.reads2 = "${params.bucket}/Fastqs/${params.run_id}*_{L001,L002}_R2_001.fastq.gz"
    }

    log.info """\

         W G S - U F L    P I P E L I N E
         ================================
         reference     : ${params.reference}
         reads R1      : ${params.reads1}
         reads R2      : ${params.reads2}
         cnv control   : ${params.cnv_control}
         dbnsfp        : ${params.dbnsfp}
         outdir        : ${params.outdir}

         """
         .stripIndent()

    reads1_ch = Channel.fromFilePairs(params.reads1)
    reads2_ch = Channel.fromFilePairs(params.reads2)

    process catLanes {

        tag "${sample_id}"
        label 'small_process'

        input:
        tuple sample_id, path(reads1) from reads1_ch
        tuple sample_id, path(reads2) from reads2_ch 

        output:
        tuple sample_id, file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz") into (read_pairs_ch1, read_pairs_ch2)

        script:
        """
        cat ${reads1[0]} ${reads1[1]} > ${sample_id}_R1.fastq.gz
        cat ${reads2[0]} ${reads2[1]} > ${sample_id}_R2.fastq.gz
        """
    }

    process fastqc {

        tag "${sample_id}"
        publishDir "${params.outdir}/${params.run_id}/${sample_id}/", mode: 'copy'
        label 'small_process'

        input:
        tuple sample_id, file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz") from read_pairs_ch1

        output:
        path "fastqc_${sample_id}_logs" into fastqc_ch

        script:
        """
        mkdir fastqc_${sample_id}_logs
        fastqc \
        -t ${task.cpus} \
        -o fastqc_${sample_id}_logs \
        -f fastq \
        ${sample_id}_R1.fastq.gz \
        ${sample_id}_R2.fastq.gz
        """
    }

    process trimReads {

        tag "${sample_id}"
        publishDir "${params.outdir}/${params.run_id}/${sample_id}/Trimmomatic", mode: 'copy'
        label 'small_process'

        input:
        tuple sample_id, file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz") from read_pairs_ch2

        output:
        tuple sample_id, file("${sample_id}_forward-paired.fastq.gz"), file("${sample_id}_reverse-paired.fastq.gz"), file("${sample_id}_forward-unpaired.fastq.gz"), file("${sample_id}_reverse-unpaired.fastq.gz") into trimmed_ch
        tuple sample_id, file("${sample_id}_trim_out.log") into trim_log_ch

        script:
        """
        TrimmomaticPE -threads ${task.cpus} \
        ${sample_id}_R1.fastq.gz \
        ${sample_id}_R2.fastq.gz \
        ${sample_id}_forward-paired.fastq.gz \
        ${sample_id}_forward-unpaired.fastq.gz \
        ${sample_id}_reverse-paired.fastq.gz \
        ${sample_id}_reverse-unpaired.fastq.gz \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
        LEADING:3 TRAILING:3 MINLEN:36 2> ${sample_id}_trim_out.log
        """
    }
}

process alignTrimmedReads {
    
    tag "${sample_id}"
    label 'alignment'

    input:
    path reference from params.reference
    path bwa_amb from params.bwa_amb
    path bwa_ann from params.bwa_ann
    path bwa_bwt from params.bwa_bwt
    path bwa_pac from params.bwa_pac
    path bwa_sa from params.bwa_sa
    tuple sample_id, file("${sample_id}_forward-paired.fastq.gz"), file("${sample_id}_reverse-paired.fastq.gz"), file("${sample_id}_forward-unpaired.fastq.gz"), file("${sample_id}_reverse-unpaired.fastq.gz") from trimmed_ch

    output:
    tuple sample_id, file("${sample_id}.sam") into aligned_ch

    script:
    """
    bwa mem -t ${task.cpus} \
    ${reference} \
    ${sample_id}_forward-paired.fastq.gz \
    ${sample_id}_reverse-paired.fastq.gz > \
    ${sample_id}.sam
    """
}

process samToBam {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/alignment", mode: 'copy'
    label 'high_mem'

    input:
    tuple sample_id, file("${sample_id}.sam") from aligned_ch

    output:
    tuple sample_id, file("${sample_id}-sort.bam"), file("${sample_id}-sort.bam.bai") into (bam_ch1, bam_ch2, bam_ch3, bam_ch4)

    script:
    """
    samtools view -@ ${task.cpus} -Sbu -o ${sample_id}.bam ${sample_id}.sam
    samtools sort -@ ${task.cpus} -o ${sample_id}-sort.bam ${sample_id}.bam
    samtools index -@ ${task.cpus} ${sample_id}-sort.bam
    """
}

process collectWgsMetrics {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/wgs_metrics", mode: 'copy'
    label 'high_mem'

    input:
    path reference from params.reference
    tuple sample_id, file("${sample_id}-sort.bam"), file("${sample_id}-sort.bam.bai") from bam_ch1

    output:
    file "${sample_id}_gatk_collect_wgs_metrics.txt" into wgs_metrics_ch

    script:
    """
    java -jar -XX:ParallelGCThreads=${task.cpus} -Xmx32g /picard.jar CollectWgsMetrics \
    -I ${sample_id}-sort.bam \
    -O ${sample_id}_gatk_collect_wgs_metrics.txt \
    -R ${reference}
    """
}

process callSNV {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'medium_process'

    input:
    path reference from params.reference
    path ref_fai from params.ref_fai
    tuple sample_id, file("${sample_id}-sort.bam"), file("${sample_id}-sort.bam.bai") from bam_ch2

    output:
    path "${sample_id}_strelka2/results/variants/variants.vcf.gz" into snv_ch
    
    script:
    """
    mkdir ${sample_id}_strelka2

    /strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
    --bam ${sample_id}-sort.bam \
    --referenceFasta ${reference} \
    --runDir ${sample_id}_strelka2

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g ${sample_id}_strelka2/runWorkflow.py

    ${sample_id}_strelka2/runWorkflow.py -m local -j ${task.cpus}
    """
}

process callCNV {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'medium_process'

    input:
    path control from params.cnv_control
    path header from params.cnv_vcf_header
    tuple sample_id, file("${sample_id}-sort.bam"), file("${sample_id}-sort.bam.bai") from bam_ch3

    output:
    tuple sample_id, file("${sample_id}_filtered_cnv.vcf") into cnv_ch
    file("${sample_id}_cnv.pdf")

    script:
    """
    callCNV.R ${sample_id} ${control}
    csvToVCF.sh ${sample_id} ${header}
    toGRanges.sh ${sample_id}
    CNVPlot.R ${sample_id}
    """
}

process expansionHunter {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/ExpansionHunter", mode: 'copy'
    label 'small_process'

    input:
    path reference from params.reference
    path ref_fai from params.ref_fai
    tuple sample_id, file("${sample_id}-sort.bam"), file("${sample_id}-sort.bam.bai") from bam_ch4

    output:
    tuple sample_id, file("${sample_id}_eh.vcf") into exp_hunt_ch
    file "${sample_id}_eh_realigned.bam" into exp_hunt_json_ch
    file "${sample_id}_eh.json"

    script:
    """
    /ExpansionHunter-v4.0.2-linux_x86_64/bin/ExpansionHunter \
    --reads ${sample_id}-sort.bam \
    --reference ${reference} \
    --variant-catalog /ExpansionHunter-v4.0.2-linux_x86_64/variant_catalog/hg19/variant_catalog.json \
    --output-prefix ${sample_id}_eh
    """
}

process mergeVCF {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'small_process'

    input:
    tuple sample_id, file("${sample_id}_filtered_cnv.vcf.gz"), path("${sample_id}_strelka2/results/variants/variants.vcf.gz") from vcf_ch

    output:
    tuple sample_id, file("${sample_id}_concat.vcf.gz"), file("${sample_id}_concat.vcf.gz.csi") into merged_ch

    script:
    """
    # change SAMPLE1 to ${sample_id}-sort in Strelka2 VCF file

    bcftools index --threads ${task.cpus} \
    ${sample_id}_strelka2/results/variants/variants.vcf.gz

    bcftools index --threads ${task.cpus} \
    ${sample_id}_filtered_cnv.vcf.gz

    bcftools concat --threads ${task.cpus} -a \
    -o ${sample_id}_concat.vcf \
    ${sample_id}_filtered_cnv.vcf.gz \
    ${sample_id}_strelka2/results/variants/variants.vcf.gz

    bgzip -@ ${task.cpus} ${sample_id}_concat.vcf

    bcftools index --threads ${task.cpus} ${sample_id}_concat.vcf.gz
    """
}

process annotateVCF {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'high_mem'

    input:
    path dbNSFP from params.dbnsfp
    path dbNSFP_tbi from params.dbnsfp_tbi
    path dbNSFP_data_types from params.dbnsfp_dt
    tuple sample_id, file("${sample_id}_concat.vcf.gz"), file("${sample_id}_concat.vcf.gz.csi") from merged_ch
    tuple sample_id, file("${sample_id}_eh.vcf.gz") from exp_hunt_ch

    output:
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz"), file("${sample_id}_eh_snpsift.vcf.gz") into ann_ch
    tuple sample_id, file("${sample_id}_snpeff_stats.csv"), file("${sample_id}_eh_snpeff_stats.csv") into snpeff_stats_ch

    script:
    """
    java -jar -XX:ParallelGCThreads=${task.cpus} -Xmx32g /snpEff/snpEff.jar -csvStats ${sample_id}_snpeff_stats.csv -v -canon hg19 ${sample_id}_concat.vcf.gz > ${sample_id}_concat_snpeff.vcf

    bgzip -@ ${task.cpus} ${sample_id}_concat_snpeff.vcf

    java -jar -XX:ParallelGCThreads=${task.cpus} -Xmx32g /snpEff/SnpSift.jar dbnsfp -v -db ${dbNSFP} ${sample_id}_concat_snpeff.vcf.gz > ${sample_id}_concat_snpsift.vcf

    bgzip -@ ${task.cpus} ${sample_id}_concat_snpsift.vcf

    java -jar -XX:ParallelGCThreads=${task.cpus} -Xmx32g /snpEff/snpEff.jar -csvStats ${sample_id}_eh_snpeff_stats.csv -v -canon hg19 ${sample_id}_eh.vcf.gz > ${sample_id}_eh_snpeff.vcf

    bgzip -@ ${task.cpus} ${sample_id}_eh_snpeff.vcf

    java -jar -XX:ParallelGCThreads=${task.cpus} -Xmx32g /snpEff/SnpSift.jar dbnsfp -v -db ${dbNSFP} ${sample_id}_eh_snpeff.vcf.gz > ${sample_id}_eh_snpsift.vcf

    bgzip -@ ${task.cpus} ${sample_id}_eh_snpsift.vcf
    """
}