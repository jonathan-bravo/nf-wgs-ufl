params.panels            = "s3://hakmonkey-genetics-lab/Pipeline/Reference/panels"
params.ataxia            = "${params.panels}/ataxia"
params.dementia          = "${params.panels}/dementia"
params.dystonia          = "${params.panels}/dystonia"
params.epilepsy          = "${params.panels}/epilepsy"
params.hsp               = "${params.panels}/hsp"
params.neuromuscular     = "${params.panels}/neuromuscular"
params.neuropathy        = "${params.panels}/neuropathy"
params.parkinsons        = "${params.panels}/parkinsons"
params.sma               = "${params.panels}/sma"

process panelAtaxia {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'
    label 'panel'

    input:
    path panel from params.ataxia
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz") from ann_ch1

    output:
    tuple sample_id, file("${sample_id}_ataxia.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel} | tr '\n' '|')

    zgrep '#' !{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_ataxia.vcf

    zcat !{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_ataxia.vcf
    '''
}

process panelDementia {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'
    label 'panel'

    input:
    path panel from params.dementia
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz") from ann_ch2

    output:
    tuple sample_id, file("${sample_id}_dementia.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel} | tr '\n' '|')

    zgrep '#' !{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_dementia.vcf

    zcat !{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_dementia.vcf
    '''
}

process panelDystonia {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'
    label 'panel'

    input:
    path panel from params.dystonia
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz") from ann_ch3

    output:
    tuple sample_id, file("${sample_id}_dystonia.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel} | tr '\n' '|')

    zgrep '#' !{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_dystonia.vcf

    zcat !{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_dystonia.vcf
    '''
}

process panelEpilepsy {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'
    label 'panel'

    input:
    path panel from params.epilepsy
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz") from ann_ch4

    output:
    tuple sample_id, file("${sample_id}_epilepsy.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel} | tr '\n' '|')

    zgrep '#' !{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_epilepsy.vcf

    zcat !{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_epilepsy.vcf
    '''
}

process panelHsp {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'
    label 'panel'

    input:
    path panel from params.hsp
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz") from ann_ch5

    output:
    tuple sample_id, file("${sample_id}_hsp.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel} | tr '\n' '|')

    zgrep '#' !{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_hsp.vcf

    zcat !{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_hsp.vcf
    '''
}

process panelNeuromuscular {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'
    label 'panel'

    input:
    path panel from params.neuromuscular
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz") from ann_ch6

    output:
    tuple sample_id, file("${sample_id}_neuromuscular.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel} | tr '\n' '|')

    zgrep '#' !{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_neuromuscular.vcf

    zcat !{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_neuromuscular.vcf
    '''
}

process panelNeuropathy {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'
    label 'panel'

    input:
    path panel from params.neuropathy
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz") from ann_ch7

    output:
    tuple sample_id, file("${sample_id}_neuropathy.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel} | tr '\n' '|')

    zgrep '#' !{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_neuropathy.vcf

    zcat !{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_neuropathy.vcf
    '''
}

process panelParkinsons {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'
    label 'panel'

    input:
    path panel from params.parkinsons
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz") from ann_ch8

    output:
    tuple sample_id, file("${sample_id}_parkinsons.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel} | tr '\n' '|')

    zgrep '#' !{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_parkinsons.vcf

    zcat !{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_parkinsons.vcf
    '''
}

process panelSma {

    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/variants", mode: 'copy'
    label 'panel'

    input:
    path panel from params.sma
    tuple sample_id, file("${sample_id}_concat_snpsift.vcf.gz") from ann_ch9

    output:
    tuple sample_id, file("${sample_id}_sma.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel} | tr '\n' '|')

    zgrep '#' !{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_sma.vcf

    zcat !{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_sma.vcf
    '''
}