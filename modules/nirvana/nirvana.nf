#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process NIRVANA {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/Nirvana", mode: 'copy'
    label 'nirvana'
    label 'medium_process'

    input:
    

    output:
    
    
    shell:
    '''
    ## installing mopst recent annotation data
    dotnet Nirvana/Downloader.dll --ga GRCh37 -o /Nirvana/Data

    ## running NIRVANA
    dotnet Nirvana/Nirvana.dll \
    -c Data/Cache/GRCh37/Both \
    --sd Data/SupplementaryAnnotation/GRCh37 \
    -r Data/References/Homo_sapiens.GRCh37.Nirvana.dat \
    -i !{sample_id}_variants.vcf.gz \
    -o !{sample_id}_Strelka2

    dotnet Nirvana/Nirvana.dll \
    -c Data/Cache/GRCh37/Both \
    --sd Data/SupplementaryAnnotation/GRCh37 \
    -r Data/References/Homo_sapiens.GRCh37.Nirvana.dat \
    -i !{sample_id}_filtered_cnv.vcf.gz \
    -o !{sample_id}_cn.MOPS

    dotnet Nirvana/Nirvana.dll \
    -c Data/Cache/GRCh37/Both \
    --sd Data/SupplementaryAnnotation/GRCh37 \
    -r Data/References/Homo_sapiens.GRCh37.Nirvana.dat \
    -i !{sample_id}_filtered_eh.vcf.gz \
    -o !{sample_id}_ExpansionHunter

    dotnet Nirvana/Nirvana.dll \
    -c Data/Cache/GRCh37/Both \
    --sd Data/SupplementaryAnnotation/GRCh37 \
    -r Data/References/Homo_sapiens.GRCh37.Nirvana.dat \
    -i !{sample_id}_diploidSV.vcf.gz \
    -o !{sample_id}_Manta

    ## splitting the Strelka2 JSON into managable pieces
    contigs=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y)

    for contig in "${contigs[@]}";
    do
        dotnet Nirvana/bin/Release/netcoreapp3.1/Jasix.dll -i !{sample_id}_Strelka2.json.gz -q ${contig} -o !{sample_id}_chr${contig}_Strelka2.json.gz;
    done

    for json in *_Strelka2.json.gz; do dotnet Nirvana/bin/Release/netcoreapp3.1/Jasix.dll -c -i ${json}; done
    '''
}