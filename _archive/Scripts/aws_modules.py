#!/usr/bin/env python3

import boto3
import multiqc
from os import system, listdir


def get_s3_resource():
    """
    """
    session = boto3.Session()
    s3 = session.resource('s3')
    return s3


def get_batch_resource():
    """
    """
    client = boto3.client('batch')
    return client


def get_samples(bucket, run_id, s3):
    """
    """
    my_bucket = s3.Bucket(bucket)
    sample_ids = []
    for objects in my_bucket.objects.filter(Prefix = f'Pipeline_Output/{run_id}/'):
        file_name = objects.key.split('/')[2]
        if file_name != 'MultiQC':
            sample_ids.append(file_name)
    return sorted(set(sample_ids))


def get_run_id(bucket, s3):
    """
    """
    my_bucket = s3.Bucket(bucket)
    run_ids = []
    for objects in my_bucket.objects.filter(Prefix = f'Pipeline_Output/'):
        run_id = objects.key.split('/')[1]
        if not run_id.startswith('_'):
            if run_id != '':
                run_ids.append(run_id)
    return sorted(set(run_ids))


def get_run_id_from_fastqs(bucket, s3, exome):
    """
    """
    my_bucket = s3.Bucket(bucket)
    if exome: prefix = 'Exome_Fastqs/'
    else: prefix = 'Fastqs/'
    run_ids = []
    for objects in my_bucket.objects.filter(Prefix = prefix):
        file_name = objects.key.split('/')
        if not file_name[1].startswith('_'):
            if file_name[1].split('-')[0] != '':
                run_ids.append(f"{file_name[1].split('-')[0]}-{file_name[1].split('-')[1]}-{file_name[1].split('-')[2]}")
    return sorted(set(run_ids))


def get_panels(bucket, s3):
    """
    """
    my_bucket = s3.Bucket(bucket)
    panels = []
    for objects in my_bucket.objects.filter(Prefix = 'Pipeline/Reference/panels/'):
        file_name = objects.key.split('/')
        if file_name[3] != '':
            panels.append(file_name[3])
    return panels

def get_qc_files(bucket, run_id, s3):
    """
    """
    system(f'mkdir {run_id}_multiqc_sample_data')
    my_bucket = s3.Bucket(bucket)
    for objects in my_bucket.objects.filter(Prefix = f'Pipeline_Output/{run_id}/'):
        qc_file = (
            '_md_metrics.txt' in objects.key
            or '_trim_out.log' in objects.key
            or '_snpeff_stats.csv' in objects.key
            or '_fastqc.html' in objects.key
            or '_fastqc.zip' in objects.key
            or '_wgs_metrics.txt' in objects.key
        )
        if qc_file:
            file_name = objects.key.split('/')[-1]
            my_bucket.download_file(objects.key, f'{run_id}_multiqc_sample_data/{file_name}')
            print(f'Downloaded: {file_name}')


def run_multiqc(run_id):
    """
    """
    multiqc.run(
        analysis_dir = f'{run_id}_multiqc_sample_data',
        filename = f'{run_id}_multiqc_report'
    )


def upload_multiqc(run_id, s3):
    """
    """
    s3.meta.client.upload_file(
        f'{run_id}_multiqc_report.html',
        'hakmonkey-genetics-lab',
        f'Pipeline_Output/{run_id}/MultiQC/{run_id}_multiqc_report.html'
    )
    for file in listdir(f'{run_id}_multiqc_report_data/'):
        s3.meta.client.upload_file(
            f'{run_id}_multiqc_report_data/{file}',
            'hakmonkey-genetics-lab',
            f'Pipeline_Output/{run_id}/MultiQC/{run_id}_multiqc_report_data/{file}'
        )
    system(f'rm -rf {run_id}_multiqc_sample_data {run_id}_multiqc_report_data {run_id}_multiqc_report.html')
    print(f'MultiQC data for {run_id} has been uploaded.')


def move_fastqs(bucket, run_id, s3):
    """
    """
    my_bucket = s3.Bucket(bucket)
    for objects in my_bucket.objects.filter(Prefix = f'Fastqs/{run_id}'):
        copy_source = {
            'Bucket': my_bucket.name,
            'Key': objects.key
        }
        fastq = objects.key[7:]
        print(f'Moving: {fastq}')
        s3.meta.client.copy(
            copy_source,
            my_bucket.name,
            f'Fastqs/_Processed/{run_id}/{fastq}',
            ExtraArgs = {
                'StorageClass': 'GLACIER_IR',
                'MetadataDirective': 'COPY'
            }
        )
        objects.delete()
        print(f'Finished moving: {fastq}')


def submit_bs_to_aws_job(client, run_id, bucket):
    client.submit_job(
        jobName=f'bs-to-aws_{run_id}',
        jobQueue='hakmonkey-bs-to-aws',
        jobDefinition='bs-to-aws-ufl-germline:1',
        containerOverrides={
            'command': [
                'bash',
                '-c',
                f'./bs-to-aws.sh -b {bucket} -r {run_id};'
            ]
        }
    )


def submit_nextflow_job(client, nextflow_command, pipeline, run_id, bucket):
    """
    """
    client.submit_job(
        jobName=f'UFL-{pipeline}_{run_id}',
        jobQueue='hakmonkey-nextflow',
        jobDefinition='nextflow-ufl-germline:10',
        containerOverrides={
            'command': [
                'bash',
                '-c',
                f'{nextflow_command}; aws s3 cp trace.txt s3://{bucket}/Pipeline_Output/_work/{run_id}/; aws s3 cp timeline.html s3://{bucket}/Pipeline_Output/_work/{run_id}/; aws s3 cp {run_id}_report.html s3://{bucket}/Pipeline_Output/_work/{run_id}/'
            ]
        }
    )



def main():
    """
    """
    #s3 = get_s3_resource()
    #client = get_batch_resource()
    #sample_ids = get_samples('hakmonkey-genetics-lab', 'LU-PD-01', s3)
    #f_run_ids = get_run_id_from_fastqs('hakmonkey-genetics-lab', s3, False)
    #e_run_ids = get_run_id_from_fastqs('hakmonkey-genetics-lab', s3, True)
    #get_qc_files('hakmonkey-genetics-lab', 'LU-PD-01', s3)
    #run_multiqc('LU-PD-01')
    #upload_multiqc('LU-PD-01', s3)
    #get_run_id('hakmonkey-genetics-lab', s3)
    #submit_nextflow_job(client, 'nextflow run /data/main.nf  -with-report LU-PD-01_report.html -work-dir s3://hakmonkey-genetics-lab/Pipeline_Output/_work/LU-PD-01 --bucket s3://hakmonkey-genetics-lab --run_id LU-PD-01 --match _{1,2}.fq.gz --one --germline --aws', 'germline', 'LU-PD-01')
    #panels = get_panels('hakmonkey-genetics-lab', s3)
    
    


if __name__ == '__main__':
    main()