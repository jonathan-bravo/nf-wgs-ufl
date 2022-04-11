#!/usr/bin/env python3

from tkinter import *
from tkinter.font import Font
from tkinter import ttk
import aws_modules as aws


def restart(root):
    root.destroy()
    root = Tk()
    GatorGenome(root)
    root.mainloop()


class GatorGenome:
    def __init__(self, root):
        self.root = root
        self.root.title('GatorGenome')
        self.root.geometry('640x480')
        self.font = Font(family = "Helvetica", size = 18)
        self.root.option_add("*TCombobox*Listbox*Font", self.font)
        self.container = Frame(root)
        self.container.pack(fill=BOTH, expand=1)
        self.canvas = Canvas(self.container, highlightthickness=0)
        self.canvas.pack(side=LEFT, fill=BOTH, expand=1)
        self.scrollbar = ttk.Scrollbar(self.container, orient=VERTICAL, command=self.canvas.yview)
        self.scrollbar.pack(side=RIGHT, fill=Y)
        self.scroll_frame = Frame(self.canvas)
        self.canvas.configure(yscrollcommand=self.scrollbar.set)
        self.scroll_frame.bind('<Configure>', lambda e: self.canvas.configure(scrollregion = self.canvas.bbox('all')))
        self.root.bind("<MouseWheel>", lambda e: self.canvas.yview_scroll(-e.delta, 'units'))
        self.canvas.create_window((0,0), window=self.scroll_frame, anchor='nw', width=640)
        self.exec_menu()

    
    def run_multiqc_for_run(self):
        s3 = aws.get_s3_resource()
        aws.get_qc_files(self.bucket_name.get(), self.run_id_options.get(), s3)
        aws.run_multiqc(self.run_id_options.get())
        aws.upload_multiqc(self.run_id_options.get(), s3)

        restart(root)


    def run_move_fastqs_to_processed(self):
        s3 = aws.get_s3_resource()
        aws.move_fastqs(self.bucket_name.get(), self.run_id_options.get(), s3)


    def submit_bs_to_aws_bach_job(self):
        client = aws.get_batch_resource(client, self.bs_run_id.get(), self.bucket_name.get())
        aws.submit_bs_to_aws_job()

    
    def submit_germline_batch_job(self):

        nf_command = f'nextflow run /data/main.nf ' \
            f'-with-report {self.run_id_options.get()}_report.html '

        if self.exec_options.get() == 'AWS':
            nf_command += f'-work-dir ' \
                f's3://{self.bucket_name.get()}/Pipeline_Output/_work/{self.run_id_options.get()} ' \
                f'--bucket s3://{self.bucket_name.get()} --aws '

        if self.pipeline_options.get() == 'Germline':
            nf_command += '--germline '
            name = 'Germline'
        elif self.pipeline_options.get() == 'Concatinate Reads':
            nf_command += '--cat_reads '
            name = 'Cat_Reads'
        elif self.pipeline_options.get() == 'Alignment':
            nf_command += '--align '
            name = 'Align'
        elif self.pipeline_options.get() == 'Variant Calling':
            nf_command += '--vcf '
            name = 'VCF'
        elif self.pipeline_options.get() == 'MultiQC':
            nf_command += '--multiqc '
            name = 'MultiQC'
        elif self.pipeline_options.get() == 'CNV Controls':
            nf_command += '--cnv_controls '
            name = 'CNV_Controls'

        if self.exome_options.get() == 'WES':
            nf_command += '--exome '
        try:
            if self.lane_options.get() == 'NextSeq: [One]':
                nf_command += '--one '
                if self.match_options.get() == 'Illumina: [_{R1,R2}_001.fastq.gz]':
                    nf_command += '--match \"_{R1,R2}_001.fastq.gz\" '
                elif self.match_options.get() == 'General: [_{1,2}.fq.gz]':
                    nf_command += '--match \"_{1,2}.fq.gz\" '
            elif self.lane_options.get() == 'NovaSeq SP/S1/S2 Flowcell: [Two]':
                nf_command += '--two '
            elif self.lane_options.get() == 'NovaSeq S4 Flowcell: [Four]':
                nf_command += '--four '
        except: pass

        nf_command += f'--run_id {self.run_id_options.get()}'

        client = aws.get_batch_resource()
        aws.submit_nextflow_job(
            client,
            nf_command,
            name,
            self.run_id_options.get(),
            self.bucket_name.get()
        )

        print(nf_command)

        restart(root)


    def submit_reporting_batch_job(self):
        sample_panel_pairs = ''
        for sample in self.samples:
            selected_indices = globals()[f'choices_{sample}'].curselection()
            for i in range(len(selected_indices)):
                sample_panel_pairs += f"{sample},{globals()[f'choices_{sample}'].get(selected_indices[i])},"

        nf_command = f'nextflow run /data/main.nf ' \
            f'-with-report {self.run_id_options.get()}_report.html '

        if self.exec_options.get() == 'AWS':
            nf_command += f'-work-dir ' \
                f's3://{self.bucket_name.get()}/Pipeline_Output/_work/{self.run_id_options.get()} ' \
                f'--bucket s3://{self.bucket_name.get()} --aws '

        nf_command += f'--report --run_id {self.run_id_options.get()} --panel_reporting_list "{sample_panel_pairs[:-1]}"'
        
        client = aws.get_batch_resource()
        aws.submit_nextflow_job(
            client,
            nf_command,
            'Reporting',
            self.run_id_options.get(),
            self.bucket_name.get()
        )

        restart(root)
        


    def get_run_ids(self):
        s3 = aws.get_s3_resource()
        output_check = (
            self.pipeline_options.get() == 'Variant Calling'
            or self.pipeline_options.get() == 'MultiQC'
            or self.pipeline_options.get() == 'Reporting'
            or self.pipeline_options.get() == 'MultiQC Run'
        )
        if output_check:
            run_ids = aws.get_run_id(self.bucket_name.get(), s3)
        elif self.exome_options.get() == 'WGS':
            run_ids = aws.get_run_id_from_fastqs(self.bucket_name.get(), s3, False)
        elif self.exome_options.get() == 'WES':
            run_ids = aws.get_run_id_from_fastqs(self.bucket_name.get(), s3, True)
        return run_ids


    def get_reporting_info(self):
        s3 = aws.get_s3_resource()
        panels = aws.get_panels(self.bucket_name.get(), s3)
        samples = aws.get_samples(self.bucket_name.get(), self.run_id_options.get(), s3)
        return (panels, samples)


    def panel_multiselect(self, e):
        chosen = self.parent_choices.curselection()
        for sample in self.samples:
            sample_list = globals()[f'choices_{sample}']
            try:
                for i in self.to_remove:
                    sample_list.selection_clear(i)
            except:
                pass
            for i in chosen:
                sample_list.selection_set(i)
        self.to_remove = chosen


    def show_from_exec(self, e):
        try: self.bucket_frame.destroy()
        except AttributeError: pass
        try: self.pipeline_frame.destroy()
        except AttributeError: pass
        try: self.exome_frame.destroy()
        except AttributeError: pass
        try: self.run_id_frame.destroy()
        except AttributeError: pass
        try: self.reporting_frame.destroy()
        except AttributeError: pass
        try: self.lanes_frame.destroy()
        except AttributeError: pass
        try: self.match_frame.destroy()
        except AttributeError: pass
        try: self.launch_frame.destroy()
        except AttributeError: pass
        try: self.multiqc_run_frame.destroy()
        except AttributeError: pass
        try: self.move_fastqs_frame.destroy()
        except AttributeError: pass
        try: self.bs_to_aws_frame.destroy()
        except AttributeError: pass
        if self.exec_options.get() == 'AWS':
            self.bucket_menu()
            self.pipeline_menu()
        elif self.exec_options.get() == 'Local':
            self.pipeline_menu()
            
            
    def show_from_pipeline(self, e):
        try: self.exome_frame.destroy()
        except AttributeError: pass
        try: self.run_id_frame.destroy()
        except AttributeError: pass
        try: self.reporting_frame.destroy()
        except AttributeError: pass
        try: self.lanes_frame.destroy()
        except AttributeError: pass
        try: self.match_frame.destroy()
        except AttributeError: pass
        try: self.launch_frame.destroy()
        except AttributeError: pass
        try: self.multiqc_run_frame.destroy()
        except AttributeError: pass
        try: self.move_fastqs_frame.destroy()
        except AttributeError: pass
        try: self.bs_to_aws_frame.destroy()
        except AttributeError: pass
        exome_check = (
            self.pipeline_options.get() == 'Germline'
            or self.pipeline_options.get() == 'Concatinate Reads'
            or self.pipeline_options.get() == 'Alignment'
            or self.pipeline_options.get() == 'Variant Calling'
            or self.pipeline_options.get() == 'MultiQC'
            or self.pipeline_options.get() == 'MultiQC Run'
            or self.pipeline_options.get() == 'Archive Fastqs'
        )
        run_id_check = (
            self.pipeline_options.get() == 'Reporting'
        )
        if exome_check:
            self.exome_menu()
        elif run_id_check:
            self.run_id_menu()
        elif self.pipeline_options.get() == 'BS-to-AWS':
            self.bs_to_aws_menu()
        elif self.pipeline_options.get() == 'CNV Control Generation':
            self.launch_menu()

        
    def show_from_exome(self, e):
        try: self.run_id_frame.destroy()
        except AttributeError: pass
        try: self.reporting_frame.destroy()
        except AttributeError: pass
        try: self.lanes_frame.destroy()
        except AttributeError: pass
        try: self.match_frame.destroy()
        except AttributeError: pass
        try: self.launch_frame.destroy()
        except AttributeError: pass
        try: self.multiqc_run_frame.destroy()
        except AttributeError: pass
        try: self.move_fastqs_frame.destroy()
        except AttributeError: pass
        self.run_id_menu()


    def show_from_run_id(self, e):
        try: self.reporting_frame.destroy()
        except AttributeError: pass
        try: self.lanes_frame.destroy()
        except AttributeError: pass
        try: self.match_frame.destroy()
        except AttributeError: pass
        try: self.launch_frame.destroy()
        except AttributeError: pass
        try: self.multiqc_run_frame.destroy()
        except AttributeError: pass
        try: self.move_fastqs_frame.destroy()
        except AttributeError: pass
        if self.pipeline_options.get() == 'Reporting':
            self.reporting_menu()
        elif self.pipeline_options.get() == 'Variant Calling' or self.pipeline_options.get() == 'MultiQC':
            self.launch_menu()
        elif self.pipeline_options.get() == 'MultiQC Run':
            self.multiqc_run_menu()
        elif self.pipeline_options.get == 'Archive Fastqs':
            self.move_fastqs_menu()
        else:
            self.lanes_menu()


    def show_from_lanes(self, e):
        try: self.reporting_frame.destroy()
        except AttributeError: pass
        try: self.match_frame.destroy()
        except AttributeError: pass
        try: self.launch_frame.destroy()
        except AttributeError: pass
        if self.lane_options.get() == 'NextSeq: [One]':
            self.match_menu()
        else:
            self.launch_menu()


    def show_from_match(self, e):
        try: self.reporting_frame.destroy()
        except AttributeError: pass
        try: self.launch_frame.destroy()
        except AttributeError: pass
        self.launch_menu()


    ## EXECUTORS ##
    def exec_menu(self):
        self.exec_frame = Frame(self.scroll_frame)
        self.exec_frame.pack(fill=BOTH, expand=1)
        Label(
            self.exec_frame,
            text='Choose where to execute pipeline:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        executors = [
            'Local',
            'AWS'
        ]
        self.exec_options = ttk.Combobox(self.exec_frame, value=executors)
        self.exec_options.current(1)
        self.exec_options.pack(pady=20)
        self.exec_options.bind("<<ComboboxSelected>>", self.show_from_exec)

    ## BUCKET ##
    def bucket_menu(self):
        self.bucket_frame = Frame(self.scroll_frame)
        self.bucket_frame.pack(fill=BOTH, expand=1)
        Label(
            self.bucket_frame,
            text='Enter bucket name for pipeline run:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        self.bucket_name = StringVar()
        self.bucket_name.set('hakmonkey-genetics-lab')
        self.bucket_entry = Entry(self.bucket_frame, textvariable=self.bucket_name)
        self.bucket_entry.pack(pady=5)

    ## PIPELINES ##
    def pipeline_menu(self):
        self.pipeline_frame = Frame(self.scroll_frame)
        self.pipeline_frame.pack(fill=BOTH, expand=1)
        Label(
            self.pipeline_frame,
            text='Choose which workflow to launch:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        pipelines = [
            'Germline',
            'Concatinate Reads',
            'Alignment',
            'Variant Calling',
            'MultiQC',
            'Reporting',
            'CNV Control Generation',
            'MultiQC Run',
            'Archive Fastqs',
            'BS-to-AWS'
        ]
        self.pipeline_options = ttk.Combobox(self.pipeline_frame, value= pipelines)
        self.pipeline_options.current(0)
        self.pipeline_options.pack(pady=20)
        self.pipeline_options.bind("<<ComboboxSelected>>", self.show_from_pipeline)

    ## EXOME ##
    def exome_menu(self):
        self.exome_frame = Frame(self.scroll_frame)
        self.exome_frame.pack(fill=BOTH, expand=1)
        Label(
            self.exome_frame,
            text='Whole genome or exome pipeline:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        read_types = [
            'WGS',
            'WES'
        ]
        self.exome_options = ttk.Combobox(self.exome_frame, value=read_types)
        self.exome_options.current(0)
        self.exome_options.pack(pady=20)
        self.exome_options.bind("<<ComboboxSelected>>", self.show_from_exome)

    # ## RUN ID ##
    def run_id_menu(self):
        self.run_id_frame = Frame(self.scroll_frame)
        self.run_id_frame.pack(fill=BOTH, expand=1)
        Label(
            self.run_id_frame,
            text='Select your run Id:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        run_ids = self.get_run_ids()
        self.run_id_options = ttk.Combobox(self.run_id_frame, values=run_ids)
        self.run_id_options.current(0)
        self.run_id_options.pack(pady=20)
        self.run_id_options.bind("<<ComboboxSelected>>", self.show_from_run_id)

    # ## LANES ##
    def lanes_menu(self):
        self.lanes_frame = Frame(self.scroll_frame)
        self.lanes_frame.pack(fill=BOTH, expand=1)
        Label(
            self.lanes_frame,
            text='How many lanes is your run:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        lanes = [
            'NextSeq: [One]',
            'NovaSeq SP/S1/S2 Flowcell: [Two]',
            'NovaSeq S4 Flowcell: [Four]'
        ]
        self.lane_options = ttk.Combobox(self.lanes_frame, values=lanes)
        self.lane_options.current(0)
        self.lane_options.pack(pady=20)
        self.lane_options.bind("<<ComboboxSelected>>", self.show_from_lanes)

    # ## MATCH ##
    def match_menu(self):
        self.match_frame = Frame(self.scroll_frame)
        self.match_frame.pack(fill=BOTH, expand=1)
        Label(
            self.match_frame,
            text='What is the match pattern:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        matches = [
            'Illumina: [_{R1,R2}_001.fastq.gz]',
            'General: [_{1,2}.fq.gz]'
        ]
        self.match_options = ttk.Combobox(self.match_frame, values=matches)
        self.match_options.current(1)
        self.match_options.pack(pady=20)
        self.match_options.bind('<<ComboboxSelected>>', self.show_from_match)

    # ## LAUNCH ##
    def launch_menu(self):
        self.launch_frame = Frame(self.scroll_frame)
        self.launch_frame.pack(fill=BOTH, expand=1)
        Label(
            self.launch_frame,
            text='Launch pipeline:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        Button(
            self.launch_frame,
            text='LAUNCH',
            command=self.submit_germline_batch_job,
            font=("Helvetica", 18)
        ).pack(anchor="s", fill=BOTH, expand=1, pady=10)

    # ## REPORTING ##
    def reporting_menu(self):
        self.reporting_frame = Frame(self.scroll_frame)
        self.reporting_frame.pack(fill=BOTH, expand=1)
        Label(
            self.reporting_frame,
            text='Select panels and launch reporting:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        panels, self.samples = self.get_reporting_info()
        Label(
            self.reporting_frame,
            text=f'Parent Selector:',
            font=("Helvetica", 12)
        ).pack(anchor="n", padx=20)
        self.parent_choices = Listbox(self.reporting_frame, selectmode='multiple', exportselection=False)
        self.parent_choices.pack(padx=20, pady=10)
        for panel in panels:
            self.parent_choices.insert(END, panel)
        self.parent_choices.insert(END,'general')
        self.parent_choices.bind('<<ListboxSelect>>', self.panel_multiselect)
        for i, sample in enumerate(self.samples):
            Label(
                self.reporting_frame,
                text=f'{sample}:',
                font=("Helvetica", 12)
            ).pack(anchor="n", padx=20)
            globals()[f'choices_{sample}'] = Listbox(self.reporting_frame, selectmode='multiple', exportselection=False)
            globals()[f'choices_{sample}'].pack(padx=20, pady=10)
            for panel in panels:
                globals()[f'choices_{sample}'].insert(END, panel)
            globals()[f'choices_{sample}'].insert(END,'general')
        Button(
            self.reporting_frame,
            text='LAUNCH',
            command=self.submit_reporting_batch_job,
            font=("Helvetica", 18)
        ).pack(anchor="s", fill=BOTH, expand=1, pady=10)

    def multiqc_run_menu(self):
        self.multiqc_run_frame = Frame(self.scroll_frame)
        self.multiqc_run_frame.pack(fill=BOTH, expand=1)
        Label(
            self.multiqc_run_frame,
            text='Launch MultiQC for whole run:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        Button(
            self.multiqc_run_frame,
            text='LAUNCH',
            command=self.run_multiqc_for_run,
            font=("Helvetica", 18)
        ).pack(anchor="s", fill=BOTH, expand=1, pady=10)

    def move_fastqs_menu(self):
        self.move_fastqs_frame = Frame(self.scroll_frame)
        self.move_fastqs_frame.pack(fill=BOTH, expand=1)
        Label(
            self.move_fastqs_frame,
            text='Launch fastq(s) transfer to _Processed:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        Button(
            self.move_fastqs_frame,
            text='LAUNCH',
            command=self.run_move_fastqs_to_processed,
            font=("Helvetica", 18)
        ).pack(anchor="s", fill=BOTH, expand=1, pady=10)

    def bs_to_aws_menu(self):
        self.bs_to_aws_frame = Frame(self.scroll_frame)
        self.bs_to_aws_frame.pack(fill=BOTH, expand=1)
        Label(
            self.bs_to_aws_frame,
            text='Input run Id to transfer from BaseSpace to AWS:',
            font=("Helvetica", 24)
        ).pack(anchor="n", pady=10, padx=20)
        self.bs_run_id = StringVar()
        self.bs_run_id.set('NQ-00-00')
        self.bs_to_aws_entry = Entry(self.bs_to_aws_frame, textvariable=self.bs_run_id)
        self.bs_to_aws_entry.pack(pady=5)
        Button(
            self.bs_to_aws_frame,
            text='LAUNCH',
            command=self.submit_bs_to_aws_bach_job,
            font=("Helvetica", 18)
        ).pack(anchor="s", fill=BOTH, expand=1, pady=10)



root = Tk()
GatorGenome(root)
root.mainloop()


# ## Might also add the bs-to-aws workflow