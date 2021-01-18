import os
import csv

pairs = []

with open('barg.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        panels = row[1].split(",")
        for panel in panels:
            pairs.append([row[0], panel])

# Will need to reverse this after running the nextflow pipeline

os.system("sed -i 's/pairs_ch = Channel.from()/pairs_ch = Channel.from(%s)/' test-tuple.nf"%pairs)

cmd = "sed -ri \"s/\\[([a-zA-Z0-9_.-]+),\\s([a-zA-Z0-9_.-]+)\\]/\\['\\1','\\2'\\]/g\" test-tuple.nf"

os.system(cmd)

