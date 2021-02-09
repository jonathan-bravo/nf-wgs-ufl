#!/usr/bin/env python3

import re
import sys

sample_id = str(sys.argv[1])
panel = str(sys.argv[2])

with open('{}_{}.info.tsv'.format(sample_id, panel), 'w') as f:
    with open('{}_{}.info.txt'.format(sample_id, panel), 'r') as fp:
        for cnt, line in enumerate(fp):
            if cnt == 0:
                header = sorted(line.split('\t'))
                for i, element in enumerate(header):
                    if "\n" in element:
                        header[i] = element[:-2]
                for i in header:
                    f.write('{}\t'.format(i))
                f.write('\n')
            else:
                sorted_list = ['.'] * len(header)
                sorted_values = sorted(line.split('\t'))
                for i in sorted_values:
                    x = re.match("([0-9a-zA-Z.-_]+=)",i)
                    for j, val in enumerate(header):
                        if val in x.groups()[0]:
                            sorted_list[j] = i
                for i, element in enumerate(sorted_list):
                    if "\n" in element:
                        sorted_list[i] = element[:-2]
                for i in sorted_list:
                    f.write('{}\t'.format(i))
                f.write('\n')