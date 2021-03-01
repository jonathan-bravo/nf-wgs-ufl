```bash
data_files=(`ls /data/`)
#${data_files[*]}

source_files=(`ls /source/`)
#${source_files[*]}

for ((i=0; i<"${#data_files[@]}"; i++)); do hap.py ${data_files[$i]} ${source_files[$i]} -f /bed/confidence.bed -o /data/${data_files[$i]}_vs_${source_files[$i]}.txt -r /reference/hg19.fa; done


docker run --rm -dt \
-v ${reference}:/reference/ \
-v ${bed}:/bed/ \
-v ${source}:/source/ \
-v ${data}:/data/ \
--name hap.py_run \
-c "" pkrusche/hap.py:latest
```


