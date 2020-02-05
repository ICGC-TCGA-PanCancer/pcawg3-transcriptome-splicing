#!/bin/bash

set -e

basedir=/cluster/project/grlab/projects/ICGC

for study in $(tail -n +2 ${basedir}/orig_data/metadata/rnaseq_metadata.tsv | cut -f 2 | sort -u)
do
    filelist=${basedir}/orig_data/metadata/alignments_${study}.txt
    echo "python /cluster/project/grlab/home/akahles/git/software/spladder/python/spladder_stats.py -o ${basedir}/alternative_splicing -s $filelist -l $study -v" #| bsub -M 51200 -J plt -W 2:00 -n 1 -R "rusage[mem=51200]" -R "span[hosts=1]"
    exit
done
