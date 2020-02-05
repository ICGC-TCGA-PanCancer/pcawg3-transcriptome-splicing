#!/bin/bash

set -e

threads=2
mem=25000
pmem=$(($mem / $threads))

datadir=/cluster/work/grlab/projects/ICGC/alignments/alignments_all_STAR

for fname in $(ls -1 ${datadir}/*.bam)
do
    if [ ! -f ${fname%bam}hdf5 ]
    then
        if [ "$1" == "logs" ]
        then
            echo "source activate python36; python $(pwd)/compress_bam.py -b $fname -p ${threads} -v y" | bsub -M ${mem} -J compr_gtex -We 2:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -n $threads -o ${fname%bam}cluster.log
        else
            echo "source activate python36; python $(pwd)/compress_bam.py -b $fname -p ${threads}" | bsub -M ${mem} -J compr_gtex -We 2:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -n $threads -o /dev/null 
        fi
    fi
done
