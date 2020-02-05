#!/bin/bash

set -e

threads=1
mem=25000
pmem=$(($mem / $threads))

THRESH=0

basedir=/cluster/work/grlab/projects/ICGC
datadir=${basedir}/alignments/alignments_all_STAR
outdir=${basedir}/junction_db/projected_alignments
junction_map=${basedir}/junction_db/icgc_junction_db.t${THRESH}.junction_map.pickle
mkdir -p ${outdir}

for fname in $(ls -1 ${datadir}/*.hdf5)
do
    fbase=$(basename $fname)
    outfname=${outdir}/${fbase%aligned.hdf5}t${THRESH}.projected.hdf5
    if [ ! -f ${outfname} ]
    then
        echo "source activate python36; python $(pwd)/collect_junction_db_step2_counts.py $junction_map $fname ${outfname}" | bsub -M ${mem} -J compr_gtex -We 2:00 -R "rusage[mem=${pmem}]" -R "span[hosts=1]" -n $threads -o /dev/null
    fi
done
