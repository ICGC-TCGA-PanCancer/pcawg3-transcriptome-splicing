#!/bin/bash

set -e

basedir=/cluster/project/grlab/projects/ICGC
anno=${basedir}/annotation/gencode.v19.annotation.hs37d5_chr.spladder.gtf
confidence=2

outdir=${basedir}/alternative_splicing_gtex
mkdir -p $outdir

if [ -z "$1" ]
then
    echo "Usage: $0 <event_type> [<threads>]"
    echo "  Where <event_type> is one of exon_skip, intron_retention,"
    echo "  alt_3prime, alt_5prime, mutex_exons, mult_exon_skip or all"
    exit 1
fi
event_type=$1
if [ "$event_type" == "all" ]
then
    event_type="alt_3prime,alt_5prime,exon_skip,intron_retention,mutex_exons,mult_exon_skip"
fi
if [ ! -z "$2" ]
then
    threads=$2
else
    threads=1
fi


if [ ! -f "$outdir/sample_list.txt" ]
then
    find /cluster/project/grlab/projects/ICGC/alignments_gtex_2015-03-29 -name \*.aligned.done | sed -e "s/done$/bam/g" > ${outdir}/sample_list.txt
fi

echo $anno
echo $outdir

### run spladder 
if [ "$3" == "local" ]
then
    python /cluster/project/grlab/home/akahles/git/software/spladder/python/spladder.py -b ${outdir}/sample_list.txt -o $outdir -a ${anno} -v y -c ${confidence}  -M merge_graphs -T y -V n -n 50 -P y -p y -t ${event_type} --parallel $threads
else
    echo "python /cluster/project/grlab/home/akahles/git/software/spladder/python/spladder.py -b ${outdir}/sample_list.txt -o $outdir -a ${anno} -v y -c ${confidence}  -M merge_graphs -T y -V n -n 50 -P y -p y -t ${event_type} --parallel $threads" | bsub -M 80000 -J icgc_as -W 120:00 -R "rusage[mem=10000]" -R "span[hosts=1]" -n $threads -o "run_as_gtex_${event_type}.log"   
fi
