#!/bin/bash

set -e

basedir=/cluster/project/grlab/projects/ICGC
anno=${basedir}/annotation/gencode.v19.annotation.hs37d5_chr.spladder.gtf
#anno=${basedir}/annotation/gencode.v19.annotation.hs37d5_chr.ENSG00000002079.gtf

outdir=${basedir}/alternative_splicing_gtex_count_only
mkdir -p $outdir

if [ ! -f "$outdir/sample_list.txt" ]
then
    find /cluster/project/grlab/projects/ICGC/alignments_gtex_2015-03-29 -name \*.aligned.done | sed -e "s/done$/bam/g" > ${outdir}/sample_list.txt
fi

echo $anno
echo $outdir

### run spladder 
python /cluster/project/grlab/home/akahles/git/software/spladder/python/spladder.py -b ${outdir}/sample_list.txt -o $outdir -a ${anno} -v y -c 3 -M merge_graphs -T y -V n -n 50 -P y -p y -t exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip
