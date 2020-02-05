#!/bin/bash

set -e

basedir=/cluster/project/grlab/projects/ICGC
anno=${basedir}/annotation/gencode.v19.annotation.hs37d5_chr.spladder.gtf

as_script=${HOME}/git/software/spladder/python/spladder.py

outdir=${basedir}/alternative_splicing_anno
mkdir -p $outdir

echo $anno
echo $outdir

bam=/cluster/project/grlab/projects/ICGC/alignments_ICGC_2015-04-25/012774c4-d4c0-11e4-93a5-c688eda0acae.aligned.bam
conf=2

### run spladder 
python ${as_script} --insert_ir=n --insert_es=n --insert_ni=n --remove_se=n --validate_sg=n -b $bam -o $outdir -a ${anno} -v y -c ${conf} -M merge_graphs -T y -n 60 -P y -p n --sparse_bam n -t exon_skip,intron_retention,alt_3prime,alt_5prime,mutex_exons,mult_exon_skip --parallel 1
