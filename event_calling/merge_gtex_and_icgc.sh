#!/bin/bash

set -e

basedir=/cluster/project/grlab/projects/ICGC

if [ -z "$1" ]
then
    echo "Usage: $0 <event_type>"
    exit 1
fi
eventtype=$1

python spladder_merge.py -v -i ${basedir}/alternative_splicing,${basedir}/alternative_splicing_gtex_count_only -o ${basedir}/alternative_splicing_gtex_icgc_merged -t $eventtype
