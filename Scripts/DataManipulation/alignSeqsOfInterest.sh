#!/bin/bash
# usage: sh alignSeqsOfInterest.sh <ref.fasta> <input.fasta> <output.fasta>

if [[ ! -d ${3%/*} ]]; then
    mkdir ${3%/*} ]]
fi
bwa mem $1 $2 | \
    samtools view -S -b - | \
    samtools bam2fq - | \
    seqtk seq -a - >$3
