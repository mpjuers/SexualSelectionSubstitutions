#!/bin/bash
# usage: sh alignSeqsOfInterest.sh <ref.fasta> <input.fasta> <output.fasta>

bwa index ${1}
bwa mem <(gunzip -c ${1}) ${2} | \
    samtools view -S -b - | \
    samtools bam2fq - | \
    seqtk seq -A - > ${3}
