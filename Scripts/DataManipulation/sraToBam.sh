#!/bin/bash
# usage bash sraToBam.sh <1: ref> <2: trait> <3: accessionlist> <4: num_cores>

tmp=$(mktemp -d Data/Scratch/XXXXXX)

species="${3%.txt}"
species="${species##*/}"
while read accession; do
    fasterq-dump "$accession" -O "$tmp" -t "$tmp"/FasterqDump -e "$4"
    bowtie2 -x Data/Scratch/Alignments/"$1".index \
        "$tmp"/"$accession"_1.fastq \
        "$tmp"/"$accession"_2.fastq \
        | samtools view -b -S - \
        | samtools sort - Data/Scratch/Alignments/Out/"$2"_"$species"_"$accession".sorted.bam
    samtools index Data/Scratch/Alignments/Out/"$2"_$species"_$accession".sorted.bam
done <$3
