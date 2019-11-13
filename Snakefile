# Number of SNPs for each trait dataset
# n_snps={"matingBehavior": 2_400_000, "aggression": 1_914_528}

import itertools
from os import listdir
import os.path as path

traits = config["traits"].keys()
accession_files = [f for f in listdir("Data/Identifiers/") if path.isfile(path.join("Data/Identifiers", f))]
species = [os.path.basename(f).rstrip(".txt") for f in accession_files]
accessions = []
traits_out = []
species_out = []
for trait in traits:
    for file in accession_files:
        with open(file, 'r') as filestream:
            for line in filestream:
                accessions.append(line)
                species_out.append(os.path.basename(file).rstrip.(".txt"))
                traits_out.append(trait)


localrules: all, fdr, clean, scratchsetup, dirsetup
configfile: "snakemakeConfig.json"


rule all:
    input:
        expand("Data/Interest/{trait}.interest.txt", trait=traits) ,
        expand("Data/Scratch/InterestSeqs/{trait}.interest.nsnps", trait=traits) ,
        expand("Data/Scratch/Alignments/{trait}.interest.aligned.fasta", trait=traits) ,
        expand("Data/Scratch/Alignments/{trait}.index.1.bt2", trait=traits) ,
        expand(
            "Data/Scratch/Alignments/Out/{trait}_{species}_{accession}.bam",
            zip, trait=traits_out, species=species_out, accession=accessions
            ),
        "Data/Scratch/Genomes/dMelRefSeq.fna.gz", 
        "Data/Scratch/Genomes/dMelRefSeq.fna"


rule fdr:
    input:
        data = "Data/SNPData/{trait}.csv"
    output:
        "Data/Interest/{trait}.interest.txt"
    params:
        n_snps = lambda wildcards: config["traits"][wildcards.trait]["n_snps"]
    shell:
        (
         "python Scripts/DataManipulation/fdrCorrection.py"
         " {params.n_snps} {input.data} {output}"
        )


rule get_seqs:
    input:
        snps = "Data/Interest/{trait}.interest.txt"
    output:
        seqs = "Data/Scratch/InterestSeqs/{trait}.interest.fasta",
        nsnps = "Data/Scratch/InterestSeqs/{trait}.interest.nsnps"
    shell:
         (
         "wc -l {input.snps} > {output.nsnps}"
          + " && python Scripts/GetData/windows.py {input.snps} {output.seqs} "
          + config["email"]
          )


rule unzipref:
    input: "Data/Scratch/Genomes/dMelRefSeq.fna.gz"
    output: "Data/Scratch/Genomes/dMelRefSeq.fna"
    shell: 
        (
         "gzip -dc {input} > {output} &&"
         " bwa index {output} &&"
         " samtools faidx {output}"
        )


rule align_interest:
    input:
        query = "Data/Scratch/InterestSeqs/{trait}.interest.fasta",
        ref = "Data/Scratch/Genomes/dMelRefSeq.fna"
    output:
        align = "Data/Scratch/Alignments/{trait}.interest.aligned.fasta"
    shadow:
        "full"
    shell:
        (
         "bash Scripts/DataManipulation/alignSeqsOfInterest.sh"
         " {input.ref} {input.query} {output.align}"
        )


rule get_dmel_genome:
    output:
        "Data/Scratch/Genomes/dMelRefSeq.fna.gz"
    script:
        "Scripts/GetData/dMelGenome.py"


rule dirsetup:
    output:
        "Logs/Cluster",
        "Data/Scratch/tmp"
    input:
        "Data/Scratch"
    shell:
        "bash Scripts/Setup/dirSetup.sh {output}"


rule scratchsetup:
    output:
        "Data/Scratch"
    shell:
        "if [[ ! -h {output} ]]; then ln -s " + config["scratchdir"]  + "{output}; fi"


rule bowtie2index:
    input:
        "Data/Scratch/Alignments/{trait}.interest.aligned.fasta"
    output:
        "Data/Scratch/Alignments/{trait}.index.1.bt2",
        "Data/Scratch/Alignments/{trait}.index.2.bt2",
        "Data/Scratch/Alignments/{trait}.index.3.bt2",
        "Data/Scratch/Alignments/{trait}.index.4.bt2",
        "Data/Scratch/Alignments/{trait}.index.rev.1.bt2",
        "Data/Scratch/Alignments/{trait}.index.rev.2.bt2"
    shell:
        "bowtie2-build {input} {trait}.index"


rule align_bowtie:
    input:
        "Data/Scratch/Alignments/{trait}.index.1.bt2",
        "Data/Scratch/Alignments/{trait}.index.2.bt2",
        "Data/Scratch/Alignments/{trait}.index.3.bt2",
        "Data/Scratch/Alignments/{trait}.index.4.bt2",
        "Data/Scratch/Alignments/{trait}.index.rev.1.bt2",
        "Data/Scratch/Alignments/{trait}.index.rev.2.bt2"
    output:
        "Data/Scratch/Alignments/Out/{trait}_{species}_{accession}.sorted.bam"
    shell:
        (
         "bash Scripts/DataManipulation/sraToBam.sh"
         + " Data/Scratch/Alignments/{trait}.fasta"
         + " {trait}"
         + " Data/Identifiers/{accession}"
         + " " + config["rules"]["align_bowtie"]["cores"]
        )


# Removes everything except initial dependencies.
rule clean:
    shell:
        (
         "rm -r"
         " Logs"
         " Data/Interest"
         " Data/InterestSeqs"
         " Data/Scratch"
        )
