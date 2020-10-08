# Number of SNPs for each trait dataset
# n_snps={"matingBehavior": 2_400_000, "aggression": 1_914_528}

import itertools
from os import listdir
import os.path as path

localrules: all, fdr, clean, scratchsetup, dirsetup
configfile: "snakemakeConfig.json"

traits = config["traits"].keys()
try:
    scratchdir = config["scratchdir"]
except KeyError:
    scratchdir = None
accession_files = [f for f in listdir("Data/Identifiers/") 
                   if f.endswith(".txt")]
species = [os.path.basename(f).rstrip(".txt") for f in accession_files]
accessions = []
traits_out = []
species_out = []

for trait in traits:
    for f in accession_files:
        accession = []
        with open("Data/Identifiers/" + f, 'r') as filestream:
            for line in filestream:
                if line.strip():
                    accession.append(line.strip())
                    species_out.append(path.basename(f).rstrip(".txt"))
                    traits_out.append(trait)
        accessions.append(accession)
accessions_flat = [a for sublist in accessions for a in sublist]


rule all:
    input:
        expand("Data/Interest/{trait}.interest.txt", trait=traits) ,
        expand("Data/Interest/{trait}.windows.csv", trait=traits) ,
        expand("Data/Scratch/InterestSeqs/{trait}.interest.nsnps", trait=traits) ,
        expand("Data/Scratch/Alignments/{trait}.index.1.bt2", trait=traits) ,
        expand(
            "Data/Scratch/Alignments/Out/{trait}_{species}_{accession}.sorted.bam",
            zip, trait=traits_out, species=species_out, accession=accessions_flat
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
         "python Scripts/DataManipulation/fdrCorrection.py"
         " {params.n_snps} {input.data} {output}"


rule get_seqs:
    input:
        snps = "Data/Interest/{trait}.interest.txt"
    output:
        seqs = "Data/Scratch/InterestSeqs/{trait}.interest.fasta",
        nsnps = "Data/Scratch/InterestSeqs/{trait}.interest.nsnps"
    shell:
        "wc -l {input.snps} > {output.nsnps}"
        " && python Scripts/GetData/windows.py {input.snps} {output.seqs} "
        + config["email"]


rule unzipref:
    input: "Data/Scratch/Genomes/dMelRefSeq.fna.gz"
    output: "Data/Scratch/Genomes/dMelRefSeq.fna"
    shell: 
        "gzip -dc {input} > {output} &&"
        " bwa index {output} &&"
        " samtools faidx {output}"


rule get_dmel_genome:
    output:
        "Data/Scratch/Genomes/dMelRefSeq.fna.gz"
    script:
        "Scripts/GetData/dMelGenome.py"


rule dirsetup:
    output:
        directory("Logs/Cluster"),
        directory("Data/Scratch/tmp")
    input:
        directory("Data/Scratch")
    shell:
        "bash Scripts/Setup/dirSetup.sh {output}"


rule scratchsetup:
    output:
        directory("Data/Scratch")
    shell:
        "python Scripts/Setup/scratchSetup.py " + str(scratchdir)


rule find_overlap:
    input:
        "Data/Interest/{trait}.interest.txt"
    output:
        "Data/Interest/{trait}.windows.csv"
    shell:
        "python Scripts/DataManipulation/overlapDetector.py"
        " -i {trait}.interest.txt"
        " -w " + str(config["window_size"]) + " > Data/Interest/{trait}.windows.csv"


# Removes everything except initial dependencies.
rule clean:
    shell:
        "rm -r"
        " Logs"
        " Data/Interest"
        " Data/InterestSeqs"
        " Data/Scratch"
