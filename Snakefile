# Number of SNPs for each trait dataset
# n_snps={"matingBehavior": 2_400_000, "aggression": 1_914_528}

import itertools
from os import listdir
import os.path as path

localrules: all, fdr, clean, dirsetup
configfile: "snakemakeConfig.json"

traits = config["traits"].keys()
datasets = config["datasets"].keys()
print(datasets)

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

rule all:
    input:
        expand("Data/SNPData/{trait}.csv", trait=traits),
        expand("Data/Interest/{trait}.interest.txt", trait=traits) ,
        expand("Data/Interest/{trait}.windows.csv", trait=traits) ,
        expand(scratchdir + "/InterestSeqs/{trait}.interest.nsnps", trait=traits) ,
        scratchdir + "Genomes/dMelRefSeq.fna.gz", 
        scratchdir + "Genomes/dMelRefSeq.fna"


rule get_dmel_genome:
    output:
        scratchdir + "Genomes/dMelRefSeq.fna.gz"
    script:
        "Scripts/GetData/dMelGenome.py"


rule unzipref:
    input: scratchdir + "Genomes/dMelRefSeq.fna.gz"
    output: scratchdir + "Genomes/dMelRefSeq.fna"
    shell: 
        "gzip -dc {input} > {output} &&"
        " bwa index {output} &&"
        " samtools faidx {output}"


rule dirsetup:
    output:
        touch(".mkdir_checkpoint")
    params:
        directory("Logs/Cluster"), 
        directory("Data/Genomes"),
        scratchdir + directory("tmp"),
        directory("Data/InterestSeqs"),
        directory("Data/Interest")
    shell:
        "python Scripts/Setup/scratchSetup.py " + scratchdir + " &&"
        " mkdir -p {params} || true"


rule fdr:
    input:
        checkpoint = ".mkdir_checkpoint",
        data = "Data/SNPData/{trait}.csv"
    output:
        "Data/Interest/{trait}.interest.txt"
    params:
        n_snps = lambda wildcards: config["traits"][wildcards.trait]["n_snps"]
    shell:
         "python Scripts/DataManipulation/fdrCorrection.py"
         " {params.n_snps} {input.data} {output}"


rule find_overlap:
    input:
        "Data/Interest/{trait}.interest.txt"
    output:
        "Data/Interest/{trait}.windows.csv"
    shell:
        "python Scripts/DataManipulation/overlapDetector.py"
        " -i {input} -w " + str(config["window_size"]) + " > {output}"


rule get_seqs:
    input:
        windows = "Data/Interest/{trait}.windows.csv",
        snps = "Data/Interest/{trait}.interest.txt"
    output:
        seqs = scratchdir + "/InterestSeqs/{trait}.interest.fasta",
        nsnps = scratchdir + "/InterestSeqs/{trait}.interest.nsnps"
    shell:
        "wc -l {input.snps} > {output.nsnps}"
        " && python Scripts/GetData/windows.py {input.windows} {output.seqs} "
        + config["email"] + str(config["window_size"])


# Removes everything except initial dependencies.
rule clean:
    shell:
        "rm -r"
        " Logs"
        " Data/Interest"
        " Data/InterestSeqs"
        " Data/Scratch"
        " .mkdir_checkpoint"
        " || true"
