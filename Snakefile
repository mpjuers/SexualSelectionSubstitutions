# Number of SNPs for each trait dataset
# n_snps={"matingBehavior": 2_400_000, "aggression": 1_914_528}


localrules: all, fdr, clean
configfile: "snakemakeConfig.json"
traits = config["traits"].keys()


rule all:
    input:
        ["Data/Interest/" + trait + ".interest.txt" for trait in traits],
        ["Data/Scratch/InterestSeqs/" + trait + ".interest.nsnps" for trait in traits],
        ["Data/Scratch/Alignments/" + trait + ".interest.aligned.fasta" for trait in traits],
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
        "Logs/Cluster"
    shell:
        "bash Scripts/Setup/dirSetup.sh {output}"


rule scratchsetup:
    output:
        "Data/Scratch"
    shell:
        "if [[ ! -h {output} ]]; then ln -s {config.scratchdir} {output}; fi"


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
