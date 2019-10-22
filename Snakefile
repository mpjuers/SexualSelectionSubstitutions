# Number of SNPs for each trait dataset
# n_snps={"matingBehavior": 2_400_000, "aggression": 1_914_528}


configfile:
    "snakemakeConfig.json"


rule all:
    input:
        ["Data/Interest/" + trait + ".interest.txt" for trait in config["traits"].keys()],
        ["Data/InterestSeqs/" + trait + ".interest.nsnps" for trait in config["traits"].keys()],
        ["Data/Alignments/" + trait + ".interest.aligned.fasta" for trait in config["traits"].keys()],
        "Data/Genomes/dMelRefSeq.fna.gz"


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
        seqs = "Data/InterestSeqs/{trait}.interest.fasta",
        nsnps = "Data/InterestSeqs/{trait}.interest.nsnps"
    shadow:
        "full"
    shell:
         ("wc -l {input.snps} > {output.nsnps}"
          + " && python Scripts/GetData/windows.py {input.snps} {output.seqs} "
          + config["email"])


rule multialign:
    input: 
        "Data/InterestSeqs/{trait}.interest.fasta"
    output:
        "Data/Alignments/{trait}.interest.aligned.fasta"
    shell:
        ("sh Scripts/DataManipulation/alignSeqsOfInterest.sh {input} {output}"
         + " " + config["muscleparams"])


rule get_dmel_genome:
    output:
        "Data/Genomes/dMelRefSeq.fna.gz"
    script:
        "Scripts/GetData/dMelGenome.py"
        
        
# Removes everything except initial dependencies.
rule clean:
    shell:
        "rm -r"
        " Logs"
        " Data/Interest"
        " Data/InterestSeqs"
        " Data/Alignments"
