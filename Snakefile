# Number of SNPs for each trait dataset
# n_snps={"matingBehavior": 2_400_000, "aggression": 1_914_528}


localrules: all, fdr, clean
configfile: "snakemakeConfig.json"
traits = config["traits"].keys()


rule all:
    input:
        ["Data/Interest/" + trait + ".interest.txt" for trait in traits],
        ["Data/InterestSeqs/" + trait + ".interest.nsnps" for trait in traits],
        ["Data/Alignments/" + trait + ".interest.aligned.fasta" for trait in traits],
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


rule align_interest:
    input:
        query = "Data/InterestSeqs/{trait}.interest.fasta",
        ref = "Data/Genomes/dMelRefSeq.fna"
    output:
        align = "Data/Alignments/{trait}.interest.aligned.fasta"
    shell:
        "sh Scripts/DataManipulation/alignSeqsOfInterest.sh"
        " {input.ref} {input.query} {output.align}"
        

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
