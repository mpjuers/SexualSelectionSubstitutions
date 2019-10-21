# Number of SNPs for each trait dataset
# n_snps={"matingBehavior": 2_400_000, "aggression": 1_914_528}


configfile:
    "snakemakeConfig.json"


rule all:
    input:
        ["Data/Interest/" + x + ".interest.txt" for x in config["traits"].keys()],


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


rule clean:
    shell:
        "rm -r Data/Interest"


rule clean_all:
    shell:
        "rm -r Data/Interest"
