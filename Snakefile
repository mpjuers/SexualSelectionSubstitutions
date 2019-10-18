# Number of SNPs for each trait dataset
n_snps={"matingBehavior": 2_400_000}


rule all:
    input:
        ["Data/Interest/" + x + ".interest.txt" for x in n_snps.keys()]


rule fdr:
    input:
        data = "Data/SNPData/{trait}.csv"
    output:
        "Data/Interest/{trait}.interest.txt"
    params:
        n_snps = lambda wildcards: n_snps[wildcards.trait]
    shell:
        "python Scripts/DataManipulation/fdrCorrection.py"
        " {params.n_snps} {input.data} {output}"


rule clean:
    shell:
        "rm -r Data/Interest"
