rule fdr:
    input:
        "Data/SNPData/{traits}.csv"
    output:
        "Data/Interest/{traits}.interest.txt"
    shell:
        "Rscript Scripts/DataManipulation/fdrCorrection.r {input}"
