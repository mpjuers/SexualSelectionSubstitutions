# SexualSelectionSubstitutions
An analysis of selection regimes in regions of inter- and intrasexual selection.

Dependencies:
<<<<<<< HEAD
<<<<<<< HEAD
    Python 3:
        BioPython
    MUMmer
=======
    Python 3 (tested with 3.7.0):
        snakemake (5.7.1)
=======
    Python 3 (tested with 3.7.0):
>>>>>>> 3aa6fff... Add more dependencies to README.md.
        BioPython (1.74)
    bwa (0.7.17)
    samtools (1.9)
    seqtk (1.3)
<<<<<<< HEAD
>>>>>>> fd7a617... Add snakemake, bwa, and samtools to dependencies.
=======
>>>>>>> 3aa6fff... Add more dependencies to README.md.

Snakefile usage:
    snakemake: Build the project.
    snakemake clean: Remove all generated files.
    snakemake target: Generate specified target file.

<<<<<<< HEAD
snakemakeConfig:
    json object structure:
        email: The user email for Entrez queries.
        traits: One key for each trait. Each entry should include the following keys.
            n_snps: The number of SNPs identified by the study. Used for FDR correction.
=======
snakemakeConfig.json:
    email: The user email for Entrez queries.
    traits: One key for each trait. Each entry should include the following keys.
        n_snps: The number of SNPs identified by the study. Used for FDR correction.
>>>>>>> fd7a617... Add snakemake, bwa, and samtools to dependencies.
