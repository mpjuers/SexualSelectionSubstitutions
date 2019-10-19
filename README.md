# SexualSelectionSubstitutions
An analysis of selection regimes in regions of inter- and intrasexual selection.

Snakefile usage:
    snakemake: Build the project.
    snakemake clean: Remove certain generated files.
    snakemake clean_all: Remove all generated files. This will restore the Data
        directory to the state of the Repository.

snakemakeConfig:
    json object structure:
        traits: One key for each trait. Each entry should include the following keys.
            n_snps: The number of SNPs identified by the study. Used for FDR correction.
