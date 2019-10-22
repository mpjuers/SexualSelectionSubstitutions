# SexualSelectionSubstitutions
An analysis of selection regimes in regions of inter- and intrasexual selection.

Snakefile usage:
    snakemake: Build the project.
    snakemake clean: Remove all generated files.

snakemakeConfig:
    json object structure:
        traits: One key for each trait. Each entry should include the following keys.
            n_snps: The number of SNPs identified by the study. Used for FDR correction.
