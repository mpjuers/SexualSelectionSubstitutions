# SexualSelectionSubstitutions
An analysis of selection regimes in regions of inter- and intrasexual selection.

Dependencies:
    Python 3 (tested with 3.7.0):
        snakemake (5.7.1)
        BioPython (1.74)
    bwa (0.7.17)
    samtools (1.9)
    seqtk (1.3)
    pandas (0.25.1)

Snakefile usage:
    snakemake: Build the project.
    snakemake clean: Remove all generated files.
    snakemake target: Generate specified target file.


From workdir on cluster:
    qsub ${PWD} Scripts/Batch/execSnakemake.batch

snakemakeConfig.json:
    email: The user email for Entrez queries.
    traits: One key for each trait. Each entry should include the following keys.
        n_snps: The number of SNPs identified by the study. Used for FDR correction.

clusterConfig.json:
    Configuration parameters for rules on a per-rule basis.
