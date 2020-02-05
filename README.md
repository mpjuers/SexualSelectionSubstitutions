# SexualSelectionSubstitutions
An analysis of selection regimes in regions of inter- and intrasexual selection.

## Rules

*Please use
[mathematica-notebook-filter](https://github.com/JP-Ellis/mathematica-notebook-filter)
if committing Mathematica notebooks for code cleanliness.*
Additionally, do not use graphics as input to functions as this makes the commit
binary and thus we can't do diffs. Use the `%` notation instead.

*Jupyter notebooks should be converted to `.py` (percent format) or `.md`
using [Jupytext](https://github.com/mwouts/jupytext).*

Directories should be PascalCase and files should be camelCase.

Python scripts should be formatted with [Black](https://github.com/psf/black).
Python scripts should be written to work with input and output arguments.
`sys.argv` is fine.

## Dependencies:

- Python 3 (tested with 3.7.0):
    - snakemake (5.7.1)
    - BioPython (1.74)
- bwa (0.7.17)
- samtools (1.9)
- seqtk (1.3)
- pandas (0.25.1)
- sra-toolkit (2.9.6)

## Snakefile usage:

- `snakemake`: Build the project.
- `snakemake clean`: Remove all generated files.
- `snakemake <target>`: Generate specified target file.


### From workdir on cluster:

`qsub -d ${PWD} Scripts/Batch/execSnakemake.batch` OR
`qsub -cwd Scripts/Batch/execSnakemake.batch`

### snakemakeConfig.json:

- email: The user email for Entrez queries.
- scratchdir: One can specify an alternate directory for many intermediate
    files on, e.g., a high-speed filesystem. This directory will be symlinked
    to `Data/Scratch` in the project's root directory. If "None" (or omitted),
    the scratch directory will be in the same location, but will not be a
    symlink.
- traits: One key for each trait. Each entry should include the following keys.
    - n_snps: The number of SNPs identified by the study. Used for FDR correction.
- 

### clusterConfig.json:

Configuration parameters for rules on a per-rule basis.

## Data Sources

- *Drosophila melanogaster*
    - BioProject url: `https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=4224886`

- *Drosophila simulans*
    - BioProject url: `https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=3616273`
