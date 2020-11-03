#!/usr/bin/env python3

# Get sequences from NCBI.
# To be called from Snakefile.
# Usage: python windows.py <infile> <outfile> <email> <window_size>

import os
import sys

from Bio import Entrez
from Bio import SeqIO
import pandas as pd


def main():
    snpfile = sys.argv[1]
    outfile = sys.argv[2]
    email = sys.argv[3]
    window_size = sys.argv[4]

    interest = pd.read_csv(
        snpfile,
        header=0,
        dtype={"chrom": str, "start": np.int64, "end": np.int64},
    )
    interest.columns = interest.columns.str.lower().str.replace(" ", "_")
    interest[["chrom", "start", "end"]] = (
        interest.iloc[:, 0]
        .str.replace("_[A-Z]{3}?", "")
        .str.replace(" ", "")
        .str.split("_", expand=True)
    )
    interest.index.rename("Index", inplace=True)
    summary = pd.read_csv("Data/dmelSummary.csv")
    summary.index.rename("Index", inplace=True)
    summary.columns = summary.columns.str.lower().str.replace(" ", "_")

    seqs = []
    Entrez.email = email
    for index, row in interest.iterrows():
        with Entrez.efetch(
            db="nucleotide",
            id=summary.loc[summary["name"] == row["chrom"], "refseq"].iat[0],
            rettype="fasta",
            strand=1,
            seq_start=row["start"],
            seq_stop=row["end"],
        ) as handle:
            seqs.append(SeqIO.read(handle, "fasta"))
    SeqIO.write(seqs, outfile, "fasta")


if __name__ == "__main__":
    main()
