#!/usr/bin/env python3

# Get sequences from NCBI.
# To be called from Snakefile.
# Usage: python windows.py <infile> <outfile> <email>

import os
import sys

from Bio import Entrez
from Bio import SeqIO
import pandas as pd


def main():
    window_size = 10_000
    snpfile = sys.argv[1]
    outfile = sys.argv[2]
    email = sys.argv[3]

    interest = pd.read_csv(snpfile, header=0)
    interest.columns = interest.columns.str.lower().str.replace(" ", "_")
    interest[["chrom", "location"]] = (interest.iloc[:, 0]
        .str.replace("_[A-Z]{3}?", "")
	.str.replace(" ", "")
	.str.split("_", expand=True))
    interest["location"] = interest["location"].astype(int)
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
            seq_start=row["location"] - window_size / 2,
            seq_stop=row["location"] + window_size / 2,
        ) as handle:
            seqs.append(SeqIO.read(handle, "fasta"))
    if not os.path.exists("Data/InterestSeqs"):
        os.makedirs("Data/InterestSeqs")
    SeqIO.write(seqs, outfile, "fasta")


if __name__ == "__main__":
    main()
