#!/usr/bin/env python3

import gzip

from Bio import SeqIO
from Bio.Alphabet import IUPAC

with gzip.open("../Data/Genomes/dMelRefSeq.gb", "rt") as handle:
    records = list(SeqIO.parse(handle, "genbank", IUPAC.ambiguous_dna))

print(records)
