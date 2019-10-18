#!/usr/bin/env python3

import sys

import pandas as pd

args = sys.argv

fdr = 0.20
n_snps = int(args[1])
data = pd.read_csv(args[2])
data.columns = (data.columns
                .str.lower()
                .str.replace(' ', '_')
                .str.replace('-', '_'))

data["bh_crit"] = pd.Series(range(1, len(data.index) + 1)) / n_snps * fdr
try:
    data["p_value"] = data.p_value.str.rstrip(r'*').astype(float)
except AttributeError:
    pass
snps_of_interest = data[data.p_value < data.bh_crit].genomic_location
snps_of_interest.to_csv(args[3], sep="\n", index=False)
