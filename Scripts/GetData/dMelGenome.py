#!/usr/bin/env python3

import os
import urllib.request as ulr

url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz"
outdir = "Data/Genomes/"
outfile = "dMelRefSeq.fna.gz"
if not os.path.isdir(outdir):
    os.makedirs(outdir)
ulr.urlretrieve(url, filename=outdir + outfile)
