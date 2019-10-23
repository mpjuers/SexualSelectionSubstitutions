#!/bin/bash

# usage: bash multialignMummer.sh <ref.fasta> <query.fasta>

nucmer -p ${2%.interest.fasta} ${1} ${2}
