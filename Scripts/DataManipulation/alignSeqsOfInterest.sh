#!/bin/bash

# usage: sh alignSeqsOfInterest.sh <input.fasta> <output.fasta> <args>
# In workflow, <args> is read from snakemakeConfig["muscleparams"].

if [ ! -d Logs ]; then
  mkdir -p Logs;
fi
logfile=${1##*/}
logfile=${logfile%.*}

muscle -log Logs/"$logfile".multialign.log -in ${1} -out ${2} ${3} 
