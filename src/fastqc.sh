#!/bin/bash

#BSUB -o o.fastq
#BSUB -e e.fastq
#BSUB -n 8
#BSUB -M 40000

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

module load FastQC-0.11.2

if [ ! -d fastqc/ ]; then
	mkdir fastqc
fi

# run fastqc

fastqc -f fastq -o fastqc/ *.fastq.gz
