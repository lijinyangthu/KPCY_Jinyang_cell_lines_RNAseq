#!/bin/bash
# bsub -n 4 -M 10000 -e e.seqtk -o o.seqtk sh src/seqtk_trim.sh

#BSUB -n 6
#BSUB -M 10000
#BSUB -e e.seqtk
#BSUB -o o.seqtk

declare -a fastq=( ` find data/fastq/ -name "*.fastq.gz" | xargs -n1 -I "{}" basename {} .fastq.gz ` )

for i in "${fastq[@]}";
do
	seqtk trimfq data/fastq/"$i".fastq.gz > data/fastq/"$i".trim.fastq
	pigz data/fastq/"$i".trim.fastq
done
