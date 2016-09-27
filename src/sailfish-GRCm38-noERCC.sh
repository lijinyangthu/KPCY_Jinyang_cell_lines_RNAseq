#!/bin/bash

set -e
set -o pipefail

# using Sailfish version 0.6.3 to pseudo-align and measure transcript abundances against mm10

#BSUB -o o.sailquant
#BSUB -e e.sailquant
#BSUB -n 12
#BSUB -M 64000
#BSUB -R "span[hosts=1]"
#BSUB -N
#BSUB -u "balli.dave@gmail.com"

# usage from main Confetti_RNAseq/ folder :
# bsub < scripts/sailfish-mm10_array.sh

# fastq files must in format: PD287_Y1_R1.fastq
#  first run the following in data/fastq/ folder : 
# $ find . -name "*_R1.fq.gz" | xargs -n1 -I "{}" basename {} _R1.fq.gz > sample_files.txt

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

module load sailfish-0.6.3

# load files IDs into array 
declare -a samples=( `find data/fastq/ -name "*_R2.trim.fastq.gz" | xargs -n1 -I "{}" basename {} _R2.trim.fastq.gz`)

for i in "${samples[@]}"
do
	sailfish quant \
	-i /home/dballi/GRCm38_genome \
	-l "T=PE:O=><:S=U" \
	-1 <(gunzip -c data/fastq/"$i"_R1.trim.fastq.gz) \
	-2 <(gunzip -c data/fastq/"$i"_R2.trim.fastq.gz) \
	-o data/output_"$i" \
	-p 32
done

for i in "${samples[@]}"
do
	if [ -e data/output_"$i"/quant.sf ]
		then
			cp data/output_"$i"/quant.sf data/output_"$i"/"$i"-sailfish.txt
			cp -r data/output_"$i"/ results/
	fi
done

exit 0
