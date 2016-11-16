#!/bin/bash

# STAR alignment for RNAseq from mouse KPCY experiments

# load module STAR v2.5.2a from PMACS modules
if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
   fi

set -e
set -o pipefail

#BSUB -e e.star
#BSUB -o o.star
#BSUB -n 10
#BSUB -M 64000
#BSUB -R "span[hosts=1]"
#BSUB -N
#BSUB -u "balli.dave@gmail.com"

if [ -f /etc/profile.d/modules.sh ]; then
   source /etc/profile.d/modules.sh
fi

module load STAR-2.5.2a

# alignment to mm10 genome
REFERENCE="/home/dballi/KPCY_Jinyang_cell_lines/src/reference"
ANNO="/home/dballi/KPCY_Jinyang_cell_lines/src/reference/genes_mm10_ercc.gtf"
DATE=`date +'%Y-%m-%d'`

# load basename of fastq files into array and run  a 'for loop' using STAR for alignment of each read pair
declare -a samples=(`find data/fastq/ -maxdepth 1 -name "*_R1.trim.fastq.gz" | xargs -n1 -I "{}" basename {} _R1.trim.fastq.gz`)

for i in "${samples[@]}"
do
   STAR \
   --runMode alignReads \
   --genomeDir "$REFERENCE" \
   --sjdbGTFfile "$ANNO" \
   --readFilesIn data/fastq/"$i""_R1.trim.fastq.gz" data/fastq/"$i""_R2.trim.fastq.gz" \
   --runThreadN 16 \
   --readFilesCommand gunzip -c \
   --outSAMtype BAM SortedByCoordinate \
   --outFileNamePrefix data/bam/"$i"-
done

# using featureCounts (from Subread package) to bin fastq reads to gene features using mm10 gtf file
featureCounts \
-T 16 \
-p \
-D 2000 \
-a "$ANNO" \
-t exon \
-g gene_id \
-o results/"$DATE"-star-counts.txt \
`find data/bam/ -iname "*.bam" -mtime -1`

exit 0
