#!/bin/bash
#SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=10
#SBATCH --time=02:05:30

PPN=10
echo $1

cd /scratch/llorenzi/CEBPa_RNA_seq/$1

#1) samtools sort and index
module load samtools/1.9
samtools sort -@ $PPN -O bam -o $1.hisat2.sorted.bam $1.hisat2.sam
samtools index -@ $PPN $1.hisat2.sorted.bam 


#2) bamCoverage to generate bigWig files
module load conda/current
conda activate deeptoolsenv

bamCoverage --numberOfProcessors $PPN --binSize 10 --normalizeUsing CPM --minMappingQuality 30 --filterRNAstrand forward --bam $1.hisat2.sorted.bam -o $1.hisat2.sorted.bam.CPM.fwd.bw


bamCoverage --numberOfProcessors $PPN --binSize 10 --normalizeUsing CPM --minMappingQuality 30 --filterRNAstrand reverse --bam $1.hisat2.sorted.bam -o $1.hisat2.sorted.bam.CPM.rev.bw

