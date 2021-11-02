#!/bin/bash

#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=8
#SBATCH --time=0-3

PPN=8

hisat2Index="/home/llorenzi/references/mouse/hisat2_index/GRCm39.primary_assembly"

sample=$1



echo $sample

module load hisat2/2.1.0

cd /scratch/llorenzi/CEBPa_RNA_seq/$sample 

hisat2 -p $PPN -x $hisat2Index -1 ${sample}_1.fq.gz -2 ${sample}_2.fq.gz -S $sample.hisat2.sam --rna-strandness RF --summary-file $sample.hisat2.summary --new-summary

#submit: for f in /scratch/llorenzi/CEBPa_RNA_seq/*-R*; do s=$(basename $f); sbatch hisat2_mm39.sh $s;done

