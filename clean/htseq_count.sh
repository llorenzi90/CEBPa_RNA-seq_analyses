#!/bin/bash

#SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=0-5

sample=$1
PPN=$2


echo -e "$PPN\t$sample"

gtf="/home/llorenzi/references/mouse/annotation/gencode/gencode.vM27.annotation.gtf"
#GENCODE annotation: gtf download link: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/gencode.vM27.basic.annotation.gtf.gz (file date: 05-May-2021 14:04, 28359457)
#downloaded may 2021
#Description:
#1. gencode.vX.annotation.{gtf,gff3}.gz:
#  Main file, gene annotation on reference chromosomes in GTF and GFF3 file formats.
#  These are the main GENCODE gene annotation files. They contain annotation (genes, 
#  transcripts, exons, start_codon, stop_codon, UTRs, CDS) on the reference chromosomes,
#  which are chr1-22, X, Y, M in human and chr1-19, X, Y, M in mouse.

cd /llorenzi/scratch/CEBPa_RNA_seq_RESEQ/$sample 

#perform htseq count 
#1) first sort the sam file by name
module load samtools/1.9 

samtools sort -@ $PPN -n -o $sample.namesorted.bam $sample.hisat2.sam

module load conda/current
conda activate htseq

#2) Run htseq-count
htseq-count -n $PPN --stranded=reverse --mode intersection-nonempty --supplementary-alignments ignore $sample.namesorted.bam $gtf

#submit: 
#PPN=8 (or another value)
#for f in /scratch/llorenzi/CEBPa_RNA_seq/*-R*;do s=$(basename $f);echo $s; sbatch -J $s.htseq_count -n $PPN htseq_count.sh $s $PPN; done

