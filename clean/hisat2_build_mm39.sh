#!/bin/bash
#SBATCH -e %x%j.err
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --ntasks=8
#SBATCH --time=00:45:30
PPN=8
genome_file="/home/llorenzi/references/mouse/GRCm39.primary_assembly.genome.fa"
#genome download link http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M27/GRCm39.primary_assembly.genome.fa.gz    
#Download date: may 2021, file date: 05-May-2021 14:04           773873008
# 15. A.primary_assembly.genome.fa.gz (where A is the current assembly name):
#   Primary assembly genome sequence fasta file (sequence region names are the same as in the GTF/GFF3 files).
# It includes reference chromosomes and scaffolds only.


module load hisat2/2.1.0

cd /home/llorenzi/references/ncbi/mouse/hisat2_index
  
hisat2-build -p $PPN $genome_file GRCm39.primary_assembly

