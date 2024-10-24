#!/bin/bash
#SBATCH -o %x%j.out
#SBATCH --mail-user=lucialorenzi90@gmail.com
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=03:05:30

module load fastqc/0.11.8
pat=$1 #example p30-LPS

echo /scratch/llorenzi/CEBPa_RNA_seq/*$pat*
for fi in /scratch/llorenzi/CEBPa_RNA_seq/*$pat*
do 
	echo $fi
	samp=$(basename $fi)
	cd $fi
	mkdir fastqc
	fastqc ${samp}_1.fq.gz -o fastqc
	fastqc ${samp}_2.fq.gz -o fastqc

done



#list_of_samples.txt:
# EV-LPS-R1
# EV-LPS-R2
# EV-LPS-R3
# EV-R1
# EV-R2
# EV-R3
# p30-LPS-R1
# p30-LPS-R2
# p30-LPS-R3
# p30-R1
# p30-R2
# p30-R3
# p42-LPS-R1
# p42-LPS-R2
# p42-LPS-R3
# p42-R1
# p42-R2
# p42-R3

#run: 
#for pat in EV p30 p42; do sbatch -J $pat.fastqc fastqc.sh $path;done
