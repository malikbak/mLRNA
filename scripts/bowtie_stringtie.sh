#!/bin/bash
echo "Forward reads: $1";
echo "Reverse reads: $2";
echo "Refernce prefix: $3";
########## Main Program #################
out=$(basename $1 _1.fastq);
echo "Alignemt using bowtie2 for $out";
bowtie2 --threads 8 -x $3 -1 $1 -2 $2 -S $out.sam
echo "Sorting the Alignment file";
samtools sort -@ 16 -o $out.bam $out.sam
