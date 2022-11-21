#!/bin/bash
echo "Sorted bam file is: $1";
echo "Annotation file gtf: $2";
########## Main Program #################
out=$(basename $1 .bam);
echo "Stringtie quantification";
stringtie $1 -l $out -p 16 -G $2.gtf -o $out.gtf
echo "make txt file for gtf"
ls SRR*.gtf > gtf_file.txt
echo "merge gtf file"
stringtie --merge -p 8 -G $2.gtf -o $2.merge.gtf gtf_file.txt
echo "Final Quantification"
stringtie -e -B -p 8 -G $2.merge.gtf -o ballgown/$out/$out.gtf $1
echo "Work Done";
