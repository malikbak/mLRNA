#!/bin/bash
STAR --genomeDir $1 --runThreadN 6 --readFilesIn $2 $3 --outFileNamePrefix $4 --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard
