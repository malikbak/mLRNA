#!/bin/bash
mkdir -p ./STAR1
STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./STAR1 --genomeFastaFiles $1 --sjdbGTFfile $2


