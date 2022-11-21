
# Instructions for Pipeline

There are some instructions for running the RNAseq pipeline. Pipeline is divided into four steps.

1- Downloading the dataset

2- Alignment and sorting the Alignment file

3- Quantification of RNAseq dataset

4- Expression analysis




## Installation

Install require sudo packages

```bash
  sudo apt install bowtie2
  sudo apt install samtools
  sudo apt install stringtie
```
Install require R packages

```bash
install.packages("BiocManager")
BiocManager::install("ballgown")
BiocManager::install("genefilter")
install.packages("dplyr")
install.packages("devtools")
devtools::install_github('alyssafrazee/RSkittleBrewer')
BiocManager::install("EnhancedVolcano")
```
## Usage 

Help command

```bash
python pipeline.oy -h
usage: pipeline.py [-h] [-as ACCESSION] [-ds DUMP] [-f FORWARD] [-r REVERSE] [-g REFERENCE] [-t GTF] [-G GENOME]
                   [-p PREFIX] [-b BAM] [-is STARI] [-os OUTS] [-md METADATA] [-pb PATH] [-c CONDITION]

optional arguments:
  -h, --help            show this help message and exit
  -as ACCESSION, --Accession ACCESSION
                        Download sra file from NCBI SRA database
  -ds DUMP, --fastq-dump DUMP
                        Convert sra file into fastq
  -f FORWARD, --Forward FORWARD
                        Forward read
  -r REVERSE, --Reverse REVERSE
                        Reverse read
  -g REFERENCE, --Reference REFERENCE
                        Prefix of Reference index
  -t GTF, --GTF GTF     GTF file for genome annotation
  -G GENOME, --Genome GENOME
                        Genome fasta file
  -p PREFIX, --Prefix PREFIX
                        Name for index, write only prefix
  -b BAM, --Bam BAM     bam file for stringtie assembly
  -is STARI, --IS STARI
                        STAR index directory
  -os OUTS, --OS OUTS   Out put file name for STAR alignment: just prefix like sample name
  -md METADATA, --metadata METADATA
                        Metadata csv file for Expression analysis
  -pb PATH, --path PATH
                        Path of ballgown folder like if have folder name ballgown just write ballgown or
                        /path/ballgown
  -c CONDITION, --condition CONDITION
                        Condition column name from csv file
```
Download the sequence from SRA database of NCBI

```bash
python pipeline.py -a SRR1182374 #write assession number
```
Convert sra file into fastq

```bash
python pipeline.py -d SRR1182374/SRR1182374.sra #sra file to dump
```
Make index of fasta Genome

```bash
python pipeline.py -G Homo_sapiens.GRCh37.dna.toplevel.fa -p human
```

Mapping of fastq file

```bash
python pipeline.py -f SRR1182374_1.fastq -r SRR1182374_2.fastq -g human
-f for Forward read 
-r for Reverse read 
-g for Prefix for Reference index
```
Quantification:

Please activate stringtie conda env before running this command.

```bash
conda activate stringtie # in current server
```

```bash
python pipeline.py -b SRR1182375.bam -t Homo_sapiens.GRCh38.107

-b for sample bam file
-t for prefix of gtf file
```
## Expression analysis

There are some precuation for Expression analysis:

1- Metadata file should be name as metadata.txt

2- All sample proceed in the same directory where main.py file present.

```bash
python main.py -md metadata.csv -pb /path/ballgown -c condition
```


if you are using it in your system where i test it then it will run fine. if you are gonna test it in another server you may get some error in running Rscript.