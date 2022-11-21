#!/usr/bin/env python
# coding: utf-8

# In[3]:


import argparse
import subprocess
import os
# Initialize parser
parser = argparse.ArgumentParser()

# Adding optional argument

parser.add_argument("-as", "--Accession", dest="accession", help = "Download sra file from NCBI SRA database")

parser.add_argument("-ds", "--fastq-dump", dest="dump", help = "Convert sra file into fastq")

parser.add_argument("-f", "--Forward", dest="Forward", help = "Forward read")

parser.add_argument("-r", "--Reverse", dest="Reverse", help = "Reverse read")

parser.add_argument("-g", "--Reference", dest="Reference", help = "Prefix of Reference index")

parser.add_argument("-t", "--GTF", dest="GTF", help = "GTF file for genome annotation")

parser.add_argument("-G", "--Genome", dest="genome", help = "Genome fasta file")

parser.add_argument("-p", "--Prefix", dest="Prefix", help = "Name for index, write only prefix")

parser.add_argument("-b", "--Bam", dest="Bam", help = "bam file for stringtie assembly")

parser.add_argument("-is", "--IS", dest="StarI", help = "STAR index directory")

parser.add_argument("-os", "--OS", dest="OutS", help = "Out put file name for STAR alignment: just prefix like sample name")

parser.add_argument("-md", "--metadata", dest="metadata", help = "Metadata csv file for Expression analysis")

parser.add_argument("-pb", "--path", dest="path", help = "Path of ballgown folder like if have folder name ballgown just write ballgown or /path/ballgown")

parser.add_argument("-c", "--condition", dest="condition", help = "Condition column name from csv file")
# Read arguments from command line
args = parser.parse_args()

if args.accession:
    subprocess.run(["./scripts/download.sh", args.accession])
elif args.dump:
    subprocess.run(["./scripts/dump.sh", args.dump])
elif args.Forward and args.Reverse and args.Reference:
    subprocess.run(["./scripts/bowtie_stringtie.sh", args.Forward, args.Reverse, args.Reference])
elif args.genome and args.Prefix:
    subprocess.run(["./scripts/indexing.sh", args.genome, args.Prefix])
elif args.genome and args.GTF:
    subprocess.run(["./scripts/STAR.sh", args.genome, args.GTF])
elif args.Bam and args.GTF:
    subprocess.run(["./scripts/stringtie.sh", args.Bam, args.GTF])
elif args.StarI and args.Forward and args.Reverse and args.OutS:
    subprocess.run(["./scripts/star_alignment.sh", args.StarI, args.Forward, args.Reverse, args.OutS])
elif args.metadata and args.path and args.condition:
    subprocess.run(["Rscript","Expression.R","-m", args.metadata, "-p", args.path, "-c", args.condition])


# In[ ]:




