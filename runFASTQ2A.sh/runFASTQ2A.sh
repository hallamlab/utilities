#!/bin/bash

#This is a bash script using a sed command to do the actual conversion of FASTQ to FASTA. It requires a list of FASTQ files to convert (even if the list is only a single .fastq file). The output are .fa (fasta) files with the same prefix. For example, IX0937_D19JLACXX_1_CCAACA.fastq would become IX0937_D19JLACXX_1_CCAACA.fa. 

#USAGE:first create the list of fastq files to be converted into fasta files:
#ls >fastqs_to_convert.list
#then launch the command like so:
#./runFASTQ2A.sh fastqs_to_convert.list

list=$1

while read line
do
	tmp=$(echo ${line%%.fastq}) #### NOTE: if the file ends in the common ".fq" change the .fastq to .fq
	fastaFile=$tmp.fa

	sed -n '1~4s/^@/>/p;2~4p' $line > $fastaFile
done<$list
