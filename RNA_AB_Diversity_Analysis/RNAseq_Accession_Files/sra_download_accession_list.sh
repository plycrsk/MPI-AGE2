#!/bin/bash

FILENAME="accession_RNA.txt"
LINES=$(cat $FILENAME)

for LINE in $LINES
do
	
echo ${LINE} "processing"
fastq-dump --gzip --skip-technical --split-3 --readids --clip ${LINE}

done
exit

