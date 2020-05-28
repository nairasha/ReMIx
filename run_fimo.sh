#!/bin/sh
if [ $# != 5 ]
then
	echo "Usage:";
	echo "1. MRE sequence";
	echo "2. full path to fasta sequence file to search MRE";
	echo "3. full path to background file";
	echo "4. Output directory";
	echo "5. Input directory that contains MRE sequence in MEME format"
else
#	set -x
#	echo `date`
	MRE=$1
	FASTA=$2
	BGFILE=$3
	OUTPUT=$4
	INPUT=$5
	
#	mkdir $OUTPUT/$MRE.fimo_out
#	out=$OUTPUT/$MRE.fimo_out

	fimo=/data5/bsi/bictools/src/meme/4.11.1/BUILD/bin/fimo

	$fimo --bgfile $BGFILE --o $OUTPUT/$MRE.fimo_out $INPUT/$MRE.fimo.input $FASTA

#	echo `date`

fi
