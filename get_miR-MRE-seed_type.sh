#!/bin/sh
if [ $# != 3 ]
then
	echo "Usage:";
	echo "1. provide target gene name";
	echo "2. output directory";
	echo "3. look up table for seed type";
else
#	set -x
#	echo `date`
	gene=$1
	output=$2
	look_up=$3


	cat $look_up | grep -w $gene > $output/$gene.tmp

	rm $output/$gene.miR_info.txt

	echo -e "TargetGene\tmiR\tseed_type\tMRE_seq_7mer-m8\tUTR_start\tUTR_end\tcontext++_score" > $output/$gene.miR_info.txt
	cat $output/$gene.tmp | awk 'BEGIN{OFS="\t"}{print $2,$5,$6,$13,$7,$8,$9}'| sort -nk7 >> $output/$gene.miR_info.txt
	rm $output/$gene.tmp
	
#	echo `date`

fi
