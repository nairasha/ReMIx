#!/bin/sh
if [ $# != 9 ]
then
	echo "Usage:";
	echo "1. TCGA sample name";
	echo "2. Subtype, like TN, ER or HER2";
	echo "3. MRE sequence";
	echo "4. Output directory (path just before the sample dir)";
	echo "5. full path to flagstat file";
	echo "6. gene name";
	echo "7. UTR length";
	echo "8. miR name";
	echo "9. seed type";
else
#	set -x
#	echo `date`
	sample=$1
	subtype=$2
	MRE=$3
	OUTPUT=$4
	flagstat=$5
	gene=$6
	UTR=$7
	miR=$8
	seed=$9
	
	if [ -f $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE_Counts ]
	then
		rm $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE_Counts
	fi
	if [ -f $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_raw_MRE_Counts ]
	then
		rm $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_raw_MRE_Counts
	fi

	value="$gene@$miR@$seed"
	echo $value >> $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE
	cat $OUTPUT/$sample/$MRE.fimo_out/fimo.txt | grep -v "#" | cut -f6 | wc -l >> $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE
	cat $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE | tr "\n" "\t" > $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE_Counts.tmp
	sed -i -e '$a\' $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE_Counts.tmp
	cat $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE_Counts.tmp | cut -f1,2 > $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_raw_MRE_Counts
	reads=`cat $flagstat | grep $sample | cut -f2`
	cat $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_raw_MRE_Counts | awk '{print ($1"\t"((10^9*$2)/('$reads'*'$UTR')))}' > $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE_Counts

	rm $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE $OUTPUT/$sample/$MRE.$gene.$miR.$seed.FIMO_MRE_Counts.tmp
	rm -rf $OUTPUT/$sample/$MRE.fimo_out/

#	echo `date`

fi
