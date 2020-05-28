#!/bin/sh
if [ $# != 1 ]
then
	echo "Usage:";
	echo "Please provide configuration file (include complete path)";
else
	set -x
	echo `date`
	config=$1

	BAM=$( cat $config | grep -w '^BAM' | cut -d '=' -f2)
	UTR_BED_dir=$( cat $config | grep -w '^3UTR_BED_DIR' | cut -d '=' -f2)
	output_dir=$( cat $config | grep -w '^OUTPUT_DIR' | cut -d '=' -f2)
	sample_name=$( cat $config | grep -w '^SAMPLE_NAME' | cut -d '=' -f2)
	samtools=$( cat $config | grep -w '^SAMTOOLS' | cut -d '=' -f2)
	tophat=$( cat $config | grep -w '^TOPHAT' | cut -d '=' -f2)
	genes=$( cat $config | grep -w '^GENES' | cut -d '=' -f2)
	line=$SGE_TASK_ID

	geneNames=$( echo $genes | tr ":" "\n" )
	i=1
	for gene in $geneNames
	do
		geneArray[$i]=$gene
		let i=i+1
	done

	gene_name=${geneArray[$line]}

	mkdir $output_dir/$sample_name
	output_dir=$output_dir/$sample_name

#### Intersect with the 3'UTR BED file of the gene - BED file to contain chr, start, end coordinates of UTR
	$samtools/samtools view -hf 2 $BAM -L $UTR_BED_dir/$gene_name.bed | $samtools/samtools view -bS - > $output_dir/$gene_name.3UTR.bam
	$samtools/samtools index $output_dir/$gene_name.3UTR.bam $output_dir/$gene_name.3UTR.bam.bai

#### Convert BAM to fastq
	$tophat/bam2fastx -q -Q -A -N -o $output_dir/$gene_name.fastq $output_dir/$gene_name.3UTR.bam

#### Convert fastq to fasta
	echo "> " $sample_name > $output_dir/$gene_name.fasta
	cat $output_dir/$gene_name.fastq |perl -e '$i=0;while(<>){if(/^\@/&&$i==0){s/^\@/\>/;print;}elsif($i==1){print;$i=-3}$i++;}' | perl -ne 'print unless (0 != $. % 2)' >> $output_dir/$gene_name.fasta

### Remove temporary files
	rm $output_dir/$sample_name.3UTR.bam $output_dir/$sample_name.fastq $output_dir/$gene_name.fastq $output_dir/$gene_name.3UTR.bam

	echo `date`

fi
