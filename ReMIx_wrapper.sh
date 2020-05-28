#!/bin/sh
if [ $# != 1 ]
then
	echo "Usage:";
	echo "Please provide configuration file (include complete path)";

else
#	set -x
	echo `date`
	config=$1

	genes=$( cat $config | grep -w '^GENES' | cut -d '=' -f2)
	utrs=$( cat $config | grep -w '^UTR_LENGTHS' | cut -d '=' -f2)
	look_up=$( cat $config | grep -w '^LOOK_UP_TABLE' | cut -d '=' -f2)
	subtype=$( cat $config | grep -w '^SUBTYPE' | cut -d '=' -f2)
	input_dir=$( cat $config | grep -w '^INPUT_DIR' | cut -d '=' -f2)
	per_sample_per_gene_BAM_fasta=$( cat $config | grep -w '^PER_GENE_MAPPED_FASTA' | cut -d '=' -f2)
	meme=$( cat $config | grep -w '^MEME' | cut -d '=' -f2)
	output=$( cat $config | grep -w '^OUTPUT_DIR' | cut -d '=' -f2)
	script_path=$( cat $config | grep -w '^SCRIPT_PATH' | cut -d '=' -f2)
	line=$SGE_TASK_ID

	geneNames=$( echo $genes | tr ":" "\n" )
        i=1
        for gene in $geneNames
        do
                geneArray[$i]=$gene
                let i=i+1
        done

        gene_name=${geneArray[$line]}

	utrLengths=$( echo $utrs | tr ":" "\n" )
        i=1
        for utr in $utrLengths
        do
                utrArray[$i]=$utr
                let i=i+1
        done

        utr_length=${utrArray[$line]}

	sh $script_path/get_miR-MRE-seed_type.sh $gene_name $output $look_up
	mkdir $output/miR_contextScore
	perl $script_path/select_miR_using_contextScore.pl $output/$gene_name.miR_info.txt $output/miR_contextScore/$gene_name.miR_contextScore.txt
	sed -i '1d' $output/miR_contextScore/$gene_name.miR_contextScore.txt	
	num_miRs=`cat $output/miR_contextScore/$gene_name.miR_contextScore.txt | wc -l`

	miR=`cat $output/miR_contextScore/$gene_name.miR_contextScore.txt | cut -f2 | tr "\n" ":"`
	seed=`cat $output/miR_contextScore/$gene_name.miR_contextScore.txt | cut -f3 | tr "\n" ":"`
	MRE=`cat $output/miR_contextScore/$gene_name.miR_contextScore.txt | cut -f4 | tr "\n" ":"`
	
	miR_names=$( echo $miR | tr ":" "\n" )
	i=1
	for tmp_miR in $miR_names
	do
		miRArray[$i]=$tmp_miR
		let i=i+1
	done

	seed_names=$( echo $seed | tr ":" "\n" )	
	i=1
	for tmp_seed in $seed_names
	do
		seedArray[$i]=$tmp_seed
		let i=i+1
	done

	MRE_names=$( echo $MRE | tr ":" "\n" )
	i=1
	for tmp_MRE in $MRE_names
	do
		MREArray[$i]=$tmp_MRE
		let i=i+1
	done
	
	for ((count=1;count<=$num_miRs;count++));
	do
		for j in `cat $input_dir/${subtype}_paired_Samples | cut -f2 | sed 's/#/./g'`
		do 
			$script_path/run_fimo.sh ${MREArray[$count]} $per_sample_per_gene_BAM_fasta/$j/$gene_name.fasta $script_path/bgfile.hg19_allchr $input_dir/$j $meme/${seedArray[$count]}/
			$script_path/get_counts_individual_MRE.sh $j $subtype ${MREArray[$count]} $input_dir $input_dir/${subtype}_paired_Samples_flagstat $gene_name $utr_length ${miRArray[$count]} ${seedArray[$count]} 
			if [ -f $input_dir/FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]} ]
			then
				rm $input_dir/FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}
				rm $input_dir/FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw
			fi
		done
		for j in `cat $input_dir/${subtype}_paired_Samples | cut -f2 | sed 's/#/./g'`
		do
			cat $input_dir/$j/${MREArray[$count]}.$gene_name.${miRArray[$count]}.${seedArray[$count]}.FIMO_MRE_Counts >> $input_dir/FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}
			cat $input_dir/$j/${MREArray[$count]}.$gene_name.${miRArray[$count]}.${seedArray[$count]}.FIMO_raw_MRE_Counts >> $input_dir/FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw
		done

		cp $input_dir/sample_names $input_dir/$gene_name.sample_names
	        paste $input_dir/$gene_name.sample_names $input_dir/FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]} | tr "#" "\t" | sort -k2 > $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}
		paste $input_dir/$gene_name.sample_names $input_dir/FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw | tr "#" "\t" | sort -k2 > $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw
		perl $script_path/summarize_FIMO_median_norm_counts.pl $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]} $input_dir/scores.summarized.$gene_name.${miRArray[$count]}.${seedArray[$count]}
		
	done

	if [ -f $output/$gene_name.$subtype.miR_info.detailed.txt ]
	then
		rm $output/$gene_name.$subtype.miR_info.detailed.txt
		rm $output/$gene_name.$subtype.miR_info.raw_detailed.txt
	fi
	
	if [ -f $output/$gene_name.$subtype.miR_info.summary.txt ]
	then
		rm $output/$gene_name.$subtype.miR_info.summary.txt
	fi

	cat $output/$gene_name.miR_info.txt >> $output/$gene_name.$subtype.miR_info.detailed.txt
	cat $output/$gene_name.miR_info.txt >> $output/$gene_name.$subtype.miR_info.raw_detailed.txt
	echo -e "target_gene\tconserved_miR\tNormal_normalized_median\tTumor_normalized_median\tt-statistic\tp-value\tnull_hypothesis" >> $output/$gene_name.$subtype.miR_info.summary.txt

	for ((count=1;count<=$num_miRs;count++));
	do
		echo -e "\n${miRArray[$count]}\nSample_ID\tSubtype_Tissue\tMRE_sequence\tMRE_frequency" >> $output/$gene_name.$subtype.miR_info.detailed.txt
		echo -e "\n${miRArray[$count]}\nSample_ID\tSubtype_Tissue\tMRE_sequence\tMRE_frequency" >> $output/$gene_name.$subtype.miR_info.raw_detailed.txt
		cat $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]} >> $output/$gene_name.$subtype.miR_info.detailed.txt
		cat $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw >> $output/$gene_name.$subtype.miR_info.raw_detailed.txt
		cat $input_dir/scores.summarized.$gene_name.${miRArray[$count]}.${seedArray[$count]} >> $output/$gene_name.$subtype.miR_info.summary.txt
	done

	echo -e "Sample_ID\tType" > $output/$gene_name.header
	awk NF $output/$gene_name.header
	echo -e "Sample_ID#Type" > $output/$gene_name.raw.header
	awk NF $output/$gene_name.raw.header

	for ((count=1;count<=$num_miRs;count++));
        do
		cat $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]} | cut -f3 | head -1 >> $output/$gene_name.header
		cat $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]} | sort -k1 | cut -f4 > $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.tmp

		cat $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw | cut -f3 | head -1 >> $output/$gene_name.raw.header
                cat $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw | sort -k1 | cut -f4 > $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw.tmp
	done
	
	cat $output/$gene_name.header | tr "\n" "\t" > $output/$gene_name.$subtype.matrix.txt
	echo -e "\n" >> $output/$gene_name.$subtype.matrix.txt
	cat $output/$gene_name.raw.header | tr "\n" "\t" > $output/$gene_name.$subtype.matrix.raw.txt
        echo -e "\n" >> $output/$gene_name.$subtype.matrix.raw.txt
	
	cat $input_dir/$gene_name.sample_names | tr "#" "\t" | sort -k1 > $input_dir/$gene_name.sample_names_sorted
	paste $input_dir/$gene_name.sample_names_sorted $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[1]}.${seedArray[1]}.tmp > $output/$gene_name.$subtype.matrix.tmp
	paste $input_dir/$gene_name.sample_names $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[1]}.${seedArray[1]}.raw.tmp > $output/$gene_name.$subtype.matrix.raw.tmp
	for ((count=2;count<=$num_miRs;count++));
        do
		paste $output/$gene_name.$subtype.matrix.tmp $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.tmp > $output/$gene_name.$subtype.matrix.tmp1
		paste $output/$gene_name.$subtype.matrix.raw.tmp $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw.tmp > $output/$gene_name.$subtype.matrix.raw.tmp1
		mv $output/$gene_name.$subtype.matrix.tmp1 $output/$gene_name.$subtype.matrix.tmp
		 mv $output/$gene_name.$subtype.matrix.raw.tmp1 $output/$gene_name.$subtype.matrix.raw.tmp
	done
	
	cat $output/$gene_name.$subtype.matrix.tmp >> $output/$gene_name.$subtype.matrix.txt
	 cat $output/$gene_name.$subtype.matrix.raw.tmp >> $output/$gene_name.$subtype.matrix.raw.txt

	rm $output/$gene_name.header $output/$gene_name.raw.header $input_dir/scores.all_samples.FIMO.$gene_name.*tmp $input_dir/$gene_name.sample_names_sorted $output/$gene_name.$subtype.matrix.tmp $input_dir/$gene_name.sample_names $output/$gene_name.$subtype.matrix.raw.tmp
	rm $output/$gene_name.miR_info.txt

	awk NF $output/$gene_name.$subtype.matrix.txt > $output/$gene_name.$subtype.matrix.tmp
	mv $output/$gene_name.$subtype.matrix.tmp $output/$gene_name.$subtype.matrix.txt
	
	awk NF $output/$gene_name.$subtype.matrix.raw.txt > $output/$gene_name.$subtype.matrix.raw.tmp
        mv $output/$gene_name.$subtype.matrix.raw.tmp $output/$gene_name.$subtype.matrix.raw.txt

	for ((count=1;count<=$num_miRs;count++));
       do
		rm $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]} 
		rm $input_dir/scores.all_samples.FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw
		rm $input_dir/scores.summarized.$gene_name.${miRArray[$count]}.${seedArray[$count]}
		rm $input_dir/FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}
		rm $input_dir/FIMO.$gene_name.${miRArray[$count]}.${seedArray[$count]}.raw
     done
	for j in `cat $input_dir/${subtype}_paired_Samples | cut -f2 | sed 's/#/./g'`
	do
		for ((count=1;count<=$num_miRs;count++));
		do
			rm $input_dir/$j/${MREArray[$count]}.$gene_name.${miRArray[$count]}.${seedArray[$count]}.FIMO_MRE
			rm $input_dir/$j/${MREArray[$count]}.$gene_name.${miRArray[$count]}.${seedArray[$count]}.FIMO_MRE_Counts
			rm $input_dir/$j/${MREArray[$count]}.$gene_name.${miRArray[$count]}.${seedArray[$count]}.FIMO_raw_MRE_Counts
		done
	done

	echo `date`

fi
