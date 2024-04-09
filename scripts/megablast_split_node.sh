


echo -e "\e[97m########################################################\n \e[38;5;210mQmatey MegaBLAST \n\e[97m########################################################\n"

export local_db=$( echo $local_db | awk '{gsub(/,/," ")}1' )

if (echo $local_db | grep -q 'nt'); then
	if [[ -z $percid ]]; then
		export percid=95
	fi
fi
if (echo $local_db | grep -q 'refseq'); then
	if [[ -z $percid ]]; then
		export percid=95
	fi
fi
if (echo $local_db | grep -q '16S') || (echo $local_db | grep -q '18S') || (echo $local_db | grep -q '28S') || (echo $local_db | grep -q 'ITS'); then
	if [[ -z $percid ]]; then
		export percid=95
	fi
fi
if (echo $local_db | grep -q '16s') || (echo $local_db | grep -q '18s') || (echo $local_db | grep -q '28s') || (echo $local_db | grep -q 'ITs'); then
	if [[ -z $percid ]]; then
		export percid=95
	fi
fi
if [[ "$blast_location" == "custom" ]]; then
	if [[ -z $percid ]]; then
		export percid=95
	fi
fi
export rpm=$((reads_per_megablast * 2))


if [[ $nodes -gt 1 ]]; then
	rm -rf /tmp/Qmatey_multi_node 2> /dev/null
	mkdir -p /tmp/Qmatey_multi_node/metagenome/alignment
fi


blast () {

  if [[ "$blast_location" =~ "local" ]]; then
    echo -e "${YELLOW}- performing a local BLAST in multi-node mode"
    cd ${projdir}/metagenome/haplotig/splitccf/splitccf_node${njob}
		for ccf in $(ls * | sort -V); do
			mv $ccf /tmp/Qmatey_multi_node/metagenome/alignment/
			cd /tmp/Qmatey_multi_node/metagenome/alignment
			awk -v pat="$ccf" -v rpm="$rpm" 'NR%rpm==1{close(pat"_subfile"pat"_"i); i++}{print > pat"_subfile"pat"_"i}' $ccf & PIDsplit2=$!
			wait $PIDsplit2
			if [[ "$(ls *subfile* | wc -l)" -gt "$threads" ]]; then
				Nsize=$(ls *subfile* | wc -l)
				nbatch=$(($Nsize / $threads))
				nsize=$(($nbatch * $threads))
				for subf in $(ls *subfile* | sort -V | head -n $nsize); do mv $subf ${subf#*_}; done
			fi
			if [[ "$(ls *subfile* | wc -l)" -le "$threads" ]]; then
				for subf in $(ls *subfile* | sort -V); do mv $subf ${subf#*_}; done
			fi
			for sub in $(ls subfile* | sort -V); do (
				if [[ "$taxids" == true ]]; then
					${Qmatey_dir}/tools/ncbi-blast-2.15.0+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
					-taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
					wait
					if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
					wait
				else
					${Qmatey_dir}/tools/ncbi-blast-2.15.0+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
					-outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
					wait
					if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
					wait
				fi
				wait
				gzip ${sub}_out.blast &&
				rm $sub )&
				if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
					wait
				fi
			done
			wait
			for subfile in *_out.blast.gz; do
				cat $subfile >> ${ccf}.blast.gz
				rm $subfile
			done
			wait
			cat ${ccf}.blast.gz >> combined_compressed_node${njob}.megablast.gz &&
			rm ${ccf}.blast.gz; rm $ccf &&
			cd ${projdir}/metagenome/haplotig/splitccf/splitccf_node${njob}
		done
		wait

		cd /tmp/Qmatey_multi_node/metagenome/alignment
		if [[ "$(ls *subfile* | wc -l)" -gt 0 ]]; then
			for subf in $(ls *subfile* | sort -V); do mv $subf ${subf#*_}; done
			for sub in $(ls subfile* | sort -V); do (
				if [[ "$taxids" == true ]]; then
					${Qmatey_dir}/tools/ncbi-blast-2.15.0+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
					-taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
					wait
					if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
					wait
				else
					${Qmatey_dir}/tools/ncbi-blast-2.15.0+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
					-outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
					wait
					if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
					wait
				fi
				wait
				gzip ${sub}_out.blast &&
				rm $sub )&
				if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
					wait
				fi
			done
			wait
			for subfile in *_out.blast.gz; do
				cat $subfile >> ${ccf}.blast.gz
				rm $subfile
			done
			wait
			cat ${ccf}.blast.gz >> combined_compressed_node${njob}.megablast.gz &&
			rm ${ccf}.blast.gz &&
			wait
		fi
		cd ${projdir}/metagenome/haplotig/splitccf/splitccf_node${njob}
		wait
    mv /tmp/Qmatey_multi_node/metagenome/alignment/combined_compressed_node${njob}.megablast.gz ${projdir}/metagenome/alignment/
  fi

  if [[ "$blast_location" =~ "custom" ]]; then
  	echo -e "${YELLOW}- performing custom BLAST in multi-node mode"
    cd ${projdir}/metagenome/haplotig/splitccf/splitccf_node${njob}
		for ccf in $(ls * | sort -V); do
			mv $ccf /tmp/Qmatey_multi_node/metagenome/alignment/
			cd /tmp/Qmatey_multi_node/metagenome/alignment
			awk -v rpm=$rpm 'NR%rpm==1{close("subfile"i); i++}{print > "subfile"i}' $ccf & PIDsplit2=$!
			wait $PIDsplit2
			for sub in $(ls subfile* | sort -V); do (
				if [[ "$taxids" == true ]]; then
					${Qmatey_dir}/tools/ncbi-blast-2.15.0+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
					-taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
					wait
					if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
					wait
				else
					${Qmatey_dir}/tools/ncbi-blast-2.15.0+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
					-outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
					wait
					if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
					wait
				fi
				wait
				gzip ${sub}_out.blast &&
				rm $sub )&
				if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
					wait
				fi
			done
			wait
			for subfile in *_out.blast.gz; do
				cat $subfile >> ${ccf}.blast.gz
				rm $subfile
			done
			wait
			cat ${ccf}.blast.gz >> combined_compressed_node${njob}.megablast.gz &&
			rm ${ccf}.blast.gz; rm $ccf &&
			cd ${projdir}/metagenome/haplotig/splitccf/splitccf_node${njob}
		done
    mv /tmp/Qmatey_multi_node/metagenome/alignment/combined_compressed_node${njob}.megablast.gz ${projdir}/metagenome/alignment/
  fi
}

time blast &>> ${projdir}/log_node${njob}.out
wait
touch ${projdir}/megablast_done_node${njob}.txt
