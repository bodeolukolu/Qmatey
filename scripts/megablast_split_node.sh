


echo -e "\e[97m########################################################\n \e[38;5;210mQmatey MegaBLAST \n\e[97m########################################################\n"

export local_db=$( echo $local_db | awk '{gsub(/,/," ")}1' )

if (echo $local_db | grep -q 'nt'); then
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
			mv $ccf /tmp/Qmatey_multi_node/metagenome/alignment/$ccf
			cd /tmp/Qmatey_multi_node/metagenome/alignment
			awk -v rpm=$rpm 'NR%rpm==1{close("subfile"i); i++}{print > "subfile"i}' $ccf & PIDsplit2=$!
			wait $PIDsplit2
			for sub in $(ls subfile* | sort -V); do (
				if [[ "$taxids" == true ]]; then
					${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/blastn -task megablast -query $sub -db "${local_db}" -num_threads 1 -perc_identity $percid -max_target_seqs $max_target \
					-qcov_hsp_perc $qcov -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out ${sub}_out.blast &&
					wait
				else
					${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/blastn -task megablast -query $sub -db "${local_db}" -num_threads 1 -perc_identity $percid -max_target_seqs $max_target \
					-qcov_hsp_perc $qcov -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out ${sub}_out.blast &&
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

  if [[ "$blast_location" =~ "custom" ]]; then
  	echo -e "${YELLOW}- performing custom BLAST in multi-node mode"
    cd ${projdir}/metagenome/haplotig/splitccf/splitccf_node${njob}
		for ccf in $(ls * | sort -V); do
			mv $ccf /tmp/Qmatey_multi_node/metagenome/alignment/$ccf
			cd /tmp/Qmatey_multi_node/metagenome/alignment
			awk -v rpm=$rpm 'NR%rpm==1{close("subfile"i); i++}{print > "subfile"i}' $ccf & PIDsplit2=$!
			wait $PIDsplit2
			for sub in $(ls subfile* | sort -V); do (
				if [[ "$taxids" == true ]]; then
					${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/blastn -task megablast -query $sub -db "${custom_db}" -num_threads 1 -perc_identity $percid -max_target_seqs $max_target \
					-qcov_hsp_perc $qcov -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out ${sub}_out.blast &&
					wait
				else
					${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/blastn -task megablast -query $sub -db "${custom_db}" -num_threads 1 -perc_identity $percid -max_target_seqs $max_target \
					-qcov_hsp_perc $qcov -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out ${sub}_out.blast &&
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
