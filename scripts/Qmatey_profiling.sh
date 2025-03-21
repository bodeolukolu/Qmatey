if R --version 2> /dev/null; then
	:
else
	module add R
fi

#cd "${projdir}"
#free=$(df ~ | awk 'NR==2{print $4}')
#required=$(du -s . | awk '{print $1}')
#required=$((required * 2))
#if [[ "$free" < "$required" ]]; then
#	echo -e "${magenta}- You might not have enough disk space. Free: $((free/1000000))G; Required(approx. 1.2x-2x the size of fastq files): $((required/1000000)G  ${white}\n"
#	read -p "- Do you wish to continue running Qmatey? " -n 1 -r
#	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
#		printf '\n'
#		echo -e "${magenta}- Exiting Qmatey ${white}\n"
#		sleep 5 && exit 0
#	fi
#fi

# Create custom database from user-provided fasta file and file mapping identifiers/header to taxids
if [[ "$blast_location" == "custom" ]]; then
	export taxids=false
	custom_db=${input_dbfasta##*/}
	custom_db=${custom_db%.f*}
	if test ! -f db_created.txt || test ! -d $custom_db; then
		if [[ -z $input_dbfasta ]] || test ! -f $input_dbfasta; then
			echo -e "${magenta}- \n- Please provide fasta file(s) required to create custom database or specify correct directory ${white}\n"
			echo -e "${magenta}- \n- Please provide files mapping fasta sequences to taxid(s) or specify correct directory ${white}\n" > ${projdir}/job_killed
			exit 0
		fi
		if [[ -z $map_taxids ]] || test ! -f $map_taxids; then
			echo -e "${magenta}- \n- Please provide files mapping fasta sequences to taxid(s) or specify correct directory ${white}\n"
			echo -e "${magenta}- \n- Please provide files mapping fasta sequences to taxid(s) or specify correct directory ${white}\n" > ${projdir}/job_killed
			exit 0
		fi
		cd "${projdir}"
		echo -e "$1 \e[31m Creating custom database"
		# cp ${input_dbfasta} $custom_db
		# cp ${map_taxids} $custom_db
		mkdir -p "$custom_db"
		# cd $custom_db
		# custom_db_dir=${custom_db}/${custom_db}
		if [[ $(file $input_dbfasta 2> /dev/null | awk -F' ' '{print $2}') == gzip ]]; then
			export title_db=${custom_db}
			export custom_db=${custom_db}
			zcat $input_dbfasta | ${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/makeblastdb -in - -out ${custom_db} -title $title_db -parse_seqids -blastdb_version 5 -taxid_map $map_taxids -dbtype nucl
			mv ${custom_db}.* ${custom_db}/ &&
			mv ${custom_db}/$input_dbfasta ${projdir}/ 2> /dev/null
			# rm ${input_dbfasta##*/} ${map_taxids##*/}
		else
			${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/makeblastdb -in $input_dbfasta -out $custom_db -title $title_db -parse_seqids -blastdb_version 5 -taxid_map $map_taxids -dbtype nucl
			mv ${custom_db}.* ${custom_db}/
			mv ${custom_db}/${input_dbfasta##*/} ${projdir}/
			# rm ${input_dbfasta##*/} ${map_taxids##*/}
		fi
		# cd "${projdir}"
		touch db_created.txt
		export custom_db=${projdir}/${custom_db}/${custom_db}
		# export custom_db=${projdir}/${custom_db}/${custom_db}/${custom_db}
	else
		export custom_db=${projdir}/${custom_db}/${custom_db}
	fi
	echo custom_db: $custom_db
fi



cd "${projdir}"
if [[ -z "$subsample_shotgun_R1" ]] || [[ "$subsample_shotgun_R1" == "true" ]]; then
	subsample_shotgun_R1=ATGCAT
fi
if [[ -z "$subsample_shotgun_R2" ]] || [[ "$subsample_shotgun_R2" == "true" ]]; then
	subsample_shotgun_R2=CATG
fi
if [[ -z "$HDsubsample" ]]; then
	HDsubsample=false
fi

if [[ -z "$simulation_lib" ]]; then
	simulation_lib=complete_digest
fi
if [[ -z "$fragment_size_range" ]]; then
	export fragment_size_range=64,600
fi
if [[ -z "$gcov" ]]; then
	export gcov=3
fi
if [[ -z "$max_read_length" ]]; then
	export max_read_length=150
fi
if [[ -z "$fastMegaBLAST" ]]; then
	export fastMegaBLAST=true
fi
if [[ -z "$taxonomic_level" ]]; then
	export taxonomic_level=strain,species,genus,family,order,class,phylum
fi
if [[ $taxonomic_level =~ strain ]]; then
	export strain_level=true
fi
if [[ $taxonomic_level =~ species ]]; then
	export species_level=true
fi
if [[ $taxonomic_level =~ genus ]]; then
	export genus_level=true
fi
if [[ $taxonomic_level =~ family ]]; then
	export family_level=true
fi
if [[ $taxonomic_level =~ order ]]; then
	export order_level=true
fi
if [[ $taxonomic_level =~ class ]]; then
	export class_level=true
fi
if [[ $taxonomic_level =~ phylum ]]; then
	export phylum_level=true
fi
if [[ "$genome_scaling" ]]; then
	export genome_scaling=false
fi
if [[ -z "$zero_inflated" ]]; then
	export zero_inflated=0.01
fi
if [[ -z "$exclude_rRNA" ]]; then
	export exclude_rRNA=false
fi
if [[ -z "$node" ]]; then
	export node=1
fi
if [[ -z $sunburst_taxlevel ]]; then
	export run_sunburst=false
else
	export run_sunburst=true
fi
if [[ -z $CCLasso_corr ]]; then
	export run_corr=false
else
	export run_corr=true
fi
if [[ -z $spearman_corr ]]; then
	export run_Spearman_corr=false
else
	export run_Spearman_corr=true
fi
if [[ "$minRD" -gt 0 ]]; then
	export zminRD=false
fi
if [[ -z $minRD ]]; then
	export minRD=1
	export zminRD=true
fi
if [[ -z $fullqlen_alignment ]]; then
	export fullqlen_alignment=false
fi
if [[ -z $burn_in_diagnostic_reads ]]; then
	export burn_in_diagnostic_reads=false
fi
if [[ -z $min_strain_uniq ]]; then
	export min_strain_uniq=1
fi
if [[ -z $maxindel ]]; then
	export maxindel=100
fi
if [[ -z $max_target ]]; then
	export max_target=5000
fi
if [[ -z "$cross_taxon_validation" ]]; then
	export cross_taxon_validation=true
fi
if [[ -z "$annotate_seq" ]]; then
	export annotate_seq=false
fi

mkdir -p "${projdir}"/tmp
export TMPDIR="${projdir}"/tmp

if [[ "$library_type" =~ "RRS" ]] || [[ "$library_type" =~ "rrs" ]] || [[ "$library_type" == "WGS" ]] || [[ "$library_type" == "wgs" ]] || [[ "$library_type" == "SHOTGUN" ]] || [[ "$library_type" == "shotgun" ]]; then
	if [[ -z $reads_per_megablast_burn_in ]]; then
		export reads_per_megablast=10000
	fi
	if [[ $reads_per_megablast_burn_in -eq 10000 ]] && [[ -z $reads_per_megablast ]]; then
		export reads_per_megablast=10000
	fi
	if [[ -z $reads_per_megablast ]]; then
		export reads_per_megablast=1000
	fi
fi

if [[ "$library_type" =~ "amplicon" ]] || [[ "$library_type" =~ "Amplicon" ]] || [[ "$library_type" =~ "AMPLICON" ]] || [[ "$library_type" =~ "16S" ]] || [[ "$library_type" =~ "16s" ]]|| [[ "$library_type" =~ "ITS" ]] || [[ "$library_type" =~ "its" ]]; then
  export min_strain_uniq=1
	export exclude_rRNA=false
	if [[ -z $reads_per_megablast_burn_in ]]; then
		export reads_per_megablast=10000
	fi
	if [[ -z $reads_per_megablast ]]; then
		export reads_per_megablast=100
	fi
fi

if [[ "$blast_location" == "custom" ]]; then
	if [[ $reads_per_megablast_burn_in -eq 10000 ]]; then
		if [[ -z $reads_per_megablast ]]; then
			export reads_per_megablast=1000
		fi
	fi
	if [[ -z $reads_per_megablast ]]; then
		export reads_per_megablast=1000
	fi
fi

if [[ -z $min_percent_sample ]]; then
	export min_percent_sample=5,10,20
fi
if [[ -z $min_pos_corr ]]; then
	export min_pos_corr=0.1,0.2,0.3
fi
if [[ -z $max_neg_corr ]]; then
	export max_neg_corr=0.1,0.2,0.3
fi
if [[ -z $sunburst_nlayers ]]; then
	export sunburst_nlayers=phylum,genus,species
fi




if [[ -z "$(ls -A ./samples)" ]]; then
	if [[ -d ./simulate_genomes ]]; then
		:
	else
		echo -e "$1 \e[31m samples folder is empty, Qmatey will exit"
		echo -e "$1 \e[31m samples folder is empty, Qmatey will exit" > ${projdir}/job_killed
		sleep 10; exit 0
	fi
else
	:
fi

totalk=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)
export loopthread=2
if [[ "$threads" -gt 1 ]]; then
	export N=$((threads/2))
	export ram1=$(($totalk/$N))
else
	export N=1 && loopthread=threads
fi
export ram1=$((ram1/1000000))
export Xmx1=-Xmx${ram1}G
export ram2=$(echo "$totalk*0.00000095" | bc)
export ram2=${ram2%.*}
export Xmx2=-Xmx${ram2}G

if [[ -z "$threads" ]]; then
	export threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		export threads=$((threads-2))
	fi
fi
if  [[ "$threads" -ge 1 ]]; then
	export loopthread=2
	export N=$(($threads/2))
else
	export N=1 && export loopthread=$threads
fi

cd ${Qmatey_dir}/tools
if test -f "rankedlineage.dmp.gz"; then
	gunzip rankedlineage.dmp.gz
fi
if [[ "$threads" -le 4 ]]; then
	export gthreads=threads
	export Xmxg=$Xmx2
	export gN=1
else
	export gthreads=4
	export gN=$(( threads / gthreads ))
	export ramg=$(( ram2 / gN ))
	export Xmxg=-Xmx${ramg}G
fi


#extract_parameters_above_for_megablast_splitrun_multi_node_mode
##################################################################################################################
#Create all necessary subdirectories for downstream processing
cd "${projdir}"
if [[ -d metagenome_ref_normalize ]]; then
	if [[ ! -d metagenome_no_normalize ]]; then
		mv "${projdir}"/metagenome_ref_normalize ${projdir}/metagenome
	fi
fi
if [[ -d metagenome_no_normalize ]]; then
	if [[ ! -d metagenome_ref_normalize ]]; then
		mv "${projdir}"/metagenome_no_normalize ${projdir}/metagenome
	fi
fi
if [[ -d metagenome_ref_normalize ]]; then
	if [[ -d metagenome_no_normalize ]]; then
		echo -e "${magenta}- Reference-normalized and non-normalized directories both exist ${white}\n"
		echo -e "${magenta}- Manually rename 1 the 2 existing metagenome-directories to metagenome ${white}\n"
		echo -e "${magenta}- Reference-normalized and non-normalized directories both exist ${white}\n" > ${projdir}/job_killed
		echo -e "${magenta}- Manually rename 1 the 2 existing metagenome-directories to metagenome ${white}\n" >> ${projdir}/job_killed
		sleep 10; exit 0
	fi
fi
if [[ ! -d metagenome ]]; then
	mkdir metagenome
fi
cd metagenome
mkdir -p haplotig
mkdir -p alignment
mkdir -p sighits
mkdir -p results
mkdir -p ./results/ref_aligned_summaries

echo -e "\e[97m########################################################\n \e[38;5;210m Generating synthetic/mock community sequences for simulation \n\e[97m########################################################\n"
simulate_reads () {
	cd "${projdir}"/simulate_genomes/
	for simdir in */ ; do
		cd "$simdir" && gunzip ./*.gz 2> /dev/null
		awk '{ sub("\r$", ""); print }' abundance.txt | awk '{gsub(/ /,"_"); gsub(/.fasta.gz/,""); gsub(/.fasta/,""); gsub(/.fna/,"");}1' | sed '$ s/$//' > abundance.tmp &&  mv abundance.tmp abundance.txt
		wait
		while IFS="" read -r p || [ -n "$p" ]; do (
			endt=$(echo "$p" | awk '{print $1}')
			cat "$(echo $p | awk '{print $2}')".fasta | awk -v endt="$endt" '{a[NR]=$0}END{for (i=0; i<endt; i++){for(k in a){print a[k]}}}' |\
			gzip > ../taxa_"$(echo $p | awk '{print $2}')".fa.gz  ) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done < abundance.txt
		wait
		cat ../taxa_*.fa.gz >> ../"${simdir%/*}".fasta.gz &&
		rm ../taxa_*.fa.gz &&
		wait
		cd ../
	done
	wait

	minfrag=${fragment_size_range%,*}
	maxfrag=${fragment_size_range#*,}
	if [[ "$simulation_lib"  =~ "complete_digest" ]] || [[ "$simulation_lib"  =~ "partial_digest" ]]; then
		if [[ "$simulation_motif_R1" ]]; then
			echo $simulation_motif_R1 | awk '{gsub(/,/,"\n");}1' | awk '{print "RE"NR"\t"$1}' > REnase_R1.txt
			RE1a=$(grep 'RE1' REnase_R1.txt | awk '{print $2}')
			RE1b=$(grep 'RE2' REnase_R1.txt | awk '{print $2}')
			RE1c=$(grep 'RE3' REnase_R1.txt | awk '{print $2}')
			RE1d=$(grep 'RE4' REnase_R1.txt | awk '{print $2}')
			if [[ -z "$RE1b" ]]; then RE1b=999; fi
			if [[ -z "$RE1c" ]]; then RE1c=999; fi
			if [[ -z "$RE1d" ]]; then RE1d=999; fi
			rm REnase_R1.txt
			if [[ "$simulation_motif_R2" ]]; then
			  echo $simulation_motif_R2 | awk '{gsub(/,/,"\n");}1' | awk '{print "RE"NR"\t"$1}' > REnase_R2.txt
			  RE2a=$(grep 'RE1' REnase_R2.txt | awk '{print $2}')
			  RE2b=$(grep 'RE2' REnase_R2.txt | awk '{print $2}')
			  RE2c=$(grep 'RE3' REnase_R2.txt | awk '{print $2}')
			  RE2d=$(grep 'RE4' REnase_R2.txt | awk '{print $2}')
			  if [[ -z "$RE2b" ]]; then RE2b=999; fi
			  if [[ -z "$RE2c" ]]; then RE2c=999; fi
			  if [[ -z "$RE2d" ]]; then RE2d=999; fi
			  rm REnase_R2.txt
			else
			  RE2a=999
			  RE2b=999
			  RE2c=999
			  RE2d=999
			fi
		fi
	fi

	for unsim in *.fasta.gz; do
		if [[ "$simulation_lib" =~ "complete_digest" ]]; then
			awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' <(zcat "${unsim}") | gzip > hold0_"${unsim}"
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' <(zcat ./hold0_"${unsim}") | awk -F"\t" '{print $2}' | awk '{gsub(/a/,"A");gsub(/c/,"C");gsub(/g/,"G");gsub(/t/,"T");}1' | shuf | gzip > hold1_"${unsim}" &&
			rm hold0_"${unsim}"
			end="$(awk '{ if ( length > L ) { L=length} }END{ print L}' <(zcat hold1_${unsim}))"
			for g in $(seq $gcov); do
				for (( gline=1; gline <= $end; gline+=100000 )); do
					awk -v pat1=$gline -v pat2=$((gline+99999)) 'NR >= pat1 && NR <= pat2' <(zcat hold1_${unsim}) > hold2_${unsim%.gz} &&
					cutpos=$(shuf -i 500000-1500000 -n1)
					while [[ "$(awk '{ if ( length > L ) { L=length} }END{ print L}' hold2_${unsim%.gz})" -gt 50000 ]] || [[ "$cutpos" -gt 50000 ]]; do
						fold -w "$cutpos" hold2_${unsim%.gz} > hold2_${unsim%.gz}.tmp &&
						mv hold2_${unsim%.gz}.tmp hold2_${unsim%.gz} &&
						cutpos=$((cutpos / 2)) &&
						cutpos=$(awk -v minfrag=5000 -v cutpos=$cutpos 'BEGIN{srand();print int(rand()*((cutpos+5000)-(cutpos-5000)))+(cutpos-5000) }')
					done
					gzip -c hold2_${unsim%.gz} >> ${unsim}.tmp && rm hold2_${unsim%.gz}
				done
				zcat ${unsim}.tmp | grep -v '>' | awk -v RE1a="$RE1a" -v RE1b="$RE1b" -v RE1c="$RE1c" -v RE1d="$RE1d" -v RE2a="$RE2a" -v RE2b="$RE2b" -v RE2c="$RE2c" -v RE2d="$RE2d" \
				'{gsub(RE1a,RE1a"\n"RE1a); gsub(RE1b,RE1b"\n"RE1b); gsub(RE1c,RE1c"\n"RE1c); gsub(RE1d,RE1d"\n"RE1d); \
				gsub(RE2a,RE2a"\n"RE2a); gsub(RE2b,RE2b"\n"RE2b); gsub(RE2c,RE2c"\n"RE2c); gsub(RE2d,RE2d"\n"RE2d); }1' > ${unsim}.tmp1.txt &&
				sleep 2
				if [[ "$HDsubsample" == true ]]; then
					grep "^$RE1a.*$RE1a$\|^$RE1b.*$RE1b$\|^$RE1c.*$RE1c$\|^$RE1d.*$RE1d$\|^$RE2a.*$RE2a$\|^$RE2b.*$RE2b$\|^$RE2c.*$RE2c$\|^$RE2d.*$RE2d$" ${unsim}.tmp1.txt 2> /dev/null >> ${unsim}.tmp.txt &&
					sleep 2
				fi
				grep "^$RE1a.*$RE2a$\|^$RE1a.*$RE2b$\|^$RE1a.*$RE2c$\|^$RE1a.*$RE2d$\|^$RE1b.*$RE2a$\|^$RE1b.*$RE2b$\|^$RE1b.*$RE2c$\|^$RE1b.*$RE2d$" ${unsim}.tmp1.txt 2> /dev/null >> ${unsim}.tmp.txt &&
				sleep 2
				grep "^$RE1c.*$RE2a$\|^$RE1c.*$RE2b$\|^$RE1c.*$RE2c$\|^$RE1c.*$RE2d$\|^$RE1d.*$RE2a$\|^$RE1d.*$RE2b$\|^$RE1d.*$RE2c$\|^$RE1d.*$RE2d$" ${unsim}.tmp1.txt 2> /dev/null >> ${unsim}.tmp.txt &&
				sleep 2
				grep "^$RE2a.*$RE1a$\|^$RE2a.*$RE1b$\|^$RE2a.*$RE1c$\|^$RE2a.*$RE1d$\|^$RE2b.*$RE1a$\|^$RE2b.*$RE1b$\|^$RE2b.*$RE1c$\|^$RE2b.*$RE1d$" ${unsim}.tmp1.txt 2> /dev/null >> ${unsim}.tmp.txt &&
				sleep 2
				grep "^$RE2c.*$RE1a$\|^$RE2c.*$RE1b$\|^$RE2c.*$RE1c$\|^$RE2c.*$RE1d$\|^$RE2d.*$RE1a$\|^$RE2d.*$RE1b$\|^$RE2d.*$RE1c$\|^$RE2d.*$RE1d$" ${unsim}.tmp1.txt 2> /dev/null >> ${unsim}.tmp.txt &&
				sleep 2
				rm "${unsim}" &&
				awk '{ print length"\t"$1}' ${unsim}.tmp.txt | awk -v minfrag=$minfrag 'BEGIN{OFS="\t"} {if ($1 >= minfrag) {print $0}}' | \
				awk -v maxfrag=$maxfrag 'BEGIN{OFS="\t"} {if ($1 <= maxfrag) {print $0}}' | awk '{print ">read"NR"_"$1"\t"$2}' | $gzip > gcov${g}_${unsim} &&
				rm -rf "*tmp*" 2> /dev/null
				wait
			done
			wait
			cat gcov*_${unsim} > ${unsim}
			rm gcov*_${unsim}
			rm hold*  2> /dev/null
		fi

		if [[ "$simulation_lib" =~ "partial_digest" ]]; then
			awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' <(zcat "${unsim}") | gzip > ./hold0_"${unsim}"
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' <(zcat ./hold0_"${unsim}") | awk -F"\t" '{print $2}' | awk '{gsub(/a/,"A");gsub(/c/,"C");gsub(/g/,"G");gsub(/t/,"T");}1' | shuf | gzip > ./hold1_"${unsim}" &&
			end="$(awk '{ if ( length > L ) { L=length} }END{ print L}' <(zcat ./hold1_${unsim}))"
			rm hold0_"${unsim}"
			for g in $(seq $gcov); do
				for (( gline=1; gline <= $end; gline+=100000 )); do
					awk -v pat1=$gline -v pat2=$((gline+99999)) 'NR >= pat1 && NR <= pat2' <(zcat hold1_${unsim}) > hold2_${unsim%.gz} &&
					cutpos=$(shuf -i 500000-1500000 -n1)
					while [[ "$(awk '{ if ( length > L ) { L=length} }END{ print L}' hold2_${unsim%.gz})" -gt $maxfrag ]] || [[ "$cutpos" -gt $maxfrag ]]; do
						fold -w "$cutpos" hold2_${unsim%.gz} > hold2_${unsim%.gz}.tmp &&
						mv hold2_${unsim%.gz}.tmp hold2_${unsim%.gz} &&
						cutpos=$((cutpos / 2)) &&
						cutpos=$(awk -v minfrag=$minfrag -v cutpos=$cutpos 'BEGIN{srand();print int(rand()*((cutpos+$minfrag)-(cutpos-$minfrag)))+(cutpos-$minfrag) }')
					done
					gzip -c hold2_${unsim%.gz} >> ${unsim}.tmp && rm hold2_${unsim%.gz}
				done
				wait

				zcat ${unsim}.tmp | grep -v '>' | sed 's/[^'"$RE1a"']*\('"$RE1a"'.*\)/\1/' | sed 's!'"$RE1a"'[^'"$RE1a"']*$!'"$RE1a"'!' | \
				sed 's/[^'"$RE1b"']*\('"$RE1b"'.*\)/\1/' | sed 's!'"$RE1b"'[^'"$RE1b"']*$!'"$RE1b"'!' | \
				sed 's/[^'"$RE1c"']*\('"$RE1c"'.*\)/\1/' | sed 's!'"$RE1c"'[^'"$RE1c"']*$!'"$RE1c"'!' | \
				sed 's/[^'"$RE1d"']*\('"$RE1d"'.*\)/\1/' | sed 's!'"$RE1d"'[^'"$RE1d"']*$!'"$RE1d"'!' | \
				sed 's/[^'"$RE2a"']*\('"$RE2a"'.*\)/\1/' | sed 's!'"$RE2a"'[^'"$RE2a"']*$!'"$RE2a"'!' | \
				sed 's/[^'"$RE2b"']*\('"$RE2b"'.*\)/\1/' | sed 's!'"$RE2b"'[^'"$RE2b"']*$!'"$RE2b"'!' | \
				sed 's/[^'"$RE2c"']*\('"$RE2c"'.*\)/\1/' | sed 's!'"$RE2c"'[^'"$RE2c"']*$!'"$RE2c"'!' | \
				sed 's/[^'"$RE2d"']*\('"$RE2d"'.*\)/\1/' | sed 's!'"$RE2d"'[^'"$RE2d"']*$!'"$RE2d"'!' | \
				awk -v maxfrag=$maxfrag '{print substr($0,1,maxfrag)}' | awk '{print length"\t"$1}' | \
				awk -v minfrag=$minfrag 'BEGIN{OFS="\t"} {if ($1 >= minfrag) {print $0}}' | awk '{print ">read"NR"_"$1"\t"$2}' | $gzip > gcov${g}_${unsim} &&
				rm *.tmp
			done
			wait
			cat gcov*_${unsim} > ${unsim}
			rm gcov*_${unsim}
			rm hold*
		fi

		if [[ "$simulation_lib" =~ "shotgun" ]]; then
			awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' <(zcat "${unsim}") | gzip > ./hold0_"${unsim}"
			awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' <(zcat ./hold0_"${unsim}") | awk -F"\t" '{print $2}' | awk '{gsub(/a/,"A");gsub(/c/,"C");gsub(/g/,"G");gsub(/t/,"T");}1' | shuf | gzip> ./hold1_"${unsim}" &&
			rm hold0_"${unsim}"
			end="$(awk '{ if ( length > L ) { L=length} }END{ print L}' <(zcat ./hold1_${unsim}))"
			for g in $(seq $gcov); do
				for (( gline=1; gline<=$end; gline+=100000 )); do
					awk -v pat1=$gline -v pat2=$((gline+99999)) 'NR >= pat1 && NR <= pat2' <(zcat hold1_${unsim}) > hold2_${unsim%.gz} &&
					cutpos=$(shuf -i 500000-1500000 -n1)
					while [[ "$(awk '{ if ( length > L ) { L=length} }END{ print L}' hold2_${unsim%.gz})" -gt "$maxfrag" ]] || [[ "$cutpos" -gt "$maxfrag" ]]; do
						fold -w "$cutpos" hold2_${unsim%.gz} > hold2_${unsim%.gz}.tmp &&
						mv hold2_${unsim%.gz}.tmp hold2_${unsim%.gz} &&
						cutpos=$((cutpos / 2)) &&
						cutpos=$(awk -v minfrag=$minfrag -v cutpos=$cutpos 'BEGIN{srand();print int(rand()*((cutpos+minfrag)-(cutpos-minfrag)))+(cutpos-minfrag) }')
					done
					gzip -c hold2_${unsim%.gz} >> ${unsim}.tmp && rm hold2_${unsim%.gz}
				done

				zcat ${unsim}.tmp | grep -v '>' | awk -v maxfrag=$maxfrag '{print substr($0,1,maxfrag)}' | awk '{print length"\t"$1}' | \
				awk -v minfrag=$minfrag 'BEGIN{OFS="\t"} {if ($1 >= minfrag) {print $0}}' | awk '{print ">read"NR"_"$1"\t"$2}' | $gzip > gcov${g}_${unsim} &&
				rm *.tmp
			done
			cat gcov*_${unsim} > ${unsim}
			rm gcov*_${unsim}
			rm hold*
		fi
	done
	wait

	for unsim in *.fasta.gz; do
		awk '{print $2}' <(zcat ${unsim}) | awk -v len=$max_read_length '{print substr($0,1,len)}' | awk '{print length"\t"$1}' | \
		awk '{print ">read"NR"_fraglength"$1"\n"$2}' | $gzip > ../samples/${unsim%.fasta.gz}_R1.fasta.gz &&
		awk '{print $2}' <(zcat ${unsim}) | tr ACGTacgt TGCAtgca | rev | awk -v len=$max_read_length '{print substr($0,1,len)}' | awk '{print length"\t"$1}' | \
		awk '{print ">read"NR"_fraglength"$1"\n"$2}' | $gzip > ../samples/${unsim%.fasta.gz}_R2.fasta.gz &&
		wait
	done
	wait

	# cd "${projdir}"/simulate_genomes/
	# for simdir in */ ; do
	# 	mkdir -p refgenome
	# 	cd "$simdir"
	# 	echo -e "Taxa\tGenome_size_(bp)\tSequenced_(bp)\tPercent_Sequenced" > ../${simdir%/*}_taxa_Seq_Genome_Cov.txt
	# 	while IFS="" read -r p || [ -n "$p" ]; do
	# 		mockfile=$(echo $p | awk '{print $2}') &&
	# 		mockfile=${mockfile}.fasta &&
	# 		cp "$mockfile" ../refgenome &&
	# 		cd ../refgenome &&
	# 		$bwa index -a bwtsw "$mockfile" &&
	# 		$samtools faidx "$mockfile" &&
	# 		$java -jar $picard CreateSequenceDictionary REFERENCE=$mockfile    OUTPUT=${mockfile%.fasta}.dict &&
	# 		$bwa mem -t $threads $mockfile ../../samples/${simdir%/*}_R1.fasta.gz ../../samples/${simdir%/*}_R2.fasta.gz | $samtools view -bS - > ../${mockfile%.fasta}.bam &&
	# 		Taxa_gzfetch=$(grep -v '>' "$mockfile" | wc -c | awk '{print $1}') &&
	# 		Sequenced_fetch=$($samtools flagstat ../${mockfile%.fasta}.bam | grep '+ 0 mapped' | awk '{print $1}') &&
	# 		Perc=$(bc <<<"scale=3; $Sequenced_fetch*100/$Taxa_gzfetch") &&
	# 		printf "${mockfile%.fasta}\t$Taxa_gzfetch\t$Sequenced_fetch\t$Perc\n" | awk -F"\t" 'BEGIN{OFS="\t"}{gsub(/^./,"0.",$4);}1' | awk -F"\t" 'BEGIN{OFS="\t"}{gsub(/^0./,"0",$4);}1' >> ../${simdir%/*}_taxa_Seq_Genome_Cov.txt &&
	# 		rm -rf ./refgenome/* &&
	# 		rm ../"${mockfile%.fasta}".bam 2> /dev/null &&
	# 		cd ../$simdir
	# 		wait
	# 	done < abundance.txt
	# 	wait
	# 	rm -rf "${projdir}"/simulate_genomes/refgenome
	# 	awk '{gsub(/\.\./,".",$4);}1' ../${simdir%/*}_taxa_Seq_Genome_Cov.txt > ../${simdir%/*}_taxa_Seq_Genome_Cov.tmp &&
	# 	:> ../${simdir%/*}_taxa_Seq_Genome_Cov.txt &&
	# 	mv ../${simdir%/*}_taxa_Seq_Genome_Cov.tmp ../${simdir%/*}_taxa_Seq_Genome_Cov.txt
	# 	cd ../
	# done
}
cd "${projdir}"
if [[ "$simulate_reads" == 1 ]]; then
	if [[ "$(ls ./samples/*.f* | wc -l)" -gt 0 ]]; then
		echo -e "${magenta}- samples directory/folder contains fasta/fastq files. Delete samples directory/folder before running simulation ${white}\n"
		echo -e "${magenta}- Qmatey will quit in 10 seconds ${white}\n"
		echo -e "${magenta}- samples directory/folder contains fasta/fastq files. Delete samples directory/folder before running simulation ${white}\n" > ${projdir}/job_killed
		sleep 10 && exit 1
	else
		echo -e "${magenta}- generating sequence reads for simulated metagenome profiling of synthetic/mock community ${white}\n"
		time simulate_reads &>> log.out
	fi
fi


#################################################################################################################
#Organize fastq files for metagenomic processing
#Combine paired-end read data for downstream analysis
echo -e "\e[97m########################################################\n \e[38;5;210mOrganizing sample fastq files \n\e[97m########################################################\n"
organize_fq_files () {
cd "${projdir}"
cd samples
if [[ "simulate_reads" == 1 ]] && [[ "$simulation_lib" == "complete_digest" ]] && [[ "$simulation_lib" == "partial_digest" ]]; then
	export library_type=qRRS
fi
if [[ "simulate_reads" == 1 ]] && [[ "$simulation_lib" == "shotgun" ]]; then
	export library_type=shotgun
fi



if test -f filename_reformatted.txt; then
	echo -e "${magenta}- \n- file names reformatting was previously performed  ${white}\n"
else
	if [[ -d "se" ]]; then
		:
	else
		mkdir -p se
	fi
	if [[ -d "pe" ]]; then
		:
	else
		mkdir -p pe
	fi
	cd "${projdir}"/samples/se
	if [[ -z "$(ls -A ../pe 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ../se 2> /dev/null)" ]]; then
			cd ../
			for i in $(ls *.f*); do (
				if [[ "$i" == *"R2.f"* ]]; then
					:
				else
					if [[ "$i" == *.R1* ]]; then
						mv "$i" ${i/.R1/}
					elif [[ "$i" == *_R1* ]]; then
						mv "$i" ${i/_R1/}
					fi
				fi ) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait
		fi
	fi


	cd "${projdir}"/samples/pe
	if [ "$(ls -A ../pe 2> /dev/null)" ]; then
		if [ -z "$(ls -A ../se 2> /dev/null)" ]; then
			echo -e "${magenta}- only paired-end reads available in pe-folder ${white}\n"
			for i in $(ls *.f*); do (
				if [[ ! "$i" =~ R2.f ]]; then
					if [[ "$i" == *".R1"* ]]; then
						cat ${i%.R1*}.* > ../${i/.R1/} && rm "${i%.R1*}".* 2> /dev/null
					elif [[ "$i" == *_R1* ]]; then
						cat ${i%_R1*}_* > ../${i/_R1/} && rm "${i%_R1*}"_* 2> /dev/null
					else
						cat ${i%.f*}.* > ../$i && rm "${i%.f*}".* 2> /dev/null
					fi
				fi ) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait
		fi
	fi


	cd "${projdir}"/samples/se
	if [ "$(ls -A ../se 2> /dev/null)" ]; then
		if [ -z "$(ls -A ../pe 2> /dev/null)" ]; then
			echo -e "${magenta}- only single-end or unpaired reads available in se-folder ${white}\n"
			for i in $(ls *.f*); do (
				if [[ ! "$i" =~ R2.f ]]; then
					if [[ "$i" == *".R1"* ]]; then
						cat ${i%.R1*}.* > ../${i/.R1/} && rm "${i%.R1*}".* 2> /dev/null
					elif [[ "$i" == *_R1* ]]; then
						cat ${i%_R1*}_* > ../${i/_R1/} && rm "${i%_R1*}"_* 2> /dev/null
					else
						cat ${i%.f*}.* > ../$i && rm "${i%.f*}".* 2> /dev/null
					fi
				fi ) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait
		fi
	fi


	cd "${projdir}"/samples/pe
	if [ "$(ls -A ../se 2> /dev/null)" ]; then
		if [ "$(ls -A ../pe 2> /dev/null)" ]; then
			for i in $(ls *R1.f*); do (
				if [[ ! "$i" =~ R2.f ]]; then
					if [[ "$i" == *.R1* ]]; then
						cat ${i%.R1*}.* ../se/${i%.R1*}.* 2> /dev/null > ../${i/.R1/} && rm "${i%.R1*}".* ../se/${i%.R1*}.* 2> /dev/null
					elif [[ "$i" == *_R1* ]]; then
						cat ${i%_R1*}_* ../se/${i%_R1*}_* 2> /dev/null > ../${i/_R1/} && rm "${i%_R1*}"_* ../se/${i%_R1*}_* 2> /dev/null
					else
						cat ${i%.f*}.* ../se/${i%.f*}.* 2> /dev/null > ../${i} && rm "${i%.f*}".* ../se/${i%.f*}.* 2> /dev/null
					fi
				fi ) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait
		fi
	fi

	cd "${projdir}"/samples/
	find . -type d -empty -delete
	sampno=$(ls -1 | wc -l)
	if [[ "$sampno" == "0" ]]; then
		echo -e "${magenta}- \n- samples folder is empty, exiting pipeline ${white}\n"
		echo -e "${magenta}- \n- samples folder is empty, exiting pipeline ${white}\n" > ${projdir}/job_killed
		exit 0
	fi
	echo filename_formatted > filename_formatted.txt
fi

if test -f flushed_reads.txt; then
	echo -e "${magenta}- \n- improved flushed ends of reads was previously performed  ${white}\n"
else
	if [[ "$library_type" =~ "RRS" ]] || [[ "$library_type" =~ "rrs" ]] || [[ "$library_type" =~ "amplicon" ]] || [[ "$library_type" =~ "Amplicon" ]] || [[ "$library_type" =~ "AMPLICON" ]] || [[ "$library_type" =~ "16S" ]] || [[ "$library_type" =~ "16s" ]]|| [[ "$library_type" =~ "ITS" ]] || [[ "$library_type" =~ "its" ]]; then
		for i in $(ls *.f*); do (
			if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
				fa_fq=$(zcat ${projdir}/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
			else
				fa_fq=$(cat ${projdir}/samples/$i | head -n1 | cut -c1-1)
			fi
			if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
				if [[ "${fa_fq}" == "@" ]]; then
					awk 'NR%2==0' <(zcat "$i") | awk 'NR%2==1' | gzip > "${i}".tmp.gz && rm "$i" && mv "${i}".tmp.gz ${i%.f*}.fasta.gz
				fi
				if [[ "${fa_fq}" == ">" ]]; then
					awk 'NR%2==0' <(zcat "$i") | gzip > "${i}".tmp.gz && mv "${i}".tmp.gz "${i}"
				fi
			else
				if [[ "${fa_fq}" == "@" ]]; then
					awk 'NR%2==0' $i | awk 'NR%2==1' | gzip > "${i}".tmp.gz && rm "$i" && mv "${i}".tmp.gz ${i%.f*}.fasta.gz
				fi
				if [[ "${fa_fq}" == ">" ]]; then
					awk 'NR%2==0' $i | gzip > "${i}".tmp.gz && rm "$i" && mv "${i}".tmp.gz "${i}".gz
				fi
			fi ) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait
		for i in $(ls *_R2.fasta.gz 2> /dev/null); do
			if test -f $i; then
				cat "$i" ${i%_R2.fasta.gz}.fasta.gz > ${i%_R2.fasta.gz}.tmp.fasta &&
				rm "$i" && mv ${i%_R2.fasta.gz}.tmp.fasta ${i%_R2.fasta.gz}.fasta.gz
				wait
			fi
		done
		for i in $(ls *.R2.fasta.gz 2> /dev/null); do
			if test -f $i; then
				cat "$i" ${i%.R2.fasta.gz}.fasta.gz > ${i%.R2.fasta.gz}.tmp.fasta.gz &&
				rm "$i" && mv ${i%.R2.fasta.gz}.tmp.fasta.gz ${i%.R2.fasta.gz}.fasta.gz
				wait
			fi
		done
		wait

		if [[ "$library_type" =~ "amplicon" ]] || [[ "$library_type" =~ "Amplicon" ]] || [[ "$library_type" =~ "AMPLICON" ]] || [[ "$library_type" =~ "16S" ]] || [[ "$library_type" =~ "16s" ]]|| [[ "$library_type" =~ "ITS" ]] || [[ "$library_type" =~ "its" ]]; then
			for i in $(ls *.f*); do (
				frag=${i%.f*}
				if [[ "$i" == *"_compressed.f"* ]]; then
					:
				else
					zcat "$i" | grep -v '>' | awk '{print substr($0,1,150)}' | awk 'length >=32' | awk -v frag="$frag" '{print ">"frag"_"NR"\n"$0}' | \
					gzip > ${i%.f*}_tmp.fasta.gz && mv ${i%.f*}_tmp.fasta.gz $i
					wait
				fi
				) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait
		fi

		if [[ "$library_type" =~ "RRS" ]] || [[ "$library_type" =~ "rrs" ]]; then
			for i in $(ls *.f*); do (
				if [[ "$i" == *"_compressed.f"* ]]; then
					:
				else
					zcat "$i" | grep -v '>' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i%.f}_length_distribution.txt
				fi
				) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait

			for lenfile in $(find . -name "*_length_distribution.txt"); do
				cat "$lenfile" >> length_distribution.txt &&
				rm "$lenfile" &&
				wait
			done
			wait
			awk '{print length($0)}' length_distribution.txt | sort -T "${projdir}"/tmp -n > tmp.txt
			mv tmp.txt length_distribution.txt
			export max_seqread_len=$(awk '{all[NR] = $0} END{print all[int(NR*0.75 - 0.5)]}' length_distribution.txt)
			rm length_distribution.txt
			wait


			for i in $(ls *.f*); do (
				frag=${i%.f*}
				if [[ "$i" == *"_compressed.f"* ]]; then
					:
				else
					zcat "$i" | grep -v '>' | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk 'length >=32' | awk -v frag="$frag" '{print ">"frag"_"NR"\n"$0}' | \
					gzip > ${i%.f*}_tmp.fasta.gz && mv ${i%.f*}_tmp.fasta.gz $i
					wait
				fi
				) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait
		fi
	fi

	if [[ "$library_type" == "WGS" ]] || [[ "$library_type" == "wgs" ]] || [[ "$library_type" == "SHOTGUN" ]] || [[ "$library_type" == "shotgun" ]]; then

		for i in $(ls *.f*); do (
			if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
				fa_fq=$(zcat ${projdir}/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
			else
				fa_fq=$(cat ${projdir}/samples/$i | head -n1 | cut -c1-1)
			fi
			if [[ $(file $i 2> /dev/null) =~ gzip ]]; then
				if [[ "${fa_fq}" == "@" ]]; then
					awk 'NR%2==0' <(zcat "$i") | awk 'NR%2==1' | gzip > "${i}".tmp.gz && rm "$i" && mv "${i}".tmp.gz ${i%.f*}.fasta.gz
				fi
				if [[ "${fa_fq}" == ">" ]]; then
					awk 'NR%2==0' <(zcat "$i") | gzip > "${i}".tmp.gz && mv "${i}".tmp.gz "${i}"
				fi
			else
				if [[ "${fa_fq}" == "@" ]]; then
					awk 'NR%2==0' $i | awk 'NR%2==1' | gzip > "${i}".tmp.gz && rm "$i" && mv "${i}".tmp.gz ${i%.f*}.fasta.gz
				fi
				if [[ "${fa_fq}" == ">" ]]; then
					awk 'NR%2==0' $i | gzip > "${i}".tmp.gz && mv "${i}".tmp.gz "${i}".gz
				fi
			fi ) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait

		# concatenate R1 and R2 reads
		for i in $(ls *_R2.fasta.gz 2> /dev/null); do
			if test -f $i; then
				cat "$i" ${i%_R2.fasta.gz}.fasta.gz > ${i%_R2.fasta.gz}.tmp.fasta.gz &&
				rm "$i" && mv ${i%_R2.fasta.gz}.tmp.fasta.gz ${i%_R2.fasta.gz}.fasta.gz
				wait
			fi
		done
		wait
		for i in $(ls *.R2.fasta.gz 2> /dev/null); do
			if test -f $i; then
				cat "$i" ${i%.R2.fasta.gz}.fasta.gz > ${i%.R2.fasta.gz}.tmp.fasta.gz &&
				rm "$i" && mv ${i%.R2.fasta.gz}.tmp.fasta.gz ${i%.R2.fasta.gz}.fasta.gz
				wait
			fi
		done
		wait

		# if subsampling of WGS reads is not false, subsample with in silico qRRS
		if [[ "$subsample_shotgun_R1" != false ]]; then
			echo $subsample_shotgun_R1 | awk '{gsub(/,/,"\n");}1' | awk '{print "RE"NR"\t"$1}' > REnase_R1.txt
			RE1a=$(grep 'RE1' REnase_R1.txt | awk '{print $2}')
			RE1b=$(grep 'RE2' REnase_R1.txt | awk '{print $2}')
			RE1c=$(grep 'RE3' REnase_R1.txt | awk '{print $2}')
			RE1d=$(grep 'RE4' REnase_R1.txt | awk '{print $2}')
			if [[ -z "$RE1b" ]]; then RE1b=999; fi
			if [[ -z "$RE1c" ]]; then RE1c=999; fi
			if [[ -z "$RE1d" ]]; then RE1d=999; fi
			rm REnase_R1.txt
			if [[ "$subsample_shotgun_R2"  != false ]]; then
				echo $subsample_shotgun_R2 | awk '{gsub(/,/,"\n");}1' | awk '{print "RE"NR"\t"$1}' > REnase_R2.txt
				RE2a=$(grep 'RE1' REnase_R2.txt | awk '{print $2}')
				RE2b=$(grep 'RE2' REnase_R2.txt | awk '{print $2}')
				RE2c=$(grep 'RE3' REnase_R2.txt | awk '{print $2}')
				RE2d=$(grep 'RE4' REnase_R2.txt | awk '{print $2}')
				if [[ -z "$RE2b" ]]; then RE2b=999; fi
				if [[ -z "$RE2c" ]]; then RE2c=999; fi
				if [[ -z "$RE2d" ]]; then RE2d=999; fi
				rm REnase_R2.txt
			else
				RE2a=999
				RE2b=999
				RE2c=999
				RE2d=999
			fi
		fi

	  for i in $(ls *.f*); do (
			frag=${i%.f*}
			if [[ "$i" == *"_compressed.f"* ]]; then
				:
			else
				if [[ -z "$shotgun_min_fragment_length" ]]; then
					shotgun_min_fragment_length=100
				fi
				if [[ "$subsample_shotgun_R1" != false ]]; then
					awk -v RE1a="$RE1a" -v RE1b="$RE1b" -v RE1c="$RE1c" -v RE1d="$RE1d" -v RE2a="$RE2a" -v RE2b="$RE2b" -v RE2c="$RE2c" -v RE2d="$RE2d" \
					'{gsub(RE1a,RE1a"\n"RE1a); gsub(RE1b,RE1b"\n"RE1b); gsub(RE1c,RE1c"\n"RE1c); gsub(RE1d,RE1d"\n"RE1d); \
					gsub(RE2a,RE2a"\n"RE2a); gsub(RE2b,RE2b"\n"RE2b); gsub(RE2c,RE2c"\n"RE2c); gsub(RE2d,RE2d"\n"RE2d); }1' <(zcat "$i") | \
					awk -v min="$shotgun_min_fragment_length" 'length >= min' > ${i%.f*}_chopped.txt

					if [[ "$HDsubsample" == true ]]; then
						grep "^$RE1a.*$RE1a$\|^$RE1b.*$RE1b$\|^$RE1c.*$RE1c$\|^$RE1d.*$RE1d$\|^$RE2a.*$RE2a$\|^$RE2b.*$RE2b$\|^$RE2c.*$RE2c$\|^$RE2d.*$RE2d$" ${i%.f*}_chopped.txt > ${i%.f*}.tmp1.txt &&
						sleep 2
					else
						touch ${i%.f*}.tmp1.txt
					fi
					grep "^$RE1a.*$RE2a$\|^$RE1a.*$RE2b$\|^$RE1a.*$RE2c$\|^$RE1a.*$RE2d$\|^$RE1b.*$RE2a$\|^$RE1b.*$RE2b$\|^$RE1b.*$RE2c$\|^$RE1b.*$RE2d$" ${i%.f*}_chopped.txt 2> /dev/null > ${i%.f*}.tmp2.txt &&
					sleep 2
					grep "^$RE1c.*$RE2a$\|^$RE1c.*$RE2b$\|^$RE1c.*$RE2c$\|^$RE1c.*$RE2d$\|^$RE1d.*$RE2a$\|^$RE1d.*$RE2b$\|^$RE1d.*$RE2c$\|^$RE1d.*$RE2d$" ${i%.f*}_chopped.txt 2> /dev/null >> ${i%.f*}.tmp2.txt &&
					sleep 2
					grep "^$RE2a.*$RE1a$\|^$RE2a.*$RE1b$\|^$RE2a.*$RE1c$\|^$RE2a.*$RE1d$\|^$RE2b.*$RE1a$\|^$RE2b.*$RE1b$\|^$RE2b.*$RE1c$\|^$RE2b.*$RE1d$" ${i%.f*}_chopped.txt 2> /dev/null >> ${i%.f*}.tmp2.txt &&
					sleep 2
					grep "^$RE2c.*$RE1a$\|^$RE2c.*$RE1b$\|^$RE2c.*$RE1c$\|^$RE2c.*$RE1d$\|^$RE2d.*$RE1a$\|^$RE2d.*$RE1b$\|^$RE2d.*$RE1c$\|^$RE2d.*$RE1d$" ${i%.f*}_chopped.txt 2> /dev/null >> ${i%.f*}.tmp2.txt &&
					sleep 2
					wait
					rm "${i%.f*}"_chopped.txt

					cat ${i%.f*}.tmp1.txt ${i%.f*}.tmp2.txt | awk 'length >=32' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.fasta.gz &&
					rm "${i%.f*}".tmp1.txt ${i%.f*}.tmp2.txt
					wait
				fi
				if [[ "$subsample_shotgun_R1" == false ]]; then
					zcat "$i" | grep -v '>' | awk 'length >=32' | awk '{print ">frag"NR"\n"$0}' | $gzip > ${i%.f*}.tmp.fasta.gz &&
					mv ${i%.f*}.tmp.fasta.gz $i
				fi
			fi ) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait

		# perform reads flushing
		if [[ "$subsample_shotgun_R1" != false ]]; then
			for i in $(ls *.f*); do (
				if [[ "$i" == *"_compressed.f"* ]]; then
					:
				else
					zcat "$i" | grep -v '>' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=1000; i++){x=int(rand()*NR) + 1; print a[x];}}' > ${i%.f}_length_distribution.txt
				fi
				) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait

			for lenfile in $(find . -name "*_length_distribution.txt"); do
				cat $lenfile >> length_distribution.txt &&
				rm "$lenfile" &&
				wait
			done
			wait
			awk '{print length($0)}' length_distribution.txt | sort -T "${projdir}"/tmp -n > tmp.txt
			mv tmp.txt length_distribution.txt
			export max_seqread_len=$(awk '{all[NR] = $0} END{print all[int(NR*0.75 - 0.5)]}' length_distribution.txt)
			rm length_distribution.txt
			wait


			for i in $(ls *.f*); do (
				frag=${i%.f*}
				if [[ "$i" == *"_compressed.f"* ]]; then
					:
				else
					zcat "$i" | grep -v '>' | awk -v max=$max_seqread_len '{print substr($0,1,max)}' | awk -v frag="$frag" '{print ">"frag"_"NR"\n"$0}' | \
					gzip > ${i%.f*}_tmp.fasta.gz && mv ${i%.f*}_tmp.fasta.gz $i
					wait
				fi
				) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			wait
		fi
	fi
	rm *fastq.gz 2> /dev/null
	find . -type d -empty -delete
	echo flushed_reads > flushed_reads.txt
fi

}
cd "${projdir}"
cd samples
if [[ -d "pe" ]]; then
	fqpass=$(find ./pe -maxdepth 1 -name '*_R1.f*' -o -name '*_R2.f*' -o -name '*.R1.f*' -o -name '*.R2.f*' | wc -l)
	fqfail=$(ls ./pe/* | wc -l)
	fqfail=$((fqfail-fqpass))
	if [[ "$fqfail" -lt 1 ]]; then
		if [[ "$fqpass" -gt 0 ]]; then
			echo -e "${YELLOW}- Qmatey is organizing sample fastq files  ${WHITE}"
			time organize_fq_files &>> ${projdir}/log.out
		fi
	else
		echo -e "${magenta}- samples' PE fastq filenames requires formatting (i.e. needs to end in "_R1.fastq" or ".R1.fastq" and "_R2.fastq" or ".R2.fastq") ${white}\n"
		echo -e "${magenta}- samples' PE fastq filenames requires formatting (i.e. needs to end in "_R1.fastq" or ".R1.fastq" and "_R2.fastq" or ".R2.fastq") ${white}\n" > ${projdir}/job_killed
		sleep 5 && exit 0
	fi
else
	echo -e "${YELLOW}- Qmatey is organizing sample fastq files  ${WHITE}"
	time organize_fq_files &>> ${projdir}/log.out
fi
total_no_samples=$(ls "${projdir}"/samples/* | wc -l)


cd "${projdir}"
if [[ -z "$(ls -A ./norm_ref 2> /dev/null)" ]]; then
	echo -e "$1 \e[31m reference/normalization genome folder is empty   ${WHITE}"
else
	cd norm_ref
	#Index ref genomes using BWA, samtool, and java dependencies -- concatenate all reference genomes into a master reference file
	file=*.dict
	if test -f $file; then
		echo -e "${YELLOW}- indexed genome already available ${WHITE}\n"
	else
		echo -e "${YELLOW}- indexing normalization reference genome ${WHITE}"
    for i in $(ls *.fa*); do
      if (file "$i" | grep -q compressed ); then
        gunzip $i
      fi
    done
    wait
    for i in $(ls *.fna*); do
      if (file "$i" | grep -q compressed ); then
        gunzip $i
      fi
    done
    wait
		for i in $(ls *.fa*); do
      if test -f "$i"; then
  			n=">${i%.fa*}_"
  			awk '{ sub("\r$",""); print}' $i | awk -v n="$n" '{gsub(/>/,n); print}' >> master_ref.fasta
      fi
		done
		for i in $(ls *.fna); do
      if test -f "$i"; then
  			n=">${i%.fa*}_"
  			awk '{ sub("\r$",""); print}' $i | awk -v n="$n" '{gsub(/>/,n); print}' >> master_ref.fasta
      fi
		done
		wait

		$bwa index -a bwtsw master_ref.fasta
		$samtools faidx master_ref.fasta
		$java -jar $picard CreateSequenceDictionary REFERENCE=master_ref.fasta    OUTPUT=master_ref.dict
	fi
fi



#################################################################################################################
echo -e "\e[97m########################################################\n \e[38;5;210m Read Compression, Normalization, and exclusion of reads from reference genomes \n\e[97m########################################################\n"
ref_norm () {
	cd "${projdir}"
	# print "sample\treads_input\tsubsetted_filtered_reads\tPercent_retained" ${projdir}/metagenome/proportion_retained_reads.txt
	# export startmotif=TTT,ATT,CTT,GTT,TAT,AAT,CAT,GAT,TCT,ACT,CCT,GCT,TGT,AGT,CGT,GGT,TTA,ATA,CTA,GTA,TAA,AAA,CAA,GAA,TCA,ACA,CCA,GCA,TGA,AGA,CGA,GGA,TTC,ATC,CTC,GTC,TAC,AAC,CAC,GAC,TCC,ACC,CCC,GCC,TGC,AGC,CGC,GGC,TTG,ATG,CTG,GTG,TAG,AAG,CAG,GAG,TCG,ACG,CCG,GCG,TGG,AGG,CGG,GGG
	if [[ -z "$(ls -A ./norm_ref/*.dict 2> /dev/null)" ]]; then
		echo -e "$1 \e[31m normalization reference folder is empty, Qmatey will not exclude any read"
		echo -e "$1 \e[31m Qmatey will use read coverage of samples for normalization"
		cd "${projdir}"/samples
		#All duplicate reads are compressed into one representative read with duplication reflected as a numeric value
		#Increases the speed of reference genome alignment -- especially if read depth is high
		rm "${projdir}"/metagenome/microbiome_coverage.txt 2> /dev/null


		for i in $(ls *.f*); do (
			if [[ "$i" == *"_compressed.f"* ]]; then
				:
			else
				if [[ "$subsample_shotgun_R1" == false ]]; then
					zcat "$i" 2> /dev/null | awk 'NR%2==0' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -v min="$minRD" '$1>=min{print $1"\t"$2}' | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | gzip > ${i%.f*}_compressed.fasta.gz
				else
					zcat "$i" 2> /dev/null | awk 'NR%2==0' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -v min="$minRD" '$1>=min{print $1"\t"$2}' | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | gzip > ${i%.f*}_compressed.fasta.gz
					wait
				fi
				wait
				if [[ "$simulate_reads" == 1 ]]; then
					:
				else
					if [[ "$library_type" =~ "amplicon" ]] || [[ "$library_type" =~ "Amplicon" ]] || [[ "$library_type" =~ "AMPLICON" ]] || [[ "$library_type" =~ "16S" ]] || [[ "$library_type" =~ "16s" ]]|| [[ "$library_type" =~ "ITS" ]] || [[ "$library_type" =~ "its" ]]; then
						minimumRD=4
						printf "minimum read depth threshold(${i%.f*}): \t$minimumRD\n" >> ${projdir}/metagenome/minimumRD_empirical.txt
						zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t");}1' | awk -v pat="$minimumRD" '$2 >= pat' | awk '{print $1"-"$2"\t"$3}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
						mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
					else
						minimumRD=$(zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t"); print $2}' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=10000; i++){x=int(rand()*NR) + 1; print a[x];}}' | \
						sort -T "${projdir}"/tmp -Vr | awk '{all[NR] = $0} END{print all[int(NR*0.25 - 0.5)]}')
						if [[ "$library_type" == "WGS" ]] || [[ "$library_type" == "wgs" ]] || [[ "$library_type" == "SHOTGUN" ]] || [[ "$library_type" == "shotgun" ]]; then
							if [[ "$subsample_shotgun_R1" == false ]]; then
								minimumRD=$((minimumRD * 1))
							fi
						fi
						printf "minimum read depth threshold(${i%.f*}): \t$minimumRD\n" >> ${projdir}/metagenome/minimumRD_empirical.txt
						zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t");}1' | awk -v pat="$minimumRD" '$2 >= pat' | awk '{print $1"-"$2"\t"$3}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
						mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
					fi
				fi
				# Nreads_input=$(zcat $i | grep '>' | wc -l)
				# Nreads_filter=$(zcat ${i%.f*}_compressed.fasta.gz | grep '>' | awk -F '-' '{s+=$2}END{print s}')
				# Percentage=$($Nreads_filter * 100)
				# awk -v samp="${i%.f*}" -v pat1="$Nreads_input" pat2="$Nreads_filter" -v pat3="$Percentage" '{print samp"\t"pat1"\t"pat2"\t"pat3/pat1}' >> ${projdir}/metagenome/proportion_retained_reads.txt
			fi ) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait
		#sample read depth is used to normalize quantification data
		echo -e "${YELLOW}- calculating a normalization factor"
		for i in $(ls *_compressed.fasta.gz); do
			zcat "$i" 2> /dev/null | grep '^>' | awk -v sample=${i%_compressed.fasta.gz} -F '-' '{s+=$2}END{print sample"\t"s}' > ${projdir}/metagenome/${i%_compressed.fasta.gz}_microbiome_coverage.txt && \
			cat ${projdir}/metagenome/${i%_compressed.fasta.gz}_microbiome_coverage.txt >> ${projdir}/metagenome/microbiome_coverage.txt
			rm "${projdir}"/metagenome/${i%_compressed.fasta.gz}_microbiome_coverage.txt
		done
		wait
		maximum=$(sort -T "${projdir}"/tmp -Vr -k2,2 ${projdir}/metagenome/microbiome_coverage.txt | awk 'NF > 0' | awk 'NR==1{print $2; exit}')
		awk -v maximum=$maximum '{print $1,maximum/$2}' ${projdir}/metagenome/microbiome_coverage.txt | cat <(printf 'Sample_ID\tNormalization_factor\n') - > ${projdir}/metagenome/coverage_normalization_factor.txt
		rm "${projdir}"/metagenome/microbiome_coverage.txt

	else
		cd "${projdir}"/samples
		#All duplicate reads are compressed into one representative read with duplication reflected as a numeric value
		#Increases the speed of reference genome alignment -- especially if read depth is high
		rm "${projdir}"/metagenome/microbiome_coverage.txt 2> /dev/null

		for i in $(ls *.f*); do (
			if [[ "$i" == *"_compressed.f"* ]]; then
				:
			else
				if [[ "$subsample_shotgun_R1" == false ]]; then
					zcat "$i" 2> /dev/null | awk 'NR%2==0' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -v min="$minRD" '$1>=min{print $1"\t"$2}' | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | gzip > ${i%.f*}_compressed.fasta.gz
				else
					zcat "$i" 2> /dev/null | awk 'NR%2==0' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -v min="$minRD" '$1>=min{print $1"\t"$2}' | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | gzip > ${i%.f*}_compressed.fasta.gz
					wait
				fi
				wait
				if [[ "$simulate_reads" == 1 ]]; then
					zcat ${i%.f*}_compressed.fasta.gz | awk '{print $1"\n"$2}' | $gzip > ${i%.f*}_compressed.tmpfasta.gz
					mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
				else
					if [[ "$library_type" =~ "amplicon" ]] || [[ "$library_type" =~ "Amplicon" ]] || [[ "$library_type" =~ "AMPLICON" ]] || [[ "$library_type" =~ "16S" ]] || [[ "$library_type" =~ "16s" ]]|| [[ "$library_type" =~ "ITS" ]] || [[ "$library_type" =~ "its" ]]; then
						minimumRD=4
						printf "minimum read depth threshold(${i%.f*}): \t$minimumRD\n" >> ${projdir}/metagenome/minimumRD_empirical.txt
						zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t");}1' | awk -v pat="$minimumRD" '$2 >= pat' | awk '{print $1"-"$2"\n"$3}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
						mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
					else
						minimumRD=$(zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t"); print $2}' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=10000; i++){x=int(rand()*NR) + 1; print a[x];}}' | \
						sort -T "${projdir}"/tmp -Vr | awk '{all[NR] = $0} END{print all[int(NR*0.25 - 0.5)]}')
						if [[ "$library_type" == "WGS" ]] || [[ "$library_type" == "wgs" ]] || [[ "$library_type" == "SHOTGUN" ]] || [[ "$library_type" == "shotgun" ]]; then
							if [[ "$subsample_shotgun_R1" == false ]]; then
								minimumRD=$((minimumRD * 1))
							fi
						fi
						printf "minimum read depth threshold(${i%.f*}): \t$minimumRD\n" >> ${projdir}/metagenome/minimumRD_empirical.txt
						zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t");}1' | awk -v pat="$minimumRD" '$2 >= pat' | awk '{print $1"-"$2"\n"$3}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
						mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
					fi
				fi
			fi ) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait

		cd "${projdir}"/samples
		if [[ "${norm_method}" == spike_host ]]; then
			#Aligning compressed sample files (read-depth accounted for) to the master reference genomes
			#Exludes all host/reference genome data from the samples -- leaving only metagenomic reads
			echo -e "${YELLOW}- aligning sample reads to normalization reference genome${WHITE}"
			for i in $(ls *_compressed.fasta.gz); do
				if test ! -f ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt; then
					cd "${projdir}"/norm_ref
					if [[ -z "$(ls "${projdir}"/metagenome/"${i%_compressed.fasta.gz}".sam 2> /dev/null)" ]]; then
						$bwa mem -t "$threads" master_ref.fasta ${projdir}/samples/$i > ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam
						wait
					fi
					cd "${projdir}"/samples
					$samtools view -S -b ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam > ${projdir}/metagenome/${i%_compressed.fasta.gz}.bam
					#$java -XX:ParallelGCThreads=$gthreads -jar $picard SortSam I= ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam O= ${projdir}/metagenome/${i%_compressed.fasta.gz}.bam SORT_ORDER=coordinate && \
					printf '\n###---'${i%.f*}'---###\n' > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt && \
					$samtools flagstat ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt &&
					rm "${projdir}"/metagenome/${i%_compressed.fasta.gz}.sam
				fi
			done
			wait

			cat ${projdir}/metagenome/results/ref_aligned_summaries/*_summ.txt > ${projdir}/metagenome/results/ref_aligned_summaries/ref_aligned_summaries_unique_reads.txt &&
			rm -r "${projdir}"/metagenome/results/ref_aligned_summaries/*_summ.txt

			cd "${projdir}"/metagenome
			rm host_coverage.txt microbiome_coverage.txt 2> /dev/null
			#Host-reference alignment coverage relative to other samples is used to normalize quantification data
			echo -e "${YELLOW}- calculating a normalization factor"
			for i in $(ls *.bam); do
				$samtools view -F 4 $i | grep -vwE "(@HD|@SQ|@PG)" | awk '{print $1}' | awk -v sample=${i%.bam} -F '-' '{s+=$2}END{print sample"\t"s}' > ${i%.bam}_host_coverage.txt && \
				cat ${i%.bam}_host_coverage.txt >> host_coverage.txt && \
				rm "${i%.bam}"_host_coverage.txt  && \
				$samtools view -f 4 $i | awk '{print $1}' | awk -v sample=${i%.bam} -F '-' '{s+=$2}END{print sample"\t"s}' > ${i%.bam}_microbiome_coverage.txt && \
				cat ${i%.bam}_microbiome_coverage.txt >> microbiome_coverage.txt && \
				rm "${i%.bam}"_microbiome_coverage.txt
			done
			wait
			awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' host_coverage.txt microbiome_coverage.txt | awk '{print $1,"\t",$2,"\t",$4,"\t",$2+$4,"\t",($2/($2+$4))*100}' | cat <(printf 'Sample_ID\t#_metagenome_reads\t#_Host_reads\t#_total_reads\tpercent_metagenome\n') - > coverage_normalize.txt &&
			maximum=$(sort -T "${projdir}"/tmp -Vr -k3,3 coverage_normalize.txt | awk '!($3 == 0)' | awk 'NF > 0' | awk 'NR==1{print $3; exit}')
			awk -v maximum=$maximum 'NR>1 && $3>0{print $1"\t"maximum/$3}' coverage_normalize.txt | cat <(printf 'Sample_ID\tNormalization_factor\n') - | \
			cat - <(awk 'NR>1 && $3==0{print $1"\t"$3}' coverage_normalize.txt) > coverage_normalization_factor.txt &&
			awk 'NR>1{print $1,"\t",($2/$4)*100,"\t",($3/$4)*100}' coverage_normalize.txt | cat <(printf 'Sample_ID\tPercent_Metagenome\tPercent_Host\n') - > ./results/metagenome_derived_perc.txt &&
			rm host_coverage.txt microbiome_coverage.txt &&
			wait
		fi

		if [[ "${norm_method}" == samples ]]; then
			#Exludes all host/reference genome data from the samples -- leaving only metagenomic reads
			echo -e "${YELLOW}- aligning sample reads to normalization reference genome${WHITE}"
			for i in $(ls *_compressed.fasta.gz); do
				if test ! -f ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt; then
					cd "${projdir}"/norm_ref
					if [[ -z "$(ls "${projdir}"/metagenome/"${i%_compressed.fasta.gz}".sam 2> /dev/null)" ]]; then
						$bwa mem -t "$threads" master_ref.fasta ${projdir}/samples/$i > ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam
						wait
					fi
					cd "${projdir}"/samples
					$samtools view -S -b ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam > ${projdir}/metagenome/${i%_compressed.fasta.gz}.bam &&
					# $java -XX:ParallelGCThreads=$gthreads -jar $picard SortSam I= ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam O= ${projdir}/metagenome/${i%_compressed.fasta.gz}.bam SORT_ORDER=coordinate && \
					printf '\n###---'${i%.f*}'---###\n' > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt && \
					$samtools flagstat ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt && \
					rm "${projdir}"/metagenome/${i%_compressed.fasta.gz}.sam &&
					wait
				fi
			done
			wait

			cat ${projdir}/metagenome/results/ref_aligned_summaries/*_summ.txt > ${projdir}/metagenome/results/ref_aligned_summaries/ref_aligned_summaries_unique_reads.txt &&
			rm -r "${projdir}"/metagenome/results/ref_aligned_summaries/*_summ.txt

			cd "${projdir}"/metagenome/
			for i in $(ls *.bam); do (
				$samtools view -f 4 $i | awk '{print $1}' | \
				awk -F '\t' 'NR==FNR{a[$0];next} $1 in a' - <(zcat ../samples/${i%.bam}_compressed.fasta.gz | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | awk '{gsub(/>/,"");}1') > ../samples/${i%.bam}_compressed.fasta &&
				if test -f ../samples/${i%.bam}_compressed.fasta; then rm "$i" ../samples/${i%.bam}_compressed.fasta.gz; fi
				wait
				gzip ../samples/${i%.bam}_compressed.fasta &&
				wait ) &
	      if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
	      wait
	      fi
			done
			wait

			cd "${projdir}"/samples/
			#sample read depth is used to normalize quantification data
			echo -e "${YELLOW}- calculating a normalization factor"
			for i in $(ls *_compressed.fasta.gz); do
				zcat "$i" 2> /dev/null | awk -v sample=${i%_compressed.fasta.gz} -F '-' '{s+=$2}END{print sample"\t"s}' > ${projdir}/metagenome/${i%_compressed.fasta.gz}_microbiome_coverage.txt && \
				cat "${projdir}"/metagenome/${i%_compressed.fasta.gz}_microbiome_coverage.txt >> "${projdir}"/metagenome/microbiome_coverage.txt &&
				rm "${projdir}"/metagenome/${i%_compressed.fasta.gz}_microbiome_coverage.txt &&
				wait
			done
			wait
			maximum=$(sort -T "${projdir}"/tmp -Vr -k2,2 ${projdir}/metagenome/microbiome_coverage.txt | awk 'NF > 0' | awk 'NR==1{print $2; exit}')
			awk -v maximum=$maximum '{print $1,maximum/$2}' ${projdir}/metagenome/microbiome_coverage.txt | cat <(printf 'Sample_ID\tNormalization_factor\n') - > ${projdir}/metagenome/coverage_normalization_factor.txt &&
			rm "${projdir}"/metagenome/microbiome_coverage.txt &&
			wait
		fi

	fi

	if [[ -z "$(ls -A ${projdir}/norm_ref/*.dict 2> /dev/null)" ]]; then
		cd "${projdir}"/samples
		echo -e "${YELLOW}- compile metagenome reads & compute relative read depth ${WHITE}"
		for i in $(ls *_compressed.fasta.gz); do
			normfactor=$( awk -v sample=${i%_compressed.fasta.gz} '$1 == sample' ${projdir}/metagenome/coverage_normalization_factor.txt | awk '{print $2}' ) && \
			awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t"); print}' <(zcat "$i" 2> /dev/null) | awk -v norm=$normfactor '{print $1"-"$2*norm"\n"$3}' | $gzip > ${projdir}/metagenome/haplotig/${i%_compressed.fasta.gz}_metagenome.fasta.gz
		done
		wait
		rm "${projdir}"/samples/*_compressed.fasta.gz

	else
		if [[ "${norm_method}" == spike_host ]]; then
			cd "${projdir}"/metagenome
			echo -e "${YELLOW}- compile metagenome reads into fasta format & compute relative read depth ${WHITE}"
			for i in $(ls *.bam); do (
				normfactor=$( awk -v sample=${i%.bam} '$1 == sample' coverage_normalization_factor.txt | awk '{print $2}' ) && \
				$samtools view -f 4 $i | awk '{print $1}' | \
				awk -F '\t' 'NR==FNR{a[$0];next} $1 in a' - <(zcat ../samples/${i%.bam}_compressed.fasta.gz | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | awk '{gsub(/>/,"");}1') | \
				awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t"); gsub(/>/,""); print}' | awk -v norm=$normfactor '{print ">"$1"-"$2*norm"\n"$3}' | $gzip > ./haplotig/${i%.bam}_metagenome.fasta.gz &&
				wait
				) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
				fi
			done
			wait
			find ./haplotig/ -size 0 -delete
			if [[ "$(ls ./haplotig/${i%.bam}_metagenome.fasta.gz | wc -l)" -gt 0 ]]; then rm *.bam "${projdir}"/samples/*_compressed.fasta.gz; fi
			wait
		fi

		if [[ "${norm_method}" == samples ]]; then
			cd "${projdir}"/samples
			echo -e "${YELLOW}- compile metagenome reads & compute relative read depth ${WHITE}"
			for i in $(ls *_compressed.fasta.gz); do (
				normfactor=$( awk -v sample=${i%_compressed.fasta.gz} '$1 == sample' ${projdir}/metagenome/coverage_normalization_factor.txt | awk '{print $2}' ) &&
				awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t"); gsub(/>/,""); print}' <(zcat "$i" 2> /dev/null) | awk -v norm=$normfactor '{print ">"$1"-"$2*norm"\n"$3}' | $gzip > ${projdir}/metagenome/haplotig/${i%_compressed.fasta.gz}_metagenome.fasta.gz
				) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
				fi
			done
			wait
			rm "${projdir}"/samples/*_compressed.fasta.gz

		fi
	fi

}
cd "${projdir}"
if [[ "$normalization" == true ]]; then
	if ls ${projdir}/metagenome/coverage_normalization_factor.txt 1> /dev/null 2>&1; then
		echo -e "${YELLOW}- coverage normalization factor already computed ${WHITE}"
		echo -e "${YELLOW}- Qmatey is skipping normalization ${WHITE}"
	else
		echo -e "${YELLOW}- Qmatey is performing normalization and file compression ${WHITE}"
		time ref_norm &>> ${projdir}/log.out
	fi
else
	echo -e "${YELLOW}- Qmatey is skipping normalization ${WHITE}"
fi

#################################################################################################################
echo -e "\e[97m########################################################\n \e[38;5;210m Read compression and exclusion of reads from reference genomes \n\e[97m########################################################\n"
no_norm () {
	cd "${projdir}"
	export startmotif=TTT,ATT,CTT,GTT,TAT,AAT,CAT,GAT,TCT,ACT,CCT,GCT,TGT,AGT,CGT,GGT,TTA,ATA,CTA,GTA,TAA,AAA,CAA,GAA,TCA,ACA,CCA,GCA,TGA,AGA,CGA,GGA,TTC,ATC,CTC,GTC,TAC,AAC,CAC,GAC,TCC,ACC,CCC,GCC,TGC,AGC,CGC,GGC,TTG,ATG,CTG,GTG,TAG,AAG,CAG,GAG,TCG,ACG,CCG,GCG,TGG,AGG,CGG,GGG
	if [[ -z "$(ls -A ./norm_ref/*.dict 2> /dev/null)" ]]; then
		echo -e "$1 \e[31m normalization reference folder is empty, Qmatey will not exclude any read"
		cd "${projdir}"/samples

		for i in $(ls *.f*); do (
			if [[ "$i" == *"_compressed.f"* ]]; then
				:
			else
				if [[ "$subsample_shotgun_R1" == false ]]; then
					zcat "$i" 2> /dev/null | awk 'NR%2==0' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -v min="$minRD" '$1>=min{print $1"\t"$2}' | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | gzip > ${i%.f*}_compressed.fasta.gz
				else
					zcat "$i" 2> /dev/null | awk 'NR%2==0' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -v min="$minRD" '$1>=min{print $1"\t"$2}' | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | gzip > ${i%.f*}_compressed.fasta.gz
					wait
				fi
				wait
				if [[ "$simulate_reads" == 1 ]]; then
					zcat ${i%.f*}_compressed.fasta.gz | awk '{print $1"\n"$2}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
					mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
				else
					if [[ "$library_type" =~ "amplicon" ]] || [[ "$library_type" =~ "Amplicon" ]] || [[ "$library_type" =~ "AMPLICON" ]] || [[ "$library_type" =~ "16S" ]] || [[ "$library_type" =~ "16s" ]]|| [[ "$library_type" =~ "ITS" ]] || [[ "$library_type" =~ "its" ]]; then
						minimumRD=4
						printf "minimum read depth threshold(${i%.f*}): \t$minimumRD\n" >> ${projdir}/metagenome/minimumRD_empirical.txt
						zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t");}1' | awk -v pat="$minimumRD" '$2 >= pat' | awk '{print $1"-"$2"\n"$3}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
						mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
					else
						minimumRD=$(zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t"); print $2}' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=10000; i++){x=int(rand()*NR) + 1; print a[x];}}' | \
						sort -T "${projdir}"/tmp -Vr | awk '{all[NR] = $0} END{print all[int(NR*0.25 - 0.5)]}')
						if [[ "$library_type" == "WGS" ]] || [[ "$library_type" == "wgs" ]] || [[ "$library_type" == "SHOTGUN" ]] || [[ "$library_type" == "shotgun" ]]; then
							if [[ "$subsample_shotgun_R1" == false ]]; then
								minimumRD=$((minimumRD * 1))
							fi
						fi
						printf "minimum read depth threshold(${i%.f*}): \t$minimumRD\n" >> ${projdir}/metagenome/minimumRD_empirical.txt
						zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t");}1' | awk -v pat="$minimumRD" '$2 >= pat' | awk '{print $1"-"$2"\n"$3}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
						mv ${i%.f*}_compressed.tmp.fasta.gz ${projdir}/metagenome/haplotig/${i%.f*}_compressed.fasta.gz
					fi
				fi
			fi ) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait
	else
		cd "${projdir}"/samples
		#All duplicate reads are compressed into one representative read with duplication reflected as a numeric value
		#Increases the spead of reference genome alignment -- especially if read depth is high
		for i in $(ls *.f*); do (
			if [[ "$i" == *"_compressed.f"* ]]; then
				:
			else
				if [[ "$subsample_shotgun_R1" == false ]]; then
					zcat "$i" 2> /dev/null | awk 'NR%2==0' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -v min="$minRD" '$1>=min{print $1"\t"$2}' | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | gzip > ${i%.f*}_compressed.fasta.gz
				else
					zcat "$i" 2> /dev/null | awk 'NR%2==0' | awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -v min="$minRD" '$1>=min{print $1"\t"$2}' | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | gzip > ${i%.f*}_compressed.fasta.gz
					wait
				fi
				wait
				if [[ "$simulate_reads" == 1 ]]; then
					zcat ${i%.f*}_compressed.fasta.gz | awk '{print $1"\n"$2}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
					mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
				else
					if [[ "$library_type" =~ "amplicon" ]] || [[ "$library_type" =~ "Amplicon" ]] || [[ "$library_type" =~ "AMPLICON" ]] || [[ "$library_type" =~ "16S" ]] || [[ "$library_type" =~ "16s" ]]|| [[ "$library_type" =~ "ITS" ]] || [[ "$library_type" =~ "its" ]]; then
						minimumRD=4
						printf "minimum read depth threshold(${i%.f*}): \t$minimumRD\n" >> ${projdir}/metagenome/minimumRD_empirical.txt
						zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t");}1' | awk -v pat="$minimumRD" '$2 >= pat' | awk '{print $1"-"$2"\n"$3}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
						mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
					else
						minimumRD=$(zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t"); print $2}' | awk 'BEGIN{srand();} {a[NR]=$0} END{for(i=1; i<=10000; i++){x=int(rand()*NR) + 1; print a[x];}}' | \
						sort -T "${projdir}"/tmp -Vr | awk '{all[NR] = $0} END{print all[int(NR*0.25 - 0.5)]}')
						if [[ "$library_type" == "WGS" ]] || [[ "$library_type" == "wgs" ]] || [[ "$library_type" == "SHOTGUN" ]] || [[ "$library_type" == "shotgun" ]]; then
							if [[ "$subsample_shotgun_R1" == false ]]; then
								minimumRD=$((minimumRD * 1))
							fi
						fi
						printf "minimum read depth threshold(${i%.f*}): \t$minimumRD\n" >> ${projdir}/metagenome/minimumRD_empirical.txt
						zcat ${i%.f*}_compressed.fasta.gz | awk '{gsub(/-/,"\t");}1' | awk -v pat="$minimumRD" '$2 >= pat' | awk '{print $1"-"$2"\n"$3}' | $gzip > ${i%.f*}_compressed.tmp.fasta.gz
						mv ${i%.f*}_compressed.tmp.fasta.gz ${i%.f*}_compressed.fasta.gz
					fi
				fi
			fi ) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait

		cd "${projdir}"/samples
		#Aligning compressed sample files (read-depth accounted for) to the master reference genomes
		#Exludes all host/reference genome data from the samples -- leaving only metagenomic reads
		echo -e "${YELLOW}- aligning sample reads to normalization reference genome${WHITE}"
		for i in $(ls *_compressed.fasta.gz); do
			if test ! -f ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt; then
				cd "${projdir}"/norm_ref/
				if [[ -z "$(ls "${projdir}"/metagenome/"${i%_compressed.fasta.gz}".sam 2> /dev/null)" ]]; then
					$bwa mem -t "$threads" master_ref.fasta ${projdir}/samples/$i > ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam
					wait
				fi
				cd "${projdir}"/samples
				$samtools view -S -b ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam > ${projdir}/metagenome/${i%_compressed.fasta.gz}.bam &&
				# $java -XX:ParallelGCThreads=$gthreads -jar $picard SortSam I= ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam O= ${projdir}/metagenome/${i%_compressed.fasta.gz}.bam SORT_ORDER=coordinate && \
				printf '\n###---'${i%.f*}'---###\n' > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt && \
				$samtools flagstat ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt && \
				rm "${projdir}"/metagenome/${i%_compressed.fasta.gz}.sam
				wait
			fi
		done
		wait

		cat "${projdir}"/metagenome/results/ref_aligned_summaries/*_summ.txt > "${projdir}"/metagenome/results/ref_aligned_summaries/ref_aligned_summaries_unique_reads.txt &&
		rm -r "${projdir}"/metagenome/results/ref_aligned_summaries/*_summ.txt &&
		wait

		cd "${projdir}"/metagenome
		echo -e "${YELLOW}- compile metagenome reads into fasta format ${WHITE}"
		for i in $(ls *.bam); do (
			if test ! -f ./haplotig/${i%.bam}_metagenome.fasta.gz; then
				$samtools view -f 4 $i | awk '{print $1}' | \
				awk -F '\t' 'NR==FNR{a[$0];next} $1 in a' - <(zcat ../samples/${i%.bam}_compressed.fasta.gz | awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | awk '{gsub(/>/,"");}1') | \
				awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t"); gsub(/>/,""); print}' | awk '{print ">"$1"-"$2"\n"$3}' | $gzip > ./haplotig/${i%.bam}_metagenome.fasta.gz &&
				wait
			fi
			) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			wait
			fi
		done
		wait
		find ./haplotig/ -size 0 -delete
		if [[ "$(ls ./haplotig/${i%.bam}_metagenome.fasta.gz | wc -l)" -gt 0 ]]; then rm *.bam "${projdir}"/samples/*_compressed.fasta.gz; fi
	fi
}
cd "${projdir}"
if [[ "$normalization" == false ]]; then
	metagout=$(ls ${projdir}/metagenome/haplotig/*metagenome.fasta.gz 2> /dev/null | wc -l)
	samplein=$(ls ${projdir}/samples/*.f 2> /dev/null | grep compressed | wc -l)
	if [[ "$metagout" -eq 0 ]]; then
		if [[ "$samplein" -eq "$metagout" ]];then
			echo -e "${YELLOW}- Qmatey is skipping normalization and only performing file compression ${WHITE}"
			time no_norm &>> ${projdir}/log.out
		else
			echo -e "${YELLOW}- Qmatey has already performed file compression ${WHITE}"
		fi
	fi
else
	echo -e "${YELLOW}- Qmatey has already performed file compression ${WHITE}"
fi

#################################################################################################################

cd "${projdir}"
if [[ "$taxids" == true ]]; then
	:> ${projdir}/metagenome/All.txids
	for i in ${projdir}/taxids/*.txids; do
		cat "$i" >> ${projdir}/metagenome/All.txids
	done

  sort -T "${projdir}"/tmp ${projdir}/metagenome/All.txids | uniq > ${projdir}/metagenome/All.tmp && mv "${projdir}"/metagenome/All.tmp ${projdir}/metagenome/All.txids
	wait
	awk 'NR>1{gsub(/\t\t/,"\tNA\t"); print}' ${Qmatey_dir}/tools/rankedlineage.dmp | awk '{gsub(/[|]/,""); print}' | awk '{gsub(/\t\t/,"\t"); print}' > ${projdir}/rankedlineage_tabdelimited.dmp &&
	awk -F'\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${projdir}/rankedlineage_tabdelimited.dmp ${projdir}/metagenome/All.txids | \
	cat <(printf "tax_id\ttaxname\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain\n") - > ${projdir}/rankedlineage_edited.dmp
	cat ${projdir}/rankedlineage_edited.dmp | grep -i 'Viruses' > ${projdir}/rankedlineage_edited_viruses.txt
	cat ${projdir}/rankedlineage_edited.dmp | grep -vi 'Viruses' > ${projdir}/rankedlineage_edited_other.txt
	cat ${projdir}/rankedlineage_edited_other.txt | grep -i 'uncultured\|unclassified\|unidentified\|Candidatus' > ${projdir}/rankedlineage_edited_other_uncultured.txt
	cat ${projdir}/rankedlineage_edited_other.txt| grep -iv 'uncultured\|unclassified\|unidentified\|Candidatus' > ${projdir}/rankedlineage_edited_other_cultured.txt
	awk -F'\t' '{print $2}' ${projdir}/rankedlineage_edited_other_cultured.txt | awk -F' ' '$1~/^[a-z]/{$0="NA"}1' | cut -f1,2 -d' ' | \
	awk '{gsub(/taxname/,"species");}1' | paste - ${projdir}/rankedlineage_edited_other_cultured.txt | \
	awk -F'\t' 'BEGIN{OFS="\t"} $4=="NA"{$4=$1}1' | awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_other_cultured_species.txt
	awk -F'\t' '{print $2}' ${projdir}/rankedlineage_edited_other_uncultured.txt | cut -f1,2,3 -d' ' | awk '{gsub(/taxname/,"species");}1' | paste - ${projdir}/rankedlineage_edited_other_uncultured.txt | \
	awk -F'\t' 'BEGIN{OFS="\t"} $4=="NA" && $1~/Candidatus/ {$4=$1}1' | \
	awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_other_uncultured_species.txt &&
	awk -F'\t' '{print $3}' ${projdir}/rankedlineage_edited_other_cultured_species.txt | cut -f1 -d' ' | \
	paste - ${projdir}/rankedlineage_edited_other_cultured_species.txt | awk -F'\t' 'BEGIN{OFS="\t"} $5=="NA"{$5=$1}1' | \
	awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_other_cultured_genus.txt
	awk -F'\t' '{print $3}' ${projdir}/rankedlineage_edited_other_uncultured_species.txt | cut -f1,2 -d' ' | \
	paste - ${projdir}/rankedlineage_edited_other_uncultured_species.txt | awk -F'\t' 'BEGIN{OFS="\t"} $5=="NA" && $1~/Candidatus/ {$5=$1}1' | \
	awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_other_uncultured_genus.txt
	cat ${projdir}/rankedlineage_edited_other_cultured_genus.txt ${projdir}/rankedlineage_edited_other_uncultured_genus.txt > ${projdir}/rankedlineage_edited_other_genus.txt
	mv "${projdir}"/rankedlineage_edited_other_genus.txt ${projdir}/rankedlineage_edited_final.txt
	rm "${projdir}"/rankedlineage_edited_other* &&
	awk -F'\t' '{print $2}' ${projdir}/rankedlineage_edited_viruses.txt | awk '{gsub(/taxname/,"species");}1' | paste - ${projdir}/rankedlineage_edited_viruses.txt | \
	awk -F'\t' 'BEGIN{OFS="\t"} $4=="NA"{$4=$1}1' | awk -F'\t' 'NR>1{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_viruses_species.txt
	cat ${projdir}/rankedlineage_edited_final.txt ${projdir}/rankedlineage_edited_viruses_species.txt > ${projdir}/temp
	mv "${projdir}"/temp ${projdir}/rankedlineage_edited.dmp
	rm "${projdir}"/rankedlineage_edited_viruses* ${projdir}/rankedlineage_edited_final.txt ${projdir}/rankedlineage_tabdelimited.dmp 2> /dev/null
else
	awk 'NR>1{gsub(/\t\t/,"\tNA\t"); print}' ${Qmatey_dir}/tools/rankedlineage.dmp | awk '{gsub(/[|]/,""); print}' | awk '{gsub(/\t\t/,"\t"); print}' | \
	cat <(printf "tax_id\ttaxname\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain\n") - > ${projdir}/rankedlineage_edited.dmp
	cat ${projdir}/rankedlineage_edited.dmp | grep -i 'Viruses' > ${projdir}/rankedlineage_edited_viruses.txt
	cat ${projdir}/rankedlineage_edited.dmp | grep -vi 'Viruses' > ${projdir}/rankedlineage_edited_other.txt
	cat ${projdir}/rankedlineage_edited_other.txt | grep -i 'uncultured\|unclassified\|unidentified\|Candidatus' > ${projdir}/rankedlineage_edited_other_uncultured.txt
	cat ${projdir}/rankedlineage_edited_other.txt| grep -iv 'uncultured\|unclassified\|unidentified\|Candidatus' > ${projdir}/rankedlineage_edited_other_cultured.txt
	awk -F'\t' '{print $2}' ${projdir}/rankedlineage_edited_other_cultured.txt | awk -F' ' '$1~/^[a-z]/{$0="NA"}1' | cut -f1,2 -d' ' | \
	awk '{gsub(/taxname/,"species");}1' | paste - ${projdir}/rankedlineage_edited_other_cultured.txt | \
	awk -F'\t' 'BEGIN{OFS="\t"} $4=="NA"{$4=$1}1' | awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_other_cultured_species.txt
	awk -F'\t' '{print $2}' ${projdir}/rankedlineage_edited_other_uncultured.txt | cut -f1,2,3 -d' ' | awk '{gsub(/taxname/,"species");}1' | paste - ${projdir}/rankedlineage_edited_other_uncultured.txt | \
	awk -F'\t' 'BEGIN{OFS="\t"} $4=="NA" && $1~/Candidatus/ {$4=$1}1' | \
	awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_other_uncultured_species.txt &&
	awk -F'\t' '{print $3}' ${projdir}/rankedlineage_edited_other_cultured_species.txt | cut -f1 -d' ' | \
	paste - ${projdir}/rankedlineage_edited_other_cultured_species.txt | awk -F'\t' 'BEGIN{OFS="\t"} $5=="NA"{$5=$1}1' | \
	awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_other_cultured_genus.txt
	awk -F'\t' '{print $3}' ${projdir}/rankedlineage_edited_other_uncultured_species.txt | cut -f1,2 -d' ' | \
	paste - ${projdir}/rankedlineage_edited_other_uncultured_species.txt | awk -F'\t' 'BEGIN{OFS="\t"} $5=="NA" && $1~/Candidatus/ {$5=$1}1' | \
	awk -F'\t' '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_other_uncultured_genus.txt
	cat ${projdir}/rankedlineage_edited_other_cultured_genus.txt ${projdir}/rankedlineage_edited_other_uncultured_genus.txt > ${projdir}/rankedlineage_edited_other_genus.txt
	mv "${projdir}"/rankedlineage_edited_other_genus.txt ${projdir}/rankedlineage_edited_final.txt
	rm "${projdir}"/rankedlineage_edited_other* &&
	awk -F'\t' '{print $2}' ${projdir}/rankedlineage_edited_viruses.txt | awk '{gsub(/taxname/,"species");}1' | paste - ${projdir}/rankedlineage_edited_viruses.txt | \
	awk -F'\t' 'BEGIN{OFS="\t"} $4=="NA"{$4=$1}1' | awk -F'\t' 'NR>1{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > ${projdir}/rankedlineage_edited_viruses_species.txt
	cat ${projdir}/rankedlineage_edited_final.txt ${projdir}/rankedlineage_edited_viruses_species.txt > ${projdir}/temp
	mv "${projdir}"/temp ${projdir}/rankedlineage_edited.dmp
	rm "${projdir}"/rankedlineage_edited_viruses* ${projdir}/rankedlineage_edited_final.txt
fi

lineagedb=${projdir}/lineage_subset.txt
if test -f $lineagedb; then
	echo -e "${YELLOW}- Compiling subset lineage file"
	if awk 'NR==1{/tax_id/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="tax_id") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > tax_id.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' tax_id.txt ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.dmp && rm tax_id.txt
	elif awk 'NR==1{/taxname/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="taxname") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > taxname.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' taxname.txt ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.dmp && rm taxname.txt
	elif awk 'NR==1{/genus/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="genus") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > genus.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$4] > 0' genus.txt ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.dmp && rm genus.txt
	elif awk 'NR==1{/family}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="family") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > family.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$5] > 0' family.txt ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.dmp && rm family.txt
	elif awk 'NR==1{/order/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="order") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > order.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$6] > 0' order.txt ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.dmp && rm order.txt
	elif awk 'NR==1{/class/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="class") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > class.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$7] > 0' class.txt ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.dmp && rm class.txt
	elif awk 'NR==1{/phylum/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="phylum") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > phylum.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$8] > 0' phylum.txt ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.dmp && rm phylum.txt
	elif awk 'NR==1{/kingdom/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="kingdom") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > kingdom.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$9] > 0' kingdom.txt ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.dmp && rm kingdom.txt
	elif awk 'NR==1{/domain/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="domain") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > domain.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$10] > 0' domain.txt ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.dmp && rm domain.txt
	else
		echo -e "${magenta}- lineage database is not formatted properly ${white}\n"
		echo -e "${magenta}- Do you wish to continue running Qmatey? ${white}\n"
		read -p "- y(YES) or n(NO) " -n 1 -r
		if [[ ! $REPLY =~ ^[Yy]$ ]]; then
			printf '\n'
			exit 0
		fi
	fi
fi

if test -f "${projdir}/remove_taxa.txids"; then
  awk 'NR==FNR {a[$1]} FNR!=NR && !($1 in a)' ${projdir}/remove_taxa.txids ${projdir}/rankedlineage_edited.dmp > ${projdir}/rankedlineage_edited.temp &&
  mv ${projdir}/rankedlineage_edited.temp ${projdir}/rankedlineage_edited.dmp &&
  wait
fi


#################################################################################################################
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey MegaBLAST \n\e[97m########################################################\n"

cd ${projdir}
local_db=$(awk '{gsub(/,/," ")}1' <<< "$local_db")
local_dblist=$local_db
if [[ "$local_db" =~ " " ]]; then
	local_db=$(awk '{gsub(/ /,"~ ~");}1' <<< $local_db | tr '~' '"')
	local_db=$(echo \"$local_db\")
	local_dblistdir=$(echo ${local_dblist} | tr ' ' '\n' | head -n1)
	local_dblistdir=${local_dblistdir%/*}
	dblists=$(awk -v dir=$local_dblistdir '{gsub(dir,"");gsub(/ /,"_");gsub(/\//,"");}1' <<< ${local_dblist})
	printf "TITLE all my databases\n" > ${local_dblistdir}/${dblists}.nal
	printf "DBLIST $local_db \n" | awk -v dir=$local_dblistdir '{gsub(dir,"");gsub(/\//,"");}1' >> ${local_dblistdir}/${dblists}.nal
	local_db=${local_dblistdir}/${dblists}
fi
echo -e "${magenta}- database: $local_dblist ${white}\n"


if (echo $local_db | grep -q 'nt'); then
	if [[ -z $percid ]]; then
		export percid=95
	fi
fi
if (echo $local_db | grep -q 'ref'); then
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
export rpm=$(($reads_per_megablast * 2))
export rpmb=$(($reads_per_megablast_burn_in * 2))



blast () {

if [[ "$fastMegaBLAST" == true ]]; then
	cd "${projdir}"/metagenome/haplotig
	if test ! -f combined_compressed_metagenomes.fasta.gz; then
		# find . -name "*.fasta.gz" | xargs zcat | grep -v '^>' | awk '{A[$1]++}END{for(i in A)print i}' | awk '{print ">"NR"\n"$1}' > combined_compressed_metagenomes.fasta
		if [[ "$(ls *_compressed.fasta.gz 2>/dev/null | wc -l)" -gt 0 ]]; then
			for i in $(ls *_compressed.fasta.gz); do (
				zcat $i | awk '{gsub(/AAAAAAAAAA$/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");}1' | \
				awk '{gsub(/CCCCCCCCCC$/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");}1' | \
				awk '{gsub(/GGGGGGGGGG$/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");}1' | \
				awk '{gsub(/TTTTTTTTTT$/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");}1' | \
				awk '{gsub(/^AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");}1' | \
				awk '{gsub(/^CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");}1' | \
				awk '{gsub(/^GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");}1' | \
				awk '{gsub(/^TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");}1' | \
				awk '{gsub(/AAAAA~/,"~");gsub(/CCCCC~/,"~");gsub(/GGGGG~/,"~");gsub(/TTTTT~/,"~");gsub(/~AAAAA/,"~");gsub(/~CCCCC/,"~");gsub(/~GGGGG/,"~");gsub(/~TTTTT/,"~");}1' | \
				awk '{gsub(/AA~/,"~");gsub(/CC~/,"~");gsub(/GG~/,"~");gsub(/TT~/,"~");gsub(/~AA/,"~");gsub(/~CC/,"~");gsub(/~GG/,"~");gsub(/~TT/,"~");gsub(/AA~/,"~");gsub(/CC~/,"~");gsub(/GG~/,"~");gsub(/TT~/,"~");gsub(/~AA/,"~");gsub(/~CC/,"~");gsub(/~GG/,"~");gsub(/~TT/,"~");}1' | \
				awk '{gsub(/A~/,"~");gsub(/C~/,"~");gsub(/G~/,"~");gsub(/T~/,"~");gsub(/~A/,"~");gsub(/~C/,"~");gsub(/~G/,"~");gsub(/~T/,"~");gsub(/~/,"");}1' | gzip > ${i}.tmp &&
				mv ${i}.tmp $i ) &
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
					wait
				fi
			done
			for i in $(ls *_compressed.fasta.gz); do
				zcat $i | grep -v '^>' | awk '{A[$1]++}END{for(i in A)print i}' | gzip >> combined_compressed_metagenomes.tmp.gz &&
				mv $i ${i%*_compressed.fasta.gz}_metagenome.fasta.gz &&
				wait
			done
		fi
		if [[ "$(ls *_metagenome.fasta.gz 2>/dev/null | wc -l)" -gt 0 ]]; then
			for i in $(ls *_metagenome.fasta.gz); do
				zcat $i | grep -v '^>' | awk '{A[$1]++}END{for(i in A)print i}' | gzip >> combined_compressed_metagenomes.tmp.gz
				wait
			done
		fi
		wait
		awk '{A[$1]++}END{for(i in A)print i}' <(zcat combined_compressed_metagenomes.tmp) | gzip > combined_compressed_metagenomes_full.tmp.gz
		rm combined_compressed_metagenomes.tmp.gz
		LC_ALL=C sort -T "${projdir}"/tmp <(zcat ombined_compressed_metagenomes_full.tmp.gz) | gzip > combined_compressed_metagenomes_tmp.fasta.gz
		wait
		rm combined_compressed_metagenomes_full.tmp.gz
		awk '{print ">"NR"\n"$1}' <(zcat combined_compressed_metagenomes_tmp.fasta) | $gzip > combined_compressed_metagenomes.fasta.gz
		wait
		rm combined_compressed_metagenomes_tmp.fasta.gz
	fi

	if [[ "$taxids" == true ]]; then
		for i in ${projdir}/taxids/*.txids; do
			cat "$i" >> ${projdir}/metagenome/All.txids
		done

		wait
		# awk -F'\t' '!seen[$4]++' ${projdir}/rankedlineage_edited.dmp | awk '{print $1}' |
		awk -F'\t' '$4!=p {if(p)print ""; p=$4; l=NR+0} NR<=l; END {print ""}' ${projdir}/rankedlineage_edited.dmp | awk '{print $1}' | \
		grep . | awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${projdir}/metagenome/All.txids - > ${projdir}/metagenome/burnin.txids
		wait
		awk -F'\t' '$4!=p {if(p)print ""; p=$4; l=NR+5} NR<=l; END {print ""}' ${projdir}/rankedlineage_edited.dmp | awk '{print $1}' | \
		grep . | awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ${projdir}/metagenome/All.txids - | \
		awk 'FNR==NR {f2[$1];next} !($0 in f2)' ${projdir}/metagenome/burnin.txids - > ${projdir}/metagenome/burnin2.txids
		wait
	fi
	if [[ -n "$(ls ${projdir}/metagenome/alignment/subfile* 2> /dev/null)" ]]; then
		mv ${projdir}/metagenome/alignment ${projdir}/metagenome/alignment_w_subfile
	fi

	if [[ "$blast_location" =~ "local" ]]; then
		echo -e "${YELLOW}- performing local BLAST"
		if [[ -z "$(ls -R ${projdir}/metagenome/alignment/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]] && [[ -z "$(ls -R ${projdir}/metagenome/alignment/cultured/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]; then

			if [[ "$reads_per_megablast_burn_in" -gt 0 ]]; then
			  if test ! -f ${projdir}/metagenome/haplotig/repseq_aln/repseq.txids ; then
			    if [[ -d repseq_aln ]]; then
			      cd repseq_aln
			    else
			      mkdir repseq_aln; cd repseq_aln
			    fi
			    if test ! -f "representative_seq_retireved.txt"; then
						longestread=$(zcat ../combined_compressed_metagenomes.fasta.gz 2> /dev/null | wc -L | awk '{print $1}')
						longestread=$((longestread*95/100))
			      awk 'NR%50000==1{close("combined_compressed_metagenomes_sub"i); i++}{print > "combined_compressed_metagenomes_sub"i}'  <(zcat ../combined_compressed_metagenomes.fasta.gz 2> /dev/null) & PIDsplit1=$!
			      wait $PIDsplit1
			      wait
			      for i in combined_compressed_metagenomes_sub*; do (
							$cd_hit_est -i "$i" -o "${i}"_clust -c 0.95 -n 10 -g 0 -M 500 -T 1 -A $longestread
			        wait
			        rm "$i" "${i}"_clust.clstr )&
			        if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
			          wait
			        fi
			      done
			      wait
			      :> combined_compressed_metagenomes_repseq.txt
			      for i in $(ls *clust); do
			        grep -v '>' $i >> combined_compressed_metagenomes_repseq.txt &&
			        wait
			      done
			      wait
			      repseq_lineN=$(wc -l combined_compressed_metagenomes_repseq.txt | awk '{print $1}')
			      cat combined_compressed_metagenomes_repseq.txt | LC_ALL=C sort -T "${projdir}"/tmp | awk '{print ">"NR"\n"$1}' | gzip > combined_compressed_metagenomes_repseq.fasta.gz &&
			      rm *clust combined_compressed_metagenomes_repseq.txt
			      wait
			      printf "extraction of representative metagenome reads completed (using cdhit to cluster sequences)\n" > representative_seq_retireved.txt &&
			      wait
			    fi
			    if [[ -d splitccf ]]; then
			      cd splitccf
			    else
			      mkdir splitccf &&
						cd splitccf &&
						mkdir alignment &&
			      awk 'NR%2000000==1{close("F"i); i++}{print > "F"i}'  <(zcat ../combined_compressed_metagenomes_repseq.fasta.gz 2> /dev/null) & PIDsplit1=$!
			      wait $PIDsplit1
			      wait
			    fi
			    if [[ -n "$(ls ./alignment/subfile* 2> /dev/null)" ]]; then
			      ls ./alignment/subfile* 2> /dev/null | xargs rm
			      mv ./alignment/F* ./
			    fi
					cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf/

			    for ccf in $(ls * | grep -v alignment | sort -V); do
			      mv $ccf ${projdir}/metagenome/haplotig/repseq_aln/splitccf/alignment/
			      cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf/alignment/
			      awk -v pat="$ccf" -v rpmb="$rpmb" 'NR%rpmb==1{close(pat"_subfile"pat"_"i); i++}{print > pat"_subfile"pat"_"i}' $ccf & PIDsplit2=$!
			      wait $PIDsplit2
						if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt "$threads" ]]; then
							Nsize=$(ls *subfile* 2> /dev/null | wc -l)
							nbatch=$(($Nsize / $threads))
							nsize=$(($nbatch * $threads))
							for subf in $(ls *subfile* 2> /dev/null | sort -V | head -n $nsize); do mv $subf ${subf#*_}; done
						fi
						if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -le "$threads" ]]; then
							for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
						fi
			      for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
							if test ! -f "${sub}_out.blast.gz"; then
				        if [[ "$taxids" == true ]]; then
				          ${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity 95 -max_target_seqs $max_target -evalue 0.01 \
				          -word_size 28 -taxidlist ${projdir}/metagenome/burnin.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
				          wait
				          if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				          wait
				        else
				          ${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity 95 -max_target_seqs $max_target -evalue 0.01 \
				          -word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
				          wait
				          if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				          wait
				        fi
								wait
								gzip "${sub}_out.blast" &&
								cat
								wait
							fi
							wait
			        rm "$sub" )&
			        if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
			          wait
			        fi
			      done
			      wait
			      rm $ccf &&
			      for subfile in *_out.blast.gz; do
			        cat "$subfile" >> ${ccf}.blast.gz &&
			        rm "$subfile" &&
			        wait
			      done
			      wait
			      cat "${ccf}".blast.gz >> ../combined_compressed.megablast.gz &&
			      rm "${ccf}".blast.gz &&
			      cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf
			      wait
			    done
			    wait

					cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf/alignment/
					if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt 0 ]]; then
						for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
						for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
							if test ! -f "${sub}_out.blast.gz"; then
								if [[ "$taxids" == true ]]; then
									${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity 95 -max_target_seqs $max_target -evalue 0.01 \
									-word_size 28 -taxidlist ${projdir}/metagenome/burnin.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
									wait
									if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
									wait
								else
									${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity 95 -max_target_seqs $max_target -evalue 0.01 \
									-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
									wait
									if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
									wait
								fi
								wait
								gzip "${sub}_out.blast" &&
								cat
								wait
							fi
							wait
							rm "$sub" )&
							if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
								wait
							fi
						done
						wait
						for subfile in *_out.blast.gz; do
							cat "$subfile" >> ${ccf}.blast.gz &&
							rm "$subfile" &&
							wait
						done
						wait
						cat "${ccf}".blast.gz >> ../combined_compressed.megablast.gz &&
						rm "${ccf}".blast.gz &&
						wait
					fi
					cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf
					wait

			    cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf &&
					zcat combined_compressed.megablast.gz | awk '{print $9}' | awk '{gsub(/;/,"\n");}1' | sort | uniq > ${projdir}/metagenome/repseq.txids &
					awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' <(zcat ../combined_compressed_metagenomes_repseq.fasta.gz) | gzip > repseq_linear.fasta.gz &&
					zcat combined_compressed.megablast.gz | awk '{print ">"$1}' | awk '!a[$1]++' | grep -f - -v <(zcat repseq_linear.fasta.gz) | \
					awk '{gsub(/\t/,"\n");}1' | gzip -f > ../combined_compressed_metagenomes_repseq.fasta.gz &&
					rm -rf repseq_linear.fasta.gz alignment &&
					wait

					if [[ "$(zcat ../combined_compressed_metagenomes_repseq.fasta.gz | head | wc -l)" -gt 1 ]]; then
						awk 'NR%2000000==1{close("F"i); i++}{print > "F"i}'  <(zcat ../combined_compressed_metagenomes_repseq.fasta.gz 2> /dev/null) & PIDsplit1=$!
						wait $PIDsplit1 &&
						mkdir alignment &&
						wait
						for ccf in $(ls * | grep -v 'alignment\|mega' | sort -V); do
							mv $ccf ${projdir}/metagenome/haplotig/repseq_aln/splitccf/alignment/
							cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf/alignment/
							awk -v pat="$ccf" -v rpmb="$rpmb" 'NR%rpmb==1{close(pat"_subfile"pat"_"i); i++}{print > pat"_subfile"pat"_"i}' $ccf & PIDsplit2=$!
							wait $PIDsplit2
							if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt "$threads" ]]; then
								Nsize=$(ls *subfile* 2> /dev/null | wc -l)
								nbatch=$(($Nsize / $threads))
								nsize=$(($nbatch * $threads))
								for subf in $(ls *subfile* 2> /dev/null | sort -V | head -n $nsize); do mv $subf ${subf#*_}; done
							fi
							if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -le "$threads" ]]; then
								for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
							fi
							for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
								if test ! -f "${sub}_out.blast.gz"; then
									if [[ "$taxids" == true ]]; then
										${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity 90 -max_target_seqs $max_target -evalue 0.01 \
										-word_size 28 -taxidlist ${projdir}/metagenome/burnin2.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
										wait
										if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
										wait
									else
										${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity 90 -max_target_seqs $max_target -evalue 0.01 \
										-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
										wait
										if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
										wait
									fi
									wait
									gzip "${sub}_out.blast" &&
									wait
								fi
								wait
								rm "$sub" )&
								if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
									wait
								fi
							done
							wait
							rm $ccf &&
							for subfile in *_out.blast.gz; do
								cat "$subfile" >> ${ccf}.blast.gz &&
								rm "$subfile" &&
								wait
							done
							wait
							cat "${ccf}".blast.gz >> ../combined_compressed.megablast.gz &&
							rm "${ccf}".blast.gz &&
							cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf
							wait
						done
						wait

						cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf/alignment/
						if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt 0 ]]; then
							for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
							for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
								if test ! -f "${sub}_out.blast.gz"; then
									if [[ "$taxids" == true ]]; then
										${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity 90 -max_target_seqs $max_target -evalue 0.01 \
										-word_size 28 -taxidlist ${projdir}/metagenome/burnin2.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
										wait
										if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
										wait
									else
										${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity 90 -max_target_seqs $max_target -evalue 0.01 \
										-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
										wait
										if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
										wait
									fi
									wait
									gzip "${sub}_out.blast" &&
									wait
								fi
								wait
								rm "$sub" )&
								if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
									wait
								fi
							done
							wait
							for subfile in *_out.blast.gz; do
								cat "$subfile" >> ${ccf}.blast.gz &&
								rm "$subfile" &&
								wait
							done
							wait
							cat "${ccf}".blast.gz >> ../combined_compressed.megablast.gz &&
							rm "${ccf}".blast.gz &&
							wait
						fi
						cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf
						wait

						cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf
						zcat combined_compressed.megablast.gz | awk '{print $9}' | awk '{gsub(/;/,"\n");}1' | sort | uniq >> ${projdir}/metagenome/repseq2.txids &&
						wait
					fi

					cd ${projdir}/metagenome/haplotig/repseq_aln/splitccf
			    cp ${projdir}/metagenome/repseq.txids ${projdir}/metagenome/haplotig/repseq_aln/
			    cd ${projdir}/metagenome/ &&
			    mv All.txids matched_and_unmatched.txids &&
					awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' repseq.txids ${projdir}/rankedlineage_edited.dmp | awk -F'\t' '{print $4}' | \
					awk -F'\t' 'NR==FNR{c[$1]++;next};c[$4] > 0' - ${projdir}/rankedlineage_edited.dmp | awk '{print $1}' > All1.txids &&
					awk -F'\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' repseq2.txids ${projdir}/rankedlineage_edited.dmp | awk -F'\t' '{print $4}' | \
					awk -F'\t' 'NR==FNR{c[$1]++;next};c[$4] > 0' - ${projdir}/rankedlineage_edited.dmp | awk '{print $1}' >> All2.txids &&
					cat All1.txids All2.txids | awk '!a[$1]++' > All.txids &&
			    rm ${projdir}/metagenome/haplotig/repseq_aln/combined_compressed_metagenomes_repseq.fasta.gz All1.txids All2.txids repseq.txids repseq2.txids &&
			    wait
			  fi
			fi

			rm -rf ${projdir}/metagenome/burnin*.txids
			cd ${projdir}/metagenome/haplotig
			mv ${projdir}/metagenome/alignment_w_subfile ${projdir}/metagenome/alignment 2> /dev/null
			if [[ -d splitccf ]]; then
				cd splitccf
			else
				mkdir splitccf; cd splitccf
				cp ../combined_compressed_metagenomes.fasta.gz ./combined_compressed_metagenomes.fasta.gz
				awk 'NR%2000000==1{close("F"i); i++}{print > "F"i}'  <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) & PIDsplit1=$!
				wait $PIDsplit1
				rm combined_compressed_metagenomes.fasta.gz
			fi
			if [[ -n "$(ls ../../alignment/subfile* 2> /dev/null)" ]]; then
				ls ../../alignment/subfile* 2> /dev/null | xargs rm
				mv ../../alignment/F* ./
			fi
			if [[ $nodes -gt 1 ]]; then
				splitnumt=$(ls F* | wc -l) && splitnum=$(($splitnumt / $nodes))
				start_split=0
				for nn in $(seq 1 "$nodes"); do
					mkdir -p splitccf_node${nn}
					end_split=$(($start_split + $splitnum))
					start_split=$(($start_split + 1))
					for snode in $(seq $start_split $end_split); do mv F${snode} ./splitccf_node${nn}/; done
					wait
					start_split=$end_split
				done
				wait
				for nn in $(seq 1 "$nodes"); do start_split=$(($start_split + 1)) && mv F${start_split} ./splitccf_node${nn}/ 2> /dev/null; done
				wait
				touch ${projdir}/multi_node_run_ready.txt

				echo -e "${YELLOW}- performing local BLAST in multi-node mode"
				cd "${projdir}"/metagenome/haplotig/splitccf/splitccf_node1
				for ccf in $(ls * | grep -v alignment | sort -V); do
					mv $ccf ${projdir}/metagenome/alignment/$ccf
					cd "${projdir}"/metagenome/alignment
					awk -v pat="$ccf" -v rpm="$rpm" 'NR%rpm==1{close(pat"_subfile"pat"_"i); i++}{print > pat"_subfile"pat"_"i}' $ccf & PIDsplit2=$!
					wait $PIDsplit2
					if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt "$threads" ]]; then
						Nsize=$(ls *subfile* 2> /dev/null | wc -l)
						nbatch=$(($Nsize / $threads))
						nsize=$(($nbatch * $threads))
						for subf in $(ls *subfile* 2> /dev/null | sort -V | head -n $nsize); do mv $subf ${subf#*_}; done
					fi
					if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -le "$threads" ]]; then
						for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
					fi
					for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
						if test ! -f "${sub}_out.blast.gz"; then
							if [[ "$taxids" == true ]]; then
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							else
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							fi
							wait
							gzip "${sub}_out.blast" &&
							wait
						fi
						wait
						rm "$sub" )&
						if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
							wait
						fi
					done
					wait
					for subfile in *_out.blast.gz; do
						cat "$subfile" >> ${ccf}.blast.gz &&
						rm "$subfile" && wait
					done
					wait
					cat ${ccf}.blast.gz >> combined_compressed_node1.megablast.gz &&
					rm "${ccf}".blast.gz; rm "$ccf" &&
					cd "${projdir}"/metagenome/haplotig/splitccf/splitccf_node1
				done
				wait

				cd "${projdir}"/metagenome/alignment
				if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt 0 ]]; then
					for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
					for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
						if test ! -f "${sub}_out.blast.gz"; then
							if [[ "$taxids" == true ]]; then
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							else
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							fi
							wait
							gzip "${sub}_out.blast" &&
							wait
						fi
						wait
						rm "$sub" )&
						if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
							wait
						fi
					done
					wait
					for subfile in *_out.blast.gz; do
						cat "$subfile" >> ${ccf}.blast.gz &&
						rm "$subfile" && wait
					done
					wait
					cat ${ccf}.blast.gz >> combined_compressed_node1.megablast.gz &&
					rm "${ccf}".blast.gz &&
					wait
				fi
				cd "${projdir}"/metagenome/haplotig/splitccf/splitccf_node1
				wait

			else
				for ccf in $(ls * | grep -v alignment | sort -V); do
					mv $ccf ../../alignment/$ccf
					cd ../../alignment
					awk -v pat="$ccf" -v rpm="$rpm" 'NR%rpm==1{close(pat"_subfile"pat"_"i); i++}{print > pat"_subfile"pat"_"i}' $ccf & PIDsplit2=$!
					wait $PIDsplit2
					if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt "$threads" ]]; then
						Nsize=$(ls *subfile* 2> /dev/null | wc -l)
						nbatch=$(($Nsize / $threads))
						nsize=$(($nbatch * $threads))
						for subf in $(ls *subfile* 2> /dev/null | sort -V | head -n $nsize); do mv $subf ${subf#*_}; done
					fi
					if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -le "$threads" ]]; then
						for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
					fi
					for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
						if test ! -f "${sub}_out.blast.gz"; then
							if [[ "$taxids" == true ]]; then
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							else
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							fi
							wait
							gzip "${sub}_out.blast" &&
							wait
						fi
						wait
						rm "$sub" )&
						if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
							wait
						fi
					done
					wait
					for subfile in *_out.blast.gz; do
						cat "$subfile" >> ${ccf}.blast.gz &&
						rm "$subfile" && wait
					done
					wait
					# zcat ${ccf}.blast.gz 2> /dev/null | awk -v percid=$percid '$3 >= $5*(percid/100) {print $0}' | $gzip >> combined_compressed.megablast.gz &&
					cat "${ccf}".blast.gz >> combined_compressed.megablast.gz &&
					rm "${ccf}".blast.gz; rm "$ccf" &&
					cd ../haplotig/splitccf/
				done
				wait

				cd "${projdir}"/metagenome/alignment
				if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt 0 ]]; then
					for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
					for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
						if test ! -f "${sub}_out.blast.gz"; then
							if [[ "$taxids" == true ]]; then
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							else
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							fi
							wait
							gzip "${sub}_out.blast" &&
							wait
						fi
						wait
						rm "$sub" )&
						if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
							wait
						fi
					done
					wait
					for subfile in *_out.blast.gz; do
						cat "$subfile" >> ${ccf}.blast.gz &&
						rm "$subfile" && wait
					done
					wait
					# zcat ${ccf}.blast.gz 2> /dev/null | awk -v percid=$percid '$3 >= $5*(percid/100) {print $0}' | $gzip >> combined_compressed.megablast.gz &&
					cat "${ccf}".blast.gz >> combined_compressed.megablast.gz &&
					rm "${ccf}".blast.gz &&
					wait
				fi
				cd ../haplotig/splitccf/

			fi
			if [[ $nodes -gt 1 ]]; then
				count_megablast_node=$(ls ${projdir}/megablast_done_node*.txt 2> /dev/null | wc -l)
				while [[ "$count_megablast_node" < $(($nodes - 1)) ]]; do
					sleep 300
				done
				wait
				cd "${projdir}"/metagenome/alignment
				cat combined_compressed_node*.megablast.gz > combined_compressed.megablast.gz
				rm combined_compressed_node*.megablast.gz
				wait
			fi
			cd "${projdir}"/metagenome/haplotig/splitccf/
			rmdir * 2> /dev/null
			rm "${projdir}"/megablast_node* ${projdir}/multi_node_run_ready.txt ${projdir}/megablast_splitrun_node_*.sh 2> /dev/null
			cd ../
			rmdir splitccf 2> /dev/null
		else
			echo -e "${YELLOW}- Primary BLAST ouput already exist"
			echo -e "${YELLOW}- Skipping BLAST and filtering hits based on defined parameters"
		fi
		wait

		cd "${projdir}"/metagenome/alignment/
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA ==  false ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -i $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -vi $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA == true ]] ; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -i $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			mkdir rRNA
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -i 'rRNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/rRNA/rRNA_combined_compressed.megablast.gz
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -vi $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi
		wait

		cd "${projdir}"/metagenome/haplotig
		for i in $(ls *metagenome.fasta.gz); do (
		  if test ! -f ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz && test ! -f ../alignment/cultured/${i%_metagenome.fasta.gz}_haplotig.megablast.gz; then
				awk '!/^$/' <(zcat "$i" 2> /dev/null) | awk -F'\t' 'ORS=NR%2?"\t":"\n"' | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz &&
				wait
				awk -F'\t' 'ORS=NR%2?"\t":"\n"' <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | \
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' - <(zcat ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $3,$1}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz &&
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&

				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/uncultured_combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/uncultured_${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&
				wait
			  # if [[ "$taxids" = true ]]; then
			  #   for taxid_files in $(ls ${projdir}/taxids/*.txids); do
			  #     taxid=${taxid_files%*.txids}
			  #     taxid=${taxid/*\/}
			  #     awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($8) in a{print $0, a[$1]}' $taxid_files <(zcat ../alignment/${i%_metagenome.fasta.gz}_haplotig_temp.megablast.gz 2> /dev/null) | $gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid${taxid}.megablast.gz
			  #   done
			  #   wait
			  #   rm ../alignment/"${i%_metagenome.fasta.gz}"_haplotig_temp.megablast.gz
			  #   find ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid*.megablast.gz | xargs cat > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz &&
			  #   rm ../alignment/"${i%_metagenome.fasta.gz}"_haplotig_taxid*.megablast.gz
			  # fi
				rm ../alignment/"${i%_metagenome.fasta.gz}"_step*.txt.gz
			fi )&
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait

	fi

	if [[ "$blast_location" =~ "remote" ]]; then
		echo -e "${YELLOW}- performing a remote BLAST"
		if [[ -z "$(ls -R ${projdir}/metagenome/alignment/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]  && [[ -z "$(ls -R ${projdir}/metagenome/alignment/cultured/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]; then
			mv ${projdir}/metagenome/alignment_w_subfile ${projdir}/metagenome/alignment 2> /dev/null
			if [[ "$taxids" == true ]]; then
				${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) -db $blast_location -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
				-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" \
				-out ../alignment/combined_compressed.megablast -remote
				wait
				if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				wait
			else
				${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) -db $blast_location -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
				-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" \
				-out ../alignment/combined_compressed.megablast -remote
				wait
				if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				wait
			fi
		else
			echo -e "${YELLOW}- BLAST ouput already exist"
			echo -e "${YELLOW}- Skipping BLAST and filtering hits based on defined parameters"
		fi
		wait

		if test -f ../alignment/combined_compressed.megablast; then
			$gzip ../alignment/combined_compressed.megablast &&
			rm ../alignment/combined_compressed.megablast
		fi

		cd "${projdir}"/metagenome/alignment/
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA ==  false ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -i $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -vi $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA == true ]] ; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -i $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			mkdir rRNA
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -i 'rRNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/rRNA/rRNA_combined_compressed.megablast.gz
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -vi $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi
		wait



		for i in $(ls *metagenome.fasta.gz); do (
			if test ! -f ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz && test ! -f ../alignment/cultured/${i%_metagenome.fasta.gz}_haplotig.megablast.gz; then
				awk '!/^$/' <(zcat "$i" 2> /dev/null) | awk -F'\t' 'ORS=NR%2?"\t":"\n"' | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz &&
				wait
				awk -F'\t' 'ORS=NR%2?"\t":"\n"' <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | \
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' - <(zcat ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $3,$1}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz &&
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&

				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/uncultured_combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/uncultured_${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&
				wait
				rm ../alignment/"${i%_metagenome.fasta.gz}"_step*.txt
			fi )&
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait
	fi

	if [[ "$blast_location" =~ "custom" ]]; then
		echo -e "${YELLOW}- performing custom BLAST"
		if [[ -z "$(ls -R ${projdir}/metagenome/alignment/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]  && [[ -z "$(ls -R ${projdir}/metagenome/alignment/cultured/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]; then
			mv ${projdir}/metagenome/alignment_w_subfile ${projdir}/metagenome/alignment 2> /dev/null
			if [[ -d splitccf ]]; then
			  cd splitccf
			else
			  mkdir splitccf; cd splitccf
			  cp ../combined_compressed_metagenomes.fasta.gz ./combined_compressed_metagenomes.fasta.gz
			  awk 'NR%2000000==1{close("F"i); i++}{print > "F"i}'  <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) & PIDsplit1=$!
			  wait $PIDsplit1
			  rm combined_compressed_metagenomes.fasta.gz
			fi
			if [[ -n "$(ls ../../alignment/subfile* 2> /dev/null)" ]]; then
			  ls ../../alignment/subfile* 2> /dev/null | grep -v .gz$ | xargs rm
			  mv ../../alignment/F* ./
			fi
			if [[ $nodes -gt 1 ]]; then
				splitnumt=$(ls F* | wc -l) && splitnum=$(($splitnumt / $nodes))
				start_split=0
				for nn in $(seq 1 "$nodes"); do
					mkdir -p splitccf_node${nn}
					end_split=$(($start_split + $splitnum))
					start_split=$(($start_split + 1))
					for snode in $(seq $start_split $end_split); do mv F${snode} ./splitccf_node${nn}/; done
					start_split=$end_split
				done
				wait
				for nn in $(seq 1 "$nodes"); do start_split=$(($start_split + 1)) && mv F${start_split} ./splitccf_node${nn}/ 2> /dev/null; done
				wait
				touch ${projdir}/multi_node_run_ready.txt

				echo -e "${YELLOW}- performing custom BLAST in multi-node mode"
				cd "${projdir}"/metagenome/haplotig/splitccf/splitccf_node1
				for ccf in $(ls * | grep -v alignment | sort -V); do
					mv $ccf ${projdir}/metagenome/alignment/$ccf
					cd "${projdir}"/metagenome/alignment
					awk -v pat="$ccf" -v rpm="$rpm" 'NR%rpm==1{close(pat"_subfile"pat"_"i); i++}{print > pat"_subfile"pat"_"i}' $ccf & PIDsplit2=$!
					wait $PIDsplit2
					if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt "$threads" ]]; then
						Nsize=$(ls *subfile* 2> /dev/null | wc -l)
						nbatch=$(($Nsize / $threads))
						nsize=$(($nbatch * $threads))
						for subf in $(ls *subfile* 2> /dev/null | sort -V | head -n $nsize); do mv $subf ${subf#*_}; done
					fi
					if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -le "$threads" ]]; then
						for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
					fi
					for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
						if test ! -f "${sub}_out.blast.gz"; then
							if [[ "$taxids" == true ]]; then
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							else
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							fi
							wait
							gzip "${sub}_out.blast" &&
							wait
						fi
						wait
						rm "$sub" )&
						if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
							wait
						fi
					done
					wait
					for subfile in *_out.blast.gz; do
						cat "$subfile" >> ${ccf}.blast.gz &&
						rm "$subfile" && wait
					done
					wait
					cat "${ccf}".blast.gz >> combined_compressed_node1.megablast.gz &&
					rm "${ccf}".blast.gz; rm "$ccf" &&
					cd "${projdir}"/metagenome/haplotig/splitccf/splitccf_node1
				done
				wait

				if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt 0 ]]; then
					for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
					for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
						if test ! -f "${sub}_out.blast.gz"; then
							if [[ "$taxids" == true ]]; then
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							else
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							fi
							wait
							gzip "${sub}_out.blast" &&
							wait
						fi
						wait
						rm "$sub" )&
						if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
							wait
						fi
					done
					wait
					for subfile in *_out.blast.gz; do
						cat "$subfile" >> ${ccf}.blast.gz &&
						rm "$subfile" && wait
					done
					wait
					cat "${ccf}".blast.gz >> combined_compressed_node1.megablast.gz &&
					rm "${ccf}".blast.gz &&
					wait
				fi
				cd "${projdir}"/metagenome/haplotig/splitccf/splitccf_node1
				wait

			else
				for ccf in $(ls * | grep -v alignment | sort -V); do
					mv $ccf ../../alignment/$ccf
					cd ../../alignment
					awk -v pat="$ccf" -v rpm="$rpm" 'NR%rpm==1{close(pat"_subfile"pat"_"i); i++}{print > pat"_subfile"pat"_"i}' $ccf & PIDsplit2=$!
					wait $PIDsplit2
					if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt "$threads" ]]; then
						Nsize=$(ls *subfile* 2> /dev/null | wc -l)
						nbatch=$(($Nsize / $threads))
						nsize=$(($nbatch * $threads))
						for subf in $(ls *subfile* 2> /dev/null | sort -V | head -n $nsize); do mv $subf ${subf#*_}; done
					fi
					if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -le "$threads" ]]; then
						for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
					fi
					for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
						if test ! -f "${sub}_out.blast.gz"; then
							if [[ "$taxids" == true ]]; then
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							else
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							fi
							wait
							gzip "${sub}_out.blast" &&
							wait
						fi
						rm "$sub" )&
						if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
							wait
						fi
					done
					wait
					for subfile in *_out.blast.gz; do
						cat "$subfile" >> ${ccf}.blast.gz &&
						rm "$subfile" && wait
					done
					wait
					# zcat ${ccf}.blast.gz 2> /dev/null | awk -v percid=$percid '$3 >= $5*(percid/100) {print $0}' | $gzip >> combined_compressed.megablast.gz &&
					cat "${ccf}".blast.gz >> combined_compressed.megablast.gz &&
					rm "${ccf}".blast.gz; rm "$ccf" &&
					cd ../haplotig/splitccf/
				done
				wait

				cd "${projdir}"/metagenome/alignment
				if [[ "$(ls *subfile* 2> /dev/null | wc -l)" -gt 0 ]]; then
					for subf in $(ls *subfile* 2> /dev/null | sort -V); do mv $subf ${subf#*_}; done
					for sub in $(ls subfile* 2> /dev/null | sort -T "${projdir}"/tmp -V); do (
						if test ! -f "${sub}_out.blast.gz"; then
							if [[ "$taxids" == true ]]; then
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							else
								${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query $sub -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
								-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out "${sub}_out.blast"
								wait
								if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
								wait
							fi
							wait
							gzip "${sub}_out.blast" &&
							wait
						fi
						rm "$sub" )&
						if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
							wait
						fi
					done
					wait
					for subfile in *_out.blast.gz; do
						cat "$subfile" >> ${ccf}.blast.gz &&
						rm "$subfile" && wait
					done
					wait
					# zcat ${ccf}.blast.gz 2> /dev/null | awk -v percid=$percid '$3 >= $5*(percid/100) {print $0}' | $gzip >> combined_compressed.megablast.gz &&
					cat "${ccf}".blast.gz >> combined_compressed.megablast.gz &&
					rm "${ccf}".blast.gz &&
					wait
				fi
				cd ../haplotig/splitccf/
				wait

			fi
			if [[ $nodes -gt 1 ]]; then
				count_megablast_node=$(ls ${projdir}/megablast_done_node*.txt 2> /dev/null | wc -l)
				while [[ "$count_megablast_node" < $(($nodes - 1)) ]]; do
					sleep 300
				done
				wait
				cd "${projdir}"/metagenome/alignment
				cat combined_compressed_node*.megablast.gz > combined_compressed.megablast.gz
				rm combined_compressed_node*.megablast.gz
				wait
			fi
			cd "${projdir}"/metagenome/haplotig/splitccf/
			rmdir * 2> /dev/null
			rm "${projdir}"/megablast_node* ${projdir}/multi_node_run_ready.txt ${projdir}/megablast_splitrun_node_*.sh 2> /dev/null
			cd ../
			rmdir splitccf 2> /dev/null
		else
			echo -e "${YELLOW}- Primary BLAST ouput already exist"
			echo -e "${YELLOW}- Skipping BLAST and filtering hits based on defined parameters"
		fi

		wait
		cd "${projdir}"/metagenome/alignment/
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA ==  false ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -i $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -vi $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA == true ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -i $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			mkdir rRNA
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -i 'rRNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/rRNA/rRNA_combined_compressed.megablast.gz
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | awk '{gsub(/77133;/,""); gsub(/340016;/,""); gsub(/175245;/,""); gsub(/137771;/,""); }1' | grep -vi $'uncultured\|unculture\|\t77133\t\|\t340016\t\|\t175245\t\|\t137771\t' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi
		wait


		cd "${projdir}"/metagenome/haplotig
		for i in $(ls *metagenome.fasta.gz); do (
		  if test ! -f ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz && test ! -f ../alignment/cultured/${i%_metagenome.fasta.gz}_haplotig.megablast.gz; then
				awk '!/^$/' <(zcat "$i" 2> /dev/null) | awk -F'\t' 'ORS=NR%2?"\t":"\n"' | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz &&
				wait
				awk -F'\t' 'ORS=NR%2?"\t":"\n"' <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | \
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' - <(zcat ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $3,$1}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz &&
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&

				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/uncultured_combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/uncultured_${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&
				wait
				rm ../alignment/"${i%_metagenome.fasta.gz}"_step*.txt.gz
			fi )&
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait

	fi
	wait
	find ../alignment/ -size 0 -delete
else
	cd "${projdir}"/metagenome/haplotig
	if test ! -f combined_compressed_metagenomes.fasta.gz; then
		# find . -name "*.fasta.gz" | xargs zcat | grep -v '^>' | awk '{A[$1]++}END{for(i in A)print i}' | awk '{print ">"NR"\n"$1}' > combined_compressed_metagenomes.fasta
		for i in $(ls *_compressed.fasta.gz); do (
			zcat $i | awk '{gsub(/AAAAAAAAAA$/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");gsub(/AAAAAAAAAA~/,"A~");}1' | \
			awk '{gsub(/CCCCCCCCCC$/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");gsub(/CCCCCCCCCC~/,"C~");}1' | \
			awk '{gsub(/GGGGGGGGGG$/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");gsub(/GGGGGGGGGG~/,"G~");}1' | \
			awk '{gsub(/TTTTTTTTTT$/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");gsub(/TTTTTTTTTT~/,"T~");}1' | \
			awk '{gsub(/^AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");gsub(/~AAAAAAAAAA/,"~A");}1' | \
			awk '{gsub(/^CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");gsub(/~CCCCCCCCCC/,"~C");}1' | \
			awk '{gsub(/^GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");gsub(/~GGGGGGGGGG/,"~G");}1' | \
			awk '{gsub(/^TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");gsub(/~TTTTTTTTTT/,"~T");}1' | \
			awk '{gsub(/AAAAA~/,"~");gsub(/CCCCC~/,"~");gsub(/GGGGG~/,"~");gsub(/TTTTT~/,"~");gsub(/~AAAAA/,"~");gsub(/~CCCCC/,"~");gsub(/~GGGGG/,"~");gsub(/~TTTTT/,"~");}1' | \
			awk '{gsub(/AA~/,"~");gsub(/CC~/,"~");gsub(/GG~/,"~");gsub(/TT~/,"~");gsub(/~AA/,"~");gsub(/~CC/,"~");gsub(/~GG/,"~");gsub(/~TT/,"~");gsub(/AA~/,"~");gsub(/CC~/,"~");gsub(/GG~/,"~");gsub(/TT~/,"~");gsub(/~AA/,"~");gsub(/~CC/,"~");gsub(/~GG/,"~");gsub(/~TT/,"~");}1' | \
			awk '{gsub(/A~/,"~");gsub(/C~/,"~");gsub(/G~/,"~");gsub(/T~/,"~");gsub(/~A/,"~");gsub(/~C/,"~");gsub(/~G/,"~");gsub(/~T/,"~");gsub(/~/,"");}1' | gzip > ${i}.tmp &&
			mv ${i}.tmp $i ) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		for i in $(ls *_compressed.fasta.gz); do
			zcat $i | grep -v '^>' | awk '{A[$1]++}END{for(i in A)print i}' | gzip >> combined_compressed_metagenomes.tmp.gz
			wait
		done
		wait
		awk '{A[$1]++}END{for(i in A)print i}' <(zcat combined_compressed_metagenomes.tmp.gz) | gzip > combined_compressed_metagenomes_full.tmp.gz
		rm combined_compressed_metagenomes.tmp.gz
		LC_ALL=C sort -T "${projdir}"/tmp <(zcat combined_compressed_metagenomes_full.tmp.gz) | gzip > combined_compressed_metagenomes_tmp.fasta.gz
		wait
		rm combined_compressed_metagenomes_full.tmp.gz
		awk '{print ">"NR"\n"$1}' <(zcat combined_compressed_metagenomes_tmp.fasta.gz) | $gzip > combined_compressed_metagenomes.fasta.gz
		wait
		rm combined_compressed_metagenomes_tmp.fasta.gz
	fi

	if [[ "$taxids" == true ]]; then
		for i in ${projdir}/taxids/*.txids; do
			cat "$i" >> ${projdir}/metagenome/All.txids
		done
		wait
	fi

	if [[ "$blast_location" =~ "local" ]]; then
		echo -e "${YELLOW}- performing local BLAST"
		if [[ -z "$(ls -R ${projdir}/metagenome/alignment/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]  && [[ -z "$(ls -R ${projdir}/metagenome/alignment/cultured/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]; then
  		if [[ "$taxids" == true ]]; then
  			${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
  			-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out combined_compressed.megablast
				wait
				if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				wait
  		else
  			${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) -db $local_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
  			-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out combined_compressed.megablast
				wait
				if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				wait
  		fi
  		wait
      $gzip combined_compressed.megablast
      wait
		else
			echo -e "${YELLOW}- Primary BLAST ouput already exist"
			echo -e "${YELLOW}- Skipping BLAST and filtering hits based on defined parameters"
		fi
		wait

		cd "${projdir}"/metagenome/alignment/
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA ==  false ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -vi 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA == true ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'uncultured\|unculture' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			mkdir rRNA
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'rRNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/rRNA/rRNA_combined_compressed.megablast.gz
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -vi 'uncultured\|unculture' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi

		wait

		cd "${projdir}"/metagenome/haplotig
		for i in $(ls *metagenome.fasta.gz); do (
		  if test ! -f ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz && test ! -f ../alignment/cultured/${i%_metagenome.fasta.gz}_haplotig.megablast.gz; then
				awk '!/^$/' <(zcat "$i" 2> /dev/null) | awk -F'\t' 'ORS=NR%2?"\t":"\n"' | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz &&
				wait
				awk -F'\t' 'ORS=NR%2?"\t":"\n"' <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | \
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' - <(zcat ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $3,$1}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz &&
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&

				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/uncultured_combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/uncultured_${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&
				wait
				rm ../alignment/"${i%_metagenome.fasta.gz}"_step*.txt.gz
			fi )&
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait
	fi

	if [[ "$blast_location" =~ "remote" ]]; then
		echo -e "${YELLOW}- performing a remote BLAST"
		if [[ -z "$(ls -R ${projdir}/metagenome/alignment/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]  && [[ -z "$(ls -R ${projdir}/metagenome/alignment/cultured/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]; then
			if [[ "$taxids" == true ]]; then
				${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) -db $blast_location -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
				-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" \
				-out ../alignment/combined_compressed.megablast -remote
				wait
				if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				wait
			else
				${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) -db $blast_location -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
				-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" \
				-out ../alignment/combined_compressed.megablast -remote
				wait
				if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				wait
			fi
		else
			echo -e "${YELLOW}- BLAST ouput already exist"
			echo -e "${YELLOW}- Skipping BLAST and filtering hits based on defined parameters"
		fi
		wait

		if test -f ../alignment/combined_compressed.megablast; then
			$gzip ../alignment/combined_compressed.megablast &&
			rm ../alignment/combined_compressed.megablast
		fi

		cd "${projdir}"/metagenome/alignment/
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA == false ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -vi 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
			zcat ../alignment/combined_compressed.megablast.gz > ../alignment/temp.megablast
			awk 'BEGIN{FS="\t";}{if(a[$1]<$3){a[$1]=$3;}}END{for(i in a){print i"\t"a[i];}}' temp.megablast | sort -T "${projdir}"/tmp -V -k1,1n | \
			awk -F'\t' 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1FS$2]++;next};c[$1FS$3] > 0' - ../alignment/temp.megablast  | $gzip > ../alignment/combined_compressed.megablast.gz
			zcat ../alignment/uncultured_combined_compressed.megablast.gz > ../alignment/uncultured_temp.megablast
			awk 'BEGIN{FS="\t";}{if(a[$1]<$3){a[$1]=$3;}}END{for(i in a){print i"\t"a[i];}}' uncultured_temp.megablast | sort -T "${projdir}"/tmp -V -k1,1n | \
			awk -F'\t' 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1FS$2]++;next};c[$1FS$3] > 0' - ../alignment/uncultured_temp.megablast  | $gzip > ../alignment/uncultured_combined_compressed.megablast.gz
		fi
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA == true ]] && [[ "$library_type" =~ "RRS" ]] || [[ "$library_type" =~ "rrs" ]] || [[ "$library_type" == "WGS" ]] || [[ "$library_type" == "wgs" ]] || [[ "$library_type" == "SHOTGUN" ]] || [[ "$library_type" == "shotgun" ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'uncultured\|unculture' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			mkdir rRNA
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'rRNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/rRNA/rRNA_combined_compressed.megablast.gz
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -vi 'uncultured\|unculture' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
			zcat ../alignment/combined_compressed.megablast.gz > ../alignment/temp.megablast
			awk 'BEGIN{FS="\t";}{if(a[$1]<$3){a[$1]=$3;}}END{for(i in a){print i"\t"a[i];}}' temp.megablast | sort -T "${projdir}"/tmp -V -k1,1n | \
			awk -F'\t' 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1FS$2]++;next};c[$1FS$3] > 0' - ../alignment/temp.megablast  | $gzip > ../alignment/combined_compressed.megablast.gz
			zcat ../alignment/uncultured_combined_compressed.megablast.gz > ../alignment/uncultured_temp.megablast
			awk 'BEGIN{FS="\t";}{if(a[$1]<$3){a[$1]=$3;}}END{for(i in a){print i"\t"a[i];}}' uncultured_temp.megablast | sort -T "${projdir}"/tmp -V -k1,1n | \
			awk -F'\t' 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1FS$2]++;next};c[$1FS$3] > 0' - ../alignment/uncultured_temp.megablast  | $gzip > ../alignment/uncultured_combined_compressed.megablast.gz
		fi
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA == true ]] && [[ "$library_type" =~ "amplicon" ]] || [[ "$library_type" =~ "Amplicon" ]] || [[ "$library_type" =~ "AMPLICON" ]] || [[ "$library_type" =~ "16S" ]] || [[ "$library_type" =~ "16s" ]]|| [[ "$library_type" =~ "ITS" ]] || [[ "$library_type" =~ "its" ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -vi 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
			zcat ../alignment/combined_compressed.megablast.gz > ../alignment/temp.megablast
			awk 'BEGIN{FS="\t";}{if(a[$1]<$3){a[$1]=$3;}}END{for(i in a){print i"\t"a[i];}}' temp.megablast | sort -T "${projdir}"/tmp -V -k1,1n | \
			awk -F'\t' 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1FS$2]++;next};c[$1FS$3] > 0' - ../alignment/temp.megablast  | $gzip > ../alignment/combined_compressed.megablast.gz
			zcat ../alignment/uncultured_combined_compressed.megablast.gz > ../alignment/uncultured_temp.megablast
			awk 'BEGIN{FS="\t";}{if(a[$1]<$3){a[$1]=$3;}}END{for(i in a){print i"\t"a[i];}}' uncultured_temp.megablast | sort -T "${projdir}"/tmp -V -k1,1n | \
			awk -F'\t' 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1FS$2]++;next};c[$1FS$3] > 0' - ../alignment/uncultured_temp.megablast  | $gzip > ../alignment/uncultured_combined_compressed.megablast.gz
		fi
		wait



		for i in $(ls *metagenome.fasta.gz); do (
			if test ! -f ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz && test ! -f ../alignment/cultured/${i%_metagenome.fasta.gz}_haplotig.megablast.gz; then
				awk '!/^$/' <(zcat "$i" 2> /dev/null) | awk -F'\t' 'ORS=NR%2?"\t":"\n"' | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz &&
				wait
				awk -F'\t' 'ORS=NR%2?"\t":"\n"' <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | \
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' - <(zcat ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $3,$1}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz &&
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&

				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/uncultured_combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/uncultured_${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&
				wait
				rm ../alignment/"${i%_metagenome.fasta.gz}"_step*.txt
			fi )&
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait
	fi

	if [[ "$blast_location" =~ "custom" ]]; then
		echo -e "${YELLOW}- performing custom BLAST"
		if [[ -z "$(ls -R ${projdir}/metagenome/alignment/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]  && [[ -z "$(ls -R ${projdir}/metagenome/alignment/cultured/ 2> /dev/null | grep combined_compressed.megablast.gz)" ]]; then
  		if [[ "$taxids" == true ]]; then
  			${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
  			-word_size 28 -taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out combined_compressed.megablast
				wait
				if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				wait
  		else
  			${Qmatey_dir}/tools/ncbi-blast-2.14.1+/bin/blastn -task megablast -query <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) -db $custom_db -num_threads 1 -perc_identity $percid -max_target_seqs $max_target -evalue 0.01 \
  			-word_size 28 -outfmt "6 qseqid sseqid length qstart qlen pident qseq sseq staxids stitle" -out combined_compressed.megablast
				wait
				if grep -qE 'Killed.*ncbi.*blastn.*megablast' ${projdir}/log.out; then printf "\nreduce parameter value for <reads_per_meagablast>, \nand then resubmit job to continue with megablast alignment\n" > ${projdir}/Megablast_killed_readme.txt; trap 'trap - SIGTERM && kill 0' SIGINT SIGTERM EXIT; fi
				wait
  		fi
  		wait
      $gzip combined_compressed.megablast
		else
			echo -e "${YELLOW}- Primary BLAST ouput already exist"
			echo -e "${YELLOW}- Skipping BLAST and filtering hits based on defined parameters"
		fi

		wait
		cd "${projdir}"/metagenome/alignment/
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA ==  false ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -vi 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA == true ]] && [[ "$library_type" =~ "RRS" ]] || [[ "$library_type" =~ "rrs" ]] || [[ "$library_type" == "WGS" ]] || [[ "$library_type" == "wgs" ]] || [[ "$library_type" == "SHOTGUN" ]] || [[ "$library_type" == "shotgun" ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'uncultured\|unculture' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
				mkdir rRNA
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'rRNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/rRNA/rRNA_combined_compressed.megablast.gz
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -vi 'uncultured\|unculture' | grep -vi 'rRNA\|ribosomal RNA\|ribosomal_RNA' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
	fi
		if test -f ${projdir}/metagenome/alignment/combined_compressed.megablast.gz && [[ $exclude_rRNA == true ]] && [[ "$library_type" =~ "amplicon" ]] || [[ "$library_type" =~ "Amplicon" ]] || [[ "$library_type" =~ "AMPLICON" ]] || [[ "$library_type" =~ "16S" ]] || [[ "$library_type" =~ "16s" ]]|| [[ "$library_type" =~ "ITS" ]] || [[ "$library_type" =~ "its" ]]; then
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -i 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/uncultured_combined_compressed.megablast.gz &&
			zcat ${projdir}/metagenome/alignment/combined_compressed.megablast.gz | grep -vi 'uncultured\|unculture' | $gzip > ${projdir}/metagenome/alignment/tmp_compressed.megablast.gz &&
			mv "${projdir}"/metagenome/alignment/tmp_compressed.megablast.gz ${projdir}/metagenome/alignment/combined_compressed.megablast.gz
		fi
		wait

		cd "${projdir}"/metagenome/haplotig
		for i in $(ls *metagenome.fasta.gz); do (
		  if test ! -f ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz && test ! -f ../alignment/cultured/${i%_metagenome.fasta.gz}_haplotig.megablast.gz; then
				awk '!/^$/' <(zcat "$i" 2> /dev/null) | awk -F'\t' 'ORS=NR%2?"\t":"\n"' | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz &&
				wait
				awk -F'\t' 'ORS=NR%2?"\t":"\n"' <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) | awk '{gsub(/-0/,""); gsub(/>/,"");}1' | \
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' - <(zcat ../alignment/${i%_metagenome.fasta.gz}_step1.txt.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $3,$1}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz &&
				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&

				awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}_step2.txt.gz) <( zcat ../alignment/uncultured_combined_compressed.megablast.gz 2> /dev/null) | \
				awk -F'\t' 'BEGIN{OFS="\t"}{print $12,$2,$3,$4,$5,$6,$7,$8,$9,$10}' | awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | gzip > ../alignment/uncultured_${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&
				wait
				rm ../alignment/"${i%_metagenome.fasta.gz}"_step*.txt.gz
			fi )&
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				wait
			fi
		done
		wait

	fi
	wait
	find ../alignment/ -size 0 -delete
fi

cd "${projdir}"/metagenome/alignment
mkdir -p uncultured
mv uncultured* ./uncultured/ 2> /dev/null
mkdir -p cultured
ls *megablast.gz 2> /dev/null | grep -v 'uncultured' | xargs mv -t ./cultured/ 2> /dev/null

printf "MegaBLAST was completed" > ${projdir}/metagenome/MegaBLAST_completed.txt

}
cd "${projdir}"
if test -f ${projdir}/metagenome/MegaBLAST_completed.txt; then
	echo -e "${YELLOW}- Qmatey has already performed NCBI MegaBLAST ${WHITE}"
else
	echo -e "${YELLOW}- Qmatey is performing sequence alignment using NCBI MegaBLAST ${WHITE}"
	time blast &>> ${projdir}/log.out
fi



#################################################################################################################
cd "${projdir}"/metagenome/alignment
if [[ "$(ls ./uncultured/*haplotig.megablast.gz | wc -l)" -lt 1 ]]; then
	for i in $(ls *haplotig.megablast.gz); do
		mv "$i" ./uncultured/uncultured_"${i}"
	done
fi
if [[ "$(ls ./cultured/*haplotig.megablast.gz | wc -l)" -lt 1 ]]; then
	for i in $(ls *haplotig.megablast.gz); do
		mv "$i" ./cultured/"${i}"
	done
fi

#################################################################################################################
strain_level() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Strain-Level classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact matching for strain-level profiling"
cd "${projdir}"/metagenome/sighits
mkdir -p sighits_strain
cd "${projdir}"/metagenome/results
mkdir -p strain_level
cd "${projdir}"/metagenome/alignment
if [[ "$fullqlen_alignment" == "true" ]]; then
	export percid=95
else
	export alnqlen=32
fi
if find ../sighits/sighits_strain/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits of at strain-level already available for each sample"
else
	for i in $(ls *_haplotig.megablast.gz); do (
		if [[ ! -f "../sighits/sighits_strain/${i%_haplotig.megablast.gz}_sighits.txt.gz" ]]; then
			if [[ "$taxids" == true ]]; then
				for emg in ${projdir}/taxids/*.txids; do
					if [[ "$fullqlen_alignment" == "false" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1] {print $0}' <(zcat "${i}" | awk '$6==100') <(zcat ${i%} | awk '$6==100') | awk '$3 >= 32 {print $0}' | awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | gzip > ${i%.gz}strain.gz
						wait
					fi
					if [[ "$fullqlen_alignment" == "true" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1] {print $0}' <(zcat "${i}" | awk '$6==100') <(zcat ${i%} | awk '$6==100') | awk '$3 == $5 {print $0}' | awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | gzip > ${i%.gz}strain.gz
						wait
					fi
					awk 'gsub(" ","_",$0)' <(zcat ${i%.gz}strain.gz) | awk -F'\t' '{print $1"___"$9}' | sort -T "${projdir}"/tmp | uniq | awk 'BEGIN{OFS="\t"}{gsub(/___/,"\t");}1' | awk '{print $1}' | \
					awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -F ' ' '{print $2"\t"$1}' | awk '$2 == 1' | awk '{print $1}' > ${i%_haplotig.megablast.gz}_exactmatch.txt
					awk 'gsub(" ","_",$0)' <(zcat ${i%.gz}strain.gz) | awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a {print; delete a[$1]}' ${i%_haplotig.megablast.gz}_exactmatch.txt - | \
					awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | awk 'BEGIN{OFS="\t"}{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | \
					awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_strain/${i%_haplotig.megablast.gz}_sighits.txt.gz &&
					rm ${i%_haplotig.megablast.gz}_exactmatch.txt ${i%.gz}strain.gz
					wait
				done
				zcat ../sighits/sighits_strain/${i%_haplotig.megablast.gz}_sighits.txt.gz | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | gzip > ../sighits/sighits_strain/${i%_haplotig.megablast.gz}_sighits.tmp.gz &&
				mv ../sighits/sighits_strain/${i%_haplotig.megablast.gz}_sighits.tmp.gz ../sighits/sighits_strain/${i%_haplotig.megablast.gz}_sighits.txt.gz
				wait
			else
				if [[ "$fullqlen_alignment" == "false" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1] {print $0}' <(zcat "$i" | awk '$6==100') <(zcat "$i" | awk '$6==100') | awk '$3 >= 32 {print $0}' | $gzip > ${i%.gz}strain.gz
				fi
				if [[ "$fullqlen_alignment" == "true" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1] {print $0}' <(zcat "$i" | awk '$6==100') <(zcat "$i" | awk '$6==100') | awk '$3 == $5 {print $0}' | $gzip > ${i%.gz}strain.gz
				fi
				awk 'gsub(" ","_",$0)1' <(zcat ${i%.gz}strain.gz) | awk -F'\t' '{print $1"___"$9}' | sort -T "${projdir}"/tmp | uniq | awk 'BEGIN{OFS="\t"}{gsub(/___/,"\t");}1' | awk '{print $1}' | \
				awk '{!seen[$0]++}END{for (i in seen) print seen[i], i}' | awk -F ' ' '{print $2"\t"$1}' | awk '$2 == 1' | awk '{print $1}' > ${i%_haplotig.megablast.gz}_exactmatch.txt
				awk 'gsub(" ","_",$0)1' <(zcat ${i%.gz}strain.gz) | awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a {print; delete a[$1]}' ${i%_haplotig.megablast.gz}_exactmatch.txt - | \
				awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | awk 'BEGIN{OFS="\t"}{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | \
				cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_strain/${i%_haplotig.megablast.gz}_sighits.txt.gz &&
				rm "${i%_haplotig.megablast.gz}"_exactmatch.txt "${i%.gz}"strain.gz
				wait
			fi
		fi ) &
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
fi

echo -e "${YELLOW}- compiling taxonomic information"
cd "${projdir}"/metagenome/sighits/sighits_strain
find . -type f -name '*_sighits.txt.gz' -exec cat {} + > sighits.txt.gz
awk -F '\t' '{print $8";"}' <(zcat sighits.txt.gz) | awk -F ';' '{print $1}' | sort -T "${projdir}"/tmp -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt.gz
awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  ${projdir}/rankedlineage_edited.dmp taxids_sighits.txt | \
awk '{gsub(/ /,"_"); print }' | awk '$2!=$3' > rankedlineage_subhits.txt
rm taxids_sighits.txt
cd "${projdir}"/metagenome/sighits/sighits_strain/
awk '{print $1}' rankedlineage_subhits.txt > strain_taxa_mean_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > strain_taxa_unique_sequences_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > strain_taxa_quantification_accuracy_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > strain_taxa_rel_quantification_accuracy_temp.txt

strain_level=strain
for i in $(ls *_sighits.txt.gz);do
	gunzip $i
	Rscript "${Qmatey_dir}/scripts/stats_summary.R" ${i%.gz} $strain_level "${Qmatey_dir}/tools/R"
	wait
	$gzip ${i%.gz}
	echo $'tax_id\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt strain_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > strain_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt strain_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > strain_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt strain_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > strain_taxa_quantification_accuracy_temp.txt
	id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdrelstderr.txt strain_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt && cat holdrelstderr2.txt > strain_taxa_rel_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	rm *stats1* *stats2* *stats3* *hold*
done
wait

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' strain_taxa_mean_temp.txt > strain_taxa_mean_temp2.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_mean_temp2.txt > ../../results/strain_level/strain_taxa_mean.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' strain_taxa_unique_sequences_temp.txt > strain_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_unique_sequences_temp2.txt > ../../results/strain_level/strain_taxa_unique_sequences.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_quantification_accuracy_temp.txt > strain_taxa_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/strain_level/strain_taxa_mean.txt strain_taxa_quantification_accuracy_temp2.txt > ../../results/strain_level/strain_taxa_quantification_accuracy.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_rel_quantification_accuracy_temp.txt > strain_taxa_rel_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/strain_level/strain_taxa_mean.txt strain_taxa_rel_quantification_accuracy_temp2.txt > ../../results/strain_level/strain_taxa_rel_quantification_accuracy.txt &&
rm *_temp*


cd "${projdir}"/metagenome/results/strain_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_strain/rankedlineage_subhits.txt strain_taxa_"${i}".txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > strain_taxainfo_"${i}".txt &&
	wait
done
wait
rm *_taxa_*


awk '{print $1}' strain_taxainfo_mean.txt > strain_taxainfo_mean_holdingtaxid.txt
awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
strain_taxainfo_mean.txt | awk '!($1="")' | awk 'BEGIN{OFS="\t"} {$1=$1};1' > strain_taxainfo_mean_holdingtaxinfo.txt
touch strain_taxainfo_mean_buildnorm.txt
for i in ../../../metagenome/haplotig/*_metagenome.fasta.gz; do
	sample=${i%_metagenome.fasta.gz}; sample=${sample##*/}
	if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
		normfactor=1
	else
		normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
	fi
	awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" strain_taxainfo_mean.txt | \
	paste - strain_taxainfo_mean_buildnorm.txt > strain_taxainfo_mean_buildnorm0.txt &&
	mv strain_taxainfo_mean_buildnorm0.txt strain_taxainfo_mean_buildnorm.txt
done
wait
paste strain_taxainfo_mean_holdingtaxid.txt strain_taxainfo_mean_buildnorm.txt > strain_taxainfo_mean_norm0.txt
paste strain_taxainfo_mean_norm0.txt strain_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > strain_taxainfo_mean_normalized.txt
rm strain_taxainfo_mean_holdingtaxid.txt strain_taxainfo_mean_buildnorm.txt strain_taxainfo_mean_holdingtaxinfo.txt strain_taxainfo_mean_norm0.txt

for i in $(ls *.txt); do
	awk '{gsub(/-/,"_"); print}' $i > ${i%.txt}.temp &&
	:> $i &&
	mv ${i%.txt}.temp $i
done
wait

Rscript "${Qmatey_dir}/scripts/cleanup_rank.R" "strain" "${Qmatey_dir}/tools/R" 2>/dev/null
wait

cd "${projdir}"/metagenome/results/
cp -r strain_level strain_level_hold
for min_strain_uniq_ematch in ${min_strain_uniq//,/ }; do
	cd ./strain_level
	# if [[ "$genome_scaling" == true ]]; then
	# 	if [[ "$subsample_shotgun_R1" == false ]]; then
	# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" strain "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "true" &>/dev/null
	# 	else
	# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" strain "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "false" &>/dev/null
	# 	fi
	# else
	# 		Rscript "${Qmatey_dir}/scripts/error_zero_inflation.R" strain "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated &>/dev/null
	# fi

	file=${projdir}/exclude_taxa.txt
	# if test -f $file; then
	# 	cat strain_taxainfo_mean.txt > strain_taxainfo_mean_filtered.txt &&
	# 	cat strain_taxainfo_unique_sequences.txt > strain_taxainfo_unique_sequences_filtered.txt &&
	# 	cat strain_taxainfo_quantification_accuracy.txt > strain_taxainfo_quantification_accuracy_filtered.txt &&
	# 	cat strain_taxainfo_rel_quantification_accuracy.txt > strain_taxainfo_rel_quantification_accuracy_filtered.txt &&
	# 	wait
	# 	while IFS="" read -r line || [ -n "$line" ]; do
	# 		while IFS="" read -r i || [ -n "$i" ]; do
	# 			awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt &&
	# 			mv ${i%.txt}_temp.txt $i &&
	# 			wait
	# 		done <(ls *filtered.txt)
	# 	done < $file
	# 	wait
	# fi


	if test -f $file; then
		echo -e "${YELLOW}- creating strain-level visualizations"
		cd "${projdir}"/metagenome/results/strain_level
		strain_level_mean=strain_taxainfo_mean_filtered.txt
		strain_level_uniq=strain_taxainfo_unique_sequences_filtered.txt
		strain_level_stderr=strain_taxainfo_quantification_accuracy_filtered.txt
		strain_level_rel_stderr=strain_taxainfo_rel_quantification_accuracy_filtered.txt

		if [[ -z $min_percent_sample ]]; then
			min_percent_sample=5,10,20
		fi


		for min_perc in ${min_percent_sample//,/ }; do (
			Rscript "${Qmatey_dir}/scripts/strain_level_boxplots.R" "$strain_level_mean" "$strain_level_uniq" "$strain_level_stderr" "$strain_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
			)&
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			 wait
			fi
		done
		wait
		mkdir -p boxplots
		mv *_files ${projdir}/metagenome/results/strain_level/boxplots/ 2> /dev/null
		mv *.html ${projdir}/metagenome/results/strain_level/boxplots/ 2> /dev/null

	else
		echo -e "${YELLOW}- creating strain-level visualizations"
		cd "${projdir}"/metagenome/results/strain_level
		strain_level_mean=strain_taxainfo_mean.txt
		strain_level_mean_norm=strain_taxainfo_mean_normalized.txt
		strain_level_uniq=strain_taxainfo_unique_sequences.txt
		strain_level_stderr=strain_taxainfo_quantification_accuracy.txt
		strain_level_rel_stderr=strain_taxainfo_rel_quantification_accuracy.txt

		if [[ -z $min_percent_sample ]]; then
			min_percent_sample=5,10,20
		fi

		for min_perc in ${min_percent_sample//,/ }; do (
			Rscript "${Qmatey_dir}/scripts/strain_level_boxplots.R" "$strain_level_mean" "$strain_level_mean_norm" "$strain_level_uniq" "$strain_level_stderr" "$strain_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
			)&
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			 wait
			fi
		done
		wait
		mkdir -p boxplots
		mv *_files ${projdir}/metagenome/results/strain_level/boxplots 2> /dev/null
		mv *.html ${projdir}/metagenome/results/strain_level/boxplots 2> /dev/null

	fi
	cd ../
	mv strain_level strain_level_minUniq_${min_strain_uniq_ematch}
	cp -r strain_level_hold strain_level
done
wait
rm -rf strain_level_hold strain_level

}
if [[ "$strain_level" == "true" ]]; then
	if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_strain_level*/strain_taxainfo* 2> /dev/null)" ]]; then
		cd "${projdir}"/metagenome/alignment
		mv ./uncultured/uncultured*megablast.gz ./ &&
		mv ./uncultured_combined_compressed.megablast.gz ./uncultured/ &&
		for i in $(ls *megablast.gz); do mv "$i" ${i#uncultured_}; done
		wait
		cd "${projdir}"
		time strain_level 2>> ${projdir}/log.out
		wait
		cd "${projdir}"/metagenome/alignment/
		for i in $(ls *megablast.gz); do mv "$i" ./uncultured/uncultured_"${i}"; done
		wait
		cd "${projdir}"/metagenome/results/ &&
		for dirc in strain_level*; do mv $dirc uncultured_${dirc}; done
		wait
		mv ../sighits/sighits_strain ../sighits/uncultured_sighits_strain
		wait
	fi
	if [[ -z "$(ls -A ${projdir}/metagenome/results/strain_level*/strain_taxainfo* 2> /dev/null)" ]]; then
		cd "${projdir}"/metagenome/alignment/cultured/
		mv *megablast* ../ &&
		mv ../combined_compressed.megablast.gz ./ &&
		cd "${projdir}"
		time strain_level 2>> ${projdir}/log.out
		wait
		mv "${projdir}"/metagenome/alignment/*megablast* ${projdir}/metagenome/alignment/cultured/
	fi
else
	echo -e "${magenta}- skipping exact matching for strain-level profiling ${white}\n"
fi

#################################################################################################################
species_level() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Species-Level classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact matching for species-level profiling"
cd "${projdir}"/metagenome/sighits
mkdir -p sighits_species
cd "${projdir}"/metagenome/results
mkdir -p species_level
cd "${projdir}"/metagenome/alignment
if find ../sighits/sighits_species/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits of at species-level already available for each sample"
else
	for i in $(ls *_haplotig.megablast.gz);do (
		if [[ ! -f "../sighits/sighits_species/${i%_haplotig.megablast.gz}_sighits.txt.gz" ]]; then
			if [[ "$taxids" == true ]]; then
			  for emg in ${projdir}/taxids/*.txids; do
					if [[ "$fullqlen_alignment" == "false" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-2 {print $0}' <(zcat "$i" | awk '$6>=99') <(zcat "$i" | awk '$6>=99') | awk '$3 >= 32 {print $0}' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_species/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
					if [[ "$fullqlen_alignment" == "true" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-2 {print $0}' <(zcat "$i" | awk '$6>=99') <(zcat "$i" | awk '$6>=99') | awk '$3 >= ($5-($5*0.01)) {print $0}' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_species/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
			  done
				zcat ../sighits/sighits_species/${i%_haplotig.megablast.gz}_sighits.txt.gz | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | gzip > ../sighits/sighits_species/${i%_haplotig.megablast.gz}_sighits.tmp.gz &&
				mv ../sighits/sighits_species/${i%_haplotig.megablast.gz}_sighits.tmp.gz ../sighits/sighits_species/${i%_haplotig.megablast.gz}_sighits.txt.gz
				wait
			else
				if [[ "$fullqlen_alignment" == "false" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-2 {print $0}' <(zcat "$i" | awk '$6>=99') <(zcat "$i" | awk '$6>=99') | awk '$3 >= 32 {print $0}' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_species/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
				if [[ "$fullqlen_alignment" == "true" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-2 {print $0}' <(zcat "$i" | awk '$6>=99') <(zcat "$i" | awk '$6>=99') | awk '$3 >= ($5-($5*0.01)) {print $0}' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_species/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
			fi
		fi ) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_species
	for i in $(ls *_sighits.txt.gz);do (
		zcat "$i" | awk -F '\t' 'NR>1{a[$10]++;b[$10]=$0}END{for(x in a)if(a[x]==1)print b[x]}' > ${i%_sighits*}_species_unique_reads.txt
		$gzip ${i%_sighits*}_species_unique_reads.txt ) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

  echo -e "${YELLOW}- compiling species-level multi-alignment algorithm"

  cd "${projdir}"/metagenome/sighits/sighits_species
  		for i in $(ls *_sighits.txt.gz);do (
				gunzip $i
				gunzip ${i%_sighits*}_species_unique_reads.txt.gz
				awk -F '\t' 'FNR==NR{a[$10]=1; next} {print $0}' ${i%_sighits*}_species_unique_reads.txt <(awk 'NR>1{print $0}' ${i%.gz}) OFS='\t' > ${i%_sighits*}_dup.txt
				rm "${i%.gz}"
				$gzip ${i%_sighits*}_species_unique_reads.txt
				$gzip ${i%_sighits*}_dup.txt
				) &
				if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
					wait
				fi
  		done
  		wait


  cd "${projdir}"/metagenome/sighits/sighits_species
  for i in $(ls *_dup.txt.gz);do (
  	awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_dup*}_taxids_dup.txt
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  cd "${projdir}"/metagenome/sighits/sighits_species

  for i in $(ls *_taxids_dup.txt); do
  	awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		cut -f 2-11 -d$'\t' > ${i%_taxids_dup*}_dup_inter.txt
  done
  wait

	for i in $(ls *_dup_inter.txt); do (
	   awk -F '\t'  '{print $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait
	rm *_taxids_dup.txt



  for i in $(ls *_dup.txt.gz);do
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <(zcat "$i") ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8}' OFS='\t' ${i%*_dup.txt.gz}_species_taxa.txt) | \
		awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18}' OFS='\t' > ${i%_dup*}_species_duplicates.txt
		$gzip ${i%_dup*}_species_duplicates.txt
  done
	wait

	for i in $(ls *_species_duplicates.txt.gz);do (
		# zcat "$i" | grep -v $'^\([^\t]*\t\)\{11\}\"NA"\t' > ${i%*_species_duplicates.txt.gz}temp.txt
		awk -F'\t' '$11 != "NA"' OFS='\t' <( zcat "$i") > ${i%*_species_duplicates.txt.gz}temp.txt
		$gzip ${i%*_species_duplicates.txt.gz}temp.txt
		rm "$i" &&
		mv ${i%*_species_duplicates.txt.gz}temp.txt.gz $i
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
	wait
  rm *_dup_inter.txt *_dup.txt.gz *_species_taxa.txt

  for i in $(ls *_species_duplicates.txt.gz);do (
    awk -F '\t' '{print $1, $10"~"$11, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, $14, $15, $16, $17, $18}' OFS='\t' <(zcat "$i") > ${i%_species_duplicates*}_species_inter.txt ) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

  for i in $(ls *_species_inter.txt);do (
    awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}' OFS='\t' $i > ${i%_species_inter*}_species_inter2.txt ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
  wait

  for i in $(ls *_species_inter2.txt);do (
    awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -T "${projdir}"/tmp -k1,1  > ${i%_species_inter2*}_duplicate_count.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_duplicate_count.txt);do (
    awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -T "${projdir}"/tmp -k1,1 > ${i%_duplicate_count*}_multialign_species_reads.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_multialign_species_reads.txt);do (
    awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_species_reads*}_species_inter2.txt | sort -T "${projdir}"/tmp -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $4, $5, $6, $7, $8, $9, $10, $11, $1, $2, $12, $13, $14, $15, $16, $17, $18}' OFS='\t' > ${i%_multialign_species_reads*}_species_OTU.txt
		) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait


  for i in $(ls *_species_unique_reads.txt.gz);do (
    awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_species_unique_reads*}_taxids_uniq.txt ) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_species
	for i in $(ls *_taxids_uniq.txt); do
	   awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		 cut -f 2-11 -d$'\t' > ${i%_taxids_uniq*}_uniq_inter.txt
	done
	wait

	for i in $(ls *_uniq_inter.txt);do (
	   awk -F '\t'  '{print $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxa.txt ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_taxids_uniq.txt


	for i in $(ls *_unique_reads.txt.gz);do (
	   paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <(zcat "$i" )) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8}' OFS='\t' "${i%*_species_unique_reads*}"_species_taxa.txt) | \
		 awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18}' OFS='\t' > "${i%_species_uniq*}"_species_unique_uncultured.txt &&
		 wait ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait
	for i in $(ls *_species_unique_uncultured.txt);do (
		awk -F '\t' '$11 != "NA"' OFS='\t' "$i" > "${i%*_species_unique_uncultured.txt}"temp.txt &&
		mv "${i%*_species_unique_uncultured.txt}"temp.txt "$i" 2> /dev/null
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
	wait
	rm *_uniq_inter.txt && rm *_species_taxa.txt 2> /dev/null

	for i in $(ls *_species_unique_uncultured.txt);do (
	   mv "$i" "${i%*_species_unique_uncultured*}"_unique_sequences.txt 2> /dev/null ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_species_unique_uncultured.txt *_species_inter.txt *_species_inter2.txt *_duplicate_count.txt *_multialign_species_reads.txt *_species_duplicates.txt.gz 2> /dev/null

	for i in $(ls *_species_OTU.txt);do (
	   cat "$i" ${i%_species_OTU*}_unique_sequences.txt > ${i%_species_OTU*}_complete_species_reads.txt ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	echo -e "${YELLOW}- compiling taxonomic information"
	for i in $(ls *_complete_species_reads.txt);do (
	   awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_species_reads*}_sighits_temp.txt
	    echo $'abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	     cat - ${i%_complete_species_reads*}_sighits_temp.txt > ${i%_complete_species_reads*}_sighits_temp2.txt
			 ) &
			 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
				 wait
			 fi
	done
	wait

	for i in $(ls *_sighits_temp2.txt);do (
	   awk -F '\t' '{gsub(/ /,"_");print}' $i > ${i%_sighits_temp2*}_sighits.txt
     $gzip ${i%_sighits_temp2*}_sighits.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_complete_species_reads.txt *_sighits_temp.txt *_unique_reads.txt.gz *_unique_sequences.txt *_sighits_temp2.txt *_species_OTU.txt
fi

cd "${projdir}"/metagenome/sighits/sighits_species
find . -type f -name '*_sighits.txt.gz' -exec cat {} + > sighits.txt.gz
awk '{print $8}' <(zcat sighits.txt.gz) | awk '{gsub(";","\n"); print}' | sort -T "${projdir}"/tmp -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt.gz
awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  ${projdir}/rankedlineage_edited.dmp taxids_sighits.txt | \
awk '{gsub(/ /,"_"); print }' | awk '{print $3, $4, $5, $6, $7, $8, $9, $10}' | awk 'NR == 1; NR > 1 {print $0}' | awk '{gsub(/ /,"\t");}1' > rankedlineage_subhits.txt
rm taxids_sighits.txt

cd "${projdir}"/metagenome/sighits/sighits_species/
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/species/' > species_taxa_mean_temp1.txt
echo -e 'species' | cat - species_taxa_mean_temp1.txt > species_taxa_mean_temp.txt && rm species_taxa_mean_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/species/'  > species_taxa_unique_sequences_temp1.txt
echo -e 'species' | cat - species_taxa_unique_sequences_temp1.txt > species_taxa_unique_sequences_temp.txt && rm species_taxa_unique_sequences_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt |  sort -T "${projdir}"/tmp -u | awk '!/species/'  > species_taxa_quantification_accuracy_temp1.txt
echo -e 'species' | cat - species_taxa_quantification_accuracy_temp1.txt > species_taxa_quantification_accuracy_temp.txt && rm species_taxa_quantification_accuracy_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/species/'  > species_taxa_rel_quantification_accuracy_temp1.txt
echo -e 'species' | cat - species_taxa_rel_quantification_accuracy_temp1.txt > species_taxa_rel_quantification_accuracy_temp.txt && rm species_taxa_rel_quantification_accuracy_temp1.txt

species_level=species
for i in $(ls *_sighits.txt.gz);do
	gunzip $i
 	Rscript ${Qmatey_dir}/scripts/stats_summary.R ${i%.gz} species "${Qmatey_dir}/tools/R"
	wait
	$gzip ${i%.gz}
  echo $'species\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
  id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdmean.txt species_taxa_mean_temp.txt > holdmean2.txt
  cat holdmean2.txt > species_taxa_mean_temp.txt
  id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holduniq_reads.txt species_taxa_unique_sequences_temp.txt > holduniq_reads2.txt
  cat holduniq_reads2.txt > species_taxa_unique_sequences_temp.txt
  id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt species_taxa_quantification_accuracy_temp.txt > holdstderr2.txt
  cat holdstderr2.txt > species_taxa_quantification_accuracy_temp.txt
  id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdrelstderr.txt species_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt
  cat holdrelstderr2.txt > species_taxa_rel_quantification_accuracy_temp.txt
	wait
  rm *stats1* *stats2* *hold*
done
wait

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' species_taxa_mean_temp.txt > species_taxa_mean_temp2.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_mean_temp2.txt > ../../results/species_level/species_taxa_mean.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' species_taxa_unique_sequences_temp.txt > species_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_unique_sequences_temp2.txt > ../../results/species_level/species_taxa_unique_sequences.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_quantification_accuracy_temp.txt > species_taxa_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/species_level/species_taxa_mean.txt species_taxa_quantification_accuracy_temp2.txt > ../../results/species_level/species_taxa_quantification_accuracy.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_rel_quantification_accuracy_temp.txt > species_taxa_rel_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/species_level/species_taxa_mean.txt species_taxa_rel_quantification_accuracy_temp2.txt > ../../results/species_level/species_taxa_rel_quantification_accuracy.txt
rm *_temp*

cd "${projdir}"/metagenome/results/species_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do
  awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_species/rankedlineage_subhits.txt species_taxa_"${i}".txt | \
  awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > species_taxainfo_"${i}".txt &&
  wait
done
wait
rm *_taxa_*

awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
species_taxainfo_mean.txt | awk 'BEGIN{OFS="\t"} {$1=$1};1' > species_taxainfo_mean_holdingtaxinfo.txt
touch species_taxainfo_mean_buildnorm.txt
for i in ../../../metagenome/haplotig/*_metagenome.fasta.gz; do
  sample=${i%_metagenome.fasta.gz}; sample=${sample##*/}
  if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
    normfactor=1
  else
    normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
  fi
awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" species_taxainfo_mean.txt | \
paste - species_taxainfo_mean_buildnorm.txt > species_taxainfo_mean_buildnorm0.txt &&
mv species_taxainfo_mean_buildnorm0.txt species_taxainfo_mean_buildnorm.txt
done
wait
paste species_taxainfo_mean_buildnorm.txt > species_taxainfo_mean_norm0.txt
paste species_taxainfo_mean_norm0.txt species_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > species_taxainfo_mean_normalized.txt
# pos=$(awk -v RS='\t' '/species/{print NR; exit}' species_taxainfo_mean_normalized.txt)
# awk -v pat="$pos" -F'\t'  'BEGIN {OFS="\t"}; {k=$pat; $pat=""; print k,$0}' species_taxainfo_mean_normalized.txt | awk 'BEGIN{FS="\t+"; OFS="\t"} {$1=$1; print}' > species_taxainfo_mean_normalized.tmp
# mv species_taxainfo_mean_normalized.tmp species_taxainfo_mean_normalized.txt
rm species_taxainfo_mean_buildnorm.txt species_taxainfo_mean_holdingtaxinfo.txt species_taxainfo_mean_norm0.txt 2> /dev/null

for i in $(ls *.txt); do
  awk '{gsub(/-/,"_"); print}' $i > ${i%.txt}.temp &&
  mv ${i%.txt}.temp $i
done
wait

Rscript "${Qmatey_dir}/scripts/cleanup_rank.R" "species" "${Qmatey_dir}/tools/R" 2>/dev/null
wait


# if [[ "$genome_scaling" == true ]]; then
# 	if [[ "$subsample_shotgun_R1" == false ]]; then
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" species "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "true" &>/dev/null
# 	else
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" species "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "false" &>/dev/null
# 	fi
# else
# 		Rscript "${Qmatey_dir}/scripts/error_zero_inflation.R" species "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated &>/dev/null
# fi

file=${projdir}/exclude_taxa.txt
# if test -f $file; then
#   cat species_taxainfo_mean.txt > species_taxainfo_mean_filtered.txt &&
#   cat species_taxainfo_unique_sequences.txt > species_taxainfo_unique_sequences_filtered.txt &&
#   cat species_taxainfo_quantification_accuracy.txt > species_taxainfo_quantification_accuracy_filtered.txt &&
#   cat species_taxainfo_rel_quantification_accuracy.txt > species_taxainfo_rel_quantification_accuracy_filtered.txt &&
#   while read -r line; do
#     for i in $( ls *filtered.txt ); do
#       awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i
#     done
# 		wait
#   done < $file
# 	wait
# fi

if test -f $file; then
	echo -e "${YELLOW}- creating species-level visualizations"
	cd "${projdir}"/metagenome/results/species_level
	species_level_mean=species_taxainfo_mean_filtered.txt
	species_level_uniq=species_taxainfo_unique_sequences_filtered.txt
	species_level_stderr=species_taxainfo_quantification_accuracy_filtered.txt
	species_level_rel_stderr=species_taxainfo_rel_quantification_accuracy_filtered.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi


	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/species_level_boxplots.R" "$species_level_mean" "$species_level_uniq" "$species_level_stderr" "$species_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/species_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/species_level/boxplots/ 2> /dev/null
else
	echo -e "${YELLOW}- creating species-level visualizations"
	cd "${projdir}"/metagenome/results/species_level
	species_level_mean=species_taxainfo_mean.txt
	species_level_mean_norm=species_taxainfo_mean_normalized.txt
	species_level_uniq=species_taxainfo_unique_sequences.txt
	species_level_stderr=species_taxainfo_quantification_accuracy.txt
	species_level_rel_stderr=species_taxainfo_rel_quantification_accuracy.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/species_level_boxplots.R" "$species_level_mean" "$species_level_mean_norm" "$species_level_uniq" "$species_level_stderr" "$species_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/species_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/species_level/boxplots/ 2> /dev/null
fi

}
if [[ "$species_level" == "true" ]]; then
	if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_species_level/species_taxainfo* 2> /dev/null)" ]]; then
		cd "${projdir}"/metagenome/alignment/
		mv ./uncultured/uncultured*megablast.gz ./ &&
		mv ./uncultured_combined_compressed.megablast.gz ./uncultured/ &&
		for i in $(ls *megablast.gz); do mv "$i" ${i#uncultured_}; done
		wait
		cd "${projdir}"
		time species_level 2>> ${projdir}/log.out
		wait
		cd "${projdir}"/metagenome/alignment/
		for i in $(ls *megablast.gz); do mv "$i" ./uncultured/uncultured_"${i}"; done
		wait
		mv "${projdir}"/metagenome/results/species_level ${projdir}/metagenome/results/uncultured_species_level
		mv "${projdir}"/metagenome/sighits/sighits_species ${projdir}/metagenome/sighits/uncultured_sighits_species
	fi
	if [[ -z "$(ls -A ${projdir}/metagenome/results/species_level/species_taxainfo* 2> /dev/null)" ]]; then
		cd "${projdir}"/metagenome/alignment/cultured/
		mv *megablast* ../
		mv ../combined_compressed.megablast.gz ./
		cd "${projdir}"
		time species_level 2>> ${projdir}/log.out
		mv "${projdir}"/metagenome/alignment/*megablast* ${projdir}/metagenome/alignment/cultured/
	fi
else
	echo -e "${magenta}- skipping exact matching for species-level profiling ${white}\n"
fi

###########################################################################################
genus_level() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Genus-Level classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact matching for genus-level profiling"
cd "${projdir}"/metagenome/sighits
mkdir -p sighits_genus
cd "${projdir}"/metagenome/results
mkdir -p genus_level
cd "${projdir}"/metagenome/alignment
if find ../sighits/sighits_genus/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits of at genus-level already available for each sample"
else
	for i in $(ls *_haplotig.megablast.gz);do (
		if [[ ! -f "../sighits/sighits_genus/${i%_haplotig.megablast.gz}_sighits.txt.gz" ]]; then
			if [[ "$taxids" == true ]]; then
			  for emg in ${projdir}/taxids/*.txids; do
					if [[ "$fullqlen_alignment" == "false" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-4 {print $0}' <(zcat "${i}" | awk '$6>=98') <(zcat "${i}" | awk '$6>=98') | awk '$3 >= 32 {print $0}' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_genus/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
					if [[ "$fullqlen_alignment" == "true" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-4 {print $0}' <(zcat "${i}" | awk '$6>=98') <(zcat "${i}" | awk '$6>=98') | awk '$3 >= ($5-($5*0.02))' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_genus/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
			  done
				zcat ../sighits/sighits_genus/${i%_haplotig.megablast.gz}_sighits.txt.gz | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | gzip > ../sighits/sighits_genus/${i%_haplotig.megablast.gz}_sighits.tmp.gz &&
				mv ../sighits/sighits_genus/${i%_haplotig.megablast.gz}_sighits.tmp.gz ../sighits/sighits_genus/${i%_haplotig.megablast.gz}_sighits.txt.gz
				wait
			else
				if [[ "$fullqlen_alignment" == "false" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-4 {print $0}' <(zcat "$i" | awk '$6>=98') <(zcat "$i" | awk '$6>=98') | awk '$3 >= 32 {print $0}' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_genus/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
				if [[ "$fullqlen_alignment" == "true" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-4 {print $0}' <(zcat "$i" | awk '$6>=98') <(zcat "$i" | awk '$6>=98') | awk '$3 >= ($5-($5*0.02))' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_genus/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
			fi
		fi
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_genus
	for i in $(ls *_sighits.txt.gz);do (
		zcat "$i" | awk -F '\t' 'NR>1{a[$10]++;b[$10]=$0}END{for(x in a)if(a[x]==1)print b[x]}' > ${i%_sighits*}_genus_unique_reads.txt
		$gzip ${i%_sighits*}_genus_unique_reads.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

  echo -e "${YELLOW}- compiling genus-level multi-alignment algorithm"

  cd "${projdir}"/metagenome/sighits/sighits_genus
  		for i in $(ls *_sighits.txt.gz);do (
				gunzip $i
				gunzip ${i%_sighits*}_genus_unique_reads.txt.gz
				awk -F '\t' 'FNR==NR{a[$10]=1; next} {print $0}' ${i%_sighits*}_genus_unique_reads.txt <(awk 'NR>1{print $0}' ${i%.gz}) OFS='\t' > ${i%_sighits*}_dup.txt
				rm "${i%.gz}"
				$gzip ${i%_sighits*}_genus_unique_reads.txt
				$gzip ${i%_sighits*}_dup.txt
				) &
				if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
					wait
				fi
  		done
  		wait


  cd "${projdir}"/metagenome/sighits/sighits_genus
  for i in $(ls *_dup.txt.gz);do (
  	awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_dup*}_taxids_dup.txt

		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

  cd "${projdir}"/metagenome/sighits/sighits_genus

  for i in $(ls *_taxids_dup.txt); do
  	awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		cut -f 2-11 -d$'\t' > ${i%_taxids_dup*}_dup_inter.txt
  done
  wait



	for i in $(ls *_dup_inter.txt);do (
	   awk -F '\t'  '{print $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_genus_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait
	rm *_taxids_dup.txt


  for i in $(ls *_dup.txt.gz);do (
    paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <(zcat "$i") ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7}' OFS='\t' ${i%*_dup.txt.gz}_genus_taxa.txt)  | \
		awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}' OFS='\t' > ${i%_dup*}_genus_duplicates.txt
		$gzip ${i%_dup*}_genus_duplicates.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

	for i in $(ls *_genus_duplicates.txt.gz);do (
     awk -F'\t' '$11 != "NA"' OFS='\t' <( zcat "$i") > ${i%*_genus_duplicates.txt.gz}temp.txt
		 $gzip ${i%*_genus_duplicates.txt.gz}temp.txt
 		 mv ${i%*_genus_duplicates.txt.gz}temp.txt.gz $i
 		) &
 		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
 			wait
 		fi
   done
 	wait
  rm *_dup_inter.txt *_dup.txt.gz *_genus_taxa.txt

  for i in $(ls *_genus_duplicates.txt.gz);do (
    awk -F '\t' '{print $1, $10"~"$11, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, $14, $15, $16, $17}' OFS='\t' <(zcat "$i") > ${i%_genus_duplicates*}_genus_inter.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

  for i in $(ls *_genus_inter.txt);do (
    awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' OFS='\t' $i > ${i%_genus_inter*}_genus_inter2.txt ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
  wait

  for i in $(ls *_genus_inter2.txt);do (
    awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -T "${projdir}"/tmp -k1,1  > ${i%_genus_inter2*}_duplicate_count.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_duplicate_count.txt);do (
    awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -T "${projdir}"/tmp -k1,1 > ${i%_duplicate_count*}_multialign_genus_reads.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_multialign_genus_reads.txt);do (
    awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_genus_reads*}_genus_inter2.txt | sort -T "${projdir}"/tmp -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $4, $5, $6, $7, $8, $9, $10, $11, $1, $2, $12, $13, $14, $15, $16, $17}' OFS='\t' > ${i%_multialign_genus_reads*}_genus_OTU.txt
		) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
  wait

  for i in $(ls *_genus_unique_reads.txt.gz);do (
    awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_genus_unique_reads*}_taxids_uniq.txt ) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait


	cd "${projdir}"/metagenome/sighits/sighits_genus
	for i in $(ls *_taxids_uniq.txt);do
	   awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		 cut -f 2-11 -d$'\t' > ${i%_taxids_uniq*}_uniq_inter.txt
	done
	wait

	for i in $(ls *_uniq_inter.txt);do (
	   awk -F '\t'  '{print $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_genus_taxa.txt ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_taxids_uniq.txt


	for i in $(ls *_unique_reads.txt.gz);do (
	   paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <(zcat "$i" )) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7}' OFS='\t' "${i%*_genus_unique_reads*}"_genus_taxa.txt) | \
		 awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17}' OFS='\t' > "${i%_genus_uniq*}"_genus_unique_uncultured.txt &&
		 wait ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait
	for i in $(ls *_genus_unique_uncultured.txt);do (
		awk -F '\t' '$11 != "NA"' OFS='\t' "$i" > "${i%*_genus_unique_uncultured.txt}"temp.txt &&
		:> $i &&
		mv "${i%*_genus_unique_uncultured.txt}"temp.txt "$i" 2> /dev/null
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
	wait
	rm *_uniq_inter.txt && rm *_genus_taxa.txt 2> /dev/null
	for i in $(ls *_genus_unique_uncultured.txt);do (
	   mv "$i" ${i%*_genus_unique_uncultured*}_unique_sequences.txt 2> /dev/null
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_genus_unique_uncultured.txt *_genus_inter.txt *_genus_inter2.txt *_duplicate_count.txt *_multialign_genus_reads.txt *_genus_duplicates.txt.gz 2> /dev/null

	for i in $(ls *_genus_OTU.txt);do (
	   cat "$i" ${i%_genus_OTU*}_unique_sequences.txt > ${i%_genus_OTU*}_complete_genus_reads.txt ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	echo -e "${YELLOW}- compiling taxonomic information"
	for i in $(ls *_complete_genus_reads.txt);do (
	   awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_genus_reads*}_sighits_temp.txt
	    echo $'abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	     cat - ${i%_complete_genus_reads*}_sighits_temp.txt > ${i%_complete_genus_reads*}_sighits_temp2.txt
			 ) &
			 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
				 wait
			 fi
	done
	wait

	for i in $(ls *_sighits_temp2.txt);do (
	   awk -F '\t' '{gsub(/ /,"_");print}' $i > ${i%_sighits_temp2*}_sighits.txt
     $gzip ${i%_sighits_temp2*}_sighits.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_complete_genus_reads.txt *_sighits_temp.txt *_unique_reads.txt.gz *_unique_sequences.txt *_sighits_temp2.txt *_genus_OTU.txt
fi

cd "${projdir}"/metagenome/sighits/sighits_genus
find . -type f -name '*_sighits.txt.gz' -exec cat {} + > sighits.txt.gz
awk '{print $8}' <(zcat sighits.txt.gz) | awk '{gsub(";","\n"); print}' | sort -T "${projdir}"/tmp -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt.gz
awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  ${projdir}/rankedlineage_edited.dmp taxids_sighits.txt | \
awk '{gsub(/ /,"_"); print }' | awk '{print $4, $5, $6, $7, $8, $9, $10}' | awk 'NR == 1; NR > 1 {print $0}' | awk '{gsub(/ /,"\t");}1' > rankedlineage_subhits.txt
rm taxids_sighits.txt

cd "${projdir}"/metagenome/sighits/sighits_genus/
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/genus/' > genus_taxa_mean_temp1.txt
echo -e 'genus' | cat - genus_taxa_mean_temp1.txt > genus_taxa_mean_temp.txt && rm genus_taxa_mean_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/genus/'  > genus_taxa_unique_sequences_temp1.txt
echo -e 'genus' | cat - genus_taxa_unique_sequences_temp1.txt > genus_taxa_unique_sequences_temp.txt && rm genus_taxa_unique_sequences_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt |  sort -T "${projdir}"/tmp -u | awk '!/genus/'  > genus_taxa_quantification_accuracy_temp1.txt
echo -e 'genus' | cat - genus_taxa_quantification_accuracy_temp1.txt > genus_taxa_quantification_accuracy_temp.txt && rm genus_taxa_quantification_accuracy_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/genus/'  > genus_taxa_rel_quantification_accuracy_temp1.txt
echo -e 'genus' | cat - genus_taxa_rel_quantification_accuracy_temp1.txt > genus_taxa_rel_quantification_accuracy_temp.txt && rm genus_taxa_rel_quantification_accuracy_temp1.txt

genus_level=genus
for i in $(ls *_sighits.txt.gz);do
gunzip $i
	Rscript ${Qmatey_dir}/scripts/stats_summary.R ${i%.gz} $genus_level "${Qmatey_dir}/tools/R"
wait
$gzip ${i%.gz}
echo $'genus\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdmean.txt genus_taxa_mean_temp.txt > holdmean2.txt
cat holdmean2.txt > genus_taxa_mean_temp.txt
id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holduniq_reads.txt genus_taxa_unique_sequences_temp.txt > holduniq_reads2.txt
cat holduniq_reads2.txt > genus_taxa_unique_sequences_temp.txt
id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt genus_taxa_quantification_accuracy_temp.txt > holdstderr2.txt
cat holdstderr2.txt > genus_taxa_quantification_accuracy_temp.txt
id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdrelstderr.txt genus_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt
cat holdrelstderr2.txt > genus_taxa_rel_quantification_accuracy_temp.txt
wait
rm *stats1* *stats2* *hold*
done
wait

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' genus_taxa_mean_temp.txt > genus_taxa_mean_temp2.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_mean_temp2.txt > ../../results/genus_level/genus_taxa_mean.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' genus_taxa_unique_sequences_temp.txt > genus_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_unique_sequences_temp2.txt > ../../results/genus_level/genus_taxa_unique_sequences.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_quantification_accuracy_temp.txt > genus_taxa_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/genus_level/genus_taxa_mean.txt genus_taxa_quantification_accuracy_temp2.txt > ../../results/genus_level/genus_taxa_quantification_accuracy.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_rel_quantification_accuracy_temp.txt > genus_taxa_rel_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/genus_level/genus_taxa_mean.txt genus_taxa_rel_quantification_accuracy_temp2.txt > ../../results/genus_level/genus_taxa_rel_quantification_accuracy.txt
rm *_temp*

cd "${projdir}"/metagenome/results/genus_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do
  awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_genus/rankedlineage_subhits.txt genus_taxa_"${i}".txt | \
  awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > genus_taxainfo_"${i}".txt &&
  wait
done
wait
rm *_taxa_*


awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
genus_taxainfo_mean.txt | awk 'BEGIN{OFS="\t"} {$1=$1};1' > genus_taxainfo_mean_holdingtaxinfo.txt
touch genus_taxainfo_mean_buildnorm.txt
for i in ../../../metagenome/haplotig/*_metagenome.fasta.gz; do
  sample=${i%_metagenome.fasta.gz}; sample=${sample##*/}
  if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
    normfactor=1
  else
    normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
  fi
  awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" genus_taxainfo_mean.txt | \
  paste - genus_taxainfo_mean_buildnorm.txt > genus_taxainfo_mean_buildnorm0.txt &&
  mv genus_taxainfo_mean_buildnorm0.txt genus_taxainfo_mean_buildnorm.txt
done
wait
paste genus_taxainfo_mean_buildnorm.txt > genus_taxainfo_mean_norm0.txt
paste genus_taxainfo_mean_norm0.txt genus_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > genus_taxainfo_mean_normalized.txt
# pos=$(awk -v RS='\t' '/genus/{print NR; exit}' genus_taxainfo_mean_normalized.txt)
# awk -v pat="$pos" -F'\t'  'BEGIN {OFS="\t"}; {k=$pat; $pat=""; print k,$0}' genus_taxainfo_mean_normalized.txt | awk 'BEGIN{FS="\t+"; OFS="\t"} {$1=$1; print}' > genus_taxainfo_mean_normalized.tmp
# mv genus_taxainfo_mean_normalized.tmp genus_taxainfo_mean_normalized.txt
rm genus_taxainfo_mean_buildnorm.txt genus_taxainfo_mean_holdingtaxinfo.txt genus_taxainfo_mean_norm0.txt 2> /dev/null

for i in $(ls *.txt); do
  awk '{gsub(/-/,"_"); print}' $i > ${i%.txt}.temp &&
  mv ${i%.txt}.temp $i
done
wait

Rscript "${Qmatey_dir}/scripts/cleanup_rank.R" "genus" "${Qmatey_dir}/tools/R" 2>/dev/null
wait


# if [[ "$genome_scaling" == true ]]; then
# 	if [[ "$subsample_shotgun_R1" == false ]]; then
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" genus "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "true" &>/dev/null
# 	else
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" genus "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "false" &>/dev/null
# 	fi
# else
# 		Rscript "${Qmatey_dir}/scripts/error_zero_inflation.R" genus "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated &>/dev/null
# fi

file=${projdir}/exclude_taxa.txt
# if test -f $file; then
#   cat genus_taxainfo_mean.txt > genus_taxainfo_mean_filtered.txt &&
#   cat genus_taxainfo_unique_sequences.txt > genus_taxainfo_unique_sequences_filtered.txt &&
#   cat genus_taxainfo_quantification_accuracy.txt > genus_taxainfo_quantification_accuracy_filtered.txt &&
#   cat genus_taxainfo_rel_quantification_accuracy.txt > genus_taxainfo_rel_quantification_accuracy_filtered.txt &&
#   while read -r line; do
#     for i in $( ls *filtered.txt ); do
#       awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i
#     done
# 		wait
#   done < $file
# 	wait
# fi

if test -f $file; then
	echo -e "${YELLOW}- creating genus-level visualizations"
	cd "${projdir}"/metagenome/results/genus_level
	genus_level_mean=genus_taxainfo_mean_filtered.txt
	genus_level_uniq=genus_taxainfo_unique_sequences_filtered.txt
	genus_level_stderr=genus_taxainfo_quantification_accuracy_filtered.txt
	genus_level_rel_stderr=genus_taxainfo_rel_quantification_accuracy_filtered.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi


	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/genus_level_boxplots.R" "$genus_level_mean" "$genus_level_uniq" "$genus_level_stderr" "$genus_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/genus_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/genus_level/boxplots/ 2> /dev/null
else
	echo -e "${YELLOW}- creating genus-level visualizations"
	cd "${projdir}"/metagenome/results/genus_level
	genus_level_mean=genus_taxainfo_mean.txt
	genus_level_mean_norm=genus_taxainfo_mean_normalized.txt
	genus_level_uniq=genus_taxainfo_unique_sequences.txt
	genus_level_stderr=genus_taxainfo_quantification_accuracy.txt
	genus_level_rel_stderr=genus_taxainfo_rel_quantification_accuracy.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/genus_level_boxplots.R" "$genus_level_mean" "$genus_level_mean_norm" "$genus_level_uniq" "$genus_level_stderr" "$genus_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/genus_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/genus_level/boxplots/ 2> /dev/null
fi

}
if [[ "$genus_level" == "true" ]] && [[ -z "$(ls -A ${projdir}/metagenome/results/genus_level/genus_taxainfo* 2> /dev/null)" ]]; then
	if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_genus_level/genus_taxainfo* 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_genus_level/genus_taxainfo* 2> /dev/null)" ]]; then
			cd "${projdir}"/metagenome/alignment/
			mv ./uncultured/uncultured*megablast.gz ./ &&
			mv ./uncultured_combined_compressed.megablast.gz ./uncultured/ &&
			for i in $(ls *megablast.gz); do mv "$i" ${i#uncultured_}; done
			wait
			cd "${projdir}"
			time genus_level 2>> ${projdir}/log.out
			wait
			cd "${projdir}"/metagenome/alignment/
			for i in $(ls *megablast.gz); do mv "$i" ./uncultured/uncultured_"${i}"; done
			wait
			mv "${projdir}"/metagenome/results/genus_level ${projdir}/metagenome/results/uncultured_genus_level
			mv "${projdir}"/metagenome/sighits/sighits_genus ${projdir}/metagenome/sighits/uncultured_sighits_genus
		fi
	fi
	if [[ -z "$(ls -A ${projdir}/metagenome/results/genus_level/genus_taxainfo* 2> /dev/null)" ]]; then
		cd "${projdir}"/metagenome/alignment/cultured/
		mv *megablast* ../
		mv ../combined_compressed.megablast.gz ./
		cd "${projdir}"
		time genus_level 2>> ${projdir}/log.out
		mv "${projdir}"/metagenome/alignment/*megablast* ${projdir}/metagenome/alignment/cultured/
	fi
else
	echo -e "${magenta}- skipping exact matching for genus-level profiling ${white}\n"
fi

############################################################################
family_level() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Family-Level classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact matching for family-level profiling"
cd "${projdir}"/metagenome/sighits
mkdir -p sighits_family
cd "${projdir}"/metagenome/results
mkdir -p family_level
cd "${projdir}"/metagenome/alignment
if find ../sighits/sighits_family/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits of at family-level already available for each sample"
else
	for i in $(ls *_haplotig.megablast.gz);do (
		if [[ ! -f "../sighits/sighits_family/${i%_haplotig.megablast.gz}_sighits.txt.gz" ]]; then
			if [[ "$taxids" == true ]]; then
			  for emg in ${projdir}/taxids/*.txids; do
					if [[ "$fullqlen_alignment" == "false" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-6 {print $0}' <(zcat "$i" | awk '$6>=97') <(zcat "$i" | awk '$6>=97') | awk '$3 >= 32 {print $0}' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_family/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
					if [[ "$fullqlen_alignment" == "true" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-6 {print $0}' <(zcat "$i" | awk '$6>=97') <(zcat "$i" | awk '$6>=97') | awk '$3 >= ($5-($5*0.03))' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_family/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
			  done
				zcat ../sighits/sighits_family/${i%_haplotig.megablast.gz}_sighits.txt.gz | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | gzip > ../sighits/sighits_family/${i%_haplotig.megablast.gz}_sighits.tmp.gz &&
				mv ../sighits/sighits_family/${i%_haplotig.megablast.gz}_sighits.tmp.gz ../sighits/sighits_family/${i%_haplotig.megablast.gz}_sighits.txt.gz
				wait
			else
				if [[ "$fullqlen_alignment" == "false" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-6 {print $0}' <(zcat "$i" | awk '$6>=97') <(zcat "$i" | awk '$6>=97') | awk '$3 >= 32 {print $0}' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_family/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
				if [[ "$fullqlen_alignment" == "true" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-6 {print $0}' <(zcat "$i" | awk '$6>=97') <(zcat "$i" | awk '$6>=97') | awk '$3 >= ($5-($5*0.03))' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_family/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
			fi
		fi
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_family
	for i in $(ls *_sighits.txt.gz);do (
		zcat "$i" | awk -F '\t' 'NR>1{a[$10]++;b[$10]=$0}END{for(x in a)if(a[x]==1)print b[x]}' > ${i%_sighits*}_family_unique_reads.txt
		$gzip ${i%_sighits*}_family_unique_reads.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

  echo -e "${YELLOW}- compiling family-level multi-alignment algorithm"

  cd "${projdir}"/metagenome/sighits/sighits_family
  		for i in $(ls *_sighits.txt.gz);do (
				gunzip $i
				gunzip ${i%_sighits*}_family_unique_reads.txt.gz
				awk -F '\t' 'FNR==NR{a[$10]=1; next} {print $0}' ${i%_sighits*}_family_unique_reads.txt <(awk 'NR>1{print $0}' ${i%.gz}) OFS='\t' > ${i%_sighits*}_dup.txt
				rm "${i%.gz}"
				$gzip ${i%_sighits*}_family_unique_reads.txt
				$gzip ${i%_sighits*}_dup.txt
				) &
				if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
					wait
				fi
  		done
  		wait


  cd "${projdir}"/metagenome/sighits/sighits_family
  for i in $(ls *_dup.txt.gz);do (
  	awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_dup*}_taxids_dup.txt

		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

  cd "${projdir}"/metagenome/sighits/sighits_family

  for i in $(ls *_taxids_dup.txt); do
  	awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		cut -f 2-11 -d$'\t' > ${i%_taxids_dup*}_dup_inter.txt
  done
  wait

	for i in $(ls *_dup_inter.txt);do (
	   awk -F '\t'  '{print $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_family_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait
	rm *_taxids_dup.txt

  for i in $(ls *_dup.txt.gz);do (
    paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <(zcat "$i") ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6}' OFS='\t' ${i%*_dup.txt.gz}_family_taxa.txt) | \
		awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' OFS='\t' > ${i%_dup*}_family_duplicates.txt
		$gzip ${i%_dup*}_family_duplicates.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait
	for i in $(ls *_family_duplicates.txt.gz);do (
		 awk -F '\t' '$11 != "NA"' OFS='\t' <(zcat "$i") > ${i%*_family_duplicates.txt.gz}temp.txt
		 $gzip ${i%*_family_duplicates.txt.gz}temp.txt
 		mv ${i%*_family_duplicates.txt.gz}temp.txt.gz $i
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
	wait

  rm *_dup_inter.txt *_dup.txt.gz *_family_taxa.txt

  for i in $(ls *_family_duplicates.txt.gz);do (
    awk -F '\t' '{print $1, $10"~"$11, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, $14, $15, $16}' OFS='\t' <(zcat "$i") > ${i%_family_duplicates*}_family_inter.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

  wait
  for i in $(ls *_family_inter.txt);do (
    awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' OFS='\t' $i > ${i%_family_inter*}_family_inter2.txt ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait


  for i in $(ls *_family_inter2.txt);do (
    awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -T "${projdir}"/tmp -k1,1  > ${i%_family_inter2*}_duplicate_count.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_duplicate_count.txt);do (
    awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -T "${projdir}"/tmp -k1,1 > ${i%_duplicate_count*}_multialign_family_reads.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_multialign_family_reads.txt);do (
    awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_family_reads*}_family_inter2.txt | sort -T "${projdir}"/tmp -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $4, $5, $6, $7, $8, $9, $10, $11, $1, $2, $12, $13, $14, $15, $16}' OFS='\t' > ${i%_multialign_family_reads*}_family_OTU.txt
    ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  wait
  for i in $(ls *_family_unique_reads.txt.gz);do (
    awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_family_unique_reads*}_taxids_uniq.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_family
	for i in $(ls *_taxids_uniq.txt); do
	   awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		 cut -f 2-11 -d$'\t' > ${i%_taxids_uniq*}_uniq_inter.txt
	done
	wait

  for i in $(ls *_uniq_inter.txt);do (
     awk -F '\t'  '{print $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_family_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait

  rm *_taxids_uniq.txt


  for i in $(ls *_unique_reads.txt.gz);do (
     paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <( zcat "$i")) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6}' OFS='\t' "${i%*_family_unique_reads*}"_family_taxa.txt) | \
		 awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16}' OFS='\t' > "${i%_family_uniq*}"_family_unique_uncultured.txt &&
		 wait ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait
	for i in $(ls *_family_unique_uncultured.txt);do (
		awk -F '\t' '$11 != "NA"' OFS='\t' "$i" > "${i%*_family_unique_uncultured.txt}"temp.txt &&
		:> "$i" &&
		mv "${i%*_family_unique_uncultured.txt}"temp.txt "$i" 2> /dev/null
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
	wait
  rm *_uniq_inter.txt && rm *_family_taxa.txt 2> /dev/null
  for i in $(ls *_family_unique_uncultured.txt);do (
     mv "$i" ${i%*_family_unique_uncultured*}_unique_sequences.txt 2> /dev/null
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait

  rm *_family_unique_uncultured.txt *_family_inter.txt *_family_inter2.txt *_duplicate_count.txt *_multialign_family_reads.txt *_family_duplicates.txt.gz 2> /dev/null

	for i in $(ls *_family_OTU.txt);do (
	   cat "$i" ${i%_family_OTU*}_unique_sequences.txt > ${i%_family_OTU*}_complete_family_reads.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	echo -e "${YELLOW}- compiling taxonomic information"
	for i in $(ls *_complete_family_reads.txt);do (
	   awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_family_reads*}_sighits_temp.txt
	    echo $'abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
	     cat - ${i%_complete_family_reads*}_sighits_temp.txt > ${i%_complete_family_reads*}_sighits_temp2.txt
			 ) &
			 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
				 wait
			 fi
	done
	wait

	for i in $(ls *_sighits_temp2.txt);do (
	   awk -F '\t' '{gsub(/ /,"_");print}' $i > ${i%_sighits_temp2*}_sighits.txt
     $gzip ${i%_sighits_temp2*}_sighits.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_complete_family_reads.txt *_sighits_temp.txt *_unique_reads.txt.gz *_unique_sequences.txt *_sighits_temp2.txt *_family_OTU.txt
fi


cd "${projdir}"/metagenome/sighits/sighits_family
find . -type f -name '*_sighits.txt.gz' -exec cat {} + > sighits.txt.gz
awk '{print $8}' <(zcat sighits.txt.gz) | awk '{gsub(";","\n"); print}' | sort -T "${projdir}"/tmp -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt.gz
awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  ${projdir}/rankedlineage_edited.dmp taxids_sighits.txt > rankedlineage_subhits_temp.txt
awk '{gsub(/ /,"_");}1' rankedlineage_subhits_temp.txt | awk '{print $5, $6, $7, $8, $9, $10}' | awk 'NR == 1; NR > 1 {print $0}' | awk '{gsub(/ /,"\t");}1' > rankedlineage_subhits.txt
rm taxids_sighits.txt rankedlineage_subhits_temp.txt

cd "${projdir}"/metagenome/sighits/sighits_family/
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/family/' > family_taxa_mean_temp1.txt
echo -e 'family' | cat - family_taxa_mean_temp1.txt > family_taxa_mean_temp.txt && rm family_taxa_mean_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/family/'  > family_taxa_unique_sequences_temp1.txt
echo -e 'family' | cat - family_taxa_unique_sequences_temp1.txt > family_taxa_unique_sequences_temp.txt && rm family_taxa_unique_sequences_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt |  sort -T "${projdir}"/tmp -u | awk '!/family/'  > family_taxa_quantification_accuracy_temp1.txt
echo -e 'family' | cat - family_taxa_quantification_accuracy_temp1.txt > family_taxa_quantification_accuracy_temp.txt && rm family_taxa_quantification_accuracy_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/family/'  > family_taxa_rel_quantification_accuracy_temp1.txt
echo -e 'family' | cat - family_taxa_rel_quantification_accuracy_temp1.txt > family_taxa_rel_quantification_accuracy_temp.txt && rm family_taxa_rel_quantification_accuracy_temp1.txt

family_level=family
for i in $(ls *_sighits.txt.gz);do
	gunzip $i
 	Rscript ${Qmatey_dir}/scripts/stats_summary.R ${i%.gz} $family_level "${Qmatey_dir}/tools/R"
	wait
	$gzip ${i%.gz}
  echo $'family\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
  id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdmean.txt family_taxa_mean_temp.txt > holdmean2.txt
  cat holdmean2.txt > family_taxa_mean_temp.txt
  id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holduniq_reads.txt family_taxa_unique_sequences_temp.txt > holduniq_reads2.txt
  cat holduniq_reads2.txt > family_taxa_unique_sequences_temp.txt
  id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt family_taxa_quantification_accuracy_temp.txt > holdstderr2.txt
  cat holdstderr2.txt > family_taxa_quantification_accuracy_temp.txt
  id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdrelstderr.txt family_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt
  cat holdrelstderr2.txt > family_taxa_rel_quantification_accuracy_temp.txt
	wait
  rm *stats1* *stats2* *hold*
done
wait

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' family_taxa_mean_temp.txt > family_taxa_mean_temp2.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_mean_temp2.txt > ../../results/family_level/family_taxa_mean.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' family_taxa_unique_sequences_temp.txt > family_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_unique_sequences_temp2.txt > ../../results/family_level/family_taxa_unique_sequences.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_quantification_accuracy_temp.txt > family_taxa_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/family_level/family_taxa_mean.txt family_taxa_quantification_accuracy_temp2.txt > ../../results/family_level/family_taxa_quantification_accuracy.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_rel_quantification_accuracy_temp.txt > family_taxa_rel_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/family_level/family_taxa_mean.txt family_taxa_rel_quantification_accuracy_temp2.txt > ../../results/family_level/family_taxa_rel_quantification_accuracy.txt
rm *_temp*

cd "${projdir}"/metagenome/results/family_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do
  awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_family/rankedlineage_subhits.txt family_taxa_"${i}".txt | \
  awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > family_taxainfo_"${i}".txt &&
  wait
done
wait
rm *_taxa_*


awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
family_taxainfo_mean.txt | awk 'BEGIN{OFS="\t"} {$1=$1};1' > family_taxainfo_mean_holdingtaxinfo.txt
touch family_taxainfo_mean_buildnorm.txt
for i in ../../../metagenome/haplotig/*_metagenome.fasta.gz; do
  sample=${i%_metagenome.fasta.gz}; sample=${sample##*/}
  if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
    normfactor=1
  else
    normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
  fi
  awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" family_taxainfo_mean.txt | \
  paste - family_taxainfo_mean_buildnorm.txt > family_taxainfo_mean_buildnorm0.txt &&
  mv family_taxainfo_mean_buildnorm0.txt family_taxainfo_mean_buildnorm.txt
done
wait
paste family_taxainfo_mean_buildnorm.txt > family_taxainfo_mean_norm0.txt
paste family_taxainfo_mean_norm0.txt family_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > family_taxainfo_mean_normalized.txt
# pos=$(awk -v RS='\t' '/family/{print NR; exit}' family_taxainfo_mean_normalized.txt)
# awk -v pat="$pos" -F'\t'  'BEGIN {OFS="\t"}; {k=$pat; $pat=""; print k,$0}' family_taxainfo_mean_normalized.txt | awk 'BEGIN{FS="\t+"; OFS="\t"} {$1=$1; print}' > family_taxainfo_mean_normalized.tmp
# mv family_taxainfo_mean_normalized.tmp family_taxainfo_mean_normalized.txt
rm family_taxainfo_mean_buildnorm.txt family_taxainfo_mean_holdingtaxinfo.txt family_taxainfo_mean_norm0.txt 2> /dev/null

for i in $(ls *.txt); do
  awk '{gsub(/-/,"_"); print}' $i > ${i%.txt}.temp &&
  mv ${i%.txt}.temp $i
done
wait


# if [[ "$genome_scaling" == true ]]; then
# 	if [[ "$subsample_shotgun_R1" == false ]]; then
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" family "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "true" &>/dev/null
# 	else
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" family "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "false" &>/dev/null
# 	fi
# else
# 		Rscript "${Qmatey_dir}/scripts/error_zero_inflation.R" family "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated &>/dev/null
# fi

file=${projdir}/exclude_taxa.txt
# if test -f $file; then
#   cat family_taxainfo_mean.txt > family_taxainfo_mean_filtered.txt &&
#   cat family_taxainfo_unique_sequences.txt > family_taxainfo_unique_sequences_filtered.txt &&
#   cat family_taxainfo_quantification_accuracy.txt > family_taxainfo_quantification_accuracy_filtered.txt &&
#   cat family_taxainfo_rel_quantification_accuracy.txt > family_taxainfo_rel_quantification_accuracy_filtered.txt &&
#   while read -r line; do
#     for i in $( ls *filtered.txt ); do
#       awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i
#     done
# 		wait
#   done < $file
# 	wait
# fi

if test -f $file; then
	echo -e "${YELLOW}- creating family-level visualizations"
	cd "${projdir}"/metagenome/results/family_level
	family_level_mean=family_taxainfo_mean_filtered.txt
	family_level_uniq=family_taxainfo_unique_sequences_filtered.txt
	family_level_stderr=family_taxainfo_quantification_accuracy_filtered.txt
	family_level_rel_stderr=family_taxainfo_rel_quantification_accuracy_filtered.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi


	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/family_level_boxplots.R" "$family_level_mean" "$family_level_uniq" "$family_level_stderr" "$family_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/family_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/family_level/boxplots/ 2> /dev/null
else
	echo -e "${YELLOW}- creating family-level visualizations"
	cd "${projdir}"/metagenome/results/family_level
	family_level_mean=family_taxainfo_mean.txt
	family_level_mean_norm=family_taxainfo_mean_normalized.txt
	family_level_uniq=family_taxainfo_unique_sequences.txt
	family_level_stderr=family_taxainfo_quantification_accuracy.txt
	family_level_rel_stderr=family_taxainfo_rel_quantification_accuracy.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/family_level_boxplots.R" "$family_level_mean" "$family_level_mean_norm" "$family_level_uniq" "$family_level_stderr" "$family_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/family_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/family_level/boxplots/ 2> /dev/null
fi

}
if [[ "$family_level" == "true" ]] && [[ -z "$(ls -A ${projdir}/metagenome/results/family_level/family_taxainfo* 2> /dev/null)" ]]; then
	if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_family_level/family_taxainfo* 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_family_level/family_taxainfo* 2> /dev/null)" ]]; then
			cd "${projdir}"/metagenome/alignment/
			mv ./uncultured/uncultured*megablast.gz ./ &&
			mv ./uncultured_combined_compressed.megablast.gz ./uncultured/ &&
			for i in $(ls *megablast.gz); do mv "$i" ${i#uncultured_}; done
			wait
			cd "${projdir}"
			time family_level 2>> ${projdir}/log.out
			wait
			cd "${projdir}"/metagenome/alignment/
			for i in $(ls *megablast.gz); do mv "$i" ./uncultured/uncultured_"${i}"; done
			wait
			mv "${projdir}"/metagenome/results/family_level ${projdir}/metagenome/results/uncultured_family_level
			mv "${projdir}"/metagenome/sighits/sighits_family ${projdir}/metagenome/sighits/uncultured_sighits_family
		fi
	fi
	if [[ -z "$(ls -A ${projdir}/metagenome/results/family_level/family_taxainfo* 2> /dev/null)" ]]; then
		cd "${projdir}"/metagenome/alignment/cultured/
		mv *megablast* ../
		mv ../combined_compressed.megablast.gz ./
		cd "${projdir}"
		time family_level 2>> ${projdir}/log.out
		mv "${projdir}"/metagenome/alignment/*megablast* ${projdir}/metagenome/alignment/cultured/
	fi
else
	echo -e "${magenta}- skipping exact matching for family-level profiling ${white}\n"
fi

#########################################################################
order_level() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Order-Level classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact matching for order-level profiling"
cd "${projdir}"/metagenome/sighits
mkdir -p sighits_order
cd "${projdir}"/metagenome/results
mkdir -p order_level
cd "${projdir}"/metagenome/alignment
if find ../sighits/sighits_order/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits of at order-level already available for each sample"
else
	for i in $(ls *_haplotig.megablast.gz);do (
		if [[ ! -f "../sighits/sighits_order/${i%_haplotig.megablast.gz}_sighits.txt.gz" ]]; then
			if [[ "$taxids" == true ]]; then
			  for emg in ${projdir}/taxids/*.txids; do
					if [[ "$fullqlen_alignment" == "false" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-8 {print $0}' <(zcat "$i" | awk '$6>=96') <(zcat "$i" | awk '$6>=96') | awk '$3 >= 32 {print $0}' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_order/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
					if [[ "$fullqlen_alignment" == "true" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-8 {print $0}' <(zcat "$i" | awk '$6>=96') <(zcat "$i" | awk '$6>=96') | awk '$3 >= ($5-($5*0.04))' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_order/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
			  done
				zcat ../sighits/sighits_order/${i%_haplotig.megablast.gz}_sighits.txt.gz | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | gzip > ../sighits/sighits_order/${i%_haplotig.megablast.gz}_sighits.tmp.gz &&
				mv ../sighits/sighits_order/${i%_haplotig.megablast.gz}_sighits.tmp.gz ../sighits/sighits_order/${i%_haplotig.megablast.gz}_sighits.txt.gz
				wait
			else
				if [[ "$fullqlen_alignment" == "false" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-8 {print $0}' <(zcat "$i" | awk '$6>=96') <(zcat "$i" | awk '$6>=96') | awk '$3 >= 32 {print $0}' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_order/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
				if [[ "$fullqlen_alignment" == "true" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-8 {print $0}' <(zcat "$i" | awk '$6>=96') <(zcat "$i" | awk '$6>=96') | awk '$3 >= ($5-($5*0.04))' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_order/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
			fi
		fi
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_order
	for i in $(ls *_sighits.txt.gz);do (
		zcat "$i" | awk -F '\t' 'NR>1{a[$10]++;b[$10]=$0}END{for(x in a)if(a[x]==1)print b[x]}' > ${i%_sighits*}_order_unique_reads.txt
		$gzip ${i%_sighits*}_order_unique_reads.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

  echo -e "${YELLOW}- compiling order-level multi-alignment algorithm"

  cd "${projdir}"/metagenome/sighits/sighits_order
  		for i in $(ls *_sighits.txt.gz);do (
				gunzip $i
				gunzip ${i%_sighits*}_order_unique_reads.txt.gz
				awk -F '\t' 'FNR==NR{a[$10]=1; next} {print $0}' ${i%_sighits*}_order_unique_reads.txt <(awk 'NR>1{print $0}' ${i%.gz}) OFS='\t' > ${i%_sighits*}_dup.txt
				rm "${i%.gz}"
				$gzip ${i%_sighits*}_order_unique_reads.txt
				$gzip ${i%_sighits*}_dup.txt
				) &
				if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
					wait
				fi
  		done
  		wait


  cd "${projdir}"/metagenome/sighits/sighits_order
  for i in $(ls *_dup.txt.gz);do (
  	awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_dup*}_taxids_dup.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait


  cd "${projdir}"/metagenome/sighits/sighits_order

  for i in $(ls *_taxids_dup.txt); do
  	awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		cut -f 2-11 -d$'\t' > ${i%_taxids_dup*}_dup_inter.txt
  done
  wait

	for i in $(ls *_dup_inter.txt);do (
	   awk -F '\t'  '{print $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_order_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait
	rm *_taxids_dup.txt

  for i in $(ls *_dup.txt.gz);do (
    paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <(zcat "$i") ) <(awk -F '\t' '{print $1, $2, $3, $4, $5}' OFS='\t' ${i%*_dup.txt.gz}_order_taxa.txt) | \
		awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' OFS='\t' > ${i%_dup*}_order_duplicates.txt
		$gzip ${i%_dup*}_order_duplicates.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait
	for i in $(ls *_order_duplicates.txt.gz);do (
		 awk -F '\t' '$11 != "NA"' OFS='\t' <(zcat "$i") > ${i%*_order_duplicates.txt.gz}temp.txt
		 $gzip ${i%*_order_duplicates.txt.gz}temp.txt
 		mv ${i%*_order_duplicates.txt.gz}temp.txt.gz $i
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
	wait
  rm *_dup_inter.txt *_dup.txt.gz *_order_taxa.txt

  for i in $(ls *_order_duplicates.txt.gz);do (
    awk -F '\t' '{print $1, $10"~"$11, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, $14, $15}' OFS='\t' <(zcat "$i") > ${i%_order_duplicates*}_order_inter.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done

  wait
  for i in $(ls *_order_inter.txt);do (
    awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' OFS='\t' $i > ${i%_order_inter*}_order_inter2.txt ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
  wait

  for i in $(ls *_order_inter2.txt);do (
    awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -T "${projdir}"/tmp -k1,1  > ${i%_order_inter2*}_duplicate_count.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_duplicate_count.txt);do (
    awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -T "${projdir}"/tmp -k1,1 > ${i%_duplicate_count*}_multialign_order_reads.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_multialign_order_reads.txt);do (
    awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_order_reads*}_order_inter2.txt | sort -T "${projdir}"/tmp -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $4, $5, $6, $7, $8, $9, $10, $11, $1, $2, $12, $13, $14, $15}' OFS='\t' > ${i%_multialign_order_reads*}_order_OTU.txt
	  ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  wait
  for i in $(ls *_order_unique_reads.txt.gz);do (
    awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_order_unique_reads*}_taxids_uniq.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait


	cd "${projdir}"/metagenome/sighits/sighits_order
	for i in $(ls *_taxids_uniq.txt);do
	   awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		 cut -f 2-11 -d$'\t' > ${i%_taxids_uniq*}_uniq_inter.txt
	done
	wait

  for i in $(ls *_uniq_inter.txt);do (
     awk -F '\t'  '{print $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_order_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait

  rm *_taxids_uniq.txt


  for i in $(ls *_unique_reads.txt.gz);do (
     paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <( zcat "$i" )) <(awk -F '\t' '{print $1, $2, $3, $4, $5}' OFS='\t' ${i%*_order_unique_reads*}_order_taxa.txt) | \
		 awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15}' OFS='\t' > ${i%_order_uniq*}_order_unique_uncultured.txt &&
		 wait ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait
	for i in $(ls *_order_unique_uncultured.txt);do (
		awk -F '\t' '$11 != "NA"' OFS='\t' "$i" > "${i%*_order_unique_uncultured.txt}"temp.txt &&
 		:> "$i" &&
		mv "${i%*_order_unique_uncultured.txt}"temp.txt "$i" 2> /dev/null
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
	wait
  rm *_uniq_inter.txt && rm *_order_taxa.txt 2> /dev/null
  for i in $(ls *_order_unique_uncultured.txt);do (
     mv "$i" "${i%*_order_unique_uncultured*}"_unique_sequences.txt 2> /dev/null
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait

	rm *_order_unique_uncultured.txt *_order_inter.txt *_order_inter2.txt *_duplicate_count.txt *_multialign_order_reads.txt *_order_duplicates.txt.gz 2> /dev/null

	for i in $(ls *_order_OTU.txt);do (
	   cat "$i" ${i%_order_OTU*}_unique_sequences.txt > ${i%_order_OTU*}_complete_order_reads.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	echo -e "${YELLOW}- compiling taxonomic information"
	for i in $(ls *_complete_order_reads.txt);do (
	   awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_order_reads*}_sighits_temp.txt
	    echo $'abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\torder\tclass\tphylum\tkingdom\tdomain' | \
	     cat - ${i%_complete_order_reads*}_sighits_temp.txt > ${i%_complete_order_reads*}_sighits_temp2.txt
			 ) &
			 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
				 wait
			 fi
	done
	wait

	for i in $(ls *_sighits_temp2.txt);do (
	   awk -F '\t' '{gsub(/ /,"_");print}' $i > ${i%_sighits_temp2*}_sighits.txt
     $gzip ${i%_sighits_temp2*}_sighits.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_complete_order_reads.txt *_sighits_temp.txt *_unique_reads.txt.gz *_unique_sequences.txt *_sighits_temp2.txt *_order_OTU.txt
fi


cd "${projdir}"/metagenome/sighits/sighits_order
find . -type f -name '*_sighits.txt.gz' -exec cat {} + > sighits.txt.gz
awk '{print $8}' <(zcat sighits.txt.gz) | awk '{gsub(";","\n"); print}' | sort -T "${projdir}"/tmp -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt.gz
awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  ${projdir}/rankedlineage_edited.dmp taxids_sighits.txt > rankedlineage_subhits_temp.txt
awk '{gsub(/ /,"_");}1' rankedlineage_subhits_temp.txt | awk '{print $6, $7, $8, $9, $10}' | awk 'NR == 1; NR > 1 {print $0}' | awk '{gsub(/ /,"\t");}1' > rankedlineage_subhits.txt
rm taxids_sighits.txt rankedlineage_subhits_temp.txt

cd "${projdir}"/metagenome/sighits/sighits_order/
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/order/' > order_taxa_mean_temp1.txt
echo -e 'order' | cat - order_taxa_mean_temp1.txt > order_taxa_mean_temp.txt && rm order_taxa_mean_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/order/'  > order_taxa_unique_sequences_temp1.txt
echo -e 'order' | cat - order_taxa_unique_sequences_temp1.txt > order_taxa_unique_sequences_temp.txt && rm order_taxa_unique_sequences_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt |  sort -T "${projdir}"/tmp -u | awk '!/order/'  > order_taxa_quantification_accuracy_temp1.txt
echo -e 'order' | cat - order_taxa_quantification_accuracy_temp1.txt > order_taxa_quantification_accuracy_temp.txt && rm order_taxa_quantification_accuracy_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/order/'  > order_taxa_rel_quantification_accuracy_temp1.txt
echo -e 'order' | cat - order_taxa_rel_quantification_accuracy_temp1.txt > order_taxa_rel_quantification_accuracy_temp.txt && rm order_taxa_rel_quantification_accuracy_temp1.txt

order_level=order
for i in $(ls *_sighits.txt.gz);do
	gunzip $i
 	Rscript ${Qmatey_dir}/scripts/stats_summary.R ${i%.gz} $order_level "${Qmatey_dir}/tools/R"
	wait
	$gzip ${i%.gz}
  echo $'order\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
  id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdmean.txt order_taxa_mean_temp.txt > holdmean2.txt
  cat holdmean2.txt > order_taxa_mean_temp.txt
  id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holduniq_reads.txt order_taxa_unique_sequences_temp.txt > holduniq_reads2.txt
  cat holduniq_reads2.txt > order_taxa_unique_sequences_temp.txt
  id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt order_taxa_quantification_accuracy_temp.txt > holdstderr2.txt
  cat holdstderr2.txt > order_taxa_quantification_accuracy_temp.txt
  id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdrelstderr.txt order_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt
  cat holdrelstderr2.txt > order_taxa_rel_quantification_accuracy_temp.txt
	wait
  rm *stats1* *stats2* *hold*
done
wait

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' order_taxa_mean_temp.txt > order_taxa_mean_temp2.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_mean_temp2.txt > ../../results/order_level/order_taxa_mean.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' order_taxa_unique_sequences_temp.txt > order_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_unique_sequences_temp2.txt > ../../results/order_level/order_taxa_unique_sequences.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_quantification_accuracy_temp.txt > order_taxa_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/order_level/order_taxa_mean.txt order_taxa_quantification_accuracy_temp2.txt > ../../results/order_level/order_taxa_quantification_accuracy.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_rel_quantification_accuracy_temp.txt > order_taxa_rel_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/order_level/order_taxa_mean.txt order_taxa_rel_quantification_accuracy_temp2.txt > ../../results/order_level/order_taxa_rel_quantification_accuracy.txt
rm *_temp*

cd "${projdir}"/metagenome/results/order_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do
  awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_order/rankedlineage_subhits.txt order_taxa_"${i}".txt | \
  awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > order_taxainfo_"${i}".txt &&
  wait
done
wait
rm *_taxa_*


awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
order_taxainfo_mean.txt | awk 'BEGIN{OFS="\t"} {$1=$1};1' > order_taxainfo_mean_holdingtaxinfo.txt
touch order_taxainfo_mean_buildnorm.txt
for i in ../../../metagenome/haplotig/*_metagenome.fasta.gz; do
  sample=${i%_metagenome.fasta.gz}; sample=${sample##*/}
  if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
    normfactor=1
  else
    normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
  fi
  awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" order_taxainfo_mean.txt | \
  paste - order_taxainfo_mean_buildnorm.txt > order_taxainfo_mean_buildnorm0.txt &&
  mv order_taxainfo_mean_buildnorm0.txt order_taxainfo_mean_buildnorm.txt
done
wait
paste order_taxainfo_mean_buildnorm.txt > order_taxainfo_mean_norm0.txt
paste order_taxainfo_mean_norm0.txt order_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > order_taxainfo_mean_normalized.txt
# pos=$(awk -v RS='\t' '/order/{print NR; exit}' order_taxainfo_mean_normalized.txt)
# awk -v pat="$pos" -F'\t'  'BEGIN {OFS="\t"}; {k=$pat; $pat=""; print k,$0}' order_taxainfo_mean_normalized.txt | awk 'BEGIN{FS="\t+"; OFS="\t"} {$1=$1; print}' > order_taxainfo_mean_normalized.tmp
# mv order_taxainfo_mean_normalized.tmp order_taxainfo_mean_normalized.txt
rm order_taxainfo_mean_buildnorm.txt order_taxainfo_mean_holdingtaxinfo.txt order_taxainfo_mean_norm0.txt 2> /dev/null

for i in $(ls *.txt); do
  awk '{gsub(/-/,"_"); print}' $i > ${i%.txt}.temp &&
  mv ${i%.txt}.temp $i
done
wait


# if [[ "$genome_scaling" == true ]]; then
# 	if [[ "$subsample_shotgun_R1" == false ]]; then
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" order "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "true" &>/dev/null
# 	else
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" order "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "false" &>/dev/null
# 	fi
# else
# 		Rscript "${Qmatey_dir}/scripts/error_zero_inflation.R" order "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated &>/dev/null
# fi

file=${projdir}/exclude_taxa.txt
# if test -f $file; then
#   cat order_taxainfo_mean.txt > order_taxainfo_mean_filtered.txt &&
#   cat order_taxainfo_unique_sequences.txt > order_taxainfo_unique_sequences_filtered.txt &&
#   cat order_taxainfo_quantification_accuracy.txt > order_taxainfo_quantification_accuracy_filtered.txt &&
#   cat order_taxainfo_rel_quantification_accuracy.txt > order_taxainfo_rel_quantification_accuracy_filtered.txt &&
#   while read -r line; do
#     for i in $( ls *filtered.txt ); do
#       awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i
#     done
# 		wait
#   done < $file
# 	wait
# fi

if test -f $file; then
	echo -e "${YELLOW}- creating order-level visualizations"
	cd "${projdir}"/metagenome/results/order_level
	order_level_mean=order_taxainfo_mean_filtered.txt
	order_level_uniq=order_taxainfo_unique_sequences_filtered.txt
	order_level_stderr=order_taxainfo_quantification_accuracy_filtered.txt
	order_level_rel_stderr=order_taxainfo_rel_quantification_accuracy_filtered.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi


	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/order_level_boxplots.R" "$order_level_mean" "$order_level_uniq" "$order_level_stderr" "$order_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/order_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/order_level/boxplots/ 2> /dev/null
else
	echo -e "${YELLOW}- creating order-level visualizations"
	cd "${projdir}"/metagenome/results/order_level
	order_level_mean=order_taxainfo_mean.txt
	order_level_mean_norm=order_taxainfo_mean_normalized.txt
	order_level_uniq=order_taxainfo_unique_sequences.txt
	order_level_stderr=order_taxainfo_quantification_accuracy.txt
	order_level_rel_stderr=order_taxainfo_rel_quantification_accuracy.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/order_level_boxplots.R" "$order_level_mean" "$order_level_mean_norm" "$order_level_uniq" "$order_level_stderr" "$order_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/order_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/order_level/boxplots/ 2> /dev/null
fi

}
if [[ "$order_level" == "true" ]] && [[ -z "$(ls -A ${projdir}/metagenome/results/order_level/order_taxainfo* 2> /dev/null)" ]]; then
	if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_order_level/order_taxainfo* 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_order_level/order_taxainfo* 2> /dev/null)" ]]; then
			cd "${projdir}"/metagenome/alignment/
			mv ./uncultured/uncultured*megablast.gz ./ &&
			mv ./uncultured_combined_compressed.megablast.gz ./uncultured/ &&
			wait
			for i in $(ls *megablast.gz); do mv "$i" ${i#uncultured_}; done
			wait
			cd "${projdir}"
			time order_level 2>> ${projdir}/log.out
			wait
			cd "${projdir}"/metagenome/alignment/
			for i in $(ls *megablast.gz); do mv "$i" ./uncultured/uncultured_"${i}"; done
			wait
			mv "${projdir}"/metagenome/results/order_level ${projdir}/metagenome/results/uncultured_order_level
			mv "${projdir}"/metagenome/sighits/sighits_order ${projdir}/metagenome/sighits/uncultured_sighits_order
		fi
	fi
	if [[ -z "$(ls -A ${projdir}/metagenome/results/order_level/order_taxainfo* 2> /dev/null)" ]]; then
		cd "${projdir}"/metagenome/alignment/cultured/
		mv *megablast* ../
		mv ../combined_compressed.megablast.gz ./
		cd "${projdir}"
		time order_level 2>> ${projdir}/log.out
		mv "${projdir}"/metagenome/alignment/*megablast* ${projdir}/metagenome/alignment/cultured/
	fi
else
	echo -e "${magenta}- skipping exact matching for order-level profiling ${white}\n"
fi

##########################################################################
##########################################################################
class_level() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Class-Level classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact matching for class-level profiling"
cd "${projdir}"/metagenome/sighits
mkdir -p sighits_class
cd "${projdir}"/metagenome/results
mkdir -p class_level
cd "${projdir}"/metagenome/alignment
if find ../sighits/sighits_class/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits of at class-level already available for each sample"
else
	for i in $(ls *_haplotig.megablast.gz);do (
		if [[ ! -f "../sighits/sighits_class/${i%_haplotig.megablast.gz}_sighits.txt.gz" ]]; then
			if [[ "$taxids" == true ]]; then
			  for emg in ${projdir}/taxids/*.txids; do
					if [[ "$fullqlen_alignment" == "false" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-10 {print $0}' <(zcat "$i" | awk '$6>=95') <(zcat "$i" | awk '$6>=95') | awk '$3 >= 32 {print $0}' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_class/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
					if [[ "$fullqlen_alignment" == "true" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-10 {print $0}' <(zcat "$i" | awk '$6>=95') <(zcat "$i" | awk '$6>=95') | awk '$3 >= ($5-($5*0.05))' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_class/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
			  done
				zcat ../sighits/sighits_class/${i%_haplotig.megablast.gz}_sighits.txt.gz | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | gzip > ../sighits/sighits_class/${i%_haplotig.megablast.gz}_sighits.tmp.gz
				mv ../sighits/sighits_class/${i%_haplotig.megablast.gz}_sighits.tmp.gz ../sighits/sighits_class/${i%_haplotig.megablast.gz}_sighits.txt.gz
				wait
			else
				if [[ "$fullqlen_alignment" == "false" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-10 {print $0}' <(zcat "$i" | awk '$6>=95') <(zcat "$i" | awk '$6>=95') | awk '$3 >= 32 {print $0}' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_class/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
				if [[ "$fullqlen_alignment" == "true" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-10 {print $0}' <(zcat "$i" | awk '$6>=95') <(zcat "$i" | awk '$6>=95') | awk '$3 >= ($5-($5*0.05))' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_class/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
			fi
		fi
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_class
	for i in $(ls *_sighits.txt.gz);do (
		zcat "$i" | awk -F '\t' 'NR>1{a[$10]++;b[$10]=$0}END{for(x in a)if(a[x]==1)print b[x]}' > ${i%_sighits*}_class_unique_reads.txt
		$gzip ${i%_sighits*}_class_unique_reads.txt ) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

  echo -e "${YELLOW}- compiling class-level multi-alignment algorithm"

  cd "${projdir}"/metagenome/sighits/sighits_class
  		for i in $(ls *_sighits.txt.gz);do (
				gunzip $i
				gunzip ${i%_sighits*}_class_unique_reads.txt.gz
				awk -F '\t' 'FNR==NR{a[$10]=1; next} {print $0}' ${i%_sighits*}_class_unique_reads.txt <(awk 'NR>1{print $0}' ${i%.gz}) OFS='\t' > ${i%_sighits*}_dup.txt
			rm "${i%.gz}"
				$gzip ${i%_sighits*}_class_unique_reads.txt
				$gzip ${i%_sighits*}_dup.txt
				) &
				if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
					wait
				fi
  		done
  		wait


  cd "${projdir}"/metagenome/sighits/sighits_class
  for i in $(ls *_dup.txt.gz);do (
  	awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_dup*}_taxids_dup.txt

		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait


  cd "${projdir}"/metagenome/sighits/sighits_class

  for i in $(ls *_taxids_dup.txt); do
  	awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		cut -f 2-11 -d$'\t' > ${i%_taxids_dup*}_dup_inter.txt
	done
  wait

	for i in $(ls *_dup_inter.txt);do (
	   awk -F '\t'  '{print $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_class_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait
	rm *_taxids_dup.txt

  for i in $(ls *_dup.txt.gz);do (
    paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <(zcat "$i") ) <(awk -F '\t' '{print $1, $2, $3, $4}' OFS='\t' ${i%*_dup.txt.gz}_class_taxa.txt) | \
		awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' OFS='\t' > ${i%_dup*}_class_duplicates.txt
		$gzip ${i%_dup*}_class_duplicates.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait
	for i in $(ls *_class_duplicates.txt.gz);do (
		 awk -F '\t' '$11 != "NA"' OFS='\t' <(zcat "$i") > ${i%*_class_duplicates.txt.gz}temp.txt
		 $gzip ${i%*_class_duplicates.txt.gz}temp.txt
 		mv ${i%*_class_duplicates.txt.gz}temp.txt.gz $i
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
	wait
  rm *_dup_inter.txt *_dup.txt.gz *_class_taxa.txt

  for i in $(ls *_class_duplicates.txt.gz);do (
    awk -F '\t' '{print $1, $10"~"$11, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13, $14}' OFS='\t' <(zcat "$i") > ${i%_class_duplicates*}_class_inter.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

  wait
  for i in $(ls *_class_inter.txt);do (
    awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' OFS='\t' $i > ${i%_class_inter*}_class_inter2.txt ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
  wait

  for i in $(ls *_class_inter2.txt);do (
    awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -T "${projdir}"/tmp -k1,1  > ${i%_class_inter2*}_duplicate_count.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_duplicate_count.txt);do (
    awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -T "${projdir}"/tmp -k1,1 > ${i%_duplicate_count*}_multialign_class_reads.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_multialign_class_reads.txt);do (
    awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_class_reads*}_class_inter2.txt | sort -T "${projdir}"/tmp -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $4, $5, $6, $7, $8, $9, $10, $11, $1, $2, $12, $13, $14}' OFS='\t' > ${i%_multialign_class_reads*}_class_OTU.txt
		 ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  wait
  for i in $(ls *_class_unique_reads.txt.gz);do (
    awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_class_unique_reads*}_taxids_uniq.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_class
	for i in $(ls *_taxids_uniq.txt); do
	   awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		 cut -f 2-11 -d$'\t' > ${i%_taxids_uniq*}_uniq_inter.txt
	done
	wait

  for i in $(ls *_uniq_inter.txt);do (
     awk -F '\t'  '{print $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_class_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait

  rm *_taxids_uniq.txt


  for i in $(ls *_unique_reads.txt.gz);do (
     paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <( zcat "$i" )) <(awk -F '\t' '{print $1, $2, $3, $4}' OFS='\t' ${i%*_class_unique_reads*}_class_taxa.txt) | \
		 awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14}' OFS='\t' > ${i%_class_uniq*}_class_unique_uncultured.txt &&
		 wait ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait
	for i in $(ls *_class_unique_uncultured.txt);do (
		awk -F '\t' '$11 != "NA"' OFS='\t' "$i" > "${i%*_class_unique_uncultured.txt}"temp.txt &&
 		:> "$i" &&
		mv "${i%*_class_unique_uncultured.txt}"temp.txt "$i" 2> /dev/null
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
  rm *_uniq_inter.txt && rm *_class_taxa.txt 2> /dev/null
  for i in $(ls *_class_unique_uncultured.txt);do (
     mv "$i" "${i%*_class_unique_uncultured*}"_unique_sequences.txt 2> /dev/null
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait

  rm *_class_unique_uncultured.txt *_class_inter.txt *_class_inter2.txt *_duplicate_count.txt *_multialign_class_reads.txt *_class_duplicates.txt.gz 2> /dev/null

	for i in $(ls *_class_OTU.txt);do (
	   cat "$i" ${i%_class_OTU*}_unique_sequences.txt > ${i%_class_OTU*}_complete_class_reads.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	echo -e "${YELLOW}- compiling taxonomic information"
	for i in $(ls *_complete_class_reads.txt);do (
	   awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_class_reads*}_sighits_temp.txt
	    echo $'abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\tclass\tphylum\tkingdom\tdomain' | \
	     cat - ${i%_complete_class_reads*}_sighits_temp.txt > ${i%_complete_class_reads*}_sighits_temp2.txt
			 ) &
			 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
				 wait
			 fi
	done
	wait

	for i in $(ls *_sighits_temp2.txt);do (
	   awk -F '\t' '{gsub(/ /,"_");print}' $i > ${i%_sighits_temp2*}_sighits.txt
     $gzip ${i%_sighits_temp2*}_sighits.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_complete_class_reads.txt *_sighits_temp.txt *_unique_reads.txt.gz *_unique_sequences.txt *_sighits_temp2.txt *_class_OTU.txt
fi


cd "${projdir}"/metagenome/sighits/sighits_class
find . -type f -name '*_sighits.txt.gz' -exec cat {} + > sighits.txt.gz
awk '{print $8}' <(zcat sighits.txt.gz) | awk '{gsub(";","\n"); print}' | sort -T "${projdir}"/tmp -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt.gz
awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  ${projdir}/rankedlineage_edited.dmp taxids_sighits.txt > rankedlineage_subhits_temp.txt
awk '{gsub(/ /,"_");}1' rankedlineage_subhits_temp.txt | awk '{print $7, $8, $9, $10}' | awk 'NR == 1; NR > 1 {print $0}' | awk '{gsub(/ /,"\t");}1' > rankedlineage_subhits.txt
rm taxids_sighits.txt rankedlineage_subhits_temp.txt

cd "${projdir}"/metagenome/sighits/sighits_class/
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/class/' > class_taxa_mean_temp1.txt
echo -e 'class' | cat - class_taxa_mean_temp1.txt > class_taxa_mean_temp.txt && rm class_taxa_mean_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/class/'  > class_taxa_unique_sequences_temp1.txt
echo -e 'class' | cat - class_taxa_unique_sequences_temp1.txt > class_taxa_unique_sequences_temp.txt && rm class_taxa_unique_sequences_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt |  sort -T "${projdir}"/tmp -u | awk '!/class/'  > class_taxa_quantification_accuracy_temp1.txt
echo -e 'class' | cat - class_taxa_quantification_accuracy_temp1.txt > class_taxa_quantification_accuracy_temp.txt && rm class_taxa_quantification_accuracy_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/class/'  > class_taxa_rel_quantification_accuracy_temp1.txt
echo -e 'class' | cat - class_taxa_rel_quantification_accuracy_temp1.txt > class_taxa_rel_quantification_accuracy_temp.txt && rm class_taxa_rel_quantification_accuracy_temp1.txt

class_level=class
for i in $(ls *_sighits.txt.gz);do
	gunzip $i
 	Rscript ${Qmatey_dir}/scripts/stats_summary.R ${i%.gz} $class_level "${Qmatey_dir}/tools/R"
	wait
	$gzip ${i%.gz}
  echo $'class\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
  id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdmean.txt class_taxa_mean_temp.txt > holdmean2.txt
  cat holdmean2.txt > class_taxa_mean_temp.txt
  id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holduniq_reads.txt class_taxa_unique_sequences_temp.txt > holduniq_reads2.txt
  cat holduniq_reads2.txt > class_taxa_unique_sequences_temp.txt
  id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt class_taxa_quantification_accuracy_temp.txt > holdstderr2.txt
  cat holdstderr2.txt > class_taxa_quantification_accuracy_temp.txt
  id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdrelstderr.txt class_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt
  cat holdrelstderr2.txt > class_taxa_rel_quantification_accuracy_temp.txt
	wait
  rm *stats1* *stats2* *hold*
done
wait

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' class_taxa_mean_temp.txt > class_taxa_mean_temp2.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_mean_temp2.txt > ../../results/class_level/class_taxa_mean.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' class_taxa_unique_sequences_temp.txt > class_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_unique_sequences_temp2.txt > ../../results/class_level/class_taxa_unique_sequences.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_quantification_accuracy_temp.txt > class_taxa_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/class_level/class_taxa_mean.txt class_taxa_quantification_accuracy_temp2.txt > ../../results/class_level/class_taxa_quantification_accuracy.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_rel_quantification_accuracy_temp.txt > class_taxa_rel_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/class_level/class_taxa_mean.txt class_taxa_rel_quantification_accuracy_temp2.txt > ../../results/class_level/class_taxa_rel_quantification_accuracy.txt
rm *_temp*

cd "${projdir}"/metagenome/results/class_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do
  awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_class/rankedlineage_subhits.txt class_taxa_"${i}".txt | \
  awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > class_taxainfo_"${i}".txt &&
  wait
done
wait
rm *_taxa_*


awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
class_taxainfo_mean.txt | awk 'BEGIN{OFS="\t"} {$1=$1};1' > class_taxainfo_mean_holdingtaxinfo.txt
touch class_taxainfo_mean_buildnorm.txt
for i in ../../../metagenome/haplotig/*_metagenome.fasta.gz; do
  sample=${i%_metagenome.fasta.gz}; sample=${sample##*/}
  if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
    normfactor=1
  else
    normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
  fi
  awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" class_taxainfo_mean.txt | \
  paste - class_taxainfo_mean_buildnorm.txt > class_taxainfo_mean_buildnorm0.txt &&
  mv class_taxainfo_mean_buildnorm0.txt class_taxainfo_mean_buildnorm.txt
done
wait
paste class_taxainfo_mean_buildnorm.txt > class_taxainfo_mean_norm0.txt
paste class_taxainfo_mean_norm0.txt class_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > class_taxainfo_mean_normalized.txt
# pos=$(awk -v RS='\t' '/class/{print NR; exit}' class_taxainfo_mean_normalized.txt)
# awk -v pat="$pos" -F'\t'  'BEGIN {OFS="\t"}; {k=$pat; $pat=""; print k,$0}' class_taxainfo_mean_normalized.txt | awk 'BEGIN{FS="\t+"; OFS="\t"} {$1=$1; print}' > class_taxainfo_mean_normalized.tmp
# mv class_taxainfo_mean_normalized.tmp class_taxainfo_mean_normalized.txt
rm class_taxainfo_mean_buildnorm.txt class_taxainfo_mean_holdingtaxinfo.txt class_taxainfo_mean_norm0.txt 2> /dev/null

for i in $(ls *.txt); do
  awk '{gsub(/-/,"_"); print}' $i > ${i%.txt}.temp &&
  mv ${i%.txt}.temp $i
done
wait


# if [[ "$genome_scaling" == true ]]; then
# 	if [[ "$subsample_shotgun_R1" == false ]]; then
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" class "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "true" &>/dev/null
# 	else
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" class "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "false" &>/dev/null
# 	fi
# else
# 		Rscript "${Qmatey_dir}/scripts/error_zero_inflation.R" class "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated &>/dev/null
# fi

file=${projdir}/exclude_taxa.txt
# if test -f $file; then
#   cat class_taxainfo_mean.txt > class_taxainfo_mean_filtered.txt &&
#   cat class_taxainfo_unique_sequences.txt > class_taxainfo_unique_sequences_filtered.txt &&
#   cat class_taxainfo_quantification_accuracy.txt > class_taxainfo_quantification_accuracy_filtered.txt &&
#   cat class_taxainfo_rel_quantification_accuracy.txt > class_taxainfo_rel_quantification_accuracy_filtered.txt &&
#   while read -r line; do
#     for i in $( ls *filtered.txt ); do
#       awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i
#     done
# 		wait
#   done < $file
# 	wait
# fi

if test -f $file; then
	echo -e "${YELLOW}- creating class-level visualizations"
	cd "${projdir}"/metagenome/results/class_level
	class_level_mean=class_taxainfo_mean_filtered.txt
	class_level_uniq=class_taxainfo_unique_sequences_filtered.txt
	class_level_stderr=class_taxainfo_quantification_accuracy_filtered.txt
	class_level_rel_stderr=class_taxainfo_rel_quantification_accuracy_filtered.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi


	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/class_level_boxplots.R" "$class_level_mean" "$class_level_uniq" "$class_level_stderr" "$class_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/class_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/class_level/boxplots/ 2> /dev/null
else
	echo -e "${YELLOW}- creating class-level visualizations"
	cd "${projdir}"/metagenome/results/class_level
	class_level_mean=class_taxainfo_mean.txt
	class_level_mean_norm=class_taxainfo_mean_normalized.txt
	class_level_uniq=class_taxainfo_unique_sequences.txt
	class_level_stderr=class_taxainfo_quantification_accuracy.txt
	class_level_rel_stderr=class_taxainfo_rel_quantification_accuracy.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/class_level_boxplots.R" "$class_level_mean" "$class_level_mean_norm" "$class_level_uniq" "$class_level_stderr" "$class_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/class_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/class_level/boxplots/ 2> /dev/null
fi

}
if [[ "$class_level" == "true" ]] && [[ -z "$(ls -A ${projdir}/metagenome/results/class_level/class_taxainfo* 2> /dev/null)" ]]; then
	if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_class_level/class_taxainfo* 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_class_level/class_taxainfo* 2> /dev/null)" ]]; then
			cd "${projdir}"/metagenome/alignment/
			mv ./uncultured/uncultured*megablast.gz ./ &&
			mv ./uncultured_combined_compressed.megablast.gz ./uncultured/ &&
			wait
			for i in $(ls *megablast.gz); do mv "$i" ${i#uncultured_}; done
			wait
			cd "${projdir}"
			time class_level 2>> ${projdir}/log.out
			wait
			cd "${projdir}"/metagenome/alignment/
			for i in $(ls *megablast.gz); do mv "$i" ./uncultured/uncultured_"${i}"; done
			wait
			mv "${projdir}"/metagenome/results/class_level ${projdir}/metagenome/results/uncultured_class_level
			mv "${projdir}"/metagenome/sighits/sighits_class ${projdir}/metagenome/sighits/uncultured_sighits_class
		fi
	fi
	if [[ -z "$(ls -A ${projdir}/metagenome/results/class_level/class_taxainfo* 2> /dev/null)" ]]; then
		cd "${projdir}"/metagenome/alignment/cultured/
		mv *megablast* ../
		mv ../combined_compressed.megablast.gz ./
		cd "${projdir}"
		time class_level 2>> ${projdir}/log.out
		mv "${projdir}"/metagenome/alignment/*megablast* ${projdir}/metagenome/alignment/cultured/
	fi
else
	echo -e "${magenta}- skipping exact matching for class-level profiling ${white}\n"
fi
##########################################################################
phylum_level() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Phylum-Level classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact matching for phylum-level profiling"
cd "${projdir}"/metagenome/sighits
mkdir -p sighits_phylum
cd "${projdir}"/metagenome/results
mkdir -p phylum_level
cd "${projdir}"/metagenome/alignment
if find ../sighits/sighits_phylum/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits of at phylum-level already available for each sample"
else
	for i in $(ls *_haplotig.megablast.gz);do (
		if [[ ! -f "../sighits/sighits_phylum/${i%_haplotig.megablast.gz}_sighits.txt.gz" ]]; then
			if [[ "$taxids" == true ]]; then
			  for emg in ${projdir}/taxids/*.txids; do
					if [[ "$fullqlen_alignment" == "false" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-10 {print $0}' <(zcat "$i" | awk '$6>=95') <(zcat "$i" | awk '$6>=95') | awk '$3 >= 32 {print $0}' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_phylum/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
					if [[ "$fullqlen_alignment" == "true" ]]; then
						awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
						next} $3 >= max[$1]-10 {print $0}' <(zcat "$i" | awk '$6>=95') <(zcat "$i" | awk '$6>=95') | awk '$3 >= ($5-($5*0.05))' | \
						awk 'NR==FNR {a[$1]++; next} $9 in a' $emg - | awk 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
						awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | awk '{gsub(" ","\t",$0);}1' | gzip >> ../sighits/sighits_phylum/${i%_haplotig.megablast.gz}_sighits.txt.gz
						wait
					fi
			  done
				zcat ../sighits/sighits_phylum/${i%_haplotig.megablast.gz}_sighits.txt.gz | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | gzip > ../sighits/sighits_phylum/${i%_haplotig.megablast.gz}_sighits.tmp.gz
				mv ../sighits/sighits_phylum/${i%_haplotig.megablast.gz}_sighits.tmp.gz ../sighits/sighits_phylum/${i%_haplotig.megablast.gz}_sighits.txt.gz
				wait
			else
				if [[ "$fullqlen_alignment" == "false" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-10 {print $0}' <(zcat "$i" | awk '$6>=95') <(zcat "$i" | awk '$6>=95') | awk '$3 >= 32 {print $0}' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_phylum/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
				if [[ "$fullqlen_alignment" == "true" ]]; then
					awk 'NR == FNR {if (FNR == 1 || $3 > max[$1]) max[$1] = $3
					next} $3 >= max[$1]-10 {print $0}' <(zcat "$i" | awk '$6>=95') <(zcat "$i" | awk '$6>=95') | awk '$3 >= ($5-($5*0.05))' | awk 'gsub(" ","_",$0)1' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
					awk '{print $2,$3,$5,$6,$7,$8,$9,$10,$11,$1}' | cat <(printf "abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\n") - | \
					awk '{gsub(" ","\t",$0);}1' | gzip > ../sighits/sighits_phylum/${i%_haplotig.megablast.gz}_sighits.txt.gz
					wait
				fi
			fi
		fi
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_phylum
	for i in $(ls *_sighits.txt.gz);do (
		zcat "$i" | awk -F '\t' 'NR>1{a[$10]++;b[$10]=$0}END{for(x in a)if(a[x]==1)print b[x]}' > ${i%_sighits*}_phylum_unique_reads.txt
		$gzip ${i%_sighits*}_phylum_unique_reads.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait

  echo -e "${YELLOW}- compiling phylum-level multi-alignment algorithm"

  cd "${projdir}"/metagenome/sighits/sighits_phylum
  		for i in $(ls *_sighits.txt.gz);do (
				gunzip "$i"
				gunzip ${i%_sighits*}_phylum_unique_reads.txt.gz
				awk -F '\t' 'FNR==NR{a[$10]=1; next} {print $0}' ${i%_sighits*}_phylum_unique_reads.txt <(awk 'NR>1{print $0}' ${i%.gz}) OFS='\t' > ${i%_sighits*}_dup.txt
				rm "${i%.gz}"
				$gzip ${i%_sighits*}_phylum_unique_reads.txt
				$gzip ${i%_sighits*}_dup.txt
				) &
				if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
					wait
				fi
  		done
  		wait


  cd "${projdir}"/metagenome/sighits/sighits_phylum
  for i in $(ls *_dup.txt.gz);do (
  	awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_dup*}_taxids_dup.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

  cd "${projdir}"/metagenome/sighits/sighits_phylum

  for i in $(ls *_taxids_dup.txt);do
  	awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		cut -f 2-11 -d$'\t' > ${i%_taxids_dup*}_dup_inter.txt
  done
  wait

	for i in $(ls *_dup_inter.txt);do (
	   awk -F '\t'  '{print $8, $9, $10}' OFS='\t' "$i" > ${i%_dup_inter*}_phylum_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait
	rm *_taxids_dup.txt

  for i in $(ls *_dup.txt.gz);do (
    paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <(zcat "$i") ) <(awk -F '\t' '{print $1, $2, $3}' OFS='\t' ${i%*_dup.txt.gz}_phylum_taxa.txt) | \
		awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' OFS='\t' > ${i%_dup*}_phylum_duplicates.txt
		$gzip ${i%_dup*}_phylum_duplicates.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait
	for i in $(ls *_phylum_duplicates.txt.gz);do (
		 awk -F '\t' '$11 != "NA"' OFS='\t' <(zcat "$i") > ${i%*_phylum_duplicates.txt.gz}temp.txt
		 $gzip ${i%*_phylum_duplicates.txt.gz}temp.txt
 		mv ${i%*_phylum_duplicates.txt.gz}temp.txt.gz "$i"
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
	wait
  rm *_dup_inter.txt *_dup.txt.gz *_phylum_taxa.txt

  for i in $(ls *_phylum_duplicates.txt.gz);do (
    awk -F '\t' '{print $1, $10"~"$11, $2, $3, $4, $5, $6, $7, $8, $9, $12, $13}' OFS='\t' <(zcat "$i") > ${i%_phylum_duplicates*}_phylum_inter.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

  wait
  for i in $(ls *_phylum_inter.txt);do (
    awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' "$i" > ${i%_phylum_inter*}_phylum_inter2.txt ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  wait

  for i in $(ls *_phylum_inter2.txt);do (
    awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' "$i" | sort -T "${projdir}"/tmp -k1,1  > ${i%_phylum_inter2*}_duplicate_count.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_duplicate_count.txt);do (
    awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' "$i" | sort -T "${projdir}"/tmp -k1,1 > ${i%_duplicate_count*}_multialign_phylum_reads.txt) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  for i in $(ls *_multialign_phylum_reads.txt);do (
    awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' "$i" ${i%_multialign_phylum_reads*}_phylum_inter2.txt | sort -T "${projdir}"/tmp -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $4, $5, $6, $7, $8, $9, $10, $11, $1, $2, $12, $13}' OFS='\t' > ${i%_multialign_phylum_reads*}_phylum_OTU.txt
	  ) &
    if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
      wait
    fi
  done
	wait

  wait
  for i in $(ls *_phylum_unique_reads.txt.gz);do (
    awk -F '\t' '{print $8";"}' <(zcat "${i}") | awk -F ';' '{print $1}' > ${i%_phylum_unique_reads*}_taxids_uniq.txt
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
  done
	wait

	cd "${projdir}"/metagenome/sighits/sighits_phylum
	for i in $(ls *_taxids_uniq.txt); do
	   awk -F'\t' 'NR==FNR{a[$1]=$0;next} {print $0, ($1 in a ? a[$1]:"NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA")}' OFS='\t' ${projdir}/rankedlineage_edited.dmp $i | \
		 cut -f 2-11 -d$'\t' > ${i%_taxids_uniq*}_uniq_inter.txt
	done
	wait

  for i in $(ls *_uniq_inter.txt);do (
     awk -F '\t'  '{print $8, $9, $10}' OFS='\t' "$i" > ${i%_uniq_inter*}_phylum_taxa.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait

  rm *_taxids_uniq.txt


  for i in $(ls *_unique_reads.txt.gz);do (
     paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' <( zcat "$i") ) <(awk -F '\t' '{print $1, $2, $3}' OFS='\t' ${i%*_phylum_unique_reads*}_phylum_taxa.txt) | \
		 awk -F'\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13}' OFS='\t' > ${i%_phylum_uniq*}_phylum_unique_uncultured.txt &&
		 wait ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait
	for i in $(ls *_phylum_unique_uncultured.txt);do (
		awk -F '\t' '$11 != "NA"' OFS='\t' "$i" > "${i%*_phylum_unique_uncultured.txt}"temp.txt &&
		:> "$i" &&
		mv "${i%*_phylum_unique_uncultured.txt}"temp.txt "$i" 2> /dev/null
		) &
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	 done
  rm *_uniq_inter.txt && rm *_phylum_taxa.txt 2> /dev/null
  for i in $(ls *_phylum_unique_uncultured.txt);do (
     mv "$i" ${i%*_phylum_unique_uncultured*}_unique_sequences.txt 2> /dev/null
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
  done
	wait

	rm *_phylum_unique_uncultured.txt *_phylum_inter.txt *_phylum_inter2.txt *_duplicate_count.txt *_multialign_phylum_reads.txt *_phylum_duplicates.txt.gz 2> /dev/null

	for i in $(ls *_phylum_OTU.txt);do (
	   cat "$i" ${i%_phylum_OTU*}_unique_sequences.txt > ${i%_phylum_OTU*}_complete_phylum_reads.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	echo -e "${YELLOW}- compiling taxonomic information"
	for i in $(ls *_complete_phylum_reads.txt);do (
	   awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' "$i" > ${i%_complete_phylum_reads*}_sighits_temp.txt
	    echo $'abundance\tsseqid\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\tphylum\tkingdom\tdomain' | \
	     cat - ${i%_complete_phylum_reads*}_sighits_temp.txt > ${i%_complete_phylum_reads*}_sighits_temp2.txt
			 ) &
			 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
				 wait
			 fi
	done
	wait

	for i in $(ls *_sighits_temp2.txt);do (
	   awk -F '\t' '{gsub(/ /,"_");print}' "$i" > ${i%_sighits_temp2*}_sighits.txt
     $gzip ${i%_sighits_temp2*}_sighits.txt
		 ) &
		 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			 wait
		 fi
	done
	wait

	rm *_complete_phylum_reads.txt *_sighits_temp.txt *_unique_reads.txt.gz *_unique_sequences.txt *_sighits_temp2.txt *_phylum_OTU.txt
fi


cd "${projdir}"/metagenome/sighits/sighits_phylum
find . -type f -name '*_sighits.txt.gz' -exec cat {} + > sighits.txt.gz
awk '{print $8}' <(zcat sighits.txt.gz) | awk '{gsub(";","\n"); print}' | sort -T "${projdir}"/tmp -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt.gz
awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  ${projdir}/rankedlineage_edited.dmp taxids_sighits.txt > rankedlineage_subhits_temp.txt
awk '{gsub(/ /,"_");}1' rankedlineage_subhits_temp.txt | awk '{print $8, $9, $10}' | awk 'NR == 1; NR > 1 {print $0}' | awk '{gsub(/ /,"\t");}1' > rankedlineage_subhits.txt
rm taxids_sighits.txt rankedlineage_subhits_temp.txt

cd "${projdir}"/metagenome/sighits/sighits_phylum/
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/phylum/' > phylum_taxa_mean_temp1.txt
echo -e 'phylum' | cat - phylum_taxa_mean_temp1.txt > phylum_taxa_mean_temp.txt && rm phylum_taxa_mean_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/phylum/'  > phylum_taxa_unique_sequences_temp1.txt
echo -e 'phylum' | cat - phylum_taxa_unique_sequences_temp1.txt > phylum_taxa_unique_sequences_temp.txt && rm phylum_taxa_unique_sequences_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt |  sort -T "${projdir}"/tmp -u | awk '!/phylum/'  > phylum_taxa_quantification_accuracy_temp1.txt
echo -e 'phylum' | cat - phylum_taxa_quantification_accuracy_temp1.txt > phylum_taxa_quantification_accuracy_temp.txt && rm phylum_taxa_quantification_accuracy_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -T "${projdir}"/tmp -u | awk '!/phylum/'  > phylum_taxa_rel_quantification_accuracy_temp1.txt
echo -e 'phylum' | cat - phylum_taxa_rel_quantification_accuracy_temp1.txt > phylum_taxa_rel_quantification_accuracy_temp.txt && rm phylum_taxa_rel_quantification_accuracy_temp1.txt

phylum_level=phylum
for i in $(ls *_sighits.txt.gz);do
	gunzip "$i"
 	Rscript ${Qmatey_dir}/scripts/stats_summary.R ${i%.gz} $phylum_level "${Qmatey_dir}/tools/R"
	wait
	$gzip ${i%.gz}
  echo $'phylum\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
  id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdmean.txt phylum_taxa_mean_temp.txt > holdmean2.txt
  cat holdmean2.txt > phylum_taxa_mean_temp.txt
  id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holduniq_reads.txt phylum_taxa_unique_sequences_temp.txt > holduniq_reads2.txt
  cat holduniq_reads2.txt > phylum_taxa_unique_sequences_temp.txt
  id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt phylum_taxa_quantification_accuracy_temp.txt > holdstderr2.txt
  cat holdstderr2.txt > phylum_taxa_quantification_accuracy_temp.txt
  id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
  awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}' holdrelstderr.txt phylum_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt
  cat holdrelstderr2.txt > phylum_taxa_rel_quantification_accuracy_temp.txt
	wait
  rm *stats1* *stats2* *hold*
done
wait

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' phylum_taxa_mean_temp.txt > phylum_taxa_mean_temp2.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_mean_temp2.txt > ../../results/phylum_level/phylum_taxa_mean.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' phylum_taxa_unique_sequences_temp.txt > phylum_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_unique_sequences_temp2.txt > ../../results/phylum_level/phylum_taxa_unique_sequences.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_quantification_accuracy_temp.txt > phylum_taxa_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/phylum_level/phylum_taxa_mean.txt phylum_taxa_quantification_accuracy_temp2.txt > ../../results/phylum_level/phylum_taxa_quantification_accuracy.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_rel_quantification_accuracy_temp.txt > phylum_taxa_rel_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/phylum_level/phylum_taxa_mean.txt phylum_taxa_rel_quantification_accuracy_temp2.txt > ../../results/phylum_level/phylum_taxa_rel_quantification_accuracy.txt
rm *_temp*

cd "${projdir}"/metagenome/results/phylum_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do
  awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_phylum/rankedlineage_subhits.txt phylum_taxa_"${i}".txt | \
  awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > phylum_taxainfo_"${i}".txt &&
  wait
done
wait
rm *_taxa_*


awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
phylum_taxainfo_mean.txt | awk 'BEGIN{OFS="\t"} {$1=$1};1' > phylum_taxainfo_mean_holdingtaxinfo.txt
touch phylum_taxainfo_mean_buildnorm.txt
for i in ../../../metagenome/haplotig/*_metagenome.fasta.gz; do
  sample=${i%_metagenome.fasta.gz}; sample=${sample##*/}
  if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
    normfactor=1
  else
    normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
  fi
  awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" phylum_taxainfo_mean.txt | \
  paste - phylum_taxainfo_mean_buildnorm.txt > phylum_taxainfo_mean_buildnorm0.txt &&
  mv phylum_taxainfo_mean_buildnorm0.txt phylum_taxainfo_mean_buildnorm.txt
done
wait
paste phylum_taxainfo_mean_buildnorm.txt > phylum_taxainfo_mean_norm0.txt
paste phylum_taxainfo_mean_norm0.txt phylum_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > phylum_taxainfo_mean_normalized.txt
# pos=$(awk -v RS='\t' '/phylum/{print NR; exit}' phylum_taxainfo_mean_normalized.txt)
# awk -v pat="$pos" -F'\t'  'BEGIN {OFS="\t"}; {k=$pat; $pat=""; print k,$0}' phylum_taxainfo_mean_normalized.txt | awk 'BEGIN{FS="\t+"; OFS="\t"} {$1=$1; print}' > phylum_taxainfo_mean_normalized.tmp
# mv phylum_taxainfo_mean_normalized.tmp phylum_taxainfo_mean_normalized.txt
rm phylum_taxainfo_mean_buildnorm.txt phylum_taxainfo_mean_holdingtaxinfo.txt phylum_taxainfo_mean_norm0.txt 2> /dev/null

for i in $(ls *.txt); do
  awk '{gsub(/-/,"_"); print}' $i > ${i%.txt}.temp &&
  mv ${i%.txt}.temp $i
done
wait


# if [[ "$genome_scaling" == true ]]; then
# 	if [[ "$subsample_shotgun_R1" == false ]]; then
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" phylum "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "true" &>/dev/null
# 	else
# 		Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" phylum "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated "false" &>/dev/null
# 	fi
# else
# 		Rscript "${Qmatey_dir}/scripts/error_zero_inflation.R" phylum "${Qmatey_dir}/tools/R" $min_strain_uniq_ematch $zero_inflated &>/dev/null
# fi

file=${projdir}/exclude_taxa.txt
# if test -f $file; then
#   cat phylum_taxainfo_mean.txt > phylum_taxainfo_mean_filtered.txt &&
#   cat phylum_taxainfo_unique_sequences.txt > phylum_taxainfo_unique_sequences_filtered.txt &&
#   cat phylum_taxainfo_quantification_accuracy.txt > phylum_taxainfo_quantification_accuracy_filtered.txt &&
#   cat phylum_taxainfo_rel_quantification_accuracy.txt > phylum_taxainfo_rel_quantification_accuracy_filtered.txt &&
#   while read -r line; do
#     for i in $( ls *filtered.txt ); do
#       awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i
#     done
# 		wait
#   done < $file
# 	wait
# fi

if test -f $file; then
	echo -e "${YELLOW}- creating phylum-level visualizations"
	cd "${projdir}"/metagenome/results/phylum_level
	phylum_level_mean=phylum_taxainfo_mean_filtered.txt
	phylum_level_uniq=phylum_taxainfo_unique_sequences_filtered.txt
	phylum_level_stderr=phylum_taxainfo_quantification_accuracy_filtered.txt
	phylum_level_rel_stderr=phylum_taxainfo_rel_quantification_accuracy_filtered.txt

	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi


	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/phylum_level_boxplots.R" "$phylum_level_mean" "$phylum_level_mean_norm" "$phylum_level_uniq" "$phylum_level_stderr" "$phylum_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/phylum_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/phylum_level/boxplots/ 2> /dev/null

else
	echo -e "${YELLOW}- creating phylum-level visualizations"
	cd "${projdir}"/metagenome/results/phylum_level
	phylum_level_mean=phylum_taxainfo_mean.txt
	phylum_level_mean_norm=phylum_taxainfo_mean_normalized.txt
	phylum_level_uniq=phylum_taxainfo_unique_sequences.txt
	phylum_level_stderr=phylum_taxainfo_quantification_accuracy.txt
	phylum_level_rel_stderr=phylum_taxainfo_rel_quantification_accuracy.txt
	if [[ -z $phylum_level_mean_norm ]]; then
		:
	else
		phylum_level_mean=phylum_taxainfo_mean_normalized.txt
	fi
	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	for min_perc in ${min_percent_sample//,/ }; do (
		Rscript "${Qmatey_dir}/scripts/phylum_level_boxplots.R" "$phylum_level_mean" "$phylum_level_uniq" "$phylum_level_stderr" "$phylum_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	wait
	mkdir -p boxplots
	mv *_files ${projdir}/metagenome/results/phylum_level/boxplots/ 2> /dev/null
	mv *.html ${projdir}/metagenome/results/phylum_level/boxplots/ 2> /dev/null
fi

}
if [[ "$phylum_level" == "true" ]] && [[ -z "$(ls -A ${projdir}/metagenome/results/phylum_level/phylum_taxainfo* 2> /dev/null)" ]]; then
	if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_phylum_level/phylum_taxainfo* 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ${projdir}/metagenome/results/uncultured_phylum_level/phylum_taxainfo* 2> /dev/null)" ]]; then
			cd "${projdir}"/metagenome/alignment/
			mv ./uncultured/uncultured*megablast.gz ./ &&
			mv ./uncultured_combined_compressed.megablast.gz ./uncultured/ &&
			wait
			for i in $(ls *megablast.gz); do mv "$i" ${i#uncultured_}; done
			wait
			cd "${projdir}"
			time phylum_level 2>> ${projdir}/log.out
			wait
			cd "${projdir}"/metagenome/alignment/
			for i in $(ls *megablast.gz); do mv "$i" ./uncultured/uncultured_"${i}"; done
			wait
			mv "${projdir}"/metagenome/results/phylum_level ${projdir}/metagenome/results/uncultured_phylum_level
			mv "${projdir}"/metagenome/sighits/sighits_phylum ${projdir}/metagenome/sighits/uncultured_sighits_phylum
		fi
	fi
	if [[ -z "$(ls -A ${projdir}/metagenome/results/phylum_level/phylum_taxainfo* 2> /dev/null)" ]]; then
		cd "${projdir}"/metagenome/alignment/cultured/
		mv *megablast* ../
		mv ../combined_compressed.megablast.gz ./
		cd "${projdir}"
		time phylum_level 2>> ${projdir}/log.out
		mv "${projdir}"/metagenome/alignment/*megablast* ${projdir}/metagenome/alignment/cultured/
	fi
else
	echo -e "${magenta}- skipping exact matching for phylum-level profiling ${white}\n"
fi


##########################################################################
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing cross-taxon validation and filtering \n\e[97m########################################################\n"
cross_taxon_validate() {
	#############################
	echo -e "${YELLOW}- performing cross-taxon validation and filtering at class level"
	cd "${projdir}"/metagenome/results/class_level
	mkdir -p validated
	Rscript "${Qmatey_dir}/scripts/cross_taxon_validation.R" "../phylum_level/phylum" "class" "phylum" "${Qmatey_dir}/tools/R" "./validated/class" &&
	mv ./validated ../class_level_validated
	wait
	#############################
	echo -e "${YELLOW}- performing cross-taxon validation and filtering at order level"
	cd "${projdir}"/metagenome/results/order_level
	mkdir -p validated
	Rscript "${Qmatey_dir}/scripts/cross_taxon_validation.R" "../class_level_validated/class" "order" "class" "${Qmatey_dir}/tools/R" "./validated/order" &&
	mv ./validated ../order_level_validated
	wait
	#############################
	echo -e "${YELLOW}- performing cross-taxon validation and filtering at family level"
	cd "${projdir}"/metagenome/results/family_level
	mkdir -p validated
	Rscript "${Qmatey_dir}/scripts/cross_taxon_validation.R" "../order_level_validated/order" "family" "order" "${Qmatey_dir}/tools/R" "./validated/family" &&
	mv ./validated ../family_level_validated
	wait
	#############################
	echo -e "${YELLOW}- performing cross-taxon validation and filtering at genus level"
	cd "${projdir}"/metagenome/results/genus_level
	mkdir -p validated
	Rscript "${Qmatey_dir}/scripts/cross_taxon_validation.R" "../family_level_validated/family" "genus" "family" "${Qmatey_dir}/tools/R" "./validated/genus" &&
	mv ./validated ../genus_level_validated
	wait
	#############################
	echo -e "${YELLOW}- performing cross-taxon validation and filtering at species level"
	cd "${projdir}"/metagenome/results/species_level
	mkdir -p validated
	Rscript "${Qmatey_dir}/scripts/cross_taxon_validation.R" "../genus_level_validated/genus" "species" "genus" "${Qmatey_dir}/tools/R" "./validated/species" &&
	mv ./validated ../species_level_validated
	wait
	#############################
	echo -e "${YELLOW}- performing cross-taxon validation and filtering at strain level (mininimum unique seqeuence = 1)"
	cd "${projdir}"/metagenome/results/strain_level_minUniq_1
	mkdir -p validated
	Rscript "${Qmatey_dir}/scripts/cross_taxon_validation.R" "../species_level_validated/species" "strain" "species" "${Qmatey_dir}/tools/R" "./validated/strain" &&
	mv ./validated ../strain_level_minUniq_1_validated
	wait
	#############################
	# echo -e "${YELLOW}- performing cross-taxon validation and filtering at strain level (mininimum unique seqeuence = 2)"
	# cd "${projdir}"/metagenome/results/strain_level_minUniq_2
	# mkdir -p validated
	# Rscript "${Qmatey_dir}/scripts/cross_taxon_validation.R" "../species_level_validated/species" "strain" "species" "${Qmatey_dir}/tools/R" "./validated/strain" &&
	# mv ./validated ../strain_level_minUniq_2_validated
	# wait
	#############################

	cd "${projdir}"/metagenome/results/
	Rscript "${Qmatey_dir}/scripts/ctv_visualization.R" "${Qmatey_dir}/tools/R" &&
	wait

}
cd "${projdir}"
if [[ "$cross_taxon_validation" == true ]]; then
	time cross_taxon_validate &>> ${projdir}/log.out
else
	echo -e "${magenta}- skipping exact matching for phylum-level profiling ${white}\n"
fi


######################################################################################################################################################
sunburst() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is creating sunburst \n\e[97m########################################################\n"

cd "${projdir}"/metagenome/results

for tsun in ${sunburst_taxlevel//,/ }; do
	echo -e "${YELLOW}- creating sunburst from ${tsun}_level taxonomic profile"
	cd "${projdir}"/metagenome/results

	if [[ "$tsun" == strain ]]; then
		for strain_minUniq in *strain_level_minUniq_*/; do
		  cd ${strain_minUniq}
		  mean=${tsun}_taxainfo_mean.txt
		  mean_norm=${tsun}_taxainfo_mean_normalized.txt
		  if [[ -z $mean_norm ]]; then
		    mean_norm=$mean
		  fi

			if [[ -z $min_percent_sample ]]; then
				min_percent_sample=5,10,20
			fi

		  for min_perc in ${min_percent_sample//,/ }; do (
		    Rscript "${Qmatey_dir}/scripts/sunburst.R" "$mean_norm" "$min_perc" "${sunburst_nlayers}" "${Qmatey_dir}/tools/R" $tsun 2>/dev/null )&
				if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				  wait
				fi
		  done
			cd ../
			wait
		done
		wait
	else
		cd "${projdir}"/metagenome/results
		for culture_type in *${tsun}_level*/; do
			cd $culture_type
			if [[ ! -d "sunburst" ]] ;then
				mean=${tsun}_taxainfo_mean.txt
				mean_norm=${tsun}_taxainfo_mean_normalized.txt
				if [[ -z $mean_norm ]]; then
					mean_norm=$mean
				fi

				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$tsun" == species ]] ; then
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/sunburst.R" "$mean_norm" "$min_perc" "${sunburst_nlayers}" "${Qmatey_dir}/tools/R" $tsun 2>/dev/null )&
						if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							wait
						fi
					done
					wait
				fi
				if [[ "$tsun" == genus ]] || [[ "$tsun" == family ]] || [[ "$tsun" == order ]] || [[ "$tsun" == class ]]; then
					sunburst_nlayers2=phylum,$tsun
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/sunburst.R" "$mean_norm" "$min_perc" "${sunburst_nlayers2}" "${Qmatey_dir}/tools/R" $tsun 2>/dev/null )&
						if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							wait
						fi
					done
					wait
				fi
				if [[ "$tsun" == phylum ]]; then
					sunburst_nlayers2=phylum
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/sunburst.R" "$mean_norm" "$min_perc" "${sunburst_nlayers2}" "${Qmatey_dir}/tools/R" $tsun 2>/dev/null )&
						if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							wait
						fi
					done
					wait
				fi
				cd ../
			fi
		done
		wait
	fi
	wait
done
wait

}
cd "${projdir}"
if [[ "$run_sunburst" == true ]]; then
	time sunburst &>> ${projdir}/log.out
fi


######################################################################################################################################################
correlogram() {
	echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is performing correlation of compositional data and creating correlogams \n\e[97m########################################################\n"
	file=${projdir}/exclude_taxa.txt


	#### strain: compositionality-corrected p-values, q-values, and Z-scores for all pairwise correlations
	######################################################################################################
	cd "${projdir}"/metagenome/results
	for strain_minUniq in strain_level_minUniq_*/; do
		if test -f $file; then
			echo -e "${YELLOW}- creating strain-level visualizations"
			cd "${projdir}"/metagenome/results/"${strain_minUniq}"
			if test ! -f *corr.txt; then
				strain_level_mean=strain_taxainfo_mean_filtered.txt
				strain_level_uniq=strain_taxainfo_unique_sequences_filtered.txt
				strain_level_stderr=strain_taxainfo_quantification_accuracy_filtered.txt
				strain_level_rel_stderr=strain_taxainfo_rel_quantification_accuracy_filtered.txt

				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && test -f $strain_level_mean; then
					if [[ "$CCLasso_corr" =~  strain || "$spearman_corr" =~  strain ]]; then
						if [[ "$CCLasso_corr" =~  strain ]] ; then corr_CCLasso=true; fi
						if [[ "$spearman_corr" =~  strain ]] ; then corr_spearman=true; fi
						for min_perc in ${min_percent_sample//,/ }; do (
							Rscript "${Qmatey_dir}/scripts/strain_level_corr.R" "$strain_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
							)&
						 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							 wait
						 fi
						done
						wait
					fi
				fi
			fi
		else
			echo -e "${YELLOW}- creating strain-level visualizations"
			cd "${projdir}"/metagenome/results/"${strain_minUniq}"
			if test ! -f *corr.txt; then
				strain_level_mean=strain_taxainfo_mean.txt
				strain_level_mean_norm=strain_taxainfo_mean_normalized.txt
				strain_level_uniq=strain_taxainfo_unique_sequences.txt
				strain_level_stderr=strain_taxainfo_quantification_accuracy.txt
				strain_level_rel_stderr=strain_taxainfo_rel_quantification_accuracy.txt
				if [[ -z $strain_level_mean_norm ]]; then
					:
				else
					strain_level_mean=strain_taxainfo_mean_normalized.txt
				fi
				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && [[ "$CCLasso_corr" =~ strain ]] && test -f $strain_level_mean; then
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/strain_level_corr.R" "$strain_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
						)&
					 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
						 wait
					 fi
					done
					wait
				fi
			fi
		fi
		wait
		mkdir -p CCLasso_correlation
		mv *corr.tiff ./CCLasso_correlation/ 2> /dev/null
		cd "${projdir}"/metagenome/results/
	done
	wait

	#### species: compositionality-corrected p-values, q-values, and Z-scores for all pairwise correlations
	######################################################################################################
	cd "${projdir}"/metagenome/results
	for speciesdir in species_level*/; do
		if test -f $file; then
			echo -e "${YELLOW}- creating species-level visualizations"
			cd "${projdir}"/metagenome/results/$speciesdir
			if test ! -f *corr.txt; then
				species_level_mean=species_taxainfo_mean_filtered.txt
				species_level_uniq=species_taxainfo_unique_sequences_filtered.txt
				species_level_stderr=species_taxainfo_quantification_accuracy_filtered.txt
				species_level_rel_stderr=species_taxainfo_rel_quantification_accuracy_filtered.txt

				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && test -f $species_level_mean; then
					if [[ "$CCLasso_corr" =~  species || "$spearman_corr" =~  species ]]; then
						if [[ "$CCLasso_corr" =~  species ]] ; then corr_CCLasso=true; fi
						if [[ "$spearman_corr" =~  species ]] ; then corr_spearman=true; fi
						for min_perc in ${min_percent_sample//,/ }; do (
							Rscript "${Qmatey_dir}/scripts/species_level_corr.R" "$species_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
							)&
						 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							 wait
						 fi
						done
						wait
					fi
				fi
			fi
		else
			echo -e "${YELLOW}- creating species-level visualizations"
			cd "${projdir}"/metagenome/results/$speciesdir
			if test ! -f *corr.txt; then
				species_level_mean=species_taxainfo_mean.txt
				species_level_mean_norm=species_taxainfo_mean_normalized.txt
				species_level_uniq=species_taxainfo_unique_sequences.txt
				species_level_stderr=species_taxainfo_quantification_accuracy.txt
				species_level_rel_stderr=species_taxainfo_rel_quantification_accuracy.txt
				if [[ -z $species_level_mean_norm ]]; then
					:
				else
					species_level_mean=species_taxainfo_mean_normalized.txt
				fi
				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && [[ "$CCLasso_corr" =~ species ]] && test -f $species_level_mean; then
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/species_level_corr.R" "$species_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
						)&
					 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
						 wait
					 fi
					done
					wait
				fi
			fi
		fi
		wait
		mkdir -p CCLasso_correlation
		mv *corr.tiff ./CCLasso_correlation/ 2> /dev/null
		cd "${projdir}"/metagenome/results
	done
	wait

	#### genus: compositionality-corrected p-values, q-values, and Z-scores for all pairwise correlations
	######################################################################################################
	cd "${projdir}"/metagenome/results
	for genusdir in genus_level*; do
		if test -f $file; then
			echo -e "${YELLOW}- creating genus-level visualizations"
			cd "${projdir}"/metagenome/results/$genusdir
			if test ! -f *corr.txt; then
				genus_level_mean=genus_taxainfo_mean_filtered.txt
				genus_level_uniq=genus_taxainfo_unique_sequences_filtered.txt
				genus_level_stderr=genus_taxainfo_quantification_accuracy_filtered.txt
				genus_level_rel_stderr=genus_taxainfo_rel_quantification_accuracy_filtered.txt

				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && test -f $genus_level_mean; then
					if [[ "$CCLasso_corr" =~  genus || "$spearman_corr" =~  genus ]]; then
						if [[ "$CCLasso_corr" =~  genus ]] ; then corr_CCLasso=true; fi
						if [[ "$spearman_corr" =~  genus ]] ; then corr_spearman=true; fi
						for min_perc in ${min_percent_sample//,/ }; do (
							Rscript "${Qmatey_dir}/scripts/genus_level_corr.R" "$genus_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
							)&
						 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							 wait
						 fi
						done
						wait
					fi
				fi
			fi
		else
			echo -e "${YELLOW}- creating genus-level visualizations"
			cd "${projdir}"/metagenome/results/$genusdir
			if test ! -f *corr.txt; then
				genus_level_mean=genus_taxainfo_mean.txt
				genus_level_mean_norm=genus_taxainfo_mean_normalized.txt
				genus_level_uniq=genus_taxainfo_unique_sequences.txt
				genus_level_stderr=genus_taxainfo_quantification_accuracy.txt
				genus_level_rel_stderr=genus_taxainfo_rel_quantification_accuracy.txt
				if [[ -z $genus_level_mean_norm ]]; then
					:
				else
					genus_level_mean=genus_taxainfo_mean_normalized.txt
				fi
				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && [[ "$CCLasso_corr" =~ genus ]] && test -f $genus_level_mean; then
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/genus_level_corr.R" "$genus_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
						)&
					 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
						 wait
					 fi
					done
					wait
				fi
			fi
		fi
		wait
		mkdir -p CCLasso_correlation
		mv *corr.tiff ./CCLasso_correlation/ 2> /dev/null
		cd "${projdir}"/metagenome/results
	done
	wait


	#### family: compositionality-corrected p-values, q-values, and Z-scores for all pairwise correlations
	######################################################################################################
	cd "${projdir}"/metagenome/results
	for familydir in family_level*/; do
		if test -f $file; then
			echo -e "${YELLOW}- creating family-level visualizations"
			cd "${projdir}"/metagenome/results/$familydir
			if test ! -f *corr.txt; then
				family_level_mean=family_taxainfo_mean_filtered.txt
				family_level_uniq=family_taxainfo_unique_sequences_filtered.txt
				family_level_stderr=family_taxainfo_quantification_accuracy_filtered.txt
				family_level_rel_stderr=family_taxainfo_rel_quantification_accuracy_filtered.txt

				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && test -f $family_level_mean; then
					if [[ "$CCLasso_corr" =~  family || "$spearman_corr" =~  family ]]; then
						if [[ "$CCLasso_corr" =~  family ]] ; then corr_CCLasso=true; fi
						if [[ "$spearman_corr" =~  family ]] ; then corr_spearman=true; fi
						for min_perc in ${min_percent_sample//,/ }; do (
							Rscript "${Qmatey_dir}/scripts/family_level_corr.R" "$family_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
							)&
						 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							 wait
						 fi
						done
						wait
					fi
				fi
			fi
		else
			echo -e "${YELLOW}- creating family-level visualizations"
			cd "${projdir}"/metagenome/results/$familydir
			if test ! -f *corr.txt; then
				family_level_mean=family_taxainfo_mean.txt
				family_level_mean_norm=family_taxainfo_mean_normalized.txt
				family_level_uniq=family_taxainfo_unique_sequences.txt
				family_level_stderr=family_taxainfo_quantification_accuracy.txt
				family_level_rel_stderr=family_taxainfo_rel_quantification_accuracy.txt
				if [[ -z $family_level_mean_norm ]]; then
					:
				else
					family_level_mean=family_taxainfo_mean_normalized.txt
				fi
				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && [[ "$CCLasso_corr" =~ family ]] && test -f $family_level_mean; then
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/family_level_corr.R" "$family_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
						)&
					 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
						 wait
					 fi
					done
					wait
				fi
			fi
		fi
		wait
		mkdir -p CCLasso_correlation
		mv *corr.tiff ./CCLasso_correlation/ 2> /dev/null
		cd "${projdir}"/metagenome/results
	done
	wait


	#### order: compositionality-corrected p-values, q-values, and Z-scores for all pairwise correlations
	######################################################################################################
	cd "${projdir}"/metagenome/results
	for orderdir in order_level*/; do
		if test -f $file; then
			echo -e "${YELLOW}- creating order-level visualizations"
			cd "${projdir}"/metagenome/results/$orderdir
			if test ! -f *corr.txt; then
				order_level_mean=order_taxainfo_mean_filtered.txt
				order_level_uniq=order_taxainfo_unique_sequences_filtered.txt
				order_level_stderr=order_taxainfo_quantification_accuracy_filtered.txt
				order_level_rel_stderr=order_taxainfo_rel_quantification_accuracy_filtered.txt

				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && test -f $order_level_mean; then
					if [[ "$CCLasso_corr" =~  order || "$spearman_corr" =~  order ]]; then
						if [[ "$CCLasso_corr" =~  order ]] ; then corr_CCLasso=true; fi
						if [[ "$spearman_corr" =~  order ]] ; then corr_spearman=true; fi
						for min_perc in ${min_percent_sample//,/ }; do (
							Rscript "${Qmatey_dir}/scripts/order_level_corr.R" "$order_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
							)&
						 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							 wait
						 fi
						done
						wait
					fi
				fi
			fi
		else
			echo -e "${YELLOW}- creating order-level visualizations"
			cd "${projdir}"/metagenome/results/$orderdir
			if test ! -f *corr.txt; then
				order_level_mean=order_taxainfo_mean.txt
				order_level_mean_norm=order_taxainfo_mean_normalized.txt
				order_level_uniq=order_taxainfo_unique_sequences.txt
				order_level_stderr=order_taxainfo_quantification_accuracy.txt
				order_level_rel_stderr=order_taxainfo_rel_quantification_accuracy.txt
				if [[ -z $order_level_mean_norm ]]; then
					:
				else
					order_level_mean=order_taxainfo_mean_normalized.txt
				fi
				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && [[ "$CCLasso_corr" =~ order ]] && test -f $order_level_mean; then
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/order_level_corr.R" "$order_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
						)&
					 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
						 wait
					 fi
					done
					wait
				fi
			fi
		fi
		wait
		mkdir -p CCLasso_correlation
		mv *corr.tiff ./CCLasso_correlation/ 2> /dev/null
		cd "${projdir}"/metagenome/results
	done
	wait


	#### class: compositionality-corrected p-values, q-values, and Z-scores for all pairwise correlations
	######################################################################################################
	cd "${projdir}"/metagenome/results
	for classdir in class_level*; do
		if test -f $file; then
			echo -e "${YELLOW}- creating class-level visualizations"
			cd "${projdir}"/metagenome/results/$classdir
			if test ! -f *corr.txt; then
				class_level_mean=class_taxainfo_mean_filtered.txt
				class_level_uniq=class_taxainfo_unique_sequences_filtered.txt
				class_level_stderr=class_taxainfo_quantification_accuracy_filtered.txt
				class_level_rel_stderr=class_taxainfo_rel_quantification_accuracy_filtered.txt

				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && [[ "$CCLasso_corr" =~ class ]] && test -f $class_level_mean; then
					if [[ "$CCLasso_corr" =~  class || "$spearman_corr" =~  class ]]; then
						if [[ "$CCLasso_corr" =~  class ]] ; then corr_CCLasso=true; fi
						if [[ "$spearman_corr" =~  class ]] ; then corr_spearman=true; fi
						for min_perc in ${min_percent_sample//,/ }; do (
							Rscript "${Qmatey_dir}/scripts/class_level_corr.R" "$class_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
							)&
						 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							 wait
						 fi
						done
						wait
					fi
				fi
			fi
		else
			echo -e "${YELLOW}- creating class-level visualizations"
			cd "${projdir}"/metagenome/results/$classdir
			if test ! -f *corr.txt; then
				class_level_mean=class_taxainfo_mean.txt
				class_level_mean_norm=class_taxainfo_mean_normalized.txt
				class_level_uniq=class_taxainfo_unique_sequences.txt
				class_level_stderr=class_taxainfo_quantification_accuracy.txt
				class_level_rel_stderr=class_taxainfo_rel_quantification_accuracy.txt
				if [[ -z $class_level_mean_norm ]]; then
					:
				else
					class_level_mean=class_taxainfo_mean_normalized.txt
				fi
				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && [[ "$CCLasso_corr" =~ class ]] && test -f $class_level_mean; then
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/class_level_corr.R" "$class_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
						)&
					 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
						 wait
					 fi
					done
					wait
				fi
			fi
		fi
		wait
		cd "${projdir}"/metagenome/results/class_level
		mkdir -p CCLasso_correlation
		mv *corr.tiff ./CCLasso_correlation/ 2> /dev/null
		cd "${projdir}"/metagenome/results
	done
	wait


	#### phylum: compositionality-corrected p-values, q-values, and Z-scores for all pairwise correlations
	######################################################################################################
	cd "${projdir}"/metagenome/results
	for phylumdir in phylum_level*/; do
		if test -f $file; then
			echo -e "${YELLOW}- creating phylum-level visualizations"
			cd "${projdir}"/metagenome/results/$phylumdir
			if test ! -f *corr.txt; then
				phylum_level_mean=phylum_taxainfo_mean_filtered.txt
				phylum_level_uniq=phylum_taxainfo_unique_sequences_filtered.txt
				phylum_level_stderr=phylum_taxainfo_quantification_accuracy_filtered.txt
				phylum_level_rel_stderr=phylum_taxainfo_rel_quantification_accuracy_filtered.txt

				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && [[ "$CCLasso_corr" =~ phylum ]] && test -f $phylum_level_mean; then
					if [[ "$CCLasso_corr" =~  phylum || "$spearman_corr" =~  phylum ]]; then
						if [[ "$CCLasso_corr" =~  phylum ]] ; then corr_CCLasso=true; fi
						if [[ "$spearman_corr" =~  phylum ]] ; then corr_spearman=true; fi
						for min_perc in ${min_percent_sample//,/ }; do (
							Rscript "${Qmatey_dir}/scripts/phylum_level_corr.R" "$phylum_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
							)&
						 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
							 wait
						 fi
						done
						wait
					fi
				fi
			fi
		else
			echo -e "${YELLOW}- creating phylum-level visualizations"
			cd "${projdir}"/metagenome/results/$phylumdir
			if test ! -f *corr.txt; then
				phylum_level_mean=phylum_taxainfo_mean.txt
				phylum_level_mean_norm=phylum_taxainfo_mean_normalized.txt
				phylum_level_uniq=phylum_taxainfo_unique_sequences.txt
				phylum_level_stderr=phylum_taxainfo_quantification_accuracy.txt
				phylum_level_rel_stderr=phylum_taxainfo_rel_quantification_accuracy.txt
				if [[ -z $phylum_level_mean_norm ]]; then
					:
				else
					phylum_level_mean=phylum_taxainfo_mean_normalized.txt
				fi
				if [[ -z $min_percent_sample ]]; then
					min_percent_sample=5,10,20
				fi

				if [[ "$total_no_samples" -ge 24 ]] && [[ "$CCLasso_corr" =~ phylum ]] && test -f $phylum_level_mean; then
					for min_perc in ${min_percent_sample//,/ }; do (
						Rscript "${Qmatey_dir}/scripts/phylum_level_corr.R" "$phylum_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "$corr_CCLasso" "$corr_spearman" "${Qmatey_dir}/tools/R" 2>/dev/null
						)&
					 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
						 wait
					 fi
					done
					wait
				fi
			fi
		fi
		wait
		cd "${projdir}"/metagenome/results/phylum_level
		mkdir -p CCLasso_correlation
		mv *corr.tiff ./CCLasso_correlation/ 2> /dev/null
		cd "${projdir}"/metagenome/results
	done
	wait
}
if [[ "$run_corr" == true || "$run_spearmancorr" == true ]]; then
	correlogram &>> ${projdir}/log.out
fi

######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${white}\n"
annotate() {
	cd "${projdir}"
	find . -depth -type d -exec rmdir {} + 2> /dev/null
	mkdir -p norm_ref

	rm "${projdir}"/rankedlineage_edited.dmp

	cd "${projdir}"
	mkdir -p ${projdir}/metagenome/results/results_uncultured
	mv "${projdir}"/metagenome/results/uncultured_* ${projdir}/metagenome/results/results_uncultured/
	mv "${projdir}"/metagenome/results/results_uncultured/ ${projdir}/metagenome/
	cd "${projdir}"/metagenome/

	if [[ "$(ls ./results/strain_level_minUniq_*_validated/strain_taxainfo_mean_normalized.txt 2> /dev/null | wc -l)" -gt 0 ]]; then
		mkdir gene_annotation_count
		taxid_genes=$(ls ./results/strain_level_minUniq_*_validated/strain_taxainfo_mean_normalized.txt | tail -n1)
		if test -f ./alignment/rRNA/rRNA_combined_compressed.megablast.gz; then
			zcat ./alignment/cultured/combined_compressed.megablast.gz | awk -F'\t' '{print $9"\t"$7"\t"$2"\t"$10}' | \
			cat - <(zcat ./alignment/rRNA/rRNA_combined_compressed.megablast.gz | awk -F'\t' '{print $9"\t"$7"\t"$2"\t"$10}') | gzip > All_compressed.megablast.gz
		else
			zcat ./alignment/cultured/combined_compressed.megablast.gz | awk -F'\t' '{print $9"\t"$7"\t"$2"\t"$10}' | gzip > All_compressed.megablast.gz
		fi
		awk 'NR>1{print $1}' $taxid_genes | sort -T "${projdir}"/tmp | uniq | grep -Fwf - <(zcat All_compressed.megablast.gz) | \
		awk -F'\t' '!seen[$1$3]++' | awk -F'\t' '!seen[$1$2]++' | \
		cat <(printf "tax_id\tsequence\tGenBank_ID\tgene_annotation\n") - | gzip > ./gene_annotation_count/combined_taxids_sequences_genes_geneID.txt.gz
		taxnamecol=$(head -n1 $taxid_genes | tr '\t' '\n' | cat -n | grep 'taxname' | awk '{print $1}')
		zcat ./gene_annotation_count/combined_taxids_sequences_genes_geneID.txt.gz | awk 'NR>1{print $1}' | sort -T "${projdir}"/tmp | uniq -c | awk '{$1=$1};1' | awk '{gsub(/ /,"\t");}1' | \
		awk -F'\t' 'NR==FNR {h[$2] = $1; next} {print $1,$2,h[$2]}' - <(awk -v taxname=$taxnamecol '{print $taxname"\t"$1}' $taxid_genes) | \
		awk '{gsub(/ /,"\t");}1' | awk 'NR>1{print $2"\t"$3"\t"$1}' | cat <(printf "tax_id\tgene_count\ttaxname\n") - > ./gene_annotation_count/combined_genes_per_taxid.txt
		rm All_compressed.megablast.gz

		for i in ./alignment/cultured/*haplotig.megablast.gz; do (
			taxid_genes=$(ls ./results/strain_level_minUniq_*/strain_taxainfo_mean_normalized.txt | tail -n1)
			awk 'NR>1{print $1}' $taxid_genes | sort -T "${projdir}"/tmp | uniq | grep -Fwf - \
			<(zcat "$i" | awk -F'\t' '{print $9"\t"$7"\t"$1"\t"$2"\t"$10}') | awk -F'\t' '!seen[$1$4]++' | awk -F'\t' '!seen[$1$2]++' | awk '{gsub(/-/,"\t",$3);}1' | awk '{$3=""}1' |\
			cat <(printf "tax_id\tsequence\tRelative_Abundance\tGenBank_ID\tgene_annotation\n") - | gzip > ${i%*haplotig.megablast.gz}taxids_sequences_genes_geneID.txt.gz
			taxnamecol=$(head -n1 $taxid_genes | tr '\t' '\n' | cat -n | grep 'taxname' | awk '{print $1}')
			zcat ${i%*haplotig.megablast.gz}taxids_sequences_genes_geneID.txt.gz | awk '{print $1}' | sort -T "${projdir}"/tmp | uniq -c | awk '{$1=$1};1' | awk '{gsub(/ /,"\t");}1' | \
			awk -F'\t' 'NR==FNR {h[$2] = $1; next} {print $1,$2,h[$2]}' - <(awk -v taxname=$taxnamecol '{print $taxname"\t"$1}' $taxid_genes) | \
			awk '{gsub(/ /,"\t");}1' | awk 'NR>1{print $2"\t"$3"\t"$1}' | cat <(printf "tax_id\tgene_count\ttaxname\n") - > ${i%*haplotig.megablast.gz}_genes_per_taxid.txt
			) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			 wait
			fi
		done
		wait
		for i in ./sighits/sighits_strain/*_sighits.txt.gz; do (
			taxid_genes=$(ls ./results/strain_level_minUniq_*/strain_taxainfo_mean_normalized.txt | tail -n1)
			awk 'NR>1{print $1}' $taxid_genes | sort -T "${projdir}"/tmp | uniq | grep -Fwf - \
			<(zcat "$i" | awk -F'\t' '{print $8"\t"$6"\t"$10"\t"$1"\t"$2"\t"$9}') | awk -F'\t' '!seen[$1$4]++' | awk -F'\t' '!seen[$1$2]++' | \
			cat <(printf "tax_id\tsequence\tseqid\tRelative_Abundance\tGenBank_ID\tgene_annotation\n") - | gzip > ${i%*_sighits.txt.gz}taxids_sequences_genes_geneID_Diagnostic.txt.gz
			taxnamecol=$(head -n1 $taxid_genes | tr '\t' '\n' | cat -n | grep 'taxname' | awk '{print $1}')
			zcat ${i%*_sighits.txt.gz}taxids_sequences_genes_geneID_Diagnostic.txt.gz | awk '{print $1}' | sort -T "${projdir}"/tmp | uniq -c | awk '{$1=$1};1' | awk '{gsub(/ /,"\t");}1' | \
			awk -F'\t' 'NR==FNR {h[$2] = $1; next} {print $1,$2,h[$2]}' - <(awk -v taxname=$taxnamecol '{print $taxname"\t"$1}' $taxid_genes) | \
			awk '{gsub(/ /,"\t");}1' | awk 'NR>1{print $2"\t"$3"\t"$1}' | cat <(printf "tax_id\tgene_count\ttaxname\n") - > ${i%*_sighits.txt.gz}_genes_per_taxid_Diagnostic.txt
			) &
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			 wait
			fi
		done
		wait
		mv ./alignment/cultured/*taxids_sequences_genes_geneID*.txt.gz ./gene_annotation_count/
		mv ./alignment/cultured/*genes_per_taxid*.txt ./gene_annotation_count/
		mv ./sighits/sighits_strain/*Diagnostic* ./gene_annotation_count/
		zcat ./gene_annotation_count/*taxids_sequences_genes_geneID_Diagnostic.txt.gz | grep -v '^tax_id' | cat <(printf "tax_id\tsequence\tseqid\tRelative_Abundance\tGenBank_ID\tgene_annotation\n") - | gzip > ./gene_annotation_count/combined_taxids_sequences_genes_geneID_Diagnostic.txt.gz

	fi
}
cd "${projdir}"
if [[ "$annotate_seq" == true ]]; then
	time annotate &>> ${projdir}/log.out
fi

if [[ "$normalization" == true ]]; then
	mv "${projdir}"/metagenome ${projdir}/metagenome_ref_normalize
fi
if [[ "$normalization" == false ]]; then
	mv "${projdir}"/metagenome ${projdir}/metagenome_no_normalize
fi

if [[ -z $samples_alt_dir ||  $samples_alt_dir == false ]]; then
 	:
else
	rm samples
fi
cd "${projdir}"
rm -rf "${projdir}"/metagenome*/sighits
rm rankedlineage_edited.dmp
touch Analysis_Complete
echo -e "\n\n${magenta}- Run Complete ${white}"
