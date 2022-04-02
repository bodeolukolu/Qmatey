if R --version; then
	:
else
	module add R
fi

#cd $projdir
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
	if [[ -z $input_dbfasta ]]; then
		echo -e "${magenta}- \n- Please provide fasta file(s) required to create custom database  ${white}\n"
		exit 0
	fi
	if [[ -z $map_taxids ]]; then
		echo -e "${magenta}- \n- Please provide files mapping fasta sequences to taxid(s)  ${white}\n"
		exit 0
	fi
	cd $projdir
	echo -e "$1 \e[31m Creating custom database"
	${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/makeblastdb -in $input_dbfasta -parse_seqids -blastdb_version 10 -taxid_map $map_taxids -title "custom_db" -dbtype nt
fi


cd $projdir
if [[ -z $min_unique_seqs ]]; then
	min_unique_seqs=1
fi
if [[ -z $maxindel ]]; then
	maxindel=100
fi
if [[ -z $min_percent_sample ]]; then
	min_percent_sample=5,10,20
fi
min_percent_sample=$( echo "$min_percent_sample" | awk '{gsub(/,/,"\n")}1' )
if [[ -z $min_pos_corr ]]; then
	min_pos_corr=0.1,0.2,0.3
fi
if [[ -z $max_neg_corr ]]; then
	max_neg_corr=0.1,0.2,0.3
fi
if [[ -z $genome_scaling ]]; then
	genome_scaling=0
fi
if [[ -z $sunburst_taxlevel ]]; then
	sunburst_taxlevel=strain
fi
if [[ -z $sunburst_nlayers ]]; then
	sunburst_nlayers=phylum,genus,species
fi




if [[ -z "$(ls -A ./samples)" ]]; then
	echo -e "$1 \e[31m samples folder is empty, Qmatey will exit"; sleep 10; exit 0
else
	:
fi

totalk=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)
loopthreads=2
if [[ "$threads" -gt 1 ]]; then
	N=$((threads/2))
	ram1=$(($totalk/$N))
else
	N=1 && loopthreads=threads
fi
ram1=$((ram1/1000000))
Xmx1=-Xmx${ram1}G
ram2=$(echo "$totalk*0.00000095" | bc)
ram2=${ram2%.*}
Xmx2=-Xmx${ram2}G

if [[ -z "$threads" ]]; then
	threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		threads=$((threads-2))
	fi
fi
if  [[ "$threads" -ge 1 ]]; then
	loopthread=2
	N=$(($threads/2))
else
	N=1 && loopthread=$threads
fi

cd ${Qmatey_dir}/tools
if test -f "rankedlineage.dmp.gz"; then
	gunzip rankedlineage.dmp.gz
fi
if [[ "$threads" -le 4 ]]; then
	gthreads=threads
	Xmxg=$Xmx2
	gN=1
else
	gthreads=4
	gN=$(( threads / gthreads ))
	ramg=$(( ram2 / gN ))
	Xmxg=-Xmx${ramg}G
fi



##################################################################################################################
#Create all necessary subdirectories for downstream processing
cd $projdir
if [[ -d metagenome_ref_normalize ]]; then
	if [[ ! -d metagenome_no_normalize ]]; then
		mv ${projdir}/metagenome_ref_normalize ${projdir}/metagenome
	fi
fi
if [[ -d metagenome_no_normalize ]]; then
	if [[ ! -d metagenome_ref_normalize ]]; then
		mv ${projdir}/metagenome_no_normalize ${projdir}/metagenome
	fi
fi
if [[ -d metagenome_ref_normalize ]]; then
	if [[ -d metagenome_no_normalize ]]; then
		echo -e "${magenta}- Reference-normalized and non-normalized directories both exist ${white}\n"
		echo -e "${magenta}- Manually rename 1 the 2 existing metagenome-directories to metagenome ${white}\n"
		sleep 10; exit 0
	fi
fi
if [[ ! -d metagenome ]]; then
	mkdir metagenome
fi
cd metagenome
mkdir haplotig
mkdir alignment
mkdir sighits
mkdir results
mkdir ./results/ref_aligned_summaries

#################################################################################################################
#Organize fastq files for metagenomic processing
#Combine paired-end read data for downstream analysis
echo -e "\e[97m########################################################\n \e[38;5;210mOrganizing sample fastq files \n\e[97m########################################################\n"
organize_fq_files () {
cd $projdir
cd samples
if [[ -d "se" ]]; then
	:
else
	mkdir se
fi
if [[ -d "pe" ]]; then
	:
else
	mkdir pe
fi

cd se
if [[ -z "$(ls -A ../pe)" ]]; then
	if [[ -z "$(ls -A ../se)" ]]; then
		cd ../
		for i in $(ls *.f* | grep -v R2.f); do
			if [[ "$i" == *.R1* ]]; then
				mv $i ${i/.R1/}
			elif [[ "$i" == *_R1* ]]; then
				mv $i ${i/_R1/}
			fi
		done
	fi
fi

cd se
if [[ -z "$(ls -A ../pe)" ]]; then
	if [[ "$(ls -A ../se)" ]]; then
		echo -e "${magenta}- only single-end reads available in se-folder ${white}\n"
		for i in $(ls *.f* ); do
			if [[ "$i" == *.R1* ]]; then
				mv $i ../${i/.R1/}
			elif [[ "$i" == *_R1* ]]; then
				mv $i ../${i/_R1/}
			else
				mv $i ../$i
			fi
		done
	fi
fi

cd ../pe
if [[ -z "$(ls -A ../se)" ]]; then
	if [[ "$(ls -A ../pe)" ]]; then
		echo -e "${magenta}- only paired-end reads available in pe-folder ${white}\n"
		for i in $(ls *R1.f*); do
			mv ${i%R1.f*}R2.f* ../
			if [[ "$i" == *.R1* ]]; then
				mv $i ../${i/.R1/}
			elif [[ "$i" == *_R1* ]]; then
				mv $i ../${i/_R1/}
			else
				echo -e "${magenta}- check paired-end filenames for proper filename format (.R1 or _R1 and .R2 or _R2)  ${white}\n"
				echo -e "${magenta}- Do you wish to continue running Qmatey? ${white}\n"
				read -p "- y(YES) or n(NO) " -n 1 -r
				if [[ ! $REPLY =~ ^[Yy]$ ]]; then
					printf '\n'
					exit 0
				fi
			fi
		done
	fi
fi

cd ../pe
if [[ "$(ls -A ../se)" ]]; then
	if [[ "$(ls -A ../pe)" ]]; then
		for i in $(ls *R1.f*); do
			mv ${i%R1.f*}R2.f* ../
			if [[ "$i" == *.R1* ]]; then
				cat $i ../se/${i%.R1.f*}* > ../${i}
				mv ../${i} ../${i/.R1/}
				rm ../pe/$i ../se/${i%.R1.f*}*
			elif [[ "$i" == *_R1* ]]; then
				cat $i ../se/${i%_R1.f*}* > ../${i}
				mv ../${i} ../${i/_R1/}
				rm ../pe/$i ../se/${i%_R1.f*}*
			else
				echo -e "${magenta}- check paired-end filenames for proper filename format (.R1 or _R1 and .R2 or _R2) ${white}\n"
				echo -e "${magenta}- Do you wish to continue running Qmatey? ${white}\n"
				read -p "- y(YES) or n(NO) " -n 1 -r
				if [[ ! $REPLY =~ ^[Yy]$ ]]; then
					printf '\n'
					exit 0
				fi
			fi
		done
	fi
fi
cd ../
mkdir hold
for i in $(ls *.f* | grep -v R2.f); do
	checkfiles=$( ls ${i%.f*}.* | wc -l )
	if [[ "$checkfiles" -gt 1 ]]; then
		cat ${i%.f*}.* > ./hold/$i
		rm -r ${i%.f*}.*
		mv ./hold/$i ./
	fi
done
sampno=$(ls -1 | wc -l)
if [[ "$sampno" == "0" ]]; then
	echo -e "${magenta}- \n- samples folder is empty, exiting pipeline ${white}\n"
	exit 0
fi
find . -type d -empty -delete
}
cd $projdir
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
		sleep 5 && exit 0
	fi
else
	cd $projdir
	echo -e "${YELLOW}- Qmatey is organizing sample fastq files  ${WHITE}"
	time organize_fq_files &>> log.out
fi
total_no_samples=$(ls "${projdir}"/samples/* | wc -l)


cd $projdir
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
		for i in *.fa*; do
			n=">${i%.fa*}_"
			awk '{ sub("\r$",""); print}' $i | awk -v n="$n" '{gsub(n,">"); print}' | awk -v n="$n" '{gsub(/>/,n); print}' >> master_ref.fasta
		done

		ncontigscaffold=$(grep '>' master_ref.fasta | wc -l)
		if [[ $ncontigscaffold -gt 1000 ]]; then
			nfakechr=$((threads/2))
			awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < master_ref.fasta > master_ref.txt
			awk 'BEGIN{srand() }
				{ lines[++d]=$0 }
					END{
					while (1){
						if (e==d) {break}
							RANDOM = int(1 + rand() * d)
							if ( RANDOM in lines  ){
							print lines[RANDOM]
							delete lines[RANDOM]
							++e
						}
					}
				}' master_ref.txt > master_ref0.fasta
			flength=$(wc -l master_ref0.fasta | awk '{print $1}'); nsplit=$(( flength / nfakechr ))
			split -a 2 -d -l $nsplit master_ref0.fasta Chr
			rm master_ref*
			for i in Chr*; do
				Nstitch=$(printf "A%.0s" $(seq 100))
				awk -v Nstitch=$Nstitch '{print $1"\t"$2}' $i | awk '{gsub(/\t/,"\n"); print $0}' > ${i}.fasta
				grep -v '>'  ${i}.fasta |  awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\t",$0);next;} {printf("%s",$0);}END{printf("\n");}' | fold -w 100 > ${i}.txt
			done
			for filename in Chr*.txt; do
				echo ">""${filename%.txt}" >> master_ref.fasta
				cat "$filename" >> master_ref.fasta
			done
			rm Chr*
		fi
			$bwa index -a bwtsw master_ref.fasta
			$samtools faidx master_ref.fasta
			$java -jar $picard CreateSequenceDictionary REFERENCE=master_ref.fasta    OUTPUT=master_ref.dict
	fi
fi

#################################################################################################################
echo -e "\e[97m########################################################\n \e[38;5;210m Read Compression, Normalization, and exclusion of reads from reference genomes \n\e[97m########################################################\n"
ref_norm () {
	cd $projdir
	if [[ -z "$(ls -A ./norm_ref/*.dict 2> /dev/null)" ]]; then
		echo -e "$1 \e[31m normalization reference folder is empty, Qmatey will not exclude any read"
		echo -e "$1 \e[31m Qmatey will use read coverage of samples for normalization"
		cd $projdir/samples
		#All duplicate reads are compressed into one representative read with duplication reflected as a numeric value
		#Increases the speed of reference genome alignment -- especially if read depth is high
		rm ${projdir}/metagenome/microbiome_coverage.txt 2> /dev/null

		for i in $(ls -S *.f* | grep -v R2.f | grep -v _compressed.f); do

			if [[ $(file $i | awk -F' ' '{print $2}') == gzip ]]; then
				fa_fq=$(zcat $projdir/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
			else
				fa_fq=$(cat $projdir/samples/$i | head -n1 | cut -c1-1)
			fi
			wait

			if test ! -f ${i%.f*}_compressed.fasta.gz && [[ "${fa_fq}" == "@" ]]; then
				if test ! -f ${i%.f*}_R2.f* && test ! -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						awk 'NR%2==0' $i | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}_R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}_R2.f* 2> /dev/null) | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}_R2.f* | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}.R2.f* 2> /dev/null) | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}.R2.f* | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
			fi
			if test ! -f ${i%.f*}_compressed.fasta.gz && [[ "${fa_fq}" == ">" ]]; then
				if test ! -f ${i%.f*}_R2.f* && test ! -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz
					else
						awk 'NR%2==0' $i | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz
					fi
					wait
				fi
				if test -f ${i%.f*}_R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}_R2.f* 2> /dev/null) | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}_R2.f* | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}.R2.f* 2> /dev/null) | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}.R2.f* | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
			fi
		done
		wait
		#sample read depth is used to normalize quantification data
		echo -e "${YELLOW}- calculating a normalization factor"
		for i in $(ls -S *_compressed.fasta.gz); do
			zcat $i 2> /dev/null | awk '{print $1}' | awk -v sample=${i%_compressed.fasta.gz} -F '-' '{s+=$2}END{print sample"\t"s}' > ${projdir}/metagenome/${i%_compressed.fasta.gz}_microbiome_coverage.txt && \
			cat ${projdir}/metagenome/${i%_compressed.fasta.gz}_microbiome_coverage.txt >> ${projdir}/metagenome/microbiome_coverage.txt
			rm ${projdir}/metagenome/${i%_compressed.fasta.gz}_microbiome_coverage.txt
		done
		wait
		maximum=$(sort -nr -k2,2 ${projdir}/metagenome/microbiome_coverage.txt | awk 'NF > 0' | awk 'NR==1{print $2; exit}')
		awk -v maximum=$maximum '{print $1,maximum/$2}' ${projdir}/metagenome/microbiome_coverage.txt | cat <(printf 'Sample_ID\tNormalization_factor\n') - > ${projdir}/metagenome/coverage_normalization_factor.txt
		rm ${projdir}/metagenome/microbiome_coverage.txt
	else
		cd $projdir/samples
		#All duplicate reads are compressed into one representative read with duplication reflected as a numeric value
		#Increases the speed of reference genome alignment -- especially if read depth is high
		rm ${projdir}/metagenome/microbiome_coverage.txt 2> /dev/null

		for i in $(ls -S *.f* | grep -v R2.f | grep -v _compressed.f); do

			if [[ $(file $i | awk -F' ' '{print $2}') == gzip ]]; then
				fa_fq=$(zcat $projdir/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
			else
				fa_fq=$(cat $projdir/samples/$i | head -n1 | cut -c1-1)
			fi
			wait

			if test ! -f ${i%.f*}_compressed.fasta.gz && [[ "${fa_fq}" == "@" ]]; then
				if test ! -f ${i%.f*}_R2.f* && test ! -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						awk 'NR%2==0' $i | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}_R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}_R2.f* 2> /dev/null) | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}_R2.f* | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						 wait
					fi
					wait
				fi
				if test -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}.R2.f* 2> /dev/null) | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}.R2.f* | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
			fi
			if test ! -f ${i%.f*}_compressed.fasta.gz && [[ "${fa_fq}" == ">" ]]; then
				if test ! -f ${i%.f*}_R2.f* && test ! -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						awk 'NR%2==0' $i | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}_R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}_R2.f* 2> /dev/null) | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}_R2.f* | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}.R2.f* 2> /dev/null) | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}.R2.f* | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
			fi
		done
		wait

		cd ${projdir}/samples
		#Aligning compressed sample files (read-depth accounted for) to the master reference genomes
		#Exludes all host/reference genome data from the samples -- leaving only metagenomic reads
		echo -e "${YELLOW}- aligning sample reads to normalization reference genome${WHITE}"
		for i in $(ls *_compressed.fasta.gz); do
			if test ! -f ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt; then
				cd $projdir/norm_ref
				$bwa mem -t "$threads" master_ref.fasta <(zcat ${projdir}/samples/$i 2> /dev/null) > ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam && \
				cd $projdir/samples
				$java -XX:ParallelGCThreads=$gthreads -jar $picard SortSam I= ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam O= ${projdir}/metagenome/${i%_compressed.fasta.gz}.bam SORT_ORDER=coordinate && \
				printf '\n###---'${i%.f*}'---###\n' > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt && \
				$samtools flagstat ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt && \
				rm ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam
			fi
		done
		wait

		cat ${projdir}/metagenome/results/ref_aligned_summaries/*_summ.txt > ${projdir}/metagenome/results/ref_aligned_summaries/ref_aligned_summaries_unique_reads.txt
		rm -r ${projdir}/metagenome/results/ref_aligned_summaries/*_summ.txt

		cd $projdir/metagenome
		rm host_coverage.txt microbiome_coverage.txt
		#Host-reference alignment coverage relative to other samples is used to normalize quantification data
		echo -e "${YELLOW}- calculating a normalization factor"
		for i in $(ls -S *.bam); do
			$samtools view -F 4 $i | grep -vwE "(@HD|@SQ|@PG)" | awk '{print $1}' | awk -v sample=${i%.bam} -F '-' '{s+=$2}END{print sample"\t"s}' > ${i%.bam}_host_coverage.txt && \
			cat ${i%.bam}_host_coverage.txt >> host_coverage.txt && \
			rm ${i%.bam}_host_coverage.txt  && \
			$samtools view -f 4 $i | grep -vwE "(@HD|@SQ|@PG)" | awk '{print $1}' | awk -v sample=${i%.bam} -F '-' '{s+=$2}END{print sample"\t"s}' > ${i%.bam}_microbiome_coverage.txt && \
			cat ${i%.bam}_microbiome_coverage.txt >> microbiome_coverage.txt && \
			rm ${i%.bam}_microbiome_coverage.txt
		done
		wait
		awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' host_coverage.txt microbiome_coverage.txt | awk '{print $1,"\t",$4,"\t",$2}' | cat <(printf 'Sample_ID\tmetagenome_reads\ttotal_reads\n') - > coverage_normalize.txt
		maximum=$(sort -nr -k2,2 coverage_normalize.txt | awk 'NF > 0' | awk 'NR==1{print $2; exit}')
		awk -v maximum=$maximum 'NR>1{print $1,maximum/$2}' coverage_normalize.txt | cat <(printf 'Sample_ID\tNormalization_factor\n') - > coverage_normalization_factor.txt
		awk 'NR>1{print $1,($3/($2+$3))*100}' coverage_normalize.txt | cat <(printf 'Sample_ID\tPercent_metagenome\n') - > ./results/metagenome_derived_perc.txt
		rm host_coverage.txt microbiome_coverage.txt
	fi

	if [[ -z "$(ls -A $projdir/norm_ref/*.dict 2> /dev/null)" ]]; then
		cd $projdir/samples
		echo -e "${YELLOW}- compile metagenome reads & compute relative read depth ${WHITE}"
		for i in $(ls -S *_compressed.fasta.gz); do
			normfactor=$( awk -v sample=${i%_compressed.fasta.gz} '$1 == sample' ${projdir}/metagenome/coverage_normalization_factor.txt | awk '{print $2}' ) && \
			awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t"); print}' <(zcat $i 2> /dev/null) | awk -v norm=$normfactor '{print $1"-"$2*norm"\n"$3}' | $gzip > ${projdir}/metagenome/haplotig/${i%_compressed.fasta.gz}_metagenome.fasta.gz
		done
		wait
		rm ${projdir}/samples/*_compressed.fasta.gz

	else
		cd $projdir/metagenome
		echo -e "${YELLOW}- compile metagenome reads into fasta format & compute relative read depth ${WHITE}"
		for i in $(ls -S *.bam); do
			normfactor=$( awk -v sample=${i%.bam} '$1 == sample' coverage_normalization_factor.txt | awk '{print $2}' ) && \
			$samtools view -f 4 $i | grep -vwE "(@HD|@SQ|@PG)" | awk '{print $1"\t"$10}' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t"); print}' | \
			awk -v norm=$normfactor '{print ">"$1"-"$2*norm"\n"$3}' | $gzip > ./haplotig/${i%.bam}_metagenome.fasta.gz
		done
		wait
		rm *.bam ${projdir}/samples/*_compressed.fasta.gz
	fi

}
cd $projdir
if [[ "$normalization" == true ]]; then
	if ls ${projdir}/metagenome/coverage_normalization_factor.txt 1> /dev/null 2>&1; then
		echo -e "${YELLOW}- coverage normalization factor already computed ${WHITE}"
		echo -e "${YELLOW}- Qmatey is skipping normalization ${WHITE}"
	else
		echo -e "${YELLOW}- Qmatey is performing normalization and file compression ${WHITE}"
		time ref_norm &>> $projdir/log.out
	fi
else
	echo -e "${YELLOW}- Qmatey is skipping normalization ${WHITE}"
fi

#################################################################################################################
echo -e "\e[97m########################################################\n \e[38;5;210m Read compression and exclusion of reads from reference genomes \n\e[97m########################################################\n"
no_norm () {
	cd $projdir
	if [[ -z "$(ls -A ./norm_ref/*.dict 2> /dev/null)" ]]; then
		echo -e "$1 \e[31m normalization reference folder is empty, Qmatey will not exclude any read"
		cd ${projdir}/samples

		for i in $(ls -S *.f* | grep -v R2.f | grep -v _compressed.f); do

			if [[ $(file $i | awk -F' ' '{print $2}') == gzip ]]; then
				fa_fq=$(zcat $projdir/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
			else
				fa_fq=$(cat $projdir/samples/$i | head -n1 | cut -c1-1)
			fi
			wait

			if test ! -f ${i%.f*}_compressed.fasta.gz && [[ "${fa_fq}" == "@" ]]; then
				if test ! -f ${i%.f*}_R2.f* && test ! -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						awk 'NR%2==0' $i | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}_R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}_R2.f*) | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}_R2.f* | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}.R2.f* 2> /dev/null) | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}.R2.f* | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
			fi
			if test ! -f ${i%.f*}_compressed.fasta.gz && [[ "${fa_fq}" == ">" ]]; then
				if test ! -f ${i%.f*}_R2.f* && test ! -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						awk 'NR%2==0' $i | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}_R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}_R2.f* 2> /dev/null) | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}_R2.f* | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}.R2.f* 2> /dev/null) | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}.R2.f* | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\t"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
			fi
		done
		wait
	else
		cd $projdir/samples
		#All duplicate reads are compressed into one representative read with duplication reflected as a numeric value
		#Increased the spead of reference genome alignment -- especially if read depth is high
		for i in $(ls -S *.f* | grep -v R2.f | grep -v _compressed.f); do

			if [[ $(file $i | awk -F' ' '{print $2}') == gzip ]]; then
				fa_fq=$(zcat $projdir/samples/$i 2> /dev/null | head -n1 | cut -c1-1)
			else
				fa_fq=$(cat $projdir/samples/$i | head -n1 | cut -c1-1)
			fi
			wait

			if test ! -f ${i%.f*}_compressed.fasta.gz && [[ "${fa_fq}" == "@" ]]; then
				if test ! -f ${i%.f*}_R2.f* && test ! -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						awk 'NR%2==0' $i | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}_R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}_R2.f* 2> /dev/null) | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}_R2.f* | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}.R2.f* 2> /dev/null) | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}.R2.f* | awk 'NR%2==0' | awk 'NR%2==1' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
			fi
			if test ! -f ${i%.f*}_compressed.fasta.gz && [[ "${fa_fq}" == ">" ]]; then
				if test ! -f ${i%.f*}_R2.f* && test ! -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						awk 'NR%2==0' $i | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}_R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}_R2.f* 2> /dev/null) | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}_R2.f* | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
				if test -f ${i%.f*}.R2.f*; then
					if gzip -t $i; then
						zcat $i 2> /dev/null | cat - <(zcat ${i%.f*}.R2.f* 2> /dev/null) | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					else
						cat $i ${i%.f*}.R2.f* | awk 'NR%2==0' | sort | uniq -c | awk -v sample=${i%.f*} '{print ">"sample"_seq"NR"-"$1"\n"$2}' | $gzip > ${i%.f*}_compressed.fasta.gz &&
						wait
					fi
					wait
				fi
			fi
		done
		wait

		cd ${projdir}/samples
		#Aligning compressed sample files (read-depth accounted for) to the master reference genomes
		#Exludes all host/reference genome data from the samples -- leaving only metagenomic reads
		echo -e "${YELLOW}- aligning sample reads to normalization reference genome${WHITE}"
		for i in $(ls *_compressed.fasta.gz); do
			if test ! -f ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt; then
				cd ${projdir}/norm_ref/
				$bwa mem -t "$threads" master_ref.fasta <(zcat ${projdir}/samples/$i 2> /dev/null) > ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam && \
				cd ${projdir}/samples
				$java -XX:ParallelGCThreads=$gthreads -jar $picard SortSam I= ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam O= ${projdir}/metagenome/${i%_compressed.fasta.gz}.bam SORT_ORDER=coordinate && \
				printf '\n###---'${i%.f*}'---###\n' > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt && \
				$samtools flagstat ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam > ${projdir}/metagenome/results/ref_aligned_summaries/${i%_compressed.fasta.gz}_summ.txt && \
				rm ${projdir}/metagenome/${i%_compressed.fasta.gz}.sam
				wait
			fi
		done
		wait

		cat ${projdir}/metagenome/results/ref_aligned_summaries/*_summ.txt > ${projdir}/metagenome/results/ref_aligned_summaries/ref_aligned_summaries_unique_reads.txt
		rm -r ${projdir}/metagenome/results/ref_aligned_summaries/*_summ.txt

		cd $projdir/metagenome
		echo -e "${YELLOW}- compile metagenome reads into fasta format ${WHITE}"
		for i in $(ls -S *.bam); do
			if test ! -f ./haplotig/${i%.bam}_metagenome.fasta.gz; then
				$samtools view -f 4 $i | grep -vwE "(@HD|@SQ|@PG)" | awk '{print $1"\t"$10}' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t"); print}' | \
				awk '{print ">"$1"-"$2"\n"$3}' | $gzip > ./haplotig/${i%.bam}_metagenome.fasta.gz
			fi
		done
		wait
		rm *.bam ${projdir}/samples/*_compressed.fasta.gz
	fi
}
cd $projdir
if [[ "$normalization" == false ]]; then
	metagout=$(ls $projdir/metagenome/haplotig/*metagenome.fasta.gz 2> /dev/null | wc -l)
	samplein=$(ls $projdir/samples/*.f 2> /dev/null | grep compressed | wc -l)
	if [[ "$samplein" -eq "$metagout" ]]; then
		echo -e "${YELLOW}- Qmatey is skipping normalization and only performing file compression ${WHITE}"
		time no_norm &>> $projdir/log.out
	else
		echo -e "${YELLOW}- Qmatey has already performed file compression ${WHITE}"
	fi
else
	echo -e "${YELLOW}- Qmatey has already performed file compression ${WHITE}"
fi

#################################################################################################################
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey MegaBLAST \n\e[97m########################################################\n"

cd ${Qmatey_dir}
local_db=$( echo $local_db | awk '{gsub(/,/," ")}1' )
percid=90
if (echo $local_db | grep -q 'nt'); then
        percid=90
fi
if (echo $local_db | grep -q '16S') || (echo $local_db | grep -q '18S') || (echo $local_db | grep -q '28S') || (echo $local_db | grep -q 'ITS'); then
        percid=97
fi
if (echo $local_db | grep -q '16s') || (echo $local_db | grep -q '18s') || (echo $local_db | grep -q '28s') || (echo $local_db | grep -q 'ITs'); then
        percid=97
fi


blast () {

cd $projdir/metagenome/haplotig
if test ! -f combined_compressed_metagenomes.fasta.gz; then
	for compfa in *.fasta.gz; do
		zcat $compfa 2> /dev/null |  awk NF | awk 'NR%2==0' | $gzip >> hold_metagenomes.fasta.gz
	done
	wait
	awk -F, '!seen[$0]++' <(zcat hold_metagenomes.fasta.gz 2> /dev/null) | awk '{print ">seq"NR"\n"$1}' | $gzip > combined_compressed_metagenomes.fasta.gz &&
	rm hold_metagenomes.fasta.gz
fi

if [[ "$taxids" == true ]]; then
	for i in ${projdir}/taxids/*.txids; do
		cat $i >> ${projdir}/metagenome/All.txids
	done
fi


if [[ "$blast_location" == "local" ]]; then
	echo -e "${YELLOW}- performing local BLAST"

	file=${projdir}/metagenome/alignment/combined_compressed.megablast.gz
	if test -f $file; then
		echo -e "${YELLOW}- Primary BLAST (nt database) ouput already exist"
		echo -e "${YELLOW}- Skipping BLAST and filtering hits based on defined parameters"
	else
		if [[ -d splitccf ]]; then rm -r splitccf; fi
		mkdir splitccf; cd splitccf
		awk -v RS=">" -v FS="\n" -v ORS="\n" -v OFS="" '$0 {$1=">"$1"\n"; print}' <(zcat ../combined_compressed_metagenomes.fasta.gz 2> /dev/null) | $gzip > combined_compressed_metagenomes.fasta.gz
		awk 'NR%2000000==1{close("F"i); i++}{print > "F"i}'  <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) & PIDsplit1=$!
		wait $PIDsplit1
		rm combined_compressed_metagenomes.fasta.gz
		for ccf in $(ls * | sort -V); do
			mv $ccf ../../alignment/$ccf
			cd ../../alignment
			awk 'NR%200==1{close("subfile"i); i++}{print > "subfile"i}' $ccf & PIDsplit2=$!
			wait $PIDsplit2
			for sub in $(ls subfile* | sort -V); do (
				if [[ "$taxids" == true ]]; then
					${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/blastn -task megablast -query $sub -db "${local_db}" -num_threads 1 -perc_identity $percid -max_target_seqs 1000000 \
					-taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qcovs pident qseq sseq staxids stitle" -out ${sub}_out.blast &&
					wait
				else
					${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/blastn -task megablast -query $sub -db "${local_db}" -num_threads 1 -perc_identity $percid -max_target_seqs 1000000 \
					-outfmt "6 qseqid sseqid length qstart qcovs pident qseq sseq staxids stitle" -out ${sub}_out.blast &&
					wait
				fi
				wait
				gzip ${sub}_out.blast )&
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
			zcat ${ccf}.blast.gz 2> /dev/null | awk '$4 <= 6{print}' | awk '$3 >= 32{print}' | $gzip >> combined_compressed.megablast.gz &&
			rm ${ccf}.blast.gz; rm $ccf &&
			cd ../haplotig/splitccf/
		done
		cd ../
		rmdir splitccf
	fi

	for i in $(ls -S *metagenome.fasta.gz); do
	  if test ! -f ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast; then
		  awk NF=NF RS= OFS=@ <(zcat $i 2> /dev/null) | awk '{gsub(/@/,"\t");gsub(/>/,"\n"); print $0}' | awk 'NR>1' | sort -k2 | $gzip > ../alignment/${i%_metagenome.fasta.gz}.txt.gz &&
			wait
		  awk -F'\t' 'ORS=NR%2?"\t":"\n"' <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) | awk -F'\t' 'BEGIN{OFS="\t"}{gsub(/>/,"")}1' | \
		  awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}.txt.gz 2> /dev/null) - | \
		  awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3}' | awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' - <( zcat ../alignment/combined_compressed.megablast.gz 2> /dev/null) | \
		  awk -F'\t' 'BEGIN{OFS="\t"}{print $14,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | \
		  awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | $gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&
			wait
		  # if [[ "$taxids" = true ]]; then
		  #   for taxid_files in $(ls ${projdir}/taxids/*.txids); do
		  #     taxid=${taxid_files%*.txids}
		  #     taxid=${taxid/*\/}
		  #     awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($8) in a{print $0, a[$1]}' $taxid_files <(zcat ../alignment/${i%_metagenome.fasta.gz}_haplotig_temp.megablast.gz 2> /dev/null) | $gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid${taxid}.megablast.gz
		  #   done
		  #   wait
		  #   rm ../alignment/${i%_metagenome.fasta.gz}_haplotig_temp.megablast.gz
		  #   find ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid*.megablast.gz | xargs cat > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz &&
		  #   rm ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid*.megablast.gz
		  # fi
			rm ../alignment/${i%_metagenome.fasta.gz}.txt.gz
		fi
	done
fi

if [[ "$blast_location" == "remote" ]]; then
	echo -e "${YELLOW}- performing a remote BLAST"
	file=${projdir}/metagenome/alignment/combined_compressed.megablast.gz
	if test -f $file; then
		echo -e "${YELLOW}- BLAST ouput already exist"
		echo -e "${YELLOW}- Skipping BLAST and filtering hits based on defined parameters"
	else
		${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/blastn -task megablast -query <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) -db "${remote_db}" -perc_identity $percid  -max_target_seqs 1000000 \
		–taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qcovs pident gapopen qseq sseq staxids stitle" \
		-out ../alignment/combined_compressed.megablast -remote &&
		wait
	fi
	wait
	cat ../alignment/combined_compressed.megablast > ../alignment/temp.megablast
	awk 'BEGIN{FS="\t";}{if(a[$1]<$3){a[$1]=$3;}}END{for(i in a){print i"\t"a[i];}}' temp.megablast | sort -V -k1,1n | \
	awk -F'\t' 'BEGIN{FS=OFS="\t"} NR==FNR{c[$1FS$2]++;next};c[$1FS$3] > 0' - ../alignment/temp.megablast > ../alignment/combined_compressed.megablast
	rm ../alignment/temp.megablast

	for i in $(ls -S *metagenome.fasta.gz); do
	  if test ! -f ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast; then
		  awk NF=NF RS= OFS=@ <(zcat $i 2> /dev/null) | awk '{gsub(/@/,"\t");gsub(/>/,"\n"); print $0}' | awk 'NR>1' | sort -k2 > ../alignment/${i%_metagenome.fasta.gz}.txt &&
			wait
		  awk -F'\t' 'ORS=NR%2?"\t":"\n"' <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) | awk -F'\t' 'BEGIN{OFS="\t"}{gsub(/>/,"")}1' | \
		  awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' ../alignment/${i%_metagenome.fasta.gz}.txt - | \
		  awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3}' | awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' - ../alignment/combined_compressed.megablast | \
		  awk -F'\t' 'BEGIN{OFS="\t"}{print $14,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | \
		  awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast  &&
			wait
		  # if [[ "$taxids" == true ]]; then
		  #   for taxid_files in $(ls ${projdir}/taxids/*.txids); do
		  #     taxid=${taxid_files%*.txids}
		  #     taxid=${taxid/*\/}
		  #     awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($8) in a{print $0, a[$1]}' $taxid_files ../alignment/${i%_metagenome.fasta.gz}_haplotig_temp.megablast > ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid${taxid}.megablast
		  #   done
		  #   wait
		  #   rm ../alignment/${i%_metagenome.fasta.gz}_haplotig_temp.megablast
		  #   find ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid*.megablast | xargs cat | $gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz &&
		  #   rm ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid*.megablast
		  # fi
			rm ../alignment/${i%_metagenome.fasta.gz}.txt
		fi
	done
fi

if [[ "$blast_location" == "custom" ]]; then
	echo -e "${YELLOW}- performing local BLAST"
	file=${projdir}/metagenome/alignment/combined_compressed.megablast.gz
	if test -f $file; then
		echo -e "${YELLOW}- Primary BLAST ouput already exist"
		echo -e "${YELLOW}- Skipping BLAST and filtering hits based on defined parameters"
	else
		if [[ -d splitccf ]]; then rm -r splitccf; fi
		mkdir splitccf; cd splitccf
		awk -v RS=">" -v FS="\n" -v ORS="\n" -v OFS="" '$0 {$1=">"$1"\n"; print}' <(zcat ../combined_compressed_metagenomes.fasta.gz 2> /dev/null) | $gzip > combined_compressed_metagenomes.fasta.gz
		awk 'NR%2000000==1{close("F"i); i++}{print > "F"i}'  <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) & PIDsplit1=$!
		wait $PIDsplit1
		rm combined_compressed_metagenomes.fasta.gz
		for ccf in $(ls * | sort -V); do
			mv $ccf ../../alignment/$ccf
			cd ../../alignment
			awk 'NR%200==1{close("subfile"i); i++}{print > "subfile"i}' $ccf & PIDsplit2=$!
			wait $PIDsplit2
			for sub in $(ls subfile* | sort -V); do (
				if [[ "$taxids" == true ]]; then
					${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/blastn -task megablast -query $sub -db "${custom_db}" -num_threads 1 -perc_identity $percid  -max_target_seqs 1000000 \
					-taxidlist ${projdir}/metagenome/All.txids -outfmt "6 qseqid sseqid length qstart qcovs pident gapopen qseq sseq staxids stitle" -out ${sub}_out.blast &&
					wait
				else
					${Qmatey_dir}/tools/ncbi-blast-2.13.0+/bin/blastn -task megablast -query $sub -db "${custom_db}" -num_threads 1 -perc_identity $percid  -max_target_seqs 1000000 \
					-outfmt "6 qseqid sseqid length qstart qcovs pident gapopen qseq sseq staxids stitle" -out ${sub}_out.blast &&
					wait
				fi
				wait
				gzip ${sub}_out.blast )&
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
			zcat ${ccf}.blast.gz 2> /dev/null | awk '$4 <= 6{print}' | awk '$3 >= 32{print}' | $gzip >> combined_compressed.megablast.gz &&
			rm ${ccf}.blast.gz; rm $ccf &&
			cd ../haplotig/splitccf/
		done
		cd ../
		rmdir splitccf
	fi
	wait

	for i in $(ls -S *metagenome.fasta.gz); do
	  if test ! -f ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast; then
		  awk NF=NF RS= OFS=@ <(zcat $i 2> /dev/null) | awk '{gsub(/@/,"\t");gsub(/>/,"\n"); print $0}' | awk 'NR>1' | sort -k2 | $gzip > ../alignment/${i%_metagenome.fasta.gz}.txt.gz &&
			wait
		  awk -F'\t' 'ORS=NR%2?"\t":"\n"' <(zcat combined_compressed_metagenomes.fasta.gz 2> /dev/null) | awk -F'\t' 'BEGIN{OFS="\t"}{gsub(/>/,"")}1' | \
		  awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$2]=$0;next} ($2) in a{print $0, a[$2]}' <(zcat ../alignment/${i%_metagenome.fasta.gz}.txt.gz 2> /dev/null) - | \
		  awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2,$3}' | awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}' - <( zcat ../alignment/combined_compressed.megablast.gz 2> /dev/null) | \
		  awk -F'\t' 'BEGIN{OFS="\t"}{print $14,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}' | \
		  awk '{gsub(/^[ \t]+|[ \t]+$/,""); print;}' | $gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz  &&
			wait
		  # if [[ "$taxids" == true ]]; then
		  #   for taxid_files in $(ls ${projdir}/taxids/*.txids); do
		  #     taxid=${taxid_files%*.txids}
		  #     taxid=${taxid/*\/}
		  #     awk 'BEGIN{OFS="\t"} NR==FNR{a[$1]=$0;next} ($8) in a{print $0, a[$1]}' $taxid_files <(zcat ../alignment/${i%_metagenome.fasta.gz}_haplotig_temp.megablast.gz 2> /dev/null) | $gzip > ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid${taxid}.megablast.gz
		  #   done
		  #   wait
		  #   rm ../alignment/${i%_metagenome.fasta.gz}_haplotig_temp.megablast.gz
		  #   find ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid*.megablast.gz | xargs cat > ../alignment/${i%_metagenome.fasta.gz}_haplotig.megablast.gz &&
		  #   rm ../alignment/${i%_metagenome.fasta.gz}_haplotig_taxid*.megablast.gz
		  # fi
			rm ../alignment/${i%_metagenome.fasta.gz}.txt.gz
		fi
	done
fi

wait
find ../alignment/ -size 0 -delete
}
cd $projdir
metagout=$(ls $projdir/metagenome/haplotig/*metagenome.fasta.gz 2> /dev/null | wc -l)
mblastout=$(ls $projdir/metagenome/alignment/*haplotig.megablast.gz 2> /dev/null | wc -l)
if [[ "$metagout" -eq "$mblastout" ]]; then
	echo -e "${YELLOW}- Qmatey has already performed ncbi megablast ${WHITE}"
else
	echo -e "${YELLOW}- Qmatey is performing sequence alignment using ncbi megablast ${WHITE}"
	time blast &>> $projdir/log.out
fi


#################################################################################################################
awk '{gsub(/\t\t/,"\tNA\t"); print}' ${Qmatey_dir}/tools/rankedlineage.dmp | awk '{gsub(/[|]/,""); print}' | awk '{gsub(/\t\t/,"\t"); print}' > ${Qmatey_dir}/tools/rankedlineage_tabdelimited.dmp
echo $'tax_id\ttaxname\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain\t' | \
cat - ${Qmatey_dir}/tools/rankedlineage_tabdelimited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp
rm ${Qmatey_dir}/tools/rankedlineage_tabdelimited.dmp

lineagedb=${projdir}/lineage_subset.txt
if test -f $lineagedb; then
	echo -e "${YELLOW}- Compiling subset lineage file"
	if awk 'NR==1{/tax_id/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="tax_id") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > tax_id.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' tax_id.txt ${Qmatey_dir}/tools/rankedlineage_edited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp && rm tax_id.txt
	elif awk 'NR==1{/taxname/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="taxname") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > taxname.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$2] > 0' taxname.txt ${Qmatey_dir}/tools/rankedlineage_edited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp && rm taxname.txt
	elif awk 'NR==1{/genus/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="genus") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > genus.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$4] > 0' genus.txt ${Qmatey_dir}/tools/rankedlineage_edited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp && rm genus.txt
	elif awk 'NR==1{/family}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="family") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > family.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$5] > 0' family.txt ${Qmatey_dir}/tools/rankedlineage_edited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp && rm family.txt
	elif awk 'NR==1{/order/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="order") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > order.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$6] > 0' order.txt ${Qmatey_dir}/tools/rankedlineage_edited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp && rm order.txt
	elif awk 'NR==1{/class/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="class") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > class.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$7] > 0' class.txt ${Qmatey_dir}/tools/rankedlineage_edited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp && rm class.txt
	elif awk 'NR==1{/phylum/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="phylum") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > phylum.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$8] > 0' phylum.txt ${Qmatey_dir}/tools/rankedlineage_edited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp && rm phylum.txt
	elif awk 'NR==1{/kingdom/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="kingdom") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > kingdom.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$9] > 0' kingdom.txt ${Qmatey_dir}/tools/rankedlineage_edited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp && rm kingdom.txt
	elif awk 'NR==1{/domain/}' lineage_subset.txt;then
		awk -F '\t' 'NR==1{for(i=1; i<=NF; i++) if ($i=="domain") {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"}' lineage_db.txt > domain.txt
		awk -F '\t' 'NR==FNR{c[$1]++;next};c[$10] > 0' domain.txt ${Qmatey_dir}/tools/rankedlineage_edited.dmp > ${Qmatey_dir}/tools/rankedlineage_edited.dmp && rm domain.txt
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


#################################################################################################################
strain_level() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Strain-Level Classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact-matching algorithm"
cd $projdir/metagenome/sighits
mkdir -p sighits_strain
cd $projdir/metagenome/results
mkdir -p strain_level
cd $projdir/metagenome/alignment
if find ../sighits/sighits_strain/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits of at strain-level already available for each sample"
else
	for i in $(ls -S *_haplotig.megablast.gz); do
		zcat $i 2> /dev/null | awk '$5==100' | awk '$3<=10' | awk -F'\t' '{print $1"___"$10}' | sort | uniq | awk 'BEGIN{OFS="\t"}{gsub(/___/,"\t"); print $1}' | \
		uniq -c | awk -F ' ' '{print $2"\t"$1}' | awk '$2 == 1' | awk '{print $1}' > ${i%_haplotig.megablast.gz}_exactmatch.txt
		awk '$5==100' $i | awk '$3==1' | awk -F'\t' 'NR==FNR {a[$1]; next} $1 in a {print; delete a[$1]}' ${i%_haplotig.megablast.gz}_exactmatch.txt - | awk 'gsub(" ","_",$0)' | \
		awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | awk 'BEGIN{print "sseqid\tabundance\tqstart\tqcovs\tpident\tqseq\tsseq\tstaxids\tstitle\tqseqid\tfqseq"}{print $0}' | \
		awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_strain/${i%.gz} &&
		rm ${i%_haplotig.megablast.gz}_exactmatch.txt
	done
	wait
	cd $projdir/metagenome/haplotig
	for i in $(ls -S *_metagenome.fasta.gz); do
		awk '{ print ; nextfile}' <(zcat ../sighits/sighits_strain/${i%_metagenome.fasta.gz}_haplotig.megablast 2> /dev/null) | head -n 1 > ../sighits/sighits_strain/${i%_metagenome.fasta.gz}_sighits.txt &&
		ls <(zcat ../sighits/sighits_strain/${i%_metagenome.fasta.gz}_haplotig.megablast 2> /dev/null) | xargs -n 1 tail -n +2 >> ../sighits/sighits_strain/${i%_metagenome.fasta.gz}_sighits.txt &&
		rm ../sighits/sighits_strain/${i%_metagenome.fasta.gz}_haplotig.megablast
	done
fi

echo -e "${YELLOW}- compiling taxonomic information"
cd $projdir/metagenome/sighits/sighits_strain
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk '{print $8}' sighits.txt | awk '{gsub(";","\n"); print}' | sort -u -n | sed -e '1s/staxids/tax_id/' > taxids_sighits.txt && rm sighits.txt
awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}'  ${Qmatey_dir}/tools/rankedlineage_edited.dmp taxids_sighits.txt | \
awk '{gsub(/ /,"_"); print }' > rankedlineage_subhits.txt
rm taxids_sighits.txt
cd $projdir/metagenome/sighits/sighits_strain/
awk '{print $1}' rankedlineage_subhits.txt > strain_taxa_mean_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > strain_taxa_unique_sequences_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > strain_taxa_quantification_accuracy_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > strain_taxa_rel_quantification_accuracy_temp.txt

strain_level=strain
for i in $(ls *_sighits.txt);do
	Rscript "${Qmatey_dir}/scripts/stats_summary.R" $i $min_unique_seqs $strain_level "${Qmatey_dir}/tools/R"
	wait
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
	awk '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$11,$12,$13,$14,$15 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt &&
	rm *stats1* *stats2* *stats3* *hold*
done

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' strain_taxa_mean_temp.txt > strain_taxa_mean_temp2.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_mean_temp2.txt > ../../results/strain_level/strain_taxa_mean.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' strain_taxa_unique_sequences_temp.txt > strain_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_unique_sequences_temp2.txt > ../../results/strain_level/strain_taxa_unique_sequences.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_quantification_accuracy_temp.txt > strain_taxa_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/strain_level/strain_taxa_mean.txt strain_taxa_quantification_accuracy_temp2.txt > ../../results/strain_level/strain_taxa_quantification_accuracy.txt
awk '{gsub(/ /,"\t"); print}' strain_taxa_rel_quantification_accuracy_temp.txt > strain_taxa_rel_quantification_accuracy_temp2.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/strain_level/strain_taxa_mean.txt strain_taxa_rel_quantification_accuracy_temp2.txt > ../../results/strain_level/strain_taxa_rel_quantification_accuracy.txt &&
rm *_temp*


cd $projdir/metagenome/results/strain_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_strain/rankedlineage_subhits.txt strain_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > strain_taxainfo_${i}.txt &&
	wait
done
rm *_taxa_*


awk '{print $1}' strain_taxainfo_mean.txt > strain_taxainfo_mean_holdingtaxid.txt
awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
strain_taxainfo_mean.txt | awk '!($1="")' | awk 'BEGIN{OFS="\t"} {$1=$1};1' > strain_taxainfo_mean_holdingtaxinfo.txt
touch strain_taxainfo_mean_buildnorm.txt
for i in $(ls -S ../../../metagenome/haplotig/*_metagenome.fasta.gz); do
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
paste strain_taxainfo_mean_holdingtaxid.txt strain_taxainfo_mean_buildnorm.txt > strain_taxainfo_mean_norm0.txt
paste strain_taxainfo_mean_norm0.txt strain_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > strain_taxainfo_mean_normalized.txt
rm strain_taxainfo_mean_holdingtaxid.txt strain_taxainfo_mean_buildnorm.txt strain_taxainfo_mean_holdingtaxinfo.txt strain_taxainfo_mean_norm0.txt

for i in *.txt; do
	awk '{gsub(/-/,"_"); print}' $i > ${i%.txt}.temp &&
	mv ${i%.txt}.temp $i
done


Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" strain $genome_scaling "${Qmatey_dir}/tools/R" &>/dev/null

file=${projdir}/exclude_taxa.txt
if test -f $file; then
	cat strain_taxainfo_mean.txt > strain_taxainfo_mean_filtered.txt &&
	cat strain_taxainfo_unique_sequences.txt > strain_taxainfo_unique_sequences_filtered.txt &&
	cat strain_taxainfo_quantification_accuracy.txt > strain_taxainfo_quantification_accuracy_filtered.txt &&
	cat strain_taxainfo_rel_quantification_accuracy.txt > strain_taxainfo_rel_quantification_accuracy_filtered.txt &&
	while read -r line; do
		for i in $( ls *filtered.txt ); do
			awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i
		done
	done < $file
fi

if test -f $file; then
	echo -e "${YELLOW}- creating strain-level visualizations"
	cd $projdir/metagenome/results/strain_level
	strain_level_mean=strain_taxainfo_mean_filtered.txt
	strain_level_uniq=strain_taxainfo_unique_sequences_filtered.txt
	strain_level_stderr=strain_taxainfo_quantification_accuracy_filtered.txt
	strain_level_rel_stderr=strain_taxainfo_rel_quantification_accuracy_filtered.txt


	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
			Rscript "${Qmatey_dir}/scripts/strain_level_corr.R" "$strain_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" 2>/dev/null
			)&
		 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			 wait
		 fi
		done
	fi
	for min_perc in $min_percent_sample; do (
		Rscript "${Qmatey_dir}/scripts/strain_level_boxplots.R" "$strain_level_mean" "$strain_level_uniq" "$strain_level_stderr" "$strain_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	mkdir boxplots
	mv *_files $projdir/metagenome/results/strain_level/boxplots
	mv *.html $projdir/metagenome/results/strain_level/boxplots
else
	echo -e "${YELLOW}- creating strain-level visualizations"
	cd $projdir/metagenome/results/strain_level
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

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
			Rscript "${Qmatey_dir}/scripts/strain_level_corr.R" "$strain_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" 2>/dev/null
			)&
		 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			 wait
		 fi
		done
	fi
	for min_perc in $min_percent_sample; do (
		Rscript "${Qmatey_dir}/scripts/strain_level_boxplots.R" "$strain_level_mean" "$strain_level_uniq" "$strain_level_stderr" "$strain_level_rel_stderr" "$min_perc" "${Qmatey_dir}/tools/R" 2>/dev/null
		)&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		 wait
		fi
	done
	mkdir boxplots
	mv *_files $projdir/metagenome/results/strain_level/boxplots
	mv *.html $projdir/metagenome/results/strain_level/boxplots
fi

}
if [[ "$strain_level" == "true" ]]; then
	time strain_level 2>> $projdir/log.out
fi


#################################################################################################################
#################################################################################################################
#################################################################################################################
species() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey Species-Level Classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact-matching algorithm"
cd $projdir/metagenome/sighits
mkdir sighits_species
cd $projdir/metagenome/results
mkdir species_level
cd $projdir/metagenome/alignment
if find ../sighits/sighits_species/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits at species-level already available for each sample"
else
	for i in $(ls -S *_haplotig.megablast);do (
		awk '$3<=10' $i | awk -F'\t' 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
		awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_species/${i} )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		  wait
		fi
	done
	wait
	if [[ "$taxids" == true ]]; then
		:
	else
		ls ../sighits/sighits_species/*_haplotig.megablast | grep -v taxid | xargs rm
	fi

	cd $projdir/metagenome/haplotig
	for i in $(ls -S *metagenome.fasta); do
		awk '{ print ; nextfile}' ../sighits/sighits_species/*${i%_metagenome.fasta}_haplotig.megablast | head -n 1 > ../sighits/sighits_species/${i%_metagenome.fasta}_sighits.txt
		ls ../sighits/sighits_species/*${i%_metagenome.fasta}_haplotig.megablast | xargs -n 1 tail -n +2 >> ../sighits/sighits_species/${i%_metagenome.fasta}_sighits.txt
		rm ../sighits/sighits_species/*${i%_metagenome.fasta}_haplotig.megablast
	done

	cd $projdir/metagenome/sighits/sighits_species
	for i in $(ls *_sighits.txt); do (
		awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_species_unique_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	echo -e "${YELLOW}- compiling species-level multi-alignment algorithm"
	cd $projdir/metagenome/sighits/sighits_species
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits*}_species_unique_reads.txt $i OFS='\t' > ${i%_sighits*}_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	cd $projdir/metagenome/sighits/sighits_species
	for i in $(ls *_dup.txt);do (
		awk -F '\t' '{print $8}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_dup_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_dup_inter.txt
	cd $projdir/metagenome/sighits/sighits_species
	for i in $(ls -S *_taxids_dup.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait

	for i in $(ls *_dup_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_taxids_dup.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_dup.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_species_duplicates_virome.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_duplicates_virome.txt);do (
		awk -F '\t' '{ if ($21!="Viruses") print $0}' $i | awk -F '\t' '!/Uncultured/' > ${i%*_species_duplicates_virome*}_species_duplicates.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt && rm *_species_duplicates_virome.txt
	for i in $(ls *_species_duplicates.txt);do (
		awk -F '\t' '{print $1, $2"~"$14, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_species_duplicates*}_species_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_species_inter.txt);do (
		awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_species_inter*}_species_inter2.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_species_inter2.txt);do (
		awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_species_inter2*}_duplicate_count.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_duplicate_count.txt);do (
		awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_species_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_multialign_species_reads.txt);do (
		awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_species_reads*}_species_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_species_reads*}_species_EMC.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_species_unique_reads.txt);do (
		awk -F '\t' '{print $11}' OFS=';' $i > ${i%_species_unique_reads*}_taxids_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_uniq_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_uniq_inter.txt
	cd $projdir/metagenome/sighits/sighits_species
	for i in $(ls *_taxids_uniq.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_uniq_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_taxids_uniq.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_unique_reads.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_species_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_species_unique_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_species_column.txt && rm *_species_taxa.txt
	for i in $(ls *_species_unique_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_species_unique_uncultured*}_unique_sequences.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_unique_uncultured.txt && rm *_species_inter.txt && rm *_species_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_species_reads.txt && rm *_species_duplicates.txt
	for i in $(ls *_species_EMC.txt);do (
		cat $i ${i%_species_EMC*}_species_unique_sequences.txt > ${i%_species_EMC*}_complete_species_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	echo -e "${YELLOW}- compiling taxonomic information"
	for i in $(ls *_complete_species_reads.txt);do (
		awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_species_reads*}_sighits_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tgapopen\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
		cat - ${i%_complete_species_reads*}_sighits_temp.txt > ${i%_complete_species_reads*}_sighits_temp2.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_sighits_temp2.txt);do (
		awk -F '\t' '{gsub(/ /,"_");print}' $i > ${i%_sighits_temp2*}_sighits.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_complete_species_reads.txt && rm *_sighits_temp.txt && rm *_species_EMC.txt && rm *_unique_reads.txt && rm *_unique_sequences.txt && rm *_sighits_temp2.txt
fi


cd $projdir/metagenome/sighits/sighits_species
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits.txt && rm sighits.txt
cd $projdir/metagenome/sighits/sighits_species/
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/species/' > species_taxa_mean_temp1.txt
echo -e 'species' | cat - species_taxa_mean_temp1.txt > species_taxa_mean_temp.txt && rm species_taxa_mean_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/species/' > species_taxa_unique_sequences_temp1.txt
echo -e 'species' | cat - species_taxa_unique_sequences_temp1.txt > species_taxa_unique_sequences_temp.txt && rm species_taxa_unique_sequences_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/species/' > species_taxa_quantification_accuracy_temp1.txt
echo -e 'species' | cat - species_taxa_quantification_accuracy_temp1.txt > species_taxa_quantification_accuracy_temp.txt && rm species_taxa_quantification_accuracy_temp1.txt
awk -F '\t' '{print $1}' rankedlineage_subhits.txt | sort -u | awk '!/species/' > species_taxa_rel_quantification_accuracy_temp1.txt
echo -e 'species' | cat - species_taxa_rel_quantification_accuracy_temp1.txt > species_taxa_rel_quantification_accuracy_temp.txt && rm species_taxa_rel_quantification_accuracy_temp1.txt


species_level=species
for i in $(ls *_sighits.txt); do
	Rscript "${Qmatey_dir}/scripts/stats_summary.R" $i $min_unique_seqs $species_level "${Qmatey_dir}/tools/R"
	echo $'species\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt species_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > species_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt species_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > species_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt species_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > species_taxa_quantification_accuracy_temp.txt
	id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdrelstderr.txt species_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt && cat holdrelstderr2.txt > species_taxa_rel_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$11,$12 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' species_taxa_mean_temp.txt > species_taxa_mean_temp2.txt && rm species_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_mean_temp2.txt > ../../results/species_level/species_taxa_mean.txt && rm species_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' species_taxa_unique_sequences_temp.txt > species_taxa_unique_sequences_temp2.txt && rm species_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_unique_sequences_temp2.txt > ../../results/species_level/species_taxa_unique_sequences.txt && rm species_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_quantification_accuracy_temp.txt > species_taxa_quantification_accuracy_temp2.txt && rm species_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/species_level/species_taxa_mean.txt species_taxa_quantification_accuracy_temp2.txt > ../../results/species_level/species_taxa_quantification_accuracy.txt && rm species_taxa_quantification_accuracy_temp2.txt
awk '{gsub(/ /,"\t"); print}' species_taxa_rel_quantification_accuracy_temp.txt > species_taxa_rel_quantification_accuracy_temp2.txt && rm species_taxa_rel_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/species_level/species_taxa_mean.txt species_taxa_rel_quantification_accuracy_temp2.txt > ../../results/species_level/species_taxa_rel_quantification_accuracy.txt && rm species_taxa_rel_quantification_accuracy_temp2.txt

cd $projdir/metagenome/results/species_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do (
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_species/rankedlineage_subhits.txt species_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > species_taxainfo_${i}.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done
rm *_taxa_*


awk '{print $1}' species_taxainfo_mean.txt > species_taxainfo_mean_holdingtaxid.txt
awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
species_taxainfo_mean.txt | awk '!($1="")' | awk 'BEGIN{OFS="\t"} {$1=$1};1' > species_taxainfo_mean_holdingtaxinfo.txt
touch species_taxainfo_mean_buildnorm.txt
for i in $(ls -S ../../../metagenome/haplotig/*_metagenome.fasta); do (
	sample=${i%_metagenome.fasta}; sample=${sample##*/}
	if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
		normfactor=1
	else
		normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
	fi
	awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" species_taxainfo_mean.txt | \
	paste - species_taxainfo_mean_buildnorm.txt > species_taxainfo_mean_buildnorm0.txt
	cat species_taxainfo_mean_buildnorm0.txt > species_taxainfo_mean_buildnorm.txt
	rm species_taxainfo_mean_buildnorm0.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done
paste species_taxainfo_mean_holdingtaxid.txt species_taxainfo_mean_buildnorm.txt > species_taxainfo_mean_norm0.txt
paste species_taxainfo_mean_norm0.txt species_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > species_taxainfo_mean_normalized.txt
rm species_taxainfo_mean_holdingtaxid.txt species_taxainfo_mean_buildnorm.txt species_taxainfo_mean_holdingtaxinfo.txt species_taxainfo_mean_norm0.txt

for i in *.txt; do (
	awk '{gsub(/-/,"_"); print}' $i > temp.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
	mv temp.txt $i
done


Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" species $genome_scaling "${Qmatey_dir}/tools/R" &>/dev/null

file=${projdir}/exclude_taxa.txt
if test -f $file; then
	cat species_taxainfo_mean.txt > species_taxainfo_mean_filtered.txt
	cat species_taxainfo_unique_sequences.txt > species_taxainfo_unique_sequences_filtered.txt
	cat species_taxainfo_quantification_accuracy.txt > species_taxainfo_quantification_accuracy_filtered.txt
	cat species_taxainfo_rel_quantification_accuracy.txt > species_taxainfo_rel_quantification_accuracy_filtered.txt
	while read -r line; do
		for i in $( ls *filtered.txt ); do (
			awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i )&
			if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			  wait
			fi
		done
	done < $file
fi


if test -f $file; then
	echo -e "${YELLOW}- creating species-level visualizations"
	cd $projdir/metagenome/results/species_level
	species_level_mean=species_taxainfo_mean_filtered.txt
	species_level_uniq=species_taxainfo_unique_sequences_filtered.txt
	species_level_stderr=species_taxainfo_quantification_accuracy_filtered.txt
	species_level_rel_stderr=species_taxainfo_rel_quantification_accuracy_filtered.txt
	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/species_level_corr.R" "$species_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/species_level_boxplots.R" "$species_level_mean" "$species_level_uniq" "$species_level_stderr" "$species_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/species_level/boxplots
	mv *.html $projdir/metagenome/results/species_level/boxplots
else
	echo -e "${YELLOW}- creating species-level visualizations"
	cd $projdir/metagenome/results/species_level
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

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/species_level_corr.R" "$species_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/species_level_boxplots.R" "$species_level_mean" "$species_level_uniq" "$species_level_stderr" "$species_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/species_level/boxplots
	mv *.html $projdir/metagenome/results/species_level/boxplots
fi
}
if [[ "$species_level" == true ]]; then
	time species 2>> $projdir/log.out
fi

#################################################################################################################
#################################################################################################################
#################################################################################################################
genus() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey Genus-Level Classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact-matching algorithm"
cd $projdir/metagenome/sighits
mkdir sighits_genus
cd $projdir/metagenome/results
mkdir genus_level
cd $projdir/metagenome/alignment
if find ../sighits/sighits_genus/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits at genus-level already available for each sample"
else
	for i in $(ls -S *_haplotig.megablast);do (
		awk '$3<=10' $i | awk -F'\t' 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
		awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_genus/${i} )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		  wait
		fi
	done
	wait

	if [[ "$taxids" == true ]]; then
		:
	else
		ls ../sighits/sighits_genus/*_haplotig.megablast | grep -v taxid | xargs rm
	fi

	cd $projdir/metagenome/haplotig
	for i in $(ls -S *metagenome.fasta); do
		awk '{ print ; nextfile}' ../sighits/sighits_genus/*${i%_metagenome.fasta}_haplotig.megablast | head -n 1 > ../sighits/sighits_genus/${i%_metagenome.fasta}_sighits.txt
		ls ../sighits/sighits_genus/*${i%_metagenome.fasta}_haplotig.megablast | xargs -n 1 tail -n +2 >> ../sighits/sighits_genus/${i%_metagenome.fasta}_sighits.txt
		rm ../sighits/sighits_genus/*${i%_metagenome.fasta}_haplotig.megablast
	done


	echo -e "${YELLOW}- compiling genus-level multi-alignment algorithm"
	cd $projdir/metagenome/sighits/sighits_genus/
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_genus_unique_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	cd $projdir/metagenome/sighits/sighits_genus
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits*}_genus_unique_reads.txt $i OFS='\t' > ${i%_sighits*}_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	cd $projdir/metagenome/sighits/sighits_genus
	for i in $(ls *_dup.txt);do (
		awk -F '\t' '{print $8}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_dup_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_taxids_dup_inter.txt

	cd $projdir/metagenome/sighits/sighits_genus
	for i in $(ls *_taxids_dup.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait

	for i in $(ls *_dup_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_taxids_dup.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_dup.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_genus_duplicates_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
	for i in $(ls *_genus_duplicates_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_genus_duplicates_uncultured*}_genus_duplicates.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_genus_duplicates_uncultured.txt
	for i in $(ls *_genus_duplicates.txt);do (
		awk -F '\t' '{print $1, $2"~"$15, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_genus_duplicates*}_genus_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_genus_inter.txt);do (
		awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_genus_inter*}_genus_inter2.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_genus_inter2.txt);do (
		awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_genus_inter2*}_duplicate_count.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_duplicate_count.txt);do (
		awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_genus_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_multialign_genus_reads.txt);do (
		awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_genus_reads*}_genus_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $2, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_genus_reads*}_genus_EMC.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		  wait
		fi
	done
	wait

	for i in $(ls *_genus_unique_reads.txt);do (
		awk -F '\t' '{print $11}' OFS=';' $i > ${i%_genus_unique_reads*}_taxids_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_uniq_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_taxids_uniq_inter.txt
	cd $projdir/metagenome/sighits/sighits_genus
	for i in $(ls *_taxids_uniq.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_uniq_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_taxids_uniq.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_unique_reads.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_genus_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_genus_unique_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_species_column.txt && rm *_species_taxa.txt

	for i in $(ls *_genus_unique_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_genus_unique_uncultured*}_unique_sequences.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_genus_unique_uncultured.txt && rm *_genus_inter.txt && rm *_genus_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_genus_reads.txt && rm *_genus_duplicates.txt
	for i in $(ls *_genus_EMC.txt);do (
		cat $i ${i%_genus_EMC*}_genus_unique_sequences.txt > ${i%_genus_EMC*}_complete_genus_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	echo -e "${YELLOW}- compiling taxonomic information"
	for i in $(ls *_complete_genus_reads.txt);do (
		awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_genus_reads*}_sighits_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tgapopen\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
		cat - ${i%_complete_genus_reads*}_sighits_temp.txt > ${i%_complete_genus_reads*}_sighits.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_complete_genus_reads.txt && rm *_sighits_temp.txt && rm *_genus_EMC.txt && rm *_unique_reads.txt && rm *_genus_unique_sequences.txt
fi

cd $projdir/metagenome/sighits/sighits_genus
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt && rm sighits.txt
awk -F '\t' 'NR==1; NR >1 {print $0 | "sort -u"}' rankedlineage_subhits_temp.txt> rankedlineage_subhits.txt && rm rankedlineage_subhits_temp.txt


cd $projdir/metagenome/sighits/sighits_genus
echo -e "${YELLOW}- quantifying the genus-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt > genus_taxa_mean_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > genus_taxa_unique_sequences_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > genus_taxa_quantification_accuracy_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > genus_taxa_rel_quantification_accuracy_temp.txt



genus_level=genus
for i in $(ls *_sighits.txt);do
	Rscript "${Qmatey_dir}/scripts/stats_summary.R" $i $min_unique_seqs $genus_level "${Qmatey_dir}/tools/R"
	echo $'genus\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt genus_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > genus_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt genus_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > genus_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt genus_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > genus_taxa_quantification_accuracy_temp.txt
	id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdrelstderr.txt genus_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt && cat holdrelstderr2.txt > genus_taxa_rel_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$5,$7,$8,$9,$10,$11 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' genus_taxa_mean_temp.txt > genus_taxa_mean_temp2.txt && rm genus_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_mean_temp2.txt > ../../results/genus_level/genus_taxa_mean.txt && rm genus_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' genus_taxa_unique_sequences_temp.txt > genus_taxa_unique_sequences_temp2.txt && rm genus_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_unique_sequences_temp2.txt > ../../results/genus_level/genus_taxa_unique_sequences.txt && rm genus_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_quantification_accuracy_temp.txt > genus_taxa_quantification_accuracy_temp2.txt && rm genus_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/genus_level/genus_taxa_mean.txt genus_taxa_quantification_accuracy_temp2.txt > genus_taxa_quantification_accuracy_temp3.txt && rm genus_taxa_quantification_accuracy_temp2.txt
sed '2,${/genus/d;}' genus_taxa_quantification_accuracy_temp3.txt > ../../results/genus_level/genus_taxa_quantification_accuracy.txt && rm genus_taxa_quantification_accuracy_temp3.txt
awk '{gsub(/ /,"\t"); print}' genus_taxa_rel_quantification_accuracy_temp.txt > genus_taxa_rel_quantification_accuracy_temp2.txt && rm genus_taxa_rel_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/genus_level/genus_taxa_mean.txt genus_taxa_rel_quantification_accuracy_temp2.txt > genus_taxa_rel_quantification_accuracy_temp3.txt && rm genus_taxa_rel_quantification_accuracy_temp2.txt
sed '2,${/genus/d;}' genus_taxa_rel_quantification_accuracy_temp3.txt > ../../results/genus_level/genus_taxa_rel_quantification_accuracy.txt && rm genus_taxa_rel_quantification_accuracy_temp3.txt

cd $projdir/metagenome/results/genus_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do (
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_genus/rankedlineage_subhits.txt genus_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > genus_taxainfo_${i}.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
	  wait
	fi
done
rm *_taxa_*


cd $projdir/metagenome/results/genus_level
awk '{print $1}' genus_taxainfo_mean.txt > genus_taxainfo_mean_holdingtaxid.txt
awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
genus_taxainfo_mean.txt | awk '!($1="")' | awk 'BEGIN{OFS="\t"} {$1=$1};1' > genus_taxainfo_mean_holdingtaxinfo.txt
touch genus_taxainfo_mean_buildnorm.txt
for i in $(ls -S ../../../metagenome/haplotig/*_metagenome.fasta); do (
	sample=${i%_metagenome.fasta}; sample=${sample##*/}
	if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
		normfactor=1
	else
		normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
	fi
	awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" genus_taxainfo_mean.txt | \
	paste - genus_taxainfo_mean_buildnorm.txt > genus_taxainfo_mean_buildnorm0.txt
	cat genus_taxainfo_mean_buildnorm0.txt > genus_taxainfo_mean_buildnorm.txt
	rm genus_taxainfo_mean_buildnorm0.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done
paste genus_taxainfo_mean_holdingtaxid.txt genus_taxainfo_mean_buildnorm.txt > genus_taxainfo_mean_norm0.txt
paste genus_taxainfo_mean_norm0.txt genus_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > genus_taxainfo_mean_normalized.txt
rm genus_taxainfo_mean_holdingtaxid.txt genus_taxainfo_mean_buildnorm.txt genus_taxainfo_mean_holdingtaxinfo.txt genus_taxainfo_mean_norm0.txt

for i in *.txt; do
	awk '{gsub(/-/,"_"); print}' $i > temp.txt
	mv temp.txt $i
done


Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" genus $genome_scaling "${Qmatey_dir}/tools/R" &>/dev/null

file=${projdir}/exclude_taxa.txt
if test -f $file; then
	cat genus_taxainfo_mean.txt > genus_taxainfo_mean_filtered.txt
	cat genus_taxainfo_unique_sequences.txt > genus_taxainfo_unique_sequences_filtered.txt
	cat genus_taxainfo_quantification_accuracy.txt > genus_taxainfo_quantification_accuracy_filtered.txt
	cat genus_taxainfo_rel_quantification_accuracy.txt > genus_taxainfo_rel_quantification_accuracy_filtered.txt
	while read -r line; do
		for i in $( ls *filtered.txt ); do (
			awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i )&
			if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			  wait
			fi
		done
	done < $file
fi

if test -f $file; then
	echo -e "${YELLOW}- creating genus-level visualizations"
	cd $projdir/metagenome/results/genus_level
	genus_level_mean=genus_taxainfo_mean_filtered.txt
	genus_level_uniq=genus_taxainfo_unique_sequences_filtered.txt
	genus_level_stderr=genus_taxainfo_quantification_accuracy_filtered.txt
	genus_level_rel_stderr=genus_taxainfo_rel_quantification_accuracy_filtered.txt
	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/genus_level_corr.R" "$genus_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/genus_level_boxplots.R" "$genus_level_mean" "$genus_level_uniq" "$genus_level_stderr" "$genus_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/genus_level/boxplots
	mv *.html $projdir/metagenome/results/genus_level/boxplots
else
	echo -e "${YELLOW}- creating genus-level visualizations"
	cd $projdir/metagenome/results/genus_level
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

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/genus_level_corr.R" "$genus_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/genus_level_boxplots.R" "$genus_level_mean" "$genus_level_uniq" "$genus_level_stderr" "$genus_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/genus_level/boxplots
	mv *.html $projdir/metagenome/results/genus_level/boxplots
fi

}
if [[ "$genus_level" == true ]]; then
	time genus 2>> $projdir/log.out
fi

#################################################################################################################
#################################################################################################################
#################################################################################################################
family() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Family-Level Classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact-matching algorithm"
cd $projdir/metagenome/sighits
mkdir sighits_family
cd $projdir/metagenome/results
mkdir family_level
cp seqcov.txt ./family_level
cd $projdir/metagenome/alignment
if find ../sighits/sighits_family/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits at family-level already available for each sample"
else
	for i in $(ls -S *_haplotig.megablast);do (
		awk '$3<=10' $i | awk -F'\t' 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
		awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_family/${i} )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		  wait
		fi
	done
	wait

	if [[ "$taxids" == true ]]; then
		:
	else
		ls ../sighits/sighits_family/*_haplotig.megablast | grep -v taxid | xargs rm
	fi

	cd $projdir/metagenome/haplotig
	for i in $(ls -S *metagenome.fasta); do
		awk '{ print ; nextfile}' ../sighits/sighits_family/*${i%_metagenome.fasta}_haplotig.megablast | head -n 1 > ../sighits/sighits_family/${i%_metagenome.fasta}_sighits.txt
		ls ../sighits/sighits_family/*${i%_metagenome.fasta}_haplotig.megablast | xargs -n 1 tail -n +2 >> ../sighits/sighits_family/${i%_metagenome.fasta}_sighits.txt
		rm ../sighits/sighits_family/*${i%_metagenome.fasta}_haplotig.megablast
	done

	echo -e "${YELLOW}- compiling family-level multi-alignment algorithm"
	cd $projdir/metagenome/sighits/sighits_family/
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_family_unique_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	cd $projdir/metagenome/sighits/sighits_family
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits*}_family_unique_reads.txt $i OFS='\t' > ${i%_sighits*}_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	cd $projdir/metagenome/sighits/sighits_family
	for i in $(ls *_dup.txt);do (
		awk -F '\t' '{print $8}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_dup_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_taxids_dup_inter.txt
	cd $projdir/metagenome/sighits/sighits_family
	for i in $(ls *_taxids_dup.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_dup_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_taxids_dup.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_dup.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_family_duplicates_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
	for i in $(ls *_family_duplicates_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_family_duplicates_uncultured*}_family_duplicates.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_family_duplicates_uncultured.txt
	for i in $(ls *_family_duplicates.txt);do (
		awk -F '\t' '{print $1, $2"~"$16, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_family_duplicates*}_family_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_family_inter.txt);do (
		awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_family_inter*}_family_inter2.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_family_inter2.txt);do (
		awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_family_inter2*}_duplicate_count.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_duplicate_count.txt);do (
		awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_family_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_multialign_family_reads.txt);do (
		awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_family_reads*}_family_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_family_reads*}_family_EMC.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		  wait
		fi
	done
	wait

	for i in $(ls *_family_unique_reads.txt);do (
		awk -F '\t' '{print $11}' OFS=';' $i > ${i%_family_unique_reads*}_taxids_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_uniq_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_taxids_uniq_inter.txt
	cd $projdir/metagenome/sighits/sighits_family
	for i in $(ls *_taxids_uniq.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_uniq_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_taxids_uniq.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_unique_reads.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_family_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_unique_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_uniq.txt && rm *_species_column.txt && rm *_species_taxa.txt
	for i in $(ls *_family_unique_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_family_unique_uncultured*}_unique_sequences.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_family_unique_uncultured.txt && rm *_family_inter.txt && rm *_family_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_family_reads.txt && rm *_family_duplicates.txt
	for i in $(ls *_family_EMC.txt);do (
		cat $i ${i%_family_EMC*}_unique_sequences.txt > ${i%_family_EMC*}_complete_family_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait


	for i in $(ls *_complete_family_reads.txt);do (
		awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_family_reads*}_sighits_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tgapopen\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
		cat - ${i%_complete_family_reads*}_sighits_temp.txt > ${i%_complete_family_reads*}_sighits.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_complete_family_reads.txt && rm *_sighits_temp.txt && rm *_family_EMC.txt && rm *_unique_reads.txt && rm *_unique_sequences.txt && rm *_species_column.txt
	rm *_species_taxa.txt
fi

cd $projdir/metagenome/sighits/sighits_family
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt
awk -F '\t' '!/family/' rankedlineage_subhits_temp.txt > rankedlineage_subhits_temp2.txt
echo $'family\torder\tclass\tphylum\tkingdom\tdomain' | cat - rankedlineage_subhits_temp2.txt > rankedlineage_subhits_temp3.txt
awk -F '\t' 'NR==1; NR >1 {print $0 | "sort -u"}' rankedlineage_subhits_temp3.txt> rankedlineage_subhits.txt
rm *temp* && rm sighits.txt


cd $projdir/metagenome/sighits/sighits_family
echo -e "${YELLOW}- quantifying the family-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt > family_taxa_mean_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > family_taxa_unique_sequences_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > family_taxa_quantification_accuracy_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > family_taxa_rel_quantification_accuracy_temp.txt


family_level=family
for i in $(ls *_sighits.txt);do
	Rscript "${Qmatey_dir}/scripts/stats_summary.R" $i $min_unique_seqs $family_level "${Qmatey_dir}/tools/R"
	echo $'family\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt family_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > family_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt family_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > family_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt family_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > family_taxa_quantification_accuracy_temp.txt
	id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdrelstderr.txt family_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt && cat holdrelstderr2.txt > family_taxa_rel_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$5,$7,$8,$9,$10 }' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' family_taxa_mean_temp.txt > family_taxa_mean_temp2.txt && rm family_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_mean_temp2.txt > ../../results/family_level/family_taxa_mean.txt && rm family_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' family_taxa_unique_sequences_temp.txt > family_taxa_unique_sequences_temp2.txt && rm family_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_unique_sequences_temp2.txt > ../../results/family_level/family_taxa_unique_sequences.txt && rm family_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_quantification_accuracy_temp.txt > family_taxa_quantification_accuracy_temp2.txt && rm family_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/family_level/family_taxa_mean.txt family_taxa_quantification_accuracy_temp2.txt > family_taxa_quantification_accuracy_temp3.txt && rm family_taxa_quantification_accuracy_temp2.txt
sed '2,${/family/d;}' family_taxa_quantification_accuracy_temp3.txt > ../../results/family_level/family_taxa_quantification_accuracy.txt && rm family_taxa_quantification_accuracy_temp3.txt
awk '{gsub(/ /,"\t"); print}' family_taxa_rel_quantification_accuracy_temp.txt > family_taxa_rel_quantification_accuracy_temp2.txt && rm family_taxa_rel_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/family_level/family_taxa_mean.txt family_taxa_rel_quantification_accuracy_temp2.txt > family_taxa_rel_quantification_accuracy_temp3.txt && rm family_taxa_rel_quantification_accuracy_temp2.txt
sed '2,${/family/d;}' family_taxa_rel_quantification_accuracy_temp3.txt > ../../results/family_level/family_taxa_rel_quantification_accuracy.txt && rm family_taxa_rel_quantification_accuracy_temp3.txt

cd $projdir/metagenome/results/family_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do (
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_family/rankedlineage_subhits.txt family_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' | awk '$0!=x;(NR==1){x=$0}' > family_taxainfo_${i}.txt
	)&
 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	 wait
 fi
done
rm *_taxa_*


cd $projdir/metagenome/results/family_level
awk '{print $1}' family_taxainfo_mean.txt > family_taxainfo_mean_holdingtaxid.txt
awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
family_taxainfo_mean.txt | awk '!($1="")' | awk 'BEGIN{OFS="\t"} {$1=$1};1' > family_taxainfo_mean_holdingtaxinfo.txt
touch family_taxainfo_mean_buildnorm.txt
for i in $(ls -S ../../../metagenome/haplotig/*_metagenome.fasta); do (
	sample=${i%_metagenome.fasta}; sample=${sample##*/}
	if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
		normfactor=1
	else
		normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
	fi
	awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" family_taxainfo_mean.txt | \
	paste - family_taxainfo_mean_buildnorm.txt > family_taxainfo_mean_buildnorm0.txt
	cat family_taxainfo_mean_buildnorm0.txt > family_taxainfo_mean_buildnorm.txt
	rm family_taxainfo_mean_buildnorm0.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done
paste family_taxainfo_mean_holdingtaxid.txt family_taxainfo_mean_buildnorm.txt > family_taxainfo_mean_norm0.txt
paste family_taxainfo_mean_norm0.txt family_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > family_taxainfo_mean_normalized.txt
rm family_taxainfo_mean_holdingtaxid.txt family_taxainfo_mean_buildnorm.txt family_taxainfo_mean_holdingtaxinfo.txt family_taxainfo_mean_norm0.txt

for i in *.txt; do (
	awk '{gsub(/-/,"_"); print}' $i > temp.txt
	mv temp.txt $i )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done


Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" family $genome_scaling "${Qmatey_dir}/tools/R" &>/dev/null

file=${projdir}/exclude_taxa.txt
if test -f $file; then
	cat family_taxainfo_mean.txt > family_taxainfo_mean_filtered.txt
	cat family_taxainfo_unique_sequences.txt > family_taxainfo_unique_sequences_filtered.txt
	cat family_taxainfo_quantification_accuracy.txt > family_taxainfo_quantification_accuracy_filtered.txt
	cat family_taxainfo_rel_quantification_accuracy.txt > family_taxainfo_rel_quantification_accuracy_filtered.txt
	while read -r line; do
		for i in $( ls *filtered.txt ); do (
			awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i )&
			if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
				wait
			fi
		done
	done > $file
fi

if test -f $file; then
	echo -e "${YELLOW}- creating family-level visualizations"
	cd $projdir/metagenome/results/family_level
	family_level_mean=family_taxainfo_mean_filtered.txt
	family_level_uniq=family_taxainfo_unique_sequences_filtered.txt
	family_level_stderr=family_taxainfo_quantification_accuracy_filtered.txt
	family_level_rel_stderr=family_taxainfo_rel_quantification_accuracy_filtered.txt
	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/family_level_corr.R" "$family_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/family_level_boxplots.R" "$family_level_mean" "$family_level_uniq" "$family_level_stderr" "$family_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/family_level/boxplots
	mv *.html $projdir/metagenome/results/family_level/boxplots
else
	echo -e "${YELLOW}- creating family-level visualizations"
	cd $projdir/metagenome/results/family_level
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

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/family_level_corr.R" "$family_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/family_level_boxplots.R" "$family_level_mean" "$family_level_uniq" "$family_level_stderr" "$family_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/family_level/boxplots
	mv *.html $projdir/metagenome/results/family_level/boxplots
fi

}
if [[ "$family_level" == true ]]; then
	time family 2>> $projdir/log.out
fi

#################################################################################################################
#################################################################################################################
#################################################################################################################
order() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Order-Level Classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact-matching algorithm"
cd $projdir/metagenome/sighits
mkdir sighits_order
cd $projdir/metagenome/results
mkdir order_level
cd $projdir/metagenome/alignment
if find ../sighits/sighits_order/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits at order-level already available for each sample"
else
	for i in $(ls -S *_haplotig.megablast);do (
		awk '$3<=10' $i | awk -F'\t' 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
		awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_order/${i} )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		  wait
		fi
	done
	wait

	if [[ "$taxids" == true ]]; then
		:
	else
		ls ../sighits/sighits_order/*_haplotig.megablast | grep -v taxid | xargs rm
	fi

	cd $projdir/metagenome/haplotig
	awk '{ print ; nextfile}' ../sighits/sighits_order/*_haplotig.megablast | head -n 1 > ../sighits/sighits_order/${i%_metagenome.fasta}_sighits.txt
	for i in $(ls -S *metagenome.fasta); do
		ls ../sighits/sighits_order/*${i%_metagenome.fasta}_haplotig.megablast | xargs -n 1 tail -n +2 >> ../sighits/sighits_order/${i%_metagenome.fasta}_sighits.txt
		rm ../sighits/sighits_order/*${i%_metagenome.fasta}_haplotig.megablast
	done

	echo -e "${YELLOW}- compiling order-level multi-alignment algorithm"
	cd $projdir/metagenome/sighits/sighits_order/
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_order_unique_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	cd $projdir/metagenome/sighits/sighits_order
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits*}_order_unique_reads.txt $i OFS='\t' > ${i%_sighits*}_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	cd $projdir/metagenome/sighits/sighits_order
	for i in $(ls *_dup.txt);do (
		awk -F '\t' '{print $8}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_dup_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_taxids_dup_inter.txt
	cd $projdir/metagenome/sighits/sighits_order
	for i in $(ls *_taxids_dup.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_dup_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_taxids_dup.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_dup.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_order_duplicates_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
	for i in $(ls *_order_duplicates_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_order_duplicates_uncultured*}_order_duplicates.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_order_duplicates_uncultured.txt
	for i in $(ls *_order_duplicates.txt);do (
		awk -F '\t' '{print $1, $2"~"$17, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_order_duplicates*}_order_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_order_inter.txt);do (
		awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_order_inter*}_order_inter2.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_order_inter2.txt);do (
		awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_order_inter2*}_duplicate_count.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_duplicate_count.txt);do (
		awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_order_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_multialign_order_reads.txt);do (
		awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_order_reads*}_order_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_order_reads*}_order_EMC.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		  wait
		fi
	done
	wait
	rm *_order_inter.txt && rm *_order_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_order_reads.txt && rm *_order_duplicates.txt

	for i in $(ls *_order_unique_reads.txt);do (
		awk -F '\t' '{print $11}' OFS=';' $i > ${i%_order_unique_reads*}_taxids_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_uniq_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_taxids_uniq_inter.txt
	cd $projdir/metagenome/sighits/sighits_order
	for i in $(ls *_taxids_uniq.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_uniq_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_taxids_uniq.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_unique_reads.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_order_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_order_unique_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_uniq.txt && rm *_species_column.txt && rm *_species_taxa.txt
	for i in $(ls *_order_unique_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_order_unique_uncultured*}_order_unique_sequences.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_order_unique_uncultured.txt &&
	for i in $(ls *_order_EMC.txt);do (
		cat $i ${i%_order_EMC*}_order_order_unique_sequences.txt > ${i%_order_EMC*}_complete_order_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_complete_order_reads.txt);do (
		awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_order_reads*}_sighits_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tgapopen\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
		cat - ${i%_complete_order_reads*}_sighits_temp.txt > ${i%_complete_order_reads*}_sighits.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_complete_order_reads.txt && rm *_sighits_temp.txt && rm *_order_EMC.txt && rm *_unique_reads.txt && rm *_unique_sequences.txt && rm *_species_column.txt && rm *_species_taxa.txt
fi

cd $projdir/metagenome/sighits/sighits_order
echo -e "${YELLOW}- compiling taxonomic information"
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $16"\t"$17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt && rm sighits.txt
awk -F '\t' '!/order/' rankedlineage_subhits_temp.txt > rankedlineage_subhits_temp2.txt && rm rankedlineage_subhits_temp.txt
echo $'order\tclass\tphylum\tkingdom\tdomain' | cat - rankedlineage_subhits_temp2.txt > rankedlineage_subhits_temp3.txt && rm rankedlineage_subhits_temp2.txt
awk -F '\t' 'NR==1; NR >1 {print $0 | "sort -u"}' rankedlineage_subhits_temp3.txt> rankedlineage_subhits.txt && rm rankedlineage_subhits_temp3.txt


cd $projdir/metagenome/sighits/sighits_order
echo -e "${YELLOW}- quantifying the order-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt > order_taxa_mean_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > order_taxa_unique_sequences_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > order_taxa_quantification_accuracy_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > order_taxa_rel_quantification_accuracy_temp.txt



order_level=order
for i in $(ls *_sighits.txt);do
	Rscript "${Qmatey_dir}/scripts/stats_summary.R" $i $min_unique_seqs $order_level "${Qmatey_dir}/tools/R"
	echo $'order\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt order_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > order_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt order_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > order_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt order_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > order_taxa_quantification_accuracy_temp.txt
	id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdrelstderr.txt order_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt && cat holdrelstderr2.txt > order_taxa_rel_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$5,$7,$8,$9}' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' order_taxa_mean_temp.txt > order_taxa_mean_temp2.txt && rm order_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_mean_temp2.txt > ../../results/order_level/order_taxa_mean.txt && rm order_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' order_taxa_unique_sequences_temp.txt > order_taxa_unique_sequences_temp2.txt && rm order_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_unique_sequences_temp2.txt > ../../results/order_level/order_taxa_unique_sequences.txt && rm order_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_quantification_accuracy_temp.txt > order_taxa_quantification_accuracy_temp2.txt && rm order_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/order_level/order_taxa_mean.txt order_taxa_quantification_accuracy_temp2.txt > order_taxa_quantification_accuracy_temp3.txt && rm order_taxa_quantification_accuracy_temp2.txt
sed '2,${/order/d;}' order_taxa_quantification_accuracy_temp3.txt > ../../results/order_level/order_taxa_quantification_accuracy.txt && rm *order_taxa_quantification_accuracy_temp3.txt
awk '{gsub(/ /,"\t"); print}' order_taxa_rel_quantification_accuracy_temp.txt > order_taxa_rel_quantification_accuracy_temp2.txt && rm order_taxa_rel_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/order_level/order_taxa_mean.txt order_taxa_rel_quantification_accuracy_temp2.txt > order_taxa_rel_quantification_accuracy_temp3.txt && rm order_taxa_rel_quantification_accuracy_temp2.txt
sed '2,${/order/d;}' order_taxa_rel_quantification_accuracy_temp3.txt > ../../results/order_level/order_taxa_rel_quantification_accuracy.txt && rm *order_taxa_rel_quantification_accuracy_temp3.txt


cd $projdir/metagenome/results/order_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do (
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_order/rankedlineage_subhits.txt order_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > order_taxainfo_${i}.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
	  wait
	fi
done
rm *_taxa_*


cd $projdir/metagenome/results/order_level
awk '{print $1}' order_taxainfo_mean.txt > order_taxainfo_mean_holdingtaxid.txt
awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
order_taxainfo_mean.txt | awk '!($1="")' | awk 'BEGIN{OFS="\t"} {$1=$1};1' > order_taxainfo_mean_holdingtaxinfo.txt
touch order_taxainfo_mean_buildnorm.txt
for i in $(ls -S ../../../metagenome/haplotig/*_metagenome.fasta); do (
	sample=${i%_metagenome.fasta}; sample=${sample##*/}
	if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
		normfactor=1
	else
		normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
	fi
	awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" order_taxainfo_mean.txt | \
	paste - order_taxainfo_mean_buildnorm.txt > order_taxainfo_mean_buildnorm0.txt
	cat order_taxainfo_mean_buildnorm0.txt > order_taxainfo_mean_buildnorm.txt
	rm order_taxainfo_mean_buildnorm0.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done
paste order_taxainfo_mean_holdingtaxid.txt order_taxainfo_mean_buildnorm.txt > order_taxainfo_mean_norm0.txt
paste order_taxainfo_mean_norm0.txt order_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > order_taxainfo_mean_normalized.txt
rm order_taxainfo_mean_holdingtaxid.txt order_taxainfo_mean_buildnorm.txt order_taxainfo_mean_holdingtaxinfo.txt order_taxainfo_mean_norm0.txt

for i in *.txt; do (
	awk '{gsub(/-/,"_"); print}' $i > temp.txt
	mv temp.txt $i )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done


Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" order $genome_scaling "${Qmatey_dir}/tools/R" &>/dev/null

file=${projdir}/exclude_taxa.txt
if test -f $file; then
	cat order_taxainfo_mean.txt > order_taxainfo_mean_filtered.txt
	cat order_taxainfo_unique_sequences.txt > order_taxainfo_unique_sequences_filtered.txt
	cat order_taxainfo_quantification_accuracy.txt > order_taxainfo_quantification_accuracy_filtered.txt
	cat order_taxainfo_rel_quantification_accuracy.txt > order_taxainfo_rel_quantification_accuracy_filtered.txt
	while read -r line; do
		for i in $( ls *filtered.txt ); do (
			awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i )&
			if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			  wait
			fi
		done
	done > $file
fi

if test -f $file; then
	echo -e "${YELLOW}- creating order-level visualizations"
	cd $projdir/metagenome/results/order_level
	order_level_mean=order_taxainfo_mean_filtered.txt
	order_level_uniq=order_taxainfo_unique_sequences_filtered.txt
	order_level_stderr=order_taxainfo_quantification_accuracy_filtered.txt
	order_level_rel_stderr=order_taxainfo_rel_quantification_accuracy_filtered.txt
	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/order_level_corr.R" "$order_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/order_level_boxplots.R" "$order_level_mean" "$order_level_uniq" "$order_level_stderr" "$order_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/order_level/boxplots
	mv *.html $projdir/metagenome/results/order_level/boxplots
else
	echo -e "${YELLOW}- creating order-level visualizations"
	cd $projdir/metagenome/results/order_level
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

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/order_level_corr.R" "$order_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/order_level_boxplots.R" "$order_level_mean" "$order_level_uniq" "$order_level_stderr" "$order_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/order_level/boxplots
	mv *.html $projdir/metagenome/results/order_level/boxplots
fi
}
if [[ "$order_level" == true ]]; then
	time order 2>> $projdir/log.out
fi

#################################################################################################################
#################################################################################################################
#################################################################################################################
class() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Class-Level Classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact-matching algorithm"
cd $projdir/metagenome/sighits
mkdir sighits_class
cd $projdir/metagenome/results
mkdir class_level
cd $projdir/metagenome/alignment
if find ../sighits/sighits_class/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits at class-level already available for each sample"
else
	for i in $(ls -S *_haplotig.megablast); do (
		awk '$3<=10' $i | awk -F'\t' 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
		awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_class/${i} )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		  wait
		fi
	done
	wait

	if [[ "$taxids" == true ]]; then
		:
	else
		ls ../sighits/sighits_class/*_haplotig.megablast | grep -v taxid | xargs rm
	fi

	cd $projdir/metagenome/haplotig
	for i in $(ls -S *metagenome.fasta); do
		awk '{ print ; nextfile}' ../sighits/sighits_class/*${i%_metagenome.fasta}_haplotig.megablast | head -n 1 > ../sighits/sighits_class/${i%_metagenome.fasta}_sighits.txt
		ls ../sighits/sighits_class/*${i%_metagenome.fasta}_haplotig.megablast | xargs -n 1 tail -n +2 >> ../sighits/sighits_class/${i%_metagenome.fasta}_sighits.txt
		rm ../sighits/sighits_class/*${i%_metagenome.fasta}_haplotig.megablast
	done

	echo -e "${YELLOW}- compiling class-level multi-alignment algorithm"
	cd $projdir/metagenome/sighits/sighits_class/
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_class_unique_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	cd $projdir/metagenome/sighits/sighits_class
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits*}_class_unique_reads.txt $i OFS='\t' > ${i%_sighits*}_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	cd $projdir/metagenome/sighits/sighits_class
	for i in $(ls *_dup.txt);do (
		awk -F '\t' '{print $8}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_dup_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_taxids_dup_inter.txt
	cd $projdir/metagenome/sighits/sighits_class
	for i in $(ls *_taxids_dup.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_dup_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_taxids_dup.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_dup.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_class_duplicates_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
	for i in $(ls *_class_duplicates_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_class_duplicates_uncultured*}_class_duplicates.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_class_duplicates_uncultured.txt
	for i in $(ls *_class_duplicates.txt);do (
		awk -F '\t' '{print $1, $2"~"$18, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_class_duplicates*}_class_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_class_inter.txt);do (
		awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_class_inter*}_class_inter2.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_class_inter2.txt);do (
		awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_class_inter2*}_duplicate_count.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_duplicate_count.txt);do (
		awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_class_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_multialign_class_reads.txt);do (
		awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_class_reads*}_class_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_class_reads*}_class_EMC.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_class_inter.txt && rm *_class_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_class_reads.txt && rm *_class_duplicates.txt
	for i in $(ls *_class_unique_reads.txt);do (
		awk -F '\t' '{print $11}' OFS=';' $i > ${i%_class_unique_reads*}_taxids_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	for i in $(ls *_uniq_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done

	rm *_taxids_uniq_inter.txt
	cd $projdir/metagenome/sighits/sighits_class
	for i in $(ls *_taxids_uniq.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_uniq_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	rm *_taxids_uniq.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_unique_reads.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_class_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_class_unique_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_uniq.txt && rm *_species_column.txt && rm *_species_taxa.txt
	for i in $(ls *_class_unique_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_class_unique_uncultured*}_class_unique_sequences.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	for i in $(ls *_class_EMC.txt);do (
		cat $i ${i%_class_EMC*}_class_class_unique_sequences.txt > ${i%_class_EMC*}_complete_class_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	wait
	for i in $(ls *_complete_class_reads.txt);do (
		awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_class_reads*}_sighits_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tgapopen\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
		cat - ${i%_complete_class_reads*}_sighits_temp.txt > ${i%_complete_class_reads*}_sighits.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		  wait
		fi
	done
	rm *_complete_class_reads.txt && rm *_sighits_temp.txt && rm *_class_EMC.txt && rm *_unique_reads.txt && rm *_unique_sequences.txt && rm *_uncultured.txt
	rm *_species_column.txt && rm *_species_taxa.txt
fi

cd $projdir/metagenome/sighits/sighits_class
echo -e "${YELLOW}- compiling taxonomic information"
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $17"\t"$18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt && rm sighits.txt
awk -F '\t' '!/class/' rankedlineage_subhits_temp.txt > rankedlineage_subhits_temp2.txt && rm rankedlineage_subhits_temp.txt
echo $'class\tphylum\tkingdom\tdomain' | cat - rankedlineage_subhits_temp2.txt > rankedlineage_subhits_temp3.txt && rm rankedlineage_subhits_temp2.txt
awk -F '\t' 'NR==1; NR >1 {print $0 | "sort -u"}' rankedlineage_subhits_temp3.txt > rankedlineage_subhits.txt && rm rankedlineage_subhits_temp3.txt

cd $projdir/metagenome/sighits/sighits_class
echo -e "${YELLOW}- quantifying the class-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt > class_taxa_mean_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > class_taxa_unique_sequences_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > class_taxa_quantification_accuracy_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > class_taxa_rel_quantification_accuracy_temp.txt



class_level=class
for i in $(ls *_sighits.txt);do
	Rscript "${Qmatey_dir}/scripts/stats_summary.R" $i $min_unique_seqs $class_level "${Qmatey_dir}/tools/R"
	echo $'class\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt class_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > class_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt class_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > class_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt class_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > class_taxa_quantification_accuracy_temp.txt
	id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdrelstderr.txt class_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt && cat holdrelstderr2.txt > class_taxa_rel_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$5,$7,$8}' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' class_taxa_mean_temp.txt > class_taxa_mean_temp2.txt && rm class_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_mean_temp2.txt > ../../results/class_level/class_taxa_mean.txt && rm class_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' class_taxa_unique_sequences_temp.txt > class_taxa_unique_sequences_temp2.txt && rm class_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_unique_sequences_temp2.txt > ../../results/class_level/class_taxa_unique_sequences.txt && rm class_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_quantification_accuracy_temp.txt > class_taxa_quantification_accuracy_temp2.txt && rm class_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/class_level/class_taxa_mean.txt class_taxa_quantification_accuracy_temp2.txt > class_taxa_quantification_accuracy_temp3.txt && rm class_taxa_quantification_accuracy_temp2.txt
sed '2,${/class/d;}' class_taxa_quantification_accuracy_temp3.txt > ../../results/class_level/class_taxa_quantification_accuracy.txt && rm class_taxa_quantification_accuracy_temp3.txt
awk '{gsub(/ /,"\t"); print}' class_taxa_rel_quantification_accuracy_temp.txt > class_taxa_rel_quantification_accuracy_temp2.txt && rm class_taxa_rel_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/class_level/class_taxa_mean.txt class_taxa_rel_quantification_accuracy_temp2.txt > class_taxa_rel_quantification_accuracy_temp3.txt && rm class_taxa_rel_quantification_accuracy_temp2.txt
sed '2,${/class/d;}' class_taxa_rel_quantification_accuracy_temp3.txt > ../../results/class_level/class_taxa_rel_quantification_accuracy.txt && rm class_taxa_rel_quantification_accuracy_temp3.txt

cd $projdir/metagenome/results/class_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do (
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_class/rankedlineage_subhits.txt class_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > class_taxainfo_${i}.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done
rm *_taxa_*


cd $projdir/metagenome/results/class_level
awk '{print $1}' class_taxainfo_mean.txt > class_taxainfo_mean_holdingtaxid.txt
awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
class_taxainfo_mean.txt | awk '!($1="")' | awk 'BEGIN{OFS="\t"} {$1=$1};1' > class_taxainfo_mean_holdingtaxinfo.txt
touch class_taxainfo_mean_buildnorm.txt
for i in $(ls -S ../../../metagenome/haplotig/*_metagenome.fasta); do (
	sample=${i%_metagenome.fasta}; sample=${sample##*/}
	if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
		normfactor=1
	else
		normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
	fi
	awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" class_taxainfo_mean.txt | \
	paste - class_taxainfo_mean_buildnorm.txt > class_taxainfo_mean_buildnorm0.txt
	cat class_taxainfo_mean_buildnorm0.txt > class_taxainfo_mean_buildnorm.txt
	rm class_taxainfo_mean_buildnorm0.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done
paste class_taxainfo_mean_holdingtaxid.txt class_taxainfo_mean_buildnorm.txt > class_taxainfo_mean_norm0.txt
paste class_taxainfo_mean_norm0.txt class_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > class_taxainfo_mean_normalized.txt
rm class_taxainfo_mean_holdingtaxid.txt class_taxainfo_mean_buildnorm.txt class_taxainfo_mean_holdingtaxinfo.txt class_taxainfo_mean_norm0.txt

for i in *.txt; do (
	awk '{gsub(/-/,"_"); print}' $i > temp.txt
	mv temp.txt $i )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
	  wait
	fi
done


Rscript "${Qmatey_dir}/scripts/phylum_level_genome_scaling.R" class $genome_scaling "${Qmatey_dir}/tools/R" &>/dev/null

file=${projdir}/exclude_taxa.txt
if test -f $file; then
	cat class_taxainfo_mean.txt > class_taxainfo_mean_filtered.txt
	cat class_taxainfo_unique_sequences.txt > class_taxainfo_unique_sequences_filtered.txt
	cat class_taxainfo_quantification_accuracy.txt > class_taxainfo_quantification_accuracy_filtered.txt
	cat class_taxainfo_rel_quantification_accuracy.txt > class_taxainfo_rel_quantification_accuracy_filtered.txt
	while read -r line; do
		for i in $( ls *filtered.txt ); do (
			awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i )&
			if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			  wait
			fi
		done
	done > $file
fi


if test -f $file; then
	echo -e "${YELLOW}- creating class-level visualizations"
	cd $projdir/metagenome/results/class_level
	class_level_mean=class_taxainfo_mean_filtered.txt
	class_level_uniq=class_taxainfo_unique_sequences_filtered.txt
	class_level_stderr=class_taxainfo_quantification_accuracy_filtered.txt
	class_level_rel_stderr=class_taxainfo_rel_quantification_accuracy_filtered.txt
	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/class_level_corr.R" "$class_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/class_level_boxplots.R" "$class_level_mean" "$class_level_uniq" "$class_level_stderr" "$class_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/class_level/boxplots
	mv *.html $projdir/metagenome/results/class_level/boxplots
else
	echo -e "${YELLOW}- creating class-level visualizations"
	cd $projdir/metagenome/results/class_level
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

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/class_level_corr.R" "$class_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/class_level_boxplots.R" "$class_level_mean" "$class_level_uniq" "$class_level_stderr" "$class_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/class_level/boxplots
	mv *.html $projdir/metagenome/results/class_level/boxplots
fi
}
if [[ "$class_level" == true ]]; then
	time class 2>> $projdir/log.out
fi

#################################################################################################################
#################################################################################################################
#################################################################################################################
phylum() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is Performing Phylum-Level Classification \n\e[97m########################################################\n"
echo -e "${YELLOW}- performing exact-matching algorithm"
cd $projdir/metagenome/sighits
mkdir sighits_phylum
cd $projdir/metagenome/results
mkdir phylum_level
cd $projdir/metagenome/alignment
if find ../sighits/sighits_phylum/ -mindepth 1 | read; then
	echo -e "${YELLOW}- significant hits at phylum-level already available for each sample"
else
	for i in $(ls -S *_haplotig.megablast);do (
		awk '$3<=10' $i | awk -F'\t' 'gsub(" ","_",$0)' | awk 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$1); print}' | \
		awk '{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' | awk 'gsub(" ","\t",$0)' > ../sighits/sighits_phylum/${i} )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			wait
		fi
	done
	wait

	if [[ "$taxids" == true ]]; then
		:
	else
		ls ../sighits/sighits_phylum/*_haplotig.megablast | grep -v taxid | xargs rm
	fi

	cd $projdir/metagenome/haplotig
	for i in $(ls -S *metagenome.fasta); do
		awk '{ print ; nextfile}' ../sighits/sighits_phylum/*${i%_metagenome.fasta}_haplotig.megablast | head -n 1 > ../sighits/sighits_phylum/${i%_metagenome.fasta}_sighits.txt
		ls ../sighits/sighits_phylum/*${i%_metagenome.fasta}_haplotig.megablast | xargs -n 1 tail -n +2 >> ../sighits/sighits_phylum/${i%_metagenome.fasta}_sighits.txt
		rm ../sighits/sighits_phylum/*${i%_metagenome.fasta}_haplotig.megablast
	done

	echo -e "${YELLOW}- compiling phylum-level multi-alignment algorithm"
	cd $projdir/metagenome/sighits/sighits_phylum/
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' '{a[$2]++;b[$2]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i > ${i%_sighits*}_phylum_unique_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	cd $projdir/metagenome/sighits/sighits_phylum
	for i in $(ls *_sighits.txt);do (
		awk -F '\t' 'FNR==NR{a[$1,$2]=1; next}  !a[$1,$2]' ${i%_sighits*}_phylum_unique_reads.txt $i OFS='\t' > ${i%_sighits*}_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
	cd $projdir/metagenome/sighits/sighits_phylum
	for i in $(ls *_dup.txt);do (
		awk -F '\t' '{print $8}' OFS=';' $i > ${i%_dup*}_taxids_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done

	for i in $(ls *_dup_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_dup_inter*}_taxids_dup.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done

	rm *_taxids_dup_inter.txt
	cd $projdir/metagenome/sighits/sighits_phylum
	for i in $(ls *_taxids_dup.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_dup*}_dup_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
	for i in $(ls *_dup_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_dup_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
	rm *_taxids_dup.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	for i in $(ls *_dup.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9}' OFS='\t' ${i%*_dup.txt}_species_taxa.txt) > ${i%_dup*}_phylum_duplicates_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	rm *_species_taxid.txt && rm *_dup_inter.txt && rm *_dup.txt && rm *_species_column.txt && rm *_species_taxa.txt
	for i in $(ls *_phylum_duplicates_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_phylum_duplicates_uncultured*}_phylum_duplicates.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	rm *_phylum_duplicates_uncultured.txt
	for i in $(ls *_phylum_duplicates.txt);do (
		awk -F '\t' '{print $1, $2"~"$19, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_phylum_duplicates*}_phylum_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
	for i in $(ls *_phylum_inter.txt);do (
		awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' $i > ${i%_phylum_inter*}_phylum_inter2.txt
		)&
	 if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		 wait
	 fi
	done
	wait
	for i in $(ls *_phylum_inter2.txt);do (
		awk -F '\t' '{dups[$1]++} END {for (num in dups) {print num}}' $i | sort -k1,1  > ${i%_phylum_inter2*}_duplicate_count.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
	for i in $(ls *_duplicate_count.txt);do (
		awk -F '~' '{a[$1]++;b[$1]=$0}END{for(x in a)if(a[x]==1)print b[x]}' $i | sort -k1,1 > ${i%_duplicate_count*}_multialign_phylum_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
	for i in $(ls *_multialign_phylum_reads.txt);do (
		awk -F '\t'  'FNR==NR {a[$1]; next}; $1 in a' $i ${i%_multialign_phylum_reads*}_phylum_inter2.txt | sort -u -k1,1 | awk 'gsub("~","\t",$0)'| awk -F '\t' '{print $3, $1, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21}' OFS='\t' > ${i%_multialign_phylum_reads*}_phylum_EMC.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			wait
		fi
	done
	wait
	rm *_phylum_inter.txt && rm *_phylum_inter2.txt && rm *_duplicate_count.txt && rm *_multialign_phylum_reads.txt && rm *_phylum_duplicates.txt
	for i in $(ls *_phylum_unique_reads.txt);do (
		awk -F '\t' '{print $11}' OFS=';' $i > ${i%_phylum_unique_reads*}_taxids_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done

	for i in $(ls *_uniq_inter.txt);do (
		awk -F ';' '{print $1}' OFS='\t' $i > ${i%_taxids_uniq_inter*}_taxids_uniq.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done

	rm *_taxids_uniq_inter.txt
	cd $projdir/metagenome/sighits/sighits_phylum
	for i in $(ls *_taxids_uniq.txt);do (
		awk -F '\t' 'NR==FNR{a[$1]=$0;next} ($1) in a{print a[$1]}' ${Qmatey_dir}/tools/rankedlineage_edited.dmp OFS='\t' $i> ${i%_taxids_uniq*}_uniq_inter.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
	for i in $(ls *_uniq_inter.txt);do (
		awk -F '\t'  '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' $i > ${i%_uniq_inter*}_species_taxid.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
	rm *_taxids_uniq.txt
	for i in $(ls *_species_taxid.txt);do (
		awk -F '\t' '{print $1}' $i | awk -F ' ' '{print $1, $2}' > ${i%_species_taxid*}_species_column.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	for i in $(ls *_species_column.txt);do (
		paste <(awk '{print $0}' OFS='\t' $i) <(awk -F '\t' '{print $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%_species_column*}_species_taxid.txt) | awk -F '\t' '{print $2, $1, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' > ${i%_species_column*}_species_taxa.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	for i in $(ls *_unique_reads.txt);do (
		paste <(awk -F '\t' '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' OFS='\t' $i ) <(awk -F '\t' '{print $2, $3, $4, $5, $6, $7, $8, $9, $10}' OFS='\t' ${i%*_phylum_unique_reads*}_species_taxa.txt) > ${i%_uniq*}_phylum_unique_uncultured.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	rm *_species_taxid.txt && rm *_uniq_inter.txt && rm *_species_column.txt && rm *_species_taxa.txt

	for i in $(ls *_phylum_unique_uncultured.txt);do (
		awk -F '\t' '!/Uncultured/' $i > ${i%*_phylum_unique_uncultured*}_phylum_unique_sequences.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done

	for i in $(ls *_phylum_EMC.txt);do (
		cat $i ${i%_phylum_EMC*}_phylum_phylum_unique_sequences.txt > ${i%_phylum_EMC*}_complete_phylum_reads.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	wait
	for i in $(ls *_complete_phylum_reads.txt);do (
		awk -F '\t' '{ for(i=1;i<=NF;i++){if(i==NF){printf("%s\n",$NF);}else {printf("%s\t",$i)}}}' $i > ${i%_complete_phylum_reads*}_sighits_temp.txt
		echo $'abundance\tqseqid\tsseqid\tlength\tmismatch\tevalue\tpident\tgapopen\tqseq\tsseq\tstaxids\tstitle\tspecies\tgenus\tfamily\torder\tclass\tphylum\tkingdom\tdomain' | \
		cat - ${i%_complete_phylum_reads*}_sighits_temp.txt > ${i%_complete_phylum_reads*}_sighits.txt )&
		if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			wait
		fi
	done
	rm *_complete_phylum_reads.txt && rm *_sighits_temp.txt && rm *_phylum_EMC.txt && rm *_unique_reads.txt && rm *_unique_sequences.txt && rm *_unique_uncultured.txt
fi

cd $projdir/metagenome/sighits/sighits_phylum
echo -e "${YELLOW}- compiling taxonomic information"
find . -type f -name '*_sighits.txt' -exec cat {} + > sighits.txt
awk -F '\t' '{print $18"\t"$19"\t"$20}' sighits.txt > rankedlineage_subhits_temp.txt && rm sighits.txt
awk -F '\t' '!/phylum/' rankedlineage_subhits_temp.txt > rankedlineage_subhits_temp2.txt && rm rankedlineage_subhits_temp.txt
echo $'phylum\tkingdom\tdomain' | cat - rankedlineage_subhits_temp2.txt > rankedlineage_subhits_temp3.txt && rm rankedlineage_subhits_temp2.txt
awk -F '\t' 'NR==1; NR >1 {print $0 | "sort -u"}' rankedlineage_subhits_temp3.txt > rankedlineage_subhits.txt && rm rankedlineage_subhits_temp3.txt

cd $projdir/metagenome/sighits/sighits_phylum
echo -e "${YELLOW}- quantifying the phylum-level taxonomy"
awk '{print $1}' rankedlineage_subhits.txt > phylum_taxa_mean_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > phylum_taxa_unique_sequences_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > phylum_taxa_quantification_accuracy_temp.txt
awk '{print $1}' rankedlineage_subhits.txt > phylum_taxa_rel_quantification_accuracy_temp.txt



phylum_level=phylum
for i in $(ls *_sighits.txt);do
	Rscript "${Qmatey_dir}/scripts/stats_summary.R" $i $min_unique_seqs $phylum_level "${Qmatey_dir}/tools/R"
	echo $'phylum\tmean\tuniq_reads\tstderr\trel_stderr' | cat - stats1.txt > stats2.txt
	id=${i%_sighits*}_mean && awk -v id=$id '{gsub(/mean/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$2}' > holdmean.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdmean.txt phylum_taxa_mean_temp.txt > holdmean2.txt && cat holdmean2.txt > phylum_taxa_mean_temp.txt
	id=${i%_sighits*}_uniq_reads && awk -v id=$id '{gsub(/uniq_reads/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$3}' > holduniq_reads.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holduniq_reads.txt phylum_taxa_unique_sequences_temp.txt > holduniq_reads2.txt && cat holduniq_reads2.txt > phylum_taxa_unique_sequences_temp.txt
	id=${i%_sighits*}_stderr && awk -v id=$id '{gsub(/stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$4}' > holdstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdstderr.txt phylum_taxa_quantification_accuracy_temp.txt > holdstderr2.txt && cat holdstderr2.txt > phylum_taxa_quantification_accuracy_temp.txt
	id=${i%_sighits*}_rel_stderr && awk -v id=$id '{gsub(/rel_stderr/,id); print }' stats2.txt | awk -F '\t' '{print $1,"\t",$5}' > holdrelstderr.txt
	awk 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){a[$1]=0}; print $0, a[$1]}'  holdrelstderr.txt phylum_taxa_rel_quantification_accuracy_temp.txt > holdrelstderr2.txt && cat holdrelstderr2.txt > phylum_taxa_rel_quantification_accuracy_temp.txt
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $1,"\t",$2,"\t",$3,"\t",$4,"\t",$5,"\t", a[$1]}'  rankedlineage_subhits.txt stats2.txt > stats3.txt
	awk '{print $1,$2,$3,$4,$5}' stats3.txt | awk '{gsub(/ /,"\t"); print }' > ${i%_sighits*}_taxastats.txt
	rm *stats1* *stats2* *stats3* *hold*
done

awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' phylum_taxa_mean_temp.txt > phylum_taxa_mean_temp2.txt && rm phylum_taxa_mean_temp.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_mean_temp2.txt > ../../results/phylum_level/phylum_taxa_mean.txt && rm phylum_taxa_mean_temp2.txt
awk 'NR==1; NR > 1 {s=0; for (i=2;i<=NF;i++) s+=$i; if (s!=0)print}' phylum_taxa_unique_sequences_temp.txt > phylum_taxa_unique_sequences_temp2.txt && rm phylum_taxa_unique_sequences_temp.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_unique_sequences_temp2.txt > ../../results/phylum_level/phylum_taxa_unique_sequences.txt && rm phylum_taxa_unique_sequences_temp2.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_quantification_accuracy_temp.txt > phylum_taxa_quantification_accuracy_temp2.txt && rm phylum_taxa_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/phylum_level/phylum_taxa_mean.txt phylum_taxa_quantification_accuracy_temp2.txt > phylum_taxa_quantification_accuracy_temp3.txt && rm phylum_taxa_quantification_accuracy_temp2.txt
sed '2,${/phylum/d;}' phylum_taxa_quantification_accuracy_temp3.txt > ../../results/phylum_level/phylum_taxa_quantification_accuracy.txt && rm phylum_taxa_quantification_accuracy_temp3.txt
awk '{gsub(/ /,"\t"); print}' phylum_taxa_rel_quantification_accuracy_temp.txt > phylum_taxa_rel_quantification_accuracy_temp2.txt && rm phylum_taxa_rel_quantification_accuracy_temp.txt
awk -F '\t' 'NR==FNR{c[$1]++;next};c[$1] > 0' ../../results/phylum_level/phylum_taxa_mean.txt phylum_taxa_rel_quantification_accuracy_temp2.txt > phylum_taxa_rel_quantification_accuracy_temp3.txt && rm phylum_taxa_rel_quantification_accuracy_temp2.txt
sed '2,${/phylum/d;}' phylum_taxa_rel_quantification_accuracy_temp3.txt > ../../results/phylum_level/phylum_taxa_rel_quantification_accuracy.txt && rm phylum_taxa_rel_quantification_accuracy_temp3.txt

cd $projdir/metagenome/results/phylum_level
for i in {mean,unique_sequences,quantification_accuracy,rel_quantification_accuracy}; do (
	awk 'NR==FNR{a[$1]=$0;next} ($1) in a{print $0, a[$1]}'  ../../sighits/sighits_phylum/rankedlineage_subhits.txt phylum_taxa_${i}.txt | \
	awk 'NR==1{for(i=1;i<=NF;i++)b[$i]++&&a[i]}{for(i in a)$i="";gsub(" +"," ")}1' | awk '{gsub(/ /,"\t"); print }' | awk '$0!=x;(NR==1){x=$0}' > phylum_taxainfo_${i}.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
		wait
	fi
done
rm *_taxa_*


cd $projdir/metagenome/results/phylum_level
awk '{print $1}' phylum_taxainfo_mean.txt > phylum_taxainfo_mean_holdingtaxid.txt
awk -F '\t' -v OFS='\t' 'FNR==1{for(i=1;i<=NF;++i)if($i!~/_mean/)k[i]=1}{n=split($0,f,FS);$0=j="";for(i=1;i<=n;++i)if(i in k)$(++j)=f[i]}1' \
phylum_taxainfo_mean.txt | awk '!($1="")' | awk 'BEGIN{OFS="\t"} {$1=$1};1' > phylum_taxainfo_mean_holdingtaxinfo.txt
touch phylum_taxainfo_mean_buildnorm.txt
for i in $(ls -S ../../../metagenome/haplotig/*_metagenome.fasta); do (
	sample=${i%_metagenome.fasta}; sample=${sample##*/}
	if [[ -z "$(ls -A ${projdir}/metagenome/coverage_normalization_factor.txt &> /dev/null)" ]]; then
		normfactor=1
	else
		normfactor=$( awk -v sample=$sample '$1 == sample' ../../coverage_normalization_factor.txt | awk '{print $2}' )
	fi
	awk -v sample=${sample}_mean -v norm=$normfactor 'BEGIN{OFS="\t"} NR==1 {for (i=1; i<=NF; i++) if ($i==sample) break} {print $i}' "$2" phylum_taxainfo_mean.txt | \
	paste - phylum_taxainfo_mean_buildnorm.txt > phylum_taxainfo_mean_buildnorm0.txt
	cat phylum_taxainfo_mean_buildnorm0.txt > phylum_taxainfo_mean_buildnorm.txt
	rm phylum_taxainfo_mean_buildnorm0.txt )&
	if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
		wait
	fi
done
paste phylum_taxainfo_mean_holdingtaxid.txt phylum_taxainfo_mean_buildnorm.txt > phylum_taxainfo_mean_norm0.txt
paste phylum_taxainfo_mean_norm0.txt phylum_taxainfo_mean_holdingtaxinfo.txt | awk '{gsub(/\t\t/,"\t"); print $0 }' > phylum_taxainfo_mean_normalized.txt
rm phylum_taxainfo_mean_holdingtaxid.txt phylum_taxainfo_mean_buildnorm.txt phylum_taxainfo_mean_holdingtaxinfo.txt phylum_taxainfo_mean_norm0.txt

for i in *.txt; do
	awk '{gsub(/-/,"_"); print}' $i > temp.txt
	mv temp.txt $i
done
rm temp.txt

file=${projdir}/exclude_taxa.txt
if test -f $file; then
	cat phylum_taxainfo_mean.txt > phylum_taxainfo_mean_filtered.txt
	cat phylum_taxainfo_unique_sequences.txt > phylum_taxainfo_unique_sequences_filtered.txt
	cat phylum_taxainfo_quantification_accuracy.txt > phylum_taxainfo_quantification_accuracy_filtered.txt
	cat phylum_taxainfo_rel_quantification_accuracy.txt > pylum_taxainfo_rel_quantification_accuracy_filtered.txt
	while read -r line; do
		for i in $( ls *filtered.txt ); do (
			awk -v line=$line '!/\tline\t/' $i > ${i%.txt}_temp.txt && mv ${i%.txt}_temp.txt $i )&
			if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
			  wait
			fi
		done
	done > $file
fi


if test -f $file; then
	echo -e "${YELLOW}- creating phylum-level visualizations"
	cd $projdir/metagenome/results/phylum_level
	phylum_level_mean=phylum_taxainfo_mean_filtered.txt
	phylum_level_uniq=phylum_taxainfo_unique_sequences_filtered.txt
	phylum_level_stderr=phylum_taxainfo_quantification_accuracy_filtered.txt
	phylum_level_rel_stderr=phylum_taxainfo_rel_quantification_accuracy_filtered.txt
	if [[ -z $min_percent_sample ]]; then
		min_percent_sample=5,10,20
	fi

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/phylum_level_corr.R" "$phylum_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/phylum_level_boxplots.R" "$phylum_level_mean" "$phylum_level_uniq" "$phylum_level_stderr" "$phylum_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/phylum_level/boxplots
	mv *.html $projdir/metagenome/results/phylum_level/boxplots
else
	echo -e "${YELLOW}- creating phylum-level visualizations"
	cd $projdir/metagenome/results/phylum_level
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

	if [[ "$total_no_samples" -ge 24 ]]; then
		for min_perc in $min_percent_sample; do (
				Rscript "${Qmatey_dir}/scripts/phylum_level_corr.R" "$phylum_level_mean" "$min_perc" "$min_pos_corr" "$max_neg_corr" "${Qmatey_dir}/tools/R" &>/dev/null
				)&
			 if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
				 wait
			 fi
		done
	fi
	Rscript "${Qmatey_dir}/scripts/phylum_level_boxplots.R" "$phylum_level_mean" "$phylum_level_uniq" "$phylum_level_stderr" "$phylum_level_rel_stderr" "$min_percent_sample" "${Qmatey_dir}/tools/R" &>/dev/null
	mkdir boxplots
	mv *_files $projdir/metagenome/results/phylum_level/boxplots
	mv *.html $projdir/metagenome/results/phylum_level/boxplots
fi
}
if [[ "$phylum_level" == true ]]; then
	time phylum 2>> $projdir/log.out
fi


######################################################################################################################################################
sunburst() {
echo -e "\e[97m########################################################\n \e[38;5;210mQmatey is creating sunburst \n\e[97m########################################################\n"

cd $projdir/metagenome/results

if [[ -z $sunburst_taxlevel ]]; then
	sunburst_taxlevel=strain
fi
sunburst_taxlevel=$( echo "$sunburst_taxlevel" | awk '{gsub(/,/,"\n")}1' )

for tsun in "$sunburst_taxlevel"; do
    cd $projdir/metagenome/results/${tsun}_level
    mean=${tsun}_taxainfo_mean.txt
    mean_norm=${tsun}_taxainfo_mean_normalized.txt
    if [[ -z $mean_norm ]]; then
      :
    else
      mean=$mean_norm
    fi
    for min_perc in $min_percent_sample; do (
      Rscript "${Qmatey_dir}/scripts/sunburst.R" "$mean" "$min_perc" "${sunburst_nlayers}" "${Qmatey_dir}/tools/R" &>/dev/null )&
			if [[ $(jobs -r -p | wc -l) -ge $gN ]]; then
			  wait
			fi
    done
done

}
cd $projdir
if [[ "$run_sunburst" == true ]]; then
	time sunburst &>> $projdir/log.out
fi

cd ${Qmatey_dir}/tools
gzip rankedlineage.dmp && rm rankedlineage_edited.dmp
if [[ "$normalization" == true ]]; then
	mv ${projdir}/metagenome ${projdir}/metagenome_ref_normalize
fi
if [[ "$normalization" == false ]]; then
	mv ${projdir}/metagenome ${projdir}/metagenome_no_normalize
fi

######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${white}\n"
cd $projdir
if [[ -z $samples_alt_dir ||  $samples_alt_dir == false ]]; then
 	:
else
	rm sammples
fi
echo -e "\n\n${magenta}- Run Complete ${white}"
