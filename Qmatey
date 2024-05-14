#!/bin/bash

magenta=$'\e[1;35m'
white=$'\e[0m'
yellow=$'\e[1;33m'
# Qmatey_dir=$(dirname "$0")
# projdir=$1
threads=""
# relpath=$(pwd)

export Qmatey_dir="$( cd -- "$(dirname "$0 ")" >/dev/null 2>&1 ; pwd -P )/"

if [[ "$1" == "" ]]; then
  echo -e "${white}\n\tUsage:\n\t\t./Qmatey/Qmatey\t\t<command>\n"
  echo -e "${yellow}\tQantitative Metagenomic Alignment and Taxonomic Exact-matching"
  echo -e "${yellow}\t\t- Exact-Matching of sequences at strain-level"
  echo -e "${yellow}\t\t- Exact-Matching of Consensus sequences (EMC) at higher taxonomic levels\n"
  echo -e "${white}\tCommand:"
  echo -e "${white}\t\t--version, -v\t\tprint software version"
  echo -e "${white}\t\t--help, -h\t\tprint help message"
  echo -e "${white}\t\tinstall\t\t\tsoftware dependencies"
  echo -e "${white}\t\tproj_dir\t\tspecify absolute or relative path to project directory\n"
  exit 0
fi
if [[ "$1" == "-v" || "$1" == "--version" ]]; then
  echo -e "${yellow}\n\tProgram:	Qmatey"
  echo -e "${yellow}\tVersion:	0.5.1"
  echo -e "${yellow}\tDescription:	Quantitative Metagenomic Alignment and Taxonomic Exact-matching"
  echo -e "${yellow}\tDevelopers:	Alison K. Adams (UTK, TN), Brandon Kristy (ORNL, TN), and Bode A. Olukolu (UTK, TN)\n"
  echo -e "${yellow}\tContact:	Bode Olukolu <bolukolu@utk.edu>\n"
  exit 0
fi
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  echo -e "${white}\n\tUsage:\n\t\t./Qmatey/Qmatey\t\t<command>\n"
  echo -e "${yellow}\tQuantitative Metagenomic Alignment and Taxonomic Exact-matching"
  echo -e "${yellow}\t\t- Exact-Matching at strain-level"
  echo -e "${yellow}\t\t- Exact-Matching of Consensus (EMC) at higher taxonomic levels\n"
  echo -e "${white}\tCommand:"
  echo -e "${white}\t\t--version, -v\t\tprint software version"
  echo -e "${white}\t\t--help, -h\t\tprint help message"
  echo -e "${white}\t\tinstall\t\t\tsoftware dependencies"
  echo -e "${white}\t\tproj_dir\t\tspecify absolute or relative path to project directory\n"
  echo -e "${magenta}- Do you want to view pipeline help? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ $REPLY =~ ^[Yy]$ ]]; then
    awk '{gsub(/\|/," | ");}1' $Qmatey_dir/misc/pipeline_help.txt | less
  fi
  exit 0
fi
if [[ "$1" == "install" ]]; then
  mkdir -p $Qmatey_dir/tools
  cd $Qmatey_dir/tools
  bash $Qmatey_dir/scripts/install.sh
  wait; exit 0
fi

projdir="$( cd -- "$(dirname "$1 ")" >/dev/null 2>&1 ; pwd -P )/"
export projdir="$( cd -- "$(dirname "$1 ")" >/dev/null 2>&1 ; pwd -P )/"
if [[ ! -d "$Qmatey_dir/tools" ]]; then
  mkdir -p $Qmatey_dir/tools
  cd $Qmatey_dir/tools
  bash $Qmatey_dir/scripts/install.sh
fi
if [[ -d "$Qmatey_dir/tools" ]] && [[ -z $(ls -A "$Qmatey_dir/tools") ]] ; then
  mkdir -p $Qmatey_dir/tools
  cd $Qmatey_dir/tools
  bash $Qmatey_dir/scripts/install.sh
fi

cd $Qmatey_dir/tools
if [ ! -d ./R/CCLasso ]; then echo -e "${magenta}- R-package  is not installed or is not properly installed\n- Consider checking dependencies and re-running install command\n\n ${white}"; fi



# if [ "${Qmatey_dir:0:1}" = "." ]; then
#   if [ "${Qmatey_dir:0:15}" = "../../../../../" ]; then
#     Qmatey_dir="${relpath%/*/*/*/*/*}${Qmatey_dir//*..}"
#   fi
#   if [ "${Qmatey_dir:0:12}" = "../../../../" ]; then
#     Qmatey_dir="${relpath%/*/*/*/*}${Qmatey_dir//*..}"
#   fi
#   if [ "${Qmatey_dir:0:9}" = "../../../" ]; then
#     Qmatey_dir="${relpath%/*/*/*}${Qmatey_dir//*..}"
#   fi
#   if [ "${Qmatey_dir:0:6}" = "../../" ]; then
#     Qmatey_dir="${relpath%/*/*}${Qmatey_dir//*..}"
#   fi
#   if [ "${Qmatey_dir:0:3}" = "../" ]; then
#     Qmatey_dir="${relpath%/*}${Qmatey_dir//*..}"
#   fi
#   if [ "${Qmatey_dir:0:2}" != ".." ]; then
#     if [ "${Qmatey_dir:0:1}" = "." ]; then
#       Qmatey_dir="${relpath}${Qmatey_dir:1}"
#     fi
#   fi
# fi
# if [ "${Qmatey_dir: -1}" != "/" ]; then
#   Qmatey_dir="${Qmatey_dir}/"
# fi
#
# if [ "${projdir:0:1}" = "." ]; then
#   if [ "${projdir:0:15}" = "../../../../../" ]; then
#     projdir="${relpath%/*/*/*/*/*}${projdir//*..}"
#   fi
#   if [ "${projdir:0:12}" = "../../../../" ]; then
#    projdir="${relpath%/*/*/*/*}${projdir//*..}"
#  fi
#   if [ "${projdir:0:9}" = "../../../" ]; then
#     projdir="${relpath%/*/*/*}${projdir//*..}"
#   fi
#   if [ "${projdir:0:6}" = "../../" ]; then
#     projdir="{relpath%/*/*}${projdir//*..}"
#   fi
#   if [ "${projdir:0:3}" = "../" ]; then
#     projdir="${relpath%/*}${projdir//*..}"
#   fi
#   if [ "${projdir:0:2}" = "./" ]; then
#     projdir="${relpath}${projdir:1}"
#   fi
#   if [ "${projdir}" = . ]; then
#     projdir="${relpath}"
#   fi
# fi
# if [ "${projdir:-1}" != "/" ]; then
#   projdir="${projdir}/"
# fi

####################################################################################################################
####################################################################################################################
####################################################################################################################

cd $projdir
echo -e "${white}\n#########################################################################################"
echo -e "${yellow}\t- Program:\tQmatey"
echo -e "${yellow}\t- Version:\t0.5.0"
echo -e "${yellow}\t- Description:\tQuantitative Metagenomic Alignment and Taxonomic Exact-Matching"
echo -e "${yellow}\t\t\t * Exact-Matching of sequences at strain-level"
echo -e "${yellow}\t\t\t * Exact-Matching of Consensus sequences (EMC) at higher taxonomic levels\n"
echo -e "${yellow}\t- Developers:\tAlison K. Adams (UTK, TN), Brandon Kristy (ORNL, TN), and Bode A. Olukolu (UTK, TN)\n"
echo -e "${yellow}\t- Contact:\tBode Olukolu <bolukolu@utk.edu>"
echo -e "${white}\n#########################################################################################"


# d=($(find . -maxdepth 1 -type d | wc -l))
# if [ $d == 3 ]; then
#   echo ""
# else
#   echo -e "${magenta}- Expecting at least 2 folders (i.e. norm_ref and samples) ${white}"
#   echo -e "${magenta}- Do you wish to continue running Qmatey? ${white}"
#   read -p "- y(YES) or n(NO)? " -n 1 -r
#   if [[ ! $REPLY =~ ^[Yy]$ ]]; then
#     exit 1
#   else
#       echo -e "${white}\n#########################################################################################"
#   fi
# fi


tr -d '\15\32' < config.sh > temp; mv temp config.sh
f=($(find . -maxdepth 1 -type f | wc -l))
if [ $f == 1 ]; then
  echo ""
else
  echo -e "${magenta}- Expecting only 1 file (i.e. config.sh) ${white}"
  echo -e "${magenta}- Do you wish to continue running Qmatey? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  else
      echo -e "${white}\n#########################################################################################"
  fi
fi

####################################################################################################################
####################################################################################################################
####################################################################################################################
cd $projdir
samples_alt_dir="$(grep samples_alt_dir config.sh)"
samples_alt_dir=${samples_alt_dir//*=}
if [[ -z $samples_alt_dir ]]; then
  samples_alt_dir==false
fi
mkdir -p samples
if [[ $samples_alt_dir == false ]]; then
	:
else
	if [ -d ./samples ] && [ "$(ls -A ./samples 2> /dev/null | wc -l)" -gt 0 ]; then
		echo -e "${magenta}- sample files exist in ./project_directory/samples/ ${white}\n"
		echo -e "${magenta}- cannot use both <alternate_samples_directory> and <./project_directory/samples/ ${white}\n"
		echo -e "${magenta}- Qmatey will quit in 5 seconds ${white}\n"
		sleep 5; exit 0
	fi
	if [ ! -d ./samples ]; then
		ln -s $samples_alt_dir .
	fi
	if [ -d ./samples ] && [ "$(ls -A ./samples 2> /dev/null | wc -l)" -eq 0 ]; then
		rmdir samples
		ln -s $samples_alt_dir .
	fi
fi

taxids="$(grep taxids config.sh | grep -v "map")"
blast_location="$(grep blast_location config.sh)"
blast_location=${blast_location//*=}
taxids=${taxids//*=}
if [[ -z $taxids ]]; then
	$taxids=true
fi
if [[ "$blast_location" != "custom" ]]; then
  if [[ $taxids == true ]] && [[ ! -d ./taxids ]]; then
  	echo -e "${magenta}- taxids directory does not exist ${white}\n"
  	echo -e "${magenta}- to subset MegaBLAST with list of taxids, provide files containing list of tax_ids ${white}\n"
  	echo -e "${magenta}- else, set taxids=false in config.sh file ${white}\n"
  	echo -e "${magenta}- Qmatey will quit in 5 seconds ${white}\n"
  	sleep 5; exit 0
  fi
  if [[ $taxids == true ]] && [[ -d ./taxids ]]; then
  	if [[ -z "$(ls -A ./taxids)" ]]; then
  		echo -e "${magenta}- taxids directory is empty ${white}\n"
  		echo -e "${magenta}- to subset MegaBLAST with list of tax_ids, provide files containing list of tax_ids ${white}\n"
  		echo -e "${magenta}- else, set taxids=false in config.sh file ${white}\n"
  		echo -e "${magenta}- Qmatey will quit in 5 seconds ${white}\n"
  		sleep 5; exit 0
  	fi
  fi
fi

cd samples
if [ -d "paired" ]; then
 mv paired pe
fi
if [ -d "single" ]; then
 mv single se
fi




# check rankedlineage. Perform update?
cd $Qmatey_dir/tools
file=rankedlineage.dmp
if test -f $file; then
  echo -e "${magenta}- we recommend updating the ranked lineage file every few weeks/months. ${white}"
  echo -e "${magenta}- warning: Updating the ranked lineage file may significantly alter the consistency of your results ${white}"
  echo -e "${magenta}- would you like to update the ranked lineage file? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	:
  else
    lineage_update() {
      rm rankedlineage.dmp
      wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
      gunzip -c new_taxdump.tar.gz | tar xf -
      rm *tar.gz
      rm citations.dmp; rm delnodes.dmp; rm division.dmp; rm fullnamelineage.dmp; rm gencode.dmp; rm host.dmp; rm merged.dmp; rm names.dmp; rm nodes.dmp
      rm taxidlineage.dmp; rm typematerial.dmp; rm typeoftype.dmp
    }
    echo -e "${white}- updating ranked lineage taxa dump ${white}"
    lineage_update &>> ./log.out
  fi
else
  lineage_install() {
    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
    tar -zxvf new_taxdump.tar.gz; rm *tar.gz
    rm citations.dmp; rm delnodes.dmp; rm division.dmp; rm fullnamelineage.dmp; rm gencode.dmp; rm host.dmp; rm merged.dmp; rm names.dmp; rm nodes.dmp
    rm taxidlineage.dmp; rm typematerial.dmp; rm typeoftype.dmp
  }
  echo -e "${white}- Performing installation of ranked lineage taxa dump ${white}"
  lineage_install &>> ./log.out
fi


# download and setup NCBI database. Perform update?
cd $projdir
awk '{gsub(/ /,""); print}' config.sh > temp
mv temp config.sh
db="$(grep -v =NA config.sh | grep -v '#' | grep _db=)"
db=${db//*=}
db=$(echo $db | awk '{gsub(/,/," "); print}')
input_db=${db%/*}

setup_database() {
  cd $input_db
  if [[ "${db/*\/}" == nt ]]; then
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*tar.gz
    for f in nt*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  if [[ "${db/*\/}" =~ refseq ]]; then
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_viroids_rep_genomes*tar.gz
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_prok_rep_genomes*tar.gz
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_euk_rep_genomes*tar.gz
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_rna*tar.gz
    for f in ref*.tar.gz; do tar -xzf "$f"; done
    rm ref*.tar.gz
  fi
  if [[ "${db/*\/}" =~ 16S ]]; then
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA*tar.gz
    for f in 16S_ribosomal_RNA*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  if [[ "${db/*\/}" =~ 18S ]]; then
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/18S_fungal_sequences*tar.gz
    for f in 18S_fungal_sequences*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  if [[ "${db/*\/}" =~ 28S ]]; then
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/28S_fungal_sequences*tar.gz
    for f in 28S_fungal_sequences*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  if [[ "${db/*\/}" =~ ITS ]]; then
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/ITS_eukaryote_sequences*tar.gz
    for f in ITS_eukaryote_sequences*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  cd $projdir
}

update_database() {
  cd $input_db
  if [[ "${db/*\/}" =~ refseq ]]; then
    rm -r ref* taxdb.btd taxdb.bti index.html 2> /dev/null
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_viroids_rep_genomes*tar.gz
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_prok_rep_genomes*tar.gz
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/ref_euk_rep_genomes*tar.gz
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/refseq_rna*tar.gz
    for f in ref*.tar.gz; do tar -xzf "$f"; done
    rm ref*.tar.gz
  fi
  if [[ "${db/*\/}" == nt ]]; then
    rm -r nt* taxdb.btd taxdb.bti index.html 2> /dev/null
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*tar.gz
    for f in nt*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  if [[ "${db/*\/}" =~ 16S ]]; then
    rm -r 16S_ribosomal_RNA* taxdb.btd taxdb.bti index.html 2> /dev/null
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA*tar.gz
    for f in 16S_ribosomal_RNA*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  if [[ "${db/*\/}" =~ 18S ]]; then
    rm -r 18S_fungal_sequences* taxdb.btd taxdb.bti index.html 2> /dev/null
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/18S_fungal_sequences*tar.gz
    for f in 18S_fungal_sequences*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  if [[ "${db/*\/}" =~ 28S ]]; then
    rm -r 28S_fungal_sequences* taxdb.btd taxdb.bti index.html 2> /dev/null
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/28S_fungal_sequences*tar.gz
    for f in 28S_fungal_sequences*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  if [[ "${db/*\/}" =~ ITS ]]; then
    rm -r ITS_eukaryote_sequences* taxdb.btd taxdb.bti index.html 2> /dev/null
    wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/ITS_eukaryote_sequences*tar.gz
    for f in ITS_eukaryote_sequences*.tar.gz; do tar -xzf "$f"; done
    rm nt*.tar.gz
  fi
  cd $projdir
}

if [[ "${db}" == "${db% *}" ]] && [[ ! -z ${db} ]]; then
  input_db=${db%/*}
  mkdir -p "$input_db"
  if [[ $(ls -A ${input_db}/* 2> /dev/null | wc -l ) -gt 0 ]]; then
    echo -e "${magenta}- $input_db database already exist. ${white}"
    echo -e "${magenta}- it is recommended to update NCBI database every few months. ${white}"
    echo -e "${magenta}- would you like to update your NCBI "${input_db/*\/}" database? ${white}"
    read -p "- y(YES) or n(NO) " -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      :
    else
      printf "\n"
      echo -e "${magenta}- downloading NCBI database(s) might take a few hours? ${white}"
      echo -e "${white}- update of NCBI "${input_db/*\/}" database will start in 60 secs. <Ctrl+C> to abort update ${white}"
      sleep 60
      wait
      echo -e "${white}- updating your NCBI "${input_db/*\/}" database ${white}"
      update_database &>> ./db_download.log
    fi
  else
    echo -e "${white}- user-defined database does not contain "${input_db/*\/}" database ${white}"
    echo -e "${magenta}- would you like to download the NCBI "${input_db/*\/}" database in $db? ${white}"
    read -p "- y(YES) or n(NO) " -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      :
    else
      printf "\n"
      echo -e "${magenta}- downloading NCBI database(s) might take a few hours? ${white}"
      echo -e "${white}- download of NCBI "${input_db/*\/}" database will start in 60 secs. <Ctrl+C> to abort update ${white}"
      sleep 60
      wait
      echo -e "${white}- downloading NCBI "${input_db/*\/}" database in user-defined directory ${white}"
      setup_database &>> ./db_download.log
    fi
  fi
else
  declare -a arr=($db)
  for i in "${arr[@]}"; do
    j=${i%/*}
    mkdir -p "$j"
  done
  for i in "${arr[@]}"; do
    j=${i%/*}
    if [[ $(ls -A ${j}/* 2> /dev/null | wc -l ) -eq 0 ]]; then
      db_miss=yes
    fi
  done
  if [[ "$db_miss" != yes ]]; then
    echo -e "${magenta}- all listed database already exist. ${white}"
    echo -e "${magenta}- it is recommended to update NCBI database every few weeks/months. ${white}"
    echo -e "${magenta}- would you like to update the listed NCBI database? ${white}"
    read -p "- y(YES) or n(NO) " -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      :
    else
      for i in "${arr[@]}"; do
        input_db=${i%/*}
        echo -e "${white}- update of NCBI "${input_db/*\/}" database will start in 60 secs. <Ctrl+C> to abort update ${white}"
        sleep 60
        wait
        echo -e "${white}- updating your NCBI "${input_db/*\/}" database ${white}"
        update_database &>> ./db_download.log
      done
    fi
  else
    echo -e "${white}- user-defined database does not contain database ${white}"
    echo -e "${magenta}- would you like to download the NCBI database? ${white}"
    read -p "- y(YES) or n(NO) " -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      :
    else
      printf "\n"
      for i in "${arr[@]}"; do
        input_db=${i%/*}
        echo -e "${white}- download of NCBI "${input_db/*\/}" database will start in 60 secs. <Ctrl+C> to abort update ${white}"
        sleep 60
        wait
        echo -e "${white}- downloading NCBI "${input_db/*\/}" database in user-defined directory ${white}"
        setup_database &>> ./db_download.log
      done
    fi
  fi
fi


cd $projdir
if [[ -d "./samples/pe" ]] || [[ -d "./samples/se" ]]; then
  if test -f "./samples/filename_reformatted.txt" || test -f "./samples/flushed_reads.txt"; then
    :
  else
    echo -e "${magenta}- Qmatey will modify your sample fastq files by concatenating <pe-R1-reads>, <pe-R2-reads>, and <se-reads> ${white}"
    echo -e "${magenta}- Do you want create a copy/backup of your sample fastq file? ${white}"
    read -p "- y(YES) or n(NO)? " -n 1 -r
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      printf "\n"
      echo -e "${magenta}- A copy of your sample fastq files is in ${projdir}/samples_backup/ ${white}"
      cp -r samples samples_backup
    else
        echo -e "${white}\n#########################################################################################"
    fi
  fi
else
  if test -f "./samples/filename_reformatted.txt" || test -f "./samples/flushed_reads.txt"; then
    :
  else
    echo -e "${magenta}- Qmatey will modify your sample fastq files ${white}"
    echo -e "${magenta}- Do you want create a copy/backup of your sample fastq file? ${white}"
    read -p "- y(YES) or n(NO)? " -n 1 -r
    if [[ $REPLY =~ ^[Yy]$ ]]; then
      printf "\n"
      echo -e "${magenta}- A copy of your sample fastq files is in ${projdir}/samples_backup/ ${white}"
      cp -r samples samples_backup
    else
        echo -e "${white}\n#########################################################################################"
    fi
  fi
fi


cd ${projdir}
normalization="$(grep normalization config.sh)"
normalization=${normalization//*=}
if [[ -z "${normalization}" ]]; then
  normalization=true
fi
if [[ "${normalization}" == true ]]; then
  echo -e "${yellow} \e[31m Options for data normalization to obtain relative abundance/expression:"
  echo -e "${yellow} \e[31m \tMethod-1: Use spiked-in standard total read coverage (synthetic oligo/coomunity)"
  echo -e "${yellow} \e[31m \tMethod-2: Use host genome total read coverage (only for endophytic community and/or fairly consistent host genomc DNA contribution)"
  echo -e "${yellow} \e[31m \tMethod-3: Use sample total read coverage (relative abundance between samples)"

  echo -e "${magenta}- would you like to use Method 1 or 2 (i.e spike-in or host genome)? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${magenta}- \tUsing sample total read coverage for normalization (Method-3) ${white}"
    echo export norm_method=samples > norm_method.sh
  else
    echo -e "${magenta}- \tUsing spike-in/host genome (Method-1/2) ${white}"
    echo export norm_method=spike_host > norm_method.sh
  fi
fi
if [[ "${normalization}" == false ]]; then
  echo -e "${yellow} \e[31m Data normalization will not be performed i.e. output will be only absolute abundance/expression"
  echo export norm_method=false > norm_method.sh
fi

mkdir -p norm_ref
if [[ "$(ls -A ./norm_ref/ 2> /dev/null | wc -l)" -eq 0 ]]; then
  if [[ "${normalization}" == true ]]; then
    echo -e "${magenta} \e[31m normalization reference folder is empty (required for normalization methods 1 and 2)"
    echo -e "${magenta} \e[31m also, Qmatey will not subtract/eliminate reads derived from host genome/transciptome\n\n"
    echo -e "${magenta} \e[31m Qmatey can explore an alternative approach (below) for normalization:\n"
    echo -e "${magenta}- would you like to use read coverage of samples for normalization (method-3, assumes normalized conc. of input DNA) ${white}"
    read -p "- y(YES) or n(NO) " -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      echo -e "${magenta}- normalization will not be performed (skipping normalization) ${white}"
      awk '{gsub(/normalization=true/,"normalization=false"); print}' config.sh > temp.txt
      mv temp.txt config.sh
    else
      echo -e "${magenta}- normalization will be perfromed based on read coverage of samples (method-3, assumes normalized conc. of input DNA) ${white}"
      awk '{gsub(/normalization=false/,"normalization=true"); print}' config.sh > temp.txt
      mv temp.txt config.sh
    fi
  fi
fi


####################################################################################################################
# check if the job was previously submitted and interrupted
cd ${projdir}
nthreads="$(grep threads config.sh)"
nthreads=${nthreads//*=}
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


file=${projdir}/metagenome*/haplotig/combined_compressed_metagenomes.fasta
if test -f $file; then
  echo -e "${magenta}- job was interrupted post-MegaBLAST,and while extracting alignments for each sample ${white}"
  echo -e "${magenta}- Do you want to continue from where this analytical step was interrupted? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${magenta}- Qmatey will delete all sample megablast output files, and will redo analysis ${white}"
    rm ${projdir}/metagenome*/alignment/*_haplotig.megablast
  else
    echo -e "${magenta}- Qmatey will continue analysis, while removing interrupted files ${white}"
    while [ $gN -gt 0 ]; do
      rm "$(ls -lht * | head -n 1 | awk '{print $NF}')"
      gN=$(($gN-1))
    done
  fi
fi

if [[ -d ${projdir}/metagenome*/haplotig/splitccf ]]; then
  echo -e "${magenta}- job interrupted during MegaBLAST alignment ${white}"
  echo -e "${magenta}- Do you want to continue from where this analytical step was interrupted? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${magenta}- Qmatey will delete all alignments outputs, and will redo analysis ${white}"
    rm -rf ${projdir}/metagenome*/splitccf
    rm ${projdir}/metagenome*/alignment/*
  else
    echo -e "${magenta}- Qmatey will continue analysis, while removing interrupted files ${white}"
    if test ! -f ${projdir}/metagenome*/alignment/F*.blast; then
      rm subfile*; mv F* ../haplotig/splitccf
    fi
  fi
fi

####################################################################################################################
####################################################################################################################

main() {
string="${Qmatey_dir}/scripts/Qmatey_internal_parameters.sh"
string2=${string//'//'/'/'}
string="${Qmatey_dir}/scripts/Qmatey_profiling.sh"
string3=${string//'//'/'/'}

cd $projdir
awk '{ sub("\r$",""); print}' config.sh > unixformat.sh
mv unixformat.sh config.sh

cd $projdir
printf "#""!""/bin/bash \n\n" > header.txt
printf "Qmatey_dir=${Qmatey_dir%}\n" > fetchdir.txt
printf "projdir=${projdir%}" >> fetchdir.txt

string="$(grep walkaway config.sh)"
string=${string//*=}
cluster="$(grep cluster config.sh)"
cluster=${cluster//*=}
nthread="$(grep threads config.sh)"
nthread=${nthread//*=}
simulation_lib="$(grep simulation_lib config.sh)"
simulation_lib=${simulation_lib//*=}
simulation_motif="$(grep simulation_motif config.sh)"
simulation_motif=${simulation_motif//*=}


if [ -z "$string" ]; then
 string=true
fi
if [ "$cluster"  == false ]; then
 unset cluster
fi


if [ -d "${projdir}/simulate_genomes/" ]; then
	echo -e "${magenta}- Do you want to perform simulation based on WGS data in directory < ${projdir}simulate_genomes >? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping simulation ${white}\n"
    echo "simulate_reads=0" | cat config.sh - > unixformat.sh
    mv unixformat.sh config.sh
	else
		printf '\n'
		if [[ -z "$simulation_lib" ]] || [[ -z "$simulation_motif" ]]; then
			echo -e "${magenta}- \n- please specify simulation parameters. simulation_lib and simulation_motif are required ${white}\n"
			echo -e "${magenta}- \n- Qmatey will exit in 5 seconds ${white}\n"
			sleep 5 && exit 1
		fi
		if [[ "$(ls ./samples/* 2> /dev/null | wc -l)" -gt 0 ]] || [[ ! -d "simulate_genomes" ]]; then
			echo -e "${magenta}- \n- samples folder needs to be empty and WGS datasets needs to be provided in directory <${projdir}/simulate_genomes> ${white}\n"
			echo -e "${magenta}- \n- Qmatey will exit in 5 seconds ${white}\n"
			sleep 5 && exit 1
		fi
    echo "simulate_reads=1" | cat config.sh - > unixformat.sh
    mv unixformat.sh config.sh
	fi
else
  echo "simulate_reads=0" | cat config.sh - > unixformat.sh
  mv unixformat.sh config.sh
fi




string_norm="$(grep normalization config.sh)"
string_norm=${string_norm//*=}
if [ -z "$string_norm"  ]; then
	echo -e "${magenta}- Expecting variable normalization: true or false> ${white}"
	echo -e "${magenta}- Do you wish to perform normalization? ${white}"
	read -p "- y(YES) or n(NO)? " -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf "\n"
		echo "normalization=false" | cat config.sh - > unixformat.sh
		mv unixformat.sh config.sh
	else
		printf "\n"
		echo "normalization=true" | cat config.sh - > unixformat.sh
		mv unixformat.sh config.sh
	fi
fi

string_sunburst="$(grep run_sunburst config.sh)"
string_sunburst=${string_sunburst//*=}
if [ -z "$string_norm"  ]; then
	echo -e "${magenta}- Expecting variable run_sunburst (visualization of diversity & relative abundance): true or false> ${white}"
	echo -e "${magenta}- Do you wish to generate visualization of diversity & relative abundance? ${white}"
	read -p "- y(YES) or n(NO)? " -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf "\n"
		echo "string_sunburst=false" | cat config.sh - > unixformat.sh
		mv unixformat.sh config.sh
	else
		printf "\n"
		echo "string_sunburst=true" | cat config.sh - > unixformat.sh
		mv unixformat.sh config.sh
	fi
fi


if [[ "$string" == false ]]; then
  if [ -z "$nthread" ]; then
    echo -e "${magenta}- Qmatey will use (total processors/cores)-2? ${white}"
    echo -e "${magenta}- Quit (ctrl + c), else Qmatey will continue in 5 seconds ${white}"
    sleep 5
  fi
  cat header.txt <(grep -v '^[[:space:]]*$' config.sh | awk '{gsub(/^/,"export "); gsub(/export #/,"#");}1') fetchdir.txt $string2 $string3 | awk '{ sub("\r$",""); print}' > Qmatey_run.sh
  rm header.txt fetchdir.txt

  echo -e "${magenta}- Qmatey is ready to submit job ${white}"
  echo -e "${magenta}- Do you want to continue? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${magenta}\n- shell script will be saved to ${projdir}Qmatey_run.sh ${white}\n"
    exit 1
  else
    bash ${projdir}Qmatey_run.sh
  fi
fi





if [[ "$string" == true ]]; then
if [ -z $cluster ]; then
  if [ -z "$nthread" ]; then
    echo -e "${magenta}- Qmatey will use (total processors/cores)-2? ${white}"
    echo -e "${magenta}- Quit (ctrl + c), else Qmatey will continue in 5 seconds ${white}"
    sleep 5
  fi
  printf "#""!""/bin/bash \n" > header.txt
  printf "walkaway=true \n\n" > walkaway.txt
  printf "Qmatey_dir=${Qmatey_dir}\n" > fetchdir.txt
  printf "projdir=${projdir}" >> fetchdir.txt
  cat header.txt <(grep -v '^[[:space:]]*$' config.sh | awk '{gsub(/^/,"export "); gsub(/export #/,"#");}1') norm_method.sh walkaway.txt fetchdir.txt $string2 $string3 | awk '{ sub("\r$",""); print}' > Qmatey_run.sh
  rm header.txt walkaway.txt fetchdir.txt norm_method.sh

  echo -e "${magenta}- Qmatey is ready to submit job ${white}"
  echo -e "${magenta}- Do you want to continue? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${magenta}\n- shell script will be saved to ${projdir}Qmatey_run.sh ${white}\n"
    exit 1
  else
    nohup bash ${projdir}Qmatey_run.sh > terminal.out 2>&1 &
  fi
fi
fi





if [[ "$string" == true ]]; then
if [[ "$cluster" == true ]]; then
  if [ -z "$nthread" ]; then
    echo -e "${magenta}- Please provide number of threads config.sh (for cluster node)? ${white}"
    echo -e "${magenta}- Qmatey will quit in 5 seconds ${white}"
    sleep 5
    exit 1
  fi
  printf "#""!""/bin/bash \n#SBATCH -c ${nthread} \n\n" > cluster_header.sh
  printf "walkaway=true \n\n" > walkaway.txt
  cat cluster_header.sh <(grep -v '^[[:space:]]*$' config.sh | awk '{gsub(/^/,"export "); gsub(/export #/,"#");}1') norm_method.sh walkaway.txt fetchdir.txt $string2 $string3 | awk '{ sub("\r$",""); print}' > Qmatey_run.sh
  rm fetchdir.txt cluster_header.sh header.txt walkaway.txt norm_method.sh

  echo -e "${magenta}- Qmatey is ready to submit job ${white}"
  echo -e "${magenta}- Do you want to continue? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo -e "${magenta}\n- shell script will be saved to ${projdir}Qmatey_run.sh ${white}\n"
    exit 1
  else
    sbatch ${projdir}Qmatey_run.sh
  fi

  cd $projdir
  n_nodes="$(grep nodes config.sh)"
  n_nodes=${n_nodes//*=}
  if [ -z "$n_nodes" ]; then
   n_nodes=1
  fi
  if [[ "$n_nodes" -gt 1 ]]; then
    while [[ ! -f "${projdir}/multi_node_run_ready.txt" ]]; do sleep 300; done
    wait
    for nn in $(seq 2 "$n_nodes"); do
      printf "njob=${nn}\n\n" | cat <(awk '/extract_parameters_above_for_megablast_splitrun_multi_node_mode/{exit} 1' ${projdir}/Qmatey_run.sh ) - | \
      cat - ${Qmatey_dir}/scripts/megablast_split_node.sh > ${projdir}/megablast_splitrun_node_${nn}.sh &&
      sbatch ${projdir}/megablast_splitrun_node_${nn}.sh
    done
    wait
  fi
fi
fi
}
cd $projdir
main
echo -e "${magenta}- For jobs running in background, monitor progress in terminal.out/slurm-xxxxxx.out ${white}"
echo -e "${magenta}- The log.out file can help with troubleshooting and when reporting a bug ${white}"
