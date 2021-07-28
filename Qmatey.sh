#!/bin/bash

magenta=$'\e[1;35m'
white=$'\e[0m'
yellow=$'\e[1;33m'
# Qmatey_dir=$(dirname "$0")
# projdir=$1
threads=""
# relpath=$(pwd)

export Qmatey_dir="$( cd -- "$(dirname "$0 ")" >/dev/null 2>&1 ; pwd -P )/"
export projdir="$( cd -- "$(dirname "$1 ")" >/dev/null 2>&1 ; pwd -P )/"

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
echo -e "${yellow}- Program:	Qmatey"
echo -e "${yellow}- Version:	0.2.5"
echo -e "${yellow}- Description:	Quantitative Metagenomic Alignment and Taxonomic Exact-matching"
echo -e "${yellow}- Contact:	Bode Olukolu <bolukolu@utk.edu>; Brandon Kristy <bkristy@vols.utk.edu>"
echo -e "${white}\n#########################################################################################"


d=($(find . -maxdepth 1 -type d | wc -l))
if [ $d == 3 ]; then
  echo ""
else
  echo -e "${magenta}- Expecting only 2 folders (i.e. norm_ref and samples) ${white}"
  echo -e "${magenta}- Do you wish to continue running Qmatey? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    exit 1
  else
      echo -e "${white}\n#########################################################################################"
  fi
fi

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
cd samples
if [ -d "paired" ]; then
 mv paired pe
fi
if [ -d "single" ]; then
 mv single se
fi


# Initial questions before running walkaway
###########################################
cd ${Qmatey_dir}/tools/

# check rankedlineage. Perform update?
file=rankedlineage.dmp*
if test -f $file; then
  echo -e "${magenta}- we recommend updating the ranked lineage file every few months. ${white}"
  echo -e "${magenta}- warning: Updating the ranked lineage file may significantly alter the consistency of your results ${white}"
  echo -e "${magenta}- would you like to update the ranked lineage file? ${white}"
  read -p "- y(YES) or n(NO) " -n 1 -r
  if [[ ! $REPLY =~ ^[Yy]$ ]]; then
  	:
  else
    lineage_update() {
      rm rankedlineage.dmp.gz && rm rankedlineage.dmp
      wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
      tar -zxvf new_taxdump.tar.gz; rm *tar.gz
      rm citations.dmp; rm delnodes.dmp; rm division.dmp; rm fullnamelineage.dmp; rm gencode.dmp; rm host.dmp; rm merged.dmp; rm names.dmp; rm nodes.dmp
      rm taxidlineage.dmp; rm typematerial.dmp; rm typeoftype.dmp
      gzip rankedlineage.dmp
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
    gzip rankedlineage.dmp
  }
  echo -e "${white}- Performing installation of ranked lineage taxa dump ${white}"
  lineage_install &>> ./log.out
fi

# download and setup NCBI database. Perform update?
cd $projdir
db="$(grep -v NA config.sh | grep _db=)"
db=${db//*=}
mkdir -p "$db"
setup_database() {
  wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*.tar.gz
  for f in nt*.tar.gz; do tar -xzf "$f"; done
  rm nt*.tar.gz
}
update_database() {
  rm -r *
  wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*.tar.gz
  for f in nt*.tar.gz; do tar -xzf "$f"; done
  rm nt*.tar.gz
}


if [ -d "$db" ]; then
  echo "db directory already exists"
  if [ "$(ls -A "$db")" ]; then
    echo -e "${magenta}- if using the NCBI database, we recommend updating it every few months. ${white}"
    echo -e "${magenta}- warning: updating might take a few hours and may significantly alter the consistency of your results ${white}"
    echo -e "${magenta}- would you like to update your NCBI nt database? ${white}"
    read -p "- y(YES) or n(NO) " -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      :
    else
      echo -e "${white}- updating your NCBI nt database ${white}"
      update_database &>> ./log.out
    fi
  else
    echo -e "${white}- user-defined database is empty ${white}"
    echo -e "${white}- downloading and setting up NCBI nt database in user-defined directory ${white}"
    cd $db
    setup_database &>> ./log.out
  fi
else
  echo -e "${white}- user-defined database does not exist ${white}"
  cd ${Qmatey_dir}/tools/
  mkdir ncbi_db; mkdir ./ncbi_db/nt
  cd ${Qmatey_dir}/tools/ncbi_db/nt
  export local_db=${Qmatey_dir}/tools/ncbi_db/nt
  if [ "$(ls -A .)" ]; then
    echo -e "${magenta}- we recommend updating the NCBI nt database every few months. ${white}"
    echo -e "${magenta}- warning: updating might take a few hours and may significantly alter the consistency of your results ${white}"
    echo -e "${magenta}- would you like to update your NCBI nt database in Qmatey tools directory? ${white}"
    read -p "- y(YES) or n(NO) " -n 1 -r
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
      :
    else
      echo -e "${white}- updating NCBI nt database within Qmatey tools directory ${white}"
      update_database &>> ./log.out
    fi
  else
    echo -e "${white}- downloading and setting up NCBI nt database within Qmatey tools directory ${white}"
    setup_database &>> ./log.out
  fi
fi




if [ -d "./samples/pe" ]; then
  echo -e "${magenta}- Qmatey will modify your sample fastq files by concatenating <pe-R1-reads>, <pe-R2-reads>, and <se-reads> ${white}"
  echo -e "${magenta}- Do you want create a copy/backup of your sample fastq file? ${white}"
  read -p "- y(YES) or n(NO)? " -n 1 -r
  if [[ $REPLY = ^[Yy]$ ]]; then
    printf "\n"
    echo -e "${magenta}- A copy of your sample fastq files is in ${projdir}/samples_backup/ ${white}"
    mkdir samples_backup
    cp -r samples samples_backup
    else
      echo -e "${white}\n#########################################################################################"
  fi
else
  :
fi


####################################################################################################################
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

if [ -z "$string" ]; then
 string=true
fi
if [ "$cluster"  == false ]; then
 unset cluster
fi
string_norm="$(grep normalization config.sh)"
string_norm=${string_norm//*=}
if [ -z "$string_norm"  ]; then
	echo -e "${magenta}- Expecting variable <normalization: true or false> ${white}"
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





if [[ "$string" == false ]]; then
  echo -e "${magenta}- Qmatey will run in walk-through mode? ${white}"
  if [ -z "$nthread" ]; then
    echo -e "${magenta}- Qmatey will use (total processors/cores)-2? ${white}"
    echo -e "${magenta}- Quit (ctrl + c), else Qmatey will continue in 5 seconds ${white}"
    sleep 5
  fi
  cat header.txt config.sh fetchdir.txt $string2 $string3 | awk '{ sub("\r$",""); print}' > Qmatey_run.sh
  rm header.txt fetchdir.txt
  bash ${projdir}Qmatey_run.sh
fi





if [[ "$string" == true ]]; then
if [ -z $cluster ]; then
  echo -e "${magenta}- Qmatey will run in walkaway mode? ${white}"
  if [ -z "$nthread" ]; then
    echo -e "${magenta}- Qmatey will use (total processors/cores)-2? ${white}"
    echo -e "${magenta}- Quit (ctrl + c), else Qmatey will continue in 5 seconds ${white}"
    sleep 5
  fi
  printf "#""!""/bin/bash \n" > header.txt
  printf "walkaway=true \n\n" > walkaway.txt
  printf "Qmatey_dir=${Qmatey_dir}\n" > fetchdir.txt
  printf "projdir=${projdir}" >> fetchdir.txt
  cat header.txt config.sh walkaway.txt fetchdir.txt $string2 $string3 | awk '{ sub("\r$",""); print}' > Qmatey_run.sh
  rm header.txt walkaway.txt fetchdir.txt
  nohup bash ${projdir}Qmatey_run.sh > terminal.out 2>&1 &
fi
fi





if [[ "$string" == true ]]; then
if [[ "$cluster" == true ]]; then
  if [ -z "$thread_node" ]; then
    echo -e "${magenta}- Please provide number of threads config.sh (for cluster node)? ${white}"
    echo -e "${magenta}- Qmatey will quit in 5 seconds ${white}"
    sleep 5
    exit 1
  fi
  echo -e "${magenta}- Qmatey will run in walkaway mode? ${white}"
  printf "#""!""/bin/bash \n#SBATCH -c ${thread_node} \n\n" > cluster_header.sh
  printf "walkaway=true \n\n" > walkaway.txt
  cat cluster_header.sh steps.txt config.sh walkaway.txt fetchdir.txt $string2 $string3 | awk '{ sub("\r$",""); print}' > Qmatey_run.sh
  rm fetchdir.txt cluster_header.sh header.txt steps.txt walkaway.txt
  sbatch ${projdir}Qmatey_run.sh
fi
fi
}
cd $projdir
main
echo -e "${magenta}- For jobs running in background, monitor progress in terminal.out/slurm-xxxxxx.out ${white}"
echo -e "${magenta}- The log.out file can help with troubleshooting and when reporting a bug ${white}"
