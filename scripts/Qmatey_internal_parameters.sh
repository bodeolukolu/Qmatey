
BLUE='\e[38;5;210m'
WHITE='\e[97m'
ORANGE='\e[38;5;210m'

######################################################################################################################################################
export Qmatey_dir=${Qmatey_dir%/*}
export projdir=${projdir%/*}

pop=${projdir%/}
export pop=${pop##*/}

######################################################################################################################################################
# tools
export bwa=${Qmatey_dir}/tools/bwa*/bwa && bwa=${bwa//'//'/'/'}
export samtools=${Qmatey_dir}/tools/samtools*/samtools && samtools=${samtools//'//'/'/'}
export picard=${Qmatey_dir}/tools/picard.jar && picard=${picard//'//'/'/'}
export java=${Qmatey_dir}/tools/jdk8*/bin/java && java=${java//'//'/'/'}
export cd_hit_est=${Qmatey_dir}/tools/cd-hit*/cd-hit-est && cdhit=${cdhit//'//'/'/'}
export blast=${Qmatey_dir}/tools/ncbi-blast*/bin/blastn && blast=${blast//'//'/'/'}
if command -v pigz &>/dev/null; then
  export gzip=pigz
else
  export gzip=gzip
fi

slurm_module=$(module --version 2> /dev/null | head -n1)

if [[ "$slurm_module" =~ "Module"  ]]; then
  samtoolsout=$($samtools --version | head -n 3)
  if [ -z "$samtoolsout" ];then
    echo -e "${white}- samtools installation within Qmatey is probably missing a dependency on host system ${white}"
    echo -e "${white}- Qmatey will use host system samtools installation ${white}"
    module add samtools
    export samtools=samtools
    $samtools --version | head -n 3
  else
    $samtools --version | head -n 3
  fi

  Rout=$(R --version 2>/dev/null | head -n 3)
  if [ -z "$Rout" ];then
  	echo -e "${white}- R not found, so loading R with module load/add on cluster  ${white}"
    module add R
  	R --version
  else
    R --version
  fi
fi

cd $Qmatey_dir/tools/R
Rscript ../../scripts/check_R_tools.R "${Qmatey_dir}/tools/R"
cd $projdir
