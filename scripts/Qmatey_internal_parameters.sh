
BLUE='\e[38;5;210m'
WHITE='\e[97m'
ORANGE='\e[38;5;210m'

######################################################################################################################################################
export Qmatey_dir=${Qmatey_dir%/*}
export projdir=${projdir%/*}



######################################################################################################################################################
# Software defined parameters
cd ${Qmatey_dir}/tools
bash ../scripts/install.sh

######################################################################################################################################################
# tools
export bwa=${Qmatey_dir}/tools/bwa*/bwa && bwa=${bwa//'//'/'/'}
export samtools=${Qmatey_dir}/tools/samtools*/samtools && samtools=${samtools//'//'/'/'}
picard=${Qmatey_dir}/tools/picard.jar && picard=${picard//'//'/'/'}
java=${Qmatey_dir}/tools/jdk8*/bin/java && java=${java//'//'/'/'}
blast=${Qmatey_dir}/tools/ncbi-blast*/bin/blastn && blast=${blast//'//'/'/'}
