#!/bin/bash
magenta=$'\e[1;35m'
blue='\e[38;5;210m'
white='\e[97m'
orange='\e[38;5;210m'


######################################################################################################################################################
main () {
  echo -e "${white}\n############################################## ${orange}\n- downloading and installing BWA ${white}\n##############################################${white}"
  wget https://sourceforge.net/projects/bio-bwa/files/latest/download &&
  tar -vxjf download*; rm download*; cd bwa*; make; cd ..
}
dirtool=bwa*
if [ -d $dirtool ]; then
  :
else
  echo -e "${white}- Performing installation of dependency (bwa) ${white}"
  main &>> ./log.out
fi



main () {
  echo -e "${white}\n############################################## ${orange}\n- downloading and installing samtools ${white}\n##############################################${white}"
  wget https://sourceforge.net/projects/samtools/files/latest/download &&
  tar -xvjf download*; rm download*; cd samtools*; make; cd ..
}
dirtool=samtools*
if [ -d $dirtool ]; then
  :
else
  echo -e "${white}- Performing installation of dependency (samtools) ${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${white}\n############################################## ${orange}\n- downloading PICARD tools ${white}\n##############################################${white}"
  wget https://github.com/broadinstitute/picard/releases/download/2.25.6/picard.jar
}
dirtool=picard*
if [ -f $dirtool ]; then
  :
else
  echo -e "${white}- Performing installation of dependency (picard tools) ${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${white}\n############################################## ${orange}\n- installing java 1.8. It will be run in manual mode ${white}\n##############################################${white}"
  wget https://github.com/adoptium/temurin8-binaries/releases/download/jdk8u322-b06/OpenJDK8U-jdk_x64_linux_hotspot_8u322b06.tar.gz
  tar -xvf OpenJDK8U-jdk_x64_linux_hotspot_8u322b06.tar.gz; rm *tar.gz
}
dirtool=jdk*
if [ -d $dirtool ]; then
  :
else
  echo -e "${white}- Performing installation of dependency (java 1.8)${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${white}\n############################################## ${orange}\n- installing latest version of BLAST ${white}\n##############################################${white}"
  wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.12.0+-x64-linux.tar.gz &&
  tar -zxvf ncbi-blast-2.12.0+-x64-linux.tar.gz; rm *tar.gz
}
dirtool=ncbi-blast*
if [ -d $dirtool ]; then
  :
else
  echo -e "${white}- Performing installation of dependency (BLAST) ${white}"
  main &>> ./log.out
fi


main () {
  echo -e "${white}\n############################################## ${orange}\n- check for R installation ${white}\n##############################################${white}"
  if R --version; then
    :
  else
    module add R
    if R --version; then
      :
    else
      echo -e "${white}- install R before proceeding ${white}"
      echo -e "${white}- dependencies for R in linux: <sudo apt install libcurl4-openssl-dev> and <sudo apt install libssl-dev>"
    fi
  fi
}
main &>> ./log.out


main () {
echo -e "${blue}\n############################################## \n- installing R-packages  ${blue}\n##############################################${white}"
  mkdir -p R
  cd R
  dirtool=plotly*; if [ ! -d $dirtool ]; then
		R -e 'install.packages("plotly", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
  wait
  dirtool=htmlwidgets; if [ ! -d $dirtool ]; then
    R -e 'install.packages("htmlwidgets", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
  wait
  dirtool=data.table; if [ ! -d $dirtool ]; then
	   R -e 'install.packages("data.table", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
  wait
  dirtool=rehsape2; if [ ! -d $dirtool ]; then
     R -e 'install.packages("reshape2", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
  wait
  dirtool=ggplot2; if [ ! -d $dirtool ]; then
		R -e 'install.packages("ggplot2", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
  wait
  dirtool=ggcorrplot; if [ ! -d $dirtool ]; then
	   R -e 'install.packages("ggcorrplot", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
  wait
  dirtool=stringi; if [ ! -d $dirtool ]; then
	   R -e 'install.packages("stringi", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
  wait
  dirtool=infotheo; if [ ! -d $dirtool ]; then
     R -e 'install.packages("infotheo", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
  wait
  dirtool=ccrepe; if [ ! -d $dirtool ]; then
    wget https://www.bioconductor.org/packages/release/bioc/src/contrib/ccrepe_1.30.0.tar.gz
    gunzip ccrepe_1.30.0.tar.gz
    R -e 'install.packages("ccrepe_1.30.0.tar", dependencies = TRUE, repos = NULL, type="source", lib="./")'
    rm ccrepe_1.30.0.tar
  fi
  dirtool=plotme; if [ ! -d $dirtool ]; then
    wget https://github.com/yogevherz/plotme/archive/refs/heads/master.tar.gz -O plotme.tar.gz
    gunzip plotme.tar.gz
    R -e 'install.packages("plotme.tar", dependencies = TRUE, repos = NULL, type="source", lib="./")'
    rm plotme.tar
  fi
  cd ../
}
echo -e "${white}- Checking and performing installation of R-packages ${white}"
main &>> ./log.out
