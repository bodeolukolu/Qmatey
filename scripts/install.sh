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
  tar -vxjf download*; rm download*; cd samtools*; make; cd ..
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
  wget https://github.com/AdoptOpenJDK/openjdk8-binaries/releases/download/jdk8u222-b10/OpenJDK8U-jre_x64_linux_hotspot_8u222b10.tar.gz &&
  tar -xvf OpenJDK8U-jre_x64_linux_hotspot_8u222b10.tar.gz; rm *tar.gz
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
  echo -e "${white}\n############################################## ${orange}\n- installing NCBI Entrez Direct ${white}\n##############################################${white}"
  cd ~
  /bin/bash
  perl -MNet::FTP -e \
    '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
     $ftp->login; $ftp->binary;
     $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
  gunzip -c edirect.tar.gz | tar xf -
  rm edirect.tar.gz
  builtin exit
  export PATH=${PATH}:$HOME/edirect >& /dev/null || setenv PATH "${PATH}:$HOME/edirect"
  sleep 5; echo "y" | ./edirect/setup.sh 2> /dev/null
}
test_edirect=""
test_edirect=$(xtract -version)
if [[ -z "$test_edirect" ]]; then
  echo -e "${white}- Performing installation of NCBI Entrex Direct (BLAST) ${white}"
  main &>> ./log.out
fi

main () {
  if R --version; then
  	:
  else
  	echo -e "${white}- install R before proceeding ${white}"
    echo -e "${white}- dependencies for R in linux: <sudo apt install libcurl4-openssl-dev> and <sudo apt install libssl-dev>"
  	echo -e "${white}- Qmatey will exit in 10 secconds ${white}"
  	sleep 5s
  	exit 1
  fi
}
main &>> ./log.out

main () {
echo -e "${blue}\n############################################## \n- installing R-packages  ${blue}\n##############################################${white}"
  mkdir ./R
  cd ./R
  dirtool=plotly*; if [ ! -d $dirtool ]; then
		R -e 'install.packages("plotly", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
  wait
  dirtool=data.table; if [ ! -d $dirtool ]; then
	R -e 'install.packages("data.table", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
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
  dirtool=.htmlwidgets; if [ ! -d $dirtool ]; then
 	R -e 'install.packages("htmlwidgets", dependencies = TRUE, repos="http://cran.r-project.org", lib="./")'
  fi
}
main &>> ./log.out
