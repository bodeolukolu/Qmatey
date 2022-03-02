#!/bin/bash




# Option-1
sudo apt install ncbi-entrez-direct

# Option-2
# cd ~
# /bin/bash
# perl -MNet::FTP -e \
#   '$ftp = new Net::FTP("ftp.ncbi.nlm.nih.gov", Passive => 1);
#    $ftp->login; $ftp->binary;
#    $ftp->get("/entrez/entrezdirect/edirect.tar.gz");'
# gunzip -c edirect.tar.gz | tar xf -
# rm edirect.tar.gz
# builtin exit
# export PATH=$PATH:$HOME/edirect >& /dev/null || setenv PATH "${PATH}:$HOME/edirect"


#################################################################################################################
echo -e "\e[97m########################################################\n \e[38;5;210m Generate taxids to subset BLAST output \n\e[97m########################################################\n"

taxids=10239,2,2157,4751,4762,6231,6843,61985,554674,30264,33342,33341,33339,30261,7509,85817,27420,7399,7041,85604,7203,7148

mkdir -p taxids
cd taxids
taxids=$( echo $taxids | awk '{gsub(/,/," ")}1' )
for get_species_taxids in ${taxids[@]}; do
	if [[ -z ${projdir}/taxids/${get_species_taxids}.txids ]]; then
		:
	else
		~/Documents/tools/Qmatey/tools/ncbi-blast-2.12.0+/bin/get_species_taxids.sh -t ${get_species_taxids} > ${get_species_taxids}.txids
	fi
done
