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
echo -e "\e[97m########################################################\n \e[38;5;210m Generate taxids to limit BLAST output to specific taxonomic groups \n\e[97m########################################################\n"

taxids (Viruses): 10239
taxids (Bacteria):  2
taxids (Archaea): 2157
taxids (Protista/Protista-like): 554915,554296,2608109,2686027,590648,2608240,2611352,2489521,2795258,2611341,2598132,2698737,660925,2683617,2686024,127916,98350,2018064,2687318,28009
taxids (Algea/Algea-like): 3027,2763,38254,2806169,3041,96475,2218517,131220
taxids (Fungi): 4751
taxids (Oomycota):  4762
taxids (Nematoda): 6231
taxids (Arthropoda): 6656

# Specify list of taxonomic groups below
taxids=10239,2,2157,4751,4762

# mkdir -p taxids
# cd taxids
taxids=$( echo $taxids | awk '{gsub(/,/," ")}1' )
for get_species_taxids in ${taxids[@]}; do
	if [[ -z ${projdir}/taxids/${get_species_taxids}.txids ]]; then
		:
	else
		~/Documents/tools/Qmatey/tools/ncbi-blast-2.14.1+/bin/get_species_taxids.sh -t ${get_species_taxids} > ${get_species_taxids}.txids
	fi
done



#################################################################################################################
echo -e "\e[97m########################################################\n \e[38;5;210m Generate taxids to remove from MegaBLAST output \n\e[97m########################################################\n"

# Specify list of taxonomic groups below
taxids=4640

# mkdir -p taxids
# cd taxids
taxids=$( echo $taxids | awk '{gsub(/,/," ")}1' )
for get_species_taxids in ${taxids[@]}; do
	if [[ -z ${projdir}/taxids/${get_species_taxids}.txids ]]; then
		:
	else
		~/Documents/tools/Qmatey/tools/ncbi-blast-2.12.0+/bin/get_species_taxids.sh -t ${get_species_taxids} > remove_taxa.txids
	fi
done
