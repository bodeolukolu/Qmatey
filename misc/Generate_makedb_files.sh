#!/bin/bash
threads=12

# Download 16S (Silva) and ITS (UNITE) fasta files
#
# https://zenodo.org/record/4587955/files/silva_species_assignment_v138.1.fa.gz?download=1 -o silva_species_assignment_v138.1.fa.gz

cat silva_species_assignment_v138.1.fa.gz unite.fa.gz > 16S_ITS.fa.gz

gunzip 16S_ITS.fa.gz
mv 16S_ITS.fa 16S_ITS.fasta

grep '>' 16S_ITS.fasta | awk '{gsub(/>/,""); print $1}' > taxaid_list.txt
while read line; do
acc=${line%%.*}
(esummary -db nuccore -id $acc | xtract -pattern DocumentSummary -element TaxId | awk -v line=$line '{print line"\t"$1}' >> map_taxaids.txt
if [[ $(jobs -r -p | wc -l) -ge $threads ]]; then
wait
fi
done < taxaid_list.txt
