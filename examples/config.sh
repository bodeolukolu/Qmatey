
# General parameters
####################################################
threads=24
walkaway=true
cluster=false

# Normalization
####################################################
normalization=true

# MegaBLAST
####################################################
blast_location=local
local_db=/media/sdd/ncbi_db/nt
remote_db=NA
input_dbfasta=NA
taxids=10239,2,2157,4751,6231,6843,61985,554674,30264,33342,33341,33339,30261,7509,85817,27420,7399,7041,85604,7203,7148
map_taxids=NA

# Taxonomic Filtering
####################################################
strain_level=true
species_level=true
genus_level=true
family_level=true
order_level=true
class_level=true
phylum_level=true
####################################################

# Filtering Parameters
####################################################
min_unique_seqs=1
min_percent_sample=5,10,20
min_pos_corr=0.1,0.2,0.3
max_neg_corr=0.1,0.2,0.3
genome_scaling=0
cross_ref=false
cross_ref_db_dir=false
####################################################
