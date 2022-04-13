
# General parameters
####################################################
threads=24
walkaway=true
cluster=false
samples_alt_dir=false


# Normalization
####################################################
normalization=true

# MegaBLAST
####################################################
blast_location=local
local_db=/media/sdd/ncbi_db/nt/nt
# local_db=/media/sdb/ncbi_db/16S/16S_ribosomal_RNA,/media/sdb/ncbi_db/18S/18S_fungal_sequences,/media/sdb/ncbi_db/28S/28S_fungal_sequences,/media/sdb/ncbi_db/ITS/ITS_eukaryote_sequences,
#local_db=/media/sdd/ncbi_db/refseq/refseq_rna,/media/sdd/ncbi_db/refseq/ref_viroids_rep_genomes,/media/sdd/ncbi_db/refseq/ref_prok_rep_genomes,/media/sdd/ncbi_db/refseq/ref_euk_rep_genomes
remote_db=NA
input_dbfasta=NA
taxids=true
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

run_corr=true
sunburst_taxlevel=strain
sunburst_nlayers=phylum,genus,species
####################################################

# Filtering Parameters
####################################################
min_unique_seqs=1
min_percent_sample=5,10,20
min_pos_corr=0.1,0.2,0.3
max_neg_corr=0.1,0.2,0.3
genome_scaling=0
####################################################
run_sunburst=true
