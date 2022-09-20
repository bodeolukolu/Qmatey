
# General_parameters
####################################################
threads=24
walkaway=true
cluster=false
samples_alt_dir=false
library_type=qRRS


# Normalization
####################################################
normalization=true

# MegaBLAST
####################################################
blast_location=local
local_db=/media/sdd/ncbi_db/nt/nt
# local_db=/media/sdd/ncbi_db/nt/nr
# local_db=/media/sdb/ncbi_db/16S/16S_ribosomal_RNA,/media/sdb/ncbi_db/18S/18S_fungal_sequences,/media/sdb/ncbi_db/28S/28S_fungal_sequences,/media/sdb/ncbi_db/ITS/ITS_eukaryote_sequences
# local_db=/media/sdd/ncbi_db/refseq/refseq_rna,/media/sdd/ncbi_db/refseq/ref_viroids_rep_genomes,/media/sdd/ncbi_db/refseq/ref_prok_rep_genomes,/media/sdd/ncbi_db/refseq/ref_euk_rep_genomes
taxids=true
input_dbfasta=NA
map_taxids=NA


# Taxonomic_Profiling_and_Filtering
####################################################
minRD=2
min_percent_sample=10,20
min_pos_corr=0.1,0.2,0.3
max_neg_corr=0.1,0.2,0.3


# Visualizations
####################################################
sunburst_taxlevel=strain,species,genus,family,order,class
sunburst_nlayers=phylum,genus,species
compositional_corr=strain,species,genus,family,order,class,phylum


# Advanced_Parameters
####################################################
nodes=1
reads_per_megablast=1000
genome_scaling=true
zero_inflated=0.01
qcov=50
exclude_rRNA=true
