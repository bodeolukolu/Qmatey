
# General_parameters
####################################################
threads=24
cluster=false
samples_alt_dir=false
library_type=qRRS
HDsubsample=false
subsample_shotgun_R1=true
subsample_shotgun_R2=true
shotgun_min_read_length=100

# simulation_parameters
####################################################
simulation_lib=complete_digest
simulation_motif_R1=ATGCAT
simulation_motif_R2=CATG
fragment_size_range=64,600
max_read_length=150
gcov=3

# Normalization
####################################################
normalization=true

# MegaBLAST
####################################################
blast_location=local
local_db=/media/sdd/ncbi_db/nt/nt
# local_db=/media/sdb/ncbi_db/16S/16S_ribosomal_RNA,/media/sdb/ncbi_db/18S/18S_fungal_sequences,/media/sdb/ncbi_db/28S/28S_fungal_sequences,/media/sdb/ncbi_db/ITS/ITS_eukaryote_sequences
# local_db=/media/sdd/ncbi_db/refseq/refseq_rna,/media/sdd/ncbi_db/refseq/ref_viroids_rep_genomes,/media/sdd/ncbi_db/refseq/ref_prok_rep_genomes,/media/sdd/ncbi_db/refseq/ref_euk_rep_genomes
taxids=true
input_dbfasta=NA
map_taxids=NA

# Taxonomic_Profiling_and_Filtering
####################################################
taxonomic_level=strain,species,genus,family,order,class,phylum
spearman_corr=strain,species,genus,family,order,class,phylum
CCLasso_corr=strain,species,genus,family,order,class,phylum
min_percent_sample=10,20
min_pos_corr=0.1,0.2,0.3
max_neg_corr=0.1,0.2,0.3

# Visualizations
####################################################
sunburst_taxlevel=strain,species,genus,family,order,class
sunburst_nlayers=phylum,genus,species


# Advanced_Parameters
####################################################
nodes=1
minRD=0
fullqlen_alignment=false
reads_per_megablast=10000
reads_per_megablast_burn_in=10000
zero_inflated=0.01
exclude_rRNA=false
annotate_seq=false
