<p align="right">
<img src="https://github.com/bodeolukolu/Qmatey/blob/master/misc/Qmatey_logo.PNG" width="273" height="162">
</p>


# Introduction
Qmatey (v0.5.1) is a quantitative metagenomic/microbiome taxonomic binning and profiling pipeline. Using the NCBI MegaBLAST, it implements a fast exact matching algorithm for strain-level profiling. For species to phylum level profiling, it implements exact-matching of consensus sequence (EMC) that is unique to each taxa (e.g. at species-level, valid hits will match uniquely to each species name; at genus-level, valid hits will match uniquely to each genus name). Qmatey can also perform simulation of mock/synthetic communities using their whole genome assemblies.


## Features
* User-friendly and fully automated
* Exact-matching at strain Level
* Exact-matching of consensus sequence (EMC) at species- to phylum-level.
* User-defined parameters for strain to phylum level taxonomic classification.
- Multiple nucleotide databases can be specified for a run.
* Input data: amplicon sequencing (e.g. 16S/ITS), shotgun/WGS, reduced representation sequencing (RRS/qRRS).
* Simulates metagenome sequencing/profiling (complete/partial digest or shotgun/WGS libraries).
* Data compression and indexing improves speed by hundreds to thousands of times.
* Speed optimization of parallelized MegaBLAST jobs.
* Relative abundance based on spike-in/host sequences (recommended).
- Minimize false positive and false negative rates.
* implements a novel "cross-taxon validation" method.
- correlation network based on CCLasso (sensitive to compositional data/multi-way interactions) and Spearman corrlation method.
* QC-plots to evaluate the predictive accuracy of profiles.
* Various visualizations.

If you have any feedback pertaining to Qmatey's operation, please contact Bode Olukolu (bolukolu@utk.edu).

## Developers
* Alison K. Adams (UTK, TN)
* Brandon Kristy (ORNL, TN)
* Bode A. Olukolu (UTK, TN)


# Contents
- [Installation](#installation)
- [Usage](#usage)
  - [Basic usage](#Basic_usage)
  - [Project directory setup](#Project_directory_setup)
  - [Prepare database for Local MegaBLAST](#Prepare_database_for_Local_MegaBLAST)
  - [Limit search to taxon](#Limit_search_to_taxon)
  - [Create custom database](#Create_custom_database)
  - [Overview of workflow](#Overview_of_workflow)
  - [Configuration](#Configuration)
- [Related Software](#Related_Software)
- [Select Article Referencing Qmatey](#Select_Article_Referencing_Qmatey)
- [Acknowledgment](#Acknowledgment)
- [Troubleshooting](#Troubleshooting)
- [Versioning](#Versioning)
- [License](#License)


## Installation
- Currently, Qmatey is only available for unix-based systems.
- Clone or download the Git repository to your desired folder.
```bash
git clone https://github.com/bodeolukolu/Qmatey.git
```
- Installation occurs automatically the first time you run the pipeline.
- To install dependencies without running a job:  
```bash
Qmatey_dir/Qmatey install
```
- With the exception of the R sotfware, all dependencies are installed to a local directory within Qmatey.



## Dependencies
* NCBI MegaBLAST
* bwa
* samtools
* picard tools
* java
* R version 3.5 or above
* R dependencies: plotme, ggplot, plotly, reshape2, and ggcorrplot.
* If you are performing a **local BLAST**, you will require an NCBI sequence database downloaded to a directory


## Usage
### Basic Usage
Before running Qmatey, make sure you have created a project directory with the required files (config.sh, *.taxids, fasta/fastq files in samples folder, and reference fasta files in norm_ref folder)
- **quality filtering**: Fastq data should be QC-filtered with high stringency to prevent false positive and false negatives (paired-end and/or single-end reads allowed)
- **normalization and non-target sequences**: provide necessary reference genome(s) or spike-in standard in FASTA format (within norm_ref folder). These allows for normalization and removal of host sequence (particularly in endophytic metagenome/microbiome communities).
- **local MegaBLAST**: obtain and decompress NCBI database or databases. Custom databases can be created within Qmatey by providing fasta and tax id files (see below).
- **remote MegaBLAST**: specify remote in config.sh file (not recommended due to slow speed).
- **edit parameters in config.sh**: template available in examples folder within Qmatey download.
- **for help**: --help or -h
- **for version**: --version or -v

From the command line, type:
```
$ bash <path-to-Qmatey-directiry/Qmatey> <path-to-project-directory>
```

### Project directory setup
A project directory should contain the following sub-directories:
- **Input Sequences**: This is where your QC-filtered sequencing data will go.
- **Configuration file**: The format of the configuration file can be taken from the tools directory of the Qmatey Repository.
- **Optional**: Normalzation Reference Genomes **if you are normalizing data**
  * This is where reference genomes for normalization will go. Reference genomes must be in **Fasta format**.

  <img src="https://github.com/bodeolukolu/Qmatey/blob/master/misc/project_dir_setup.PNG" width="739" height="426">


### Prepare database for Local MegaBLAST
If necessary, install the lftp tool. For more details about tool: https://lftp.yar.ru/lftp-man.html

```
sudo apt update
sudo apt-get install lftp
```
To obtain an NCBI sequencing database, go to the FTP site ftp://ftp.ncbi.nlm.nih.gov/blast/db/.
Create a database directory and select the sequencing database from the FTP site you wish to obtain using the following commands:

```
mkdir database_directory
cd database_directory
wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*.tar.gz
```
Next, uncompress the database files and delete the original .tar.gz files:
```
for f in nt*.tar.gz; do tar xzf "$f"; done
rm nt*.tar.gz
```

Your database directory should now have the desired, uncompressed database files, and ready to use.
To perform megablast with multiple databases, provide uncompressed database files in the same directory.


# Limit search to taxon
To limit search of taxon in database, list of all taxids of interest within specific taxonomic groups (*.taxids files in examples folder) can be specified in MegaBLAST search.
-**basis for excludinng taxa**: large multicellular organisms (plants and animals) are often host rather than players within metagenomic communities.
-**mitigating false negatives due to endophytic microbial sequences assembled as host sequence**: organisms, particularly multicellular ones, often present a problem for exact-matchign for metagenomic studies. These genome sequences/assemblies often contain contaminating sequences from endophytic/epiphytic organisms. To avoid false negatives, exact matching is only performed among individuals within higher order taxonomic group (e.g. same phylum).

```
-**pre-selected taxids (Viruses)**: 10239
-**pre-selected taxids (Bacteria)**:  2
-**pre-selected taxids (Archaea)**: 2157
-**pre-selected taxids (Protista/Protista-like)**: 554915,554296,2608109,2686027,590648,2608240,2611352,2489521,2795258,2611341,2598132,2698737,660925,2683617,2686024,127916,98350,2018064,2687318,28009
-**pre-selected taxids (Algea/Algea-like)**:3027,2763,38254,2806169,3041,96475,2218517,131220
-**pre-selected taxids (Fungi)**: 4751
-**pre-selected taxids (Oomycota)**:  4762
-**pre-selected taxids (Nematoda)**: 6231
-**pre-selected taxids (Arthropoda)**: 6656
```

Exlcuded taxon
```
Opisthokonta (except unicellular organisms, Fungi, Oomycota, Nematoda, Arthropda)
Viridiplantae (except unicellular organisms, Algea)
```
To make generate a custom list of taxids.list, use the NCBI entrez-direct and BLAST tool get_species_taxids.sh as shown in script included with the Qmatey download (misc folder).


### Create custom database
Provide fasta and map_taxids files. Qmatey will use NCBI BLAST tool, makeblastdb, and these files to create the custom database. Example files are provided with the Qmatey download (examples) folder.


### Overview of workflow
<img src="https://github.com/bodeolukolu/Qmatey/blob/master/misc/Qmatey_workflow.PNG" width="876" height="301">

### configuration
Using a text editor, save a file containing any of the following variables as 'config.sh' file and include it in your project directory.

**General parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|threads|na|number of cores/processors|integer|Optional|
|cluster|false|run on compute cluster node (default: slurm) or workstation|true or false|Optional|
|samples_alt_dir|false|links samples in separate directory to project directory|true or false|Optional|
|lib_type|RRS|original input data: RRS (reduced representation sequence e.g. GBS), WGS (shotgun whole genome sequence), or 16S/ITS/amplicon|string|required|
|HDsubsample|false|high-density subsampling uses reads with all combination of motifs flanking reads|string|Optional|
|subsample_shotgun_R1|true|subsample reads based on? REnase motifs, true(ATGCAT) or false. specify different REnase instead of true|string|Optional|
|subsample_shotgun_R2|true|subsample reads based on? REnase motifs, true(CATG) or false. specify different REnase instead of true|string|Optional|
|shotgun_min_read_length|100|minimum read length after subsampling|string|Optional|

**Simulation parameters**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|simulation_lib|complete_digest|generate sequence reads in silico (complete_digest, partial_digest, or shotgun)|string|Optional|
|simulation_motif_R1|ATGCAT|genome fragmentation based on: (i) REnase site motif(s) or (ii) random |string|Optional|
|simulation_motif_R2|CATG|genome fragmentation based on: (i) REnase site motif(s) or (ii) random |string|Optional|
|fragment_size_range|64,600|minimum and maximum genomic fragment size (comma-separated)|string|Optional|
|max_read_length|150|maximum read length for each of R1 and/or R2 reads |string|Optional|
|gcov|3|whole genome sequencing coverage |string|Optional|


**Normalzation**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|normalization|na|to compute relative abundance using various methods|true or false|Optional|


**MegaBLAST**

|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|blast_location|na|local,remote, or custom|string|required|
|local_db|na|directory: ./ncbi_db/nt/nt |comma delimited string(s)|required|
|taxids|na|files mapping fasta sequences to taxid(s)|true or false|Optional|
|input_dbfasta|na|provide fasta file if custom database|string|Optional|
|map_taxids|na|provide files mapping fasta sequences to taxid(s) if custom database|string|Optional|


**Taxonomic_Profiling_and_Filtering**
|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|taxonomic_level|na|strain,species,genus,family,order,class,phylum|string|Optional|
|min_percent_sample|5,10,20|percentage of missing hits per sample allowed|comma delimited integer(s)|Optional|
|min_pos_corr|0.1,0.2,0.3|correlation coefficient threshold(s)|comma delimited decimal number(s)|Optional|
|max_neg_corr|0.1,0.2,0.3|correlation coefficient threshold(s)|comma delimited decimal number(s)|Optional|
|spearman_corr|na|Options: strain,species,genus,family,order,class,phylum|comma delimited string(s). Faster than CCLasso|Optional|
|CCLasso_corr|na|Options: strain,species,genus,family,order,class,phylum|comma delimited string(s). Analysis might take hours to days|Optional|


**Visualizations**
|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|sunburst_taxlevel|na|strain,species,genus,family,order,class|comma delimited string(s)|Optional|
|sunburst_nlayers|na|phylum,genus,species shown in sunburst|comma delimited string(s)|Optional|




**Advanced parameters**
|Variable      |Default       |Usage         |Input         |required/Optional|
|:-------------|:-------------|:-------------|:-------------|:----------------|
|nodes|1|number of nodes|integer|Optional|
|minRD|0|threshold excludes reads with erroneous base call(s); default (0=auto) is determined based on read depth distribution|integer|Optional|
|fullqlen_alignment|false|full query length should be consumed/aligned, but allowing for upto 5% mismatch|true or false|Optional|
|reads_per_megablast|na|batch size of reads/thread. Default: 1000/10000 (reads_per_megablast_burn_in=0/10000) for RRS/WGS and 100 for 16S/ITS/amplicon data|integer|Optional|
|reads_per_megablast_burn_in|10000|reduces database size by performing burn-in MegaBLAST run. Default=0 for custom databases (i.e. no burn-in) |integer|Optional|
|zero_inflated|0.01|exclude samples with proportion of zero taxa <= value|integer|Optional|
|exclude_rRNA|true|exclude rRNA for qRRS/WGS data|true or false|Optional|
|annotate_seq|false|generates gene ids and their abundance|true or false|Optional|

**Note: na indicates that variable is user-defined or intuitively hard-coded/computed.*

Below is an example of a configuration file:

**config.sh**
```
### General_parameters
####################################################
threads=24
cluster=false
samples_alt_dir=false
library_type=qRRS
HDsubsample=false
subsample_shotgun_R1=true
subsample_shotgun_R2=true
shotgun_min_read_length=50

# simulation_parameters
####################################################
simulation_lib=complete_digest
fragment_size_range=50,600
max_read_length=150
gcov=3


### Normalization
####################################################
normalization=true

### MegaBLAST
####################################################
blast_location=local
local_db=/media/sdd/ncbi_db/nt/nt
local_db=/media/sdd/ncbi_db/nt/nr
local_db=/media/sdb/ncbi_db/16S/16S_ribosomal_RNA,/media/sdb/ncbi_db/18S/18S_fungal_sequences,/media/sdb/ncbi_db/28S/28S_fungal_sequences,/media/sdb/ncbi_db/ITS/ITS_eukaryote_sequences
local_db=/media/sdd/ncbi_db/refseq/refseq_rna,/media/sdd/ncbi_db/refseq/ref_viroids_rep_genomes,/media/sdd/ncbi_db/refseq/ref_prok_rep_genomes,/media/sdd/ncbi_db/refseq/ref_euk_rep_genomes
taxids=true
input_dbfasta=NA
map_taxids=NA


### Taxonomic_Profiling_and_Filtering
####################################################
taxonomic_level=strain,species,genus,family,order,class,phylum
spearman_corr=strain,species,genus,family,order,class,phylum
CCLasso_corr=strain,species,genus,family,order,class,phylum
min_percent_sample=10,20
min_pos_corr=0.1,0.2,0.3
max_neg_corr=0.1,0.2,0.3


### Visualizations
####################################################
sunburst_taxlevel=strain,species,genus,family,order,class
sunburst_nlayers=phylum,genus,species



### Advanced_Parameters
####################################################
nodes=1
minRD=0
fullqlen_alignment=false
reads_per_megablast=10000
reads_per_megablast_burn_in=10000
zero_inflated=0.01
exclude_rRNA=false
annotate_seq=false
```


## Related Software
- [ngsComposer: Empirical Base-call error-filtering and read preprocessing pipeline.](https://github.com/bodeolukolu/ngsComposer)
- [GBSapp: automated pipeline for variant calling and filtering.](https://github.com/bodeolukolu/GBSapp)



## Select Article Referencing Qmatey
1.	Adams AK, Kristy BD, Gorman M, Balint-Kurti P, Yencho GC, and Olukolu BA. (2023) Qmatey: An automated pipeline for fast exact matching-based alignment and strain-level taxonomic binning and profiling of metagenomes. Briefings in Bioinformatics. 24 (6): bbad351.

2. Adams A, Landry D., Sykes V, Rickman T, Cham A, Timling A, Kelly H, McBeath JH, and Olukolu BA (2024) Maize Kernel-Associated Metagenomes Reveal Potential Microbe-Microbe Interactions that Underlie Fusarium Ear Rot Disease. Phytobiomes. doi.org/10.1094/PBIOMES-07-23-0074-R


## Troubleshooting
**Pre-Installation of R:**<br />
```
    - To view R version (or to check if installed), from the terminal type:
          $ python --version

    - For Ubuntu, install R using apt:
          $ sudo apt update && sudo apt upgrade
          $ sudo apt install r-base

    - For macOS, install using homebrew:
          brew install r
**If samtools doesn't install properly:**<br />
```
While the installation of samtools are automated, the installation requires some dependencies. Consider typing the commands below in terminal:
  $ sudo apt-get update
  $ sudo apt-get install gcc
  $ sudo apt-get install make
  $ sudo apt-get install libbz2-dev
  $ sudo apt-get install zlib1g-dev
  $ sudo apt-get install libncurses5-dev
  $ sudo apt-get install libncursesw5-dev
  $ sudo apt-get install liblzma-dev
  $ sudo apt-get install libcurl4-gnutls-dev
  $ sudo apt-get install libssl-dev
```
**Problem with amount of memory and/or processors/cores specified:**<br />
```
- This might be due to specifing values greater than available resources
- Re-submit job with appropriate values or modify header of the Qmatey_run.sh batch file.
- If using compute cluster managers other than SLURM, header of Qmatey_run.sh batch can also be modified to fit the syntax of the cluster manager been used.
```
**Problem with amount of memory and/or processors/cores specified:**<br />
```
- If BLAST search fails after updating NCBI database, consider updating the list of taxids in the taxids directory.
- You can use script in the Generate_taxids.sh file to update the taxid list.
```
## Versioning
Versioning will follow major.minor.patch <a href="https://semver.org">semantic versioning format</a>.

## License
<a href="https://github.com/bodeolukolu/Qmatey/blob/master/LICENSE">Apache License Version 2.0</a>
