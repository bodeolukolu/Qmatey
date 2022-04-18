<p align="right">
<img src="https://github.com/bodeolukolu/Qmatey/blob/master/misc/Qmatey_logo.PNG" width="273" height="162">
</p>

# Introduction
Qmatey (0.40) is a quantitative metagenomic/microbiome profiling pipeline. Using the NCBI MegaBLAST, it implements a fast exact-matching algorithm for strain-level profiling. For species-level to phylum-level profiling, it implements exact-matching of consensus sequence that is unique to each taxa (e.g. at species-level, valid hits will match uniquely to each specific epitet).

## Features
* Exact-matching (and exact-matching of consensus) sequence.
* User-friendly and fully automated: “walk-away” and “walk-through” mode
* Input data: whole genome shotgun sequencing (WGS), reduced representation sequencing (RRS/qRRS), and 16S/ITS amplicon sequencing.
* User-defined parameters for strain- to phylum-level taxonomic identification and quantification.
* Data compression and indexing (reads of all samples into single file) improves speed (avoids alignment of the same read multiple times).
* Optimization of parallelized MegaBLAST jobs.
* Allows for various types of normalization.
* QC-plots to evaluate the predictive accuracy of profiles.
* visualization of results


Qmatey is undergoing beta testing. If you have any feedback pertaining to Qmatey's operation, please contact Bode Olukolu (bolukolu@utk.edu).

Developers:	Brandon Kristy (ORNL, TN), Alison K. Adams (UTK, TN), and Bode A. Olukolu (UTK, TN)



## Installation
Clone or download the git repository to a desired location

```
git clone https://github.com/bodeolukolu/Qmatey.git
```

## Dependencies
* NCBI MegaBLAST
* bwa
* samtools
* picard tools
* java
* R version 3.5 or above
* R dependencies (dependencies of packages not listed): plotme, ggplot, plotly, reshape2, ggcorrplot, ccrepe
* If you are performing a **local BLAST**, you will require an NCBI sequencing database compiled into one directory

## Setting Up a Project Directory
A project directory should contain the following sub-directories:
* Input Sequences
  * This is where your QC-filtered sequencing data will go.
* Configuration file
  * The format of the configuration file can be taken from the tools directory of the Qmatey Repository.
* **Optional**: Normalzation Reference Genomes **if you are normalizing data**
  * This is where reference genomes for normalization will go. Reference genomes must be in **Fasta format**.
## Preparing A Database Directory for a Local BLAST
If necessary, install the lftp tool to navigate NCBI's FTP site:
```
sudo apt-get install lftp
```
To obtain an NCBI sequencing database, go to the FTP site ftp://ftp.ncbi.nlm.nih.gov/blast/db/.
Create a database directory and select the sequencing database from the FTP site you wish to obtain using the following commands:
```
mkdir database_directory
cd database_directory
wget ftp ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*.tar.gz
```
In the above code, I am extracting the nt.00.tar.gz nucleotide database file into my database directory. A complete sequencing database will require an extensive amount of space.

Next, uncompress the database files and remove them with the following commands:
```
for f in nt*.tar.gz; do tar xzf "$f"; done
rm nt*.tar.gz
```

Your database directory should now have the desired, uncompressed database files.


# Limit search to taxa levels
use arguement: -taxids ### (e.g. -taxids 2 for bacteria)
Viruses		10239
Bacteria	2
Archaea		2157
Fungi		4751
Ecdysozoa	1206794
Animalia/Metazoan	33208
viridiplantae		33090

The taxids.list for the taxa levels listed above are alreadyt pre-packaged into Qmatey but if you want custom taxids.list, follow the instructions below from NCBI:
# Limiting a BLAST search with a high-level taxonomic node
$ get_species_taxids.sh -n Enterobacterales
Taxid: 91347
 rank: order
 division: enterobacteria
 scientific name: Enterobacterales
 common name:
1 matches found
$ get_species_taxids.sh -t 91347 > 91347.txids
$ blastn –db nt –query QUERY –taxidlist 91347.txids –outfmt 7 –out OUTPUT.tab
Go to:
# Limiting a BLAST search with a species-level taxonomic node
$ blastn –db nt –query QUERY –taxids 9606 –outfmt 7 –out OUTPUT.tab
# To use the get_spcies_taxids.sh, you will need to install Entrez Direct (EDirect). Copy and paste the commands below into your terminal
```
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
  ./edirect/setup.sh
```
# Additionally, one may use the -negative_taxids and -negative_taxidlist options to exclude sequences by taxid


## Configuration

Variable | Usage | Input
-------------- | ------------------------------------------------------------------- | -----
data_type      |  the type of data in your input_dir                                 | FASTQ or FASTA
input_dir      | the path to QC-filtered sequencing data                             | /path to directory/
norm_ref_dir        | the path to the normalization reference genomes                              | /path to directory/
threads        | the maximum number of subprocesses that can run simultaneously      | integer
tool_dir       | the path to Qmatey's tools                                          | /path to downloaded Qmatey repository/tools
strain_level   | An option for strain-level taxonomic analysis                       | TRUE or FALSE
species_level  | An option for species-level taxonomic analysis                      | TRUE or FALSE
genus_level    | An option for genus-level taxonomic analysis                        | TRUE or FALSE
family_level   | An option for family-level taxonomic analysis                       | TRUE or FALSE
order_level    | An option for order-level taxonomic analysis                        | TRUE or FALSE
class_level    | An option for class-level taxonomic analysis                        | TRUE or FALSE
phylum_level   | An option for phylum-level taxonomic analysis                       | True or FALSE
blast_location | An option to perform BLAST locally or remotely                      | LOCAL or REMOTE
local_db_dir   | the path to a local NCBI sequencing database on your desktop        | /path to directory/database name or NA
remote_db_dir  | the NCBI database for remote BLAST performance                    | e.g. nt, 16s, nr, etc. or NA
normalization  | An option for reference-based normalization                         | TRUE OR FALSE

## Usage
Before running Qmatey, make sure you have:
* Created a project directory with all the required subdirectories
* Have QC-filtered sequencing data in the appropriate input directory
* **for normalization**: Obtained necessary reference genome(s) in .fasta format in the appropriate normalization reference genome directory
* **for a local BLAST**: Obtained an NCBI database directory and have the uncompressed files in one directory
* **for a remote BLAST**: identified an NCBI database directory online
* Have a correctly edited configuration file within the project directory

From the command line, type:
```
$ bash <path to github repository>/Qmatey_v0.1.sh <path to project directory>/Qmatey.config
```
# License
<a href="https://github.com/tararickman/metagenome/blob/add-license-1/LICENSE"> GNU General Public License v3.0
