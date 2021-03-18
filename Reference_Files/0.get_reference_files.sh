#!/bin/bash

wget "http://huttenhower.sph.harvard.edu/humann2_data/chocophlan/full_chocophlan.v296_201901.tar.gz"
tar -xvzf full_chocophlan.v296_201901.tar.gz
rm full_chocophlan.v296_201901.tar.gz

wget "http://huttenhower.sph.harvard.edu/humann2_data/uniprot/uniref_annotated/uniref90_annotated_v201901.tar.gz"
tar -xvzf uniref90_annotated_v201901.tar.gz
rm uniref90_annotated_v201901.tar.gz

wget "http://huttenhower.sph.harvard.edu/humann2_data/full_mapping_v201901.tar.gz"
tar -xvzf full_mapping_v201901.tar.gz
rm full_mapping_v201901.tar.gz

wget "https://zenodo.org/record/3957592/files/mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2"

wget "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.bowtie_index.tar.gz"
tar -xvzf GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.bowtie_index.tar.gz
rm GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna.bowtie_index.tar.gz

wget "https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz"
gunzip SILVA_132_SSURef_Nr99_tax_silva.fasta.gz

echo "Reference databases and files for running SLURM_Shotgun_Metagenomic_Pipeline.sh have been downloaded, would recommend also running script 0.prep_SILVA_132_SSURef_Nr99_tax_silva.sh if wanting to do 16S based taxonomic classification option to prep SILVA reference FASTA for use with RDP classifier"

