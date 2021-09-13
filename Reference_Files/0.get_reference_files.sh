#!/bin/bash

wget "http://huttenhower.sph.harvard.edu/humann_data/chocophlan/full_chocophlan.v296_201901b.tar.gz"
tar -xvzf full_chocophlan.v296_201901b.tar.gz
mkdir full_chocophlan.v296_201901b
mv g__* full_chocophlan.v296_201901b/
rm full_chocophlan.v296_201901b.tar.gz

wget "http://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz"
tar -xvzf uniref90_annotated_v201901b_full.tar.gz 
mkdir uniref90_annotated_v201901b_full
mv uniref90_201901b_full.dmnd uniref90_annotated_v201901b_full//
rm uniref90_annotated_v201901b_full.tar.gz 

wget "http://huttenhower.sph.harvard.edu/humann_data/full_mapping_v201901b.tar.gz"
tar -xvzf full_mapping_v201901b.tar.gz
wget "https://huttenhower.sph.harvard.edu/humann_data/map_uniprot_uniref.tsv.gz"
mkdir full_mapping_v201901b
mv map_* full_mapping_v201901b/
mv *tol-lca.dat.bz2 full_mapping_v201901b/
rm full_mapping_v201901b.tar.gz

wget "https://zenodo.org/record/3957592/files/mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2"

wget "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/reference/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_genomic.fna.gz"
gunzip GCA_000001405.28_GRCh38.p13_genomic.fna.gz

echo " *** Reference databases for running SLURM_Shotgun_Metagenomic_Pipeline.sh have been downloaded ***"

