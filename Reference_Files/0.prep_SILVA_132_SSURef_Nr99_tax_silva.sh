#!/bin/bash
#Prep SILVA_132_SSURef_Nr99_tax_silva.fasta for use with shotgun metagenomic pipeline

#Download file
wget "https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_Nr99_tax_silva.fasta.gz"
gunzip SILVA_132_SSURef_Nr99_tax_silva.fasta.gz

#Substitute Uracils for Thymines
sed -i '/^>/!s/U/T/g' SILVA_132_SSURef_Nr99_tax_silva.fasta

#Remove unneccessary newlines in sequences
awk '!/^>/ { printf "%s", $0; n = "\n" }
     /^>/ { print n $0; n = "" }
     END { printf "%s", n }' SILVA_132_SSURef_Nr99_tax_silva.fasta > temp.fasta
mv temp.fasta SILVA_132_SSURef_Nr99_tax_silva.fasta

#Replace spaces within the taxonomy lineages with an underscore
sed -i -e '/^>/s/ /_/g' -e '/^>/s/_/ /1' SILVA_132_SSURef_Nr99_tax_silva.fasta

#Replace backslashes with underscore
sed -i 's#/#_#g' SILVA_132_SSURef_Nr99_tax_silva.fasta

#Replace asterisks with underscore
sed -i 's/\*/_/g' SILVA_132_SSURef_Nr99_tax_silva.fasta

#Remove any brackets, parentheses, and greater than/less than signs
sed -i 's/\[//g' SILVA_132_SSURef_Nr99_tax_silva.fasta
sed -i 's/]//g' SILVA_132_SSURef_Nr99_tax_silva.fasta
sed -i 's/(//g' SILVA_132_SSURef_Nr99_tax_silva.fasta
sed -i 's/)//g' SILVA_132_SSURef_Nr99_tax_silva.fasta

#Remove any single quotations
sed -i "s/'//g" SILVA_132_SSURef_Nr99_tax_silva.fasta

#Replace commas with underscores
sed -i 's/,/_/g' SILVA_132_SSURef_Nr99_tax_silva.fasta

#Replace hyphens with underscores
sed -i 's/-/_/g' SILVA_132_SSURef_Nr99_tax_silva.fasta

#Remove Eukaryota, Chloroplast, and Mitochondrial sequences
grep -A 1 'Eukaryota' SILVA_132_SSURef_Nr99_tax_silva.fasta | awk '$1 != "--"' | grep -Fvxf - SILVA_132_SSURef_Nr99_tax_silva.fasta > temp.fasta
mv temp.fasta SILVA_132_SSURef_Nr99_tax_silva.fasta

grep -A 1 'Chloroplast' SILVA_132_SSURef_Nr99_tax_silva.fasta | awk '$1 != "--"' | grep -Fvxf - SILVA_132_SSURef_Nr99_tax_silva.fasta > temp.fasta
mv temp.fasta SILVA_132_SSURef_Nr99_tax_silva.fasta

grep -A 1 'Mitochondria' SILVA_132_SSURef_Nr99_tax_silva.fasta | awk '$1 != "--"' | grep -Fvxf - SILVA_132_SSURef_Nr99_tax_silva.fasta > temp.fasta
mv temp.fasta SILVA_132_SSURef_Nr99_tax_silva.fasta

#Fix some divergent lineage issues
sed -i 's/Bacillales;Family_XI/Bacillales;Bacillales_Family_XI/g' SILVA_132_SSURef_Nr99_tax_silva.fasta
sed -i 's/Clostridiales;Family_XI/Clostridiales;Clostridiales_Family_XI/g' SILVA_132_SSURef_Nr99_tax_silva.fasta

#Replace Unknown, uncultured, unidentified, uncharacterized, and metagenome entries with custom labels to be more uniform
split -l 48800 SILVA_132_SSURef_Nr99_tax_silva.fasta SPLIT

echo '#!/bin/bash' > bash_script.sh

echo "grep 'Unknown_Order;Unknown_Family' \$1 | awk '{print \$2}' | \\" >> bash_script.sh
echo "while read line; do" >> bash_script.sh
echo "  new_line=\$(echo \$line | sed 's/;/ /g' | awk '{\$4 = \$3\"_order\"; \$5 = \$3\"_family\"};1' | sed 's/ /;/g')" >> bash_script.sh
echo "  sed -i \"s/\$line/\$new_line/\" \$1" >> bash_script.sh
echo "done" >> bash_script.sh

echo "grep 'Unknown_Family' \$1 | awk '{print \$2}' | \\" >> bash_script.sh
echo "while read line; do" >> bash_script.sh
echo "  new_line=\$(echo \$line | sed 's/;/ /g' | awk '{\$5 = \$4\"_family\"};1' | sed 's/ /;/g')" >> bash_script.sh
echo "  sed -i \"s/\$line/\$new_line/\" \$1" >> bash_script.sh
echo "done" >> bash_script.sh

echo "grep 'uncultured' \$1 | awk '{print \$2}' | \\" >> bash_script.sh
echo "while read line; do" >> bash_script.sh
echo "  new_line=\$(echo \$line | sed 's/;/ /g' | \\" >> bash_script.sh
echo "              awk '\$1 ~ /uncultured/ {\$1 = \"Unknown_kingdom\"; \$2 = \"Unknown_phylum\"; \$3 = \"Unknown_class\"; \$4 = \"Unknown_order\"; \$5 = \"Unknown_family\"; \$6 = \"Unknown_genus\"; \$7 = \"Unknown_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$2 ~ /uncultured/ {\$2 = \$1\"_phylum\"; \$3 = \$1\"_class\"; \$4 = \$1\"_order\"; \$5 = \$1\"_family\"; \$6 = \$1\"_genus\"; \$7 = \$1\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$3 ~ /uncultured/ {\$3 = \$2\"_class\"; \$4 = \$2\"_order\"; \$5 = \$2\"_family\"; \$6 = \$2\"_genus\"; \$7 = \$2\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$4 ~ /uncultured/ {\$4 = \$3\"_order\"; \$5 = \$3\"_family\"; \$6 = \$3\"_genus\"; \$7 = \$3\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$5 ~ /uncultured/ {\$5 = \$4\"_family\"; \$6 = \$4\"_genus\"; \$7 = \$4\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$6 ~ /uncultured/ {\$6 = \$5\"_genus\"; \$7 = \$5\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$7 ~ /uncultured/ {\$7 = \$6\"_species\"};1' | sed 's/ /;/g')" >> bash_script.sh
echo "  sed -i \"s/\$line/\$new_line/\" \$1" >> bash_script.sh
echo "done" >> bash_script.sh

echo "grep 'unidentified' \$1 | awk '{print \$2}' | \\" >> bash_script.sh
echo "while read line; do" >> bash_script.sh
echo "  new_line=\$(echo \$line | sed 's/;/ /g' | \\" >> bash_script.sh
echo "              awk '\$1 ~ /unidentified/ {\$1 = \"Unknown_kingdom\"; \$2 = \"Unknown_phylum\"; \$3 = \"Unknown_class\"; \$4 = \"Unknown_order\"; \$5 = \"Unknown_family\"; \$6 = \"Unknown_genus\"; \$7 = \"Unknown_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$2 ~ /unidentified/ {\$2 = \$1\"_phylum\"; \$3 = \$1\"_class\"; \$4 = \$1\"_order\"; \$5 = \$1\"_family\"; \$6 = \$1\"_genus\"; \$7 = \$1\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$3 ~ /unidentified/ {\$3 = \$2\"_class\"; \$4 = \$2\"_order\"; \$5 = \$2\"_family\"; \$6 = \$2\"_genus\"; \$7 = \$2\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$4 ~ /unidentified/ {\$4 = \$3\"_order\"; \$5 = \$3\"_family\"; \$6 = \$3\"_genus\"; \$7 = \$3\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$5 ~ /unidentified/ {\$5 = \$4\"_family\"; \$6 = \$4\"_genus\"; \$7 = \$4\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$6 ~ /unidentified/ {\$6 = \$5\"_genus\"; \$7 = \$5\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$7 ~ /unidentified/ {\$7 = \$6\"_species\"};1' | sed 's/ /;/g')" >> bash_script.sh
echo "  sed -i \"s/\$line/\$new_line/\" \$1" >> bash_script.sh
echo "done" >> bash_script.sh

echo "grep 'uncharacterized' \$1 | awk '{print \$2}' | \\" >> bash_script.sh
echo "while read line; do" >> bash_script.sh
echo "  new_line=\$(echo \$line | sed 's/;/ /g' | \\" >> bash_script.sh
echo "              awk '\$1 ~ /uncharacterized/ {\$1 = \"Unknown_kingdom\"; \$2 = \"Unknown_phylum\"; \$3 = \"Unknown_class\"; \$4 = \"Unknown_order\"; \$5 = \"Unknown_family\"; \$6 = \"Unknown_genus\"; \$7 = \"Unknown_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$2 ~ /uncharacterized/ {\$2 = \$1\"_phylum\"; \$3 = \$1\"_class\"; \$4 = \$1\"_order\"; \$5 = \$1\"_family\"; \$6 = \$1\"_genus\"; \$7 = \$1\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$3 ~ /uncharacterized/ {\$3 = \$2\"_class\"; \$4 = \$2\"_order\"; \$5 = \$2\"_family\"; \$6 = \$2\"_genus\"; \$7 = \$2\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$4 ~ /uncharacterized/ {\$4 = \$3\"_order\"; \$5 = \$3\"_family\"; \$6 = \$3\"_genus\"; \$7 = \$3\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$5 ~ /uncharacterized/ {\$5 = \$4\"_family\"; \$6 = \$4\"_genus\"; \$7 = \$4\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$6 ~ /uncharacterized/ {\$6 = \$5\"_genus\"; \$7 = \$5\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$7 ~ /uncharacterized/ {\$7 = \$6\"_species\"};1' | sed 's/ /;/g')" >> bash_script.sh
echo "  sed -i \"s/\$line/\$new_line/\" \$1" >> bash_script.sh
echo "done" >> bash_script.sh

echo "grep 'metagenome' \$1 | awk '{print \$2}' | \\" >> bash_script.sh
echo "while read line; do" >> bash_script.sh
echo "  new_line=\$(echo \$line | sed 's/;/ /g' | \\" >> bash_script.sh
echo "              awk '\$1 ~ /metagenome/ {\$1 = \"Unknown_kingdom\"; \$2 = \"Unknown_phylum\"; \$3 = \"Unknown_class\"; \$4 = \"Unknown_order\"; \$5 = \"Unknown_family\"; \$6 = \"Unknown_genus\"; \$7 = \"Unknown_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$2 ~ /metagenome/ {\$2 = \$1\"_phylum\"; \$3 = \$1\"_class\"; \$4 = \$1\"_order\"; \$5 = \$1\"_family\"; \$6 = \$1\"_genus\"; \$7 = \$1\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$3 ~ /metagenome/ {\$3 = \$2\"_class\"; \$4 = \$2\"_order\"; \$5 = \$2\"_family\"; \$6 = \$2\"_genus\"; \$7 = \$2\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$4 ~ /metagenome/ {\$4 = \$3\"_order\"; \$5 = \$3\"_family\"; \$6 = \$3\"_genus\"; \$7 = \$3\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$5 ~ /metagenome/ {\$5 = \$4\"_family\"; \$6 = \$4\"_genus\"; \$7 = \$4\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$6 ~ /metagenome/ {\$6 = \$5\"_genus\"; \$7 = \$5\"_species\"};1' | \\" >> bash_script.sh
echo "              awk '\$7 ~ /metagenome/ {\$7 = \$6\"_species\"};1' | sed 's/ /;/g')" >> bash_script.sh
echo "  sed -i \"s/\$line/\$new_line/\" \$1" >> bash_script.sh
echo "done" >> bash_script.sh
chmod +x bash_script.sh

touch job_ids.txt
for file in SPLIT*; do
  sbatch --partition=short \
  --job-name=${file} \
  --error=${file}.err \
  --output=${file}.out \
  --time=12:00:00 \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem-per-cpu=64000 \
  ./bash_script.sh $file | awk '{print $4}' >> job_ids.txt
done

while squeue | awk '{print $1}' | grep -q -f job_ids.txt; do
  :
done

rm SPLIT*.err
rm SPLIT*.out
cat SPLIT* > SILVA_132_SSURef_Nr99_tax_silva.fasta
rm SPLIT*
rm bash_script.sh
rm job_ids.txt

#Remove sequences with taxonomy less than 7 fields
grep '^>' SILVA_132_SSURef_Nr99_tax_silva.fasta | awk '{print $2}' > temp.txt
echo "n_fields <- count.fields('temp.txt', sep=';')" > Rfunc.R
echo "taxa <- read.table('temp.txt', col.names=paste0('V', seq_len(max(n_fields))), fill=T, header=F, stringsAsFactors=F, comment.char='', sep=';', na.strings='')" >> Rfunc.R
echo "taxa <- taxa[rowSums(is.na(taxa)) > 0,]" >> Rfunc.R
echo "taxa[is.na(taxa)] <- ''" >> Rfunc.R
echo "write.table(taxa, 'temp.txt', sep=';', row.names=F, quote=F, col.names=F)" >> Rfunc.R
Rscript Rfunc.R
rm Rfunc.R
sed -i 's/\(;\)*$//' temp.txt
grep -A 1 -f temp.txt SILVA_132_SSURef_Nr99_tax_silva.fasta | awk '$1 != "--"' | grep -Fvxf - SILVA_132_SSURef_Nr99_tax_silva.fasta > temp.fasta
mv temp.fasta SILVA_132_SSURef_Nr99_tax_silva.fasta
rm temp.txt

#Run cutadapt to filter small sequences and extract V4 region of silva_nr_v132_train_set.fa
source /home/wallenz/miniconda3/etc/profile.d/conda.sh; conda activate base

#Remove small sequences
cutadapt -m 20 \
         -o temp.fasta \
         SILVA_132_SSURef_Nr99_tax_silva.fasta \
         > SILVA_132_SSURef_Nr99_tax_silva.log
mv temp.fasta SILVA_132_SSURef_Nr99_tax_silva.fasta

#Run cutadapt to trim for first primer sequence 515F
cutadapt -g GTGCCAGCNGCCGCGGTAA \
         -o temp.fasta \
         SILVA_132_SSURef_Nr99_tax_silva.fasta \
         > SILVA_132_SSURef_Nr99_tax_silva_V4.log

#Run cutadapt to trim for second primer sequence 806R
cutadapt -a ATTAGANACCCNNGTAGTCC \
         -o SILVA_132_SSURef_Nr99_tax_silva_V4.fasta \
         -m 20 \
         temp.fasta \
         >> SILVA_132_SSURef_Nr99_tax_silva_V4.log

rm temp.fasta

