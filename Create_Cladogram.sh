#!/bin/bash
set -e

##############################################################
# Create cladogram of most abundant clades                   #
# by Zachary D Wallen                                        #
# Last updated: 18 March 2021                                #
#                                                            #
# Description: This is a wrapper program that wraps          #
# GraPhlAn program to produce a cladogram summarizing the    #
# top most abundant clades detected by MetaPhlAn for a       #
# dataset.                                                   #
#                                                            #
# Required programs and databases:                           #
#    GraPhlAn:   For generating the cladogram. Will also need#
#                the export2graphlan.py conversion script.   #
#                                                            #
# Usage:                                                     #
# Create_Cladogram.sh -i input_metaphlan_rel_abun_table \    #
#                     -o output_cladogram \                  #
#                     -s "search_string"                     #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -i    (Required) A taxa by sample relative abundance   #
#           table created by running MetaPhlAn with "rel_ab" #
#           option on individual sample sequences, then      #
#           merging individual sample relative abundance     #
#           tables with merge_metaphlan_tables.py script.    #
#     -o    (Required) Name for the outputted cladogram.     #
#     -s    (Optional) Should the relative abundance table be#
#           subsetted for specific samples? If so, supply a  #
#           quoted pattern in the sample IDs here that can be#
#           used to extract the desired sample data columns. #
##############################################################

echo " "
echo "##############################################################"
echo "# Create cladogram of most abundant clades                   #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 18 March 2021                                #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hi:o:s:" opt; do
  case $opt in
    h)
    echo " Description: This is a wrapper program that wraps          "
    echo " GraPhlAn program to produce a cladogram summarizing the    "
    echo " top most abundant clades detected by MetaPhlAn for a       "
    echo " dataset.                                                   "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    GraPhlAn:   For generating the cladogram. Will also need"
    echo "                the export2graphlan.py conversion script.   "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " Create_Cladogram.sh -i input_metaphlan_rel_abun_table \    "
    echo "                     -o output_cladogram \                  "
    echo '                     -s "search_string"                     '
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -i    (Required) A taxa by sample relative abundance   "
    echo "           table created by running MetaPhlAn with 'rel_ab' "
    echo "           option on individual sample sequences, then      "
    echo "           merging individual sample relative abundance     "
    echo "           tables with merge_metaphlan_tables.py script.    "
    echo "     -o    (Required) Name for the outputted cladogram.     "
    echo "     -s    (Optional) Should the relative abundance table be"
    echo "           subsetted for specific samples? If so, supply a  "
    echo "           quoted pattern in the sample IDs here that can be"
    echo "           used to extract the desired sample data columns. "
    exit 0
    ;;
    i) IN_FILE="$OPTARG"
    ;;
    o) OUT_FILE="$OPTARG"
    ;;
    s) SUBSET="$OPTARG"
    ;;
    \?) echo "Invalid option: $OPTARG" 1>&2
        exit 1
    ;;
    :) echo "Invalid option: $OPTARG requires an argument" 1>&2
       exit 1
    ;;
  esac
done

# Check that valid arguments were entered

# -i
if [[ -z "$IN_FILE" ]]; then
  echo "ERROR: Argument -i is required, please supply a MetaPhlAn relative abundance table"
  exit 1
fi
if [[ -d "$IN_FILE" ]]; then
  echo "ERROR: Argument -i should be a single file, but a directory has been supplied, please supply a MetaPhlAn relative abundance table"
  exit 1
fi
if [[ ! -f "$IN_FILE" ]]; then
  echo "ERROR: File given to -i could not be found"
  exit 1
fi

# -o
if [[ -z "$OUT_FILE" ]]; then
  echo "ERROR: Argument -o is required, please supply a name for the output cladogram"
  exit 1
fi
if [[ -f "$OUT_FILE" ]]; then
  echo "WARNING: a file has been found that matches the filename given to -o, this file will be overwritten"
fi

# -s
if grep -q "$SUBSET" $IN_FILE; then
  :
else
  echo "ERROR: Pattern given to -s resulted in no hits in the input file, does the pattern exist in the sample column names?"
  exit 1
fi

##### SUBSET RELATIVE ABUNDANCE TABLE #####
if [[ ! -z "$SUBSET" ]]; then
  echo "rel_ab <- read.table('$IN_FILE', header=T, stringsAsFactors=F)" > Rfunc.R
  echo "rel_ab_filt <- rel_ab[,c(1,2,grep('$SUBSET',colnames(rel_ab)))]" >> Rfunc.R
  echo "rel_ab_filt <- rel_ab_filt[rowSums(rel_ab_filt[,3:ncol(rel_ab_filt)] > 0) > 0,]" >> Rfunc.R
  echo "write.table(rel_ab_filt, 'temp.txt', sep='\t', quote=F, row.names=F)" >> Rfunc.R
else
  echo "rel_ab <- read.table('$IN_FILE', header=T, stringsAsFactors=F)" > Rfunc.R
  echo "rel_ab_filt <- rel_ab[rowSums(rel_ab[,3:ncol(rel_ab)] > 0) > 0,]" >> Rfunc.R
  echo "write.table(rel_ab_filt, 'temp.txt', sep='\t', quote=F, row.names=F)" >> Rfunc.R
fi
Rscript Rfunc.R
rm Rfunc.R

##### CREATE CLADOGRAM FOR TOP 10% MOST ABUNDANT CLADES #####

#Get N for 10% of total taxa detected
TAXA_N=$(wc -l temp.txt | awk '{printf "%.0f\n", ($1-2)*0.1}')

#Export annotation and tree files from metaphlan results
export2graphlan.py -i temp.txt \
--tree tree.txt \
--annotation annot.txt \
--most_abundant $TAXA_N \
--abundance_threshold 1 \
--least_biomarkers 10 \
--annotations 6 \
--external_annotations 7 \
--min_clade_size 1 \
--min_font_size 2 \
--max_font_size 2

#Get number of background color annotations
ANNOT_N=$(grep "annotation_background_color" annot.txt | wc -l)

#Create file with number of colors that match how many background color annotations were found
echo "col.df <- as.data.frame(rainbow($ANNOT_N))" > Rfunc.R
echo "write.table(col.df, 'annot_color.txt', sep='\t', quote=F, row.names=F, col.names=F)" >> Rfunc.R
Rscript Rfunc.R
rm Rfunc.R

#Add generated colors to the annotation file for background colors and corresponding clade marker colors
for i in $(seq 1 $ANNOT_N)
do
  ANNOT=$(grep "annotation_background_color" annot.txt | sed -ne ${i}p)
  MARKER_ANNOT=$(echo $ANNOT | awk '{print $1,"clade_marker_color",$3}' OFS='\t')
  COLOR=$(cat annot_color.txt | sed -ne ${i}p)
  NEW_ANNOT=$(echo $ANNOT | awk -v color=$COLOR '{print $1,$2,color}' OFS='\t')
  NEW_MARKER_ANNOT=$(echo $MARKER_ANNOT | awk -v color=$COLOR '{print $1,$2,color}' OFS='\t')

  sed -i "s/${ANNOT}/${NEW_ANNOT}/" annot.txt
  sed -i "s/${MARKER_ANNOT}/${NEW_MARKER_ANNOT}/" annot.txt
done

#Remove clade marker colors if they don't have a corresponding background
grep "clade_marker_color" annot.txt | while read line
do
  ANNOT=$(echo $line | awk '{print $1,"annotation_background_color",$3}' OFS='\t')

  if grep -q "$ANNOT" annot.txt
  then
    :
  else
    sed -i "s/${line}//" annot.txt
  fi
done

#Add parameters to annotation file to rotate outer annotations 90 degrees to avoid running into each other
grep "g__.*annotation_font_size" annot.txt | while read line
do
  CLADE=$(echo $line | awk '{print $1}')
  sed -i "/${line}/a ${CLADE}\tannotation_rotation\t90" annot.txt
done

#Remove Order annotations
grep -E '\b[[:upper:]]+\b' annot.txt | while read line
do
  sed -i "s/${line}//" annot.txt
done

#Plot graph
graphlan_annotate.py --annot annot.txt tree.txt annot.xml

graphlan.py --dpi 500 annot.xml ${OUT_FILE}.png --size 7 --pad 0.1 --external_legends

#Clean up
rm temp.txt
rm tree.txt
rm annot.txt
rm annot_color.txt
rm annot.xml

