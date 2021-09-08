#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 8 Sep 2021                                   #
#                                                            #
# Description: Create tracking report of sequence read       #
# numbers through pipeline.                                  #
#                                                            #
# Usage:                                                     #
# ./6.Pipeline_Report.sh -o output_dir                       #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -o    (Required) Pipeline output directory.            #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 8 Sep 2021                                   #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":ho:" opt; do
  case $opt in
    h)
    echo " Description: Create tracking report of sequence read       "
    echo " numbers through pipeline.                                  "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./6.Pipeline_Report.sh -o output_dir                       "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -o    (Required) Pipeline output directory.            "
    echo " "
    exit 0
    ;;
    o) RESULTS_DIR=$(echo $OPTARG | sed 's#/$##')
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

# -o
if [[ -z "$RESULTS_DIR" ]]; then
  echo "ERROR: Argument -o is required, please supply path to pipeline output directory"
  exit 1
fi
if [[ ! -d "$RESULTS_DIR" ]]; then
  echo "ERROR: Argument -o should be a directory, please supply path to pipeline output directory"
  exit 1
fi

###### CREATE TRACKING REPORT #####

SECONDS=0
echo " "
echo "*** Creating tracking report of sequences through pipeline ***"
echo " "

echo "Sample Input Quality_controlled Decontaminated Profiled" | sed 's/ /\t/g' > ${RESULTS_DIR}/5.Pipeline_Report.txt

ls -l -d ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*_Profiles | awk '{print $NF}' | awk -F"/" '{print $NF}' | awk -F"_Profiles" '{print $1}' | \
while read sample; do
  echo "Processing $sample"
  BEGIN_SEQ=$(grep 'Input:' ${RESULTS_DIR}/2.Quality_Controlled_Sequences/${sample}.log | awk '{print $2}')
  QC_SEQ=$(grep 'Result:' ${RESULTS_DIR}/2.Quality_Controlled_Sequences/${sample}.log | awk '{print $2}')
  DECONTAM_SEQ=$(( $(grep 'READ COUNT: final pair1' ${RESULTS_DIR}/3.Decontaminated_Sequences/${sample}_kneaddata.log | awk '{print $NF}' | awk -F'.' '{print $1}') + \
               $(grep 'READ COUNT: final pair2' ${RESULTS_DIR}/3.Decontaminated_Sequences/${sample}_kneaddata.log | awk '{print $NF}' | awk -F'.' '{print $1}') ))
  PROFILED_SEQ=$(( $DECONTAM_SEQ - \
                   $(grep 'UNMAPPED' ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${sample}_Profiles/${sample}_genefamilies.tsv | awk '{print $NF}' | awk -F'.' '{print $1}') 
                 ))
  echo "$sample $BEGIN_SEQ $QC_SEQ $DECONTAM_SEQ $PROFILED_SEQ" | sed 's/ /\t/g' >> ${RESULTS_DIR}/5.Pipeline_Report.txt
done

echo " "
echo "Pipeline sequence tracking report complete"
echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
echo " "
#################################################
