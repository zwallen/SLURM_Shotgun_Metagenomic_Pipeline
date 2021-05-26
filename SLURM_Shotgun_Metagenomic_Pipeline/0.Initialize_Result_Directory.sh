#!/bin/bash
set -e

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 25 May 2021                                  #
#                                                            #
# Description: Initialize the directory that pipeline results#
# will be placed into.                                       #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#                                                            #
# Usage:                                                     #
# 0.Initialize_Result_Directory.sh -o output_dir             #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -o    (Required) Directory to put output of pipeline   #
#           into. NOTE: make sure output directory is in an  #
#           area that has plenty of data storage space       #
#           available if processing large datasets.          #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 25 May 2021                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":ho:" opt; do
  case $opt in
    h)
    echo " Description: Initialize the directory that pipeline results"
    echo " will be placed into.                                       "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " 0.Initialize_Result_Directory.sh -o output_dir             "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -o    (Required) Directory to put output of pipeline   "
    echo "           into. NOTE: make sure output directory is in an  "
    echo "           area that has plenty of data storage space       "
    echo "           available if processing large datasets.          "
    echo " "
    exit 0
    ;;
    o) OUT_DIR=$(echo $OPTARG | sed 's#/$##')
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
if [[ -z "$OUT_DIR" ]]; then
  echo "ERROR: Argument -o is required, please supply an output directory"
  exit 1
fi
if [[ ! -d "$OUT_DIR" ]]; then
  echo "ERROR: Argument -o should be a directory, please supply an output directory"
  exit 1
fi

###### CREATE DIRECTORY FOR PIPELINE OUTPUT #####
DATE=$(date | awk '{print $3"_"$2"_"$6}')
if [ -d "${OUT_DIR}/Metagenomic_Pipeline_${DATE}" ]
then
	:
else
	mkdir ${OUT_DIR}/Metagenomic_Pipeline_${DATE}
fi

echo " "
echo "Directory for shotgun metagenomic pipeline output: Metagenomic_Pipeline_${DATE}"
echo " "
RESULTS_DIR="${OUT_DIR}/Metagenomic_Pipeline_${DATE}"
