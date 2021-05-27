#!/bin/bash
set -e

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 27 May 2021                                  #
#                                                            #
# Description: Merge paired-end reads using BBMerge.         #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    BBMerge:    For merging paired-end reads.               #
#                                                            #
# Usage:                                                     #
# ./2.Merge_PE_Reads.sh -i input_seqs_dir \                  #
#                       -o output_dir \                      #
#                       -p 'commands; to; load; programs' \  #
#                       -f notificationEmail@forFailures.edu #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -i    (Required) Directory that contains the raw       #
#           fastq files to be processed. Sequences must have #
#           file extensions .fastq OR .fq,                   #
#           and can be gzipped or not.                       #
#     -o    (Required) Path to pipeline result directory.    #
#     -p    (Required) Single quoted string that contains    #
#           commands to load all the necessary programs      #
#           needed to run pipeline steps (e.g. activating    #
#           conda environments, loading modules, adding to   #
#           PATH, etc.).                                     #
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 27 May 2021                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hi:o:p:f:" opt; do
  case $opt in
    h)
    echo " Description: Merge paired-end reads using BBMerge.         "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    BBMerge:    For merging paired-end reads.               "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./2.Merge_PE_Reads.sh -i input_seqs_dir \                  "
    echo "                       -o output_dir \                      "
    echo "                       -p 'commands; to; load; programs' \  "
    echo "                       -f notificationEmail@forFailures.edu "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -i    (Required) Directory that contains the raw       "
    echo "           fastq files to be processed. Sequences must have "
    echo "           file extensions .fastq OR .fq,                   "
    echo "           and can be gzipped or not.                       "
    echo "     -o    (Required) Path to pipeline result directory.    "
    echo "     -p    (Required) Single quoted string that contains    "
    echo "           commands to load all the necessary programs      "
    echo "           needed to run pipeline steps (e.g. activating    "
    echo "           conda environments, loading modules, adding to   "
    echo "           PATH, etc.).                                     "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo " "
    exit 0
    ;;
    i) SEQ_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    o) RESULTS_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    p) PROG_LOAD="$OPTARG"
    ;;
    f) FAIL_EMAIL="$OPTARG"
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
if [[ -z "$SEQ_DIR" ]]; then
  echo "ERROR: Argument -i is required, please supply a directory with input fastq files"
  exit 1
fi
if [[ ! -d "$SEQ_DIR" ]]; then
  echo "ERROR: Argument -i should be a directory, please supply a directory with input fastq files"
  exit 1
fi
if ls -l $SEQ_DIR | grep -q ".fastq.gz"; then
  SEQ_EXT=fastq.gz
elif ls -l $SEQ_DIR | grep -q ".fastq"; then
  SEQ_EXT=fastq
elif ls -l $SEQ_DIR | grep -q ".fq.gz"; then
  SEQ_EXT=fq.gz
elif ls -l $SEQ_DIR | grep -q ".fq"; then
  SEQ_EXT=fq
else
  echo "ERROR: Sequences in input directory should have file extension of either .fastq[.gz] OR .fq[.gz]"
  exit 1
fi
found_1=$(ls -l $SEQ_DIR | grep -q "_R1_001")
found_2=$(ls -l $SEQ_DIR | grep -q "_R2_001")
if [[ -n "$found_1" ]] && [[ -n "$found_2" ]]; then
  echo "ERROR: Sequences in input directory expected to be paired-end Illumina sequence fastq files whose file names contain the strings '_R1_001' and '_R2_001'"
  exit 1
else
  :
fi

# -o
if [[ -z "$RESULTS_DIR" ]]; then
  echo "ERROR: Argument -o is required, please supply an output directory"
  exit 1
fi
if [[ ! -d "$RESULTS_DIR" ]]; then
  echo "ERROR: Argument -o should be a directory, please supply an output directory"
  exit 1
fi

# -p
if [[ -z "$PROG_LOAD" ]]; then
  echo "ERROR: Argument -p is required, please supply a single quoted string of commands needed to load required programs (can be an empty string ' ' if none required)"
  exit 1
fi

# -f
if [[ -z "$FAIL_EMAIL" ]]; then
  echo "ERROR: Argument -f is required, please supply an email that can be notified upon failure of any jobs ran during the pipeline"
  exit 1
elif echo $FAIL_EMAIL | grep -q -v '@'; then
  echo "ERROR: Argument -f requires a valid email, please give an email in the form of xxxx@xxxx.xxx"
fi

######### MERGE PAIRED READS WITH BBMERGE #######

  SECONDS=0
  echo " "
  echo "*** Running BBMerge to perform paired read merging ***"
  echo " "
  
  #Create directory for output
  if [ -d "${RESULTS_DIR}/1.Merged_Paired_End_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/1.Merged_Paired_End_Sequences
	  mkdir ${RESULTS_DIR}/1.Merged_Paired_End_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/1.Merged_Paired_End_Sequences/0.Output
  fi
  
  ##### Run BBMerge #####
  #Create script for running program and submit
  echo '#!/bin/bash' > bash_script.sh
  echo "#SBATCH --partition=short" >> bash_script.sh
  echo "#SBATCH --job-name=Merge" >> bash_script.sh
  echo "#SBATCH --error=${RESULTS_DIR}/1.Merged_Paired_End_Sequences/0.ErrorOut/Merge_%A_%a.err" >> bash_script.sh
  echo "#SBATCH --output=${RESULTS_DIR}/1.Merged_Paired_End_Sequences/0.Output/Merge_%A_%a.out" >> bash_script.sh
  echo "#SBATCH --time=12:00:00" >> bash_script.sh
  echo "#SBATCH --ntasks=1" >> bash_script.sh
  echo "#SBATCH --cpus-per-task=5" >> bash_script.sh
  echo "#SBATCH --mem-per-cpu=32000" >> bash_script.sh
  echo "#SBATCH --mail-type=FAIL" >> bash_script.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> bash_script.sh
  echo "#SBATCH --array=1-$(ls -l ${SEQ_DIR}/*R1_001.${SEQ_EXT} | wc -l)" >> bash_script.sh
  echo "#SBATCH --wait" >> bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE1=\$(ls ${SEQ_DIR}/*R1_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
  echo "FILE2=\$(ls ${SEQ_DIR}/*R2_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1_001' '{print \$1}')" >> bash_script.sh
  echo "bbmerge.sh in1=\$FILE1 \\" >> bash_script.sh
  echo "in2=\$FILE2 \\" >> bash_script.sh
  echo "out=${RESULTS_DIR}/1.Merged_Paired_End_Sequences/\${FILE_NAME}.fastq.gz \\" >> bash_script.sh
  echo "rem k=31 iterations=5 extend2=20 ecct t=5 -Xmx160g \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/1.Merged_Paired_End_Sequences/\${FILE_NAME}.log 2>&1"  >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh
  
  rm bash_script.sh
  
  #Signal jobs have ended
  echo "Merging of paired end reads with BBMerge complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "

#################################################
