#!/bin/bash
set -e

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 29 May 2021                                  #
#                                                            #
# Description: Remove host sequence reads using KneadData.   #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    KneadData:  For removing host contamination from wgs    #
#                reads. Requires Bowtie2 database file to    #
#                map reads against.                          #
#                                                            #
# Usage:                                                     #
# ./4.Remove_Host_Reads.sh -o output_dir \                   #
#                        -p 'commands; to; load; programs' \ #
#                        -r path/to/host/ref/files/dir \     #
#                        -f notificationEmail@forFailures.edu#
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -o    (Required) Path to pipeline result directory.    #
#     -p    (Required) Single quoted string that contains    #
#           commands to load all the necessary programs      #
#           needed to run pipeline steps (e.g. activating    #
#           conda environments, loading modules, adding to   #
#           PATH, etc.).                                     #
#     -r    (Required) Path to directory of host genome      #
#           Bowtie2 indexed reference files (.bt2 files).    #
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 29 May 2021                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":ho:p:r:f:" opt; do
  case $opt in
    h)
    echo " Description: Merge paired-end reads using BBMerge.         "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    KneadData:  For removing host contamination from wgs    "
    echo "                reads. Requires Bowtie2 database file to    "
    echo "                map reads against.                          "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./4.Remove_Host_Reads.sh -o output_dir \                   "
    echo "                        -p 'commands; to; load; programs' \ "
    echo "                        -r path/to/host/ref/files/dir \     "
    echo "                        -f notificationEmail@forFailures.edu"
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -o    (Required) Path to pipeline result directory.    "
    echo "     -p    (Required) Single quoted string that contains    "
    echo "           commands to load all the necessary programs      "
    echo "           needed to run pipeline steps (e.g. activating    "
    echo "           conda environments, loading modules, adding to   "
    echo "           PATH, etc.).                                     "
    echo "     -r    (Required) Path to directory of host genome      "
    echo "           Bowtie2 indexed reference files (.bt2 files).    "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo " "
    exit 0
    ;;
    o) RESULTS_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    p) PROG_LOAD="$OPTARG"
    ;;
    r) HOST_REF=$(echo $OPTARG | sed 's#/$##')
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

# -r
if [[ -z "$HOST_REF" ]]; then
  echo "ERROR: Argument -r is required, please supply a directory with host genome Bowtie2 indexed reference sequences"
  exit 1
fi
if [[ ! -d "$HOST_REF" ]]; then
  echo "ERROR: Argument -r should be a directory, please supply a directory with host genome Bowtie2 indexed reference sequences"
  exit 1
fi
if ls -l $HOST_REF | grep -q ".bt2"; then
  :
else
  echo "ERROR: Expecting host genome Bowtie2 indexed reference files to have extension .bt2, but none found"
  exit 1
fi

# -f
if [[ -z "$FAIL_EMAIL" ]]; then
  echo "ERROR: Argument -f is required, please supply an email that can be notified upon failure of any jobs ran during the pipeline"
  exit 1
elif echo $FAIL_EMAIL | grep -q -v '@'; then
  echo "ERROR: Argument -f requires a valid email, please give an email in the form of xxxx@xxxx.xxx"
fi

##### REMOVAL OF HOST READS WITH KNEADDATA ######

  SECONDS=0
  echo " "
  echo "*** Running KneadData to perform removal of contaminant host reads ***"
  echo " "
  
  #Create directory for quality controlled sequences
  if [ -d "${RESULTS_DIR}/3.Decontaminated_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences
	  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences/0.Output
  fi
  
  ##### Run KneadData #####
  #Create script for running program and submit
  echo '#!/bin/bash' > bash_script.sh
  echo "#SBATCH --partition=short" >> bash_script.sh
  echo "#SBATCH --job-name=Decontam" >> bash_script.sh
  echo "#SBATCH --error=${RESULTS_DIR}/3.Decontaminated_Sequences/0.ErrorOut/Decontam_%A_%a.err" >> bash_script.sh
  echo "#SBATCH --output=${RESULTS_DIR}/3.Decontaminated_Sequences/0.Output/Decontam_%A_%a.out" >> bash_script.sh
  echo "#SBATCH --time=12:00:00" >> bash_script.sh
  echo "#SBATCH --ntasks=1" >> bash_script.sh
  echo "#SBATCH --cpus-per-task=2" >> bash_script.sh
  echo "#SBATCH --mem-per-cpu=32000" >> bash_script.sh
  echo "#SBATCH --mail-type=FAIL" >> bash_script.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> bash_script.sh
  echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*.fastq.gz | wc -l)" >> bash_script.sh
  echo "#SBATCH --wait" >> bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE=\$(ls ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> bash_script.sh
  echo "kneaddata --input \$FILE \\" >> bash_script.sh
  echo "--output ${RESULTS_DIR}/3.Decontaminated_Sequences \\" >> bash_script.sh
  echo "--output-prefix \$FILE_NAME \\" >> bash_script.sh
  echo "--log ${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}_kneaddata.log \\" >> bash_script.sh
  echo "--reference-db $HOST_REF \\" >> bash_script.sh
  echo "--bypass-trim \\" >> bash_script.sh
  echo "--threads 2 \\" >> bash_script.sh
  echo "--verbose \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh
  
  rm bash_script.sh
  
  ##### Gzip output #####
  echo "Compressing KneadData output..."
  echo " "
  
  #Launch gzip for each file
  touch job_ids.txt
  for file in ${RESULTS_DIR}/3.Decontaminated_Sequences/*.fastq; do
    sbatch --partition=short \
           --job-name=Gzip \
	   --error=${RESULTS_DIR}/3.Decontaminated_Sequences/0.ErrorOut/Gzip.err \
	   --output=${RESULTS_DIR}/3.Decontaminated_Sequences/0.ErrorOut/Gzip.out \
	   --time=12:00:00 \
	   --ntasks=1 \
	   --cpus-per-task=1 \
	   --mem-per-cpu=16000 \
	   --mail-type=FAIL \
	   --mail-user=${FAIL_EMAIL} \
	   --wrap="gzip $file" | \
	   awk '{print $4}' >> job_ids.txt
	   sleep 1
  done
  
  #Hold script here until all jobs are completed
  while :
  do
    if squeue -u $USER 2>&1 | grep -q -f job_ids.txt; then
      sleep 1m
      :
    elif squeue -u $USER 2>&1 | grep -q "slurm_load_jobs error"; then
      sleep 5m
      :
    else
      break
    fi
  done
  rm job_ids.txt
  
  #Move extracted host sequences to their own directory
  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences/Extracted_Host_Sequences
  mv ${RESULTS_DIR}/3.Decontaminated_Sequences/*contam* ${RESULTS_DIR}/3.Decontaminated_Sequences/Extracted_Host_Sequences/
  
  echo "Done"
  echo " "
  
  #Signal jobs have ended
  echo "Running of KneadData complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################
