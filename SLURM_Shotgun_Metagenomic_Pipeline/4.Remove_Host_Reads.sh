#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 6 Sep 2021                                   #
#                                                            #
# Description: Remove host sequence reads using KneadData.   #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    BBMap/                                                  #
#    BBSplit:    For removing host contamination from wgs    #
#                reads. Requires FASTA genome reference file #
#                to map reads against.                       #
#                                                            #
# Usage:                                                     #
# ./4.Remove_Host_Reads.sh [-m] \                            #
#                        -o output_dir \                     #
#                        -p 'commands; to; load; programs' \ #
#                        -r path/to/host/ref/file.fa \       #
#                        -f notificationEmail@forFailures.edu#
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -m    (Required) Are input fastq files merged? Add this#
#           parameter if so. Will cause error otherwise.     #
#     -o    (Required) Path to pipeline result directory.    #
#     -p    (Required) Single quoted string that contains    #
#           commands to load all the necessary programs      #
#           needed to run pipeline steps (e.g. activating    #
#           conda environments, loading modules, adding to   #
#           PATH, etc.).                                     #
#     -r    (Required) Path to reference genome file of host.#
#           Should be in FASTA format, and uncompressed.     #
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 6 Sep 2021                                   #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hmo:p:r:f:" opt; do
  case $opt in
    h)
    echo " Description: Merge paired-end reads using BBMerge.         "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    BBMap/                                                  "
    echo "    BBSplit:    For removing host contamination from wgs    "
    echo "                reads. Requires FASTA genome reference file "
    echo "                to map reads against.                       "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./4.Remove_Host_Reads.sh [-m] \                            "
    echo "                        -o output_dir \                     "
    echo "                        -p 'commands; to; load; programs' \ "
    echo "                        -r path/to/host/ref/file.fa \       "
    echo "                        -f notificationEmail@forFailures.edu"
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -m    (Required) Are input fastq files merged? Add this"
    echo "           parameter if so. Will cause error otherwise.     "
    echo "     -o    (Required) Path to pipeline result directory.    "
    echo "     -p    (Required) Single quoted string that contains    "
    echo "           commands to load all the necessary programs      "
    echo "           needed to run pipeline steps (e.g. activating    "
    echo "           conda environments, loading modules, adding to   "
    echo "           PATH, etc.).                                     "
    echo "     -r    (Required) Path to reference genome file of host."
    echo "           Should be in FASTA format, and uncompressed.     "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo " "
    exit 0
    ;;
    m) MERGE=1
    ;;
    o) RESULTS_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    p) PROG_LOAD="$OPTARG"
    ;;
    r) HOST_REF="$OPTARG"
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
  echo "ERROR: Argument -r is required, please supply a host genome reference file in FASTA format"
  exit 1
fi
if [[ -d "$HOST_REF" ]]; then
  echo "ERROR: Argument -r should be a single FASTA file, please supply a host genome reference file in FASTA format"
  exit 1
fi
if ls -l $HOST_REF | grep -q ".fa"; then
  :
elif ls -l $HOST_REF | grep -q ".fna"; then
  :
elif ls -l $HOST_REF | grep -q ".fasta"; then
  :
else
  echo "ERROR: Expecting host genome reference file to be in FASTA format with extension '.fa', '.fna', or '.fasta'"
  exit 1
fi

# -f
if [[ -z "$FAIL_EMAIL" ]]; then
  echo "ERROR: Argument -f is required, please supply an email that can be notified upon failure of any jobs ran during the pipeline"
  exit 1
elif echo $FAIL_EMAIL | grep -q -v '@'; then
  echo "ERROR: Argument -f requires a valid email, please give an email in the form of xxxx@xxxx.xxx"
fi

##### REMOVAL OF HOST READS WITH BBMAP/BBSPLIT ######

  SECONDS=0
  echo " "
  echo "*** Running BBMap/BBSplit to perform removal of contaminant host reads ***"
  echo " "
  
  #Create directory
  if [ -d "${RESULTS_DIR}/3.Decontaminated_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences
	  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences/0.Output
  fi
  
  ##### Run BBMap/BBSplit #####
  #Create script for running program and submit
  echo '#!/bin/bash' > bash_script.sh
  echo "#SBATCH --partition=short" >> bash_script.sh
  echo "#SBATCH --job-name=Decontam" >> bash_script.sh
  echo "#SBATCH --error=${RESULTS_DIR}/3.Decontaminated_Sequences/0.ErrorOut/Decontam_%A_%a.err" >> bash_script.sh
  echo "#SBATCH --output=${RESULTS_DIR}/3.Decontaminated_Sequences/0.Output/Decontam_%A_%a.out" >> bash_script.sh
  echo "#SBATCH --time=12:00:00" >> bash_script.sh
  echo "#SBATCH --ntasks=1" >> bash_script.sh
  echo "#SBATCH --cpus-per-task=10" >> bash_script.sh
  echo "#SBATCH --mem-per-cpu=16000" >> bash_script.sh
  echo "#SBATCH --mail-type=FAIL" >> bash_script.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> bash_script.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*.fastq.gz | wc -l)" >> bash_script.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*R1_001.fastq.gz | wc -l)" >> bash_script.sh
  fi
  echo "#SBATCH --wait" >> bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> bash_script.sh
    echo "bbsplit.sh in=\${FILE} \\" >> bash_script.sh
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*R1_001.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*R2_001.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1_001' '{print \$1}')" >> bash_script.sh
    echo "bbsplit.sh in1=\${FILE1} \\" >> bash_script.sh
    echo "in2=\${FILE2} \\" >> bash_script.sh
  fi
  echo "path=${RESULTS_DIR}/3.Decontaminated_Sequences \\" >> bash_script.sh
  echo "ref=${HOST_REF} \\" >> bash_script.sh
  echo "outu1=${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}_R1_001.fastq \\" >> bash_script.sh
  echo "outu2=${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}_R2_001.fastq \\" >> bash_script.sh
  echo "basename=${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}.%_contam_#.fastq \\" >> bash_script.sh
  echo "t=10 -Xmx160g \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh > /dev/null
  
  rm bash_script.sh
  rm -r ${RESULTS_DIR}/3.Decontaminated_Sequences/ref
  
  ##### Gzip output #####
  echo "Compressing BBMap/BBSplit output..."
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
	   sleep 30s
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
  echo "Running of BBMap/BBSplit complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
#################################################
