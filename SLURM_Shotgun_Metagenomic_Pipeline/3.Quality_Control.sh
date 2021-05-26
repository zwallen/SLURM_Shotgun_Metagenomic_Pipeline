#!/bin/bash
set -e

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 25 May 2021                                  #
#                                                            #
# Description: Merge paired-end reads using BBMerge.         #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    BBDuk:      For adapter and quality trimming of raw wgs #
#                reads. Also can remove PhiX sequences.      #
#                                                            #
# Usage:                                                     #
# ./3.Quality_Control.sh -o output_dir \                     #
#                        -p 'commands; to; load; programs' \ #
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
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 25 May 2021                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":ho:p:f:" opt; do
  case $opt in
    h)
    echo " Description: Merge paired-end reads using BBMerge.         "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    BBDuk:      For adapter and quality trimming of raw wgs "
    echo "                reads. Also can remove PhiX sequences.      "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./3.Quality_Control.sh -o output_dir \                     "
    echo "                        -p 'commands; to; load; programs' \ "
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
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo " "
    exit 0
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

################# QC WITH BBDUK #################

  SECONDS=0
  echo " "
  echo "*** Running BBDuk to perform adapter/quality trimming and filtering on input fastq files ***"
  echo " "
  
  #Create directory for output
  if [ -d "${RESULTS_DIR}/2.Quality_Controlled_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/2.Quality_Controlled_Sequences
	  mkdir ${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.Output
  fi
  
  ##### Run BBDuk #####
  #Create shell script for running program
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$1 | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> bash_script.sh
  echo "FILE=\$(echo \$1 | awk -F '/' '{print \$NF}')" >> bash_script.sh
  echo "bbduk.sh in=\$1 \\" >> bash_script.sh
  echo "out=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE} \\" >> bash_script.sh
  echo "stats=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}_stats.txt \\" >> bash_script.sh
  echo "ftm=5 qtrim=rl trimq=25 minlen=50 ref=adapters,phix -Xmx64g \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}.log 2>&1"  >> bash_script.sh
  chmod +x bash_script.sh
  
  #For every sequence file submit job and grab job IDs
  touch job_ids.txt
  for file in ${RESULTS_DIR}/1.Merged_Paired_End_Sequences/*.fastq.gz; do
    FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '.fastq.gz' '{print $1}')
    
    sbatch --partition=short \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.ErrorOut/${FILE_NAME}.err \
    --output=${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.Output/${FILE_NAME}.out \
    --time=12:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script.sh $file | awk '{print $4}' >> job_ids.txt
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
  rm bash_script.sh
  rm job_ids.txt
  
  #Signal jobs have ended
  echo "Quality control with BBDuk complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "

#################################################
