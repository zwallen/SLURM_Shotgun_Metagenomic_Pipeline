#!/bin/bash
set -e

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 18 March 2021                                #
#                                                            #
# Description: This is a wrapper program that wraps various  #
# programs to process raw paired-end whole genome shotgun    #
# metagenomic sequences. The end product of this pipeline    #
# is taxonomic and functional (gene and pathway) abundances  #
# that are ready for further statistical analyses.           #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    FastQC:     For performing initial quality reports.     #
#    BBDuk:      For adapter and quality trimming of raw wgs #
#                reads. Also can remove PhiX sequences.      #
#    KneadData:  For removing host contamination from wgs    #
#                reads. Requires Bowtie2 database file to    #
#                map reads against.                          #
#    HUMAnN:     For generating taxonomic, gene family, and  #
#                pathway abundances.                         #
#    ChocoPhlAn: Database used for taxonomic profiling. Can  #
#                be any of the ChocoPhlAn databases          #
#                downloaded using humann_databases utility   #
#                program.                                    #
#    UniRef:     Database used for functional profiling. Can #
#                be any of the UniRef databases downloaded   #
#                using humann_databases utility program.     #
#    Markers:    File with ChocoPhlAn GeneIDs and marker     #
#                taxonomies. Used for replacing GeneIDs with #
#                clade marker taxonomy lineages in the       #
#                normalized abundance tables from MetaPhlAn. #
#                                                            #
# Optional programs and databases:                           #
#    SortMeRNA:  For extracting 16S rRNA gene sequences from #
#                wgs reads. Only needed if running the       #
#                additional 16S based classification         #
#                pipeline. Requires a 16S FASTA reference    #
#                file to use for extracting 16S sequences.   #
#    RDP classifier:                                         #
#                For classifying extracted 16S rRNA gene     #
#                sequences. Only needed if running the       #
#                additional 16S based classification         #
#                pipeline. Requires a 16S FASTA reference    #
#                to use for training the classifier and      #
#                classifying sequences.                      #
#                                                            #
# Usage:                                                     #
# SLURM_Shotgun_Metagenomic_Pipeline.sh -i input_seqs_dir \  #
#                    -o output_dir \                         #
#                    -p 'commands; to; load; programs' \     #
#                    -r path/to/host/ref/files/dir \         #
#                    -c path/to/chocophlan/dir \             #
#                    -u path/to/uniref/dir \                 #
#                    -m path/to/clade/marker/info/file \     #
#                    -f notificationEmail@forFailures.edu \  #
#                    [additional options]                    #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -i    (Required) Directory that contains the raw       #
#           fastq files to be processed. Sequences must have #
#           file extensions .fastq OR .fq,                   #
#           and can be gzipped or not.                       #
#     -o    (Required) Directory to put output of pipeline   #
#           into. NOTE: make sure output directory is in an  #
#           area that has plenty of data storage space       #
#           available if processing large datasets.          #
#     -p    (Required) Single quoted string that contains    #
#           commands to load all the necessary programs      #
#           needed to run pipeline steps (e.g. activating    #
#           conda environments, loading modules, adding to   #
#           PATH, etc.).                                     #
#     -r    (Required) Path to directory of host genome      #
#           Bowtie2 indexed reference files (.bt2 files).    #
#     -c    (Required) Path to ChocoPhlAn database directory.#
#     -u    (Required) Path to UniRef90 database directory.  #
#     -m    (Required) Path to clade marker info file        #
#           mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2    #
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
#     -a    (Optional) Perform additional 16S rRNA gene      #
#           based taxonomy profiling.                        #
#     -e    (Optional) 16S rRNA DNA reference FASTA file to  #
#           use for extraction of 16S rRNA DNA sequences.    #
#           Should be unzipped FASTA file with extension     #
#           fasta, fa, or fna.                               #
#     -t    (Optional) 16S rRNA DNA reference FASTA file to  #
#           use for training the RDP classifier and          #
#           classifying extracted 16S sequences. Should be   #
#           unzipped FASTA file with extension fasta, fa, or #
#           fna. Should be of the following format:          #
#                                                            #
#     >Seq_ID Kingdom;Phylum;Class;Order;Family;Genus;Species#
#     ACTGAAACTGCCGTTCAAAGCTTCGCGCGCTTTCCCGGGCGCGATATACGCGCGC#
#     AAACTGGGGTCGCGAA                                       #
#                                                            #
#     -l    (Optional) Location of RDP classifier source dir #
#           (rdp_classifier_x.xx). Safest to give full path  #
#           to the dir, but a partial path should also work. #
#     -s    (Optional) Skip certain steps in the pipeline if #
#           need be. Provide a comma separated list of steps #
#           that you wish to skip in the pipeline. List may  #
#           have the values: fastqc, bbduk, kneaddata,       #
#           humann.                                          #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 18 March 2021                                #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hi:o:p:r:c:u:m:f:ae:t:l:s:" opt; do
  case $opt in
    h)
    echo " Description: This is a wrapper program that wraps various  "
    echo " programs to process raw paired-end whole genome shotgun    "
    echo " metagenomic sequences. The end product of this pipeline    "
    echo " is taxonomic and functional (gene and pathway) abundances  "
    echo " that are ready for further statistical analyses.           "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    FastQC:     For performing initial quality reports.     "
    echo "    BBDuk:      For adapter and quality trimming of raw wgs "
    echo "                reads. Also can remove PhiX sequences.      "
    echo "    KneadData:  For removing host contamination from wgs    "
    echo "                reads. Requires Bowtie2 database file to    "
    echo "                map reads against.                          "
    echo "    HUMAnN:     For generating taxonomic, gene family, and  "
    echo "                pathway abundances.                         "
    echo "    ChocoPhlAn: Database used for taxonomic profiling. Can  "
    echo "                be any of the ChocoPhlAn databases          "
    echo "                downloaded using humann_databases utility   "
    echo "                program.                                    "
    echo "    UniRef:     Database used for functional profiling. Can "
    echo "                be any of the UniRef databases downloaded   "
    echo "                using humann_databases utility program.     "
    echo "    Markers:    File with ChocoPhlAn GeneIDs and marker     "
    echo "                taxonomies. Used for replacing GeneIDs with "
    echo "                clade marker taxonomy lineages in the       "
    echo "                normalized abundance tables from MetaPhlAn. "
    echo "                                                            "
    echo " Optional programs and databases:                           "
    echo "    SortMeRNA:  For extracting 16S rRNA gene sequences from "
    echo "                wgs reads. Only needed if running the       "
    echo "                additional 16S based classification         "
    echo "                pipeline. Requires a 16S FASTA reference    "
    echo "                file to use for extracting 16S sequences.   "
    echo "    RDP classifier:                                         "
    echo "                For classifying extracted 16S rRNA gene     "
    echo "                sequences. Only needed if running the       "
    echo "                additional 16S based classification         "
    echo "                pipeline. Requires a 16S FASTA reference    "
    echo "                to use for training the classifier and      "
    echo "                classifying sequences.                      "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " SLURM_Shotgun_Metagenomic_Pipeline.sh -i input_seqs_dir \  "
    echo "                    -o output_dir \                         "
    echo "                    -p 'commands; to; load; programs' \     "
    echo "                    -r path/to/host/ref/files/dir \         "
    echo "                    -c path/to/chocophlan/dir \             "
    echo "                    -u path/to/uniref/dir \                 "
    echo "                    -m path/to/clade/marker/info/file \     "
    echo "                    -f notificationEmail@forFailures.edu \  "
    echo "                    [additional options]                    "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -i    (Required) Directory that contains the raw       "
    echo "           fastq files to be processed. Sequences must have "
    echo "           file extensions .fastq OR .fq,                   "
    echo "           and can be gzipped or not.                       "
    echo "     -o    (Required) Directory to put output of pipeline   "
    echo "           into. NOTE: make sure output directory is in an  "
    echo "           area that has plenty of data storage space       "
    echo "           available if processing large datasets.          "
    echo "     -p    (Required) Single quoted string that contains    "
    echo "           commands to load all the necessary programs      "
    echo "           needed to run pipeline steps (e.g. activating    "
    echo "           conda environments, loading modules, adding to   "
    echo "           PATH, etc.).                                     "
    echo "     -r    (Required) Path to directory of host genome      "
    echo "           Bowtie2 indexed reference files (.bt2 files).    "
    echo "     -c    (Required) Path to ChocoPhlAn database directory."
    echo "     -u    (Required) Path to UniRef90 database directory.  "
    echo "     -m    (Required) Path to clade marker info file        "
    echo "           mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2    "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo "     -a    (Optional) Perform additional 16S rRNA gene      "
    echo "           based taxonomy profiling.                        "
    echo "     -e    (Optional) 16S rRNA DNA reference FASTA file to  "
    echo "           use for extraction of 16S rRNA DNA sequences.    "
    echo "           Should be unzipped FASTA file with extension     "
    echo "           fasta, fa, or fna.                               "
    echo "     -t    (Optional) 16S rRNA DNA reference FASTA file to  "
    echo "           use for training the RDP classifier and          "
    echo "           classifying extracted 16S sequences. Should be   "
    echo "           unzipped FASTA file with extension fasta, fa, or "
    echo "           fna. Should be of the following format:          "
    echo "                                                            "
    echo "     >Seq_ID Kingdom;Phylum;Class;Order;Family;Genus;Species"
    echo "     ACTGAAACTGCCGTTCAAAGCTTCGCGCGCTTTCCCGGGCGCGATATACGCGCGC"
    echo "     AAACTGGGGTCGCGAA                                       "
    echo "                                                            "
    echo "     -l    (Optional) Location of RDP classifier source dir "
    echo "           (rdp_classifier_x.xx). Safest to give full path  "
    echo "           to the dir, but a partial path should also work. "
    echo "     -s    (Optional) Skip certain steps in the pipeline if "
    echo "           need be. Provide a comma separated list of steps "
    echo "           that you wish to skip in the pipeline. List may  "
    echo "           have the values: fastqc, bbduk, kneaddata,       "
    echo "           humann.                                          "
    echo " "
    exit 0
    ;;
    i) SEQ_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    o) OUT_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    p) PROG_LOAD="$OPTARG"
    ;;
    r) HOST_REF=$(echo $OPTARG | sed 's#/$##')
    ;;
    c) CHOCO=$(echo $OPTARG | sed 's#/$##')
    ;;
    u) UNIREF=$(echo $OPTARG | sed 's#/$##')
    ;;
    m) MARKERS="$OPTARG"
    ;;
    f) FAIL_EMAIL="$OPTARG"
    ;;
    a) TAX_16S=1
    ;;
    e) EXTRACT_16S="$OPTARG"
    ;;
    t) TRAIN_16S="$OPTARG"
    ;;
    l) CLASS_LOC=$(echo $OPTARG | sed 's#/$##')
    ;;
    s) SKIP="$OPTARG"
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
if [[ -z "$OUT_DIR" ]]; then
  echo "ERROR: Argument -o is required, please supply an output directory"
  exit 1
fi
if [[ ! -d "$OUT_DIR" ]]; then
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

# -c
if [[ -z "$CHOCO" ]]; then
  echo "ERROR: Argument -c is required, please supply path to ChocoPhlAn database directory"
  exit 1
fi
if [[ ! -d "$CHOCO" ]]; then
  echo "ERROR: Argument -c should be a directory, please supply path to ChocoPhlAn database directory"
  exit 1
fi

# -u
if [[ -z "$UNIREF" ]]; then
  echo "ERROR: Argument -u is required, please supply path to UniRef database directory"
  exit 1
fi
if [[ ! -d "$UNIREF" ]]; then
  echo "ERROR: Argument -u should be a directory, please supply path to UniRef database directory"
  exit 1
fi

# -m
if [[ -z "$MARKERS" ]]; then
  echo "ERROR: Argument -m is required, please supply path to the clade marker info file mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2"
  exit 1
fi
if [[ -d "$MARKERS" ]]; then
  echo "ERROR: Argument -m should be the path to a single file, not a directory, please supply path to the clade marker info file mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2"
  exit 1
fi
if echo $MARKERS | grep -q -v "mpa_v30_CHOCOPhlAn_201901_marker_info\.txt\.bz2"; then
  echo "ERROR: path given to -m does not contain the file name mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2, please supply the clade marker info file mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2 to this argument"
  exit 1
fi

# -f
if [[ -z "$FAIL_EMAIL" ]]; then
  echo "ERROR: Argument -f is required, please supply an email that can be notified upon failure of any jobs ran during the pipeline"
  exit 1
elif echo $FAIL_EMAIL | grep -q -v '@'; then
  echo "ERROR: Argument -f requires a valid email, please give an email in the form of xxxx@xxxx.xxx"
fi

# -a
if [[ ! -z "$TAX_16S" ]]; then
  if [[ -z "$EXTRACT_16S" ]]; then
    echo "ERROR: when performing additional 16S-based taxonomic profiling, parameter -e is required, please supply a fasta reference file to use for extracting 16S sequences"
    exit 1
  elif [[ -z "$TRAIN_16S" ]]; then
    echo "ERROR: when performing additional 16S-based taxonomic profiling, parameter -t is required, please supply a fasta reference file to use for training and classifying with RDP classifier"
    exit 1
  elif [[ -z "$CLASS_LOC" ]]; then
    echo "ERROR: when performing additional 16S-based taxonomic profiling, parameter -l is required, please supply the location of the classifier.jar file for running the RDP classifier"
    exit 1
  else
    :
  fi
fi

# -e
if [[ ! -z "$EXTRACT_16S" ]]; then
  if [[ -z "$TAX_16S" ]]; then
    echo "WARNING: fasta reference file supplied for extracting 16S sequences, but flag -a is missing, therefore, this parameter will be ignored"
  elif echo $EXTRACT_16S | grep -q ".fasta"; then
    :
  elif echo $EXTRACT_16S | grep -q ".fa"; then
    :
  elif echo $EXTRACT_16S | grep -q ".fna"; then
    :
  else
    echo "ERROR: expecting reference file for extracting 16S sequences to be in fasta format with extension .fasta, .fa, or .fna, but none found"
    exit 1
  fi
fi

# -t
if [[ ! -z "$TRAIN_16S" ]]; then
  if [[ -z "$TAX_16S" ]]; then
    echo "WARNING: fasta reference file supplied for classifying 16S sequences, but flag -a is missing, therefore, this parameter will be ignored"
  elif echo $TRAIN_16S | grep -q ".fasta"; then
    :
  elif echo $TRAIN_16S | grep -q ".fa"; then
    :
  elif echo $TRAIN_16S | grep -q ".fna"; then
    :
  else
    echo "ERROR: expecting reference file for classifying 16S sequences to be in fasta format with extension .fasta, .fa, or .fna, but none found"
    exit 1
  fi
fi

# -l
if [[ ! -z "$CLASS_LOC" ]]; then
  if [[ -z "$TAX_16S" ]]; then
    echo "WARNING: classifier.jar location supplied for classifying 16S sequences, but flag -a is missing, therefore, this parameter will be ignored"
  elif [[ ! -d "$CLASS_LOC" ]]; then
    echo "ERROR: Argument -l should be a directory that contains the source files for RDP classifier, please supply the directory that contains these files"
    exit 1
  else
    :
  fi
fi

# -s
if [[ ! -z "$SKIP" ]]; then
  if echo $SKIP | grep -q "fastqc"; then
    :
  elif echo $SKIP | grep -q "bbduk"; then
    :
  elif echo $SKIP | grep -q "kneaddata"; then
    :
  elif echo $SKIP | grep -q "humann"; then
    :
  elif echo $SKIP | grep -q "graphlan"; then
    :
  else
    echo "ERROR: Invalid argument given to -s, please specify one or more of: fastqc_initial,bbmerge,bbduk,kneaddata,fastqc_postqc,humann,graphlan"
    exit 1
  fi
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

############# INITIAL FASTQC REPORT #############
if echo $SKIP | grep -q "fastqc"; then
  echo " "
  echo "*** Skipping running of initial FastQC report generation on input fastq files ***"
  echo " "
else
  SECONDS=0
  echo " "
  echo "*** Running FastQC on input fastq files ***"
  echo " "
  
  #Create directory for output
  if [ -d "${RESULTS_DIR}/1.FastQC_Initial_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output
  fi
  
  ##### Run FastQC #####
  #Create shell script for running program
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$1 | awk -F '/' '{print \$NF}' | awk -F '_' '{print \$1,\$2,\$3}' OFS='_')" >> bash_script.sh
  echo "fastqc \$1 \$2 -o ${RESULTS_DIR}/1.FastQC_Initial_Reports \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/1.FastQC_Initial_Reports/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  chmod +x bash_script.sh
  
  #For every sequence file submit job and grab job IDs
  touch job_ids.txt
  for file in ${SEQ_DIR}/*R1_001.${SEQ_EXT}; do
    FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '_' '{print $1,$2,$3}' OFS='_')
    
    sbatch --partition=express \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut/${FILE_NAME}.err \
    --output=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output/${FILE_NAME}.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script.sh ${SEQ_DIR}/${FILE_NAME}_R1_001.${SEQ_EXT} ${SEQ_DIR}/${FILE_NAME}_R2_001.${SEQ_EXT} | \
    awk '{print $4}' >> job_ids.txt
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
  echo "FastQC reports complete"
  echo " "
  
  ##### Get mean Q score per file #####
  echo "Grabbing average quality score per sequence file..."
  echo " "
  echo 'Filename Mean SD' > ${RESULTS_DIR}/1.FastQC_Initial_Reports/mean_Q_per_file.txt
  for file in ${RESULTS_DIR}/1.FastQC_Initial_Reports/*fastqc.zip
  do
    DIR_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '.zip' '{print $1}')
    unzip -qq $file
    sed -n '/>>Per sequence quality scores/,/>>END_MODULE/p' ${DIR_NAME}/fastqc_data.txt | \
    sed '1,2d' | sed '$d' > q_vals.txt
    echo "q_vals <- read.table('q_vals.txt')" > Rfunc.R
    echo "cat(paste(round(mean(rep(q_vals[,1],q_vals[,2])),2), round(sd(rep(q_vals[,1],q_vals[,2])),2)), '\n')" >> Rfunc.R
    echo $(echo $file | awk -F '/' '{print $NF}' | awk -F '_fastqc.zip' '{print $1}') $(Rscript --vanilla Rfunc.R) >> ${RESULTS_DIR}/1.FastQC_Initial_Reports/mean_Q_per_file.txt
    rm -r ${DIR_NAME}
    rm Rfunc.R
    rm q_vals.txt
  done
  echo "Done"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

################# QC WITH BBDUK #################
if echo $SKIP | grep -q "bbduk"; then
  echo " "
  echo "*** Skipping running of BBDuk for adapter/quality trimming and filtering of input fastq files ***"
  echo " "
else
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
  echo "FILE_NAME=\$(echo \$1 | awk -F '/' '{print \$NF}' | awk -F '_' '{print \$1,\$2,\$3}' OFS='_')" >> bash_script.sh
  echo "FILE1=\$(echo \$1 | awk -F '/' '{print \$NF}')" >> bash_script.sh
  echo "FILE2=\$(echo \$2 | awk -F '/' '{print \$NF}')" >> bash_script.sh
  echo "bbduk.sh in=\$1 \\" >> bash_script.sh
  echo "in2=\$2 \\" >> bash_script.sh
  echo "out=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE1} \\" >> bash_script.sh
  echo "out2=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE2} \\" >> bash_script.sh
  echo "stats=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}_stats.txt \\" >> bash_script.sh
  echo "ftm=5 tpe tbo qtrim=rl trimq=25 minlen=50 ref=adapters,phix -Xmx64000m \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}.log 2>&1"  >> bash_script.sh
  chmod +x bash_script.sh
  
  #For every pair of paired-end sequence files submit job and grab job IDs
  touch job_ids.txt
  for file in ${SEQ_DIR}/*R1_001.${SEQ_EXT}; do
    FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '_' '{print $1,$2,$3}' OFS='_')
    
    sbatch --partition=express \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.ErrorOut/${FILE_NAME}.err \
    --output=${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.Output/${FILE_NAME}.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script.sh ${SEQ_DIR}/${FILE_NAME}_R1_001.${SEQ_EXT} ${SEQ_DIR}/${FILE_NAME}_R2_001.${SEQ_EXT} | \
    awk '{print $4}' >> job_ids.txt
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
fi
#################################################

##### REMOVAL OF HOST READS WITH KNEADDATA ######
if echo $SKIP | grep -q "kneaddata"; then
  echo " "
  echo "*** Skipping running of KneadData for removal of contaminant host reads ***"
  echo " "
else
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
  #Create shell script for running program
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$1 | awk -F '/' '{print \$NF}' | awk -F '_' '{print \$1,\$2,\$3}' OFS='_')" >> bash_script.sh
  echo "kneaddata --input \$1 \\" >> bash_script.sh
  echo "--input \$2 \\" >> bash_script.sh
  echo "--output ${RESULTS_DIR}/3.Decontaminated_Sequences \\" >> bash_script.sh
  echo "--output-prefix \$FILE_NAME \\" >> bash_script.sh
  echo "--log ${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}_kneaddata.log \\" >> bash_script.sh
  echo "--reference-db $HOST_REF \\" >> bash_script.sh
  echo "--bypass-trim \\" >> bash_script.sh
  echo "--threads 5 \\" >> bash_script.sh
  echo "--verbose \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  chmod +x bash_script.sh
  
  #For every quality controlled pair of sequence files submit job and grab job IDs
  touch job_ids.txt
  for file in ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*R1_001.${SEQ_EXT}; do
    FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '_' '{print $1,$2,$3}' OFS='_')
    
    sbatch --partition=express \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/3.Decontaminated_Sequences/0.ErrorOut/${FILE_NAME}.err \
    --output=${RESULTS_DIR}/3.Decontaminated_Sequences/0.Output/${FILE_NAME}.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=5 \
    --mem-per-cpu=32000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script.sh ${RESULTS_DIR}/2.Quality_Controlled_Sequences/${FILE_NAME}_R1_001.${SEQ_EXT} \
    ${RESULTS_DIR}/2.Quality_Controlled_Sequences/${FILE_NAME}_R2_001.${SEQ_EXT} | \
    awk '{print $4}' >> job_ids.txt
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
  rm bash_script.sh
  
  rm ${RESULTS_DIR}/3.Decontaminated_Sequences/*unmatched*fastq
  
  ##### Gzip output #####
  echo "Compressing KneadData output..."
  echo " "
  
  #Create shell script for running program
  echo '#!/bin/bash' > bash_script.sh
  echo "gzip \$1" >> bash_script.sh
  
  #For every fastq file submit job and grab job IDs
  touch job_ids.txt
  for file in ${RESULTS_DIR}/3.Decontaminated_Sequences/*.fastq; do
    FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '_' '{print $1,$2,$3}' OFS='_')
    
    sbatch --partition=express \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/3.Decontaminated_Sequences/0.ErrorOut/${FILE_NAME}_gzip.err \
    --output=${RESULTS_DIR}/3.Decontaminated_Sequences/0.Output/${FILE_NAME}_gzip.out \
    --time=2:00:00 \
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
  rm job_ids.txt
  rm bash_script.sh
  
  echo "Done"
  echo " "
  
  #Signal jobs have ended
  echo "Running of KneadData complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
  
  #Move extracted host sequences to their own directory
  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences/Extracted_Host_Sequences
  mv ${RESULTS_DIR}/3.Decontaminated_Sequences/*contam*fastq.gz ${RESULTS_DIR}/3.Decontaminated_Sequences/Extracted_Host_Sequences/
fi
#################################################

###### TAXONOMIC AND FUNCTIONAL PROFILING #######
if echo $SKIP | grep -q "humann"; then
  echo " "
  echo "*** Skipping running of HUMAnN workflow for performing taxonomic and functional profiling ***"
  echo " "
else
  SECONDS=0
  echo " "
  echo "*** Running HUMAnN workflow for performing taxonomic and functional profiling ***"
  echo " "
  
  #Create directory for profiling output
  if [ -d "${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling
	  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut
	  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output
  fi
  
  ##### Run HUMAnN workflow #####
  #Create shell script for running program
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$1 | awk -F '/' '{print \$NF}' | awk -F '_' '{print \$1,\$2,\$3}' OFS='_')" >> bash_script.sh
  echo "zcat \$1 \$2 > ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.temp.fastq" >> bash_script.sh
  echo "humann --input ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.temp.fastq \\" >> bash_script.sh
  echo "--output ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling \\" >> bash_script.sh
  echo "--output-basename \$FILE_NAME \\" >> bash_script.sh
  echo "--metaphlan-options '-t rel_ab --add_viruses' \\" >> bash_script.sh
  echo "--nucleotide-database $CHOCO \\" >> bash_script.sh
  echo "--protein-database $UNIREF \\" >> bash_script.sh
  echo "--prescreen-threshold 0.001 \\" >> bash_script.sh
  echo "--threads 5 \\" >> bash_script.sh
  echo "--verbose \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.temp.fastq" >> bash_script.sh
  chmod +x bash_script.sh
  
  #For every pair of decontaminated sequence files submit job and grab job IDs
  touch job_ids.txt
  for file in ${RESULTS_DIR}/3.Decontaminated_Sequences/*paired_1.fastq.gz; do
    FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '_' '{print $1,$2,$3}' OFS='_')
    
    sbatch --partition=medium \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/${FILE_NAME}_humann.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/${FILE_NAME}_humann.out \
    --time=24:00:00 \
    --ntasks=1 \
    --cpus-per-task=5 \
    --mem-per-cpu=32000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script.sh ${RESULTS_DIR}/3.Decontaminated_Sequences/${FILE_NAME}_paired_1.fastq.gz \
    ${RESULTS_DIR}/3.Decontaminated_Sequences/${FILE_NAME}_paired_2.fastq.gz | \
    awk '{print $4}' >> job_ids.txt
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
  rm bash_script.sh
  
  ##### Run MetaPhlAn again to get normalized counts #####
  #Create shell script for running program
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$1 | awk -F '/' '{print \$NF}' | awk -F '_' '{print \$1,\$2,\$3}' OFS='_')" >> bash_script.sh
  echo "metaphlan --input_type bowtie2out \\" >> bash_script.sh
  echo "--add_viruses \\" >> bash_script.sh
  echo "-t marker_ab_table \\" >> bash_script.sh
  echo "\${1}/\${FILE_NAME}_metaphlan_bowtie2.txt \\" >> bash_script.sh
  echo "\${1}/\${FILE_NAME}_metaphlan_norm_abun_table.tsv \\" >> bash_script.sh
  echo ">> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  echo "join -j 1 -o 1.3,1.2,2.2 \\" >> bash_script.sh
  echo "<(sort -k1,1 <(paste <(bzcat $MARKERS | awk -F\"\t\" '{print \$1}') \\" >> bash_script.sh
  echo "                     <(bzcat $MARKERS | awk -F\"__\" '{print \$1}') \\" >> bash_script.sh
  echo "                     <(bzcat $MARKERS | awk -F\"('taxon': '|'})\" '{print \$2}'))) \\" >> bash_script.sh
  echo "<(sort -k1,1 \${1}/\${FILE_NAME}_metaphlan_norm_abun_table.tsv) | \\" >> bash_script.sh
  echo "sed '1s/^/#clade_name NCBI_tax_id normalized_abundance\n/' | sed 's/ /\t/g' | \\" >> bash_script.sh
  echo "awk 'FNR==1{print;next}{val=\$3;\$3=\"~\";a[\$0]+=val}"'!'"b[\$0]++{c[++count]=\$0}END{for(i=1;i<=count;i++){sub(\"~\",a[c[i]],c[i]);print c[i]}}' OFS='\t' | \\" >> bash_script.sh
  echo "cat <(grep '^#' \${1}/\${FILE_NAME}_metaphlan_norm_abun_table.tsv) - > \${FILE_NAME}_temp.txt" >> bash_script.sh
  echo "mv \${FILE_NAME}_temp.txt \${1}/\${FILE_NAME}_metaphlan_norm_abun_table.tsv" >> bash_script.sh
  chmod +x bash_script.sh
  
  #For every metaphlan output submit job and grab job IDs
  touch job_ids.txt
  for dir in ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*humann_temp; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_' '{print $1,$2,$3}' OFS='_')
    
    sbatch --partition=express \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/${FILE_NAME}_metaphlan.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/${FILE_NAME}_metaphlan.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script.sh $dir | awk '{print $4}' >> job_ids.txt
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
  rm bash_script.sh
  
  ##### Clean up #####
  echo "Cleaning things up... merging per sample tables and removing unneeded files"
  echo " "
  for dir in ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*humann_temp; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_' '{print $1,$2,$3}' OFS='_')
    
    mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_Profiles
    
    mv ${dir}/${FILE_NAME}_metaphlan_bugs_list.tsv \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_Profiles/${FILE_NAME}_metaphlan_rel_abun_table.tsv
    
    mv ${dir}/${FILE_NAME}_metaphlan_norm_abun_table.tsv \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_Profiles/${FILE_NAME}_metaphlan_norm_abun_table.tsv
    
    mv ${dir}/${FILE_NAME}.log \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_Profiles/${FILE_NAME}_humann.log
    
    rm -r ${dir}
    
    mv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}*.* \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_Profiles/
  done
  
  ##### Merge tables #####
  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables
  
  for dir in ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*_Profiles; do
    cp ${dir}/*metaphlan_rel_abun_table.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/
    cp ${dir}/*metaphlan_norm_abun_table.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/
    cp ${dir}/*genefamilies.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/
    cp ${dir}/*pathabundance.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/
    cp ${dir}/*pathcoverage.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/
  done
  
  touch job_ids.txt
  
  echo '#!/bin/bash' > bash_script_1.sh
  echo "$PROG_LOAD" >> bash_script_1.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_rel_abun_table.tsv \\" >> bash_script_1.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_rel_abun_table.tsv" >> bash_script_1.sh
  echo "sed -i '2s/_metaphlan_rel_abun_table//g' ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_rel_abun_table.tsv" >> bash_script_1.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_rel_abun_table.tsv" >> bash_script_1.sh
  
  sbatch --partition=express \
    --job-name=metaphlan_rel_abun_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/metaphlan_rel_abun_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/metaphlan_rel_abun_merge.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script_1.sh | awk '{print $4}' >> job_ids.txt
  
  echo '#!/bin/bash' > bash_script_2.sh
  echo "$PROG_LOAD" >> bash_script_2.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_norm_abun_table.tsv \\" >> bash_script_2.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_norm_abun_table.tsv" >> bash_script_2.sh
  echo "sed -i '2s/_metaphlan_norm_abun_table//g' ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_norm_abun_table.tsv" >> bash_script_2.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_norm_abun_table.tsv" >> bash_script_2.sh
  
  sbatch --partition=express \
    --job-name=metaphlan_norm_abun_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/metaphlan_norm_abun_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/metaphlan_norm_abun_merge.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script_2.sh | awk '{print $4}' >> job_ids.txt
    
  echo '#!/bin/bash' > bash_script_3.sh
  echo "$PROG_LOAD" >> bash_script_3.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/ \\" >> bash_script_3.sh
  echo "--output ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/humann_genefamilies.tsv \\" >> bash_script_3.sh
  echo "--file_name genefamilies.tsv" >> bash_script_3.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*genefamilies.tsv" >> bash_script_3.sh
  
  sbatch --partition=express \
    --job-name=genefamilies_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/genefamilies_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/genefamilies_merge.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script_3.sh | awk '{print $4}' >> job_ids.txt
  
  echo '#!/bin/bash' > bash_script_4.sh
  echo "$PROG_LOAD" >> bash_script_4.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/ \\" >> bash_script_4.sh
  echo "--output ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/humann_pathabundance.tsv \\" >> bash_script_4.sh
  echo "--file_name pathabundance.tsv" >> bash_script_4.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*pathabundance.tsv" >> bash_script_4.sh
  
  sbatch --partition=express \
    --job-name=pathabundance_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/pathabundance_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/pathabundance_merge.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script_4.sh | awk '{print $4}' >> job_ids.txt

  echo '#!/bin/bash' > bash_script_5.sh
  echo "$PROG_LOAD" >> bash_script_5.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/ \\" >> bash_script_5.sh
  echo "--output ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/humann_pathcoverage.tsv \\" >> bash_script_5.sh
  echo "--file_name pathcoverage.tsv" >> bash_script_5.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*pathcoverage.tsv" >> bash_script_5.sh
  
  sbatch --partition=express \
    --job-name=pathcoverage_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/pathcoverage_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/pathcoverage_merge.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script_5.sh | awk '{print $4}' >> job_ids.txt
  
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
  rm bash_script_*.sh
  
  #Signal workflow has completed
  echo "Running of HUMAnN workflow complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi

##### 16S based taxonomic profiling #####
if [[ -z "$TAX_16S" ]]; then
  :
else
  SECONDS=0
  echo " "
  echo "*** Running 16S based taxonomic profiling ***"
  echo " "

  #Create directory for profiling output
  if [ -d "${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling
	  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut
	  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output
  fi
  if [ -d "${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling" ]
  then
         :
  else
         mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling
  fi
  
  ##### Run SortMeRNA #####
  echo "Running SortMeRNA for extracting 16S rRNA gene sequences specified in reference file..."
  echo " "
  
  #Create shell script for running program
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$1 | awk -F '/' '{print \$NF}' | awk -F '_' '{print \$1,\$2,\$3}' OFS='_')" >> bash_script.sh
  echo "if [[ -d ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${FILE_NAME}_16S ]]; then" >> bash_script.sh
  echo "  rm -r ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${FILE_NAME}_16S" >> bash_script.sh
  echo "  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${FILE_NAME}_16S" >> bash_script.sh
  echo "else" >> bash_script.sh
  echo "  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${FILE_NAME}_16S" >> bash_script.sh
  echo "fi" >> bash_script.sh
  echo "zcat \$1 \$2 > ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${FILE_NAME}.temp.fastq" >> bash_script.sh
  echo "sortmerna --ref $EXTRACT_16S \\" >> bash_script.sh
  echo "--reads ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${FILE_NAME}.temp.fastq \\" >> bash_script.sh
  echo "--workdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${FILE_NAME}_16S/ \\" >> bash_script.sh
  echo "--fastx True \\" >> bash_script.sh
  echo "--threads 5 \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${FILE_NAME}_16S/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${FILE_NAME}.temp.fastq" >> bash_script.sh
  chmod +x bash_script.sh
  
  #For every pair of decontaminasted sequence files submit job and grab job IDs
  touch job_ids.txt
  for file in ${RESULTS_DIR}/3.Decontaminated_Sequences/*paired_1.fastq.gz; do
    FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '_' '{print $1,$2,$3}' OFS='_')
    
    sbatch --partition=short \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/${FILE_NAME}_sortmerna.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/${FILE_NAME}_sortmerna.out \
    --time=12:00:00 \
    --ntasks=1 \
    --cpus-per-task=5 \
    --mem-per-cpu=32000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script.sh ${RESULTS_DIR}/3.Decontaminated_Sequences/${FILE_NAME}_paired_1.fastq.gz \
    ${RESULTS_DIR}/3.Decontaminated_Sequences/${FILE_NAME}_paired_2.fastq.gz | \
    awk '{print $4}' >> job_ids.txt
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
  rm bash_script.sh
  
  #Signal jobs have ended
  echo "Done"
  echo " "
  
  ##### Training RDP classifier ####
  echo "Training the RDP classifier with provided reference FASTA file..."
  echo " "
  
  #Create script needed for creating files for training classifier
  #NOTE: following python scripts are not original, but provided by RDP staff
  #and re-implemented here to avoid having a dependency of obtaining the
  #scripts before running pipeline
  echo "import sys, string" > lineage2taxTrain.py
  echo "if not len(sys.argv) == 2:" >> lineage2taxTrain.py
  echo '        print "lineage2taxTrain.py taxonomyFile"' >> lineage2taxTrain.py
  echo "        sys.exit()" >> lineage2taxTrain.py
  echo " " >> lineage2taxTrain.py
  echo "f = open(sys.argv[1], 'r').readlines()" >> lineage2taxTrain.py
  echo "header = f[0]" >> lineage2taxTrain.py
  echo "header = f[0].strip().split('\t')[1:]#header: list of ranks" >> lineage2taxTrain.py
  echo "hash = {}#taxon name-id map" >> lineage2taxTrain.py
  echo "ranks = {}#column number-rank map" >> lineage2taxTrain.py
  echo "lineages = []#list of unique lineages" >> lineage2taxTrain.py
  echo " " >> lineage2taxTrain.py
  echo 'hash = {"Root":0}#initiate root rank taxon id map' >> lineage2taxTrain.py
  echo "for i in range(len(header)):" >> lineage2taxTrain.py
  echo "        name = header[i]" >> lineage2taxTrain.py
  echo "        ranks[i] = name" >> lineage2taxTrain.py
  echo "root = ['0', 'Root', '-1', '0', 'rootrank']#root rank info" >> lineage2taxTrain.py
  echo "print string.join(root, '*')" >> lineage2taxTrain.py
  echo "ID = 0 #taxon id" >> lineage2taxTrain.py
  echo "for line in f[1:]:" >> lineage2taxTrain.py
  echo "        cols = line.strip().split('\t')[1:]" >> lineage2taxTrain.py
  echo "        for i in range(len(cols)):#iterate each column" >> lineage2taxTrain.py
  echo "                name = []" >> lineage2taxTrain.py
  echo "                for node in cols[:i + 1]:" >> lineage2taxTrain.py
  echo "                        node = node.strip()" >> lineage2taxTrain.py
  echo "                        if not node in ('-', ''):" >> lineage2taxTrain.py
  echo "                                name.append(node)" >> lineage2taxTrain.py
  echo "                pName = string.join(name[:-1], ';')" >> lineage2taxTrain.py
  echo "                if not name in lineages:" >> lineage2taxTrain.py
  echo "                        lineages.append(name)" >> lineage2taxTrain.py
  echo "                depth = len(name)" >> lineage2taxTrain.py
  echo "                name = string.join(name, ';')" >> lineage2taxTrain.py
  echo "                if name in hash.keys():#already seen this lineage" >> lineage2taxTrain.py
  echo "                        continue" >> lineage2taxTrain.py
  echo "                try:" >> lineage2taxTrain.py
  echo "                        rank = ranks[i]" >> lineage2taxTrain.py
  echo "                except KeyError:" >> lineage2taxTrain.py
  echo "                        print cols" >> lineage2taxTrain.py
  echo "                        sys.exit()" >> lineage2taxTrain.py
  echo "                if i == 0:" >> lineage2taxTrain.py
  echo "                        pName = 'Root'" >> lineage2taxTrain.py
  echo "                pID = hash[pName]#parent taxid" >> lineage2taxTrain.py
  echo "                ID += 1" >> lineage2taxTrain.py
  echo "                hash[name] = ID #add name-id to the map" >> lineage2taxTrain.py
  echo "                out = ['%s'%ID, name.split(';')[-1], '%s'%pID, '%s'%depth, rank]" >> lineage2taxTrain.py
  echo "                print string.join(out, '*')" >> lineage2taxTrain.py
  chmod +x lineage2taxTrain.py
  
  echo "import sys, string" > addFullLineage.py
  echo 'if len(sys.argv) != 3:' >> addFullLineage.py
  echo "        print 'addFullLineage.py taxonomyFile fastaFile'" >> addFullLineage.py
  echo "        sys.exit()" >> addFullLineage.py
  echo "f1 = open(sys.argv[1], 'r').readlines()" >> addFullLineage.py
  echo "hash = {} #lineage map" >> addFullLineage.py
  echo "for line in f1[1:]:" >> addFullLineage.py
  echo "        line = line.strip()" >> addFullLineage.py
  echo "        cols = line.strip().split('\t')" >> addFullLineage.py
  echo "        lineage = ['Root']" >> addFullLineage.py
  echo "        for node in cols[1:]:" >> addFullLineage.py
  echo "                node = node.strip()" >> addFullLineage.py
  echo "                if not (node == '-' or node == ''):" >> addFullLineage.py
  echo "                        lineage.append(node)" >> addFullLineage.py
  echo "        ID = cols[0]" >> addFullLineage.py
  echo "        lineage = string.join(lineage, ';').strip()" >> addFullLineage.py
  echo "        hash[ID] = lineage" >> addFullLineage.py
  echo "f2 = open(sys.argv[2], 'r').readlines()" >> addFullLineage.py
  echo "for line in f2:" >> addFullLineage.py
  echo "        line = line.strip()" >> addFullLineage.py
  echo "        if line == '':" >> addFullLineage.py
  echo "                continue" >> addFullLineage.py
  echo "        if line[0] == '>':" >> addFullLineage.py
  echo "                ID = line.strip().split()[0].replace('>', '')" >> addFullLineage.py
  echo "                try:" >> addFullLineage.py
  echo "                        lineage = hash[ID]" >> addFullLineage.py
  echo "                except KeyError:" >> addFullLineage.py
  echo "                        print ID, 'not in taxonomy file'" >> addFullLineage.py
  echo "                        sys.exit()" >> addFullLineage.py
  echo "                print '>' + ID + '\t' + lineage" >> addFullLineage.py
  echo "        else:" >> addFullLineage.py
  echo "                        print line.strip()" >> addFullLineage.py
  chmod +x addFullLineage.py
  
  #Create shell script for running RDP training
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "echo 'Removing duplicate sequences in provided reference file...'" >> bash_script.sh
  echo "echo ' '" >> bash_script.sh
  echo "awk 'BEGIN {i = 1;} { if (\$1 ~ /^>/) { tmp = h[i]; h[i] = \$0; } else if ("'!'"a[\$1]) { s[i] = \$0; a[\$1] = \"1\"; i++; } else { h[i] = tmp; } } END { for (j = 1; j < i; j++) { print h[j]; print s[j]; } }' $TRAIN_16S > ref.fasta" >> bash_script.sh
  echo "echo 'Sequence removal complete'" >> bash_script.sh
  echo "echo ' '" >> bash_script.sh
  echo " " >> bash_script.sh
  echo "echo 'Creating taxonomy file for training classifier...'" >> bash_script.sh
  echo "echo ' '" >> bash_script.sh
  echo "python2 lineage2taxTrain.py \\" >> bash_script.sh
  echo "<(cat ref.fasta | grep '>' | sed 's/>//' | sed 's/ /\t/' | sed 's/;/\t/g' | \\" >> bash_script.sh
  echo "sed '1s/^/Seq_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' | \\" >> bash_script.sh
  echo "awk '\$2 == \"\"{\$2 = \"Unknown_kingdom\"; \$3 = \"Unknown_phylum\"; \$4 = \"Unknown_class\"; \$5 = \"Unknown_order\"; \$6 = \"Unknown_family\"; \$7 = \"Unknown_genus\"; \$8 = \"Unknown_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$3 == \"\"{\$3 = \$2\"_phylum\"; \$4 = \$2\"_class\"; \$5 = \$2\"_order\"; \$6 = \$2\"_family\"; \$7 = \$2\"_genus\"; \$8 = \$2\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$4 == \"\"{\$4 = \$3\"_class\"; \$5 = \$3\"_order\"; \$6 = \$3\"_family\"; \$7 = \$3\"_genus\"; \$8 = \$3\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$5 == \"\"{\$5 = \$4\"_order\"; \$6 = \$4\"_family\"; \$7 = \$4\"_genus\"; \$8 = \$4\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$6 == \"\"{\$6 = \$5\"_family\"; \$7 = \$5\"_genus\"; \$8 = \$5\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$7 == \"\"{\$7 = \$6\"_genus\"; \$8 = \$6\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$8 == \"\"{\$8 = \$7\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' OFS=\"\t\") \\" >> bash_script.sh
  echo "> train.tax" >> bash_script.sh
  echo "echo 'Done'" >> bash_script.sh
  echo "echo ' '" >> bash_script.sh
  echo "echo 'Finalizing FASTA for training...'" >> bash_script.sh
  echo "echo ' '" >> bash_script.sh
  echo "python2 addFullLineage.py \\" >> bash_script.sh
  echo "<(cat ref.fasta | grep '>' | sed 's/>//' | sed 's/ /\t/' | sed 's/;/\t/g' | \\" >> bash_script.sh
  echo "sed '1s/^/Seq_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' | \\" >> bash_script.sh
  echo "awk '\$2 == \"\"{\$2 = \"Unknown_kingdom\"; \$3 = \"Unknown_phylum\"; \$4 = \"Unknown_class\"; \$5 = \"Unknown_order\"; \$6 = \"Unknown_family\"; \$7 = \"Unknown_genus\"; \$8 = \"Unknown_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$3 == \"\"{\$3 = \$2\"_phylum\"; \$4 = \$2\"_class\"; \$5 = \$2\"_order\"; \$6 = \$2\"_family\"; \$7 = \$2\"_genus\"; \$8 = \$2\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$4 == \"\"{\$4 = \$3\"_class\"; \$5 = \$3\"_order\"; \$6 = \$3\"_family\"; \$7 = \$3\"_genus\"; \$8 = \$3\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$5 == \"\"{\$5 = \$4\"_order\"; \$6 = \$4\"_family\"; \$7 = \$4\"_genus\"; \$8 = \$4\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$6 == \"\"{\$6 = \$5\"_family\"; \$7 = \$5\"_genus\"; \$8 = \$5\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$7 == \"\"{\$7 = \$6\"_genus\"; \$8 = \$6\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '\$8 == \"\"{\$8 = \$7\"_species\"};1' OFS=\"\t\" | \\" >> bash_script.sh
  echo "awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' OFS=\"\t\") \\" >> bash_script.sh
  echo "ref.fasta > train.fasta" >> bash_script.sh
  echo "echo 'Done'" >> bash_script.sh
  echo "echo ' '" >> bash_script.sh
  echo "echo 'Training classifier using supplied reference FASTA file...'" >> bash_script.sh
  echo "echo ' '" >> bash_script.sh
  echo "if [ -d \"${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/RDP_training_files\" ]" >> bash_script.sh
  echo "then" >> bash_script.sh
  echo "  rm -r ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/RDP_training_files" >> bash_script.sh
  echo "else" >> bash_script.sh
  echo "  :" >> bash_script.sh
  echo "fi" >> bash_script.sh
  echo "java -Xmx64g -jar ${CLASS_LOC}/dist/classifier.jar train \\" >> bash_script.sh
  echo "-o ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/RDP_training_files \\" >> bash_script.sh
  echo "-s train.fasta \\" >> bash_script.sh
  echo "-t train.tax" >> bash_script.sh
  echo "cp ${CLASS_LOC}/src/data/classifier/16srrna/rRNAClassifier.properties \\" >> bash_script.sh
  echo "${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/RDP_training_files/" >> bash_script.sh
  echo "echo 'Classifier training complete'" >> bash_script.sh
  echo "echo ' '" >> bash_script.sh
  chmod +x bash_script.sh
  
  #Submit job
  touch job_ids.txt
  sbatch --partition=short \
  --job-name=RDP_train \
  --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/RDP_train.err \
  --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/RDP_train.out \
  --time=12:00:00 \
  --ntasks=1 \
  --cpus-per-task=1 \
  --mem-per-cpu=64000 \
  --mail-type=FAIL \
  --mail-user=${FAIL_EMAIL} \
  ./bash_script.sh | awk '{print $4}' >> job_ids.txt
  
  #Hold script here until job completed
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
  rm lineage2taxTrain.py
  rm addFullLineage.py
  rm ref.fasta
  rm train.fasta
  rm train.tax
  rm bash_script.sh
  rm job_ids.txt
  
  #Signal jobs have ended
  echo "Done"
  echo " "
  
  ##### Classifying 16S sequences ####
  echo "Classifying extracted 16S sequences..."
  echo " "
  
  #Create shell script for running classification
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "java -Xmx64g -jar ${CLASS_LOC}/dist/classifier.jar classify \\" >> bash_script.sh
  echo "-c 0.8 \\" >> bash_script.sh
  echo "-f allrank \\" >> bash_script.sh
  echo "-t ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/RDP_training_files/rRNAClassifier.properties \\" >> bash_script.sh
  echo "-o \${1}/out/aligned_classified.txt \\" >> bash_script.sh
  echo "-h \${1}/out/aligned_classified_counts.txt \\" >> bash_script.sh
  echo "\${1}/out/aligned.fastq" >> bash_script.sh
  chmod +x bash_script.sh
  
  #Submit jobs for each directory of extracted sequences
  touch job_ids.txt
  for dir in ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/*_16S; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_' '{print $1,$2,$3}' OFS='_')
    
    sbatch --partition=express \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/${FILE_NAME}_classify.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/${FILE_NAME}_classify.out \
    --time=2:00:00 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu=64000 \
    --mail-type=FAIL \
    --mail-user=${FAIL_EMAIL} \
    ./bash_script.sh ${dir} | awk '{print $4}' >> job_ids.txt
  done
  
  #Hold script here until job completed
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
  echo "Done"
  echo " "
  
  ##### Create abundance tables from classification #####
  echo "Constructing abundance tables from classifications..."
  echo " "
  
  #Create shell script for running construction of abundance tables
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "for dir in ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/*_16S; do" >> bash_script.sh
  echo "  ID=\$(echo \$dir | awk -F '/' '{print \$NF}' | awk -F '_' '{print \$1}' OFS='_')" >> bash_script.sh
  #echo "  cat \\" >> bash_script.sh
  #echo "  <(awk '\$4 == \"Species\"{print \$2, \$5}' \${dir}/out/aligned_classified_counts.txt | \\" >> bash_script.sh
  #echo "  sed 's/Root;rootrank;//' | sed 's/;Kingdom;/;/' | sed 's/;Phylum;/;/' | sed 's/;Class;/;/' | \\" >> bash_script.sh
  #echo "  sed 's/;Order;/;/' | sed 's/;Family;/;/' | sed 's/;Genus;/;/' | sed 's/;Species;//' | \\" >> bash_script.sh
  #echo "  sed \"1s/^/Taxa \$ID\n/\") \\" >> bash_script.sh
  #echo "  <(awk '/unclassified/{print \$2, \$4}' \${dir}/out/aligned_classified_counts.txt | \\" >> bash_script.sh
  #echo "  sed 's/Root;rootrank;//' | sed 's/;Kingdom;/;/' | sed 's/;Phylum;/;/' | sed 's/;Class;/;/' | \\" >> bash_script.sh
  #echo "  sed 's/;Order;/;/' | sed 's/;Family;/;/' | sed 's/;Genus;/;/' | sed 's/;;//') > \${ID}_temp.txt" >> bash_script.sh
  #echo "  sed 's/;/ /g' \${ID}_temp.txt | \\" >> bash_script.sh
  #echo "  awk '\$1 ~ /unclassified/{\$8 = \$2; \$1 = \"Unclassified_kingdom\"; \$2 = \"Unclassified_phylum\"; \$3 = \"Unclassified_class\"; \$4 = \"Unclassified_order\"; \$5 = \"Unclassified_family\"; \$6 = \"Unclassified_genus\"; \$7 = \"Unclassified_species\"};1' | \\" >> bash_script.sh
  #echo "  awk '\$2 ~ /unclassified/{\$8 = \$3; \$2 = \$1\"_phylum\"; \$3 = \$1\"_class\"; \$4 = \$1\"_order\"; \$5 = \$1\"_family\"; \$6 = \$1\"_genus\"; \$7 = \$1\"_species\"};1' | \\" >> bash_script.sh
  #echo "  awk '\$3 ~ /unclassified/{\$8 = \$4; \$3 = \$2\"_class\"; \$4 = \$2\"_order\"; \$5 = \$2\"_family\"; \$6 = \$2\"_genus\"; \$7 = \$2\"_species\"};1' | \\" >> bash_script.sh
  #echo "  awk '\$4 ~ /unclassified/{\$8 = \$5; \$4 = \$3\"_order\"; \$5 = \$3\"_family\"; \$6 = \$3\"_genus\"; \$7 = \$3\"_species\"};1' | \\" >> bash_script.sh
  #echo "  awk '\$5 ~ /unclassified/{\$8 = \$6; \$5 = \$4\"_family\"; \$6 = \$4\"_genus\"; \$7 = \$4\"_species\"};1' | \\" >> bash_script.sh
  #echo "  awk '\$6 ~ /unclassified/{\$8 = \$7; \$6 = \$5\"_genus\"; \$7 = \$5\"_species\"};1' | \\" >> bash_script.sh
  #echo "  awk '\$7 ~ /unclassified/{\$8 = \$8; \$7 = \$6\"_species\"};1' | sed 's/ /;/g' | sed 's/\(.*\);/\1 /' | \\" >> bash_script.sh
  #echo "  awk 'FNR==1{print;next}{val=\$2;\$2=\"~\";a[\$0]+=val}"'!'"b[\$0]++{c[++count]=\$0}END{for(i=1;i<=count;i++){sub(\"~\",a[c[i]],c[i]);print c[i]}}' \\" >> bash_script.sh
  #echo "  > ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${ID}_species.txt" >> bash_script.sh
  echo "  cat \\" >> bash_script.sh
  echo "  <(awk '\$4 == \"Genus\"{print \$2, \$5}' \${dir}/out/aligned_classified_counts.txt | \\" >> bash_script.sh
  echo "  sed 's/Root;rootrank;//' | sed 's/;Kingdom;/;/' | sed 's/;Phylum;/;/' | sed 's/;Class;/;/' | \\" >> bash_script.sh
  echo "  sed 's/;Order;/;/' | sed 's/;Family;/;/' | sed 's/;Genus;//' | \\" >> bash_script.sh
  echo "  sed \"1s/^/Taxa \$ID\n/\") \\" >> bash_script.sh
  echo "  <(awk '/unclassified/{print \$2, \$4}' \${dir}/out/aligned_classified_counts.txt | \\" >> bash_script.sh
  echo "  sed 's/Root;rootrank;//' | sed 's/;Kingdom;/;/' | sed 's/;Phylum;/;/' | sed 's/;Class;/;/' | \\" >> bash_script.sh
  echo "  sed 's/;Order;/;/' | sed 's/;Family;/;/' | sed 's/;;//') > \${ID}_temp.txt" >> bash_script.sh
  echo "  sed 's/;/ /g' \${ID}_temp.txt | \\" >> bash_script.sh
  echo "  awk '\$1 ~ /unclassified/{\$7 = \$2; \$1 = \"Unclassified_kingdom\"; \$2 = \"Unclassified_phylum\"; \$3 = \"Unclassified_class\"; \$4 = \"Unclassified_order\"; \$5 = \"Unclassified_family\"; \$6 = \"Unclassified_genus\"};1' | \\" >> bash_script.sh
  echo "  awk '\$2 ~ /unclassified/{\$7 = \$3; \$2 = \$1\"_phylum\"; \$3 = \$1\"_class\"; \$4 = \$1\"_order\"; \$5 = \$1\"_family\"; \$6 = \$1\"_genus\"};1' | \\" >> bash_script.sh
  echo "  awk '\$3 ~ /unclassified/{\$7 = \$4; \$3 = \$2\"_class\"; \$4 = \$2\"_order\"; \$5 = \$2\"_family\"; \$6 = \$2\"_genus\"};1' | \\" >> bash_script.sh
  echo "  awk '\$4 ~ /unclassified/{\$7 = \$5; \$4 = \$3\"_order\"; \$5 = \$3\"_family\"; \$6 = \$3\"_genus\"};1' | \\" >> bash_script.sh
  echo "  awk '\$5 ~ /unclassified/{\$7 = \$6; \$5 = \$4\"_family\"; \$6 = \$4\"_genus\"};1' | \\" >> bash_script.sh
  echo "  awk '\$6 ~ /unclassified/{\$7 = \$7; \$6 = \$5\"_genus\"};1' | sed 's/ /;/g' | sed 's/\(.*\);/\1 /' | \\" >> bash_script.sh
  echo "  awk 'FNR==1{print;next}{val=\$2;\$2=\"~\";a[\$0]+=val}"'!'"b[\$0]++{c[++count]=\$0}END{for(i=1;i<=count;i++){sub(\"~\",a[c[i]],c[i]);print c[i]}}' \\" >> bash_script.sh
  echo "  > ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/\${ID}_genus.txt" >> bash_script.sh
  echo "  rm \${ID}_temp.txt" >> bash_script.sh
  echo "done" >> bash_script.sh
  echo "echo \"path <- '${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/'\" > Rfunc.R" >> bash_script.sh
  #echo "echo \"species_files <- list.files(path, 'species.txt')\" >> Rfunc.R" >> bash_script.sh
  #echo "echo \"merged_species <- read.table(paste(path, species_files[1], sep=''), header=T, stringsAsFactors=F, comment.char='')\" >> Rfunc.R" >> bash_script.sh
  #echo "echo \"for (i in 2:length(species_files)){\" >> Rfunc.R" >> bash_script.sh
  #echo "echo \"  file <- read.table(paste(path, species_files[i], sep=''), header=T, stringsAsFactors=F, comment.char='')\" >> Rfunc.R" >> bash_script.sh
  #echo "echo \"  merged_species <- merge(merged_species, file, by='Taxa', all=T)\" >> Rfunc.R" >> bash_script.sh
  #echo "echo \"}\" >> Rfunc.R" >> bash_script.sh
  #echo "echo \"merged_species[is.na(merged_species)] <- 0\" >> Rfunc.R" >> bash_script.sh
  #echo "echo \"write.table(merged_species, paste(path, 'RDP_species_counts.txt', sep=''), row.names=F, quote=F, sep='\t')\" >> Rfunc.R" >> bash_script.sh
  echo "echo \"genus_files <- list.files(path, 'genus.txt')\" >> Rfunc.R" >> bash_script.sh
  echo "echo \"merged_genus <- read.table(paste(path, genus_files[1], sep=''), header=T, stringsAsFactors=F, comment.char='')\" >> Rfunc.R" >> bash_script.sh
  echo "echo \"for (i in 2:length(genus_files)){\" >> Rfunc.R" >> bash_script.sh
  echo "echo \"  file <- read.table(paste(path, genus_files[i], sep=''), header=T, stringsAsFactors=F, comment.char='')\" >> Rfunc.R" >> bash_script.sh
  echo "echo \"  merged_genus <- merge(merged_genus, file, by='Taxa', all=T)\" >> Rfunc.R" >> bash_script.sh
  echo "echo \"}\" >> Rfunc.R" >> bash_script.sh
  echo "echo \"merged_genus[is.na(merged_genus)] <- 0\" >> Rfunc.R" >> bash_script.sh
  echo "echo \"write.table(merged_genus, paste(path, 'RDP_genera_counts.txt', sep=''), row.names=F, quote=F, sep='\t')\" >> Rfunc.R" >> bash_script.sh
  echo "Rscript --vanilla Rfunc.R" >> bash_script.sh
  echo "rm Rfunc.R" >> bash_script.sh
  #echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/*_species.txt" >> bash_script.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/16S_Based_Taxa_Profiling/*_genus.txt" >> bash_script.sh
  chmod +x bash_script.sh
  
  #Run script
  ./bash_script.sh
  rm bash_script.sh
  
  #Signal script and 16S workflow has ended
  echo "Done"
  echo " "
  echo "Running of 16S-based taxonomy classification complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

echo "*** Metagenomics pipeline complete ***"
