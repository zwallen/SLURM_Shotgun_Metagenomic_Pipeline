#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 8 June 2021                                  #
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
#    R base:     For performing functions in pipeline script.#
#    FastQC:     For performing initial quality reports.     #
#    BBMerge:    For merging paired-end reads.               #
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
# Usage:                                                     #
# SLURM_Shotgun_Metagenomic_Pipeline.sh -i input_seqs_dir \  #
#                    -o output_dir \                         #
#                    -p 'commands; to; load; programs' \     #
#                    -r path/to/host/ref/files/dir \         #
#                    -c path/to/chocophlan/dir \             #
#                    -u path/to/uniref/dir \                 #
#                    -t path/to/clade/marker/info/file \     #
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
#     -t    (Required) Path to clade marker info file        #
#           mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2    #
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
#     -m    (Optional) Merge paired-end reads before         #
#           performing the pipeline using BBMerge.           #
#     -a    (Optional) Path to adapters.fa file that comes   #
#           packaged with BBMerge and BBDuk. Required when   #
#           merging reads.                                   #
#     -s    (Optional) Skip certain steps in the pipeline if #
#           need be. Provide a comma separated list of steps #
#           that you wish to skip in the pipeline. List may  #
#           have the values: fastqc, bbduk, kneaddata, humann#
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 8 June 2021                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hi:o:p:r:c:u:t:f:ma:s:" opt; do
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
    echo "    R base:     For performing functions in pipeline script."
    echo "    FastQC:     For performing initial quality reports.     "
    echo "    BBMerge:    For merging paired-end reads.               "
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
    echo " Usage:                                                     "
    echo " SLURM_Shotgun_Metagenomic_Pipeline.sh -i input_seqs_dir \  "
    echo "                    -o output_dir \                         "
    echo "                    -p 'commands; to; load; programs' \     "
    echo "                    -r path/to/host/ref/files/dir \         "
    echo "                    -c path/to/chocophlan/dir \             "
    echo "                    -u path/to/uniref/dir \                 "
    echo "                    -t path/to/clade/marker/info/file \     "
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
    echo "     -t    (Required) Path to clade marker info file        "
    echo "           mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2    "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo "     -m    (Optional) Merge paired-end reads before         "
    echo "           performing the pipeline using BBMerge.           "
    echo "     -a    (Optional) Path to adapters.fa file that comes   "
    echo "           packaged with BBMerge and BBDuk. Required when   "
    echo "           merging reads.                                   "
    echo "     -s    (Optional) Skip certain steps in the pipeline if "
    echo "           need be. Provide a comma separated list of steps "
    echo "           that you wish to skip in the pipeline. List may  "
    echo "           have the values: fastqc, bbduk, kneaddata, humann"
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
    t) MARKERS="$OPTARG"
    ;;
    f) FAIL_EMAIL="$OPTARG"
    ;;
    m) MERGE=1
    ;;
    a) ADAPTERS="$OPTARG"
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

# -t
if [[ -z "$MARKERS" ]]; then
  echo "ERROR: Argument -t is required, please supply path to the clade marker info file mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2"
  exit 1
fi
if [[ -d "$MARKERS" ]]; then
  echo "ERROR: Argument -t should be the path to a single file, not a directory, please supply path to the clade marker info file mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2"
  exit 1
fi
if echo $MARKERS | grep -q -v "mpa_v30_CHOCOPhlAn_201901_marker_info\.txt\.bz2"; then
  echo "ERROR: path given to -t does not contain the file name mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2, please supply the clade marker info file mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2 to this argument"
  exit 1
fi

# -f
if [[ -z "$FAIL_EMAIL" ]]; then
  echo "ERROR: Argument -f is required, please supply an email that can be notified upon failure of any jobs ran during the pipeline"
  exit 1
elif echo $FAIL_EMAIL | grep -q -v '@'; then
  echo "ERROR: Argument -f requires a valid email, please give an email in the form of xxxx@xxxx.xxx"
fi

# -m
if [[ ! -z "$MERGE" ]]; then
  if [[ -z "$ADAPTERS" ]]; then
    echo "ERROR: when specifying the -m parameter, the -a parameter must also be specified"
    exit 1
  fi
fi

# -a
if [[ ! -z "$ADAPTERS" ]]; then
  if [[ -z "$MERGE" ]]; then
    echo "ERROR: when specifying the -a parameter, the -m parameter must also be specified"
    exit 1
  fi
  if [[ -d "$ADAPTERS" ]]; then
    echo "ERROR: Argument -a should be the path to a single file, not a directory, please supply path to the adapters.fa file that comes with BBMerge and BBDuk"
    exit 1
  fi
  if echo $ADAPTERS | grep -q -v "adapters\.fa"; then
    echo "ERROR: path given to -o does not contain the file name adapters.fa, please supply the adapters.fa file that comes with BBMerge and BBDuk to this argument"
    exit 1
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
  else
    echo "ERROR: Invalid argument given to -s, please specify one or more of: fastqc,bbduk,kneaddata,humann"
    exit 1
  fi
fi

# Load programs
$PROG_LOAD

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
  
  #Create directory
  if [ -d "${RESULTS_DIR}/0.FastQC_Initial_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/0.FastQC_Initial_Reports
	  mkdir ${RESULTS_DIR}/0.FastQC_Initial_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/0.FastQC_Initial_Reports/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Running FastQC on input fastq files ***"
  echo " "
  
  #Create directory for output
  if [ -d "${RESULTS_DIR}/0.FastQC_Initial_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/0.FastQC_Initial_Reports
	  mkdir ${RESULTS_DIR}/0.FastQC_Initial_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/0.FastQC_Initial_Reports/0.Output
  fi
  
  ##### Run FastQC #####
  #Create script for running program and submit
  echo '#!/bin/bash' > bash_script.sh
  echo "#SBATCH --partition=express" >> bash_script.sh
  echo "#SBATCH --job-name=FastQC" >> bash_script.sh
  echo "#SBATCH --error=${RESULTS_DIR}/0.FastQC_Initial_Reports/0.ErrorOut/FastQC_%A_%a.err" >> bash_script.sh
  echo "#SBATCH --output=${RESULTS_DIR}/0.FastQC_Initial_Reports/0.Output/FastQC_%A_%a.out" >> bash_script.sh
  echo "#SBATCH --time=2:00:00" >> bash_script.sh
  echo "#SBATCH --ntasks=1" >> bash_script.sh
  echo "#SBATCH --cpus-per-task=1" >> bash_script.sh
  echo "#SBATCH --mem-per-cpu=32000" >> bash_script.sh
  echo "#SBATCH --mail-type=FAIL" >> bash_script.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> bash_script.sh
  echo "#SBATCH --array=1-$(ls -l ${SEQ_DIR}/*R1_001.${SEQ_EXT} | wc -l)" >> bash_script.sh
  echo "#SBATCH --wait" >> bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE1=\$(ls ${SEQ_DIR}/*R1_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
  echo "FILE2=\$(ls ${SEQ_DIR}/*R2_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1_001' '{print \$1}')" >> bash_script.sh
  echo "fastqc \$FILE1 \$FILE2 -d ${RESULTS_DIR}/0.FastQC_Initial_Reports -o ${RESULTS_DIR}/0.FastQC_Initial_Reports \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/0.FastQC_Initial_Reports/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh > /dev/null
  
  rm bash_script.sh
  
  #Signal jobs have ended
  echo "FastQC reports complete"
  echo " "
  
  ##### Get mean Q score per file #####
  echo "Grabbing average quality score per sequence file..."
  echo " "
  echo 'Filename Mean SD' > ${RESULTS_DIR}/0.FastQC_Initial_Reports/mean_Q_per_file.txt
  for file in ${RESULTS_DIR}/0.FastQC_Initial_Reports/*fastqc.zip
  do
    DIR_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '.zip' '{print $1}')
    unzip -qq $file
    sed -n '/>>Per sequence quality scores/,/>>END_MODULE/p' ${DIR_NAME}/fastqc_data.txt | \
    sed '1,2d' | sed '$d' > q_vals.txt
    echo "q_vals <- read.table('q_vals.txt')" > Rfunc.R
    echo "cat(paste(round(mean(rep(q_vals[,1],q_vals[,2])),2), round(sd(rep(q_vals[,1],q_vals[,2])),2)), '\n')" >> Rfunc.R
    echo $(echo $file | awk -F '/' '{print $NF}' | awk -F '_fastqc.zip' '{print $1}') $(Rscript --vanilla Rfunc.R) >> ${RESULTS_DIR}/0.FastQC_Initial_Reports/mean_Q_per_file.txt
    rm -r ${DIR_NAME}
    rm Rfunc.R
    rm q_vals.txt
  done
  echo "Done"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

######### MERGE PAIRED READS WITH BBMERGE #######
if [[ -z "$MERGE" ]]; then
  echo " "
  echo "*** Skipping running of BBMerge for paired read merging ***"
  echo "*** Making empty directory to keep numbering consistent ***"
  echo " "
  
  #Create directory
  if [ -d "${RESULTS_DIR}/1.Merged_Paired_End_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/1.Merged_Paired_End_Sequences
	  mkdir ${RESULTS_DIR}/1.Merged_Paired_End_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/1.Merged_Paired_End_Sequences/0.Output
  fi
else
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
  echo "#SBATCH --cpus-per-task=2" >> bash_script.sh
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
  echo "adapters=${ADAPTERS} \\" >> bash_script.sh
  echo "rem iterations=5 extend2=20 ecct t=2 -Xmx64g \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/1.Merged_Paired_End_Sequences/\${FILE_NAME}.log 2>&1"  >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh > /dev/null
  
  rm bash_script.sh
  
  #Signal jobs have ended
  echo "Merging of paired end reads with BBMerge complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

################# QC WITH BBDUK #################
if echo $SKIP | grep -q "bbduk"; then
  echo " "
  echo "*** Skipping running of BBDuk for adapter/quality trimming and filtering of input fastq files ***"
  echo " "
  
  #Create directory
  if [ -d "${RESULTS_DIR}/2.Quality_Controlled_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/2.Quality_Controlled_Sequences
	  mkdir ${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.Output
  fi
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
  #Create script for running program and submit
  echo '#!/bin/bash' > bash_script.sh
  echo "#SBATCH --partition=short" >> bash_script.sh
  echo "#SBATCH --job-name=QC" >> bash_script.sh
  echo "#SBATCH --error=${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.ErrorOut/QC_%A_%a.err" >> bash_script.sh
  echo "#SBATCH --output=${RESULTS_DIR}/2.Quality_Controlled_Sequences/0.Output/QC_%A_%a.out" >> bash_script.sh
  echo "#SBATCH --time=12:00:00" >> bash_script.sh
  echo "#SBATCH --ntasks=1" >> bash_script.sh
  echo "#SBATCH --cpus-per-task=1" >> bash_script.sh
  echo "#SBATCH --mem-per-cpu=32000" >> bash_script.sh
  echo "#SBATCH --mail-type=FAIL" >> bash_script.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> bash_script.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/1.Merged_Paired_End_Sequences/*.fastq.gz | wc -l)" >> bash_script.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${SEQ_DIR}/*R1_001.${SEQ_EXT} | wc -l)" >> bash_script.sh
  fi
  echo "#SBATCH --wait" >> bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/1.Merged_Paired_End_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> bash_script.sh
    echo "bbduk.sh in=\$FILE \\" >> bash_script.sh
    echo "out=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}.fastq.gz \\" >> bash_script.sh
    echo "stats=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}_stats.txt \\" >> bash_script.sh
    echo "ftm=5 qtrim=rl trimq=25 minlen=50 ref=adapters,phix -Xmx32g \\" >> bash_script.sh
  else
    echo "FILE1=\$(ls ${SEQ_DIR}/*R1_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE2=\$(ls ${SEQ_DIR}/*R2_001.${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1_001' '{print \$1}')" >> bash_script.sh
    echo "bbduk.sh in=\$FILE1 \\" >> bash_script.sh
    echo "in2=\$FILE2 \\" >> bash_script.sh
    echo "out=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}_R1_001.fastq.gz \\" >> bash_script.sh
    echo "out2=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}_R2_001.fastq.gz \\" >> bash_script.sh
    echo "stats=${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}_stats.txt \\" >> bash_script.sh
    echo "ftm=5 tpe tbo qtrim=rl trimq=25 minlen=50 ref=adapters,phix -Xmx32g \\" >> bash_script.sh
  fi
  echo "> ${RESULTS_DIR}/2.Quality_Controlled_Sequences/\${FILE_NAME}.log 2>&1"  >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh > /dev/null
  
  rm bash_script.sh
  
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
  
  #Create directory
  if [ -d "${RESULTS_DIR}/3.Decontaminated_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences
	  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/3.Decontaminated_Sequences/0.Output
  fi
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
    echo "kneaddata --input \$FILE \\" >> bash_script.sh
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*R1_001.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/2.Quality_Controlled_Sequences/*R2_001.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1_001' '{print \$1}')" >> bash_script.sh
    echo "kneaddata --input \$FILE1 \\" >> bash_script.sh
    echo "--input \$FILE2 \\" >> bash_script.sh
  fi
  echo "--output ${RESULTS_DIR}/3.Decontaminated_Sequences \\" >> bash_script.sh
  echo "--output-prefix \$FILE_NAME \\" >> bash_script.sh
  echo "--log ${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}_kneaddata.log \\" >> bash_script.sh
  echo "--reference-db $HOST_REF \\" >> bash_script.sh
  echo "--bypass-trim \\" >> bash_script.sh
  echo "--threads 2 \\" >> bash_script.sh
  echo "--verbose \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/3.Decontaminated_Sequences/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh > /dev/null
  
  rm bash_script.sh
  
  if [[ -z "$MERGE" ]]; then
    rm ${RESULTS_DIR}/3.Decontaminated_Sequences/*unmatched*fastq
  fi
  
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

###### TAXONOMIC AND FUNCTIONAL PROFILING #######
if echo $SKIP | grep -q "humann"; then
  echo " "
  echo "*** Skipping running of HUMAnN/MetaPhlAn workflow for performing taxonomic and functional profiling ***"
  echo " "
else
  SECONDS=0
  echo " "
  echo "*** Running HUMAnN/MetaPhlAn workflow for performing taxonomic and functional profiling ***"
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
  #Create script for running program and submit
  echo '#!/bin/bash' > bash_script.sh
  echo "#SBATCH --partition=medium" >> bash_script.sh
  echo "#SBATCH --job-name=Profiling" >> bash_script.sh
  echo "#SBATCH --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/Profiling_%A_%a.err" >> bash_script.sh
  echo "#SBATCH --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/Profiling_%A_%a.out" >> bash_script.sh
  echo "#SBATCH --time=50:00:00" >> bash_script.sh
  echo "#SBATCH --ntasks=1" >> bash_script.sh
  echo "#SBATCH --cpus-per-task=4" >> bash_script.sh
  echo "#SBATCH --mem-per-cpu=32000" >> bash_script.sh
  echo "#SBATCH --mail-type=FAIL" >> bash_script.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> bash_script.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/3.Decontaminated_Sequences/*.fastq.gz | wc -l)" >> bash_script.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/3.Decontaminated_Sequences/*paired_1.fastq.gz | wc -l)" >> bash_script.sh
  fi
  echo "#SBATCH --wait" >> bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/3.Decontaminated_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> bash_script.sh
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/3.Decontaminated_Sequences/*paired_1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/3.Decontaminated_Sequences/*paired_2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_paired_1' '{print \$1}')" >> bash_script.sh
    echo "zcat \$FILE1 \$FILE2 > ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.temp.fastq" >> bash_script.sh
    echo "FILE=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.temp.fastq" >> bash_script.sh
  fi
  echo "humann --input \$FILE \\" >> bash_script.sh
  echo "--output ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling \\" >> bash_script.sh
  echo "--output-basename \$FILE_NAME \\" >> bash_script.sh
  echo "--metaphlan-options '-t rel_ab --add_viruses' \\" >> bash_script.sh
  echo "--nucleotide-database $CHOCO \\" >> bash_script.sh
  echo "--protein-database $UNIREF \\" >> bash_script.sh
  echo "--prescreen-threshold 0.001 \\" >> bash_script.sh
  echo "--threads 4 \\" >> bash_script.sh
  echo "--verbose \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.temp.fastq" >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh > /dev/null
  
  rm bash_script.sh
  
  ##### Run MetaPhlAn again to get normalized counts #####
  #Create script for running program and submit
  echo '#!/bin/bash' > bash_script.sh
  echo "#SBATCH --partition=short" >> bash_script.sh
  echo "#SBATCH --job-name=Norm_Abun" >> bash_script.sh
  echo "#SBATCH --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/Norm_Abun_%A_%a.err" >> bash_script.sh
  echo "#SBATCH --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/Norm_Abun_%A_%a.out" >> bash_script.sh
  echo "#SBATCH --time=12:00:00" >> bash_script.sh
  echo "#SBATCH --ntasks=1" >> bash_script.sh
  echo "#SBATCH --cpus-per-task=1" >> bash_script.sh
  echo "#SBATCH --mem-per-cpu=32000" >> bash_script.sh
  echo "#SBATCH --mail-type=FAIL" >> bash_script.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> bash_script.sh
  echo "#SBATCH --array=1-$(ls -d -l ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*humann_temp | wc -l)" >> bash_script.sh
  echo "#SBATCH --wait" >> bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "DIR=\$(ls -d ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*humann_temp | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$DIR | awk -F '/' '{print \$NF}' | awk -F '_humann_temp' '{print \$1}')" >> bash_script.sh
  echo "metaphlan --input_type bowtie2out \\" >> bash_script.sh
  echo "--add_viruses \\" >> bash_script.sh
  echo "-t marker_ab_table \\" >> bash_script.sh
  echo "\${DIR}/\${FILE_NAME}_metaphlan_bowtie2.txt \\" >> bash_script.sh
  echo "\${DIR}/\${FILE_NAME}_metaphlan_norm_abun_table.tsv \\" >> bash_script.sh
  echo ">> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  echo "join -j 1 -o 1.3,1.2,2.2 \\" >> bash_script.sh
  echo "<(sort -k1,1 <(paste <(bzcat $MARKERS | awk -F\"\t\" '{print \$1}') \\" >> bash_script.sh
  echo "                     <(bzcat $MARKERS | awk -F\"__\" '{print \$1}') \\" >> bash_script.sh
  echo "                     <(bzcat $MARKERS | awk -F\"('taxon': '|'})\" '{print \$2}'))) \\" >> bash_script.sh
  echo "<(sort -k1,1 \${DIR}/\${FILE_NAME}_metaphlan_norm_abun_table.tsv) | \\" >> bash_script.sh
  echo "sed '1s/^/#clade_name NCBI_tax_id normalized_abundance\n/' | sed 's/ /\t/g' | \\" >> bash_script.sh
  echo "awk 'FNR==1{print;next}{val=\$3;\$3=\"~\";a[\$0]+=val}"'!'"b[\$0]++{c[++count]=\$0}END{for(i=1;i<=count;i++){sub(\"~\",a[c[i]],c[i]);print c[i]}}' OFS='\t' | \\" >> bash_script.sh
  echo "cat <(grep '^#' \${DIR}/\${FILE_NAME}_metaphlan_norm_abun_table.tsv) - > \${FILE_NAME}_temp.txt" >> bash_script.sh
  echo "mv \${FILE_NAME}_temp.txt \${DIR}/\${FILE_NAME}_metaphlan_norm_abun_table.tsv" >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh > /dev/null
  
  rm bash_script.sh
  
  ##### Clean up #####
  echo "Cleaning things up... merging per sample tables and removing unneeded files"
  echo " "
  for dir in ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*humann_temp; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_humann_temp' '{print $1}')
    
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
  
  #Create individual scripts for running table merging
  echo '#!/bin/bash' > bash_script_1.sh
  echo "$PROG_LOAD" >> bash_script_1.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_rel_abun_table.tsv \\" >> bash_script_1.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_rel_abun_table.tsv" >> bash_script_1.sh
  echo "sed -i '2s/_metaphlan_rel_abun_table//g' ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_rel_abun_table.tsv" >> bash_script_1.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_rel_abun_table.tsv" >> bash_script_1.sh
  chmod +x bash_script_1.sh
  
  echo '#!/bin/bash' > bash_script_2.sh
  echo "$PROG_LOAD" >> bash_script_2.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_norm_abun_table.tsv \\" >> bash_script_2.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_norm_abun_table.tsv" >> bash_script_2.sh
  echo "sed -i '2s/_metaphlan_norm_abun_table//g' ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_norm_abun_table.tsv" >> bash_script_2.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_norm_abun_table.tsv" >> bash_script_2.sh
  chmod +x bash_script_2.sh
  
  echo '#!/bin/bash' > bash_script_3.sh
  echo "$PROG_LOAD" >> bash_script_3.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/ \\" >> bash_script_3.sh
  echo "--output ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/humann_genefamilies.tsv \\" >> bash_script_3.sh
  echo "--file_name genefamilies.tsv" >> bash_script_3.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*genefamilies.tsv" >> bash_script_3.sh
  chmod +x bash_script_3.sh
  
  echo '#!/bin/bash' > bash_script_4.sh
  echo "$PROG_LOAD" >> bash_script_4.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/ \\" >> bash_script_4.sh
  echo "--output ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/humann_pathabundance.tsv \\" >> bash_script_4.sh
  echo "--file_name pathabundance.tsv" >> bash_script_4.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*pathabundance.tsv" >> bash_script_4.sh
  chmod +x bash_script_4.sh

  echo '#!/bin/bash' > bash_script_5.sh
  echo "$PROG_LOAD" >> bash_script_5.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/ \\" >> bash_script_5.sh
  echo "--output ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/humann_pathcoverage.tsv \\" >> bash_script_5.sh
  echo "--file_name pathcoverage.tsv" >> bash_script_5.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*pathcoverage.tsv" >> bash_script_5.sh
  chmod +x bash_script_5.sh
  
  #Create script for running and submitting individual scripts
  echo '#!/bin/bash' > bash_script.sh
  echo "#SBATCH --partition=short" >> bash_script.sh
  echo "#SBATCH --job-name=Merge_tables" >> bash_script.sh
  echo "#SBATCH --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/Merge_tables_%A_%a.err" >> bash_script.sh
  echo "#SBATCH --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/Merge_tables_%A_%a.out" >> bash_script.sh
  echo "#SBATCH --time=12:00:00" >> bash_script.sh
  echo "#SBATCH --ntasks=1" >> bash_script.sh
  echo "#SBATCH --cpus-per-task=1" >> bash_script.sh
  echo "#SBATCH --mem-per-cpu=32000" >> bash_script.sh
  echo "#SBATCH --mail-type=FAIL" >> bash_script.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> bash_script.sh
  echo "#SBATCH --array=1-5" >> bash_script.sh
  echo "#SBATCH --wait" >> bash_script.sh
  echo "./bash_script_\${SLURM_ARRAY_TASK_ID}.sh" >> bash_script.sh
  chmod +x bash_script.sh
  
  sbatch bash_script.sh > /dev/null
  
  rm bash_script*sh
  
  #Signal workflow has completed
  echo "Running of HUMAnN workflow complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

echo "*** Metagenomics pipeline complete ***"
