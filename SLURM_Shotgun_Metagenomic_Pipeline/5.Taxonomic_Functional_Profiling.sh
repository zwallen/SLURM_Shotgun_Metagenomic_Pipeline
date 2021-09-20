#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 8 Sep 2021                                   #
#                                                            #
# Description: Perform taxonomic and functional profiling    #
# using HUMAnN/MetaPhlAn workflow.                           #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
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
# ./5.Taxonomic_Functional_Profiling.sh [-m] \               #
#                        -o output_dir \                     #
#                        -p 'commands; to; load; programs' \ #
#                        -c path/to/chocophlan/dir \         #
#                        -u path/to/uniref/dir \             #
#                        -t path/to/clade/marker/info/file \ #
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
#     -c    (Required) Path to ChocoPhlAn database directory.#
#     -u    (Required) Path to UniRef90 database directory.  #
#     -t    (Required) Path to clade marker info file        #
#           mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2    #
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# by Zachary D Wallen                                        #"
echo "# Last updated: 8 Sep 2021                                   #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hmo:p:c:u:t:f:" opt; do
  case $opt in
    h)
    echo " Description: Merge paired-end reads using BBMerge.         "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
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
    echo " ./5.Taxonomic_Functional_Profiling.sh [-m] \               "
    echo "                        -o output_dir \                     "
    echo "                        -p 'commands; to; load; programs' \ "
    echo "                        -c path/to/chocophlan/dir \         "
    echo "                        -u path/to/uniref/dir \             "
    echo "                        -t path/to/clade/marker/info/file \ "
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
    echo "     -c    (Required) Path to ChocoPhlAn database directory."
    echo "     -u    (Required) Path to UniRef90 database directory.  "
    echo "     -t    (Required) Path to clade marker info file        "
    echo "           mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2    "
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
    c) CHOCO=$(echo $OPTARG | sed 's#/$##')
    ;;
    u) UNIREF=$(echo $OPTARG | sed 's#/$##')
    ;;
    t) MARKERS="$OPTARG"
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

# Load programs
$PROG_LOAD

###### TAXONOMIC AND FUNCTIONAL PROFILING #######

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
  echo "#SBATCH --cpus-per-task=10" >> bash_script.sh
  echo "#SBATCH --mem-per-cpu=16000" >> bash_script.sh
  echo "#SBATCH --mail-type=FAIL" >> bash_script.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> bash_script.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/3.Decontaminated_Sequences/*.fastq.gz | wc -l)" >> bash_script.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/3.Decontaminated_Sequences/*R1_001.fastq.gz | wc -l)" >> bash_script.sh
  fi
  echo "#SBATCH --wait" >> bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/3.Decontaminated_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> bash_script.sh
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/3.Decontaminated_Sequences/*R1_001.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/3.Decontaminated_Sequences/*R2_001.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> bash_script.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1_001' '{print \$1}')" >> bash_script.sh
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
  echo "--threads 10 \\" >> bash_script.sh
  echo "--verbose \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.temp.fastq" >> bash_script.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}_humann_temp/*_bowtie_*" >> bash_script.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}_humann_temp/*_diamond_*" >> bash_script.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}_humann_temp/*_custom_chocophlan_database.ffn" >> bash_script.sh
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
  echo "cat <(sed -n 1p \${DIR}/\${FILE_NAME}_metaphlan_norm_abun_table.tsv) - > \${FILE_NAME}_temp.txt" >> bash_script.sh
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
    
    mv ${dir}/${FILE_NAME}_metaphlan_bowtie2.txt \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_Profiles/${FILE_NAME}_metaphlan_bowtie2.txt
    
    mv ${dir}/${FILE_NAME}.log \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_Profiles/${FILE_NAME}_humann.log
    
    rm -r ${dir}
    
    mv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}*.* \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_Profiles/
  done
  
  ##### Merge tables #####
  mkdir ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables
  
  for dir in ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*_Profiles; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_Profiles' '{print $1}')

    cp ${dir}/${FILE_NAME}_metaphlan_rel_abun_table.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/

    cp ${dir}/${FILE_NAME}_metaphlan_norm_abun_table.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/

    sed -i 1d ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_metaphlan_norm_abun_table.tsv
    cat <(head -n 3 ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_metaphlan_rel_abun_table.tsv) \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_metaphlan_norm_abun_table.tsv > \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/temp.txt
    mv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/temp.txt \
    ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/${FILE_NAME}_metaphlan_norm_abun_table.tsv

    cp ${dir}/${FILE_NAME}_genefamilies.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/

    cp ${dir}/${FILE_NAME}_pathabundance.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/

    cp ${dir}/${FILE_NAME}_pathcoverage.tsv ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/
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
#################################################
