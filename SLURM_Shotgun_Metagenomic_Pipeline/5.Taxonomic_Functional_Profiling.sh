#!/bin/bash
set -e

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 25 May 2021                                  #
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
# ./5.Taxonomic_Functional_Profiling.sh -o output_dir \      #
#                        -p 'commands; to; load; programs' \ #
#                        -c path/to/chocophlan/dir \         #
#                        -u path/to/uniref/dir \             #
#                        -m path/to/clade/marker/info/file \ #
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
#     -c    (Required) Path to ChocoPhlAn database directory.#
#     -u    (Required) Path to UniRef90 database directory.  #
#     -m    (Required) Path to clade marker info file        #
#           mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2    #
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
while getopts ":ho:p:c:u:m:f:" opt; do
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
    echo " ./5.Taxonomic_Functional_Profiling.sh -o output_dir \      "
    echo "                        -p 'commands; to; load; programs' \ "
    echo "                        -c path/to/chocophlan/dir \         "
    echo "                        -u path/to/uniref/dir \             "
    echo "                        -m path/to/clade/marker/info/file \ "
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
    echo "     -c    (Required) Path to ChocoPhlAn database directory."
    echo "     -u    (Required) Path to UniRef90 database directory.  "
    echo "     -m    (Required) Path to clade marker info file        "
    echo "           mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2    "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo " "
    exit 0
    ;;
    o) RESULTS_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    p) PROG_LOAD="$OPTARG"
    ;;
    c) CHOCO=$(echo $OPTARG | sed 's#/$##')
    ;;
    u) UNIREF=$(echo $OPTARG | sed 's#/$##')
    ;;
    m) MARKERS="$OPTARG"
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
  #Create shell script for running program
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$1 | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> bash_script.sh
  echo "humann --input \$1 \\" >> bash_script.sh
  echo "--output ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling \\" >> bash_script.sh
  echo "--output-basename \$FILE_NAME \\" >> bash_script.sh
  echo "--metaphlan-options '-t rel_ab --add_viruses' \\" >> bash_script.sh
  echo "--nucleotide-database $CHOCO \\" >> bash_script.sh
  echo "--protein-database $UNIREF \\" >> bash_script.sh
  echo "--prescreen-threshold 0.001 \\" >> bash_script.sh
  echo "--threads 5 \\" >> bash_script.sh
  echo "--verbose \\" >> bash_script.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/\${FILE_NAME}.log 2>&1" >> bash_script.sh
  chmod +x bash_script.sh
  
  #For every decontaminated sequence files submit job and grab job IDs
  touch job_ids.txt
  for file in ${RESULTS_DIR}/3.Decontaminated_Sequences/*.fastq.gz; do
    FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '.fastq.gz' '{print $1}')
    
    sbatch --partition=medium \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/${FILE_NAME}_humann.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/${FILE_NAME}_humann.out \
    --time=50:00:00 \
    --ntasks=1 \
    --cpus-per-task=5 \
    --mem-per-cpu=32000 \
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
  
  ##### Run MetaPhlAn again to get normalized counts #####
  #Create shell script for running program
  echo '#!/bin/bash' > bash_script.sh
  echo "$PROG_LOAD" >> bash_script.sh
  echo "FILE_NAME=\$(echo \$1 | awk -F '/' '{print \$NF}' | awk -F '_humann_temp' '{print \$1}')" >> bash_script.sh
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
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_humann_temp' '{print $1}')
    
    sbatch --partition=short \
    --job-name=${FILE_NAME} \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/${FILE_NAME}_metaphlan.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/${FILE_NAME}_metaphlan.out \
    --time=12:00:00 \
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
  
  touch job_ids.txt
  
  echo '#!/bin/bash' > bash_script_1.sh
  echo "$PROG_LOAD" >> bash_script_1.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_rel_abun_table.tsv \\" >> bash_script_1.sh
  echo "> ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_rel_abun_table.tsv" >> bash_script_1.sh
  echo "sed -i '2s/_metaphlan_rel_abun_table//g' ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/Merged_Sample_Tables/metaphlan_rel_abun_table.tsv" >> bash_script_1.sh
  echo "rm ${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/*metaphlan_rel_abun_table.tsv" >> bash_script_1.sh
  
  sbatch --partition=short \
    --job-name=metaphlan_rel_abun_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/metaphlan_rel_abun_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/metaphlan_rel_abun_merge.out \
    --time=12:00:00 \
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
  
  sbatch --partition=short \
    --job-name=metaphlan_norm_abun_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/metaphlan_norm_abun_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/metaphlan_norm_abun_merge.out \
    --time=12:00:00 \
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
  
  sbatch --partition=short \
    --job-name=genefamilies_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/genefamilies_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/genefamilies_merge.out \
    --time=12:00:00 \
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
  
  sbatch --partition=short \
    --job-name=pathabundance_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/pathabundance_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/pathabundance_merge.out \
    --time=12:00:00 \
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
  
  sbatch --partition=short \
    --job-name=pathcoverage_merge \
    --error=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.ErrorOut/pathcoverage_merge.err \
    --output=${RESULTS_DIR}/4.Taxonomic_and_Functional_Profiling/0.Output/pathcoverage_merge.out \
    --time=12:00:00 \
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
#################################################
