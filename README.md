This github repository houses a wrapper program (`SLURM_Shotgun_Metagenomic_Pipeline.sh`) for processing shotgun metagenomic sequences (derived from Illumina paired-end whole genome shotgun sequencing) to taxonomic (relative abundances and normalized abundance counts) and functional (gene and pathway) profiles on a high performance computing cluster using a SLURM scheduling system. The overall pipeline includes performing an intitial quality assessment on raw sequences using FastQC [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/], sequence processing using BBMap suite of tools [https://sourceforge.net/projects/bbmap/] including merging of paired-end reads using BBMerge (optional), adapter removal and quality trimming/filtering using BBDuk, removal of host contaminant reads using BBSplit/BBMap, and lastly taxonomic and functional profiling using MetaPhlAn/HUMAnN workflow [https://huttenhower.sph.harvard.edu/humann/]. A convenience script is also located in the repository that wraps GraPhlAn [https://huttenhower.sph.harvard.edu/graphlan/] for generating a cladogram of top most abundant clades detected by MetaPhlAn (`Create_Cladogram.sh`).

The following gives an overview of the overall structure of the repository:

## Directory tree for repository
```
SLURM_Shotgun_Metagenomic_Pipeline
|
|-- Environment -- Directory that contains a .yml file to create a conda environment with all the necessary packages for running the pipeline scripts.
|                  Optional to build this environment, provided for convenience to ensure proper programs and versions are downloaded.
|
|-- Reference_Files -- Directory that contains the shell script 0.get_reference_files.sh.
|                      Running 0.get_reference_files.sh will download all the necessary and most up to date reference files and databases
|                      used in the pipeline.
|
|-- SLURM_Shotgun_Metagenomic_Pipeline -- Directory that contains separate shell scripts for each step of the pipeline.
|                                         Instead of running the whole pipeline at once using the main wrapper script (SLURM_Shotgun_Metagenomic_Pipeline.sh),
|                                         one can run it in chunks using these scripts. Useful for when the dataset is large, and will not finish running
|                                         with using the main wrapper script, or if the main wrapper script failed or had to be stopped at a certain step.
|                                         Each pipeline chunk comes with its own example job script for submitting it to a SLURM scheduler.
|
|-- Create_Cladogram.sh -- Convenience script for creating a cladogram with GraPhlAn. Takes the MetaPhlAn relative abundance table generated
|                          during running of HUMAnN as input.
|
|-- SLURM_Shotgun_Metagenomic_Pipeline.job -- An example sbatch job script for submitting the SLURM_Shotgun_Metagenomic_Pipeline.sh to a
|                                             SLURM scheduling system on a high performance computing cluster.
|
|-- SLURM_Shotgun_Metagenomic_Pipeline.sh -- The wrapper program that runs the metagenomic pipeline.

```
## Important notes about the pipeline program

### Submission of internal sbatch jobs
The `SLURM_Shotgun_Metagenomic_Pipeline.sh` script will internally submit jobs for each step of the pipeline. For each job step, a `bash_script.sh` file is created in the current directory that contains the code for the currently running step, and is deleted once the step completes. Default partitions and time-limits for these jobs have been chosen based on previous experience, but there is always a chance that jobs might take longer if sequence files being processed are larger than what the program was tested with. This pipeline was created and tested using 372 paired sequence files (744 total) derived from paired-end 150bp sequencing with a target of 40M total reads per sample (gzipped sequence file sizes ranged from 1.4-4.5G). If needing to modify the partitions and/or time-limits, one can simply pop open the pipeline script and modify where needed.

### Required programs/databases and parameter descriptions
For descriptions of required programs/databases and parameters for `SLURM_Shotgun_Metagenomic_Pipeline.sh` or `Create_Cladogram.sh` scripts, run the respective script with parameter `-h`.

### Separate shell scripts for individual pipeline steps
The directory `SLURM_Shotgun_Metagenomic_Pipeline/` contains separate shell scripts that can perform each step of the pipeline individually, and are numbered in the sequence they are performed in the wrapper script `SLURM_Shotgun_Metagenomic_Pipeline.sh`. These have been provided as an alternative to running the whole pipeline in one shot. Using these can be useful if the dataset being processed is very large, and will not complete in time using the one shot `SLURM_Shotgun_Metagenomic_Pipeline.sh` wrapper script. Additionally, these are useful if the main wrapper script has failed, or stopped, at a certain step, and you want to continue with the pipeline without having to re-run the `SLURM_Shotgun_Metagenomic_Pipeline.sh` script.

### Generating a report for tracking number of sequences through pipeline
The `SLURM_Shotgun_Metagenomic_Pipeline/` that houses separate shell scripts for each step of the pipeline also houses a script that can be used to generate a report detailing how many sequences were outputted from each step of the pipeline (`6.Pipeline_Report.sh`).

```
./SLURM_Shotgun_Metagenomic_Pipeline.sh -h

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# by Zachary D Wallen                                        #
# Last updated: 8 Sep 2021                                   #
##############################################################
 
 Description: This is a wrapper program that wraps various  
 programs to process raw paired-end whole genome shotgun    
 metagenomic sequences. The end product of this pipeline    
 is taxonomic and functional (gene and pathway) abundances  
 that are ready for further statistical analyses.           
                                                            
 Required programs and databases:                           
    SLURM:      Program is designed to work with a SLURM    
                high performance computing cluster          
                scheduling system.                          
    R base:     For performing functions in pipeline script.
    FastQC:     For performing initial quality reports.     
    BBMerge:    For merging paired-end reads.               
    BBDuk:      For adapter and quality trimming of raw wgs 
                reads. Also can remove PhiX sequences.      
    BBMap/                                                  
    BBSplit:    For removing host contamination from wgs    
                reads. Requires FASTA genome reference file 
                to map reads against.                       
    HUMAnN:     For generating taxonomic, gene family, and  
                pathway abundances.                         
    ChocoPhlAn: Database used for taxonomic profiling. Can  
                be any of the ChocoPhlAn databases          
                downloaded using humann_databases utility   
                program.                                    
    UniRef:     Database used for functional profiling. Can 
                be any of the UniRef databases downloaded   
                using humann_databases utility program.     
    Markers:    File with ChocoPhlAn GeneIDs and marker     
                taxonomies. Used for replacing GeneIDs with 
                clade marker taxonomy lineages in the       
                normalized abundance tables from MetaPhlAn. 
                                                            
 Usage:                                                     
 SLURM_Shotgun_Metagenomic_Pipeline.sh -i input_seqs_dir \  
                    -o output_dir \                         
                    -p 'commands; to; load; programs' \     
                    -r path/to/host/ref/file.fa \           
                    -c path/to/chocophlan/dir \             
                    -u path/to/uniref/dir \                 
                    -t path/to/clade/marker/info/file \     
                    -f notificationEmail@forFailures.edu \  
                    [additional options]                    
                                                            
 Parameters:                                                
     -h    Print the parameter list below then exit.        
     -i    (Required) Directory that contains the raw       
           fastq files to be processed. Sequences must have 
           file extensions .fastq OR .fq,                   
           and can be gzipped or not.                       
     -o    (Required) Directory to put output of pipeline   
           into. NOTE: make sure output directory is in an  
           area that has plenty of data storage space       
           available if processing large datasets.          
     -p    (Required) Single quoted string that contains    
           commands to load all the necessary programs      
           needed to run pipeline steps (e.g. activating    
           conda environments, loading modules, adding to   
           PATH, etc.).                                     
     -r    (Required) Path to reference genome file of host.
           Should be in FASTA format, and uncompressed.     
     -c    (Required) Path to ChocoPhlAn database directory.
     -u    (Required) Path to UniRef90 database directory.  
     -t    (Required) Path to clade marker info file        
           mpa_v30_CHOCOPhlAn_201901_marker_info.txt.bz2    
     -f    (Required) E-mail to send notifications to upon  
           failure of any jobs.                             
     -m    (Optional) Merge paired-end reads before         
           performing the pipeline using BBMerge.           
     -a    (Optional) Path to adapters.fa file that comes   
           packaged with BBMerge and BBDuk. Required when   
           merging reads.                                   
     -s    (Optional) Skip certain steps in the pipeline if 
           need be. Provide a comma separated list of steps 
           that you wish to skip in the pipeline. List may  
           have the values: fastqc, QC, decontam, humann
```
```
./Create_Cladogram.sh -h

##############################################################
# Create cladogram of most abundant clades                   #
# by Zachary D Wallen                                        #
# Last updated: 18 March 2021                                #
##############################################################
 
 Description: This is a wrapper program that wraps          
 GraPhlAn program to produce a cladogram summarizing the    
 top most abundant clades detected by MetaPhlAn for a       
 dataset.                                                   
                                                            
 Required programs and databases:                           
    GraPhlAn:   For generating the cladogram. Will also need
                the export2graphlan.py conversion script.   
                                                            
 Usage:                                                     
 Create_Cladogram.sh -i input_metaphlan_rel_abun_table \    
                     -o output_cladogram \                  
                     -s "search_string"                     
                                                            
 Parameters:                                                
     -h    Print the parameter list below then exit.        
     -i    (Required) A taxa by sample relative abundance   
           table created by running MetaPhlAn with 'rel_ab' 
           option on individual sample sequences, then      
           merging individual sample relative abundance     
           tables with merge_metaphlan_tables.py script.    
     -o    (Required) Name for the outputted cladogram.     
     -s    (Optional) Should the relative abundance table be
           subsetted for specific samples? If so, supply a  
           quoted pattern in the sample IDs here that can be
           used to extract the desired sample data columns. 
```
