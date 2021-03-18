This github repository houses scripts used in the maunscript 
**Wallen, ZD. Comparison study of sixteen differential abundance methods using two large Parkinson disease gut microbiome datasets.**

The following README gives an overview of the overall structure of the repository, and important notes on how to run the scripts.

## Directory tree for repository
```
Wallen_DAMethodCompare_2021
|
|-- Dataset1_Scripts -- houses scripts used to run differential abundance tests and create visualizations for Dataset 1
|                       each differential abundance method has its own sub-directory
|
|-- Dataset2_Scripts -- houses scripts used to run differential abundance tests and create visualizations for Dataset 2
|                       each differential abundance method has its own sub-directory
|
|-- Joint_Analyses_Scripts -- houses scripts used to generate heatmaps when both datasets are involved (when accounting for replication of signals)
|  
|-- PhyloseqObjects -- phyloseq objects used in analyses that were reported in manuscript
|   |
|   |-- Dataset1 -- phyloseq object for Dataset 1
|   |
|   |-- Dataset2 -- phyloseq object for Dataset 2
|
|-- Script_Output -- directory to house output from scripts in Dataset1_Scripts, Dataset2_Scripts, and Joint_Analyses_Scripts directories
|   |
|   |-- Dataset1_Output -- output from scripts in Dataset1_Scripts directory
|   |
|   |-- Dataset2_Output -- output from scripts in Dataset2_Scripts directory
|   |
|   |-- Joint_Analyses_Output -- output from scripts in Joint_Analyses_Scripts directory
|
|-- Support_Files -- files called upon by certain scripts
```

## Important notes about this repository

#### The repository is structured to work out of the box
Once downloaded, all scripts should be able to be run as is, without any modification to the repository structure or scripts themselves.

#### Phyloseq objects used in the manuscript are included in this repository
As stated in the directory tree, phyloseq objects used in the manuscript for datasets 1 and 2 are located in the `PhyloseqObjects/` directory. Running scripts in `Dataset1_Scripts/`, `Dataset2_Scripts/`, and `Joint_Analyses_Scripts/` directories using these phyloseq objects should give same results as reported in the manuscript. 

#### Results for analyses utilizing manuscript phyloseq objects are included in this repository
Results that are outputted by analyses scripts and that were reported in the manuscript are located in the `Script_Output/` directory. Any results that are not found in that directory should be in the Supplementary Material of the manuscript.

*WARNING: As the scripts are currently set up, the results that are included in the `Script_Output/` repository will be overwritten if re-running the analyses scripts without moving or renaming the result files first.*

#### There are two types of scripts stored in this repository
Scripts with the extension `.R` are R scripts written in R programming language used to perform differential abundance analyses and generate plots. Scripts that have a `.job` extension are shell scripts that were used to submit `.R` scripts to a SLURM scheduling system on a high performance computing cluster. Each `.R` script should have an accompanying `.job` script.
