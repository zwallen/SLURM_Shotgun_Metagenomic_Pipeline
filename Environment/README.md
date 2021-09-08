This directory contains a `.yml` file that can be used for building an environment with `conda` with all the necessary packages for running the `SLURM_Shotgun_Metagenomic_Pipeline`. To create the environment (titled `shotgun_pipeline`) use the following command:
```
conda env create -f shotgun_pipeline.yml
```
and once the environment is built, activate it with
```
conda activate shotgun_pipeline
```
and make sure you move the `resources` directory found in this directory to the `shotgun_pipeline/bin` directory as these reference files are needed when running BBDuk.

If this environment is built, then `conda activate shotgun_pipeline` should suffice for the `-p` parameter in the `SLURM_Shotgun_Metagenomic_Pipeline` scripts.
