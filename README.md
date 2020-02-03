This repository contains workflows to estimate the abundance of ORFs in a 
metagenome using protein assemblies derived from the metagenome. 

To run the snakefiles in this repository, install snakemake-minimal

```
conda create -n aa snakemake-minimal=5.9.1
conda activate aa
snakemake -s control.snakefile
```
