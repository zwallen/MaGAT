This directory contains a `.yml` file that can be used for building an environment with `conda` with all the necessary packages for running `MaGAT`, `Meta_MaGAT`, and the helper scripts. To create the environment (titled `MaGAT`) use the following command:
```
$ conda env create -f MaGAT.yml
```
and once the environment is built, activate it with
```
$ conda activate MaGAT
```
The only thing that still needs to be done once the conda environment is built is to install the `phyloseq` R package
```
> BiocManager::install('phyloseq')
```
