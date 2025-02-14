This directory contains a `.yml` file that can be used for building an environment with `conda` with all the necessary packages for running `GaMbIT`, `Meta_GaMbIT`, and the helper scripts. To create the environment (titled `GaMbIT`) use the following command:
```
$ conda env create -f GaMbIT.yml
```
and once the environment is built, activate it with
```
$ conda activate GaMbIT
```
The only thing that still needs to be done once the conda environment is built is to install the `phyloseq` R package
```
> BiocManager::install('phyloseq')
```
`METASOFT` and `METAL` are not available to download via conda, so the most up to date `METASOFT` java file (filename `Metasoft.jar`) and the `METAL` binary for Linux (filename `metal`) are provided in this directory. Two additional files needed to run `METASOFT` (`plink2metasoft_modified.py` and `HanEskinPvalueTable.txt`) are also provided here.