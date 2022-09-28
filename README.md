#### Microbiome and Genetics Analysis Tool (MaGAT)
MaGAT is a wrapper program that helps to facilitate the analysis of microbiome and host genetic data. Analyses that can be run include association of microbiome data with host genetic variants with or without interaction from additional variables. Further descriptions of wrapper programs are below:

- `MaGAT.sh` - This the main wrapper program that wraps different R packages and PLINK 2 [https://www.cog-genomics.org/plink/2.0/] to perform association analyses between microbiome count data (e.g., that derived from 16S rRNA gene amplicon or shotgun metagenomic sequencing) and host genetic variation with or without covariate adjustment and/or interaction with other variables. This script does not introduce anything novel, rather, it brings together existing methods and their programs to provide a convenient and user friendly way to run complex analyses involving micorbiome and host genetic data. Run `./MaGAT.sh -h` for further description and list of parameters.

- `Meta_MaGAT.sh` - This is a script that wraps two programs for performing meta-analyses (METASOFT and METAL). It takes as input results from running `MaHGAT.sh` on two or more datasets. For most meta-analyses, the script implements the program METASOFT [http://genetics.cs.ucla.edu/meta_jemdoc/index.html] for performing fixed- and random-effects meta-analysis between results for two or more datasets. For meta-analyses involving joint test results (see PLINK 2 documentation [https://www.cog-genomics.org/plink/2.0/assoc#glm]), METAL [https://genome.sph.umich.edu/wiki/METAL] is used to perform a sample weighted Z-score-based p-value meta-analysis as no test statistics are given for joint tests with PLINK 2 output. Run `./Meta_MaGAT.sh -h` for list of parameters.

- `LocusZoom_MaGAT.sh` - This is a convenience script that takes as input the output from `MaGAT.sh` or `Meta_MaGAT.sh` and formats it for uploading to LocusZoom [https://my.locuszoom.org/].

