#!/bin/bash
set -e

####################################################################
# Microbiome and Genetics Analysis Tool (MaGAT)                    #
# by Zachary D Wallen                                              #
# Last updated: 8 June 2021                                        #
#                                                                  #
# Description: This is a program that wraps different R packages   #
# and PLINK 2 to perform association analyses between microbial    #
# feature count data and genome-wide genotypes.                    #
#                                                                  #
# Input data consists of the following:                            #
#                                                                  #
# (1) A phyloseq object (saved as RDS file) that contains a table  #
# with feature counts per sample in the otu_table() slot,          #
# taxonomic assignments (if applicable) in the tax_table() slot,   #
# and any sample metadata that one wants to include in the analysis#
# in the sample_data() slot.                                       #
#                                                                  #
# (2) Genome-wide genotypes in PLINK binary (.bed, .bim, .fam)     #
# format. Sample names of the phyloseq object and IIDs of the PLINK#
# files should contain the same IDs (they do not have to be in the #
# same order) for the program to work properly. Alternatively,     #
# genotypes may be supplied as a directory of VCF files. Sample IDs#
# for VCF files should be in the format "FID_IID".                 #
#                                                                  #
# WARNING: Files created during analysis will take up a lot of     #
# storage space if a large number of variants and features are     #
# analyzed, make sure you have plenty of free storage where the    #
# output directory is located.                                     #
#                                                                  #
# WARNING: When using VCFs, VCFs will be converted to PLINK format #
# when running this script. The symbol "_" will be taken as an ID  #
# delimiter to split VCF sample IDs into FIDs and IIDs. Please     #
# ensure VCF sample IDs are in the form "FID_IID" and that the     #
# sample names of the input phyloseq object matches the IID portion#
# of the VCF ID.                                                   #
#                                                                  #
# WARNING: If using the -c and/or -m filter for SNPs, be sure to   #
# only include VCF files in the directory that contain these SNPs, #
# else script will fail.                                           #
#                                                                  #
# Usage: ./MaGAT.sh -i phyloseq_object.rds \                      #
#                    -g genotype_file_prefix OR -d vcf_directory \ #
#                    -o output_dir \                               #
#                    [additional options]                          #
#                                                                  #
# Parameters:                                                      #
#     -h    Print the parameter list below then exit.              #
#     -i    (Required) Phyloseq object in RDS format, saved from R #
#           with the writeRDS function. Make sure tax_table() slot #
#           is made up of 7 columns (Kingdom, Phylum, Class, Order,#
#           Family, Genus, Species) if wanting to collapse to      #
#           to certain taxonomic levels and any unclassified levels#
#           have 'NA' for the entry.                               #
#     -g    (Required) File name pre-fix for host genome-wide      #
#           genotypes in PLINK binary format (files with extension #
#           .bim, .bam, .fam). Do not include file extensions,     #
#           only pre-fix.                                          #
#     -d    (Required if using VCFs) Directory containing VCF files#
#           with genotype data. Files can be broken up by          #
#           chromosome and gzipped compressed. All files should    #
#           contain the same samples.                              #
#     -o    (Required) Directory to where the output files will be #
#           written to.                                            #
#     -l    (Optional) The desired taxonomic level at which to run #
#           the analysis. One of: Kingdom, Phylum, Class, Order,   #
#           Family, Genus, Species. If this parameter is not given #
#           the analysis will be performed directly on the counts  #
#           as they are supplied to the script.                    #
#     -t    (Optional) The desired transformation/normalization to #
#           perform on the feature count data. One of:             #
#           CLR (default), log_TSS, or none.                       #
#     -p    (Optional) Pick specific features that will be in the  #
#           analysis. Associations will only be performed for the  #
#           specified features. Must be valid feature names that   #
#           would be returned by the taxa_names() function when ran#
#           on the phyloseq object, or if using the -l flag, valid #
#           taxonomic names of the level in which the phyloseq     #
#           object was collapsed to. Provide feature names in a    #
#           comma separated list with no spaces. Parameters -f and #
#           -u will be ignored if this is specified.               #
#     -k    (Optional) File listing sample names/IIDs to be keep in#
#           the analysis. Samples not included in this list will be#
#           excluded from all pre-processing steps and analysis.   #
#           Each sample name/IID should be listed one ID per line. #
#           Default is to include all samples.                     #
#     -f    (Optional) The minimum proportion of samples a taxon   #
#           must be detected in to be included in the analysis.    #
#           Default is 0.1.                                        #
#     -u    (Optional) Should unclassified taxa be removed?        #
#           Inclusion of this flag will remove taxa unclassified   #
#           at the designated level (from -l flag), while leaving  #
#           this flag out will keep them in. Only applicable when  #
#           a taxonomic level has been designated for -l, will be  #
#           ignored otherwise.                                     #
#     -v    (Optional) Column names of additional variables to be  #
#           included in the analysis. Must be present in the       #
#           sample_data() slot of the phyloseq object. Enter as a  #
#           comma separated list with no spaces.                   #
#     -e    (Optional) Exchange microbiome data for another        #
#           phenotype that has been listed among the variables     #
#           specified with the -v parameter. Listing a variable    #
#           name here will cause that variable to become the       #
#           phenotype and feature count data to be included as     #
#           predictor variables in the GWAS model. Microbial       #
#           features will be referred to as 'FEATURE' for the      #
#           variable name in the results.                          #
#     -c    (Optional) Chromosomes to be included in the           #
#           analysis. Default is to include all autosomes          #
#           (chr 1-22). Input chromosomes as comma-separated       #
#           list and/or ranges with no spaces (e.g. 1-4,22,23,25). #
#     -b    (Optional) Base pair range of SNPs to be included in   #
#           the analysis. Only applicable when '-c' is a single    #
#           chromosome. Should be in the format START-END where    #
#           START is the first base pair position and END is the   #
#           last base pair position. All SNPs included by default. #
#     -a    (Optional) The minor allele frequency treshold that    #
#           SNPs must pass in order to be included in the analysis.#
#           Default is 0.01.                                       #
#     -q    (Optional) Filter SNPs based on an INFO key contained  #
#           in VCF file (e.g. imputation quality score R2). Input  #
#           to this parameter must be in the following form:       #
#           "<key name> <operator> <value>" (e.g. "R2 >= 0.3"),    #
#           and must be a valid key name and value type, else PLINK#
#           will error out. This parameter only applies when       #
#           specifying a directory of VCF files as input with the  #
#           -d parameter. All SNPs that do not satisfy the criteria#
#           will be removed from analysis.                         #
#     -s    (Optional) How should SNP genotypes be coded in the    #
#           model, additive (ADD), dominant (DOM), recessive (REC).#
#           Specify one of: ADD (default), DOM, REC.               #
#     -x    (Optional) Variable to perform interaction with SNPs.  #
#           The default behavior is to only test main effects of   #
#           SNPs and additional covariates (specified with -v), but#
#           if a variable name is supplied here, an interaction    #
#           term between SNP and the variable will be added to the #
#           models. If parameter -e has been given use variable    #
#           name 'FEATURE' when wanting to include microbial       #
#           features in the interaction.                           #
#     -j    (Optional) Perform a 2df joint test of SNP main effect #
#           and its interaction with variable specified in -x.     #
#           Might improve power over pure interaction test if there#
#           is a modest SNP marginal effect.                       #
#     -r    (Optional) Remove results for covariates. By default,  #
#           covariates are not removed from PLINK output files, but#
#           applying this flag will remove them from outputs.      #
#     -R    (Optional) Keep R scripts generated during pre-        #
#           processing of data and plot creation. Might be useful  #
#           to keep for record, or for debugging. If this flag is  #
#           absent, all R scripts are removed at end of analysis.  #
####################################################################

echo " "
echo "####################################################################"
echo "# Microbiome and Genetics Analysis Tool (MaGAT)                    #"
echo "# by Zachary D Wallen                                              #"
echo "# Last updated: 8 June 2021                                        #"
echo "####################################################################"
echo " "

# Argument parsing
while getopts ":hi:g:d:o:l:t:p:k:f:uv:e:c:b:a:q:s:x:jprR" opt; do
  case $opt in
    h)
    echo " Description: This is a program that wraps different R packages   "
    echo " and PLINK 2 to perform association analyses between microbial    "
    echo " feature count data and genome-wide genotypes.                    "
    echo "                                                                  "
    echo " Input data consists of the following:                            "
    echo "                                                                  "
    echo " (1) A phyloseq object (saved as RDS file) that contains a table  "
    echo " with feature counts per sample in the otu_table() slot,          "
    echo " taxonomic assignments (if applicable) in the tax_table() slot,   "
    echo " and any sample metadata that one wants to include in the analysis"
    echo " in the sample_data() slot.                                       "
    echo "                                                                  "
    echo " (2) Genome-wide genotypes in PLINK binary (.bed, .bim, .fam)     "
    echo " format. Sample names of the phyloseq object and IIDs of the PLINK"
    echo " files should contain the same IDs (they do not have to be in the "
    echo " same order) for the program to work properly. Alternatively,     "
    echo " genotypes may be supplied as a directory of VCF files. Sample IDs"
    echo " for VCF files should be in the format 'FID_IID'.                 "
    echo "                                                                  "
    echo " WARNING: Files created during analysis will take up a lot of     "
    echo " storage space if a large number of variants and features are     "
    echo " analyzed, make sure you have plenty of free storage where the    "
    echo " output directory is located.                                     "
    echo "                                                                  "
    echo " WARNING: When using VCFs, VCFs will be converted to PLINK format "
    echo " when running this script. The symbol '_' will be taken as an ID  "
    echo " delimiter to split VCF sample IDs into FIDs and IIDs. Please     "
    echo " ensure VCF sample IDs are in the form 'FID_IID' and that the     "
    echo " sample names of the input phyloseq object matches the IID portion"
    echo " of the VCF ID.                                                   "
    echo "                                                                  "
    echo " WARNING: If using the -c and/or -m filter for SNPs, be sure to   "
    echo " only include VCF files in the directory that contain these SNPs, "
    echo " else script will fail.                                           "
    echo "                                                                  "
    echo " Usage: ./MaGAT.sh -i phyloseq_object.rds \                      "
    echo "                    -g genotype_file_prefix OR -d vcf_directory \ "
    echo "                    -o output_dir \                               "
    echo "                    [additional options]                          "
    echo "                                                                  "
    echo " Parameters:                                                      "
    echo "     -h    Print the parameter list below then exit.              "
    echo "     -i    (Required) Phyloseq object in RDS format, saved from R "
    echo "           with the writeRDS function. Make sure tax_table() slot "
    echo "           is made up of 7 columns (Kingdom, Phylum, Class, Order,"
    echo "           Family, Genus, Species) if wanting to collapse to      "
    echo "           to certain taxonomic levels and any unclassified levels"
    echo "           have 'NA' for the entry.                               "
    echo "     -g    (Required) File name pre-fix for host genome-wide      "
    echo "           genotypes in PLINK binary format (files with extension "
    echo "           .bim, .bam, .fam). Do not include file extensions,     "
    echo "           only pre-fix.                                          "
    echo "     -d    (Required if using VCFs) Directory containing VCF files"
    echo "           with genotype data. Files can be broken up by          "
    echo "           chromosome and gzipped compressed. All files should    "
    echo "           contain the same samples.                              "
    echo "     -o    (Required) Directory to where the output files will be "
    echo "           written to.                                            "
    echo "     -l    (Optional) The desired taxonomic level at which to run "
    echo "           the analysis. One of: Kingdom, Phylum, Class, Order,   "
    echo "           Family, Genus, Species. If this parameter is not given "
    echo "           the analysis will be performed directly on the counts  "
    echo "           as they are supplied to the script.                    "
    echo "     -t    (Optional) The desired transformation/normalization to "
    echo "           perform on the feature count data. One of:             "
    echo "           CLR (default), log_TSS, or none.                       "
    echo "     -p    (Optional) Pick specific features that will be in the  "
    echo "           analysis. Associations will only be performed for the  "
    echo "           specified features. Must be valid feature names that   "
    echo "           would be returned by the taxa_names() function when ran"
    echo "           on the phyloseq object, or if using the -l flag, valid "
    echo "           taxonomic names of the level in which the phyloseq     "
    echo "           object was collapsed to. Provide feature names in a    "
    echo "           comma separated list with no spaces. Parameters -f and "
    echo "           -u will be ignored if this is specified.               "
    echo "     -k    (Optional) File listing sample names/IIDs to be keep in"
    echo "           the analysis. Samples not included in this list will be"
    echo "           excluded from all pre-processing steps and analysis.   "
    echo "           Each sample name/IID should be listed one ID per line. "
    echo "           Default is to include all samples.                     "
    echo "     -f    (Optional) The minimum proportion of samples a taxon   "
    echo "           must be detected in to be included in the analysis.    "
    echo "           Default is 0.1.                                        "
    echo "     -u    (Optional) Should unclassified taxa be removed?        "
    echo "           Inclusion of this flag will remove taxa unclassified   "
    echo "           at the designated level (from -l flag), while leaving  "
    echo "           this flag out will keep them in. Only applicable when  "
    echo "           a taxonomic level has been designated for -l, will be  "
    echo "           ignored otherwise.                                     "
    echo "     -v    (Optional) Column names of additional variables to be  "
    echo "           included in the analysis. Must be present in the       "
    echo "           sample_data() slot of the phyloseq object. Enter as a  "
    echo "           comma separated list with no spaces.                   "
    echo "     -e    (Optional) Exchange microbiome data for another        "
    echo "           phenotype that has been listed among the variables     "
    echo "           specified with the -v parameter. Listing a variable    "
    echo "           name here will cause that variable to become the       "
    echo "           phenotype and feature count data to be included as     "
    echo "           predictor variables in the GWAS model. Microbial       "
    echo "           features will be referred to as 'FEATURE' for the      "
    echo "           variable name in the results.                          "
    echo "     -c    (Optional) Chromosomes to be included in the           "
    echo "           analysis. Default is to include all autosomes          "
    echo "           (chr 1-22). Input chromosomes as comma-separated       "
    echo "           list and/or ranges with no spaces (e.g. 1-4,22,23,25). "
    echo "     -b    (Optional) Base pair range of SNPs to be included in   "
    echo "           the analysis. Only applicable when '-c' is a single    "
    echo "           chromosome. Should be in the format START-END where    "
    echo "           START is the first base pair position and END is the   "
    echo "           last base pair position. All SNPs included by default. "
    echo "     -a    (Optional) The minor allele frequency treshold that    "
    echo "           SNPs must pass in order to be included in the analysis."
    echo "           Default is 0.01.                                       "
    echo "     -q    (Optional) Filter SNPs based on an INFO key contained  "
    echo "           in VCF file (e.g. imputation quality score R2). Input  "
    echo "           to this parameter must be in the following form:       "
    echo "           '<key name> <operator> <value>' (e.g. 'R2 >= 0.3'),    "
    echo "           and must be a valid key name and value type, else PLINK"
    echo "           will error out. This parameter only applies when       "
    echo "           specifying a directory of VCF files as input with the  "
    echo "           -d parameter. All SNPs that do not satisfy the criteria"
    echo "           will be removed from analysis.                         "
    echo "     -s    (Optional) How should SNP genotypes be coded in the    "
    echo "           model, additive (ADD), dominant (DOM), recessive (REC)."
    echo "           Specify one of: ADD (default), DOM, REC.               "
    echo "     -x    (Optional) Variable to perform interaction with SNPs.  "
    echo "           The default behavior is to only test main effects of   "
    echo "           SNPs and additional covariates (specified with -v), but"
    echo "           if a variable name is supplied here, an interaction    "
    echo "           term between SNP and the variable will be added to the "
    echo "           models. If parameter -e has been given use variable    "
    echo "           name 'FEATURE' when wanting to include microbial       "
    echo "           features in the interaction.                           "
    echo "     -j    (Optional) Perform a 2df joint test of SNP main effect "
    echo "           and its interaction with variable specified in -x.     "
    echo "           Might improve power over pure interaction test if there"
    echo "           is a modest SNP marginal effect.                       "
    echo "     -r    (Optional) Remove results for covariates. By default,  "
    echo "           covariates are not removed from PLINK output files, but"
    echo "           applying this flag will remove them from outputs.      "
    echo "     -R    (Optional) Keep R scripts generated during pre-        "
    echo "           processing of data and plot creation. Might be useful  "
    echo "           to keep for record, or for debugging. If this flag is  "
    echo "           absent, all R scripts are removed at end of analysis.  "
    echo " "
    exit 0
    ;;
    i) PHYLO_OBJ="$OPTARG"
    ;;
    g) GENOS="$OPTARG"
    ;;
    d) DOSAGE=$(echo $OPTARG | sed 's#/$##')
    ;;
    o) OUT_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    l) TAXA_LVL="$OPTARG"
    ;;
    t) TRANSF="$OPTARG"
    ;;
    p) FEAT="$OPTARG"
    ;;
    k) KEEP_SAMPS="$OPTARG"
    ;;
    f) FILTER="$OPTARG"
    ;;
    u) UNCLASS=1
    ;;
    v) VARS="$OPTARG"
    ;;
    e) SWAP="$OPTARG"
    ;;
    c) CHR="$OPTARG"
    ;;
    b) RANGE="$OPTARG"
    ;;
    a) MAF="$OPTARG"
    ;;
    q) INFO="$OPTARG"
    ;;
    s) SNP_MOD="$OPTARG"
    ;;
    x) IXN="$OPTARG"
    ;;
    j) JOINT_TEST=1
    ;;
    r) REM_COV=1
    ;;
    R) KEEP_R=1
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
if [[ -z "$PHYLO_OBJ" ]]; then
  echo "ERROR: Argument -i is required, please supply a phyloseq object"
  exit 1
fi
if [[ ! $PHYLO_OBJ =~ ".rds" ]]; then
  echo "ERROR: Extension '.rds' not found for phyloseq object, argument -i requires a phyloseq object saved as an RDS file"
  exit 1
fi

# -g & -d
if [[ -z "$GENOS" ]] && [[ -z "$DOSAGE" ]]; then
  echo "ERROR: Must supply either PLINK genotype file pre-fix (with -g), or directory with VCF files of dosage data (with -d)"
  exit 1
elif [[ ! -z "$GENOS" ]] && [[ ! -z "$DOSAGE" ]]; then
  echo "ERROR: Cannot specify both PLINK genotypes (with -g) and VCF files with dosage data (with -d)"
  exit 1
elif [[ ! -z "$GENOS" ]] && [[ -z "$DOSAGE" ]]; then
  GENO_FILES=$(ls ${GENOS}*)
  if [[ ! $GENO_FILES =~ ".bed" ]] || [[ ! $GENO_FILES =~ ".bim" ]] || [[ ! $GENO_FILES =~ ".fam" ]]; then
    echo "ERROR: Argument -g requires a file name pre-fix associated with .bed, .bim, and .fam files, make sure files have these extensions"
    exit 1
  fi
elif [[ -z "$GENOS" ]] && [[ ! -z "$DOSAGE" ]]; then
  if [[ ! -d "$DOSAGE" ]]; then
    echo "ERROR: Argument -d requires a directory to be listed, please supply a directory with VCF files"
    exit 1
  fi
  DOSE_FILES=$(ls ${DOSAGE}/*)
  if [[ ! $DOSE_FILES =~ ".vcf" ]]; then
    echo "ERROR: Argument -d requires a directory containing VCF files, expecting files with .vcf or .vcf.gz extensions"
    exit 1
  fi
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

# -l
if [[ ! -z "$TAXA_LVL" ]]; then
  ARG_LIST="Kingdom Phylum Class Order Family Genus Species"
  if echo $ARG_LIST | grep -q "$TAXA_LVL"; then
    :
  else
    echo "ERROR: Invalid argument given to -l, please specify one of: Kingdom, Phylum, Class, Order, Family, Genus, Species"
    exit 1
  fi
fi

# -t
if [[ ! -z "$TRANSF" ]]; then
  ARG_LIST="CLR log_TSS none"
  if echo $ARG_LIST | grep -q "$TRANSF"; then
    :
  else
    echo "ERROR: Invalid argument given to -t, please specify one of: CLR, log_TSS, none"
    exit 1
  fi
fi
if [[ -z "$TRANSF" ]]; then
  TRANSF="CLR"
fi

# -p
# Not an easy way to double-check this so let phyloseq give the error if its invalid

# -k
if [[ ! -z "$KEEP_SAMPS" ]]; then
  if [[ ! -f "$KEEP_SAMPS" ]]; then
    echo "ERROR: Argument -k should be a file, please supply an file with sample names/IIDs"
    exit 1
  fi
  if [[ $(awk '{print NF}' $KEEP_SAMPS | uniq | wc -l) -gt 1 ]] || \
     [[ $(awk '{print NF}' $KEEP_SAMPS | uniq) -gt 1 ]]; then
    echo "ERROR: File given to -k should only have one column with one sample name/IID per row"
    exit 1
  fi
fi

# -f
if [[ ! -z "$FILTER" ]]; then
  if [[ $(echo "$FILTER <= 1" | bc) -eq 1 ]]; then
    :
  else
    echo "ERROR: Invalid value given to -f, should be between 0 and 1"
    exit 1
  fi
  if [[ $(echo "$FILTER >= 0" | bc) -eq 1 ]]; then
    :
  else
    echo "ERROR: Invalid value given to -f, should be between 0 and 1"
    exit 1
  fi
fi

# -u
# Not an easy way to double-check this so let phyloseq give the error if its invalid

# -v
if [[ ! -z "$VARS" ]] && [[ $(echo $VARS | awk '{print NF}') > 1 ]]; then
  echo "ERROR: Invalid input give to -v, should be a comma separated list of variable names with no spaces"
  exit 1
fi

# -e
if [[ ! -z "$SWAP" ]]; then
  if [[ -z "$VARS" ]]; then
    echo "ERROR: Covariates must be given with -v in order to use parameter -e"
    exit 1
  fi
  if [[ "$SWAP" == "FEATURE" ]] || (echo $VARS | grep -q "$SWAP"); then
    :
  else
    echo "ERROR: Argument given to -e is not a covariate listed in -v or FEATURE"
    exit 1
  fi
fi

# -c
# Not an easy way to double-check this so let PLINK give the error if its invalid

# -b
if [[ ! -z "$RANGE" ]]; then
  if [[ -z "$CHR" ]]; then
    echo "ERROR: A chromosome must be specified with -c when specifying a base pair range to include in analysis"
    exit 1
  elif echo $CHR | grep -q "," || echo $CHR | grep -q "-"; then
    echo "ERROR: Only one chromosome can be specified when using a base pair range"
    exit 1
  fi
  if echo $RANGE | grep -q " "; then
    echo "ERROR: Do not include spaces when specifying a base pair range"
    exit 1
  fi
  if echo $RANGE | grep -q "-"; then
    :
  else
    echo "ERROR: START and END base pair positions for the base pair position range should be separated by a '-' with no spaces"
    exit 1
  fi
fi

# -a
if [[ ! -z "$MAF" ]]; then
  if [[ $(echo "$MAF <= 0.5" | bc) -eq 1 ]]; then
    :
  else
    echo "ERROR: Invalid value given to -f, should be between 0 and 0.5"
    exit 1
  fi
  if [[ $(echo "$MAF >= 0" | bc) -eq 1 ]]; then
    :
  else
    echo "ERROR: Invalid value given to -f, should be between 0 and 0.5"
    exit 1
  fi
fi

# -q
if [[ ! -z "$INFO" ]]; then
  if [[ -z "$DOSAGE" ]]; then
    echo "ERROR: -q parameter only valid when a directory of VCF files is specified as input using the -d parameter"
    exit 1
  fi
fi

# -s
if [[ ! -z "$SNP_MOD" ]]; then
  ARG_LIST="ADD DOM REC"
  if echo $ARG_LIST | grep -q "$SNP_MOD"; then
    :
  else
    echo "ERROR: Invalid argument given to -s, please specify one of: ADD, DOM, REC"
    exit 1
  fi
fi

# -x
if [[ ! -z "$IXN" ]]; then
  if [[ "$IXN" == "$SWAP" ]]; then
    echo "ERROR: Cannot list same variable for both -x and -e parameters"
    exit 1
  elif [[ ! -z "$SWAP" ]] && [[ "$IXN" == "FEATURE" ]]; then
    :
  elif echo $VARS | grep -q "$IXN"; then
    :
  else
    echo "ERROR: Invalid input give to -x, should be a variable name that is included in the variable list supplied to the -v parameter"
    exit 1
  fi
fi

# -j
if [[ ! -z "$JOINT_TEST" ]]; then
  if [[ ! -z "$JOINT_TEST" ]] && [[ -z "$IXN" ]]; then
    echo "ERROR: 2df joint test cannot be specified if a variable name has not been supplied to perform an interaction with (using the -x parameter)."
  else
    :
  fi
fi

# If directory with VCF files given, make temporary .fam file to use for phenotype and covariate file creation
if [[ ! -z "$DOSAGE" ]]; then
  echo " "
  echo "*** VCF files given instead of PLINK binary files, making temporary .fam file from smallest sized VCF file in supplied directory to use in downstream processes ***"
  echo " "
  if ls $OUT_DIR | grep -q "tEmPoRaRy"; then
    rm ${OUT_DIR}/*tEmPoRaRy*
  fi
  file=$(ls -l ${DOSAGE}/*vcf* | sort -g -k5,5 | awk 'NR == 1{print $NF}')
  plink2 --vcf $file --id-delim --make-just-fam --out ${OUT_DIR}/tEmPoRaRy.plink
  FAM_FILE="${OUT_DIR}/tEmPoRaRy.plink.fam"
fi

############# START PRE-PROCESS OF PHYLOSEQ DATA #############
echo " "
echo "*** Parameters for phyloseq data pre-processing ***"
echo " "

# Initiate R script for data pre-processing
echo "### Pre-processing of phyloseq data ###" > ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "suppressMessages(library(phyloseq))" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "# Read in phyloseq object" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "ps <- readRDS('${PHYLO_OBJ}')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "# Make sure phyloseq object is in sample x feature orientation" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "if (taxa_are_rows(ps)){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "    ps <- phyloseq(t(otu_table(ps)), sample_data(ps), tax_table(ps))" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "cat('\n','Summary of input phyloseq object:','\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "ps" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R

# Add syntax for collapsing taxa to desired taxonomic level
if [[ -z "$TAXA_LVL" ]]; then
  echo "- Taxa level for analysis: input level"
  echo " "
else
  echo "- Taxa level for analysis: $TAXA_LVL"
  echo " "
  echo "# Collapse phyloseq object to ${TAXA_LVL} level" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps <- tax_glom(ps, taxrank = '${TAXA_LVL}', NArm=F)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "cat('\n','Summary after collapsing phyloseq object to ${TAXA_LVL} level:','\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
fi

# Add syntax to replace taxa_names() with more meaningful labels if phyloseq object was collapsed
if [[ ! -z "$TAXA_LVL" ]]; then
  if [[ "$TAXA_LVL" == "Species" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "for (i in 1:length(taxa_names(ps))){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:7])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  if (length(taxonomy.sub) == 7){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('s__',taxonomy.sub[7],sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 6){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('g__',taxonomy.sub[6],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 5){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('f__',taxonomy.sub[5],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 4){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('o__',taxonomy.sub[4],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 3){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else{" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  elif [[ "$TAXA_LVL" == "Genus" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "for (i in 1:length(taxa_names(ps))){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:6])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  if (length(taxonomy.sub) == 6){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('g__',taxonomy.sub[6],sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 5){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('f__',taxonomy.sub[5],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 4){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('o__',taxonomy.sub[4],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 3){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else{" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  elif [[ "$TAXA_LVL" == "Family" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "for (i in 1:length(taxa_names(ps))){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:5])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  if (length(taxonomy.sub) == 5){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('f__',taxonomy.sub[5],sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 4){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('o__',taxonomy.sub[4],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 3){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else{" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  elif [[ "$TAXA_LVL" == "Order" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "for (i in 1:length(taxa_names(ps))){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:4])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  if (length(taxonomy.sub) == 4){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('o__',taxonomy.sub[4],sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 3){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else{" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  elif [[ "$TAXA_LVL" == "Class" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "for (i in 1:length(taxa_names(ps))){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:3])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  if (length(taxonomy.sub) == 3){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else{" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  elif [[ "$TAXA_LVL" == "Phylum" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "for (i in 1:length(taxa_names(ps))){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:2])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  if (length(taxonomy.sub) == 2){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else{" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  elif [[ "$TAXA_LVL" == "Kingdom" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "for (i in 1:length(taxa_names(ps))){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  if (length(taxonomy.sub) == 1){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],sep='')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }else{" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "  }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  fi
fi

# Add syntax to extract certain samples if requested
if [[ ! -z "$KEEP_SAMPS" ]]; then
  echo "# Extract specific samples to use in pre-processing and analysis" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "samp.list <- read.table('$KEEP_SAMPS')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps <- prune_samples(samp.list[,1], ps)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps <- filter_taxa(prune_samples(samp.list[,1], ps), function(x){sum(x>0)>0}, TRUE)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "cat('\n','Summary after removing samples not found int ${KEEP_SAMPS} file:','\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
fi

# Add syntax for feature data transformation
echo "- Transformation for feature count data: ${TRANSF}"
echo " "
echo "# Perform normalization/transformation on feature count data" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
if [[ "$TRANSF" = "CLR" ]]; then
  echo "ps.t <- transform_sample_counts(ps, function(x){log(x+1)-mean(log(x+1))})" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
elif [[ "$TRANSF" = "log_TSS" ]]; then
  echo "log.trans <- function(x) {" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "  y <- replace(x, x == 0, min(x[x>0]) / 2)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "  return(log(y))" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps.t <- transform_sample_counts(ps, function(x){x/sum(x)})" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "otu_table(ps.t) <- otu_table(apply(otu_table(ps.t), 2, log.trans), taxa_are_rows=F)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
elif [[ "$TRANSF" = "none" ]]; then
  echo "ps.t <- ps" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
fi

# Add syntax for subsetting phyloseq object for specific feature specified by user
if [[ ! -z "$FEAT" ]]; then
  echo "- Subsetting microbiome data for specific feature: $FEAT"
  echo " "
  echo "# Subset microbiome data for specified feature: $FEAT" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "target.feat <- strsplit('${FEAT}', ',')[[1]]" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  if [[ ! -z "$TAXA_LVL" ]]; then
    echo "ps.t <- subset_taxa(ps.t, ${TAXA_LVL} %in% target.feat)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  else
    echo "ps.t <- prune_taxa(target.feat, ps.t)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  fi
    echo "cat('\n','Summary after subsetting data for requested feature(s):', '\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "ps.t" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
fi

# Add syntax for filtering microbiome count data
if [[ ! -z "$FEAT" ]]; then
  :
elif [[ ! -z "$FILTER" ]]; then
  echo "- Proportion of samples feature must be detected in to be included in this analysis: $FILTER"
  echo " "
  echo "# Filter out features that were detected below minimum proportion of samples equal to ${FILTER}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "filt_feat <- taxa_names(filter_taxa(ps, function(x){sum(x > 0) >= (${FILTER}*length(x))}, TRUE))" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps.t <- subset_taxa(ps.t, taxa_names(ps.t) %in% filt_feat)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "cat('\n','Summary after filtering out features that were detected below minimum proportion of samples equal to ${FILTER}:', '\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps.t" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
else
  FILTER=0.1
  echo "- Proportion of samples feature must be detected in to be included in this analysis: $FILTER"
  echo " "
  echo "# Filter out features that were detected below minimum proportion of samples equal to ${FILTER}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "filt_feat <- taxa_names(filter_taxa(ps, function(x){sum(x > 0) >= (${FILTER}*length(x))}, TRUE))" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps.t <- subset_taxa(ps.t, taxa_names(ps.t) %in% filt_feat)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "cat('\n','Summary after filtering out features that were detected below minimum proportion of samples equal to ${FILTER}:', '\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "ps.t" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
fi

# Add syntax for removing unclassified taxa
if [[ ! -z "$TAXA_LVL" ]]; then
  if [[ ! -z "$FEAT" ]]; then
    :
  elif [[ ! -z "$UNCLASS" ]]; then
    echo "- Removing unclassified features at the $TAXA_LVL level"
    echo " "
    echo "# Remove unclassifed taxa at the ${TAXA_LVL}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "ps.t <- subset_taxa(ps.t, "'!'"is.na(${TAXA_LVL}))" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "cat('\n','Summary after removing unclassified taxa at the ${TAXA_LVL} level:', '\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo "ps.t" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
    echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  else
    echo "- Keeping unclassified features at the $TAXA_LVL level in the analysis"
    echo " "
  fi
fi

# Add syntax to create microbiome phenotype file for input to PLINK
if [[ ! -z "$GENOS" ]]; then
  echo "# Read in PLINK fam file" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "fam.file <- read.table('${GENOS}.fam', comment.char='', stringsAsFactors=F)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
elif [[ ! -z "$DOSAGE" ]]; then
  echo "# Read in PLINK fam file" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "fam.file <- read.table('${FAM_FILE}', comment.char='', stringsAsFactors=F)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
fi
echo "cat('\n','Number of samples found in genotype files:', nrow(fam.file), '\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "# Extract microbiome data from phyloseq object" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "feat.df <- data.frame(otu_table(ps.t))" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "cat('\n','Number of samples found in microbiome data:', nrow(feat.df), '\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "# Find overlapping samples between microbome and genotype data" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "fam.file.filt <- fam.file[fam.file[,2] %in% rownames(feat.df),]" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "feat.df.filt <- feat.df[rownames(feat.df) %in% fam.file[,2],,drop=F]" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "if (nrow(feat.df.filt) == 0){ stop('No overlapping samples found between microbiome and genotype data. Are sample names of phyloseq object and IIDs of genotype data concordant?')}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "cat('\n','Number of samples that overlap between phyloseq and genotype data and will be included in phenotype/covariate files:', nrow(feat.df.filt), '\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "# Order microbiome and genotype samples the same so the correct FID and IIDs are added to the microbiome phenotype file" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "feat.df.filt <- feat.df.filt[order(rownames(feat.df.filt)),,drop=F]" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "fam.file.filt <- fam.file.filt[order(fam.file.filt[,2]),]" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo 'if (!identical(rownames(feat.df.filt), fam.file.filt[,2])){ stop("Sample IDs between microbiome and genotype data do not match, even after finding overlaps and ordering the same.")}' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "# Create microbiome phenotype file" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "pheno.file <- data.frame(FID=fam.file.filt[,1], IID=fam.file.filt[,2], feat.df.filt)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo "write.table(pheno.file, '${OUT_DIR}/phenotype_file.txt', row.names=F, quote=F, sep='\t')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R

# Add syntax to create covariate file with requested variables
if [[ -z "$VARS" ]]; then
  echo "- No additional variables requested to be included in analysis. Linear models will only include the SNP variable as predictor"
  echo " "
else
  echo "- Additional variables to be included in analysis: $(echo $VARS | sed 's/,/ /g')"
  echo " "
  echo "# Check to make sure additional variables provided are in sample data" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "vars <- strsplit('${VARS}', ',')[[1]]" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo 'if (!sum(vars %in% colnames(sample_data(ps.t))) == length(vars)){ stop("Additional variables provided to be included in model were not all found in the phyloseq object sample data.")}' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "# Extract desired variables from phyloseq object sample data" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "samp.df <- data.frame(sample_data(ps.t)[,vars])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "cat('\n','Additional variables requested were found in phyloseq object sample data and will be included in covariate file:', colnames(samp.df), '\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "# Find overlapping samples between microbome and genotype data" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "samp.df.filt <- samp.df[rownames(samp.df) %in% fam.file.filt[,2],,drop=F]" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "# Order sample data the same as genotype data so the correct FID and IIDs are added to the covariate file" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "samp.df.filt <- samp.df.filt[order(rownames(samp.df.filt)),,drop=F]" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo 'if (!identical(rownames(samp.df.filt), fam.file.filt[,2])){ stop("Sample IDs between sample data and genotype data do not match, even after finding overlaps and ordering the same.")}' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "# Create covariate file with additional variables to be included in models" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "covar.file <- data.frame(FID=fam.file.filt[,1], IID=fam.file.filt[,2], samp.df.filt)" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "# Dummy-code any categorical variables so PLINK works correctly" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "for (i in 3:ncol(covar.file)){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo '  if (!is.numeric(covar.file[,i]) | length(table(covar.file[,i]))<=2){' >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "    cat('\n','WARNING: While creating covariate file for PLINK, variable [',colnames(covar.file)[i],'] was detected as being categorical. It will be dummy-coded to 2 and 1 for analysis. Please check covariate_file.txt to ensure this was done correctly.','\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "    levels <- names(table(covar.file[,i]))" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "    if (length(levels) == 2){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "      covar.file[,i] <- gsub(levels[1], '2', covar.file[,i])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "      covar.file[,i] <- gsub(levels[2], '1', covar.file[,i])" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "    }else if (length(levels) > 2){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "      stop('Variable [',colnames(covar.file)[i],'] was detected as a categorical variable with >2 levels. Automatic dummy-coding of categorical variables with >2 levels is unsupported at this time. Please recode the variables to numeric values manually and try again.')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "    }else if (length(levels) < 2){" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "      stop('Please check variable [',colnames(covar.file)[i],']. It was detected as categorical with <2 levels. PLINK will not perform correctly if given a variable with only 1 level.')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "    }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "  }" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "}" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo "write.table(covar.file, '${OUT_DIR}/covariate_file.txt', row.names=F, quote=F, sep='\t')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
  echo " " >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R
fi

# Signal end of pre-processing
echo "cat('\n', 'Pre-processing of phyloseq data, and creation of phenotype/covariate files complete.', '\n')" >> ${OUT_DIR}/Pre-Process_Phyloseq_Data.R

# Run Pre-Process_Phyloseq_Data.R script to perform the data pre-processing and generation of phenotype and covariate file for PLINK
echo "*** Performing phyloseq data pre-processing ***"
Rscript ${OUT_DIR}/Pre-Process_Phyloseq_Data.R

# If parameter given to make a covariate as phenotype, create new phenotype/covariate files with covariate as phenotype and add microbiome features as covariates
if [[ ! -z "$SWAP" ]]; then
  echo "pheno.file <- read.table('${OUT_DIR}/phenotype_file.txt', header=T, stringsAsFactors=F, comment.char='')" > ${OUT_DIR}/Swap_phenotype.R
  echo "covar.file <- read.table('${OUT_DIR}/covariate_file.txt', header=T, stringsAsFactors=F, comment.char='')" >> ${OUT_DIR}/Swap_phenotype.R
  echo "if (!identical(pheno.file[,c('FID','IID')], covar.file[,c('FID','IID')])){ stop('ERROR: FID and IID in phenotype_file.txt and covariate_file.txt do not match, cannot make phenotype swap.')}" >> ${OUT_DIR}/Swap_phenotype.R
  echo "pheno.file.new <- cbind(pheno.file[,c('FID','IID')], covar.file[,'${SWAP}', drop=F])" >> ${OUT_DIR}/Swap_phenotype.R
  echo "write.table(pheno.file.new, '${OUT_DIR}/phenotype_file.txt', row.names=F, quote=F, sep='\t')" >> ${OUT_DIR}/Swap_phenotype.R
  echo "covar.file.new <- cbind(pheno.file, covar.file[,which(!colnames(covar.file) %in% c('FID','IID','${SWAP}'))])" >> ${OUT_DIR}/Swap_phenotype.R
  echo "write.table(covar.file.new, '${OUT_DIR}/covariate_file.txt', row.names=F, quote=F, sep='\t')" >> ${OUT_DIR}/Swap_phenotype.R
  echo "for (i in 3:ncol(pheno.file)){" >> ${OUT_DIR}/Swap_phenotype.R
  echo "  covar.file.new <- cbind(pheno.file[,c(1,2,i)], covar.file[,which(!colnames(covar.file) %in% c('FID','IID','${SWAP}'))])" >> ${OUT_DIR}/Swap_phenotype.R
  echo "  colnames(covar.file.new)[3] <- 'FEATURE'" >> ${OUT_DIR}/Swap_phenotype.R
  echo "  write.table(covar.file.new, paste('${OUT_DIR}/',colnames(pheno.file)[i],'_cov_file_tEmPoRaRy.txt',sep=''), row.names=F, quote=F, sep='\t')" >> ${OUT_DIR}/Swap_phenotype.R
  echo "}" >> ${OUT_DIR}/Swap_phenotype.R
  
  echo " "
  echo "*** Swapping phenotypes, variable ${SWAP} will be used as phenotype while microbiome data will be included as predictors for analyses ***"
  Rscript ${OUT_DIR}/Swap_phenotype.R
fi

# Grab phenotypes for use later
PHENOS=$(awk 'NR == 1{$1=$2=""; print $0}' ${OUT_DIR}/phenotype_file.txt)

# Remove any previous results if they exist
for pheno in $PHENOS
do
  if [[ ! -z "$(ls ${OUT_DIR}/*${pheno}.glm* 2>/dev/null)" ]]; then
      rm ${OUT_DIR}/*${pheno}.glm*
  fi
done

############# START ASSOCIATION ANALYSIS #############

echo " "
echo "*** Performing genome-wide association analyses with each feature ***"
echo " "

# Apply default chromosome option if one was not supplied
if [[ -z "$CHR" ]]; then
  CHR=1-22
fi
echo "- Running analysis for chromosomes $CHR"
echo " "

# Get start and end base pair positions
if [[ ! -z "$RANGE" ]]; then
  echo "- Running analysis for base pair range $RANGE"
  echo " "
  START_BP=$(echo $RANGE | awk -F"-" '{print $1}')
  END_BP=$(echo $RANGE | awk -F"-" '{print $2}')
fi

# Apply default MAF threshold option if one was not supplied
if [[ -z "$MAF" ]]; then
  MAF=0.01
fi
echo "- MAF threshold SNPs must pass to be included in analysis is $MAF"
echo " "

# State what info field is being filtered by if given
if [[ ! -z "$INFO" ]]; then
  echo "- SNPs with $INFO will be included in the analysis"
  echo " "
fi

# Create parameter input for interaction if specified
if [[ ! -z "$IXN" ]]; then
  IXN_PARAM="interaction"
  if [[ -z "$SWAP" ]]; then
    N_VAR=$(echo $(($(awk '{print NF}' ${OUT_DIR}/covariate_file.txt | sort -nu | tail -n 1)-1)))
    VAR_N=$(($(awk -v RS='\t' "/$IXN/{print NR}" ${OUT_DIR}/covariate_file.txt)-1))
    IXN_N=$((($N_VAR-1) + $VAR_N))
    PARAM=$(echo $(seq $N_VAR) $IXN_N)
    echo "- Adding interaction term between SNP and $IXN in the linear models"
    echo " "
  elif [[ ! -z "$SWAP" ]]; then
    file=$(ls -l ${OUT_DIR}/*cov_file_tEmPoRaRy* | awk 'NR == 1{print $NF}')
    N_VAR=$(echo $(($(awk '{print NF}' $file | sort -nu | tail -n 1)-1)))
    VAR_N=$(($(awk -v RS='\t' "/$IXN/{print NR}" $file)-1))
    IXN_N=$((($N_VAR-1) + $VAR_N))
    PARAM=$(echo $(seq $N_VAR) $IXN_N)
    echo "- Adding interaction term between SNP and $IXN in the linear models"
    echo " "
  fi
fi

# Apply default SNP model option if one was not supplied
if [[ -z "$SNP_MOD" ]]; then
  SNP_MOD=ADD
fi

# Get PLINK parameter for SNP model
if [[ "$SNP_MOD" = "REC" ]]; then
  SNP_PARAM="recessive"
elif [[ "$SNP_MOD" = "DOM" ]]; then
  SNP_PARAM="dominant"
fi

# Return the model to be run in PLINK
echo "- Running analysis using PLINK with following linear model:"
if [[ -z "$SWAP" ]]; then
  if [[ ! -z "$VARS" ]] && [[ -z "$IXN" ]]; then
    echo "  Feature~SNP(${SNP_MOD})+$(head -n1  ${OUT_DIR}/covariate_file.txt | cut -f 3- | sed $'s/\t/+/g')"
    echo " "
  elif [[ ! -z "$VARS" ]] && [[ ! -z "$IXN" ]]; then
    echo "  Feature~SNP(${SNP_MOD})+$(head -n1  ${OUT_DIR}/covariate_file.txt | cut -f 3- | sed $'s/\t/+/g')+SNP:${IXN}"
    echo " "
  else
    echo "  Feature~SNP(${SNP_MOD})"
    echo " "
  fi
fi
if [[ ! -z "$SWAP" ]]; then
  file=$(ls -l ${OUT_DIR}/*cov_file_tEmPoRaRy* | awk 'NR == 1{print $NF}')
  if [[ ! -z "$VARS" ]] && [[ -z "$IXN" ]]; then
    echo "  $(head -n1 ${OUT_DIR}/phenotype_file.txt | cut -f 3-)~SNP(${SNP_MOD})+$(head -n1 $file | cut -f 3- | sed $'s/\t/+/g')"
    echo " "
  elif [[ ! -z "$VARS" ]] && [[ ! -z "$IXN" ]]; then
    echo "  $(head -n1 ${OUT_DIR}/phenotype_file.txt | cut -f 3-)~SNP(${SNP_MOD})+$(head -n1 $file | cut -f 3- | sed $'s/\t/+/g')+SNP:${IXN}"
    echo " "
  fi
fi

# Create test input for 2df joint test if specified to do so
if [[ ! -z "$JOINT_TEST" ]]; then
  TESTS="1 $(echo $PARAM | awk '{print NF}')"
  echo "- Performing 2df joint test for SNP main and SNP interaction effect between following models:"
  if [[ -z "$SWAP" ]]; then
    echo "  Feature~SNP(${SNP_MOD})+$(head -n1  ${OUT_DIR}/covariate_file.txt | cut -f 3- | sed $'s/\t/+/g')+SNP:${IXN}"
    echo "  Feature~$(head -n1  ${OUT_DIR}/covariate_file.txt | cut -f 3- | sed $'s/\t/+/g')"
    echo " "
  elif [[ ! -z "$SWAP" ]]; then
    file=$(ls -l ${OUT_DIR}/*cov_file_tEmPoRaRy* | awk 'NR == 1{print $NF}')
    echo "  $(head -n1 ${OUT_DIR}/phenotype_file.txt | cut -f 3-)~SNP(${SNP_MOD})+$(head -n1 $file | cut -f 3- | sed $'s/\t/+/g')+SNP:${IXN}"
    echo "  $(head -n1 ${OUT_DIR}/phenotype_file.txt | cut -f 3-)~$(head -n1 $file | cut -f 3- | sed $'s/\t/+/g')"
    echo " "
  fi
fi

# Run GWAS with microbiome data as phenotype
if [[ -z $SWAP ]]; then

  # If variables supplied, grab names of quantitative one for standardization during analysis
  if [[ ! -z "$VARS" ]]; then
    covars=$(sed -n 1p ${OUT_DIR}/covariate_file.txt | cut -f 3-)
    for i in $covars;
    do
      levels=$(cat ${OUT_DIR}/covariate_file.txt | \
      awk -v COVAR="$i" \
        'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} NR>1 {print $(f[COVAR])}' | \
      sort | uniq | wc -l)
      if [[ "$levels" -gt 2 ]]; then
        quant_covars=$(echo "$quant_covars $i")
      fi
    done
    echo " - Following variables were detected as quantitative and will be standardized in analysis:"
    echo "  $quant_covars"
    echo " "
  fi
  
  # Build PLINK command to run
    echo "plink2 \\" > ${OUT_DIR}/Run_PLINK.sh
  if [[ ! -z "$GENOS" ]]; then
    echo "--bfile $GENOS \\" >> ${OUT_DIR}/Run_PLINK.sh
  fi
  if [[ ! -z "$DOSAGE" ]]; then
    echo "--vcf \${1} dosage=DS \\" >> ${OUT_DIR}/Run_PLINK.sh
    echo "--id-delim _ \\" >> ${OUT_DIR}/Run_PLINK.sh
  fi
    echo "--pheno ${OUT_DIR}/phenotype_file.txt \\" >> ${OUT_DIR}/Run_PLINK.sh
  if [[ ! -z "$VARS" ]]; then
    echo "--covar ${OUT_DIR}/covariate_file.txt \\" >> ${OUT_DIR}/Run_PLINK.sh
    echo "--variance-standardize $quant_covars \\" >> ${OUT_DIR}/Run_PLINK.sh
  fi
    echo "--chr $CHR \\" >> ${OUT_DIR}/Run_PLINK.sh
  if [[ ! -z "$RANGE" ]]; then
    echo "--from-bp $START_BP \\" >> ${OUT_DIR}/Run_PLINK.sh
    echo "--to-bp $END_BP \\" >> ${OUT_DIR}/Run_PLINK.sh
  fi
    echo "--maf $MAF \\" >> ${OUT_DIR}/Run_PLINK.sh
  if [[ ! -z "$DOSAGE" ]] && [[ ! -z "$INFO" ]]; then
    echo "--extract-if-info $INFO \\" >> ${OUT_DIR}/Run_PLINK.sh
  fi
    echo "--glm $SNP_PARAM $IXN_PARAM \\" >> ${OUT_DIR}/Run_PLINK.sh
    echo "--ci 0.95 \\" >> ${OUT_DIR}/Run_PLINK.sh
  if [[ ! -z "$IXN" ]]; then
    echo "--parameters $PARAM \\" >> ${OUT_DIR}/Run_PLINK.sh
  fi
  if [[ ! -z "$JOINT_TEST" ]]; then
    echo "--tests $TESTS \\"
  fi
  if [[ ! -z "$GENOS" ]]; then
    echo "--out ${OUT_DIR}/" >> ${OUT_DIR}/Run_PLINK.sh
  fi
  if [[ ! -z "$DOSAGE" ]]; then
    echo "--out ${OUT_DIR}/\${2}" >> ${OUT_DIR}/Run_PLINK.sh
  fi
  chmod +x ${OUT_DIR}/Run_PLINK.sh
  
  # Run PLINK command
  if [[ ! -z "$GENOS" ]]; then
    ./${OUT_DIR}/Run_PLINK.sh
    for pheno in $PHENOS
    do
      assoc_file=$(ls ${OUT_DIR}/${pheno}.glm* | awk -F'/' '{print $NF}')
      sed '1s/#//' ${OUT_DIR}/${assoc_file} | sed '/^#/d' > ${OUT_DIR}/${assoc_file}.tEmPoRaRy
      rm ${OUT_DIR}/${assoc_file}
      mv ${OUT_DIR}/${assoc_file}.tEmPoRaRy ${OUT_DIR}/${assoc_file}
      rm ${OUT_DIR}/${pheno}*log
    done
  elif [[ ! -z "$DOSAGE" ]]; then
    for DOSE_FILE in ${DOSAGE}/*vcf*
    do
      echo " "
      echo " "
      OUT_FILE=$(echo $DOSE_FILE | awk -F'.vcf' '{print $1}' | awk -F'/' '{print $NF}')
      ./${OUT_DIR}/Run_PLINK.sh $DOSE_FILE $OUT_FILE
      rm ${OUT_DIR}/${OUT_FILE}.log
    done
    echo " "
    echo "Merging result files of the same feature..."
    for pheno in $PHENOS
    do
      assoc=$(ls -l ${OUT_DIR}/*${pheno}.glm* | grep ${pheno}.glm | head -n1 | awk '{print $NF}' | awk -F'.' '{print $NF}')
      cat ${OUT_DIR}/*${pheno}.glm.${assoc} | sed '1s/#//' | sed '/^#/d' > ${OUT_DIR}/${pheno}.glm.${assoc}.tEmPoRaRy
      rm ${OUT_DIR}/*${pheno}.glm.${assoc}
      mv ${OUT_DIR}/${pheno}.glm.${assoc}.tEmPoRaRy ${OUT_DIR}/${pheno}.glm.${assoc}
      echo "$pheno done."
    done
  fi
fi

# Run GWAS for swapped phenotype if applicable
if [[ ! -z $SWAP ]]; then
  for cov_file in ${OUT_DIR}/*cov_file_tEmPoRaRy*
  do
    # Get feature name that is now a covariate
    feature=$(echo $cov_file | awk -F'_cov_file_tEmPoRaRy' '{print $1}' | awk -F'/' '{print $NF}')

    # Grab names of quantitative one for standardization during analysis
    covars=$(sed -n 1p $cov_file | cut -f 3-)
    for i in $covars;
    do
      levels=$(cat $cov_file | \
      awk -v COVAR="$i" \
        'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} NR>1 {print $(f[COVAR])}' | \
      sort | uniq | wc -l)
      if [[ "$levels" -gt 2 ]]; then
        quant_covars=$(echo "$quant_covars $i")
      fi
    done
    echo " - Following variables were detected as quantitative and will be standardized in analysis:"
    echo "  $quant_covars"
    echo " "

    # Build PLINK command to run
      echo "plink2 \\" > ${OUT_DIR}/Run_PLINK.sh
    if [[ ! -z "$GENOS" ]]; then
      echo "--bfile $GENOS \\" >> ${OUT_DIR}/Run_PLINK.sh
    fi
    if [[ ! -z "$DOSAGE" ]]; then
      echo "--vcf \${1} dosage=DS \\" >> ${OUT_DIR}/Run_PLINK.sh
      echo "--id-delim _ \\" >> ${OUT_DIR}/Run_PLINK.sh
    fi
      echo "--pheno ${OUT_DIR}/phenotype_file.txt \\" >> ${OUT_DIR}/Run_PLINK.sh
      echo "--covar ${cov_file} \\" >> ${OUT_DIR}/Run_PLINK.sh
      echo "--variance-standardize $quant_covars \\" >> ${OUT_DIR}/Run_PLINK.sh
      echo "--chr $CHR \\" >> ${OUT_DIR}/Run_PLINK.sh
    if [[ ! -z "$RANGE" ]]; then
      echo "--from-bp $START_BP \\" >> ${OUT_DIR}/Run_PLINK.sh
      echo "--to-bp $END_BP \\" >> ${OUT_DIR}/Run_PLINK.sh
    fi
      echo "--maf $MAF \\" >> ${OUT_DIR}/Run_PLINK.sh
    if [[ ! -z "$DOSAGE" ]] && [[ ! -z "$INFO" ]]; then
      echo "--extract-if-info $INFO \\" >> ${OUT_DIR}/Run_PLINK.sh
    fi
      echo "--glm $SNP_PARAM $IXN_PARAM \\" >> ${OUT_DIR}/Run_PLINK.sh
      echo "--ci 0.95 \\" >> ${OUT_DIR}/Run_PLINK.sh
    if [[ ! -z "$IXN" ]]; then
      echo "--parameters $PARAM \\" >> ${OUT_DIR}/Run_PLINK.sh
    fi
   if [[ ! -z "$JOINT_TEST" ]]; then
      echo "--tests $TESTS \\"
    fi
    if [[ ! -z "$GENOS" ]]; then
      echo "--out ${OUT_DIR}/${feature}" >> ${OUT_DIR}/Run_PLINK.sh
    fi
    if [[ ! -z "$DOSAGE" ]]; then
      echo "--out ${OUT_DIR}/\${2}.${feature}" >> ${OUT_DIR}/Run_PLINK.sh
    fi
    chmod +x ${OUT_DIR}/Run_PLINK.sh
  
    # Run PLINK command
    if [[ ! -z "$GENOS" ]]; then
      ./${OUT_DIR}/Run_PLINK.sh
      rm ${OUT_DIR}/${feature}.log
      pattern=$(echo ${feature}.${PHENOS} | sed 's/ //')
      assoc_file=$(ls ${OUT_DIR}/${pattern}* | awk -F'/' '{print $NF}')
      sed '1s/#//' ${OUT_DIR}/${assoc_file} | sed '/^#/d' > ${OUT_DIR}/${assoc_file}.tEmPoRaRy
      rm ${OUT_DIR}/${assoc_file}
      mv ${OUT_DIR}/${assoc_file}.tEmPoRaRy ${OUT_DIR}/${assoc_file}
    elif [[ ! -z "$DOSAGE" ]]; then
      for DOSE_FILE in ${DOSAGE}/*vcf*
      do
        echo " "
        echo " "
        OUT_FILE=$(echo $DOSE_FILE | awk -F'.vcf' '{print $1}' | awk -F'/' '{print $NF}')
        ./${OUT_DIR}/Run_PLINK.sh $DOSE_FILE $OUT_FILE
        rm ${OUT_DIR}/${OUT_FILE}.${feature}.log
      done
      pattern=$(echo ${feature}.${PHENOS} | sed 's/ //')
      assoc=$(ls -l ${OUT_DIR}/*${pattern}.glm* | grep $pattern.glm | head -n1 | awk '{print $NF}' | awk -F'.' '{print $NF}')
      cat ${OUT_DIR}/*${pattern}.glm.${assoc} | sed '1s/#//' | sed '/^#/d' > ${OUT_DIR}/${pattern}.glm.${assoc}.tEmPoRaRy
      rm ${OUT_DIR}/*${pattern}.glm.${assoc}
      mv ${OUT_DIR}/${pattern}.glm.${assoc}.tEmPoRaRy ${OUT_DIR}/${pattern}.glm.${assoc}
    fi
    echo " "
    echo "GWAS for ${feature} done"
  done
fi

echo " "
echo "PLINK run complete."
echo " "

# Remove results for covariates if specified
if [[ ! -z $REM_COV ]]; then
  if [[ ! -z "$IXN" ]] && [[ -z "$JOINT_TEST" ]]; then
    echo " "
    echo "Subsetting PLINK result files for SNP(${SNP_MOD})x${IXN} variable results..."
    echo " "
    for pheno in $PHENOS
    do
      FILE=$(ls ${OUT_DIR}/*${pheno}.glm*)
      cat $FILE | \
      awk -v IXN="$IXN" -v SNP_MOD="$SNP_MOD" \
        'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || $(f["TEST"])==SNP_MOD"x"IXN)' > ${FILE}.mod
      mv ${FILE}.mod $FILE
    done
  elif [[ ! -z "$IXN" ]] && [[ ! -z "$JOINT_TEST" ]]; then
    echo " "
    echo "Subsetting PLINK result files for 2df joint test results..."
    echo " "
    for pheno in $PHENOS
    do
      FILE=$(ls ${OUT_DIR}/*${pheno}.glm*)
      cat $FILE | \
      awk -v IXN="$IXN" -v SNP_MOD="$SNP_MOD" \
        'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || $(f["TEST"])==SNP_MOD || $(f["TEST"])==SNP_MOD"x"IXN || $(f["TEST"])=="USER_2DF")' > ${FILE}.mod
      mv ${FILE}.mod $FILE
    done
  else
    echo " "
    echo "Subsetting PLINK result files for SNP(${SNP_MOD}) variable results..."
    echo " "
    for pheno in $PHENOS
    do
      FILE=$(ls ${OUT_DIR}/*${pheno}.glm*)
      cat $FILE | \
      awk -v SNP_MOD="$SNP_MOD" \
        'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || $(f["TEST"])==SNP_MOD)' > ${FILE}.mod
      mv ${FILE}.mod $FILE
    done
  fi
else
  echo " "
  echo "Removal of covariates from result file(s) not requested, results for all variables will be kept."
  echo " "
fi

# Clean up Rscripts
if [[ -z "$KEEP_R" ]]; then
  rm ${OUT_DIR}/*.R
elif [[ ! -z "$KEEP_R" ]]; then
  if [[ ! -d "${OUT_DIR}/helper_R_scripts" ]]; then
    mkdir ${OUT_DIR}/helper_R_scripts
  fi
  mv ${OUT_DIR}/*.R ${OUT_DIR}/helper_R_scripts/
fi

# Gzip files
echo " "
echo "Compressing PLINK result files..."
echo " "
for pheno in $PHENOS
do
  gzip -f ${OUT_DIR}/*${pheno}.glm*
done

# Remove any temporary files
rm ${OUT_DIR}/*tEmPoRaRy*

echo " "
echo "*** Running of MaGAT has completed ***"
echo " "
