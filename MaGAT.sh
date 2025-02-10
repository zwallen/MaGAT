#!/bin/bash
set -e

####################################################################
# Microbiome and Genetics Analysis Tool (MaGAT) v 1.0.0            #
# by Zachary D Wallen                                              #
# Last updated: 9 Feb 2025                                         #
#                                                                  #
# Description: This is a program that performs association         #
# analyses between microbial feature count data and genome-wide or #
# targeted genotypes with and without adjustment for covariates.   #
####################################################################

echo " "
echo "####################################################################"
echo "# Microbiome and Genetics Analysis Tool (MaGAT) v 1.0.0            #"
echo "# by Zachary D Wallen                                              #"
echo "# Last updated: 9 Feb 2025                                         #"
echo "####################################################################"
echo " "

# Argument parsing
while [[ $# -gt 0 ]]; do
  case "$1" in
  --version)
    echo "MaGAT v 1.0.0"
    exit 0
    ;;
  --help)
    echo " Description: This is a program that performs association             "
    echo " analyses between microbial feature count data and genome-wide or     "
    echo " targeted genotypes with and without adjustment for covariates.       "
    echo "                                                                      "
    echo " Program requirements:                                                "
    echo "    - R base                                                          "
    echo "    - Phyloseq R package                                              "
    echo "    - PLINK v 2                                                       "
    echo "                                                                      "
    echo " Input data consists of the following:                                "
    echo "                                                                      "
    echo " (1) A phyloseq object (saved as RDS file) that contains a table      "
    echo " with feature counts per sample in the otu_table() slot,              "
    echo " taxonomic assignments (if applicable) in the tax_table() slot,       "
    echo " and any sample metadata that one wants to include in the analysis    "
    echo " in the sample_data() slot.                                           "
    echo "                                                                      "
    echo " (2) Genotypes in PLINK binary (.bed, .bim, .fam) format. Sample      "
    echo " names of the phyloseq object and IIDs of the PLINK files should      "
    echo " contain the same IDs (they do not have to be in the same order) for  "
    echo " the program to work properly. Alternatively, genotypes may be        "
    echo " supplied as a directory of VCF files. Sample IDs for VCF files should"
    echo " be in the format 'FID_IID'.                                          "
    echo "                                                                      "
    echo " WARNING: Files created during analysis will take up a lot of         "
    echo " storage space if a large number of variants and features are         "
    echo " analyzed, make sure you have plenty of free storage where the        "
    echo " output directory is located.                                         "
    echo "                                                                      "
    echo " WARNING: When using VCFs, VCFs will be converted to PLINK format     "
    echo " when running this script. The symbol '_' will be taken as an ID      "
    echo " delimiter to split VCF sample IDs into FIDs and IIDs. Please         "
    echo " ensure VCF sample IDs are in the form 'FID_IID' and that the         "
    echo " sample names of the input phyloseq object matches the IID portion    "
    echo " of the VCF ID.                                                       "
    echo "                                                                      "
    echo " Usage: ./MaGAT.sh --phyloseq-object phyloseq_object.rds \            "
    echo "                    --genotype-file-prefix genotype_file_prefix OR    "
    echo "                    --vcf-directory vcf_directory \                   "
    echo "                    --output-dir output_dir \                         "
    echo "                    [additional options]                              "
    echo "                                                                      "
    echo " Parameters:                                                          "
    echo "     --help                 Print the parameter list below then exit. "
    echo "     --version              Print version of the workflow.            "
    echo "     --phyloseq-object      (Required) Phyloseq object in RDS format. "
    echo "                              Make sure tax_table() slot is made up of"
    echo "                              7 columns (Kingdom, Phylum, Class,      "
    echo "                              Order, Family, Genus, Species) if       "
    echo "                              wanting to collapse to certain taxonomic"
    echo "                              levels. Any unclassified levels should  "
    echo "                              be 'NA' for the entry.                  "
    echo "     --genotype-file-prefix (Required) File name prefix for host      "
    echo "                              genotypes in PLINK binary format.       "
    echo "     --vcf-directory        (Required) Directory containing VCF files "
    echo "                              of genotypes. Files can be broken up by "
    echo "                              chromosome and gzipped compressed. All  "
    echo "                              files should contain the same samples.  "
    echo "     --output-dir           (Required) Directory to where the output  "
    echo "                              files will be written to.               "
    echo "     --taxonomic-level      (Optional) The desired taxonomic level at "
    echo "                              which to run the analysis. One of:      "
    echo "                              Kingdom, Phylum, Class, Order, Family,  "
    echo "                              Genus, Species. If this parameter is not"
    echo "                              given the analysis will be performed    "
    echo "                              directly on the counts as they are      "
    echo "                              supplied to the script.                 "
    echo "     --transformation       (Optional) The desired transformation and "
    echo "                              normalization to perform on the feature "
    echo "                              count data. One of: CLR (default),      "
    echo "                              log_TSS, or none.                       "
    echo "     --features             (Optional) Pick specific features that    "
    echo "                              will be in the analysis. Associations   "
    echo "                              will only be performed for the specified"
    echo "                              features. Must be valid feature names   "
    echo "                              that would be returned by the           "
    echo "                              taxa_names() function when ran on the   "
    echo "                              phyloseq object, or if using the        "
    echo "                              --taxonomic-level flag, valid taxonomic "
    echo "                              names of the level in which the phyloseq"
    echo "                              object was collapsed to. Provide feature"
    echo "                              names in a comma separated list with no "
    echo "                              spaces. Parameters --features and       "
    echo "                              --unclassified will be ignored if this  "
    echo "                              is specified.                           "
    echo "     --keep-samples         (Optional) File listing sample names/IIDs "
    echo "                              to be keep in the analysis. Samples not "
    echo "                              included in this list will be excluded  "
    echo "                              from all pre-processing steps and       "
    echo "                              analysis. Each sample name/IID should be"
    echo "                              listed one ID per line. Default is to   "
    echo "                              include all samples.                    "
    echo "     --min-sample-prop      (Optional) The minimum proportion of      "
    echo "                              samples a taxon must be detected in to  "
    echo "                              be included in the analysis. Default is "
    echo "                              10% of samples (i.e., value of 0.1).    "
    echo "     --unclassified         (Optional) Should unclassified taxa be    "
    echo "                              removed? Inclusion of this flag will    "
    echo "                              remove taxa unclassified at the         "
    echo "                              designated level, while leaving this    "
    echo "                              flag out will keep them in. Only        "
    echo "                              applicable when a taxonomic level has   "
    echo "                              been designated with --taxonomic-level, "
    echo "                              will be ignored otherwise.              "
    echo "     --variables            (Optional) Column names of additional     "
    echo "                              variables to be included in the         "
    echo "                              analysis. Must be present in the        "
    echo "                              sample_data() slot of the phyloseq      "
    echo "                              object. Enter as a comma separated list "
    echo "                              with no spaces.                         "
    echo "     --swap                 (Optional) Swap microbiome data for       "
    echo "                              another phenotype that has been listed  "
    echo "                              among the variables specified with the  "
    echo "                              --variables parameter. Listing a        "
    echo "                              variable name here will cause it to     "
    echo "                              become the phenotype and feature count  "
    echo "                              data to be included as predictor        "
    echo "                              variables in the statistical model.     "
    echo "                              Microbial features will be referred to  "
    echo "                              as 'FEATURE' for the variable name in   "
    echo "                              the results.                            "
    echo "     --maf                  (Optional) The minor allele frequency     "
    echo "                              threshold that SNPs must pass in order  "
    echo "                              to be included in the analysis.         "
    echo "                              Default is 0.1.                         "
    echo "     --info                 (Optional) Filter SNPs based on an INFO   "
    echo "                              key contained in VCF file. Input to this"
    echo "                              parameter must be in the following form:"
    echo "                              '\"<key name> <operator> <value>\"'     "
    echo "                              (e.g. '\"R2 >= 0.3\"'), and must be a   "
    echo "                              valid key name and value type, else     "
    echo "                              PLINK will error out. This parameter    "
    echo "                              only applies when specifying a directory"
    echo "                              of VCF files as input with the          "
    echo "                              --vcf-directory parameter. All SNPs that"
    echo "                              do not satisfy the criteria will be     "
    echo "                              removed from analysis.                  "
    echo "     --pca                  (Optional) Calculate genetic principal    "
    echo "                              components (PCs) to be used as          "
    echo "                              covariates in the statistical model.    "
    echo "                              PCs will be generated using SNPs that   "
    echo "                              survive filtering specified with        "
    echo "                              parameters --maf and --info. The first  "
    echo "                              10 PCs will be used as covariates.      "
    echo "     --chromosomes          (Optional) Chromosomes to be included in  "
    echo "                              the analysis. Default is to include all "
    echo "                              chromosomes provided. Input chromosomes "
    echo "                              as comma-separated list and/or ranges   "
    echo "                              with no spaces (e.g. 1-4,22,23,25).     "
    echo "     --base-pair-range      (Optional) Base pair range of SNPs to be  "
    echo "                              included in the analysis. Only valid    "
    echo "                              when a single chromosome is specified   "
    echo "                              with --chromosomes. Should be in the    "
    echo "                              format START-END where START is the     "
    echo "                              first base pair position and END is the "
    echo "                              last base pair position. All SNPs are   "
    echo "                              included by default.                    "
    echo "     --snp-model            (Optional) How should SNP genotypes be    "
    echo "                              coded in the model, additive (ADD),     "
    echo "                              dominant (DOM), recessive (REC). Specify"
    echo "                              one of: ADD (default), DOM, REC.        "
    echo "     --interaction-var      (Optional) Variable to test interaction   "
    echo "                              with SNPs.                              "
    echo "     --joint-test           (Optional) Perform a 2df joint test of SNP"
    echo "                              main effect and its interaction with the"
    echo "                              variable specified in --interaction-var."
    echo "                              Might improve power over a pure test of "
    echo "                              interaction if there is a modest SNP    "
    echo "                              marginal effect.                        "
    echo "     --remove-covariates    (Optional) Remove results for covariates. "
    echo "                              By default, covariates are not removed  "
    echo "                              from PLINK output files, but applying   "
    echo "                              this flag will remove them from outputs."
    echo "     --keep-helper-scripts  (Optional) Keep helper R and shell scripts"
    echo "                              generated during processing of data.    "
    echo "                              Might be useful to keep for records or  "
    echo "                              debugging. If this flag is absent, all  "
    echo "                              helper scripts are removed at end of    "
    echo "                              analysis.                               "
    echo " "
    exit 0
    ;;
  --phyloseq-object)
    PHYLO_OBJ="$2"
    shift 2
    ;;
  --genotype-file-prefix)
    GENOS="$2"
    shift 2
    ;;
  --vcf-directory)
    DOSAGE=$(echo $2 | sed 's#/$##')
    shift 2
    ;;
  --output-dir)
    OUT_DIR=$(echo $2 | sed 's#/$##')
    shift 2
    ;;
  --taxonomic-level)
    TAXA_LVL="$2"
    shift 2
    ;;
  --transformation)
    TRANSF="$2"
    shift 2
    ;;
  --features)
    FEAT="$2"
    shift 2
    ;;
  --keep-samples)
    KEEP_SAMPS="$2"
    shift 2
    ;;
  --min-sample-prop)
    FILTER="$2"
    shift 2
    ;;
  --unclassified)
    UNCLASS=1
    shift 2
    ;;
  --variables)
    VARS="$2"
    shift 2
    ;;
  --swap)
    SWAP="$2"
    shift 2
    ;;
  --maf)
    MAF="$2"
    shift 2
    ;;
  --info)
    INFO="$2"
    shift 2
    ;;
  --pca)
    PCA=1
    shift 2
    ;;
  --chromosomes)
    CHR="$2"
    shift 2
    ;;
  --base-pair-range)
    RANGE="$2"
    shift 2
    ;;
  --snp-model)
    SNP_MOD="$2"
    shift 2
    ;;
  --interaction-variable)
    IXN="$2"
    shift 2
    ;;
  --joint-test)
    JOINT_TEST=1
    shift 2
    ;;
  --remove-covariates)
    REM_COV=1
    shift 2
    ;;
  --keep-helper-scripts)
    KEEP_R=1
    shift 2
    ;;
  *)
    echo "Invalid option: $1" 1>&2
    exit 1
    ;;
  esac
done

# Check that valid arguments were entered, and set defaults where needed

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
  if [[ $(awk '{print NF}' $KEEP_SAMPS | uniq | wc -l) -gt 1 ]] ||
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
if [[ -z "$FILTER" ]]; then
  FILTER=0.1
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
if [[ -z "$MAF" ]]; then
  MAF=0.1
fi

# -q
if [[ ! -z "$INFO" ]]; then
  if [[ -z "$DOSAGE" ]]; then
    echo "ERROR: -q parameter only valid when a directory of VCF files is specified as input using the -d parameter"
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
  START_BP=$(echo $RANGE | awk -F"-" '{print $1}')
  END_BP=$(echo $RANGE | awk -F"-" '{print $2}')
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
if [[ -z "$SNP_MOD" ]]; then
  SNP_MOD=ADD
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
echo " "
echo "*** Parameters for phyloseq data pre-processing ***"
echo " "

# Initiate R script for data pre-processing
echo "### Pre-processing of phyloseq data ###" >${OUT_DIR}/Pre_process_phyloseq_data.R
echo "suppressMessages(library(phyloseq))" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "# Read in phyloseq object" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "ps <- readRDS('${PHYLO_OBJ}')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "# Make sure phyloseq object is in sample x feature orientation" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "if (taxa_are_rows(ps)){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "    ps <- phyloseq(t(otu_table(ps)), sample_data(ps), tax_table(ps))" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "cat('\n','Summary of input phyloseq object:','\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "ps" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R

# Add syntax for collapsing taxa to desired taxonomic level
if [[ -z "$TAXA_LVL" ]]; then
  echo "- Taxa level for analysis: input level"
  echo " "
else
  echo "- Taxa level for analysis: $TAXA_LVL"
  echo " "
  echo "# Collapse phyloseq object to ${TAXA_LVL} level" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "ps <- tax_glom(ps, taxrank = '${TAXA_LVL}', NArm=F)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "cat('\n','Summary after collapsing phyloseq object to ${TAXA_LVL} level:','\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "ps" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
fi

# Add syntax to replace taxa_names() with more meaningful labels if phyloseq object was collapsed
if [[ ! -z "$TAXA_LVL" ]]; then
  if [[ "$TAXA_LVL" == "Species" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "for (i in 1:length(taxa_names(ps))){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:7])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  if (length(taxonomy.sub) == 7){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('s__',taxonomy.sub[7],sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 6){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('g__',taxonomy.sub[6],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 5){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('f__',taxonomy.sub[5],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 4){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('o__',taxonomy.sub[4],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 3){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else{" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  elif [[ "$TAXA_LVL" == "Genus" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "for (i in 1:length(taxa_names(ps))){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:6])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  if (length(taxonomy.sub) == 6){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('g__',taxonomy.sub[6],sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 5){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('f__',taxonomy.sub[5],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 4){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('o__',taxonomy.sub[4],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 3){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else{" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  elif [[ "$TAXA_LVL" == "Family" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "for (i in 1:length(taxa_names(ps))){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:5])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  if (length(taxonomy.sub) == 5){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('f__',taxonomy.sub[5],sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 4){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('o__',taxonomy.sub[4],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 3){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else{" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  elif [[ "$TAXA_LVL" == "Order" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "for (i in 1:length(taxa_names(ps))){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:4])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  if (length(taxonomy.sub) == 4){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('o__',taxonomy.sub[4],sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 3){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else{" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  elif [[ "$TAXA_LVL" == "Class" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "for (i in 1:length(taxa_names(ps))){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:3])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  if (length(taxonomy.sub) == 3){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('c__',taxonomy.sub[3],sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 2){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else{" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  elif [[ "$TAXA_LVL" == "Phylum" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "for (i in 1:length(taxa_names(ps))){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1:2])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  if (length(taxonomy.sub) == 2){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('p__',taxonomy.sub[2],sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else if (length(taxonomy.sub) == 1){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],'_unclass',sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else{" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  elif [[ "$TAXA_LVL" == "Kingdom" ]]; then
    echo "# Replace current labels with shorter taxonomy names" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "for (i in 1:length(taxa_names(ps))){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  taxonomy <- as.vector(tax_table(ps)[i, 1])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo '  taxonomy.sub <- taxonomy[!is.na(taxonomy)]' >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  if (length(taxonomy.sub) == 1){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- paste('k__',taxonomy.sub[1],sep='')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }else{" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "    taxa_names(ps)[i] <- 'unclassified'" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "  }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  fi
fi

# Add syntax to extract certain samples if requested
if [[ ! -z "$KEEP_SAMPS" ]]; then
  echo "# Extract specific samples to use in pre-processing and analysis" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "samp.list <- read.table('$KEEP_SAMPS')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "ps <- prune_samples(samp.list[,1], ps)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "ps <- filter_taxa(prune_samples(samp.list[,1], ps), function(x){sum(x>0)>0}, TRUE)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "cat('\n','Summary after removing samples not found in ${KEEP_SAMPS} file:','\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "ps" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
fi

# Add syntax for feature data transformation
echo "- Transformation for feature count data: ${TRANSF}"
echo " "
echo "# Perform normalization/transformation on feature count data" >>${OUT_DIR}/Pre_process_phyloseq_data.R
if [[ "$TRANSF" = "CLR" ]]; then
  echo "ps.t <- transform_sample_counts(ps, function(x){log(x+1)-mean(log(x+1))})" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
elif [[ "$TRANSF" = "log_TSS" ]]; then
  echo "log.trans <- function(x) {" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "  y <- replace(x, x == 0, min(x[x>0]) / 2)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "  return(log(y))" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "ps.t <- transform_sample_counts(ps, function(x){x/sum(x)})" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "otu_table(ps.t) <- otu_table(apply(otu_table(ps.t), 2, log.trans), taxa_are_rows=F)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
elif [[ "$TRANSF" = "none" ]]; then
  echo "ps.t <- ps" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
else
  echo "ps.t <- transform_sample_counts(ps, function(x){log(x+1)-mean(log(x+1))})" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
fi

# Add syntax for subsetting phyloseq object for specific feature specified by user
if [[ ! -z "$FEAT" ]]; then
  echo "- Subsetting microbiome data for specific feature: $FEAT"
  echo " "
  echo "# Subset microbiome data for specified feature: $FEAT" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "target.feat <- strsplit('${FEAT}', ',')[[1]]" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  if [[ ! -z "$TAXA_LVL" ]]; then
    echo "ps.t <- subset_taxa(ps.t, ${TAXA_LVL} %in% target.feat)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  else
    echo "ps.t <- prune_taxa(target.feat, ps.t)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  fi
  echo "cat('\n','Summary after subsetting data for requested feature(s):', '\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "ps.t" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
fi

# Add syntax for filtering microbiome count data
if [[ ! -z "$FEAT" ]]; then
  :
else
  echo "- Proportion of samples feature must be detected in to be included in this analysis: $FILTER"
  echo " "
  echo "# Filter out features that were detected below minimum proportion of samples equal to ${FILTER}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "filt_feat <- taxa_names(filter_taxa(ps, function(x){sum(x > 0) >= (${FILTER}*length(x))}, TRUE))" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "ps.t <- subset_taxa(ps.t, taxa_names(ps.t) %in% filt_feat)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "cat('\n','Summary after filtering out features that were detected below minimum proportion of samples equal to ${FILTER}:', '\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "ps.t" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
fi

# Add syntax for removing unclassified taxa
if [[ ! -z "$TAXA_LVL" ]]; then
  if [[ ! -z "$FEAT" ]]; then
    :
  elif [[ ! -z "$UNCLASS" ]]; then
    echo "- Removing unclassified features at the $TAXA_LVL level"
    echo " "
    echo "# Remove unclassifed taxa at the ${TAXA_LVL}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "ps.t <- subset_taxa(ps.t, "'!'"is.na(${TAXA_LVL}))" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "cat('\n','Summary after removing unclassified taxa at the ${TAXA_LVL} level:', '\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo "ps.t" >>${OUT_DIR}/Pre_process_phyloseq_data.R
    echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  else
    echo "- Keeping unclassified features at the $TAXA_LVL level in the analysis"
    echo " "
  fi
fi

# Add syntax to create microbiome phenotype file for input to PLINK
if [[ ! -z "$GENOS" ]]; then
  echo "# Read in PLINK fam file" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "fam.file <- read.table('${GENOS}.fam', comment.char='', stringsAsFactors=F)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
elif [[ ! -z "$DOSAGE" ]]; then
  echo "# Read in PLINK fam file" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "fam.file <- read.table('${FAM_FILE}', comment.char='', stringsAsFactors=F)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
fi
echo "cat('\n','Number of samples found in genotype files:', nrow(fam.file), '\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "# Extract microbiome data from phyloseq object" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "feat.df <- data.frame(otu_table(ps.t))" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "cat('\n','Number of samples found in microbiome data:', nrow(feat.df), '\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "# Find overlapping samples between microbome and genotype data" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "fam.file.filt <- fam.file[fam.file[,2] %in% rownames(feat.df),]" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "feat.df.filt <- feat.df[rownames(feat.df) %in% fam.file[,2],,drop=F]" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "if (nrow(feat.df.filt) == 0){ stop('No overlapping samples found between microbiome and genotype data. Are sample names of phyloseq object and IIDs of genotype data concordant?')}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "cat('\n','Number of samples that overlap between phyloseq and genotype data and will be included in phenotype/covariate files:', nrow(feat.df.filt), '\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "# Order microbiome and genotype samples the same so the correct FID and IIDs are added to the microbiome phenotype file" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "feat.df.filt <- feat.df.filt[order(rownames(feat.df.filt)),,drop=F]" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "fam.file.filt <- fam.file.filt[order(fam.file.filt[,2]),]" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo 'if (!identical(rownames(feat.df.filt), fam.file.filt[,2])){ stop("Sample IDs between microbiome and genotype data do not match, even after finding overlaps and ordering the same.")}' >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "# Create microbiome phenotype file" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "pheno.file <- data.frame(FID=fam.file.filt[,1], IID=fam.file.filt[,2], feat.df.filt)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo "write.table(pheno.file, '${OUT_DIR}/phenotype_file.txt', row.names=F, quote=F, sep='\t')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R

# Add syntax to create covariate file with requested variables
if [[ -z "$VARS" ]]; then
  echo "- No additional variables requested to be included in analysis. Linear models will only include the SNP variable as predictor"
  echo " "
else
  echo "- Additional variables to be included in analysis: $(echo $VARS | sed 's/,/ /g')"
  echo " "
  echo "# Check to make sure additional variables provided are in sample data" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "vars <- strsplit('${VARS}', ',')[[1]]" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo 'if (!sum(vars %in% colnames(sample_data(ps.t))) == length(vars)){ stop("Additional variables provided to be included in model were not all found in the phyloseq object sample data.")}' >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "# Extract desired variables from phyloseq object sample data" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "samp.df <- data.frame(sample_data(ps.t)[,vars])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "cat('\n','Additional variables requested were found in phyloseq object sample data and will be included in covariate file:', colnames(samp.df), '\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "# Find overlapping samples between microbome and genotype data" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "samp.df.filt <- samp.df[rownames(samp.df) %in% fam.file.filt[,2],,drop=F]" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "# Order sample data the same as genotype data so the correct FID and IIDs are added to the covariate file" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "samp.df.filt <- samp.df.filt[order(rownames(samp.df.filt)),,drop=F]" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo 'if (!identical(rownames(samp.df.filt), fam.file.filt[,2])){ stop("Sample IDs between sample data and genotype data do not match, even after finding overlaps and ordering the same.")}' >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "# Create covariate file with additional variables to be included in models" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "covar.file <- data.frame(FID=fam.file.filt[,1], IID=fam.file.filt[,2], samp.df.filt)" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "# Dummy-code any categorical variables so PLINK works correctly" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "for (i in 3:ncol(covar.file)){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo '  if (!is.numeric(covar.file[,i]) | length(table(covar.file[,i]))<=2){' >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "    cat('\n','WARNING: While creating covariate file for PLINK, variable [',colnames(covar.file)[i],'] was detected as being categorical. It will be dummy-coded to 2 and 1 for analysis. Please check covariate_file.txt to ensure this was done correctly.','\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "    levels <- names(table(covar.file[,i]))" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "    if (length(levels) == 2){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "      covar.file[,i] <- gsub(levels[1], '2', covar.file[,i])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "      covar.file[,i] <- gsub(levels[2], '1', covar.file[,i])" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "    }else if (length(levels) > 2){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "      stop('Variable [',colnames(covar.file)[i],'] was detected as a categorical variable with >2 levels. Automatic dummy-coding of categorical variables with >2 levels is unsupported at this time. Please recode the variables to numeric values manually and try again.')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "    }else if (length(levels) < 2){" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "      stop('Please check variable [',colnames(covar.file)[i],']. It was detected as categorical with <2 levels. PLINK will not perform correctly if given a variable with only 1 level.')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "    }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "  }" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "}" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo "write.table(covar.file, '${OUT_DIR}/covariate_file.txt', row.names=F, quote=F, sep='\t')" >>${OUT_DIR}/Pre_process_phyloseq_data.R
  echo " " >>${OUT_DIR}/Pre_process_phyloseq_data.R
fi

# Signal end of pre-processing
echo "cat('\n', 'Pre-processing of phyloseq data, and creation of phenotype/covariate files complete.', '\n')" >>${OUT_DIR}/Pre_process_phyloseq_data.R

# Run Pre_process_phyloseq_data.R script to perform the data pre-processing and generation of phenotype and covariate file for PLINK
echo " "
echo " "
echo "*** Performing phyloseq data pre-processing ***"
Rscript ${OUT_DIR}/Pre_process_phyloseq_data.R

# If parameter given to make a covariate as phenotype, create new phenotype/covariate files with covariate as phenotype and add microbiome features as covariates
if [[ ! -z "$SWAP" ]]; then
  echo "### Swapping out covariate as phenotype ###" >${OUT_DIR}/Swap_phenotype.R
  echo "pheno.file <- read.table('${OUT_DIR}/phenotype_file.txt', header=T, stringsAsFactors=F, comment.char='')" >>${OUT_DIR}/Swap_phenotype.R
  echo "covar.file <- read.table('${OUT_DIR}/covariate_file.txt', header=T, stringsAsFactors=F, comment.char='')" >>${OUT_DIR}/Swap_phenotype.R
  echo 'if (!identical(pheno.file[,c("FID","IID")], covar.file[,c("FID","IID")])){ stop("ERROR: FID and IID in phenotype_file.txt and covariate_file.txt do not match, cannot make phenotype swap.")}' >>${OUT_DIR}/Swap_phenotype.R
  echo "pheno.file.new <- cbind(pheno.file[,c('FID','IID')], covar.file[,'${SWAP}', drop=F])" >>${OUT_DIR}/Swap_phenotype.R
  echo "write.table(pheno.file.new, '${OUT_DIR}/phenotype_file.txt', row.names=F, quote=F, sep='\t')" >>${OUT_DIR}/Swap_phenotype.R
  echo "covar.file.new <- cbind(pheno.file, covar.file[,which("'!'"colnames(covar.file) %in% c('FID','IID','${SWAP}'))])" >>${OUT_DIR}/Swap_phenotype.R
  echo "write.table(covar.file.new, '${OUT_DIR}/covariate_file.txt', row.names=F, quote=F, sep='\t')" >>${OUT_DIR}/Swap_phenotype.R
  echo "for (i in 3:ncol(pheno.file)){" >>${OUT_DIR}/Swap_phenotype.R
  echo "  covar.file.new <- cbind(pheno.file[,c(1,2,i)], covar.file[,which("'!'"colnames(covar.file) %in% c('FID','IID','${SWAP}'))])" >>${OUT_DIR}/Swap_phenotype.R
  echo "  colnames(covar.file.new)[3] <- 'FEATURE'" >>${OUT_DIR}/Swap_phenotype.R
  echo "  write.table(covar.file.new, paste('${OUT_DIR}/',colnames(pheno.file)[i],'_cov_file_tEmPoRaRy.txt',sep=''), row.names=F, quote=F, sep='\t')" >>${OUT_DIR}/Swap_phenotype.R
  echo "}" >>${OUT_DIR}/Swap_phenotype.R

  echo " "
  echo " "
  echo "*** Swapping phenotypes, variable ${SWAP} will be used as phenotype while microbiome data will be included as predictors for analyses ***"
  Rscript ${OUT_DIR}/Swap_phenotype.R
fi

# Grab phenotypes for use later
PHENOS=$(awk 'NR == 1{$1=$2=""; print $0}' ${OUT_DIR}/phenotype_file.txt)

# Create list of subjects included in analysis for later
awk 'NR > 1{print $1,$2}' ${OUT_DIR}/covariate_file.txt >${OUT_DIR}/tEmPoRaRy.samp_list.txt

############# END PRE-PROCESS OF PHYLOSEQ DATA #############

############# START PCA #############

if [[ ! -z "$PCA" ]]; then
  echo " "
  echo " "
  echo "*** Calculating genetic PCs to use as covariates in analysis ***"
  echo " "

  # Build PLINK commands to run LD pruning
  echo "### LD prune SNPs for PCA ###" >${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "# Prune for independent SNPs (1st pass)" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "plink2 \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  if [[ ! -z "$GENOS" ]]; then
    echo "--bfile $GENOS \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  fi
  if [[ ! -z "$DOSAGE" ]]; then
    echo "--vcf \${1} dosage=DS \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
    echo "--id-delim _ \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  fi
  echo "--keep ${OUT_DIR}/tEmPoRaRy.samp_list.txt \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "--maf $MAF \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  if [[ ! -z "$DOSAGE" ]] && [[ ! -z "$INFO" ]]; then
    echo "--extract-if-info $INFO \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  fi
  echo "--indep-pairwise 50 5 0.2 \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "--out ${OUT_DIR}/ld_prune_1" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo " " >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "plink2 \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  if [[ ! -z "$GENOS" ]]; then
    echo "--bfile $GENOS \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  fi
  if [[ ! -z "$DOSAGE" ]]; then
    echo "--vcf \${1} dosage=DS \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
    echo "--id-delim _ \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  fi
  echo "--keep ${OUT_DIR}/tEmPoRaRy.samp_list.txt \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "--extract ${OUT_DIR}/ld_prune_1.prune.in \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "--make-pgen \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "--out ${OUT_DIR}/ld_prune_1" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo " " >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "# Prune for independent SNPs (2nd pass)" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "plink2 --pfile ${OUT_DIR}/ld_prune_1 \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "--indep-pairwise 138 5 0.2 \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "--out ${OUT_DIR}/ld_prune_2" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo " " >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "plink2 --pfile ${OUT_DIR}/ld_prune_1 \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "--extract ${OUT_DIR}/ld_prune_2.prune.in \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  echo "--make-pgen \\" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  if [[ ! -z "$GENOS" ]]; then
    echo "--out ${OUT_DIR}/ld_prune_2" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  fi
  if [[ ! -z "$DOSAGE" ]]; then
    echo "--out ${OUT_DIR}/ld_prune_2.\${2}" >>${OUT_DIR}/Run_PLINK_LD_prune.sh
  fi
  chmod +x ${OUT_DIR}/Run_PLINK_LD_prune.sh

  # Run PLINK LD pruning and PCA
  echo " "
  echo "Begin running PLINK commands for LD pruning and PCA..."
  echo " "
  if [[ ! -z "$GENOS" ]]; then
    ${OUT_DIR}/Run_PLINK_LD_prune.sh
    plink2 --pfile ${OUT_DIR}/ld_prune_2 \
      --pca \
      --out ${OUT_DIR}/PCs
  fi
  if [[ ! -z "$DOSAGE" ]]; then
    i=0
    for DOSE_FILE in ${DOSAGE}/*vcf*; do
      echo " "
      echo " "
      i=$(($i + 1))
      ${OUT_DIR}/Run_PLINK_LD_prune.sh $DOSE_FILE $i
      if [[ $i -eq 2 ]]; then
        j=$(($i - 1))
        plink2 --pfile ${OUT_DIR}/ld_prune_2.${i} \
          --pmerge ${OUT_DIR}/ld_prune_2.${j} \
          --make-pgen \
          --out ${OUT_DIR}/ld_prune_2
      elif [[ $i -gt 2 ]]; then
        plink2 --pfile ${OUT_DIR}/ld_prune_2.${i} \
          --pmerge ${OUT_DIR}/ld_prune_2 \
          --make-pgen \
          --out ${OUT_DIR}/ld_prune_2
      fi
    done
    plink2 --pfile ${OUT_DIR}/ld_prune_2 \
      --pca \
      --out ${OUT_DIR}/PCs
  fi
  rm ${OUT_DIR}/ld_prune*

  # Add PCs with high % variation explained to covariate file
  echo "### Extract genetic PCs with high % variation explained ###" >${OUT_DIR}/Add_PCs.R
  echo "pcs <- read.table('${OUT_DIR}/PCs.eigenvec', header=T, comment.char='')" >>${OUT_DIR}/Add_PCs.R
  echo "eigenvals <- read.table('${OUT_DIR}/PCs.eigenval', comment.char='')" >>${OUT_DIR}/Add_PCs.R
  echo "eigenvals <- cbind(paste('PC',1:nrow(eigenvals), sep=''), eigenvals)" >>${OUT_DIR}/Add_PCs.R
  echo "covars <- read.table('${OUT_DIR}/covariate_file.txt', header=T, comment.char='')" >>${OUT_DIR}/Add_PCs.R
  echo "covars <- covars[covars[,2] %in% pcs[,2],]" >>${OUT_DIR}/Add_PCs.R
  echo 'if (!identical(covars[order(covars[,2]),2], pcs[order(pcs[,2]),2])){' >>${OUT_DIR}/Add_PCs.R
  echo "  stop('IIDs not matching between phenotype file and PC file even after subsetting and sorting')" >>${OUT_DIR}/Add_PCs.R
  echo "}else{" >>${OUT_DIR}/Add_PCs.R
  echo "  covars <- cbind(covars[order(covars[,2]),], pcs[order(pcs[,2]), c(eigenvals[eigenvals[,2]>mean(eigenvals[,2])+(2*sd(eigenvals[,2])),1]),FALSE])" >>${OUT_DIR}/Add_PCs.R
  echo "}" >>${OUT_DIR}/Add_PCs.R
  echo "write.table(covars, '${OUT_DIR}/covariate_file.txt', row.names=F, quote=F, sep=' ')" >>${OUT_DIR}/Add_PCs.R
  echo "cat('\n','Genetic PCs that will be included in analysis as covariates:',c(eigenvals[eigenvals[,2]>mean(eigenvals[,2])+(2*sd(eigenvals[,2])),1]),'\n')" >>${OUT_DIR}/Add_PCs.R
  Rscript ${OUT_DIR}/Add_PCs.R
fi

############# END PCA #############

############# START ASSOCIATION ANALYSIS #############

# Remove any previous results if they exist
for pheno in $PHENOS; do
  if [[ ! -z "$(ls ${OUT_DIR}/*${pheno}.glm* 2>/dev/null)" ]]; then
    rm ${OUT_DIR}/*${pheno}.glm*
  fi
done

echo " "
echo " "
echo "*** Performing association analyses with each feature ***"
echo " "
if [[ ! -z "$CHR" ]]; then
  echo "- Running analysis for chromosomes $CHR"
  echo " "
else
  echo "- Running analysis for all chromosomes provided"
  echo " "
fi
if [[ ! -z "$RANGE" ]]; then
  echo "- Running analysis for base pair range $RANGE"
  echo " "
else
  echo "- No base pair range provided, will run analysis for all SNP locations"
  echo " "
fi
echo "- MAF threshold SNPs must pass to be included in analysis is $MAF"
echo " "
if [[ ! -z "$INFO" ]]; then
  echo "- SNPs with $INFO will be included in the analysis"
  echo " "
fi

# Create parameter input for interaction if specified
if [[ ! -z "$IXN" ]]; then
  IXN_PARAM="interaction"
  if [[ -z "$SWAP" ]]; then
    N_VAR=$(echo $(($(awk '{print NF}' ${OUT_DIR}/covariate_file.txt | sort -nu | tail -n 1) - 1)))
    VAR_N=$(($(awk -v RS='\t' "/$IXN/{print NR}" ${OUT_DIR}/covariate_file.txt) - 1))
    IXN_N=$((($N_VAR - 1) + $VAR_N))
    PARAM=$(echo $(seq $N_VAR) $IXN_N)
    echo "- Adding interaction term between SNP and $IXN in the linear models"
    echo " "
  elif [[ ! -z "$SWAP" ]]; then
    file=$(ls -l ${OUT_DIR}/*cov_file_tEmPoRaRy* | awk 'NR == 1{print $NF}')
    N_VAR=$(echo $(($(awk '{print NF}' $file | sort -nu | tail -n 1) - 1)))
    VAR_N=$(($(awk -v RS='\t' "/$IXN/{print NR}" $file) - 1))
    IXN_N=$((($N_VAR - 1) + $VAR_N))
    PARAM=$(echo $(seq $N_VAR) $IXN_N)
    echo "- Adding interaction term between SNP and $IXN in the linear models"
    echo " "
  fi
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
    echo "  Feature~SNP(${SNP_MOD})+$(head -n1 ${OUT_DIR}/covariate_file.txt | cut -f 3- | sed $'s/\t/+/g')"
    echo " "
  elif [[ ! -z "$VARS" ]] && [[ ! -z "$IXN" ]]; then
    echo "  Feature~SNP(${SNP_MOD})+$(head -n1 ${OUT_DIR}/covariate_file.txt | cut -f 3- | sed $'s/\t/+/g')+SNP:${IXN}"
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
    echo "  Feature~SNP(${SNP_MOD})+$(head -n1 ${OUT_DIR}/covariate_file.txt | cut -f 3- | sed $'s/\t/+/g')+SNP:${IXN}"
    echo "  Feature~$(head -n1 ${OUT_DIR}/covariate_file.txt | cut -f 3- | sed $'s/\t/+/g')"
    echo " "
  elif [[ ! -z "$SWAP" ]]; then
    file=$(ls -l ${OUT_DIR}/*cov_file_tEmPoRaRy* | awk 'NR == 1{print $NF}')
    echo "  $(head -n1 ${OUT_DIR}/phenotype_file.txt | cut -f 3-)~SNP(${SNP_MOD})+$(head -n1 $file | cut -f 3- | sed $'s/\t/+/g')+SNP:${IXN}"
    echo "  $(head -n1 ${OUT_DIR}/phenotype_file.txt | cut -f 3-)~$(head -n1 $file | cut -f 3- | sed $'s/\t/+/g')"
    echo " "
  fi
fi

#### BEGIN PLINK ASSOCIATION ANALYSIS ####

echo " "
echo "Begin running PLINK commands for association analysis..."
echo " "

# Run association analysis with microbiome data as phenotype
if [[ -z $SWAP ]]; then

  # If variables supplied, grab names of quantitative one for standardization during analysis
  if [[ ! -z "$VARS" ]]; then
    covars=$(sed -n 1p ${OUT_DIR}/covariate_file.txt | cut -f 3-)
    for i in $covars; do
      levels=$(cat ${OUT_DIR}/covariate_file.txt |
        awk -v COVAR="$i" \
          'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} NR>1 {print $(f[COVAR])}' |
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
  echo "### Run PLINK GLM ###" >${OUT_DIR}/Run_PLINK_GLM.sh
  echo "plink2 \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  if [[ ! -z "$GENOS" ]]; then
    echo "--bfile $GENOS \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  fi
  if [[ ! -z "$DOSAGE" ]]; then
    echo "--vcf \${1} dosage=DS \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    echo "--id-delim _ \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  fi
  echo "--keep ${OUT_DIR}/tEmPoRaRy.samp_list.txt \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  echo "--pheno ${OUT_DIR}/phenotype_file.txt \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  if [[ ! -z "$VARS" ]]; then
    echo "--covar ${OUT_DIR}/covariate_file.txt \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    echo "--variance-standardize $quant_covars \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  fi
  if [[ ! -z "$CHR" ]]; then
    echo "--chr $CHR \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  fi
  if [[ ! -z "$RANGE" ]]; then
    echo "--from-bp $START_BP \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    echo "--to-bp $END_BP \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  fi
  echo "--maf $MAF \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  if [[ ! -z "$DOSAGE" ]] && [[ ! -z "$INFO" ]]; then
    echo "--extract-if-info $INFO \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  fi
  echo "--glm $SNP_PARAM $IXN_PARAM \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  echo "--ci 0.95 \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  if [[ ! -z "$IXN" ]]; then
    echo "--parameters $PARAM \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
  fi
  if [[ ! -z "$JOINT_TEST" ]]; then
    echo "--tests $TESTS \\"
  fi
  if [[ ! -z "$GENOS" ]]; then
    echo "--out ${OUT_DIR}/" >>${OUT_DIR}/Run_PLINK_GLM.sh
  fi
  if [[ ! -z "$DOSAGE" ]]; then
    echo "--out ${OUT_DIR}/\${2}" >>${OUT_DIR}/Run_PLINK_GLM.sh
  fi
  chmod +x ${OUT_DIR}/Run_PLINK_GLM.sh

  # Run PLINK command
  if [[ ! -z "$GENOS" ]]; then
    ${OUT_DIR}/Run_PLINK_GLM.sh
    for pheno in $PHENOS; do
      assoc_file=$(ls ${OUT_DIR}/${pheno}.glm* | awk -F'/' '{print $NF}')
      sed '1s/#//' ${OUT_DIR}/${assoc_file} | sed '/^#/d' >${OUT_DIR}/${assoc_file}.tEmPoRaRy
      rm ${OUT_DIR}/${assoc_file}
      mv ${OUT_DIR}/${assoc_file}.tEmPoRaRy ${OUT_DIR}/${assoc_file}
      rm ${OUT_DIR}/${pheno}*log
    done
  elif [[ ! -z "$DOSAGE" ]]; then
    for DOSE_FILE in ${DOSAGE}/*vcf*; do
      echo " "
      echo " "
      OUT_FILE=$(echo $DOSE_FILE | awk -F'.vcf' '{print $1}' | awk -F'/' '{print $NF}')
      ${OUT_DIR}/Run_PLINK_GLM.sh $DOSE_FILE $OUT_FILE
      rm ${OUT_DIR}/${OUT_FILE}.log
    done
    echo " "
    echo "Merging result files of the same feature..."
    for pheno in $PHENOS; do
      assoc=$(ls -l ${OUT_DIR}/*${pheno}.glm* | grep ${pheno}.glm | head -n1 | awk '{print $NF}' | awk -F'.' '{print $NF}')
      cat ${OUT_DIR}/*${pheno}.glm.${assoc} | sed '1s/#//' | sed '/^#/d' >${OUT_DIR}/${pheno}.glm.${assoc}.tEmPoRaRy
      rm ${OUT_DIR}/*${pheno}.glm.${assoc}
      mv ${OUT_DIR}/${pheno}.glm.${assoc}.tEmPoRaRy ${OUT_DIR}/${pheno}.glm.${assoc}
      echo "$pheno done."
    done
  fi
fi

# Run association analysis for swapped phenotype if applicable
if [[ ! -z $SWAP ]]; then
  for cov_file in ${OUT_DIR}/*cov_file_tEmPoRaRy*; do
    # Get feature name that is now a covariate
    feature=$(echo $cov_file | awk -F'_cov_file_tEmPoRaRy' '{print $1}' | awk -F'/' '{print $NF}')

    # Grab names of quantitative one for standardization during analysis
    covars=$(sed -n 1p $cov_file | cut -f 3-)
    for i in $covars; do
      levels=$(cat $cov_file |
        awk -v COVAR="$i" \
          'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} NR>1 {print $(f[COVAR])}' |
        sort | uniq | wc -l)
      if [[ "$levels" -gt 2 ]]; then
        quant_covars=$(echo "$quant_covars $i")
      fi
    done
    echo " - Following variables were detected as quantitative and will be standardized in analysis:"
    echo "  $quant_covars"
    echo " "

    # Build PLINK command to run
    echo "### Run PLINK GLM ###" >${OUT_DIR}/Run_PLINK_GLM.sh
    echo "plink2 \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    if [[ ! -z "$GENOS" ]]; then
      echo "--bfile $GENOS \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    fi
    if [[ ! -z "$DOSAGE" ]]; then
      echo "--vcf \${1} dosage=DS \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
      echo "--id-delim _ \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    fi
    echo "--keep ${OUT_DIR}/tEmPoRaRy.samp_list.txt \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    echo "--pheno ${OUT_DIR}/phenotype_file.txt \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    echo "--covar ${cov_file} \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    echo "--variance-standardize $quant_covars \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    if [[ ! -z "$CHR" ]]; then
      echo "--chr $CHR \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    fi
    if [[ ! -z "$RANGE" ]]; then
      echo "--from-bp $START_BP \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
      echo "--to-bp $END_BP \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    fi
    echo "--maf $MAF \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    if [[ ! -z "$DOSAGE" ]] && [[ ! -z "$INFO" ]]; then
      echo "--extract-if-info $INFO \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    fi
    echo "--glm $SNP_PARAM $IXN_PARAM \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    echo "--ci 0.95 \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    if [[ ! -z "$IXN" ]]; then
      echo "--parameters $PARAM \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    fi
    if [[ ! -z "$JOINT_TEST" ]]; then
      echo "--tests $TESTS \\" >>${OUT_DIR}/Run_PLINK_GLM.sh
    fi
    if [[ ! -z "$GENOS" ]]; then
      echo "--out ${OUT_DIR}/${feature}" >>${OUT_DIR}/Run_PLINK_GLM.sh
    fi
    if [[ ! -z "$DOSAGE" ]]; then
      echo "--out ${OUT_DIR}/\${2}.${feature}" >>${OUT_DIR}/Run_PLINK_GLM.sh
    fi
    chmod +x ${OUT_DIR}/Run_PLINK_GLM.sh

    # Run PLINK command
    if [[ ! -z "$GENOS" ]]; then
      ${OUT_DIR}/Run_PLINK_GLM.sh
      rm ${OUT_DIR}/${feature}.log
      pattern=$(echo ${feature}.${PHENOS} | sed 's/ //')
      assoc_file=$(ls ${OUT_DIR}/${pattern}* | awk -F'/' '{print $NF}')
      sed '1s/#//' ${OUT_DIR}/${assoc_file} | sed '/^#/d' >${OUT_DIR}/${assoc_file}.tEmPoRaRy
      rm ${OUT_DIR}/${assoc_file}
      mv ${OUT_DIR}/${assoc_file}.tEmPoRaRy ${OUT_DIR}/${assoc_file}
    elif [[ ! -z "$DOSAGE" ]]; then
      for DOSE_FILE in ${DOSAGE}/*vcf*; do
        echo " "
        echo " "
        OUT_FILE=$(echo $DOSE_FILE | awk -F'.vcf' '{print $1}' | awk -F'/' '{print $NF}')
        ${OUT_DIR}/Run_PLINK_GLM.sh $DOSE_FILE $OUT_FILE
        rm ${OUT_DIR}/${OUT_FILE}.${feature}.log
      done
      pattern=$(echo ${feature}.${PHENOS} | sed 's/ //')
      assoc=$(ls -l ${OUT_DIR}/*${pattern}.glm* | grep $pattern.glm | head -n1 | awk '{print $NF}' | awk -F'.' '{print $NF}')
      cat ${OUT_DIR}/*${pattern}.glm.${assoc} | sed '1s/#//' | sed '/^#/d' >${OUT_DIR}/${pattern}.glm.${assoc}.tEmPoRaRy
      rm ${OUT_DIR}/*${pattern}.glm.${assoc}
      mv ${OUT_DIR}/${pattern}.glm.${assoc}.tEmPoRaRy ${OUT_DIR}/${pattern}.glm.${assoc}
    fi
    echo " "
    echo "Association analysis for ${feature} done"
  done
fi

echo " "
echo "PLINK association analysis runs complete."
echo " "

#### END PLINK ANALYSES ####

############# END ASSOCIATION ANALYSIS #############

############# START CLEAN UP #############

# Remove results for covariates if specified
if [[ ! -z $REM_COV ]]; then
  if [[ ! -z "$IXN" ]] && [[ -z "$JOINT_TEST" ]]; then
    echo " "
    echo "Subsetting PLINK result files for SNP(${SNP_MOD})x${IXN} and main effects results..."
    echo " "
    for pheno in $PHENOS; do
      FILE=$(ls ${OUT_DIR}/*${pheno}.glm*)
      cat $FILE |
        awk -v IXN="$IXN" -v SNP_MOD="$SNP_MOD" \
          'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || $(f["TEST"])==SNP_MOD || $(f["TEST"])==IXN || $(f["TEST"])==SNP_MOD"x"IXN)' >${FILE}.mod
      mv ${FILE}.mod $FILE
    done
  elif [[ ! -z "$IXN" ]] && [[ ! -z "$JOINT_TEST" ]]; then
    echo " "
    echo "Subsetting PLINK result files for 2df joint test, SNP(${SNP_MOD})x${IXN}, and main effects results..."
    echo " "
    for pheno in $PHENOS; do
      FILE=$(ls ${OUT_DIR}/*${pheno}.glm*)
      cat $FILE |
        awk -v IXN="$IXN" -v SNP_MOD="$SNP_MOD" \
          'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || $(f["TEST"])==SNP_MOD || $(f["TEST"])==IXN || $(f["TEST"])==SNP_MOD"x"IXN || $(f["TEST"])=="USER_2DF")' >${FILE}.mod
      mv ${FILE}.mod $FILE
    done
  else
    echo " "
    echo "Subsetting PLINK result files for SNP(${SNP_MOD}) results..."
    echo " "
    for pheno in $PHENOS; do
      FILE=$(ls ${OUT_DIR}/*${pheno}.glm*)
      cat $FILE |
        awk -v SNP_MOD="$SNP_MOD" \
          'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || $(f["TEST"])==SNP_MOD)' >${FILE}.mod
      mv ${FILE}.mod $FILE
    done
  fi
else
  echo " "
  echo "Removal of covariates from result file(s) not requested, results for all variables will be kept."
  echo " "
fi

# Clean up files generated during workflow
if [[ ! -d "${OUT_DIR}/0.MaGAT_generated_files" ]]; then
  mkdir ${OUT_DIR}/0.MaGAT_generated_files
fi
mv ${OUT_DIR}/phenotype_file.txt ${OUT_DIR}/0.MaGAT_generated_files/
mv ${OUT_DIR}/covariate_file.txt ${OUT_DIR}/0.MaGAT_generated_files/
if [[ ! -z "$PCA" ]]; then
  mv ${OUT_DIR}/PCs.* ${OUT_DIR}/0.MaGAT_generated_files/
fi

# Clean up helper scripts
if [[ -z "$KEEP_R" ]]; then
  rm ${OUT_DIR}/*.R
  rm ${OUT_DIR}/*.sh
elif [[ ! -z "$KEEP_R" ]]; then
  if [[ ! -d "${OUT_DIR}/0.MaGAT_helper_scripts" ]]; then
    mkdir ${OUT_DIR}/0.MaGAT_helper_scripts
  fi
  mv ${OUT_DIR}/*.R ${OUT_DIR}/0.MaGAT_helper_scripts/
  mv ${OUT_DIR}/*.sh ${OUT_DIR}/0.MaGAT_helper_scripts/
fi

# Gzip files
echo " "
echo "Compressing PLINK result files..."
echo " "
for pheno in $PHENOS; do
  gzip -f ${OUT_DIR}/*${pheno}.glm*
done

# Remove any temporary files
rm ${OUT_DIR}/*tEmPoRaRy*

############# END CLEAN UP #############

echo " "
echo " "
echo "*** Running of MaGAT has completed ***"
echo " "
