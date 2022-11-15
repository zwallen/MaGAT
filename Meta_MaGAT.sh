#!/bin/bash
set -e

######################################################################
# Wrapper program to perform meta-analysis between two (or more)     #
# datasets analyzed using MaGAT.                                     #
# by Zachary D Wallen                                                #
# Last updated: 8 June 2021                                          #
#                                                                    #
# Usage: ./Meta_MaGAT.sh -i list_of_directories_with_results \       #
#                         -o output_directory \                      #
#                         -d directory_of_needed_programs \          #
#                         -v variable_to_meta_analyze                #
#                                                                    #
# Parameters:                                                        #
#     -h    Print the parameter list below then exit.                #
#     -i    (Required) List of directories that contains the PLINK   #
#           resuts from running MaGAT. Each directory should contain #
#           results from different datasets. Should be a comma       #
#           separated list of directory names, no spaces, no '/'     #
#           after the directory name. Only files with matching file  #
#           names will be meta-analyzed. SNPs will be matched up     #
#           based on 'ID' column in PLINK result file, so make sure  #
#           SNP IDs are harmonized between datasets.                 #
#     -o    (Required) Directory to place output.                    #
#     -d    (Required) Path to directory that has the METASOFT java  #
#           program and related files (PLINK format converter script #
#           and Han & Eskin Pvalue table file) and the METAL binary. #
#     -v    (Required) The variable to perform meta-analysis for.    #
#           Must be valid entry that is found in the 'TEST' column of#
#           the result files. If USER_2DF variable specified for     #
#           meta-analysis of 2df joint test results, a sample        #
#           weighted Z-score-based p-value meta-analysis will be     #
#           performed as implemented in the program METAL (as joint  #
#           tests do not result in estimates and their standard      #
#           errors, therefore, meta-analysis must be performed on p  #
#           values). When performing meta-analysis on USER_2DF       #
#           results, direction of effect will be taken from the      #
#           interaction estimate.                                    #
######################################################################

# Argument parsing
while getopts ":hi:o:d:v:" opt; do
  case $opt in
    h)
    echo " "
    echo " Usage: ./Meta_MaGAT.sh -i list_of_directories_with_results \       "
    echo "                         -o output_directory \                      "
    echo "                         -d directory_of_needed_programs \          "
    echo "                         -v variable_to_meta_analyze                "
    echo "                                                                    "
    echo " Parameters:                                                        "
    echo "     -h    Print the parameter list below then exit.                "
    echo "     -i    (Required) List of directories that contains the PLINK   "
    echo "           resuts from running MaGAT. Each directory should contain "
    echo "           results from different datasets. Should be a comma       "
    echo "           separated list of directory names, no spaces, no '/'     "
    echo "           after the directory name. Only files with matching file  "
    echo "           names will be meta-analyzed. SNPs will be matched up     "
    echo "           based on 'ID' column in PLINK result file, so make sure  "
    echo "           SNP IDs are harmonized between datasets.                 "
    echo "     -o    (Required) Directory to place output.                    "
    echo "     -d    (Required) Path to directory that has the METASOFT java  "
    echo "           program and related files (PLINK format converter script "
    echo "           and Han & Eskin Pvalue table file) and the METAL binary. "
    echo "     -v    (Required) The variable to perform meta-analysis for.    "
    echo "           Must be valid entry that is found in the 'TEST' column of"
    echo "           the result files. If USER_2DF variable specified for     "
    echo "           meta-analysis of 2df joint test results, a sample        "
    echo "           weighted Z-score-based p-value meta-analysis will be     "
    echo "           performed as implemented in the program METAL (as joint  "
    echo "           tests do not result in estimates and their standard      "
    echo "           errors, therefore, meta-analysis must be performed on p  "
    echo "           values). When performing meta-analysis on USER_2DF       "
    echo "           results, direction of effect will be taken from the      "
    echo "           interaction estimate.                                    "
    echo " "
    exit 0
    ;;
    i) IN_DIR="$OPTARG"
    ;;
    o) OUT_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    d) META_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    v) VAR="$OPTARG"
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
if [[ -z "$IN_DIR" ]]; then
  echo "ERROR: Argument -i is required, please supply a list of directories with PLINK results from running MaGAT"
  exit 1
fi
if echo $IN_DIR | grep -q ","; then
  :
else
  echo "ERROR: Invalid input given to -i, should be a comma separated list of directory names with no spaces"
  exit 1
fi
if echo $IN_DIR | grep -q " "; then
  echo "ERROR: Invalid input given to -i, should be a comma separated list of directory names with no spaces"
  exit 1
else
  :
fi

# -o
if [[ -z "$OUT_DIR" ]]; then
  echo "ERROR: Argument -o is required, please supply an output directory for meta-analysis results"
  exit 1
fi
if [[ ! -d "$OUT_DIR" ]]; then
  echo "ERROR: Argument -o should be a directory, please supply a directory for output"
  exit 1
fi

# -d
if [[ -z "$META_DIR" ]]; then
  echo "ERROR: Argument -d is required, please supply the directory housing Metasoft.jar, plink2metasoft_modified.py, HanEskinPvalueTable.txt files"
  exit 1
fi
if [[ ! -d "$META_DIR" ]]; then
  echo "ERROR: Argument -d should be a directory, please supply a directory housing Metasoft.jar, plink2metasoft_modified.py, HanEskinPvalueTable.txt files"
  exit 1
fi
if ls $META_DIR | grep "Metasoft.jar"; then
  :
else
  echo "ERROR: Expecting file Metasoft.jar to be in directory specified with -d, but not found"
fi
if ls $META_DIR | grep "plink2metasoft_modified.py"; then
  :
else
  echo "ERROR: Expecting file plink2metasoft_modified.py to be in directory specified with -d, but not found"
fi
if ls $META_DIR | grep "HanEskinPvalueTable.txt"; then
  :
else
  echo "ERROR: Expecting file HanEskinPvalueTable.txt to be in directory specified with -d, but not found"
fi
if ls $META_DIR | grep "metal"; then
  :
else
  echo "ERROR: Expecting METAL binary called 'metal' to be in directory specified with -d, but not found"
fi

# -v
if [[ -z "$VAR" ]]; then
  echo "ERROR: Argument -v is required, please supply a variable to perform meta-analysis for"
  exit 1
fi

############# START META-ANALYSIS #############

echo " "
echo "*** Performing meta-analysis ***"
echo " "

N_DIR=$(echo $IN_DIR | awk -F"," '{print NF}')
echo "- Meta-analysis will be performed for $N_DIR datasets"
echo " "
echo "- Meta-analysis will be performed for variable $VAR"
echo " "

# Create list of files in common between directories
echo "dirs <- strsplit('${IN_DIR}', ',')[[1]]" > ${OUT_DIR}/Rscript.R
echo "dir_files <- c()" >> ${OUT_DIR}/Rscript.R
echo "for (i in 1:length(dirs)){" >> ${OUT_DIR}/Rscript.R
echo "  dir_files <- append(dir_files, list.files(dirs[i], pattern='.glm.'))" >> ${OUT_DIR}/Rscript.R
echo "}" >> ${OUT_DIR}/Rscript.R
echo "count_df <- data.frame(table(dir_files))" >> ${OUT_DIR}/Rscript.R
echo "common_files <- count_df[count_df[,2] == length(dirs),1]" >> ${OUT_DIR}/Rscript.R
echo "write.table(common_files, '${OUT_DIR}/files_in_common_between_datasets.txt', quote=F, row.names=F, col.names=F)" >> ${OUT_DIR}/Rscript.R
Rscript --vanilla ${OUT_DIR}/Rscript.R
rm ${OUT_DIR}/Rscript.R

N_FILES=$(wc -l ${OUT_DIR}/files_in_common_between_datasets.txt | awk '{print $1}')
echo "- Meta-analysis will be performed for $N_FILES result files found in common between directories"
echo " "

if [[ "$VAR" == "USER_2DF" ]]; then
  echo "- A sample weighted Z-score-based p-value meta-analysis will be performed using METAL"
  echo " "
else
  echo "- An inverse variance weighted meta-analysis with both fixed- and random-effects models will be performed using METASOFT"
  echo " "
fi
 
# Perform meta-analysis for each in common file between datasets
cat ${OUT_DIR}/files_in_common_between_datasets.txt | while read line
do
  SECONDS=0
  echo " "
  echo "Performing meta-analysis for files:"
  echo $IN_DIR | sed 's/,/\n/g' | awk -v line="$line" '$1=$1"/"line'
  echo " "
  
  # Subset files for variable results wanting to meta-analyze
  echo "Subsetting files for variable ${VAR}..."
  echo " "
  echo $IN_DIR | sed 's/,/\n/g' | awk -v line="$line" '$1=$1"/"line' | while read line_2
  do
    if [[ $file =~ ".gz" ]] && [[ "$VAR" == "USER_2DF" ]]; then
      snp_check=$(diff <(zcat $line_2 | awk -v VAR="$VAR" 'NR!=1 && $7 == VAR{print $3}') <(zcat $line_2 | grep 'ADDx' | awk '{print $3}'))
      if [[ ! -z "$snp_check" ]]; then
        echo "ERROR: In order to run meta-analysis for USER_2DF, results should be ordered by SNP"
        exit 1
      fi
      zcat $line_2 | awk -v VAR="$VAR" '$7 == VAR{print $1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14}' OFS="\t" | paste - <(zcat $line_2 | grep 'ADDx' | awk '{print $9}') | sed "1s/^/CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tSE\tL95\tU95\tT_OR_F_STAT\tP\tBETA\n/" > ${line_2}.tEmPoRaRy
    elif [[ "$VAR" == "USER_2DF" ]]; then
      snp_check=$(diff <(cat $line_2 | awk -v VAR="$VAR" 'NR!=1 && $7 == VAR{print $3}') <(cat $line_2 | grep 'ADDx' | awk '{print $3}'))
      if [[ ! -z "$snp_check" ]]; then
        echo "ERROR: In order to run meta-analysis for USER_2DF, results should be ordered by SNP"
        exit 1
      fi
      cat $line_2 | awk -v VAR="$VAR" '$7 == VAR{print $1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13,$14}' OFS="\t" | paste - <(cat $line_2 | grep 'ADDx' | awk '{print $9}') | sed "1s/^/CHROM\tPOS\tID\tREF\tALT\tA1\tTEST\tOBS_CT\tSE\tL95\tU95\tT_OR_F_STAT\tP\tBETA\n/" > ${line_2}.tEmPoRaRy
    elif [[ $file =~ ".gz" ]]; then
      zcat $line_2 | awk -v VAR="$VAR" 'NR==1 || $7 == VAR' > ${line_2}.tEmPoRaRy
    else
      cat $line_2 | awk -v VAR="$VAR" 'NR==1 || $7 == VAR' > ${line_2}.tEmPoRaRy
    fi
  done
  FILES=$(echo $IN_DIR | sed 's/,/\n/g' | awk -v line="$line" '$1=$1"/"line".tEmPoRaRy"')
  
  # Create an A2 column in temporary results to use in meta-analysis
  for file in $FILES
  do
    awk '{print $1,$2,$3,$4,$5,$6,"A2",$7,$8,$9,$10,$11,$12,$13,$14}' $file | \
    awk '{if($6==$4) $7=$5; else if ($6==$5) $7=$4; else $7="A2"; print $0}' OFS='\t' > temp.file
    if [[ "$(grep "A2" temp.file | wc -l)" -gt 1 ]]; then
      echo "ERROR: A1 allele does not match allele listed in REF or ALT columns for $(grep "A2" temp.file | awk '{print $3}')"
      exit 1
    fi
    mv temp.file $file
  done
    
  # Run meta-analysis workflows
  
  if [[ "$VAR" == "USER_2DF" ]]; then
    # Create script to run with METAL
    echo "SCHEME SAMPLESIZE" > ${OUT_DIR}/metal_script.txt
    echo "SEPARATOR TAB" >> ${OUT_DIR}/metal_script.txt
    echo "MARKER ID" >> ${OUT_DIR}/metal_script.txt
    echo "ALLELE A1 A2" >> ${OUT_DIR}/metal_script.txt
    echo "EFFECT BETA" >> ${OUT_DIR}/metal_script.txt
    echo "PVALUE P" >> ${OUT_DIR}/metal_script.txt
    echo "WEIGHT OBS_CT" >> ${OUT_DIR}/metal_script.txt
    echo $FILES | sed 's/ /\n/g' | while read line_3
    do
      echo "PROCESS $line_3" >> ${OUT_DIR}/metal_script.txt
    done
    if [[ $file =~ ".gz" ]]; then
      file_name=$(echo $line | awk -F'.gz' '{print $1}')
    else
      file_name=$(echo $line)
    fi
    echo "OUTFILE ${OUT_DIR}/${file_name}. .meta" >> ${OUT_DIR}/metal_script.txt
    echo "ANALYZE HETEROGENEITY" >> ${OUT_DIR}/metal_script.txt
    echo "Script created for running METAL meta-analysis, saved as metal_script.txt"
    echo " "
    
    # Run meta-analysis using METAL
    echo "Running meta-analysis with METAL..."
    echo " "
    
    ${META_DIR}/metal ${OUT_DIR}/metal_script.txt > ${OUT_DIR}/${file_name}.meta.log
    
    # Filter out results for only one dataset and clean up
    filt=$(( $N_DIR - 1 ))
    paste ${OUT_DIR}/${file_name}.1.meta <(awk '{print $7}' ${OUT_DIR}/${file_name}.*.meta | awk -F"?" '{print NF-1}') | awk -v filt="$filt" '$NR == 1 || $NF < filt{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11}' OFS="\t" > ${OUT_DIR}/${file_name}.meta
    rm ${OUT_DIR}/${file_name}.1.meta
    mv ${OUT_DIR}/${file_name}.1.meta.info ${OUT_DIR}/${file_name}.meta.info
    rm ${OUT_DIR}/metal_script.txt
    gzip -f ${OUT_DIR}/${file_name}.meta
  else
    # Format result files to be used in METASOFT
    echo "Formatting result files for METASOFT program..."
    echo " "
    
    python ${META_DIR}/plink2metasoft_modified.py ${OUT_DIR}/meta_input $FILES
    
    echo "Running of plink2metasoft script complete"
    echo " "
    
    # Run meta-analysis using METASOFT
    echo "Running meta-analysis with METASOFT..."
    echo " "
  
    if [[ $file =~ ".gz" ]]; then
      file_name=$(echo $line | awk -F'.gz' '{print $1}')
    else
      file_name=$(echo $line)
    fi
    java -jar ${META_DIR}/Metasoft.jar -input ${OUT_DIR}/meta_input.meta \
                                       -output ${OUT_DIR}/${file_name}.meta \
                                       -log ${OUT_DIR}/${file_name}.log \
                                       -pvalue_table ${META_DIR}/HanEskinPvalueTable.txt \
                                       -mvalue \
                                       -mvalue_p_thres 0.05 \
                                       -verbose
  
    echo "Meta-analysis complete. Adding betas of individual datasets to results..."
    echo " "
  
    # Add beta of individual datasets, reformat headers, and filter for SNPs present in >= 2 studies
    snp_check=$(diff <(awk 'NR!=1{print $1}' ${OUT_DIR}/${file_name}.meta) <(awk '{print $1}' ${OUT_DIR}/meta_input.meta))
    if [[ ! -z "$snp_check" ]]; then
      echo "ERROR: Something went wrong, SNP IDs should be the same between meta-analysis input and results"
      exit 1
    fi
    ind_dat_cols=$(echo $(for i in $(seq 1 $N_DIR); do echo 'PVALUE_STUDY'$i; done ; \
                  for i in $(seq 1 $N_DIR); do echo 'MVALUE_STUDY'$i; done ; \
                  for i in $(seq 1 $N_DIR); do echo 'BETA_STUDY'$i; echo 'STD'$i; done) | \
                  sed 's/ /\t/g')
    sed '1d' ${OUT_DIR}/${file_name}.meta | \
    paste - <(cut -f2- ${OUT_DIR}/meta_input.meta | sed 's/ /\t/g') | \
    sed 's/\t\t/\t/g' | \
    sed "1s/^/SNP_ID\tN_STUDY\tPVALUE_FE\tBETA_FE\tSTD_FE\tPVALUE_RE\tBETA_RE\tSTD_RE\tPVALUE_RE2\tSTAT1_RE2\tSTAT2_RE2\tPVALUE_BE\tI_SQUARE\tQ\tPVALUE_Q\tTAU_SQUARE\t$ind_dat_cols\n/" \
    > ${OUT_DIR}/${file_name}.meta.tmp
    awk '$2 > 1' ${OUT_DIR}/${file_name}.meta.tmp > ${OUT_DIR}/${file_name}.meta
  
    # Clean up current round
    rm ${OUT_DIR}/meta_input*
    rm ${OUT_DIR}/${file_name}.meta.tmp
    gzip -f ${OUT_DIR}/${file_name}.meta
  fi
  
  # Clean up temporary files
  echo $IN_DIR | sed 's/,/\n/g' | awk -v line="$line" '$1=$1"/"line' | while read line_2
  do
    rm ${line_2}.tEmPoRaRy
  done
  echo "Done."
  echo " "
  echo "$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds elapsed."
  echo " "
done
rm ${OUT_DIR}/files_in_common_between_datasets.txt
echo " "
echo "*** Meta-analyses complete ***"
echo " "



