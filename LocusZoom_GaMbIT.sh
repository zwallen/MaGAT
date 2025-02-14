#!/bin/bash
set -e

######################################################################
# Helper script to produce a formatted file for upload to LocusZoom  #
# (https://my.locuszoom.org/) from GaMbIT or Meta_GaMbIT results.    #
# Note: For Meta_GaMbIT results, only will work if SNP IDs are in the#
# format CHR:POS:REF:ALT                                             #
# Last updated: 12 May 2021                                          #
#                                                                    #
# Usage ./LocusZoom_GaMbIT.sh -i result_file \                       #
#                             -t type_of_results \                   #
#                             -v variable OR -m model \              #
#                             -o out_file_name                       #
#                                                                    #
# Parameters:                                                        #
#     -h    Print the parameter list below then exit.                #
#     -i    (Required) Result file from running GaMbIT or Meta_GaMbIT#
#           to format for LocusZoom.                                 #
#     -t    (Required) The type of results being submitted to the    #
#           script. One of: GaMbIT or Meta_GaMbIT.                   #
#     -v    (Required) The variable to extract results for. Must be  #
#           valid entry that is found in the 'TEST' column of the    #
#           result files. Only used when -t is GaMbIT, will be       #
#           ignored otherwise.                                       #
#     -m    (Required) Meta-analysis model to extract results for.   #
#           One of: FE, RE, RE2, FE_RE, FE_RE2, Pvalue. FE_RE/FE_RE2 #
#           options combines FE and RE/RE2 models by taking the FE   #
#           results unless significant heterogeneity was detected    #
#           (Cochran's Q P < 0.1), then will take the RE or RE2      #
#           model. If RE2 chosen, no effect size will be outputted   #
#           to the formatted file. Option Pvalue must be chosen if   #
#           results are from sample weighted z-score based p-value   #
#           meta-analysis from meta-analyzing 2df joint test results.#
#           Only valid when -t is Meta_GaMbIT and will be ignored    #
#           otherwise.                                               #
#     -o    (Required) Name for the output file.                     #
######################################################################

# Argument parsing
while getopts ":hi:t:v:m:o:" opt; do
  case $opt in
    h)
    echo " "
    echo " Usage ./LocusZoom_GaMbIT.sh -i result_file \                       "
    echo "                             -t type_of_results \                   "
    echo "                             -v variable OR -m model \              "
    echo "                             -o out_file_name                       "
    echo "                                                                    "
    echo " Parameters:                                                        "
    echo "     -h    Print the parameter list below then exit.                "
    echo "     -i    (Required) Result file from running GaMbIT or Meta_GaMbIT"
    echo "           to format for LocusZoom.                                 "
    echo "     -t    (Required) The type of results being submitted to the    "
    echo "           script. One of: GaMbIT or Meta_GaMbIT.                   "
    echo "     -v    (Required) The variable to extract results for. Must be  "
    echo "           valid entry that is found in the 'TEST' column of the    "
    echo "           result files. Only used when -t is GaMbIT, will be       "
    echo "           ignored otherwise.                                       "
    echo "     -m    (Required) Meta-analysis model to extract results for.   "
    echo "           One of: FE, RE, RE2, FE_RE, FE_RE2, Pvalue. FE_RE/FE_RE2 "
    echo "           options combines FE and RE/RE2 models by taking the FE   "
    echo "           results unless significant heterogeneity was detected    "
    echo "           (Cochran's Q P < 0.1), then will take the RE or RE2      "
    echo "           model. If RE2 chosen, no effect size will be outputted   "
    echo "           to the formatted file. Option Pvalue must be chosen if   "
    echo "           results are from sample weighted z-score based p value   "
    echo "           meta-analysis from meta-analyzing 2df joint test results."
    echo "           Only valid when -t is Meta_GaMbIT and will be ignored    "
    echo "           otherwise.                                               "
    echo "     -o    (Required) Name for the output file.                     "
    exit 0
    ;;
    i) IN_FILE="$OPTARG"
    ;;
    t) TYPE="$OPTARG"
    ;;
    v) VAR="$OPTARG"
    ;;
    m) MODEL="$OPTARG"
    ;;
    o) OUT_FILE="$OPTARG"
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
if [[ -z "$IN_FILE" ]]; then
  echo "Argument -i is required, please supply a result file from running GaMbIT or Meta_GaMbIT."
  exit 1
fi

# -o
if [[ -z "$OUT_FILE" ]]; then
  echo "Argument -o is required, please supply a name for the outputted file."
  exit 1
fi

# -t
if [[ ! -z "$TYPE" ]]; then
  ARG_LIST="GaMbIT Meta_GaMbIT"
  if echo $ARG_LIST | grep -q "$TYPE"; then
    :
  else
    echo "Invalid argument given to -t, please specify one of: GaMbIT or Meta_GaMbIT"
    exit 1
  fi
elif [[ -z "$TYPE" ]]; then
  echo "Argument -t is required, please specify if results being given to script are from running GaMbIT or Meta_GaMbIT"
  exit 1
fi

# -v
if [[ "$TYPE" == "GaMbIT" ]]; then
  if [[ -z "$VAR" ]]; then
    echo "Argument -v is required, please supply a variable to extract results for."
    exit 1
  fi
fi

# -m
if [[ "$TYPE" == "Meta_GaMbIT" ]]; then
  if [[ ! -z "$MODEL" ]]; then
    ARG_LIST="FE RE RE2 FE_RE FE_RE2 Pvalue"
    if echo $ARG_LIST | grep -q "$MODEL"; then
      :
    else
      echo "Invalid argument given to -m, please specify one of: FE, RE, RE2, FE_RE, FE_RE2, or Pvalue"
      exit 1
    fi
  elif [[ -z "$MODEL" ]]; then
    echo "Argument -m is required, please specify if meta-analysis model to extract results for."
    exit 1
  fi
fi

############# START FORMATTING RESULTS #############
echo " "
echo "*** Formatting $TYPE results for upload to LocusZoom ***"
echo " "

if [[ "$TYPE" == "GaMbIT" ]]; then
  echo "*** Formatted LocusZoom file will contain results for $VAR ***"
  echo " "
  echo "Formatting results for file $IN_FILE..."
  echo " "
  if [[ $IN_FILE =~ ".gz" ]]; then
    zcat $IN_FILE | \
    awk -v VAR="$VAR" \
      'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || $(f["TEST"])==VAR)' OFS='\t' | \
    gzip > ${OUT_FILE}.gz
  else
    awk -v VAR="$VAR" \
      'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || $(f["TEST"])==VAR)' OFS='\t' | \
    gzip > ${OUT_FILE}.gz
  fi
elif [[ "$TYPE" == "Meta_GaMbIT" ]]; then
  echo "*** Formatted LocusZoom file will contain results for $MODEL meta-analysis ***"
  echo " "
  echo " "
  echo "Formatting results for file $IN_FILE..."
  echo " "
  if [[ "$MODEL" == "FE" ]]; then
    echo "results <- read.table('${IN_FILE}', header=TRUE, stringsAsFactors=FALSE, comment.char='')" > Rscript.R
    echo 'format.res <- data.frame(stringr::str_split_fixed(results$SNP_ID, pattern=":", 4))' >> Rscript.R
    echo "colnames(format.res) <- c('CHROM', 'POS', 'REF', 'ALT')" >> Rscript.R
    echo "format.res[,1] <- gsub('chr', '', format.res[,1])" >> Rscript.R
    echo 'format.res <- cbind(format.res, data.frame(ID=results$SNP_ID, P=results$PVALUE_FE, BETA=results$BETA_FE, SE=results$STD_FE))' >> Rscript.R
    echo "write.table(format.res, '${OUT_FILE}', row.names=FALSE, quote=FALSE, sep='\t')" >> Rscript.R
    Rscript Rscript.R
    rm Rscript.R
    gzip -f $OUT_FILE
  elif [[ "$MODEL" == "RE" ]]; then
    echo "results <- read.table('${IN_FILE}', header=TRUE, stringsAsFactors=FALSE, comment.char='')" > Rscript.R
    echo 'format.res <- data.frame(stringr::str_split_fixed(results$SNP_ID, pattern=":", 4))' >> Rscript.R
    echo "colnames(format.res) <- c('CHROM', 'POS', 'REF', 'ALT')" >> Rscript.R
    echo "format.res[,1] <- gsub('chr', '', format.res[,1])" >> Rscript.R
    echo 'format.res <- cbind(format.res, data.frame(ID=results$SNP_ID, P=results$PVALUE_RE, BETA=results$BETA_RE, SE=results$STD_RE))' >> Rscript.R
    echo "write.table(format.res, '${OUT_FILE}', row.names=FALSE, quote=FALSE, sep='\t')" >> Rscript.R
    Rscript Rscript.R
    rm Rscript.R
    gzip -f $OUT_FILE
  elif [[ "$MODEL" == "RE2" ]]; then
    echo "results <- read.table('${IN_FILE}', header=TRUE, stringsAsFactors=FALSE, comment.char='')" > Rscript.R
    echo 'format.res <- data.frame(stringr::str_split_fixed(results$SNP_ID, pattern=":", 4))' >> Rscript.R
    echo "colnames(format.res) <- c('CHROM', 'POS', 'REF', 'ALT')" >> Rscript.R
    echo "format.res[,1] <- gsub('chr', '', format.res[,1])" >> Rscript.R
    echo 'format.res <- cbind(format.res, data.frame(ID=results$SNP_ID, P=results$PVALUE_RE2))' >> Rscript.R
    echo "write.table(format.res, '${OUT_FILE}', row.names=FALSE, quote=FALSE, sep='\t')" >> Rscript.R
    Rscript Rscript.R
    rm Rscript.R
    gzip -f $OUT_FILE
  elif [[ "$MODEL" == "FE_RE" ]]; then
    echo "results <- read.table('${IN_FILE}', header=TRUE, stringsAsFactors=FALSE, comment.char='')" > Rscript.R
    echo 'format.res <- data.frame(stringr::str_split_fixed(results$SNP_ID, pattern=":", 4))' >> Rscript.R
    echo "colnames(format.res) <- c('CHROM', 'POS', 'REF', 'ALT')" >> Rscript.R
    echo "format.res[,1] <- gsub('chr', '', format.res[,1])" >> Rscript.R
    echo 'format.res <- cbind(format.res, data.frame(ID=results$SNP_ID, P=results$PVALUE_FE, BETA=results$BETA_FE, SE=results$STD_FE))' >> Rscript.R
    echo 'format.res$P[results$PVALUE_Q < 0.1] <- results[results$PVALUE_Q < 0.1,"PVALUE_RE"]' >> Rscript.R
    echo 'format.res$BETA[results$PVALUE_Q < 0.1] <- results[results$PVALUE_Q < 0.1,"BETA_RE"]' >> Rscript.R
    echo 'format.res$SE[results$PVALUE_Q < 0.1] <- results[results$PVALUE_Q < 0.1,"STD_RE"]' >> Rscript.R
    echo "write.table(format.res, '${OUT_FILE}', row.names=FALSE, quote=FALSE, sep='\t')" >> Rscript.R
    Rscript Rscript.R
    rm Rscript.R
    gzip -f $OUT_FILE
  elif [[ "$MODEL" == "FE_RE2" ]]; then
    echo "results <- read.table('${IN_FILE}', header=TRUE, stringsAsFactors=FALSE, comment.char='')" > Rscript.R
    echo 'format.res <- data.frame(stringr::str_split_fixed(results$SNP_ID, pattern=":", 4))' >> Rscript.R
    echo "colnames(format.res) <- c('CHROM', 'POS', 'REF', 'ALT')" >> Rscript.R
    echo "format.res[,1] <- gsub('chr', '', format.res[,1])" >> Rscript.R
    echo 'format.res <- cbind(format.res, data.frame(ID=results$SNP_ID, P=results$PVALUE_FE))' >> Rscript.R
    echo 'format.res$P[results$PVALUE_Q < 0.1] <- results[results$PVALUE_Q < 0.1,"PVALUE_RE2"]' >> Rscript.R
    echo "write.table(format.res, '${OUT_FILE}', row.names=FALSE, quote=FALSE, sep='\t')" >> Rscript.R
    Rscript Rscript.R
    rm Rscript.R
    gzip -f $OUT_FILE
  elif [[ "$MODEL" == "Pvalue" ]]; then
    echo "results <- read.table('${IN_FILE}', header=TRUE, stringsAsFactors=FALSE, comment.char='')" > Rscript.R
    echo 'format.res <- data.frame(stringr::str_split_fixed(results$MarkerName, pattern=":", 4))' >> Rscript.R
    echo "colnames(format.res) <- c('CHROM', 'POS', 'REF', 'ALT')" >> Rscript.R
    echo "format.res[,1] <- gsub('chr', '', format.res[,1])" >> Rscript.R
    echo 'format.res <- cbind(format.res, data.frame(ID=results$MarkerName, P=results$P.value, Zscore=results$Zscore))' >> Rscript.R
    echo "write.table(format.res, '${OUT_FILE}.tmp', row.names=FALSE, quote=FALSE, sep='\t')" >> Rscript.R
    Rscript Rscript.R
    rm Rscript.R
    sed '1d' ${OUT_FILE}.tmp | sort -k1,1 -g -k2,2 > $OUT_FILE
    sed -i "1s/^/$(head -n1 ${OUT_FILE}.tmp)\n/" $OUT_FILE
    rm ${OUT_FILE}.tmp
    gzip -f $OUT_FILE
  fi
fi

echo "*** Done ***"
echo " "
