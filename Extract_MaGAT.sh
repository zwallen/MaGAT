#!/bin/bash
set -e

######################################################################
# Helper script to extract association results at a paricular pvalue #
# threshold (e.g. P<5E-8 for genome-wide hits) from PLINK results    #
# that were outputted after running MaGAT or meta-analysis results   #
# after running Meta_MaGAT.                                          #
# Last updated: 3 Nov 2022                                           #
#                                                                    #
# Usage: ./Extract_Results.sh -i directory_containing_results \      #
#                             -t type_of_results \                   #
#                             -v variable \                          #
#                             -p pvalue_treshold                     #
#                                                                    #
# Parameters:                                                        #
#     -h    Print the parameter list below then exit.                #
#     -i    (Required) Directory that contains the PLINK resuts from #
#           running MaGAT or meta-analysis results from running      #
#           Meta_MaGAT. The outputted file with extracted results    #
#           will be placed in this directory as well.                #
#     -t    (Required) The type of results being submitted to the    #
#           script. One of: MaGAT or Meta_MaGAT.                     #
#     -v    (Required) The variable to extract results for. Must be  #
#           valid entry that is found in the 'TEST' column of the    #
#           result files. Only used when -t is MaGAT, will be        #
#           ignored otherwise.                                       #
#     -p    (Required) The Pvalue threshold used for extracting      #
#           results. All results with Pvalue less than this threshold#
#           will be extracted. Should be value between 0 and 1. No   #
#           scientific notation, just float values.                  #
######################################################################

# Argument parsing
while getopts ":hi:t:v:p:" opt; do
  case $opt in
    h)
    echo " "
    echo " Usage: ./Extract_Results.sh -i directory_containing_results \      "
    echo "                             -t type_of_results \                   "
    echo "                             -v variable OR -m model \              "
    echo "                             -p pvalue_treshold                     "
    echo "                                                                    "
    echo " Parameters:                                                        "
    echo "     -h    Print the parameter list below then exit.                "
    echo "     -i    (Required) Directory that contains the PLINK resuts from "
    echo "           running MaGAT or meta-analysis results from running      "
    echo "           Meta_MaGAT. The outputted file with extracted results    "
    echo "           will be placed in this directory as well.                "
    echo "     -t    (Required) The type of results being submitted to the    "
    echo "           script. One of: MaGAT or Meta_MaGAT.                     "
    echo "     -v    (Required) The variable to extract results for. Must be  "
    echo "           valid entry that is found in the 'TEST' column of the    "
    echo "           result files. Only used when -t is MaGAT, will be        "
    echo "           ignored otherwise.                                       "
    echo "     -p    (Required) The Pvalue threshold used for extracting      "
    echo "           results. All results with Pvalue less than this threshold"
    echo "           will be extracted. Should be value between 0 and 1. No   "
    echo "           scientific notation, just float values.                  "
    echo " "
    exit 0
    ;;
    i) IN_OUT_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    t) TYPE="$OPTARG"
    ;;
    v) VAR="$OPTARG"
    ;;
    p) P_THRESH="$OPTARG"
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
if [[ -z "$IN_OUT_DIR" ]]; then
  echo "Argument -i is required, please supply a directory with PLINK results from running MaGAT."
  exit 1
fi
if [[ ! -d "$IN_OUT_DIR" ]]; then
  echo "Argument -i should be a directory, please supply a directory with PLINK results from running MaGAT."
  exit 1
fi

# -t
if [[ ! -z "$TYPE" ]]; then
  ARG_LIST="MaGAT Meta_MaGAT"
  if echo $ARG_LIST | grep -q "$TYPE"; then
    :
  else
    echo "Invalid argument given to -t, please specify one of: MaGAT or Meta_MaGAT"
    exit 1
  fi
elif [[ -z "$TYPE" ]]; then
  echo "Argument -t is required, please specify if results being given to script are from running MaGAT or Meta_MaGAT"
  exit 1
fi

# -v
if [[ "$TYPE" == "MaGAT" ]]; then
  if [[ -z "$VAR" ]]; then
    echo "Argument -v is required, please supply a variable to extract results for."
    exit 1
  fi
fi

# -p
if [[ ! -z "$P_THRESH" ]]; then
  if [[ $(echo "$P_THRESH <= 1" | bc) -eq 1 ]]; then
    :
  else
    echo "Invalid value given to -p, should be between 0 and 1"
    exit 1
  fi
  if [[ $(echo "$P_THRESH >= 0" | bc) -eq 1 ]]; then
    :
  else
    echo "Invalid value given to -p, should be between 0 and 1"
    exit 1
  fi
elif [[ -z "$P_THRESH" ]]; then
  echo "Argument -p is required, please specify what pvalue threshold to filter results by"
  exit 1
fi

############# START PROCESS OF EXTRACTING RESULTS #############
echo " "
echo "*** Extracting results ***"
echo " "

if [[ "$TYPE" == "MaGAT" ]]; then
  echo "*** Input result files are specified as MaGAT output ***"
  echo " "
  echo "*** Results reaching P < $P_THRESH will be extracted for variable $VAR ***"
  echo " "
  for file in ${IN_OUT_DIR}/*.glm.*
  do
    if [[ $file =~ ".gz" ]]; then
      echo "Extracting results from file $(echo ${file} | awk -F'/' '{print $NF}')..."
      FEAT_NAME=$(echo $file | awk -F"/" '{print $NF}' | awk -F".glm." '{print $1}')
      zcat $file | \
      awk -v VAR="$VAR" -v P_THRESH="$P_THRESH" \
        'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || ($(f["TEST"])==VAR && $(f["P"])<P_THRESH))' | \
      sed "2,\$s/$/\t${FEAT_NAME}/" | sed '1s/$/\tFEATURE/' > ${file}.tEmPoRaRy
    else
      echo "Extracting results from file $(echo ${file} | awk -F'/' '{print $NF}')..."
      FEAT_NAME=$(echo $file | awk -F"/" '{print $NF}' | awk -F".glm.linear" '{print $1}')
      awk -v VAR="$VAR" -v P_THRESH="$P_THRESH" \
        'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || ($(f["TEST"])==VAR && $(f["P"])<P_THRESH))' | \
      sed "2,\$s/$/\t${FEAT_NAME}/" | sed '1s/$/\tFEATURE/' > ${file}.tEmPoRaRy
    fi
  done
  echo " "
  echo "Merging extracted results for each taxon into one file..."
  echo " "
  pcol=$(awk 'NR==1{for(i=1;i<=NF;i++){if($i=="P"){print i;exit}}}' $(ls ${IN_OUT_DIR}/*.tEmPoRaRy | sed -n 1p))
  cat ${IN_OUT_DIR}/*.tEmPoRaRy | sort -g -k"$pcol" -u > ${IN_OUT_DIR}/Extracted_Results.txt
elif [[ "$TYPE" == "Meta_MaGAT" ]]; then
  echo "*** Input result files are specified as Meta_MaGAT output ***"
  echo " "
  echo "*** Results reaching P < $P_THRESH will be extracted for meta-analysis results ***"
  echo " "
  for file in ${IN_OUT_DIR}/*glm*meta.gz
  do
    echo "Extracting results from file $(echo ${file} | awk -F'/' '{print $NF}')..."
    FEAT_NAME=$(echo $file | awk -F"/" '{print $NF}' | awk -F".glm." '{print $1}')
    if zcat $file | head -n1 | grep -q "Zscore"; then
      zcat $file | \
      awk -v P_THRESH="$P_THRESH" \
        'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || $(f["P-value"])<P_THRESH)' | \
      sed "2,\$s/$/\t${FEAT_NAME}/" | sed '1s/$/\tFEATURE/' > ${file}.tEmPoRaRy
    else
      zcat $file | \
      awk -v P_THRESH="$P_THRESH" \
        'NR==1 {for (i=1; i<=NF; i++) {f[$i] = i}} (NR==1 || ($(f["PVALUE_FE"])<P_THRESH || $(f["PVALUE_RE"])<P_THRESH))' | \
      sed "2,\$s/$/\t${FEAT_NAME}/" | sed '1s/$/\tFEATURE/' > ${file}.tEmPoRaRy
    fi
  done
  echo " "
  echo "Merging extracted results for each taxon into one file..."
  echo " "
  if grep -q "Zscore" ${IN_OUT_DIR}/*.tEmPoRaRy; then
    pcol=$(awk 'NR==1{for(i=1;i<=NF;i++){if($i=="P-value"){print i;exit}}}' $(ls ${IN_OUT_DIR}/*.tEmPoRaRy | sed -n 1p))
    cat ${IN_OUT_DIR}/*.tEmPoRaRy | sort -g -k"$pcol" -u > ${IN_OUT_DIR}/Extracted_Results.txt
  else
    pcol=$(awk 'NR==1{for(i=1;i<=NF;i++){if($i=="PVALUE_FE"){print i;exit}}}' $(ls ${IN_OUT_DIR}/*.tEmPoRaRy | sed -n 1p))
    cat ${IN_OUT_DIR}/*.tEmPoRaRy | sort -g -k"$pcol" -u > ${IN_OUT_DIR}/Extracted_Results.txt
  fi
fi

rm ${IN_OUT_DIR}/*.tEmPoRaRy

echo "*** Done ***"
echo " "
