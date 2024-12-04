#!/bin/bash

LINE=$1

echo ${LINE}

MYDELS=$(wc -l ${LINE}.deletions.merged.bed | awk '{print $1}')

MYAMBDELS=$(cat ${LINE}.bad.deletions | awk '{print $1"\t"$2"\t"$3}' | sort | uniq -c | awk '{print $2"\t"$3"\t"$4"\t"$1}' | wc -l)

MYAMBDELSGT1=$(cat ${LINE}.bad.deletions | awk '{print $1"\t"$2"\t"$3}' | sort | uniq -c | awk '$1>1 {print $2"\t"$3"\t"$4"\t"$1}' | wc -l)

MYAMBDELSGT100=$(cat ${LINE}.bad.deletions | awk '{print $1"\t"$2"\t"$3}' | sort | uniq -c | awk '$1>100 {print $2"\t"$3"\t"$4"\t"$1}' | wc -l)

#echo -e "NAMLINE\tDELCOUNT\tAMBDELCOUNT\tAMBDELCOUNTGT1\tAMBDELCOUNTGT100"
echo -e ${LINE}"\t"$MYDELS"\t"$MYAMBDELS"\t"$MYAMBDELSGT1"\t"$MYAMBDELSGT100 >> myambdelresultsperline.tsv
