#!/bin/bash

MYTSV=$1

echo "Print ${MYTSV} total lines"
zcat ${MYTSV} | wc -l | awk '{print $1}'

echo "Print ${MYTSV} deletions"
zcat ${MYTSV} | awk '{print $3}' | awk -F"," '{print $1"\t"$2}' | awk '{sum+=$1} END {print sum}'

echo "Print ${MYTSV} missing data"
zcat ${MYTSV} | awk '{print $3}' | awk -F"," '{print $1"\t"$2}' | awk '{sum+=$2} END {print sum}'

echo "Print ${MYTSV} deletion and missing data"
zcat ${MYTSV} | awk '{print $3}' | awk '/1,1/ {sum++} END {print sum}'
