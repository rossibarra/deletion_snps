#!/bin/bash

PATH=/gpool/bin/bedtools2/bin:$PATH

#
MYGENOME=$1
#
MYDELS=$2

bedtools subtract -a ${MYGENOME} -b ${MYDELS} > $(basename ${MYDELS} .deletions.bed).nondeletions.bed
