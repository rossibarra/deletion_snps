#!/bin/bash

# B97biallelic.snpsNindels.for.SFS.tsv
MYFILE=$1

tail -n +2 ${MYFILE} | awk '{print $1"\t"$2"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24"\t"$25"\t"$26"\t"$27"\t"$28"\t"$29"\t"$30"\t"$31}' > $(basename ${MYFILE} .tsv).bed
#cat B97biallelic.snpsNindels.for.SFS.headers
