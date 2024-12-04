#!/bin/bash

PATH=/gpool/bin/bedtools2/bin:$PATH

for i in $(ls *.deletions.bed)
  do
     sort -k1,1V -k2,2n $i > $(basename $i .bed).sorted.bed
     bedtools merge -i $(basename $i .bed).sorted.bed > $(basename $i .bed).merged.bed
done
