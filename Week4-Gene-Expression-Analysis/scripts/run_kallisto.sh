#!/bin/bash

# Set variables used in each kallisto run
PUBLIC=/home/linux/ieng6/cs185s/public/week4
GTF=${PUBLIC}/gencode.vM11.primary_assembly.annotation.gtf
KINDEX=${PUBLIC}/mm10_kallisto

# Do a separate kallisto run for each dataset
for prefix in FL_Rep1 FL_Rep2 HL_Rep1 HL_Rep2 MB_Rep1 MB_Rep2
do
    mkdir -p ~/week4/${prefix}
    kallisto quant -t 3 -b 100 \
	-o $prefix --gtf $GTF -i $KINDEX \
	${PUBLIC}/${prefix}_chr5_1.fq.gz ${PUBLIC}/${prefix}_chr5_2.fq.gz
done
