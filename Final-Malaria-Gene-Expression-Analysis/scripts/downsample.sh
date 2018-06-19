#!/bin/bash

cat GSE77499_Dd2CQ_Dd2_3D7_gene_exp.diff | awk '($5=="Dd2CQ" && $6=="Dd2" && $13<=0.05)' > original_only_dd2.diff