#!/bin/bash

paste SRR3107187/abundance.tsv SRR3107188/abundance.tsv | cut -f 5,10 | grep -v tpm | awk '(!($1==0 && $2==0))' | datamash ppearson 1:2