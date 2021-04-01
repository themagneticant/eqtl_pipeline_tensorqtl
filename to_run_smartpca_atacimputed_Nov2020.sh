#!/bin/sh

cd /lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/scripts

module load nixpkgs/16.09
module load gcc/7.3.0
module load eigensoft/7.2.1

smartpca -p par_file_atacseqimputed_smartPCA.txt > smartpca_atac_ohs_1MbeqtlsKeep.log
