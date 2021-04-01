#!/bin/sh

module load StdEnv/2020 gcc/9.3.0
module load python
module load arrow/2.0.0

python -c "import pyarrow"

cd /lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl

source /lustre03/project/6032391/mjfave/CPTP/ATAC/caqtl/env_tensor_2/bin/activate

python ./Pyscript_tensorqtl_caqtl_exp1.py
