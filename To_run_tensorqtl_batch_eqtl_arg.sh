#!/bin/sh



 cd projects/ctb-awadalla/GROUP/sc_rnaseq/tensorqtl/

 module load StdEnv/2020 gcc/9.3.0
 module load python
 module load arrow/2.0.0

 source /lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/env_tensorqtl/bin/activate

# python Pyscript_tensorqtl_eqtl_multiplecelltype_all_and_subgroup.py $1
 python Pyscript_tensorqtl_eqtl_multiplecelltype_all_and_subgroup_testPerm.py $1

# salloc --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac

# sbatch --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac
