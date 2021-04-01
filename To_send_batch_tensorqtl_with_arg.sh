#!/bin/sh



cd /lustre03/project/6032391/GROUP/sc_rnaseq/tensorqtl/scripts


sbatch --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac To_run_tensorqtl_batch_eqtl_arg.sh Bcells
sbatch --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac To_run_tensorqtl_batch_eqtl_arg.sh Mono
sbatch --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac To_run_tensorqtl_batch_eqtl_arg.sh CD4
sbatch --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac To_run_tensorqtl_batch_eqtl_arg.sh CD8
sbatch --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac To_run_tensorqtl_batch_eqtl_arg.sh DC
sbatch --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac To_run_tensorqtl_batch_eqtl_arg.sh NK
sbatch --time=2:0:0 --ntasks=1 --cpus-per-task=40  --mem=96Gb --gres=gpu:1 --account=def-awadalla-ac To_run_tensorqtl_batch_eqtl_arg.sh bulk


