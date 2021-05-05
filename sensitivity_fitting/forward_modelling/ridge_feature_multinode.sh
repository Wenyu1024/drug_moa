#!/bin/bash -l
#SBATCH --job-name=r_furrr
#SBATCH --account=project_2003466
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt
#SBATCH --time=72:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH -p large
#SBATCH --mem-per-cpu=5000

# Load r-env-singularity
module load r-env-singularity/4.0.3

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
     sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2003466/glmnet_modelling_cluster/tmp/" >> ~/.Renviron

# Run the R script

srun singularity_wrapper exec RMPISNOW --no-save --slave -f ridge_feature_multinode.R
seff $SLURM_JOBID
