#!/bin/bash -l
#SBATCH --job-name=r_furrr
#SBATCH --account=project_2003466
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH -p small
#SBATCH --mem-per-cpu=1000

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
srun singularity_wrapper exec Rscript --no-save --slave supervised_target_pred.R
seff $SLURM_JOBID