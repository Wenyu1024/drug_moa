#!/bin/bash -l
#SBATCH --job-name=run_array
#SBATCH --account=project_2003466
#SBATCH --output=output_%j.txt
#SBATCH --error=error_%j.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH -p small
#SBATCH --cpus-per-task=5      
#SBATCH --mem-per-cpu=50000
#SBATCH --array=1-365


# module purge
module load r-env-singularity/4.0.3 

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
    sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2003466/glmnet_modelling_cluster/tmp" >> ~/.Renviron

srun singularity_wrapper exec Rscript --no-save /projappl/project_2003466/drug_moa/sensitivity_fitting/run_ridge_array.R ${SLURM_ARRAY_TASK_ID}
seff $SLURM_JOBID

