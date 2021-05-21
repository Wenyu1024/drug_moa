#!/bin/bash -l
#SBATCH --job-name=run_array
#SBATCH --account=project_2003466
#SBATCH --output=./message/output_%j.txt
#SBATCH --error=./message/error_%j.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH -p small
#SBATCH --cpus-per-task=10      
#SBATCH --mem-per-cpu=10000
#SBATCH --array=1-3


# module purge
module load r-env-singularity/4.0.3 

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
    sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
fi

# Specify a temp folder path
# echo "TMPDIR=/scratch/project_2003466/glmnet_modelling_cluster/tmp/" >> ~/.Renviron

srun singularity_wrapper exec Rscript --no-save /projappl/project_2003466/drug_moa/predicting_prism/run_ridge_array_prism.R ${SLURM_ARRAY_TASK_ID}
seff $SLURM_JOBID


