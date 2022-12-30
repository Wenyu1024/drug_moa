#!/bin/bash -l
#SBATCH --job-name=run_multicore
#SBATCH --account=project_2003466
#SBATCH --output=./message/output_run_multisession.txt
#SBATCH --error=./message/error_run_multisession.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH -p small
#SBATCH --cpus-per-task=40  
#SBATCH --mem-per-cpu=1000


# module purge
module load r-env-singularity

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
    sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
fi

# Specify a temp folder path
# echo "TMPDIR=/scratch/project_2003466/glmnet_modelling_cluster/tmp/" >> ~/.Renviron

srun singularity_wrapper exec Rscript --no-save ./run_null_dist.R 


