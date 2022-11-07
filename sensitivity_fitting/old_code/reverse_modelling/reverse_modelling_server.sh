#!/bin/bash -l
#SBATCH --job-name=reverse_modelling
#SBATCH --account=project_2003466
#SBATCH --output=output_1node.txt
#SBATCH --error=error_1node.txt
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p small
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=1000

# Load r-env-singularity
module load r-env-singularity/4.0.2

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
     sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/project_2003466/glmnet_modelling_cluster/tmp" >> ~/.Renviron

# Match thread and core numbers
# echo "OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK" >> ~/.Renviron

# Run the R script

srun singularity_wrapper exec Rscript --no-save --slave reverse_modelling_server.R
seff $SLURM_JOBID
