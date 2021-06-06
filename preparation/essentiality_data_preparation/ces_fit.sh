#!/bin/bash
#SBATCH --job-name=r_serial
#SBATCH --account=project_2003466
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=4000

# Load r-env-singularity
module load r-env-singularity/4.0.2

# Clean up .Renviron file in home directory
# if test -f ~/.Renviron; then
#     sed -i '/TMPDIR/d' ~/.Renviron
#     sed -i '/OMP_NUM_THREADS/d' ~/.Renviron
# fi

# Specify a temp folder path
# echo "TMPDIR=/scratch/<project>" >> ~/.Renviron

# Run the R script
srun singularity_wrapper exec Rscript --no-save fit_ces_and_impute.R
seff $SLURM_JOBID
