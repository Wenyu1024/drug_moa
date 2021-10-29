#!/bin/bash -l
#SBATCH --job-name=r_parallel
#SBATCH --account=project_2003466
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=small
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=9000

# Load r-env-singularity
module load r-env-singularity

# Clean up .Renviron file in home directory
if test -f ~/.Renviron; then
    sed -i '/TMPDIR/d' ~/.Renviron
fi

# Specify a temp folder path
echo "TMPDIR=/scratch/<project>" >> ~/.Renviron

# Run the R script
srun singularity_wrapper exec Rscript --no-save /projappl/project_2003466/drug_moa/unsupervised_target_exploration/pathway_enrichment.R
seff $SLURM_JOBID