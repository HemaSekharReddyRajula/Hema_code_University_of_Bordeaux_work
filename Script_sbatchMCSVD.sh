#!/bin/bash

################################ Slurm options #################################
### Job name
#SBATCH --job-name=MCSVD
### Max run time "hours:minutes:seconds"
#SBATCH --time=200:00:00
### Requirements nodes/servers (default: 1)
#SBATCH --nodes=1
### Requirements cpu/core/task (default: 1)
#SBATCH --ntasks-per-node=1
### Requirements memory (default: 12.5GB per cpu, requesting 25 GB)
#SBATCH --mem-per-cpu=25GB

# Useful information to print
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job Id:' $SLURM_JOB_ID
echo 'Directory:' $(pwd)
# Detail Information:
scontrol show job $SLURM_JOB_ID
echo '########################################'

# load required software modules
module load R/4.1.0

cd /mnt/cbib/INSERM_U897/HEMA/PRS_cSVD/WMH_Single_trait_best_model_generate_score_in_UKBiobank/WMH_Residuals_UKBiobank

# Run the analysis
Rscript scriptMCSVDv2.2_popQ.R
