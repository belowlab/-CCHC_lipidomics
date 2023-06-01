#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=24G
#SBATCH --array=1-830
#SBATCH --output=/data100t1/home/wanying/CCHC/lipidomics/code/fastGWA_lip_species_adj_sex_age_pc_job_slurm_%A_%a.out

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

sed -n "${SLURM_ARRAY_TASK_ID}p" < /data100t1/home/wanying/CCHC/lipidomics/code/fastGWA_run.txt | bash
