#!/bin/bash
#SBATCH -J MSA_Stats   # Job name
#SBATCH -o MSA_Stats.o%j # Name of output file
#SBATCH -e MSA_Stats.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --array=5

module load anaconda
source activate /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/venvs/binnacle_env
input_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Core_Gene_Proteins_MSA/
output_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/MSA_Stats/
mkdir ${output_dir}
mkdir ${output_dir}Core_Genes/

ls ${input_dir} > Core_Genes.txt
f=`head -n ${SLURM_ARRAY_TASK_ID} MSA_Resubmit.txt | tail -n 1`

python Calculate_MSA_Stats.py ${input_dir}${f} ${output_dir}Core_Genes/${f}