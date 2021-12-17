#!/bin/bash
#SBATCH -J HMM_Build # Job name
#SBATCH -o HMM_Build.o%j # Name of output file
#SBATCH -e HMM_Build.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=large
#SBATCH --mem=36gb
#SBATCH --array=1-54

module load hmmer

data_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Accessory_Genes_Proteins_MSA/
out_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Accessory_Genes_Proteins_MSA_HMM/
mkdir ${out_dir}

ls ${data_dir} > Proteins.txt

s=`head -n ${SLURM_ARRAY_TASK_ID} Prot_Resubmit.txt | tail -n 1`
out_hmm=${out_dir}${s::-4}.hmm
hmmbuild ${out_hmm} ${data_dir}${s}