#!/bin/bash
#SBATCH -J HMM_Scan # Job name
#SBATCH -o HMM_Scan.o%j # Name of output file
#SBATCH -e HMM_Scan.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=18:00:00
#SBATCH --qos=large
#SBATCH --mem=36gb
#SBATCH --array=1-62

module load hmmer

data_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Untagged_Accessory_Genes_Proteins/
out_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Untagged_Accessory_Genes_Proteins_HMM_Alignment/
hmm_path=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Accessory_Gene_HMM_Profile_DB/Accessory_Genes.hmm
mkdir ${out_dir}

ls ${data_dir} > Proteins.txt

s=`head -n ${SLURM_ARRAY_TASK_ID} Prot_Resubmit.txt | tail -n 1`
tbl_out=${out_dir}${s::-4}.tbl
out=${out_dir}${s::-4}.out

hmmscan --tblout ${tbl_out} -o ${out}  ${hmm_path} ${data_dir}${s}