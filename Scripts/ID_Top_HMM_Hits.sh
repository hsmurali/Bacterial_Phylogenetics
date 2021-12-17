#!/bin/bash
#SBATCH -J Filter_HMM_Hits # Job name
#SBATCH -o Filter_HMM_Hits.o%j # Name of output file
#SBATCH -e Filter_HMM_Hits.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=large
#SBATCH --mem=36gb
#SBATCH --array=2-25

module load anaconda
source activate /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/venvs/binnacle_env

filepath=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Untagged_Accessory_Genes_Proteins_HMM_Alignment/out/
tophit_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Untagged_Accessory_Genes_Proteins_HMM_Alignments_Tophits/
mkdir ${tophit_dir}
untagged_gene_seq=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Untagged_Accessory_Genes_Proteins/

ls ${filepath} | grep "^GCF" > HMM.txt
f=`head -n ${SLURM_ARRAY_TASK_ID} HMM-Resubmit.txt | tail -n 1`

python Parse_HMMER_Outputs.py ${filepath} ${tophit_dir} ${untagged_gene_seq} ${f}

