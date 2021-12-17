#!/bin/bash
#SBATCH -J MAFFT_Core # Job name
#SBATCH -o MAFFT_Core.o%j # Name of output file
#SBATCH -e MAFFT_Core.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=0-18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --cpus-per-task=8
#SBATCH --array=2-4

module load mafft/7.471

data_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Accessory_Genes_Proteins/
out_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Accessory_Genes_Proteins_MSA/
mkdir ${out_dir}

ls ${data_dir} > Proteins.txt

s=`head -n ${SLURM_ARRAY_TASK_ID} Prot_Resubmit.txt | tail -n 1`
mafft-linsi --adjustdirection --anysymbol  --thread 8 ${data_dir}${s} > ${out_dir}${s}