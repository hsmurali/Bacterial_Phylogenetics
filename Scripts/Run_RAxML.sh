#!/bin/bash
#SBATCH -J RAXML_Core # Job name
#SBATCH -o RAXML_Core.o%j # Name of output file
#SBATCH -e RAXML_Core.e%j # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=18:00:00
#SBATCH --qos=throughput
#SBATCH --mem=36gb
#SBATCH --cpus-per-task=16
#SBATCH --array=56-61

module load raxml

data_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Core_Gene_Nucleotides_MSA/
out_dir=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Core_Gene_RAXML_Tree/
mkdir ${out_dir}
cd ${out_dir}

ls ${data_dir} | grep ".fna$" > MSA.txt

s=`head -n ${SLURM_ARRAY_TASK_ID} MSA.txt | tail -n 1`
out_tree_name=T_${s::-4}
echo ${out_tree_name}
raxml -f a -m GTRGAMMA -p 12345 -x 12345 -T 16 -s ${data_dir}${s} -n ${out_tree_name} -# 20  --print-identical-sequences