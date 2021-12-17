#!/bin/bash
#SBATCH -J Calculate_RF_Dist_Mat   # Job name
#SBATCH -o Calculate_RF_Dist_Mat_MP.o # Name of output file
#SBATCH -e Calculate_RF_Dist_Mat_MP.e # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=5-00:00:00
#SBATCH --qos=large
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=64

module load anaconda
source activate /fs/cbcb-software/RedHat-7-x86_64/users/hsmurali/venvs/binnacle_env
data=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Core_Gene_RAXML_Tree/
out=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Core_Gene_RF_Matrix_MP.txt

python Calculate_RF_Distance_Matrix_MP.py ${data} ${out}