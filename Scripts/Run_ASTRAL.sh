#!/bin/bash
#SBATCH -J ASTRAL_Core   # Job name
#SBATCH -o ASTRAL_Core.o # Name of output file
#SBATCH -e ASTRAL_Core.e # Name of error file
#SBATCH --mail-user=hsmurali@terpmail.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=5-00:00:00
#SBATCH --qos=large
#SBATCH --mem=128gb
#SBATCH --cpus-per-task=16

input=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Genes_Trees.tree
output=/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Core_Species_Trees.tree
module load java

java -D"java.library.path=ASTRAL-MP/lib/" -jar ASTRAL-MP/astral.5.15.4.jar -i ${input} -o ${output}