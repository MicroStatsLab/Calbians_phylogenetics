#!/bin/bash
#SBATCH --account=def-acgerste
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=5:00:00
#SBATCH --mem=0 # all memory on the node
#SBATCH --mail-user=adamubua@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=MSA
#SBATCH --output=%x-%j.out

cat *.fasta >> All_samples_alignment.fa
#python phylip.py
#python create_phylip.py
#nano create_phylip_2.py
