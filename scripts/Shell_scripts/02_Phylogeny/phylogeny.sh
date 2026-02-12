#!/bin/bash
#SBATCH --account=def-acgerste
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --time=268:00:00
#SBATCH --partition=genlm
#SBATCH --mem=1200G # all memory on the node
#SBATCH --mail-user=adamubua@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=All_sample_DP_H_phylogeny
#SBATCH --output=%x-%j.out

./FastTreeDbl -fastest -nt -gtr All_samples_alignment.fa > Phylogeny_for_all_isolates.nwk
