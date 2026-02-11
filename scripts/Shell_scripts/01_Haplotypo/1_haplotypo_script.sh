#!/bin/bash
#SBATCH --account=def-acgerste
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --time=45:00:00
#SBATCH --mem=0
#SBATCH --mail-user=adamubua@myumanitoba.ca
#SBATCH --mail-type=ALL
#SBATCH --job-name=HSC_2_haplotypo
#SBATCH --output=%x-%j.out

module load python bwa bcftools picard bcftools freebayes gatk vcflib
source env/bin/activate
#pip install pyvcf
cat samples.txt | parallel 'mkdir /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output'
cat samples.txt | parallel 'python haplotypo/bin/3.5/haplotypo1.py -o /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output -thr 32 -c 30 -amb 0 -hapA hapA.fasta -hapB hapB.fasta -idA {}_hapA -idB {}_hapB -f1 /home/abdul/scratch/2023_C.albicans/HSC_2018/data_in/raw_data/renamed/{}_R1.fastq.gz -f2 /home/abdul/scratch/2023_C.albicans/HSC_2018/data_in/raw_data/renamed/{}_R2.fastq.gz -caller bcftools' 
cat samples.txt | parallel 'python haplotypo/bin/3.5/VCFcorr_alleles1.py -A /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output/{}_hapA.pass.snp.vcf -B /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output/{}_hapB.pass.snp.vcf -fastaA hapA.fasta -fastaB hapB.fasta -cA /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output/{}_hapA.corrected_amb0.vcf -cB /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output/{}_hapB.corrected_amb0.vcf -amb 0'
cat samples.txt | parallel 'python haplotypo/bin/3.5/haplomaker.py -o /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output/ -hapA hapA.fasta -hapB hapB.fasta -corrA /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output/{}_hapA.corrected_amb0.vcf -corrB /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output/{}_hapB.corrected_amb0.vcf'
#cat samplessamples.txt | parallel 'python /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/scripts/haplotypo/bin/3.5/haplomaker.py -o fasta_files -hapA hapA.fasta -hapB hapB.fasta -corrA /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output/{}_hapA.corrected_amb0.vcf -corrB /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/{}_output/{}_hapB.corrected_amb0.vcf'
#python haplotypo/bin/3.5/haplotypo1.py -o /home/abdul/scratch/2023_C.albicans/Haplotypo_analyses/data_out/SRR6669933_output -thr 25 -c 30 -amb 0 -hapA hapA.fasta -hapB hapB.fasta -idA SRR6669933_hapA -idB SRR6669933_hapB -f1 /home/abdul/scratch/2023_C.albicans/roppers/data_in/renamed/SRR6669933_1.fastq.gz -f2 /home/abdul/scratch/2023_C.albicans/roppers/data_in/renamed/SRR6669933_2.fastq.gz -caller bcftools
