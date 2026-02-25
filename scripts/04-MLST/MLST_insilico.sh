
#Install the stringLST
pip install stringMLST
#Download the C. albicans database
stringMLST.py --getMLST -P candida/calb --species 'Candida albicans'

#Run the analyses for the samples with names in Samples.txt. Ensure the extentions are the same
cat samples.txt | parallel -j 10 'stringMLST.py --predict -P candida/calb -1 {}_R1.fastq.gz -2 {}_R2.fastq.gz -r'
