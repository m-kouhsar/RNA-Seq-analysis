#!/bin/bash 
#SBATCH -A Research_Project1 # research project to submit under. 
#SBATCH --export=ALL # export all environment variables to the batch job. 
#SBATCH -D . # set working directory to . 
#SBATCH -p mrcq # submit to the parallel test queue 
#SBATCH --time=24:00:00 # Maximum wall time for the job. 
#SBATCH --nodes=1 # specify number of nodes. 
#SBATCH --ntasks-per-node=16 # specify number of processors. 
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address 

RefDir=./ref
OutDir=${RefDir}/gencode.v38.transcripts_index
RawDir=./RNASeq_fastq

module load Salmon

for i in  ${RawDir}/*R1*
do

sample=$i 
stem=${sample::-22}
out=${OutDir}/${stem:72}

salmon quant -i ${RefDir}/gencode.v38.transcripts_index \
	-l A -1 ${stem}_R1_001_fastp.fastq.gz -2 ${stem}_R2_001_fastp.fastq.gz \
	-p 16 --validateMappings \
	--geneMap ${RefDir}/salmon_gene_transcript_v38.txt \
	-o $out \
done
echo "Salmon processes have been done!"
