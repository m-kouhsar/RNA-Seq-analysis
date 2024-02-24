#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=200:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
##########################################################################################################

InDir=./RNA-Seq
OutDir=./STAR.Results
GenomeIndexDir=./GenomeIndex

##########################################################################################################

module load STAR
mkdir -p $OutDir
mkdir -p ${OutDir}/temp

for i in ${InDir}/*R1*.fastq.gz
do
    R1=$i
    R2=${i/R1/R2}
    sample_name=$(echo $i| cut -d'/' -f 7)
    sample_name=$(echo ${sample_name%R1*.fastq.gz})
    OutDir=${OutDir}/${sample_name}
    TempDir=${OutDir}/temp/${sample_name}

    #echo "Working on $sample_name ..."

    STAR --runThreadN 16 \
         --quantMode GeneCounts \
         --genomeDir ${GenomeIndexDir} \
         --outSAMtype None \
         --outFileNamePrefix ${OutDir}. \
         --readFilesIn $R1 $R2 \
         --readFilesCommand zcat \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMultimapNmax 20 \
         --genomeLoad LoadAndKeep \
         --limitBAMsortRAM 50000000000 \
         --outFilterScoreMinOverLread 0.5 \
         --outFilterMatchNminOverLread 0.5 \
         --outTmpDir $TempDir \

done
