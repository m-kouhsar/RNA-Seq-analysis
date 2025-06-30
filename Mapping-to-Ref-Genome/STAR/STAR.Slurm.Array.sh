#!/bin/bash
#SBATCH -A  Research_Project-MRC164847                  # research project to submit under.
#SBATCH --export=ALL                                    # export all environment variables to the batch job.
#SBATCH -D .                                            # set working directory to .
#SBATCH -p mrcq                                         # submit to the parallel test queue
#SBATCH --time=10:00:00                                  # Maximum wall time for the job
#SBATCH --nodes=1                                       # specify number of nodes.
#SBATCH --ntasks-per-node=16                            # specify number of processors.
#SBATCH --mail-type=END                                 # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk              # email address
#SBATCH --array=0-15
#SBATCH --job-name=STAR

################################################################################################################
DataDir="/lustre/projects/Research_Project-T111004/Project_11556/X0052/11_cuta_trimmed"
GenomeIndexDir="/lustre/projects/Research_Project-T111004/Morteza/GenomeIndex/STAR/GRCh38.p14.PRI.GENECODE.R48"
OutDir="/lustre/projects/Research_Project-T111004/Morteza/Ascribed_Study/Results/STAR"

################################################################################################################

echo Job started on:
date -u

module load STAR

mkdir -p $OutDir
mkdir -p ${OutDir}/temp

samples=($(ls ${DataDir}/*R1*.fastq.gz))
Num_samp=${#samples[@]}
window_size=$(( Num_samp / SLURM_ARRAY_TASK_COUNT + 1 ))

lower=$(( SLURM_ARRAY_TASK_ID * window_size ))
next=$(( SLURM_ARRAY_TASK_ID + 1 ))
upper=$(( next * window_size ))

if [ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MAX" ]; then
    upper=$Num_samp
fi

char="/"
index=$(echo "${DataDir}" | awk -F"${char}" '{print NF-1}')
index=$(( index + 2 ))

for((i=$lower;i<$upper;i++))
do

    name=${samples[$i]}
    R1=$name
    R2=${name/R1/R2}
    
    sample_name=$(echo $name| cut -d'/' -f $index)
    sample_name=$(echo ${sample_name%R1*.fastq.gz})
    temp_dir=${OutDir}/temp/${sample_name}
	OutDir1=${OutDir}/${sample_name}

    echo start sample $sample_name
    
    STAR --runThreadN 16 \
         --quantMode GeneCounts \
         --genomeDir $GenomeIndexDir \
         --outSAMtype None \
         --outFileNamePrefix ${OutDir1} \
         --readFilesIn $R1 $R2 \
         --readFilesCommand zcat \
         --alignSJoverhangMin 8 \
         --alignSJDBoverhangMin 1 \
         --outFilterMultimapNmax 20 \
         --genomeLoad LoadAndKeep \
         --limitBAMsortRAM 50000000000 \
         --outFilterScoreMinOverLread 0.5 \
         --outFilterMatchNminOverLread 0.5 \
         --outTmpDir $temp_dir \

done

echo Job finished on:
date -u
