#!/bin/bash
#SBATCH -A Research_Project-MRC164847                   # research project to submit under.
#SBATCH --export=ALL                                    # export all environment variables to the batch job.
#SBATCH -D .                                            # set working directory to .
#SBATCH -p mrcq                                         # submit to the parallel test queue
#SBATCH --time=24:00:00                                  # Maximum wall time for the job
#SBATCH --nodes=1                                       # specify number of nodes.
#SBATCH --ntasks-per-node=16                            # specify number of processors.
#SBATCH --mail-type=END                                 # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk              # email address
#SBATCH --job-name="kallisto"
#SBATCH --output=%x.%j.out
#SBATCH --array=0-15

##############################################################################################################
##############################################################################################################
fastq_dir=/lustre/projects/Research_Project-191391/Project_10202/11_fastp_trimmed                                    # Directory contains fastq files
GenomeIndex=/lustre/projects/Research_Project-191391/Morteza/kallisto/Ref/GRCh38_kallistoIndex.idx               # Genome Index file
result_dir=/lustre/projects/Research_Project-191391/Morteza/kallisto/Results/BDR                                                # Output files prefix
Threads=15

##############################################################################################################
##############################################################################################################
echo Job started on:
date -u
echo -e
fastq_files=(${fastq_dir}/*R1*.fastq.gz)

Num_samp=${#fastq_files[@]}

denom_2=$(( SLURM_ARRAY_TASK_COUNT / 2 ))
window_size=$(( ( Num_samp + denom_2 ) / SLURM_ARRAY_TASK_COUNT ))

lower=$(( SLURM_ARRAY_TASK_ID * window_size ))

fastq_files1=(${fastq_files[@]:${lower}:${window_size}})

if [ "$SLURM_ARRAY_TASK_ID" -eq "$SLURM_ARRAY_TASK_MAX" ]
then
    fastq_files1=(${fastq_files[@]:$lower})
fi

echo Output directory: $result_dir
echo Fastq files directory: $fastq_dir
echo Genome index file: $GenomeIndex
echo Number of samples: $Num_samp
echo Start array index: $SLURM_ARRAY_TASK_MIN
echo End array index : $SLURM_ARRAY_TASK_MAX
echo numer of arrays: $SLURM_ARRAY_TASK_COUNT
echo current array index: $SLURM_ARRAY_TASK_ID
echo Number of samples in current array: ${#fastq_files1[@]}
echo -e
echo "##########################################################################"
echo -e

char="/"
index=$(echo "${fastq_dir}" | awk -F"${char}" '{print NF-1}')
index=$(( index + 2 ))

if [ ${#fastq_files1[@]} != 0 ]
then
    j=0
    for i in ${fastq_files1[@]}
    do 
        name=$i
        R1=$name
        R2=${name/R1/R2}
        
        sample_name=$(echo $name| cut -d'/' -f $index)
        sample_name=$(echo ${sample_name%R1*.fastq.gz})

        result_dir1=${result_dir}/${sample_name}
        mkdir -p $result_dir1

        echo working on sample ${j}/${#fastq_files1[@]}: $sample_name

        kallisto quant -i "$GenomeIndex" -o "$result_dir1" -t $Threads -b 100 "$R1" "$R2"

        j=$(( j + 1 ))
        echo -e
    done 
else
    echo "There is no sample in this array!"
fi

echo Job finished on:
date -u
