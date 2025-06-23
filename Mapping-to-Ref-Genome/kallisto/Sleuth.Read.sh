#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=Sleuth.Read.%j.out

#################### Input Arguments ################################################################################

# kallisto_res_dir: Directory of kallisto results to read by sleuth
# pheno_file: is a csv file which must contains the following columns:
#             sample: Samples ID
#             path: path of the kallisto resuls (inside kallisto_res_dir) folder for each sample
# target_map_file: (Optional) Target mapping file contains information (eg. gene ID and type) about the transcripts. It must be a tsv file contains "target_id" column
# RemoveOutliers: Do you want to remove outlier samples from the data? (set it to 'yes' or 'no')
# OutPrefix: Results files/images prefix (can contains a directory)
# ScriptDir: Directory of all Scripts related to this analysis 

#######################################################################################################################

kallisto_res_dir="/lustre/projects/Research_Project-191391/Morteza/kallisto/BDR/"
pheno_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/BDR.Phenotype.Labeled.csv"
target_map_file="/lustre/projects/Research_Project-191391/Morteza/kallisto/Human.Transcript.Ensemble.GRCh38.txt"
factor_variables="Phenotype,Gender"  
numeric_variables="Age,RIN"
RemoveOutliers=Yes
OutPrefix="/lustre/projects/Research_Project-191391/Morteza/kallisto/Results/BDR"
ScriptDir="/lustre/projects/Research_Project-191391/Morteza/github/RNA-Seq-analysis/Mapping-to-Ref-Genome/kallisto/"

#######################################################################################################################

Rscript ${ScriptDir}/Sleuth.Read.R "$kallisto_res_dir" "$pheno_file" "$target_map_file" "$RemoveOutliers" "$OutPrefix" "$ScriptDir" 
