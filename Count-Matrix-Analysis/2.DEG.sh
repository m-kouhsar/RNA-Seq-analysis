#!/bin/bash
#SBATCH -A Research_Project-MRC164847 # research project to submit under.
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=1:00:00 # Maximum wall time for the job
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=m.kouhsar@exeter.ac.uk # email address
#SBATCH --output=%x.%j.out

#########################################################################

method="DESeq2"
counts_file="./Raw/Ascribed.Count.tsv"
pheno_file="./Raw/Ascribed.Pheno.csv"
trait_var="Group"
numeric_batch="Age,RIN"
factor_batch="Sex"
outliers="16_A008"
gFilter_min_count=5
gFilter_min_prop=0.75
runSVA=no
nSV=4
OutPrefix="./Results/Ascribed"

ScriptDir="/lustre/projects/Research_Project-T111004/Morteza/github/RNA_Seq_analysis/Count-Matrix-Analysis/"

########################################################################

if [[ "$method" == "DESeq2" ]]; then
    Rscript $ScriptDir/DESeq2.R $counts_file $pheno_file $trait_var $numeric_batch $factor_batch $outliers $gFilter_min_count $gFilter_min_prop $runSVA $nSV $OutPrefix
fi
if [[ "$method" == "edgeR" ]]; then
    Rscript $ScriptDir/edgeR.R $counts_file $pheno_file $trait_var $numeric_batch $factor_batch $outliers $gFilter_min_count $gFilter_min_prop $runSVA $nSV $OutPrefix
fi
if [[ "$method" == "limma" ]]; then
    Rscript $ScriptDir/limma.R $counts_file $pheno_file $trait_var $numeric_batch $factor_batch $outliers $gFilter_min_count $gFilter_min_prop $runSVA $nSV $OutPrefix
fi
echo "All done!"
