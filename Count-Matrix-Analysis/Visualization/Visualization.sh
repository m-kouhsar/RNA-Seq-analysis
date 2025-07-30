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

results_file="/lustre/projects/Research_Project-T111004/Morteza/Ascribed_Study/Results/Ascribed.DESeq2.Confusio.Non_Conf.csv"
gene_col="X"
logFC_col="logFC"
pvalue_col="PValue"
pvalue_cut="0.01"
QQ_title="QQ Plot"
Volcano_title="Volcano Plot"
OutPrefix="Results/Ascribed.DESeq2.Confusio.Non_Conf"

ScriptDir=/lustre/projects/Research_Project-T111004/Morteza/github/RNA_Seq_analysis/Count-Matrix-Analysis/Visualization

###########################################################################

Rscript $ScriptDir/Visualization.R "$results_file" "$gene_col" "$logFC_col" "$pvalue_col" "$pvalue_cut" "$QQ_title" "$Volcano_title" "$OutPrefix"

echo "All done!"
