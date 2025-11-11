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

results_file="./Results/NIH.limma.AD.C.csv"
gene_col="Gene"
logFC_col="logFC"
pvalue_col="PValue"
pvalue_cut="0.01"
logFC_cut="1"
QQ_title="QQ Plot"
Volcano_title="Volcano Plot"
OutPrefix="./Results/NIH.limma"

ScriptDir="./RNA-Seq-analysis/Count-Matrix-Analysis"

###########################################################################

Rscript $ScriptDir/Visualization.R "$results_file" "$gene_col" "$logFC_col" "$pvalue_col" "$pvalue_cut" "$logFC_cut" "$QQ_title" "$Volcano_title" "$OutPrefix"

echo "All done!"
