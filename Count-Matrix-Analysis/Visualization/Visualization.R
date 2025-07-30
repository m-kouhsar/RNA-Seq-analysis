
message("Loading required libraries...")

suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(qqman)))
suppressMessages(suppressWarnings(library(EnhancedVolcano)))

#########################################################
args = commandArgs(T)

results_file = args[1]
gene_col = args[2]
logFC_col = args[3]
pvalue_col = args[4]
pvalue_cut = as.numeric(trimws(args[5]))
logFC_cut = as.numeric(trimws(args[6]))
QQ.title = args[7]
Volcano.title = args[8]
OutPrefix = args[9]

cat("##########################################################################\n")
message("Input arguments:")
message("        DEG results file: ",results_file)
message("        Gene IDs column name: ",gene_col)
message("        Log fold change column name: ",logFC_col)
message("        P-values column name: ",pvalue_col)
message("        P-value threshold: ",pvalue_cut)
message("        logFC threshold: ",logFC_cut)
message("        QQ plot title: ",QQ.title)
message("        Volcano plot title:",Volcano.title)
message("        Output files prefix: ",OutPrefix)
cat("##########################################################################\n")

########################################################
message("Reading the data...")
results.table <- read.csv(results_file , stringsAsFactors = F)

chisq <- qchisq(1-results.table[,pvalue_col],1)
inflation <- median(chisq, na.rm = T)/qchisq(0.5,1)

message("Generating plots...")
tiff(filename = paste0(OutPrefix , ".QQ.tif") , res = 300 , units = "in" , height = 8 , width = 8)
qq(results.table[,pvalue_col], main=QQ.title)
text(x = 0.5,y = (par("usr")[4]-0.2),
     label = bquote(lambda == .(round(inflation, 2))),
     adj = c(0, 1),
     cex = 1)
graphics.off()

tiff(filename = paste0(OutPrefix , ".volcano.tif") , res = 300 , units = "in" , height = 8 , width = 8)

EnhancedVolcano(toptable = results.table,lab = results.table[,gene_col] , x = logFC_col , y = pvalue_col , 
                pCutoff = pvalue_cut,title = Volcano.title,FCcutoff = logFC_cut,
                legendPosition = "right",drawConnectors=T)

graphics.off()
