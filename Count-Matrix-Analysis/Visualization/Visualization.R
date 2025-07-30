
message("Loading required libraries...")

suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(qqman)))

#########################################################
args = commandArgs(T)

results_file = args[1]
gene_col = args[2]
logFC_col = args[3]
pvalue_col = args[4]
pvalue_cut = args[5]
QQ.title = args[6]
Volcano.title = args[7]
OutPrefix = args[8]

########################################################

results.table <- read.csv(results_file , stringsAsFactors = F)

chisq <- qchisq(1-results.table[,pvalue_col],1)
inflation <- median(chisq, na.rm = T)/qchisq(0.5,1)

tiff(filename = paste0(OutPrefix , ".QQ.tif") , res = 300 , units = "in" , height = 8 , width = 8)
qq(results.table[,pvalue_col], main=QQ.title)
text(x = 0.5,y = (par("usr")[4]-0.2),
     label = bquote(lambda == .(round(inflation, 2))),
     adj = c(0, 1),
     cex = 1)
graphics.off()

tiff(filename = paste0(OutPrefix , ".volcano.tif") , res = 300 , units = "in" , height = 8 , width = 8)

EnhancedVolcano(toptable = results.table,lab = gene_col , x = logFC_col , y = pvalue_col , pCutoff = pvalue_cut)

graphics.off()
