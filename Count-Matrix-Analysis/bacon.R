message("Loading requiered packages...")
library(bacon)
library(qqman)
############################################################
calculate_lambda <- function(pvals) {
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals < 1]
  chisq <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq) / qchisq(0.5, df = 1)
  return(lambda)
}
##########################################################
args <- commandArgs(T)
sumstat_file <- trimws(args[1])
effectsizes_col <- trimws(args[2])
standarderrors_col <- trimws(args[3])
OutPrefix <- trimws(args[4])

message("Input arguments:")
message("        Summary statistics file: ",sumstat_file)
message("        Effect Size column: ",effectsizes_col)
message("        Standar Errors column: ",standarderrors_col)
message("        Output files Prefix: ",OutPrefix)

############################################################
message("Reading input data...")

sumstat <- read.table(sumstat_file , stringsAsFactors = F , header = T)

if(!(effectsizes_col %in% colnames(sumstat))){
  stop("There is no column named ", effectsizes_col, " in summary statistic file.")
}
if(!(standarderrors_col %in% colnames(sumstat))){
  stop("There is no column named ", standarderrors_col, " in summary statistic file.")
}

sumstat <- sumstat[!is.na(sumstat[,standarderrors_col]),]

message("Applying bacon...")
bc_obj = bacon(effectsizes = sumstat[,effectsizes_col] , standarderrors = sumstat[,standarderrors_col])

sumstat$bacon_test_stat = as.numeric(bacon::tstat(bc_obj , corrected = T))
sumstat$bacon_pval = as.numeric(bacon::pval(bc_obj , corrected = T))

message("Saving results...")
pdf(file = paste0(OutPrefix , ".bacon.pdf"))

qq(sumstat$pval, main="QQ Plot - Uncorrected Data")
text(x = 0.5,y = (par("usr")[4]-0.2),
     label = bquote(lambda == .(round(calculate_lambda(sumstat$pval), 2))),
     adj = c(0, 1),
     cex = 1)

qq(sumstat$bacon_pval, main="QQ Plot - Corrected with Bacon")
text(x = 0.5,y = (par("usr")[4]-0.2),
     label = bquote(lambda == .(round(calculate_lambda(sumstat$bacon_pval), 2))),
     adj = c(0, 1),
     cex = 1)

graphics.off()

write.csv(sumstat , file = paste0(OutPrefix , ".bacon.csv") )



