library(bacon)
library(qqman)

calculate_lambda <- function(pvals) {
  pvals <- pvals[!is.na(pvals) & pvals > 0 & pvals < 1]
  chisq <- qchisq(1 - pvals, df = 1)
  lambda <- median(chisq) / qchisq(0.5, df = 1)
  return(lambda)
}


sumstat.all_file <- "Results/UKBBN.PC10.sleuth.DEG.tsv"
OutPrefix <- "Results/bacon/UKBBN.lnc.PC10"
Plot_title <- "UKBBN"

dir.create(dirname(OutPrefix), recursive = T , showWarnings = F)

sumstat.all <- read.table(sumstat.all_file , stringsAsFactors = F , header = T)
sumstat.all <- sumstat.all[!is.na(sumstat.all$se_b),]
sumstat.lnc <- sumstat.all[sumstat.all$biotype == "lncRNA" , ]

bc.lnc = bacon(effectsizes = sumstat.lnc$b , standarderrors = sumstat.lnc$se_b)

sumstat.lnc$bacon_test_stat = as.numeric(bacon::tstat(bc.lnc , corrected = T))
sumstat.lnc$bacon_pval = as.numeric(bacon::pval(bc.lnc , corrected = T))

pdf(file = paste0(OutPrefix , ".bacon.pdf"))

qq(sumstat.lnc$pval, main=paste0(Plot_title , " (lncRNAs, uncorrected data)"))
text(x = 0.5,y = (par("usr")[4]-0.2),
     label = bquote(lambda == .(round(calculate_lambda(sumstat.lnc$pval), 2))),
     adj = c(0, 1),
     cex = 1)

qq(sumstat.lnc$bacon_pval, main=paste0(Plot_title , " (lncRNAs, corrected with bacon)"))
text(x = 0.5,y = (par("usr")[4]-0.2),
     label = bquote(lambda == .(round(calculate_lambda(sumstat.lnc$bacon_pval), 2))),
     adj = c(0, 1),
     cex = 1)

graphics.off()

save(bc.lnc ,sumstat.all , sumstat.lnc , file = paste0(OutPrefix , ".bacon.rdat") )



