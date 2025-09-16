library(stringr)
library(reshape)

read.star.log <- function(log_file){
  
  if(!file.exists(log_file)){
    stop("File not found: ",log_file)
  }
  
  line_vector <- readLines(log_file)
  
  result <- list()
  j=0
  
  for(i in 1:length(line_vector)){
    
    line_ <- trimws(line_vector[i])
    
    if(line_ != "" & str_detect(line_ , pattern = " [|]\t")){
      j = j+1
      result[[j]] <- str_split_1(line_ , pattern = " [|]\t")
    }
  }
  
  return(do.call(rbind , result))
}

args <- commandArgs(T)

STAR.results.dir <- args[1]    # The directory contains STAR results
OutPrefix <- args[2]           # This script Output files prefix

if(is.na(OutPrefix)){
  OutPrefix = ""
}else{
  OutPrefix = paste0(OutPrefix , ".")
}

if(is.na(STAR.results.dir)){
  STAR.results.dir = "."
}

message("STAR results directory: ", STAR.results.dir)
message("Count matrix files prefix: ", OutPrefix)
cat("\n")

message("Detecting expression files...")

files.expr <- list.files(path = STAR.results.dir , pattern = "ReadsPerGene.out.tab",
                     all.files = F,full.names = T, recursive = T,ignore.case = F)
files.log <- list.files(path = STAR.results.dir , pattern = "Log.final.out",
                         all.files = F,full.names = T, recursive = T,ignore.case = F)

results.files = data.frame(ExprFile = basename(files.expr),LogFile = basename(files.log),Path =dirname( files.expr ))

if(nrow(results.files) == 0){
  stop("No expression file were detected. Check the results directory: ",STAR.results.dir)
}else{
  message(nrow(results.files), " expression files were detected:")
  cat("\n")
}

counts.unstarnded = vector(mode = "list" , length = nrow(results.files))
counts.forward = vector(mode = "list" , length = nrow(results.files))
counts.reverse = vector(mode = "list" , length = nrow(results.files))
log_list = vector(mode = "list" , length = nrow(results.files))

for (i in 1:nrow(results.files)) {
  
  message("Working on expression file ",i,"...")
  
  expr <- read.delim(paste0(results.files$Path[i] , "/",results.files$ExprFile[i]), header = F)
  names(expr) <- c("ID","unstarnded","forward","reverse")
  
  log_ <- read.star.log(log_file = paste0(results.files$Path[i] , "/",results.files$LogFile[i]))
  expr.log <- reshape::melt(expr[1:4 , ],id.vars = "ID")
  expr.log$ID <- paste(expr.log[,1],expr.log[,2] , sep = "_")
  log_ <- rbind.data.frame(log_ , cbind(expr.log[,1],expr.log[,3]))
  
  expr.unstarnded <- expr[c(-1:-4),2,drop = FALSE]
  expr.forward <- expr[c(-1:-4),3,drop = FALSE]
  expr.reverse <- expr[c(-1:-4),4,drop = FALSE]
  expr.reverse <- expr[c(-1:-4),4,drop = FALSE]
  log_value <- log_[,2,drop = FALSE ]
  
  names(expr.unstarnded)  = names(expr.forward)  = names(expr.reverse)  = names(log_value)  = 
                                                str_remove(results.files$ExprFile[i] , pattern = "ReadsPerGene.out.tab")
  
  rownames(expr.unstarnded)  = rownames(expr.forward)  = rownames(expr.reverse) = expr$ID[c(-1:-4)]
  
  counts.unstarnded[[i]] <- expr.unstarnded
  counts.forward[[i]] <- expr.forward
  counts.reverse[[i]] <- expr.reverse
  log_list[[i]] <- log_value[]
}

if(!dir.exists(dirname(OutPrefix))){
  dir.create(dirname(OutPrefix))
}

message("Writing merged results...")
write.csv(results.files , file = paste0(OutPrefix , "STAR.ResultsFiles.csv"), row.names = F)

counts.unstarnded.merged = do.call(cbind.data.frame, counts.unstarnded)
counts.forward.merged = do.call(cbind.data.frame, counts.forward)
counts.reverse.merged = do.call(cbind.data.frame, counts.reverse)
logs.merged = do.call(cbind.data.frame, list(ID = log_[,1],log_list))

write.table(counts.unstarnded.merged , file = paste0(OutPrefix,"STAR.count.unstranded.tsv"),quote = F , sep = "\t" , row.names = T , col.names = T)
write.table(counts.forward.merged , file = paste0(OutPrefix,"STAR.count.forward.tsv"),quote = F , sep = "\t" , row.names = T , col.names = T)
write.table(counts.reverse.merged , file = paste0(OutPrefix,"STAR.count.reverse.tsv"),quote = F , sep = "\t" , row.names = T , col.names = T)
write.csv(logs.merged , file = paste0(OutPrefix,"STAR.report.csv"),row.names = F)


message("ALL DONE!")
