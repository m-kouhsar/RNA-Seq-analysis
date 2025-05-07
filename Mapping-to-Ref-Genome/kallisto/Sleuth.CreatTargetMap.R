suppressMessages(library(biomaRt))
suppressMessages(library(dplyr))
mart <- biomaRt::useMart(biomart = "ensembl",
                         dataset = "hsapiens_gene_ensembl",
                         host = "https://www.ensembl.org")
ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version","chromosome_name","start_position","end_position",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id)
write.table(ttg , file = "Human.Transcript.Ensemble.GRCh38.txt", row.names = F , col.names = T , sep = "\t",quote = F)
