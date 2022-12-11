library(dplyr)

ht=read.table("htseq_OUT", sep="\t", header=FALSE)
colnames(ht) =c("hgnc_symbol", "count")
feature =read.table("Featurelength", sep="\t", header=TRUE)
feature <- feature %>% mutate(length=end_position - start_position)
new<-inner_join(ht, feature)
new<-new %>%  mutate (RPK= (count/length))
#count all RPK 
rpk_sum = sum(new$RPK)
rpk_scaling = rpk_sum/1000000
#print (rpk_scaling)
new <- new %>% mutate(TPM= RPK/rpk_scaling)
write.csv(new, file="htseq_OUT_withTPM.csv", quote=FALSE, row.names=FALSE)

op1=read.table("NEW_OUTPUT",sep="\t", header=TRUE)
colnames(op1)
colnames(op1)[19]= "hgnc_symbol"

final = inner_join(op1, new, by=c("hgnc_symbol"))
colnames(final)[2]="HLA_wild"
colnames(final)[11]="HLA_mutant"
colnames(final)[28]="htseq_count"
colnames(final)[30:34]=c("Gene_start_position","Gene_end_position", "Gene_length", "Gene_RPKM","Gene_TPM")
write.csv(final, file="Final_output_NetMHCPan.csv", quote=FALSE, row.names=FALSE)
