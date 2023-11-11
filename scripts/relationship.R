setwd("/storage/shihongjunLab/liulifeng/project/03_lxx_220912_2000_mus_wgs/exon/3.analysis_raw_vcf/relationship_GCTA")
library(reshape2)
tmp <- read.table(gzfile("counts_le_03_mut.root.gcta.grm.gz"), header = F, stringsAsFactors = F)
ids <- read.table("counts_le_03_mut.root.gcta.grm.id", header = F, stringsAsFactors = F)
tmp <- tmp[,c(1,2,4)]
result_matrix <- acast(tmp, V1~V2, value.var="V4", drop = F)

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

result_full <- makeSymm(result_matrix)
diag(result_full) <- 2
result_df <- as.data.frame(result_full)
row.names(result_df) <- ids$V2
colnames(result_df) <- ids$V2
write.table(result_df, file = "gcta.kinship.txt", row.names = T, col.names = NA, sep = "\t", quote = F)

gcta<- read.table("gcta.kinship.txt")
library("pheatmap")
pdf("heapmap_of_gcta_kinship.pdf",width=10,height=10)
pheatmap(gcta, fontsize_row = 1, fontsize_col = 1, show_rownames=T,show_colnames=T)
dev.off()

