
library(qvalue)
a = "sifg4g_stat_exclude_gene.20230822.xls"
data_ <- read.table(a,header=TRUE,sep="\t",check.names=F,quote="")


for (i in 1:length(rownames(data_))){
  
  matrix_d <- matrix(c(data_[i,]$`disease_num`, 598 - data_[i,]$`disease_num`, 
                  data_[i,]$`control_num`, 532 - data_[i,]$`control_num`), nrow = 2)
  
  matrix_c <- matrix(c(data_[i,]$`control_num`, 532 - data_[i,]$`control_num`,
                       data_[i,]$`disease_num`, 598 - data_[i,]$`disease_num`), nrow = 2)
  
 
  case_binom.test <- binom.test(data_[i,]$`disease_num`,  598, p = data_[i,]$`org_prob`,
                                alternative = c("greater"), conf.level = 0.95)
  
  data_[i, 'case_binom.test'] <-  case_binom.test$p.value
  
  control_binom.test <- binom.test(data_[i,]$`control_num`, 532, 
                                   p = data_[i,]$`org_prob`, alternative = c("greater"),
                                   conf.level = 0.95)
  
  data_[i, 'control_binom.test'] <-  control_binom.test$p.value
}

new_data <- data_[data_[["mut_mun"]] >= 6,]
other_data <- data_[data_[["mut_mun"]] < 6,]


### binom.test
new_data[['case_binom.test.pvalues.adj']] <- 
  p.adjust(p=new_data[['case_binom.test']], method = "fdr", n = length(rownames(new_data)))

new_data[['control_binom.test.pvalues.adj']] <- 
  p.adjust(p=new_data[['control_binom.test']], method = "fdr", n = length(rownames(new_data)))


case_binom_qvobj <- qvalue(p = new_data[['case_binom.test']], pi0=NULL)
new_data[["case_binom.test.qvalues"]] <- case_binom_qvobj$qvalues


control_binom_qvobj <- qvalue(p = new_data[['control_binom.test']], pi0=NULL)
new_data[["control_binom.test.qvalues"]] <- control_binom_qvobj$qvalues


file_ = "sifg4g_stat_exclude_gene.20230822.add_qvalues.txt"
write.table(x = new_data, file = file_, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(x = other_data, file = file_, sep = "\t", quote = FALSE, append = TRUE, row.names = FALSE)
