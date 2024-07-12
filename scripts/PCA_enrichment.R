library(DescTools)
library("ade4")

profile.table = "F:/02.R_project/013.human_fisher_test/mgi_20240710/last_level/random_p"
group.file = "F:/02.R_project/013.human_fisher_test/mgi_20240710/g.txt"
source('F:/02.R_project/013.human_fisher_test/mgi_20240710/labels2colors.R')

X=read.table(profile.table,header=TRUE,row.names=1,sep="\t",check.names=F,quote="")
group=read.table(group.file,header=F,row.names=1,check.names=F,quote="")
W = X[,rownames(group)]
X <- subset(W, mouse_candidate_genes <0.01)

x <-  as.data.frame(lapply(X, function(x) ifelse(x<0.01, 0, 1)), stringsAsFactors = FALSE)     
rownames(x)  <- rownames(X) 
X <- x

color_list = group2corlor(group)
sample_colors = color_list[[1]]
group_colors  = color_list[[2]]
group_names = color_list[[3]]
group = color_list[[4]]

Xdist=dist(t(X))
X.dudi=dudi.pca(t(X),center=T,scale=T,scan=F)

case <- as.numeric(X.dudi$li[1,])
random <- as.matrix(X.dudi$li[-1,])
result <- HotellingsT2Test(random, mu = case)
print(result)  # output p-value

len=c()
con=X.dudi$eig/sum(X.dudi$eig)*100
con=round(con,2)
par(mar=c(4.1,5.1,4.1,2.1))
tem = 1
xli=c(min(X.dudi$li[1])-tem,max(X.dudi$li[1])+tem)
yli=c(min(X.dudi$li[2])-tem,max(X.dudi$li[2])+tem)

plot(X.dudi$li,col=sample_colors,pch=19,
     xlab=paste("PCA1(",con[1],"%)",sep=""),
     ylab=paste("PCA2(",con[2],"%)",sep=""),
     cex=0.5,cex.axis=1.5,cex.lab=1.5,xlim = xli,ylim = yli
)

c <- c('case-enriched geneset', '10,000 random genesets')
legend("topright", legend = c, col = group_colors, pch=19, cex = 0.9)
legend("bottomright",paste0("p < 2.2e-16"))
box(which = "figure",col = "white",lwd = 6)


# enrichment annalysis
p_all_file <-  "F:/02.R_project/013.human_fisher_test/mgi_20240710/last_level/p_mouse_candidate_genes_rich.txt"
p_all <- read.table(p_all_file,header=TRUE,row.names=1,sep="\t",check.names=F,quote="")
p_all <- subset(p_all, Pvalue <0.01)
# 
p_all <- p_all[order(-p_all$Pvalue), ]
p_all <- tail(p_all, n=20)
p_all$des <- factor(p_all$des, levels = p_all$des)


ggplot(p_all, aes(x=rich, y=des)) +
  geom_point(aes( colour = -log10( Pvalue ), size= count)  ) +
  # geom_hline(yintercept =  c(1, 2, 3, 4), linetype = "dashed", color = "#CCE4F7")+
  scale_y_discrete(limits=p_all$des)+
  labs(title = "Top 20 Mammalian Phenotype Ontology", 
       colour = "-log10(pvalue)",
       size = "Gene Number")+
  xlab("Rich Factor")+
  # scale_color_gradient(low = "#FEE5D9", high = "#894C38") + 
  scale_color_gradient(low = 'blue', high = 'red') + 
  xlim(range(p_all$rich)) +
  scale_size_continuous(range = c(2, 6)) +
  theme(plot.title = element_text(family = "serif", face = "bold", color = "black", size = 16, hjust = 0.5),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12, family = "serif"),
        axis.title.x = element_text(size = 12, family = "serif", face = "bold"),
        axis.text.x = element_text(size = 10, family = "serif"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        legend.position = "right",
        legend.title = element_text(size = 10, family = "serif", face = "bold", hjust = 0.2),
        legend.text = element_text(hjust = 0.8),
        # legend.key.size = unit(1.5, "lines"),
        legend.box.just = "center", 
        
        panel.background = element_blank(), 
        panel.grid.major = element_line(size = 0.05, colour = "gray"), 
        # panel.grid.minor.y = element_line(size=0.05, colour="gray"), 
        # panel.grid.minor.x = element_line(size=0.05, colour="gray")
        # plot.background = element_blank()  
  ) +
  guides(color = guide_colorbar(order = 2), size = guide_legend(order = 1))
  