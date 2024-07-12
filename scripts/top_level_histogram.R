library(ggplot2)
profile.table = "F:/02.R_project/013.human_fisher_test/mgi_20240710/top_level/ramdom_p_0.05.txt"
group.file = "F:/02.R_project/013.human_fisher_test/mgi_20240710/top_level/g.txt"

X=read.table(profile.table,header=TRUE,row.names=1,sep="\t",check.names=F,quote="")
group=read.table(group.file,header=F,row.names=1,check.names=F,quote="")
X = X[,rownames(group)]

mp_level1 <- c("MP:0003631", "MP:0005385")
mp_name <- c("nervous system phenotype", 
             "cardiovascular system phenotype")

p_value <- numeric(length(mp_level1))
for (i in seq(length(mp_level1))) {
  mp <- X[(rownames(X) %in% mp_level1[i]), , drop = FALSE]
  target_number <- t(mp)[1,]
  sorted_numbers <- sort(as.numeric(mp))
  position <- match(target_number, sorted_numbers)
  percentage <- (position / length(sorted_numbers))
  p_value[i] = percentage
}


#### ===========================================================================
m <- 1
MP_0003631 <- X[(rownames(X) %in% mp_level1[m]), , drop = FALSE]
line <- t(MP_0003631)[1,]
data <- data.frame(
  value = t(MP_0003631)[-1,]
)
p1 <- ggplot(data, aes(x = value)) + #  binwidth = 0.5,bins=15,
  geom_histogram(aes(y = ..count.. / sum(..count..)), bins=20, fill = "black", color = "black") +
  theme_bw()+geom_vline(xintercept = line, color = "red", linetype = "solid", size = 0.5) +
  annotate("text", x = 0.5, y = 0.18, label = mp_name[m], 
           color = "black", size = 5, 
           fontface = "bold", family = "serif") +
  annotate("text", x = 0.5, y = 0.17, label = paste0("p = ", format(p_value[m], scientific = TRUE, digits = 3)), 
           color = "black", size = 5, family = "serif") +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        legend.title = element_blank()) +
  
  scale_y_continuous(limits = c(0,0.2),
                     expand = expansion(mult = c(0,0)),
                     breaks = seq(0,0.2,0.05))
p1

#### ===========================================================================
m <- 2
MP_0005385 <- X[(rownames(X) %in% mp_level1[m]), , drop = FALSE]
line <- t(MP_0005385)[1,]
data <- data.frame(
  value = t(MP_0005385)[-1,]
)
p2 <- ggplot(data, aes(x = value)) +
  geom_histogram(aes(y = ..count.. / sum(..count..)), bins=20, fill = "black", color = "black") +
  theme_bw()+geom_vline(xintercept = line, color = "red", linetype = "solid", size = 0.5) +
  annotate("text", x = 0.5, y = 0.18, label = mp_name[m], 
           color = "black", size = 5, 
           fontface = "bold", family = "serif") +
  annotate("text", x = 0.5, y = 0.17, label = paste0("p = ", format(p_value[m], scientific = TRUE, digits = 3)), 
           color = "black", size = 5, family = "serif") +
  theme(panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15),
        legend.title = element_blank()) +
  
  scale_y_continuous(limits = c(0,0.2),
                     expand = expansion(mult = c(0,0)),
                     breaks = seq(0,0.2,0.05))

p2