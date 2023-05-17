setwd("/home/vjimenez/Documentos/rotita_170523/bracken")
dir()

a <- read.csv("INS030506723_S1.bracken.txt", header=T, sep="\t")
a1 <- a[order(-a$fraction_total_reads), ]
head(a1)
a2 <- data.frame(sample=rep("A",length(a$name)), a1)
head(a2)

library(tidyr)
a3 <- separate(a2,"sample",c("sample","delete"), sep="_")
head(a3)


## reading the bracken output ## 
r <- dir()
r1 <- data.frame()
r2 <- 0
r3 <- 0 
r4 <- 0 
r5 <- 0 
r6 <- 0 
r7 <- 0 

for (i in 1:length(r)){
r1 <- read.csv(r[i], header=T, sep="\t")
#r2 <- rep(gsub("\\_S.*","",r[i],fixed=F),length(r1$name))#
r2 <- rep(gsub(".bracken.txt","",r[i],fixed=F),length(r1$name))
r3 <- append(r2,r3)
r4 <- append(r1$name, r4)
r5 <- append(r1$new_est_reads, r5)
r6 <- append(r1$fraction_total_reads, r6)
}

unique(r3)

r8 <- data.frame(sample=r3,species=r4,reads=r5,percentage=r6)
head(r8)
write.csv(r8, "Viroma_total.csv",row.names=F)

## the plot for abundances ##

abundances <- function(data,percentage,title,level){
data <- data
percentage <- as.numeric(percentage)
title <- title
level <- level
percentage2 <- as.numeric(percentage)*100

r9 <- data[data$percentage >= percentage , ]
print(r9)

library(ggplot2)
library(plotly)
pt <- ggplot(r9,aes(x=r9$sample, y=r9$percentage, fill=r9$species)) + 
  geom_bar(stat="identity") + xlab("sample") + ylab("Porcentaje") + theme_minimal() + theme(legend.position = 'bottom') + 
  geom_text(aes(label = paste0(r9$species,":",r9$reads)), position = position_stack(vjust = 0.5), colour = "black", size = 2) + 
  scale_fill_discrete(name = paste0(title,":",level)) + theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle(paste0(title,"_at_",level,"_level_",as.character(percentage2),"%_abundance"))

print(pt)
ggplotly(pt)

write.csv(r9, paste0(title,"_",level,"_",as.character(percentage),".csv"),row.names=F)
}

## test ## 

pdf(file ="Viroma_pos_rotavirusA.pdf",width = 35, height = 15)
abundances(data=r8,percentage="0.01",title="Virome",level="species")
dev.off()

pos <- unique(r8[r8$species %in% "Rotavirus A", 1])
length(pos)

unique(r8$sample)
r10 <- r8[r8$sample %in% pos, ]
dim(r10)
dim(r8)

pdf(file ="Viroma_pos_rotavirusA.pdf",width = 35, height = 15)
abundances(data=r10,percentage="0.05",title="Virome_positives_rotA",level="species")
dev.off()


r11 <- r8[r8$sample %in% c("INS030518523_S4","INS030518223_S3","INS030517723_S2","INS030517423_S1"), ]
abundances(data=r11,percentage="0.05",title="test",level="species")
