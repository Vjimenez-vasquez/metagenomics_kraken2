library(tidyr)
library(ggplot2)
library(plotly)
library(fossil)
library(tidyr)

### FUNCTION 1 : ESTIMATE FREQUENCIES AND PERCENTAGES ##

freqs <- function(title,remove,frequency){

title <- title
remove <- remove
frequency <- as.numeric(frequency)
r <- dir(); r1 <- data.frame() ; r2 <- 0 ; r3 <- 0 ; r4 <- 0 ; r5 <- 0 ; r6 <- 0 ; r7 <- 0 ; 

for (i in 1:length(r)){
r1 <- read.csv(r[i], header=T, sep="\t")
#considering "INS030518523_S4.bracken.txt" as an example name#
#r2 <- rep(gsub("\\_S.*","",r[i],fixed=F),length(r1$name))#
r2 <- rep(gsub(".bracken.txt","",r[i],fixed=F),length(r1$name))
r3 <- append(r2,r3)
r4 <- append(r1$name, r4)
r5 <- append(r1$new_est_reads, r5)
r6 <- append(r1$fraction_total_reads, r6)
}

#check the items# 
r8 <- data.frame(sample=r3,species=r4,reads=r5,frequency=r6)
r8 <- r8[!r8$species %in% remove & r8$frequency > frequency , ]

#estimate percentages#
s1 <- 0 ; s2 <- 0 ; s3 <- 0 ; s4 <- 0 ; 
for (i in unique(r8$sample)){
s1 <- r8[r8$sample %in% i, 3]
s2 <- sum(r8[r8$sample %in% i, 3])
s3 <- (s1/s2)*100
s4 <- append(s4,s3)
}

s5 <- s4[2:length(s4)]
r8$percentage <- s5
head(r8)

write.csv(r8,paste0(title,".csv"),row.names=F)

r8 <- return(r8)
}


### FUNCTION 2 : BARPLOTS ##

abundances <- function(data,percentage,title,level){
data <- data
percentage <- as.numeric(percentage)
title <- title
level <- level
percentage2 <- as.numeric(percentage)

r9 <- data[data$percentage >= percentage , ]
print(r9)

pt <- ggplot(r9,aes(x=r9$sample, y=r9$percentage, fill=r9$species)) + 
  geom_bar(stat="identity") + xlab("sample") + ylab("Porcentaje") + theme_minimal() + theme(legend.position = 'bottom') + 
  geom_text(aes(label = paste0(r9$species,":",r9$reads)), position = position_stack(vjust = 0.5), colour = "black", size = 2) + 
  scale_fill_discrete(name = paste0(title,":",level)) + theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle(paste0(title,"_at_",level,"_level_>_",as.character(percentage2),"%_(abundance)"))

write.csv(r9, paste0(title,"_",level,"_",as.character(percentage),".csv"),row.names=F)

pu <- ggplotly(pt)
print(pt)
pu

}

abundances_sp <- function(data,percentage,title,level){
data <- data
percentage <- as.numeric(percentage)
title <- title
level <- level
percentage2 <- as.numeric(percentage)

r9 <- data[data$percentage >= percentage , ]
print(r9)

pv <- ggplot(r9,aes(x=reorder(r9$species,r9$percentage), y=r9$percentage, fill=r9$sample)) + 
  geom_bar(stat="identity") + xlab("species") + ylab("Porcentaje") + theme_minimal() + theme(legend.position = 'bottom') + 
  geom_text(aes(label = paste0(r9$sample,":",r9$reads)), position = position_stack(vjust = 0.5), colour = "black", size = 2) + 
  scale_fill_discrete(name = paste0(title,":",level)) + theme(axis.text.x = element_text(angle = 45)) + 
  ggtitle(paste0(title,"_at_",level,"_level_>_",as.character(percentage2),"%_(abundance)"))

pw <- ggplotly(pv)
print(pv)
pw

}

### FUNCTION 3 : HEATMAPS ##
metaheat <- function(data,reads,title){

data <- data 
reads <- reads
title <- title

sort(unique(data$species))
sp_reads <- aggregate(data$reads, by=list(data$species), FUN=sum)
sp_reads <- sp_reads[order(sp_reads$x),]
names(sp_reads) <- c("species","abundance")
data2 <- merge(data,sp_reads, by="species", all.x=F)
dim(data)
dim(data2)
head(data2)

q <- data.frame(sample=data$sample, species=data$species, abundance=as.numeric(data$percentage)) 
head(q)
dim(q)
df <- create.matrix(q, tax.name = "sample",locality = "species", abund.col = "abundance", abund = TRUE)
class(df)
table <- as.data.frame(df)
dim(table)
head(table)

write.csv(table,"complete_table.csv", row.names=T, quot=F)

head(data)
data3 <- data2[data2$reads > reads, ]
heatmap_plot1 <- ggplot(data3, aes(x = reorder(species,-abundance), y = sample, fill = percentage)) +
		geom_tile(color = "black") +
  		geom_text(aes(label = round(reads,2)), color = "#333333", size = 1) +
		scale_fill_gradient2(low = "gray", high = "blue",
    	      mid="red", midpoint = (max(data3$percentage)/2), limits = c(0, max(data3$percentage))) +
		theme_bw() + 
		coord_fixed() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		scale_x_discrete(name=paste0("Viral species (> ",reads,"reads)")) +
	      scale_y_discrete(name="samples") +
		ggtitle(title)
print(heatmap_plot1)

heatmap_plot2 <- ggplotly(heatmap_plot1)
print(heatmap_plot1)
heatmap_plot2

}

