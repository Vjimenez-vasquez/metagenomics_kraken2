# metagenomics_kraken2
a collection of bash and R codes for metagenomics analysis using kraken2, bracken and ggplot
designed by Victor Jimenez-Vasquez (vr.jimenez.vs@gmail.com)
# Kraken2 databases availabe at : https://benlangmead.github.io/aws-indexes/k2
![metagenomics_kraken2](https://github.com/Vjimenez-vasquez/metagenomics_kraken2/assets/89874227/3cfc3e29-98cd-40ca-8124-b51606ba3c2e)
# Bracken is availabe at :  https://github.com/jenniferlu717/Bracken
 
# step 1 : The bash code for INS
```r
# path to BRACKEN: /media/ins-bio/DATA01/data_base_download/Bracken-2.7/./bracken
# path to KRAKEN VIRUS DATABASE: /media/ins-bio/DATA01/data_base_download/KRAKENVIRDB

# path to KRAKEN DATABASE (104) :/home/administrador/Documentos/KRAKENPlusDB
# path to KRAKENDB (104) :/home/administrador/Documentos/KRKDB2
# path to BRACKEN (104) : /home/administrador/Documentos/Bracken-2.7

# 1 # fastqc #
fastqc -t 25 *
mkdir fastqc ; 
mv *.html *.zip fastqc/ ; 
ls -lh ; 

# 2 # kraken viral #
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz
kraken2 --paired --use-names --gzip-compressed --db /media/ins-bio/DATA01/data_base_download/KRAKENVIRDB/ --threads 28 $r1 $r2 --report ${prefix}_report.txt --output ${prefix}_kraken2.out ;
done ;
rm *.fastq.gz_report.txt ; 
mkdir kraken_out ;
mv *.out  kraken_out/ ;
mkdir kraken_txt ; 
mv *.txt kraken_txt/ ;  
cd kraken_txt/ ; 
ls -lh ; 

# 3 # bracken #
for r1 in *_report.txt
do
prefix=$(basename $r1 _report.txt)
/media/ins-bio/DATA01/data_base_download/Bracken-2.7/./bracken -d /media/ins-bio/DATA01/data_base_download/KRAKENVIRDB/ -i $r1 -o ${prefix}.bracken.txt -l S ; 
done ; 
mkdir species_report ; 
mv *_species.txt species_report/ ;
mkdir bracken_species_abundances ;  
mv /media/ins-bio/DATA01/data_base_download/Bracken-2.7/*.txt bracken_species_abundances/ ; 
mv *.bracken.txt bracken_species_abundances/ ;
ls ;

# 4 # PAVIAN #
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port=5000)
```

# step 2 : R code for plot
```r
# 1 # set the working directory containing bracken output files (ej: INS030518523_S4.bracken.txt) #

setwd("/home/vjimenez/Documentos/rotita_170523/bracken")
dir()
library(tidyr)
#install.packages("tidyr")#

# 2 # reading the bracken output ##

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
#considering "INS030518523_S4.bracken.txt" as an example name#
#r2 <- rep(gsub("\\_S.*","",r[i],fixed=F),length(r1$name))#
r2 <- rep(gsub(".bracken.txt","",r[i],fixed=F),length(r1$name))
r3 <- append(r2,r3)
r4 <- append(r1$name, r4)
r5 <- append(r1$new_est_reads, r5)
r6 <- append(r1$fraction_total_reads, r6)
}

#check the items# 
unique(r3)

r8 <- data.frame(sample=r3,species=r4,reads=r5,percentage=r6)
head(r8)
write.csv(r8, "Viroma_total.csv",row.names=F)

# 3 # the plot for abundances ##

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

write.csv(r9, paste0(title,"_",level,"_",as.character(percentage),".csv"),row.names=F)

pu <- ggplotly(pt)
print(pt)
pu

}

# 4 : test the code ## 

pdf(file ="Viroma_1.pdf",width = 35, height = 15)
abundances(data=r8,percentage="0.01",title="Virome",level="species")
dev.off()
```

# step 3 : USAGE #
```r
# data : abundance data frame obtained in "step-2"
# percentage : minimun abundance percentage (in frequence units) to consider ej. 1% = 0.01 , 0.1% = 0.001 , 50% = 0.5
# title : a given prefix to include in the title. ej: if you use "Virome" word, the final title will include the "level" and the "percentage" to obtain the final title : "Virome_at_species_level_1%_abundance"
# level : taxonomic level of the identification. ej: "species"
```

## REMOVING "BeAn 58058 virus" ##
```r
head(r8)
t <- r8[!r8$species %in% "BeAn 58058 virus", ]
head(t)
dim(r8)

t1 <- 0 ; t2 <- 0 ; t3 <- 0 ; t4 <- 0 ;
for (i in unique(t$sample)){

t1 <- as.numeric(t[t$sample %in% i , 3])
t2 <- sum(t[t$sample %in% i , 3])
t3 <- t1/t2
t4 <- append(t4,t3)
}

t$percentage <- t4[2:length(t4)]
head(t)

pdf(file ="Viroma_1.pdf",width = 35, height = 15)
abundances(data=r8,percentage="0.01",title="Virome (suero)",level="species")
dev.off()

pdf(file ="Viroma_2_clean.pdf",width = 35, height = 15)
abundances(data=t,percentage="0.01",title="Virome (suero)",level="species")
dev.off()
```

# step 4 : HEATMAP #
```r
#install.packages("fossil")#
#install.packages("tidyr")#
#install.packages("ggplot2")#

library(ggplot2)
library(fossil)
library(tidyr)

data <- read.csv("Viroma_total.csv", header=TRUE)
dim(data)
head(data)
sort(unique(data$species))

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
data <- data[data$reads > 1000, ]
heatmap_plot1 <- ggplot(data, aes(x = species, y = sample, fill = percentage)) +
		geom_tile(color = "black") +
  		geom_text(aes(label = round(reads,2)), color = "#333333", size = 1.2) +
		scale_fill_gradient2(low = "gray", high = "blue",
    	      mid="red", midpoint = (max(data$percentage)/2), limits = c(0, max(data$percentage))) +
		theme_bw() + 
		coord_fixed() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		scale_x_discrete(name="Microbial species") +
	      scale_y_discrete(name="samples") +
		ggtitle("Microbial Diversity")
heatmap_plot1
```

# step 5 : HEATMAP (BeAn removed) #
```r
#install.packages("fossil")#
#install.packages("tidyr")#
#install.packages("ggplot2")#

library(ggplot2)
library(fossil)
library(tidyr)

data <- read.csv("Virome (suero)_species_0.01.csv", header=TRUE)
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
data3 <- data2[data2$reads > 200, ]
heatmap_plot1 <- ggplot(data3, aes(x = reorder(species,-abundance), y = sample, fill = percentage)) +
		geom_tile(color = "black") +
  		geom_text(aes(label = round(reads,2)), color = "#333333", size = 1) +
		scale_fill_gradient2(low = "gray", high = "blue",
    	      mid="red", midpoint = (max(data3$percentage)/2), limits = c(0, max(data3$percentage))) +
		theme_bw() + 
		coord_fixed() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		scale_x_discrete(name="Viral species (> 200 reads)") +
	      scale_y_discrete(name="samples") +
		ggtitle("Viral Diversity (suero)")
heatmap_plot1

heatmap_plot2 <- ggplotly(heatmap_plot1)
print(heatmap_plot1)
heatmap_plot2
```
# step 6 : METAGENOMICS_COMMANDS.R examples #
```r
## directories ##

# setwd("C:/Users/USUARIO/Documents/cietrop/analisis_2/hisopado") #
# setwd("C:/Users/USUARIO/Documents/cietrop/analisis_2/suero") #

dir()

source("C:/Users/USUARIO/Documents/cietrop/analisis_2/metagenomic_commands.R")

## TOTAL DATA ##
data1 <- freqs("Viroma_total", "all","0")
head(data1)
dim(data1)

pdf(file ="Viroma_total.pdf",width = 35, height = 15)
abundances(data=data1,percentage="1",title="Viroma_total",level="species")
dev.off()

pdf(file ="Viroma_clean.pdf",width = 35, height = 15)
abundances_sp(data=data1,percentage="1",title="Viroma_clean",level="species")
dev.off()

## REMOVING BeAn 58058 virus ## 
data2 <- freqs("Viroma_clean", "BeAn 58058 virus", "0")
head(data2)

pdf(file ="Viroma_clean.pdf",width = 35, height = 15)
abundances(data=data2,percentage="0",title="Viroma_clean",level="species")
dev.off()

pdf(file ="Viroma_clean.pdf",width = 35, height = 15)
abundances_sp(data=data2,percentage="0",title="Viroma_clean",level="species")
dev.off()

## REMOVING BeAn 58058 virus and species > 10% abundance ## 

pdf(file ="Viroma_clean_2.pdf",width = 35, height = 15)
abundances(data=data2,percentage="30",title="Viroma_clean",level="species")
dev.off()

pdf(file ="Viroma_clean_2.pdf",width = 35, height = 15)
abundances_sp(data=data2,percentage="30",title="Viroma_clean",level="species")
dev.off()

## HEATMAP ## 
metaheat(data2,1000,"Virome 1000 reads (hisopado)")
```
