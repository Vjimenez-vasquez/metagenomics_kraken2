# metagenomics_kraken2
a collection of bash and R codes for metagenomics analysis using kraken2, bracken and ggplot
designed by Victor Jimenez-Vasquez (vr.jimenez.vs@gmail.com)
# Kraken2 databases availabe at : https://benlangmead.github.io/aws-indexes/k2
```

![metagenomics_kraken2](https://github.com/Vjimenez-vasquez/metagenomics_kraken2/assets/89874227/3cfc3e29-98cd-40ca-8124-b51606ba3c2e)

# step 1 : The bash code for INS
```r
# path to BRACKEN: /media/ins-bio/DATA01/data_base_download/Bracken-2.7/./bracken
# path to KRAKEN VIRUS DATABASE: /media/ins-bio/DATA01/data_base_download/KRAKENVIRDB

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
