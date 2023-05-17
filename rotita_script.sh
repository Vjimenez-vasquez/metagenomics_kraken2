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
kraken2 --paired --use-names --gzip-compressed --db /home/administrador/Documentos/KRAKENVIRDB/ --threads 28 $r1 $r2 --report ${prefix}_report.txt --output ${prefix}_kraken2.out ;
done ;
rm *.fastq.gz_report.txt ; 
mkdir kraken_out ;
mv *.out  kraken_out/ ;
mkdir kraken_txt ; 
mv *.txt kraken_txt/ ;  
cd kraken_txt/ ; 
ls -lh ; 

# 2.1 # bracken #
for r1 in *_report.txt
do
prefix=$(basename $r1 _report.txt)
/home/administrador/Documentos/rotita_170523/Bracken-master/./bracken -d /home/administrador/Documentos/KRAKENVIRDB/ -i $r1 -o ${prefix}.bracken.txt -l S ; 
done ; 
mkdir species_report ; 
mv *_species.txt species_report/ ;
mkdir bracken_species_abundances ;  
mv /home/administrador/Documentos/rotita_170523/Bracken-master/*.txt bracken_species_abundances/ ; 
mv *.bracken.txt bracken_species_abundances/ ;
ls ; 


# 3 # PAVIAN #
if (!require(remotes)) { install.packages("remotes") }
remotes::install_github("fbreitwieser/pavian")
pavian::runApp(port=5000)

# 4 # trimming #
for r1 in *fastq.gz
do
prefix=$(basename $r1 _L001_R1_001.fastq.gz)
r2=${prefix}_L001_R2_001.fastq.gz
java -jar trimmomatic-0.39.jar PE -threads 25 $r1 $r2 ${prefix}_f_paired.fq.gz ${prefix}_f_unpaired.fq.gz ${prefix}_r_paired.fq.gz ${prefix}_r_unpaired.fq.gz SLIDINGWINDOW:4:20 MINLEN:35 ;
done ; 

mkdir trimm ; 
mv *_paired.fq.gz trimm/ ; 
rm *_unpaired.fq.gz ;
cd trimm/ ;
fastqc *.gz -t 15 ; 
mkdir fastqc ;
mv *.html *.zip fastqc/ ; 
ls -lh ;

#5 # de-novo assembly #
for r1 in *fq.gz
do
prefix=$(basename $r1 _f_paired.fq.gz)
r2=${prefix}_r_paired.fq.gz
metaspades --pe1-1 $r1 --pe1-2 $r2 -t 25 -o ${prefix}_spades ;
#spades --pe1-1 $r1 --pe1-2 $r2 --careful -t 25 -o ${prefix}_spades ;#
mv ${prefix}_spades/scaffolds.fasta ${prefix}_spades/${prefix}_spades_scaffolds.fasta ;
cp ${prefix}_spades/${prefix}_spades_scaffolds.fasta . ;
done ;
rm -r *fq.gz_spades ;
ls -lh ; 

# install.packages("seqinr") #
library(seqinr) ;
r <- dir() ;
head <- gsub("_spades_scaffolds.fasta","",r) ;
a <- 0 ;
for (i in 1:length(head)){ ;
  a <- read.fasta(r[i]) ;
  names(a) <- paste0(rep(head[i],length(a)),"_",names(a)) ;
  write.fasta(a,names(a), file.out=paste0(head[i],".fas")) ;
} ;
q("no") ;

cat *.fas > contigs.fasta

cd .. ;
mkdir scaffolds ; 
mv *scaffolds.fasta script_2.R scaffolds/ ;
cd scaffolds/ ;

#6 # looping mapping and estimate abundance #
for r1 in *fq.gz
do
prefix=$(basename $r1 _f_paired.fq.gz)
r2=${prefix}_r_paired.fq.gz
r3=${prefix}.fas

bwa index $r3 ;
bwa mem -t 25 $r3 $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 25 -bS -T $r3 ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 25 -n ${prefix}_unoa.bam -o ${prefix}_count.bam ;
samtools index -@ 25 ${prefix}_count.bam ;

rm ${prefix}_uno.sam ${prefix}_unoa.bam ;
done ;

for p1 in *_count.bam
do
prefix=$(basename $p1 _count.bam)
samtools view -@ 25 $p1 | cut -f1,3 | sort | uniq | cut -f2 | sort | uniq -c > ${prefix}_readcounts.txt ;
done ;
rm *.amb *.ann *.bwt *.fai *.pac *.sa ;
cat *_readcounts.txt > counts.txt ;
grep "NODE" counts.txt > counts2.txt ;
ls -lh ; 

# 7 # next #
## download DB : https://ccb.jhu.edu/software/kraken/ ##
kraken --db /home/administrador/Documentos/kraken/minikraken_20171019_8GB/ --fasta-input all.fasta --threads 25 --unclassified-out unclassified --classified-out classified --output iih_files ; 
kraken-translate --db /home/administrador/Documentos/kraken/minikraken_20171019_8GB/ iih_files > iih_files.labels.csv ; 

# 8 # CAT/BAT #
CAT prepare --fresh -d /home/administrador/Documentos/CAT_prepare_20210107/2021-01-07_CAT_database2/ -t /home/administrador/Documentos/CAT_prepare_20210107/2021-01-07_taxonomy2/ ; 
CAT contigs -c all.fasta -d /home/administrador/Documentos/CAT_prepare_20210107/2021-01-07_CAT_database -t /home/administrador/Documentos/CAT_prepare_20210107/2021-01-07_taxonomy -o all2 ;


