**same treatment for *A. thaliana*, *A. lyrata* and *C. rubella* all six tissues; *B. rapa* only for pollen and stigma**  
[Data](../blob/main/Data.md) 1. RNA_seq raw data  
Some scripts contain the analysis for all the data and some scripts only show the treatment for one.  

##### a. Use Trimmomatic to trim all the raw sequencing data

```
export trim_data=$HOME/dup/trim_data. 
parallel java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 $raw_data/Atha/AT_fastq/AT{}_R1.fq $raw_data/Atha/AT_fastq/AT{}_R2.fq $trim_data/Atha/AT{}_R1_paired.fq $trim_data/Atha/AT{}_R1_unpaired.fq $trim_data/Atha/AT{}_R2_paired.fq $trim_data/Atha/AT{}_R2_unpaired.fq ILLUMINACLIP:Trimmomatic-0.35/adapters/TruSeq2-PE.fa:2:30:10:1:true  LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50 ::: 1_leaf 2_flower 3_stem 4_root 5_silque
parallel java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 $raw_data/Atha/Ath{}_170628_NextSeq_R1.fastq $raw_data/Atha/Ath{}_170628_NextSeq_R2.fastq $trim_data/Atha/AT{}_R1_paired.fq $trim_data/Atha/AT{}_R1_unpaired.fq $trim_data/Atha/AT{}_R2_paired.fq $trim_data/Atha/AT{}_R2_unpaired.fq ILLUMINACLIP:Trimmomatic-0.35/adapters/TruSeq2-PE.fa:2:30:10:1:true  LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50 ::: PolA PolB PolC StigA StigB StigC
parallel java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 $raw_data/Brapa/Bra{}_150922_NextSeq_R1.fastq $raw_data/Brapa/Bra{}_150922_NextSeq_R2.fastq $trim_data/Brapa/BR{}_R1_paired.fq $trim_data/Brapa/BR{}_R1_unpaired.fq $trim_data/Brapa/BR{}_R2_paired.fq $trim_data/Brapa/BR{}_R2_unpaired.fq ILLUMINACLIP:Trimmomatic-0.35/adapters/TruSeq3-PE.fa:2:30:10:1:true  LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:50 ::: StigA StigB StigC
```	
##### b. Use seqkit check sequencing read overall status before and after trim 

```
seqkit stats dup/trim_data/Atha/*.fq > dup/report/Atha_trim.txt
seqkit stats dup/raw_data/Atha/AT_fastq/*.fq dup/raw_data/Atha/Ath*_170628_NextSeq_*.fastq > dup/report/Atha_raw.txt
```
##### c. Use fastqc check sequencing read quality before and after trim 

```
FastQC/fastqc -o dup/fastqc/trim -f fastq dup/trim_data/Atha/AT*_paired.fq 
FastQC/fastqc -o dup/fastqc/raw -f fastq dup/raw_data/Atha/Ath*_R*.fastq
```
 
##### d. Use STAR to construct the reference for mapping

```
ml star/2.6.1d
export ref_genome=$HOME/dup/ref/rna_seq
export index=$HOME/dup/index
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $index/Alyr_100 --genomeFastaFiles $ref_genome/Alyrata_384_v1.fa --sjdbGTFfile $ref_genome/Alyrata_384_v2.1.gene_exons.gff3 --sjdbGTFtagExonParentTranscript Parent
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $index/Atha_100 --genomeFastaFiles $ref_genome/Athaliana_447_TAIR10.fa --sjdbGTFfile $ref_genome/Athaliana_447_Araport11.gene_exons.gff3 --sjdbGTFtagExonParentTranscript Parent
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $index/Crub_100 --genomeFastaFiles $ref_genome/Crubella_474_v1.fa --sjdbGTFfile $ref_genome/Crubella_474_v1.1.gene_exons.gff3 --sjdbGTFtagExonParentTranscript Parent
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir $index/Brap_100 --genomeFastaFiles $ref_genome/BrapaFPsc_277_v1.fa --sjdbGTFfile $ref_genome/BrapaFPsc_277_v1.3.gene_exons.gff3 --sjdbGTFtagExonParentTranscript Parent
```
##### e. Use STAR to align the raw data to the index reference 

```
ml star/2.6.1d
index=$HOME/dup/index
trim=$HOME/dup/trim_data
parallel STAR --runThreadN 8 --runMode alignReads --genomeDir $index/Alyr_100 --readFilesIn $trim/Alyr/AL{}_R1_paired.fq $trim/Alyr/AL{}_R2_paired.fq --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 --outFileNamePrefix dup/star/Alyr/AL{} ::: 1_leaf 2_flower 3_stem 4_root 5_silque PA PB PC SA SB SC
parallel STAR --runThreadN 8 --runMode alignReads --genomeDir $index/Atha_100 --readFilesIn $trim/Atha/AT{}_R1_paired.fq $trim/Atha/AT{}_R2_paired.fq --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 --outFileNamePrefix dup/star/Atha/AT{} ::: 1_leaf 2_flower 3_stem 4_root 5_silque PolA PolB PolC StigA StigB StigC
parallel STAR --runThreadN 8 --runMode alignReads --genomeDir $index/Crub_100 --readFilesIn $trim/Crub/CR{}_R1_paired.fq $trim/Crub/CR{}_R2_paired.fq --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 --outFileNamePrefix dup/star/Crub/CR{} ::: 1_leaf 2_flower 3_stem 4_root 5_silque PA PB PC SA SB SC
parallel STAR --runThreadN 8 --runMode alignReads --genomeDir $index/Brap_100 --readFilesIn $trim_data/Brapa/BR{}_R1_paired.fq $trim_data/Brapa/BR{}_R2_paired.fq --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 8 --outFileNamePrefix dup/star/Brap/BR{} ::: PolA_150922 PolB_150922 PolC_150922
```

##### f. Use featurecount from subread to count the reads and generate the raw expression matrix

```
featurecount=$HOME/subread/bin
ref_genome=$HOME/dup/ref/rna_seq
atha=$HOME/dup/star/Atha
alyr=$HOME/dup/star/Alyr
crub=$HOME/dup/star/Crub
featurecount/featureCounts -T 8  -p -a $ref_genome/Athaliana_447_Araport11.gene_exons.gff3 -g Parent -o dup/count/Atha_counts.txt $atha/AT1_leafAligned.sortedByCoord.out.bam $atha/AT2_flowerAligned.sortedByCoord.out.bam $atha/AT3_stemAligned.sortedByCoord.out.bam $atha/AT4_rootAligned.sortedByCoord.out.bam $atha/ATPol?Aligned.sortedByCoord.out.bam $atha/ATStig?Aligned.sortedByCoord.out.bam
featurecount/featureCounts -T 8 -p -a $ref_genome/Alyrata_384_v2.1.gene_exons.gff3 -g Parent -o dup/count/Alyr_counts.txt $alyr/AL1_leafAligned.sortedByCoord.out.bam $alyr/AL2_flowerAligned.sortedByCoord.out.bam $alyr/AL3_stemAligned.sortedByCoord.out.bam $alyr/AL4_rootAligned.sortedByCoord.out.bam $alyr/ALP?Aligned.sortedByCoord.out.bam $alyr/ALS?Aligned.sortedByCoord.out.bam
featurecount/featureCounts -T 8 -p -a $ref_genome/Crubella_474_v1.1.gene_exons.gff3 -g Parent -o dup/count/Crub_counts.txt $crub/CR1_leafAligned.sortedByCoord.out.bam $crub/CR2_flowerAligned.sortedByCoord.out.bam $crub/CR3_stemAligned.sortedByCoord.out.bam $crub/CR4_rootAligned.sortedByCoord.out.bam $crub/CRP?Aligned.sortedByCoord.out.bam $crub/CRS?Aligned.sortedByCoord.out.bam
featurecount/featureCounts -T 8 -p -a $ref_genome/BrapaFPsc_277_v1.3.gene_exons.gff3 -g Parent -o dup/count/Brap_counts_2.txt $brap/BRPol?Aligned.sortedByCoord.out.bam $brap/BRStig?Aligned.sortedByCoord.out.bam $brap/BRPol?_150922Aligned.sortedByCoord.out.bam
```

##### g. Use R to generate the normalized TPM value and draw the heatmap for the studied homologs 
Normalize the raw expression value to TPM value 
R script
```
alyr <- read.table("Alyr_counts.txt",header=TRUE, sep="\t", row.names=1 )
alyr_metadata <- alyr[,1:5]
alyr_countdata <- alyr[,6:ncol(alyr)]
colnames(alyr_countdata) <- c('AL1_leaf','AL2_flower', 'AL3_stem', 'AL4_root', 'ALPA', 'ALPB', 'ALPC', 'ALSA', 'ALSB', 'ALSC')
prefix<-"couts"
alyr_kb <- alyr_metadata$Length / 1000
alyr_rpk <- alyr_countdata / alyr_kb
alyr_tpm <- t(t(alyr_rpk)/colSums(alyr_rpk) * 1000000)
write.csv(alyr_tpm,paste0(prefix,"_tpm.csv"))

```
Draw heatmap.  

```
library(devtools)
library(ComplexHeatmap)
library(circlize)
library(gridtext)
library(dplyr)
expression_raw <- read.table("heatmap_expression.matrix")
colnames(expression_raw) <- c('Gene ID','leaf','flower','stem','root','Pollen_1', 'Pollen_2', 'Pollen_3', 'Stigma_1', 'Stigma_2', 'Stigma_3')
expression_raw <- mutate(expression_raw,pollen=(Pollen_1 + Pollen_2 + Pollen_3)/3)
expression_raw <- mutate(expression_raw,Stigma=(Stigma_1 + Stigma_2 + Stigma_3)/3)
expression <- select(expression_raw,-Pollen_1, -Pollen_2, -Pollen_3, -Stigma_1, -Stigma_2, -Stigma_3)
row.names(expression) <- expression[['Gene ID']]
expression <- expression[2:7]
logplusone <- function(x) {log10(x[1] + 1)}
expression[,] <- as.data.frame(lapply(expression[,], FUN = function(x) {sapply(x, FUN = logplusone)}))
col_fun <- colorRamp2(seq(0, 4, length=4), rev(rainbow(4)), space="RGB")
expression <- as.matrix(expression)
expression.ht <- Heatmap(expression, name = "log10(TPM value + 1)", col = col_fun, rect_gp = gpar(col = "white", lwd = 2), cluster_rows = TRUE,cluster_columns = FALSE, 
                 row_names_side = "left", column_names_side = "bottom",column_names_gp = gpar(fontsize = 12), row_dend_width = unit(2, "cm"))
expression.ht
library(GetoptLong)
pdf(qq("dup_heatmap_3.pdf"), width = 10, height = 10)
draw(expression.ht, column_title_gp=gpar(fontsize = 16))
dev.off()
```

##### h. Principal component analysis (PCA) 
R script
```
library(ggplot2)
library(ggrepel)
expression_pca <- select(expression_raw,-Pollen_1, -Pollen_2, -Pollen_3, -Stigma_1, -Stigma_2, -Stigma_3)
row.names(expression_pca) <- expression_pca[['Gene ID']]
expression_pca <- expression_pca[2:7]
head(expression_pca)
pca1 <- prcomp(expression_pca[],center = TRUE,scale. = TRUE)
df1 <- pca1$x 
df1 <- as.data.frame(df1)
summ1 <- summary(pca1)
summ1
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = rownames(expression_pca)))+ 
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c('orange', 'purple','blue','black','green','yellow','red'))+
  scale_colour_manual(values = c('orange', 'purple','blue','black','green','yellow','red'))+
  theme(axis.text = element_text(size = 14, face='bold'), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca1
filename.pca.4="dup_pca_4.pdf"
ggsave(file=filename.pca.4, p.pca1, width=10, height=10, units="in")
```

##### i. Extracte the normalized expression data of the two studied *A. thaliana* duplicate genes *AT5G12950* and *AT5G12960*, which contain 2,843 samples from the previous 106 RNA-seq studies and were analyzed by Genevestigator.   
R script
Draw the heatmap.

```
expression_vestigator <- read.csv("genevestgator.csv", header = T)
row.names(expression_vestigator) <- expression_vestigator[['Sample']]
expression_vestigator <- select(expression_vestigator,AT5G12950, AT5G12960)
logplusone <- function(x) {log10(x[1] + 1)}
expression_vestigator[,] <- as.data.frame(lapply(expression_vestigator[,], FUN = function(x) {sapply(x, FUN = logplusone)}))
summary(expression_vestigator)
col_fun <- colorRamp2(seq(0, 4, length=4), rev(rainbow(4)), space="RGB")
expression_vestigator <- as.matrix(expression_vestigator)

pollen <- read.csv("pollen_id.csv", header = F)
pollen_l = rownames(expression_vestigator) %in% pollen$V1
summary(pollen_l)
expression_vestigator.list <- Heatmap(expression_vestigator,name = 'gene expression', col = col_fun, cluster_rows = TRUE,cluster_columns = FALSE, 
                                    show_row_names = FALSE,show_column_names = FALSE, use_raster = TRUE, raster_by_magick = TRUE) + 
                              Heatmap(pollen_l + 0, name = "pollen", col = c("0" = "black", "1" = "Orange"), 
                                      show_heatmap_legend = TRUE, width = unit(5, "mm"),show_column_names = FALSE)
expression_vestigator.list.ht = draw(expression_vestigator.list, main_heatmap = "gene expression")
expression_vestigator.list.ht
library(GetoptLong)
library(ggplot2)
filename.expression_vestigator.ht="dup_heatmap_vestigator_1.pdf"
ggsave(file=filename.expression_vestigator.ht, expression_vestigator.list.ht, width=10, height=10, units="in")
```

##### i. Extracte the normalized expression data of the two studied *A. thaliana* duplicate genes *AT5G12950* and *AT5G12960*, which involved within 79 organs and developmental stages from [Trava](http://travadb.org/).   
R script
Draw the heatmap.

```
trava_raw <- read.csv("trava.csv", header = T)
row.names(trava_raw) <- trava_raw[['Sample']]
trava <- select(trava_raw,-Sample)
row.names(trava) <- trava_raw[['Sample']]
logplusone <- function(x) {log2(x[1] + 1)}
trava[,] <- as.data.frame(lapply(trava[,], FUN = function(x) {sapply(x, FUN = logplusone)}))
col_fun <- colorRamp2(seq(0, 16, length=4), rev(rainbow(4)), space="RGB")
trava <- as.matrix(trava)
trava.ht <- Heatmap(trava, name = "TMM value", col = col_fun, rect_gp = gpar(col = "white", lwd = 2), cluster_rows = TRUE,cluster_columns = FALSE, 
                         row_names_side = "right", show_column_names = FALSE, row_names_gp = gpar(fontsize = 8),column_names_side = "bottom")
trava.ht
pdf(qq("dup_heatmap_trava.pdf"), width = 10, height = 10)
draw(trava.ht, column_title_gp=gpar(fontsize = 16))
```

##### k. Download the Proteomic data, do the statistical analysis and draw the plot
Use [ATHENA](http://athena.proteomics.wzw.tum.de/) to obtain the protein expression of AT5G12950 and AT5G12960 among 30 tissues.
Draw the plot and do the statistical analysis via GraphPad Prism. 
