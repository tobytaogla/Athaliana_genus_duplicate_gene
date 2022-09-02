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

##### g. Use R command to generate the normalized TPM value and draw the heatmap for the studied homologs 
Normalize the raw expression value to TPM value 

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


##### h. Extracte the normalized expression data of the two studied *A. thaliana* duplicate genes *AT5G12950* and *AT5G12960*, which contain 2,843 samples from the previous 106 RNA-seq studies and werer analyzed by Genevestigator.   
Draw the heatmap.


##### i. Download the Proteomic data, do the statistical analysis and draw the plot
Use [ATHENA](http://athena.proteomics.wzw.tum.de/) to obtain the protein expression of AT5G12950 and AT5G12960 among 30 tissues.
Draw the plot and do the statistical analysis via GraphPad Prism. 

