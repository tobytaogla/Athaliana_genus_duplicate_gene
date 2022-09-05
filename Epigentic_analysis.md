##### a. Download the raw sequencing data as below described.
Raw data downloaded from the [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra)  

Chip-Seq data:  
H3K27ac  (leaf SRR1509474, floral bud SRR1509479). 
H3K4me1  (leaf SRR1509477, floral bud SRR1509476). 
H3K27me3 (leaf SRR1509472, floral bud SRR1509478).  

MNase-Seq data: (leaf SRR1536110, floral bud SRR1536143)

##### b. Check the raw data and trimmed data quality from Chip-seq by fastqc and seqkit  
```
FastQC/fastqc -o dup/fastqc/epi_data/ -f fastq dup/epi_data/raw_data/*.fastq 
seqkit stats dup/epi_data/raw_data/*.fastq > dup/report/epi_hist_raw.txt
FastQC/fastqc -o dup/fastqc/epi_data/trim -f fastq dup/epi_data/trim/*.fastq 

```
##### c. Build the bowtie2 index for Athaliana_447_TAIR10.fa.
```
ml bowtie2/2.3.4.3
bowtie2-build -f dup/ref/rna_seq/Athaliana_447_TAIR10.fa dup/index/bowtie/Atha
```
##### d. Use fastp to remove the adaptor

```
ml fastp/0.19.4
parallel fastp -i dup/epi_data/raw_data/{}_1.fastq -I dup/epi_data/raw_data/{}_2.fastq -o dup/epi_data/trim/{}_1.fastq -O dup/epi_data/trim/{}_2.fastq ::: buds_H3K27ac buds_H3K27me3 buds_H3K4me1 leaf_H3K27ac leaf_H3K27me3 leaf_H3K4me1 
```
##### e. use bowtie to align each sequencing data to the *A.thaliana* genome, then sort the sam file and convert it to bam file.  
```
ml bowtie2/2.3.4.3
PATH=$PATH:~/samtools-1.9/
parallel bowtie2 -x dup/index/bowtie/Atha -1 dup/epi_data/trim/{}_1.fastq -2 dup/epi_data/trim/{}_2.fastq -S dup/chip_seq/{}.sam 2'>'dup/chip_seq/{}.align.log ::: buds_H3K27me3 buds_H3K4me1 leaf_H3K27ac leaf_H3K27me3 leaf_H3K4me1 buds_H3K27ac
parallel samtools-1.9/samtools sort dup/chip_seq/{}.sam -@ 8 -O bam -o dup/chip_seq/{}_sorted.bam ::: buds_H3K27me3 buds_H3K4me1 leaf_H3K27ac leaf_H3K27me3 leaf_H3K4me1 buds_H3K27ac
```
##### f. use sambamba to file the duplicate and unmatched reads
```
parallel sambamba view -h -t 8 -f bam -F "[XS] == null and not unmapped  and not duplicate" {}_sorted.bam '>' {}_sorted_filter.bam ::: buds_H3K27me3 buds_H3K4me1 leaf_H3K27ac leaf_H3K27me3 leaf_H3K4me1 buds_H3K27ac
```
##### g. use Deeptools bamCoverage to normalize and convert it to bw file for visulization
```
bamCoverage -b buds_H3K27ac_sorted_filter.bam -o visual_filter_bpm/buds_H3K27ac_sorted_filter.bw --binSize 20 --smoothLength 60 --extendReads 150 --centerRead --normalizeUsing BPM --effectiveGenomeSize 120000000 2>visual_filter_bpm/buds_H3K27ac_sorted_filter.log
```
##### h. same analysis for the MNase data
```
FastQC/fastqc -o dup/fastqc/mnase/ -f fastq dup/mnase/raw_data/*.fastq
seqkit stats dup/mnase/raw_data/*.fastq > dup/report/mnase_raw.txt
parallel bowtie2 --threads 8 -x dup/index/bowtie/Atha -1 dup/mnase/raw_data/{}_1.fastq -2 dup/mnase/raw_data/{}_2.fastq -S dup/mnase/bam/{}.sam 2'>'dup/mnase/bam/{}_align.log ::: bud_mnase leaf_mnase 
parallel samtools-1.9/samtools sort dup/mnase/bam/{}.sam -@ 8 -O bam -o dup/mnase/bam/{}.bam ::: bud_mnase leaf_mnase
parallel sambamba view -h -t 8 -f bam -F "[XS] == null and not unmapped  and not duplicate" dup/mnase/bam/{}.bam '>' dup/mnase/bam/{}_sorted_filter.bam ::: bud_mnase leaf_mnase
bamCoverage -b dup/mnase/bam/bud_mnase_sorted_filter.bam -o dup/mnase/visual_bam/bud_mnase_sorted_filter.bw --binSize 20 --smoothLength 60 --extendReads 150 --centerRead --normalizeUsing BPM --effectiveGenomeSize 120000000 2>dup/mnase/visual_bam/leaf_mnase.lognnn
bamCoverage -b dup/mnase/bam/leaf_mnase_sorted_filter.bam -o dup/mnase/visual_bam/leaf_mnase_sorted_filter.bw --binSize 20 --smoothLength 60 --extendReads 150 --centerRead --normalizeUsing BPM --effectiveGenomeSize 120000000 2>dup/mnase/visual_bam/leaf_mnase.lognnn
```
##### i. visualize the read coverage in IGV
Load the genome file with its index file, all the .bw file and gene.gff3 file into the IGV (RNA-seq, Chip-seq, MNase-seq), find the studied gene coordinates (*AT5G12950*-*AT5G12960* region: Chr5:4,089,817-4,102,332), keep the image for two types of read coverage (withou log transform and with)
##### h. EXtract the promoter sequence (1,500bp upstream) of 4 Brassicales species, use [Plant CARE](http://bioinformatics.psb.ugent.be/webtools/plantcare/html/) to predict their *cis*-activate elements (CAREs), and further use GSDS to draw the figure.
```
perl gene_family_demo/script/get_gene_weizhi.pl -in1 dup/structure/atha_lst.txt -in2 dup/structure/Athaliana_447_Araport11.gene.gff3 -out dup/promoter/atha_location
perl gene_family_demo/script/get_promoter.pl dup/data/Athaliana_447_TAIR10.fa dup/promoter/atha_location dup/promoter/atha_promoter.fa
```
Then find the unique cares and gene-specific cares. 
