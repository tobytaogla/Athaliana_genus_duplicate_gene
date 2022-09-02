###### h. visualize the read coverage in IGV

A. index the reference genome for IGV and bedGraphToBigWig
2/2/22, jobtd1, dup/index/sam_index/

```
export PATH=$PATH:$HOME/samtools-1.9/
parallel samtools-1.9/samtools faidx dup/index/sam_index/{}/{}.fa ::: Alyr Atha Brap Crub
```

B. Turn each STAR aligned bam file into bedGraph coverage, file having the .bg extension

2/2/22, jobtd1, dup/star/visual/ 
same for *A. lyrata*, *A. thaliana*, *C. rubella* and *B. rapa*, the generated file in the abvoe procedures

```
ml bedtools/2.25.0 
parallel bedtools genomecov -ibam dup/star/Alyr/AL{}Aligned.sortedByCoord.out.bam -split -bg ">" dup/star/visual/Alyr/AL{}.bg ::: 1_leaf 2_flower 3_stem 4_root PA PB PC SA SB SC
```
C. Sort all the .bg files
2/2/22, local mac, /Desktop/gene_family/rna_seq_visul

```
parallel sort -k1,1 -k2,2n Atha/AT{}.bg ">" Atha/AT{}_1.bg ::: 1_leaf 2_flower 3_stem 4_root PolA PolB PolC StigA StigB StigC
```
D. Convert each bedGraph coverage into bigWig coverage. The files will have the .bw extension.
2/2/22, local mac, /Desktop/gene_family/rna_seq_visul

```
parallel bedGraphToBigWig Alyr/AL{}_1.bg index/Alyr.fa.fai Alyr/AL{}.bw ::: 1_leaf 2_flower 3_stem 4_root PA PB PC SA SB SC
```
E. Load the genome file with its index file, all the .bw file and gene.gff3 file into the IGV, find the studied gene coordinates, keep the image for two types of read coverage (withou log transform and with)
2/3/22, local mac, /Desktop/gene_family/rna_seq_visul/image
