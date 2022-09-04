##### a. Based on the information from [OrthoDB](https://www.orthodb.org/?query=Q84W43) and [Phytozome gene family 123840092](https://phytozome-next.jgi.doe.gov/report/family/5264/123840092), get the genomic data from [Data 3](../blob/main/Data.md), do the following analysis. 

##### b. Use diamond to make the blastp referene   
```
ml load diamond/2.0.6
export protein_db=$HOME/dup/ref/protein_db
cd $protein_db
parallel diamond makedb -p 8 --in {}.protein.fa -d {} ::: Alyrata_384_v2.1 Athaliana_167_TAIR10 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Cpapaya_113_ASGPBv0.4 Crubella_474_v1.1 Esalsugineum_173_v1.0 Osativa_323_v7.0 RSA_r1.0 Bstricta_278_v1.2 Cgrandiflora_266_v1.1 Graimondii_221_v2.1 Tcacao_233_v1.1
```
##### c.  Use diamond to blast the *AT5G12960*(SSN1)& *AT5G12950*(SSN2) protein fasta file seperatelly to the reference, for finding the orthologs. Try to use diamond to blast the SSN1&SSN2 protein fasta file seperatelly to the reference with only one parameter evalue < 1e-3  
```
ml load diamond/2.0.6
export protein_db=$HOME/dup/ref/protein_db
export duplicate=$HOME/dup/ref/duplicate
export orth_1=$HOME/dup/ref/blastp_orth_1
parallel diamond blastp -p 8 -q $duplicate/SSN1_prot.fasta -d $protein_db/{}.dmnd --outfmt 6 -e 1e-3 -o $orth_1/SSN1_blastp_{}_1.outfmt6 ::: Alyrata_384_v2.1 Athaliana_167_TAIR10 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Cpapaya_113_ASGPBv0.4 Crubella_474_v1.1 Esalsugineum_173_v1.0 Osativa_323_v7.0 RSA_r1.0 Bstricta_278_v1.2 Cgrandiflora_266_v1.1 Graimondii_221_v2.1 Tcacao_233_v1.1
parallel diamond blastp -p 8 -q $duplicate/SSN2_prot.fasta -d $protein_db/{}.dmnd --outfmt 6 -e 1e-3 -o $orth_1/SSN2_blastp_{}_1.outfmt6 ::: Alyrata_384_v2.1 Athaliana_167_TAIR10 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Cpapaya_113_ASGPBv0.4 Crubella_474_v1.1 Esalsugineum_173_v1.0 Osativa_323_v7.0 RSA_r1.0 Bstricta_278_v1.2 Cgrandiflora_266_v1.1 Graimondii_221_v2.1 Tcacao_233_v1.1
```
##### d. Extract the identified homologous protein sequence, reciprocal blastp *A.thaliana* using the finding orthologs, to verify the results from above blastp results using SSN1 and SSN2  
```
other=$HOME/dup/ref/other_protein_seq
protein_db=$HOME/dup/ref/protein_db
orth_2=$HOME/dup/ref/blastp_orth_2
ml load diamond/2.0.6
parallel diamond blastp -p 8 -q $other/{}_protein.fa -d $protein_db/Athaliana_167_TAIR10.dmnd --outfmt 6 -e 1e-3 -o $orth_2/{}_blastp_Ath.outfmt6 ::: AL6G23050 AL6G23540 AL6G23554 Brara.B00465 Brara.C00538 Carub.0006s1113 Bostr.2902s0031 Cagra.1535s0011 Gorai.001G140300 Thecc.09G226900 Thhalv10012662m
```
##### e. Use perl script extractSeq.pl to extract the corresponding homologous sequence, based on the diamond blast result, combine them as one fasta file for mafft alignment  
```
perl scripts/extractSeq.pl dup/ref/multiple_align_list/Alyr_lst.txt dup/ref/protein_db/Alyrata_384_v2.1.protein.fa >> dup/ref/duplicate/alignment_protein_1.fasta
```
##### f. use mafft to build the protein sequence alignment  
```
ml load mafft
mafft --localpair --maxiterate 1000 dup/ref/duplicate/alignment_protein_1.fasta > dup/ref/duplicate/alignment.fa
``` 
##### g. Use iqtree to construct the Phylogenetic tree  
```
export PATH=$PATH:$HOME/iqtree-1.6.12-Linux/bin/
cd dup/ref/duplicate
iqtree -s alignment.fa -nt 8 -m MFP -mtree -b 1000 
```
##### i. Use the output file alignment.fa.treefile in [iTOL](https://itol.embl.de/.)  
##### j. Extract exon and intron information from the gff file of each identified homologs, use [GSDS](http://gsds.gao-lab.org/) to show the exon/intron structure.
```
perl gene_family_demo/script/get_gene_exon_from_gff.pl -in1 dup/structure/atha_lst.txt -in2 dup/structure/Athaliana_447_Araport11.gff3 -out dup/strucure/atha_gene_exon_info.gff
```
