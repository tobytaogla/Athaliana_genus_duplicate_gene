##### a.blastp to search homologs
###### 1. Based on the information from [OrthoDB](https://www.orthodb.org/?query=Q84W43) and [Phytozome gene family 123840092](https://phytozome-next.jgi.doe.gov/report/family/5264/123840092), get the genomic data from Phytozome, do the following analysis. 

###### b. Use diamond to make the blastp referene 
8/18/21, 10/27/21, jobtd2, protein reference index in dup/ref/protein_db/

```
ml load diamond/2.0.6
export protein_db=$HOME/dup/ref/protein_db
cd $protein_db
parallel diamond makedb -p 8 --in {}.protein.fa -d {} ::: Alyrata_384_v2.1 Athaliana_167_TAIR10 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Cpapaya_113_ASGPBv0.4 Crubella_474_v1.1 Esalsugineum_173_v1.0 Osativa_323_v7.0 RSA_r1.0 Bstricta_278_v1.2 Cgrandiflora_266_v1.1 Graimondii_221_v2.1 Tcacao_233_v1.1
```

###### c. Use diamond to blast the SSN1&SSN2 protein fasta file seperatelly to the reference, for finding the orthologs, identity >50%, query-coverage > 60%, evalue < 1e-3
8/18/21, jobtd2, blast result outfmt6 in dup/ref/blastp_orth

```
ml load diamond/2.0.6
export protein_db=$HOME/dup/ref/protein_db
export duplicate=$HOME/dup/ref/duplicate
export orth=$HOME/dup/ref/blastp_orth
parallel diamond blastp -p 8 -q $duplicate/SSN1_prot.fasta -d $protein_db/{}.dmnd --outfmt 6 --query-cover 60 --id 50 -e 1e-3 -o $orth/SSN1_blastp_{}.outfmt6 ::: Alyrata_384_v2.1 Athaliana_167_TAIR10 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Cpapaya_113_ASGPBv0.4 Crubella_474_v1.1 Esalsugineum_173_v1.0 Osativa_323_v7.0 RSA_r1.0
parallel diamond blastp -p 8 -q $duplicate/SSN2_prot.fasta -d $protein_db/{}.dmnd --outfmt 6 --query-cover 60 --id 50 -e 1e-3 -o $orth/SSN2_blastp_{}.outfmt6 ::: Alyrata_384_v2.1 Athaliana_167_TAIR10 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Cpapaya_113_ASGPBv0.4 Crubella_474_v1.1 Esalsugineum_173_v1.0 Osativa_323_v7.0 RSA_r1.0

```
###### d. Accoding to [OrthoDB](https://www.orthodb.org/?query=Q84W43) of SSN1 restult at the Embryophyta level, *Oryza sativa* Japonica Group has three homologous. In the step b search, one gene (*LOC_Os06g41030.1*) is missed. Its identity is 48.2%, lower than the parameter settings. Try to use diamond to blast the SSN1&SSN2 protein fasta file seperatelly to the reference with only one parameter evalue < 1e-3.
8/21/21, 10/27/21, jobtd2, blast result outfmt6 in dup/ref/blastp\_orth\_1

```
ml load diamond/2.0.6
export orth_1=$HOME/dup/ref/blastp_orth_1
parallel diamond blastp -p 8 -q $duplicate/SSN1_prot.fasta -d $protein_db/{}.dmnd --outfmt 6 -e 1e-3 -o $orth_1/SSN1_blastp_{}_1.outfmt6 ::: Alyrata_384_v2.1 Athaliana_167_TAIR10 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Cpapaya_113_ASGPBv0.4 Crubella_474_v1.1 Esalsugineum_173_v1.0 Osativa_323_v7.0 RSA_r1.0 Bstricta_278_v1.2 Cgrandiflora_266_v1.1 Graimondii_221_v2.1 Tcacao_233_v1.1
parallel diamond blastp -p 8 -q $duplicate/SSN2_prot.fasta -d $protein_db/{}.dmnd --outfmt 6 -e 1e-3 -o $orth_1/SSN2_blastp_{}_1.outfmt6 ::: Alyrata_384_v2.1 Athaliana_167_TAIR10 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Cpapaya_113_ASGPBv0.4 Crubella_474_v1.1 Esalsugineum_173_v1.0 Osativa_323_v7.0 RSA_r1.0 Bstricta_278_v1.2 Cgrandiflora_266_v1.1 Graimondii_221_v2.1 Tcacao_233_v1.1
```

1/2/22 jobtd1, download *A. halleri* protein reference, use the same analysis as above, but not find any homologs 
1/3/22, jobtd1, same command for Tcacao_523

```
parallel diamond blastp -p 8 -q $duplicate/{}_prot.fasta -d $protein_db/Tcacao_523_v2.1.dmnd --outfmt 6 -e 1e-3 -o $orth/{}_blastp_Tcacao_523_v2.1.outfmt6 ::: SSN1 SSN2
```

######e. Use perl script extractSeq.pl kept in $HOME/scripts/ to extract the corresponding homologous sequence, based on the diamond blast result. Combine them as one fasta file for mafft alignment
8/24/21, jobtd2.   
dup/ref/duplicate/raw\_15\_protein.fa (containing all the homologous protein sequences from selected 8 Embryophyta species);    
dup/ref/duplicate/raw\_8\_protein.fa (containing all the homologous protein sequences from selected 4 Brassicales species);    
10/27/21, jobtd2.   
dup/ref/duplicate/synteny\_protein.fa (containing all the protein sequence listed as the gene family in [Phytozome gene family 123840092](https://phytozome-next.jgi.doe.gov/report/family/5264/123840092) and *Oryza sativa* Japonica.

```
1. perl scripts/extractSeq.pl dup/ref/multiple_align_list/Alyr_lst.txt dup/ref/protein_db/Alyrata_384_v2.1.protein.fa >> dup/ref/duplicate/raw_15_protein.fa
```
######f. try use three differnt settings of mafft to build the alignment
8/24/21, jobtd2, finally use the --localpair --maxiterate 1000.   
15\_protein for raw\_15\_protein.fa; 8\_protein for raw\_8\_protein.fa; kept in folder dup/ref/duplicate/

```
ml load mafft
mafft --localpair --maxiterate 1000 dup/ref/duplicate/alignment_protein_1.fasta > dup/ref/duplicate/alignment.fa
#mafft --dash --localpair --maxiterate 1000 dup/ref/duplicate/alignment_protein_1.fasta > dup/ref/duplicate/alignment_dash
#mafft --auto dup/ref/duplicate/alignment_protein_1.fasta > dup/ref/duplicate/alignment_auto.fa
``` 
10/27/21, jobtd2, for dup/ref/duplicate/synteny/synteny\_protein.fa file

```
ml load mafft
cd dup/ref/duplicate
mafft --localpair --maxiterate 1000 synteny_protein.fa > synteny/synteny_alignment.fa
```

######g. Use iqtree to construct the Phylogenetic tree
8/25/21, jobtd2, dup/ref/duplicate/8\-align for raw\_8\_protein.fa, dup/ref/duplicate/15\-align for raw\_15\_protein.fa. alignment.fa.treefile is used for tree figure draw in iTOL <https://itol.embl.de/>.

```
export PATH=$PATH:$HOME/iqtree-1.6.12-Linux/bin/
cd dup/ref/duplicate
iqtree -s alignment.fa -nt 8 -m MFP -mtree -b 1000 
iqtree -s 8_aln.fa -nt 8 -m MFP -mtree -b 1000
```
10/27/21, jobtd2, dup/ref/duplicate/synteny for synteny_protein.fa tree fig draw

######h. Reciprocal blastp A.tha using the finding orthologs, to verify the results from above blastp results using SSN1 and SSN2
10/29/21, jobtd2, kept in /dup/ref/blastp\_orth\_2

```
other=$HOME/dup/ref/other_protein_seq
protein_db=$HOME/dup/ref/protein_db
orth_2=$HOME/dup/ref/blastp_orth_2
ml load diamond/2.0.6
parallel diamond blastp -p 8 -q $other/{}_protein.fa -d $protein_db/Athaliana_167_TAIR10.dmnd --outfmt 6 -e 1e-3 -o $orth_2/{}_blastp_Ath.outfmt6 ::: AL6G23050 AL6G23540 AL6G23554 Brara.B00465 Brara.C00538 Carub.0006s1113
```
1/25/22, jobtd2, kept in /dup/ref/blastp\_orth\_2

```
parallel diamond blastp -p 8 -q $other/{}_protein.fa -d $protein_db/Athaliana_167_TAIR10.dmnd --outfmt 6 -e 1e-3 -o $orth_2/{}_blastp_Ath.outfmt6 ::: Bostr.2902s0031 Cagra.1535s0011 Gorai.001G140300 Thecc.09G226900 Thhalv10012662m

