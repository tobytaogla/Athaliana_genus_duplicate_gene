##### a. Sequence divergent test using Codeml 
Use the protein sequence (*A. lyrate*: *AL6G23540*, *AL6G23554* and *AL6G23050*; *A. thaliana*: *AT5G12950* and *AT5G12960*; *B. stricta*: *Bostr.2902s0031*; *C. grandiflora*: *Cagra.1535s0011*; *C. rubella*: *Carub.006s1113*; *E. salsugineum*: *Thhalv10012662m*; *G. raimondii: Gorai.001G140300*;) to make the alignment by Mafft   
1. Atha_alignment.fa: (AT5G12960, AT5G12950, Carub.0006s1113, Thhalv10012662m, Gorai.001G140300, Thecc1EG040567t1, Bostr.2902s0031)  
2. Alyr_050_alignment.fa: (AL6G23050, AL6G23540, Carub.0006s1113, Thhalv10012662m, Gorai.001G140300, Thecc1EG040567t1, Bostr.2902s0031)  
3. Alyr_554_alignment.fa: (AL6G23554, AL6G23540, Carub.0006s1113, Thhalv10012662m, Gorai.001G140300, Thecc1EG040567t1, Bostr.2902s0031)  
```
ml load mafft
cd dup/kaks/paml/Atha_gap_test/
mafft --localpair --maxiterate 1000 synteny_protein.fa > Atha_alignment.fa
```
Use pal2nal to make the cds alignment, based on the protein sequence alignment and cds information.  
Before the analysis, extract the corresponding cds fasta information: Atha_cds.fasta, Alyr_cds.fa, Alyr_cds.fa 
```
pal2nal.v14/pal2nal.pl dup/kaks/paml/Atha_gap_test/Atha_alignment.fa dup/kaks/paml/Atha_gap_test/Atha_cds.fasta -output paml > dup/kaks/paml/Atha_gap_test/Atha_cds_alignment
pal2nal.v14/pal2nal.pl dup/kaks/paml/Alyr_050/Alyr_050_alignment.fa dup/kaks/paml/Alyr_050/Alyr_050_cds.fa -output paml > dup/kaks/paml/Alyr_050/Alyr_050_cds_alignment
pal2nal.v14/pal2nal.pl dup/kaks/paml/Alyr_554/Alyr_554_alignment.fa dup/kaks/paml/Alyr_554/Alyr_554_cds.fa -output paml > dup/kaks/paml/Alyr_554/Alyr_554_cds_alignment
```
Use multiple protein sequence alignment to construct the phylogenetic tree by iqtree.

```
cd dup/kaks/paml/Alyr_050/
export PATH=$PATH:$HOME/iqtree-1.6.12-Linux/bin/
iqtree -s Alyr_050_alignment.fa -nt 8 -m MFP -mtree -b 1000
```
Use codeml analysis from PAML do two-ratio model (model = 2) and one ratio model (model = 0). 
The phylogenetic tree usd from above step, as:
1. AT (AT5G12960 #1, AT5G12950 #2, (Carub.0006s1113, ((Thhalv10012662m, (Gorai.001G140300, Thecc1EG040567t1)), Bostr.2902s0031))); 
2. Alyr_050 (AL6G23050 #1, AL6G23540 #2, (Carub.0006s1113, ((Thhalv10012662m, (Gorai.001G140300, Thecc1EG040567t1)), Bostr.2902s0031))); 
3. Alyr_554 (AL6G23554 #1, AL6G23540 #2, (Carub.0006s1113, ((Thhalv10012662m, (Gorai.001G140300, Thecc1EG040567t1)), Bostr.2902s0031))).  
**two-ratio model (model = 2, NSsites = 0)**
```
cd dup/kaks/paml/Atha_gap_test/
$HOME/paml4.9j/bin/codeml Atha_free_ratio_test.ctl
```
**one ratio model (model = 0,NSsites = 0)**
```
cd dup/kaks/paml/Atha_gap_test/
$HOME/paml4.9j/bin/codeml Atha_one_ratio.ctl
```
Use R to caculater the p-value of the LRT statistic in a Chi-square distribution with two degrees of freedom. LRT statistic can be found in the output of codeml.
This is to test which model is best. 
```
pchisq(2*(-9366.060210+9379.108334),15-13,lower.tail=FALSE) #for A. tha
pchisq(2*(-9160.620827+9174.424730),15-13,lower.tail=FALSE) #for A. lyr 554
pchisq(2*(-9124.648202+9133.751786),15-13,lower.tail=FALSE) #for A. lyr 050
```

##### b. Pairwise KaKs value computation using KaKs_Caculator2.0
First two steps are same as the analysis a  
Use the protein sequence (*A. lyrate*: *AL6G23540*, *AL6G23554* and *AL6G23050*; *A. thaliana*: *AT5G12950* and *AT5G12960*; *C. rubella*: *Carub.006s1113*;) to separately make the alignment by Mafft, same as the previous analyses.
As: (*AT5G12950* and *AT5G12960*), (*AT5G12960*, *Carub.006s1113*), (*AL6G23540*, *AL6G23554*), (*AL6G23554*, *AL6G23050*), (*AL6G23050*, *Carub.006s1113*)
Use pal2nal to make the cds alignment, based on the protein sequence alignment and cds information, output clustal format (-output clustal), remove stop codons before alignment
```
pal2nal.v14/pal2nal.pl dup/kaks/caculator/caculator_alignment.fa dup/kaks/caculator/caculator_cds.fa -output clustal > dup/kaks/caculator/gap/caculator_cds_gap_aln
```
Using AXTConvertor from KaKs_Caculator2.0 to convert alignment file to axt file; then use KaKs_Caculator2.0 to compute the kaks value pairewiely
```
KaKs_Calculator2.0/src/AXTConvertor dup/kaks/caculator/caculator_cds_aln dup/kaks/caculator/caculator_cds_aln.axt
KaKs_Calculator2.0/bin/Linux/KaKs_Calculator -i dup/kaks/caculator/caculator_cds_aln.axt -o dup/kaks/caculator/caculator_cds_aln.kaks
```
##### c. Slide-window pairwise comparison KaKs analysis using KaKs_Caculator2.0
First two steps are same as the analysis a   
Same five seperate pairwise KaKs analyses as (AL: *AL6G23540* and *AL6G23554*; AL\_050: *AL6G23540* and *AL6G23050* ; AL\_CR: *AL6G23540* and *Carub.006s1113*; AT: *AT5G12950* and *AT5G12960*; AT\_CR: *AT5G12950* and *Carub.006s1113*)
Use the protein sequence to make the alignment by Mafft
Use pal2nal to make the cds alignment, based on the protein sequence alignment and cds information, output clustal format (-output clustal), remove stop codons before alignment 
Using AXTConvertor from KaKs_Caculator2.0 to convert alignment file to axt file; 
Split the axt file to the slide_window file, use the 210 as the window length, 6 as the step length

``` 
cd KaKs_Calculator2.0/bin/Linux/ 
java split ~/Downloads/AT_cds/AT_cds_gap_aln.axt 210 6
```
######e. use KaKs_Caculator2.0 to compute the kaks value for each group

```
KaKs_Calculator2.0/bin/Linux/KaKs_Calculator -i dup/kaks/caculator_slide/gap/AT_cds_gap_alnsplit_210_6.axt -o dup/kaks/caculator_slide/gap/AT_cds_gap_alnsplit_210_6.kaks
```
######e. extract the KaKs values and corresponding region, draw the figure with R script kaks.R

##### D. Population genetics analysis for testing selection, including tajama.d, Fu\_Li\_D, Fu\_Li\_F, Fay\_Wu and MK test.  
aligned files kept in grid dup/pop/
analysis file kept in local ~/Desktop/gene_family/popgenome/
raw and clean files in local ~/Desktop/gene_family/dnaSP/
######a. Download the Pseudogenomes of the two genes from <https://1001genomes.org/tools.html>, then manually check the sequences of all the 1,135 accssions, make them with the same length  
file name pseudogenomes_AT5G12950 and pseudogenomes_AT5G12960 (~/Desktop/gene_family/dnaSP/)

######b. use the local script to modify the raw date, (scripts in ~/Desktop/gene_family/script/) 

b.1. modify the name of eahc sequence by 1135_fasta_name_motify.py

```
python ../script/1135_fasta_name_motify.py pseudogenomes_AT5G12950.fa > pseudogenomes_AT5G12950_modify.fa
python ../script/1135_fasta_name_motify.py pseudogenomes_AT5G12960.fa >pseudogenomes_AT5G12960_modify.fa
```
b.2. As 14 accessions have the high impact mutation, need to filter these accessions
1135_fasta_filter.py

```
python ../script/1135_fasta_filter.py pseudogenomes_AT5G12960_modify.fa AT5G12960_stop_accession.txt > pseudogenomes_AT5G12960_modify_2.fa
``` 

b.3. remove the "/n", make each sequence in one line

```
awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' pseudogenomes_AT5G12960_modify_2.fa > pseudogenomes_AT5G12960_modify_3.fa
```

b.4. manually align the sequence, then extract the CDS sequence
1135_fasta_AT5G12960_cds_extract.py
1135_fasta_AT5G12950_cds_extract.py 

```
python ../script/1135_fasta_AT5G12960_cds_extract.py pseudogenomes_AT5G12960_modify_3.fa > pseudogenomes_AT5G12960_cds.fa
python ../script/1135_fasta_AT5G12950_cds_extract.py pseudogenomes_AT5G12950_modify.fa > pseudogenomes_AT5G12950_cds.fa
```

######c. add the some reference cds sequence for certain analyses, then align them
 
Only population sequences for tajama.d, Fu\_Li\_D and Fu\_Li\_F: pseudogenomes_AT5G12950_cds.fa & pseudogenomes_AT5G12960_cds.fa   

Population sequences with three reference sequences (Crubella, Esalsugineum, Bstricta) for MK test: AT5G12950_out_cds_aln.fa & AT5G12960_out_cds_aln.fa   

Population sequences with the orthologs (AT5G12950 with AL6G23540, and AT5G12960 with AL6G23554) for Fay\_Wu test:
AT5G12950_cds_out_aly540_aln.fa & AT5G12960_cds_out_aly554_aln.fa

These above files first use Mafft to align them, same as the analysis A and B, jobtd1 and jobtd2, 5/8/22. 


######d. Use R script Popgenome.R to caculate the statitics, then use tajima_draw.R to visualize the result.  

