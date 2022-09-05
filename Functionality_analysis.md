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
Same five seperate pairwise KaKs analyses as (AL: *AL6G23540* and *AL6G23554*; AL\_050: *AL6G23540* and *AL6G23050* ; AL\_CR: *AL6G23540* and *Carub.006s1113*; AT: *AT5G12950* and *AT5G12960*; AT\_CR: *AT5G12950* and *Carub.006s1113*);  
Use the protein sequence to make the alignment by Mafft;  
Use pal2nal to make the cds alignment, based on the protein sequence alignment and cds information, output clustal format (-output clustal), remove stop codons before alignment;   
Using AXTConvertor from KaKs_Caculator2.0 to convert alignment file to axt file;  
Split the axt file to the slide_window file, use the 210 as the window length, 6 as the step length

``` 
cd KaKs_Calculator2.0/bin/Linux/ 
java split ~/Downloads/AT_cds/AT_cds_gap_aln.axt 210 6
```
Use KaKs_Caculator2.0 to compute the kaks value for each slide-window pairwise comparison

```
KaKs_Calculator2.0/bin/Linux/KaKs_Calculator -i dup/kaks/caculator_slide/gap/AT_cds_gap_alnsplit_210_6.axt -o dup/kaks/caculator_slide/gap/AT_cds_gap_alnsplit_210_6.kaks
```
Extract the KaKs values and corresponding region, draw the figure with R.
