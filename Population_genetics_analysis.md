##### a. Population genetics analysis for testing selection, including tajama.d, Fu\_Li\_D, Fu\_Li\_F, Fay\_Wu.  
Download the Pseudogenomes of the two genes (*AT5G12950* and *AT5G12960*) from [*Arabidopsis* 1001 Genomes Project](https://1001genomes.org/tools.html), then manually check the sequences of all the 1,135 accssions, make them with the same length;    

Use the python script to modify the raw date, (scripts in ~/Desktop/gene_family/script/);   
###### 1. modify the name of each sequence among 1,135 accessions by 1135_fasta_name_motify.py

```
python ../script/1135_fasta_name_motify.py pseudogenomes_AT5G12950.fa > pseudogenomes_AT5G12950_modify.fa
python ../script/1135_fasta_name_motify.py pseudogenomes_AT5G12960.fa >pseudogenomes_AT5G12960_modify.fa
```
1135_fasta_name_motify.py
```
import sys
x = sys.argv[1]
file = open(x, "r")
lines = file.readlines()
for line in lines:
    if ">" in line:
       strings = line.split("|")[3]
       strings = ">" + strings
       print(strings)
    else:
       print(line.rstrip())
```
###### 2. As 14 accessions in *AT5G12960* have the high impact mutation, filter these accessions
```
python ../script/1135_fasta_filter.py pseudogenomes_AT5G12960_modify.fa AT5G12960_stop_accession.txt > pseudogenomes_AT5G12960_modify_2.fa
``` 
1135_fasta_filter.py
```
from Bio import SeqIO
import sys

ffile = SeqIO.parse(sys.argv[1], "fasta")
header_set = set(line.strip() for line in open(sys.argv[2]))

for seq_record in ffile:
    try:
        header_set.remove(seq_record.name)
    except KeyError:
        print(seq_record.format("fasta"))
        continue
if len(header_set) != 0:
    print(len(header_set),'of the headers from list were not identified in the input fasta file.', file=sys.stderr)
```

###### 3. remove the "/n" from the above generated filter file, make each sequence in one line

```
awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' pseudogenomes_AT5G12960_modify_2.fa > pseudogenomes_AT5G12960_modify_3.fa
```

###### 4. manually align the sequence, then extract the CDS sequence
1135_fasta_AT5G12960_cds_extract.py
1135_fasta_AT5G12950_cds_extract.py 

```
python ../script/1135_fasta_AT5G12960_cds_extract.py pseudogenomes_AT5G12960_modify_3.fa > pseudogenomes_AT5G12960_cds.fa
python ../script/1135_fasta_AT5G12950_cds_extract.py pseudogenomes_AT5G12950_modify.fa > pseudogenomes_AT5G12950_cds.fa
```

###### 5. add the some reference cds sequence for certain analyses, then align them by Mafft.
 
Only population sequences for tajama.d, Fu\_Li\_D and Fu\_Li\_F: pseudogenomes_AT5G12950_cds.fa & pseudogenomes_AT5G12960_cds.fa   

Population sequences with three reference sequences (Crubella, Esalsugineum, Bstricta) for MK test: AT5G12950_out_cds_aln.fa & AT5G12960_out_cds_aln.fa   

Population sequences with the orthologs (AT5G12950 with AL6G23540, and AT5G12960 with AL6G23554) for Fay\_Wu test:
AT5G12950_cds_out_aly540_aln.fa & AT5G12960_cds_out_aly554_aln.fa

###### 7. Use R script Popgenome.R to caculate the statitics, then use tajima_draw.R to visualize the result.  

```
library(PopGenome)
#/Users/feng/Desktop/gene_family/popgenome/analysis
AT5G12950_cds.class <- readData("AT5G12950_cds",include.unknown = TRUE)
AT5G12960_cds.class <- readData("AT5G12960_cds",include.unknown = TRUE)
AT5G12950_540_cds.class <-readData("AT5G12950_cds_540/",include.unknown = TRUE)
AT5G12960_554_cds.class <-readData("AT5G12960_cds_554//",include.unknown = TRUE)
AT5G12950_mk.class <- readData("AT5G12950_cds_mk/",include.unknown = TRUE)
AT5G12960_mk.class <- readData("AT5G12960_cds_mk/",include.unknown = TRUE)

#tajama.d, Fu_Li_D and Fu_Li_F
get.sum.data(AT5G12950_cds.class)
AT5G12950_cds.class <- neutrality.stats(AT5G12950_cds.class)
get.neutrality(AT5G12950_cds.class)[[1]]
AT5G12950_cds.class@region.data@synonymous
sum(AT5G12950_cds.class@region.data@synonymous[[1]]==1) #the number of synonymous
sum(AT5G12950_cds.class@region.data@synonymous[[1]]==0) #the number of non-synonymous 

get.sum.data(AT5G12960_cds.class)
AT5G12960_cds.class <- neutrality.stats(AT5G12960_cds.class)
get.neutrality(AT5G12960_cds.class)[[1]]
AT5G12960_cds.class@region.data@synonymous
sum(AT5G12960_cds.class@region.data@synonymous[[1]]==1)
sum(AT5G12960_cds.class@region.data@synonymous[[1]]==0)


#MK test
AT5G12950_mk.class <- set.outgroup(AT5G12950_mk.class, c("Crubella", "Bstricta","Esalguineum"))
AT5G12950_mk.class <- MKT(AT5G12950_mk.class,list(1:1135,1136:1138))
get.MKT(AT5G12950_mk.class)
AT5G12950.mk.fish <- data.frame("poly" = c(47, 33), "div" = c(32, 24), row.names = c("non-syn", "syn"), stringsAsFactors = FALSE )
fisher.test(AT5G12950.mk.fish)

AT5G12960_mk.class <- set.outgroup(AT5G12960_mk.class, c("Crubella", "Bstricta","Esalguineum"))
AT5G12960_mk.class <- MKT(AT5G12960_mk.class,list(1:1121,1122:1124))
get.MKT(AT5G12960_mk.class)
AT5G12960.mk.fish <- data.frame("poly" = c(74, 16), "div" = c(65, 14), row.names = c("non-syn", "syn"), stringsAsFactors = FALSE )
fisher.test(AT5G12960.mk.fish)

#Fay-wu for cds with 540 or 554
AT5G12950_540_cds.class <- neutrality.stats(AT5G12950_540_cds.class,detail=TRUE)
get.neutrality(AT5G12950_540_cds.class)[[1]]
AT5G12950_540_cds.class <- set.outgroup(AT5G12950_540_cds.class, c("AL6G23540"))
AT5G12950_540_cds.class@region.data@outgroup
AT5G12950_540_cds.class <- neutrality.stats(AT5G12950_540_cds.class,detail=TRUE)
get.neutrality(AT5G12950_540_cds.class)[[1]]

AT5G12960_554_cds.class <- neutrality.stats(AT5G12960_554_cds.class,detail=TRUE)
get.neutrality(AT5G12960_554_cds.class)[[1]]
AT5G12960_554_cds.class <- set.outgroup(AT5G12960_554_cds.class, c("AL6G23554"))
AT5G12960_554_cds.class@region.data@outgroup
AT5G12960_554_cds.class <- neutrality.stats(AT5G12960_554_cds.class,detail=TRUE)
get.neutrality(AT5G12960_554_cds.class)[[1]]
```
