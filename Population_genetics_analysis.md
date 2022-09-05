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

###### 5. add the some reference cds sequence for certain analyses, then align them
 
Only population sequences for tajama.d, Fu\_Li\_D and Fu\_Li\_F: pseudogenomes_AT5G12950_cds.fa & pseudogenomes_AT5G12960_cds.fa   

Population sequences with three reference sequences (Crubella, Esalsugineum, Bstricta) for MK test: AT5G12950_out_cds_aln.fa & AT5G12960_out_cds_aln.fa   

Population sequences with the orthologs (AT5G12950 with AL6G23540, and AT5G12960 with AL6G23554) for Fay\_Wu test:
AT5G12950_cds_out_aly540_aln.fa & AT5G12960_cds_out_aly554_aln.fa

These above files first use Mafft to align them, same as the analysis A and B, jobtd1 and jobtd2, 5/8/22. 


######d. Use R script Popgenome.R to caculate the statitics, then use tajima_draw.R to visualize the result.  


