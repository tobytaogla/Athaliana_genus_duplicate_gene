##### a. download [Glyco_hydro_127 (PF07944)](https://pfam.xfam.org/family/Glyco_hydro_127#tabview=tab6) for analysis.  

##### b. search the Glyco_hydro_127 domain in the selected 4 Brassicales species. 
```
parallel hmmer-3.3/bin/hmmsearch --domtblout dup/domain/glyco_{}_hmm_out.txt --cut_tc dup/data/Glyco_hydro_127.hmm dup/ref/protein_db/{}.protein.fa ::: Alyrata_384_v2.1 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Crubella_474_v1.1 
```
###### c. Build the species-specific model in tested species from step a, and repeat the hmmer search using this model

1.extract the domain sequence from the step a obtained file after hmmer search
11/25/21, jobtd1, kept in dup/domain
Same command for the other species

```
perl scripts/domain_xulie.pl dup/domain/glyco_Alyrata_384_v2.1_hmm_out.txt dup/ref/protein_db/Alyrata_384_v2.1.protein.fa dup/domain/glyco_Alyr_domain.fa 1
```

2.Combine these domain sequence in one file (dup/domain/all_domain/all_domain.aln), aligned them and construct the phylogenetic tree; use the aligned file to build species-specific model, and repeat the hmmer search, to further identify the studied genes are the only one containing the Glyco_hydro_127 domain.  
11/26/21. jobtd1, kept in dup/domain

```
ml load mafft
export PATH=$PATH:$HOME/iqtree-1.6.12-Linux/bin/
cd dup/domain
mafft --localpair --maxiterate 1000 dup/domain/all_domain.fa > dup/domain/all_domain.aln
iqtree -s all_domain_bif.aln -nt 8 -m MFP -mtree -b 1000
```
11/27/21, jobtd1, kept in dup/domain

```
hmmer-3.3/bin/hmmbuild dup/domain/all_domain.hmm dup/domain/all_domain/all_domain.aln 
parallel hmmer-3.3/bin/hmmsearch --domtblout dup/domain/domian_{}_hmm_out.txt dup/domain/all_domain.hmm dup/ref/protein_db/{}.protein.fa ::: Alyrata_384_v2.1 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Crubella_474_v1.1 
```
###### c. extract the domain sequence from the bacteria *Bidobacterium longicum*, make the alignment among all the extracted domain seuqneces, use iqtree to construct the Phylogenetic tree based on the domain alignment
11/26/21, jobtd1, kept in dup/domain 

```
hmmer-3.3/bin/hmmsearch --domtblout dup/domain/glyco_hyba1_hmm_out.txt --cut_tc dup/data/Glyco_hydro_127.hmm dup/domain/HYBA1_BIFL2_prot.fa
perl scripts/domain_xulie.pl dup/domain/glyco_hyba1_hmm_out.txt dup/domain/HYBA1_BIFL2_prot.fa dup/domain/HYBA1_domain.fa 1
ml load mafft
export PATH=$PATH:$HOME/iqtree-1.6.12-Linux/bin/
cd dup/domain
mafft --localpair --maxiterate 1000 dup/domain/all_domain_bif.fa > dup/domain/all_domain_bif.aln
iqtree -s all_domain_bif.aln -nt 8 -m MFP -mtree -b 1000
```

###### d. verify the function domain on three website 
Using the extracted protein sequences from selected 4 Brassicales species, check the results from SMART, NCBI CDD, Pfam. Further confirmed that the Glyco_hydro_127.hmm is the only functional motif contained in all these homologs. 
SMART (http://smart.embl.de/)
CDD (https://www.ncbi.nlm.nih.gov/cdd/) 
Pfam (https://pfam.xfam.org/) 
