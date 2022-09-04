##### a. download [Glyco_hydro_127 (GH127, PF07944)](https://pfam.xfam.org/family/Glyco_hydro_127#tabview=tab6) for analysis.  

##### b. search the Glyco_hydro_127 domain in the selected 4 Brassicales species. 
```
parallel hmmer-3.3/bin/hmmsearch --domtblout dup/domain/glyco_{}_hmm_out.txt --cut_tc dup/data/Glyco_hydro_127.hmm dup/ref/protein_db/{}.protein.fa ::: Alyrata_384_v2.1 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Crubella_474_v1.1 
```
###### c. Build the Brassicales-specific model in tested species from step a, and repeat the hmmer search using this model

extract the domain sequence from the step a obtained file after hmmer search， same command for the other tested species
```
perl scripts/domain_xulie.pl dup/domain/glyco_Alyrata_384_v2.1_hmm_out.txt dup/ref/protein_db/Alyrata_384_v2.1.protein.fa dup/domain/glyco_Alyr_domain.fa 1
```

Combine these domain sequence in one file (all_domain.aln), build species-specific model, and repeat the hmmer search, to further identify the studied genes are the only one containing the Glyco_hydro_127 domain.  
```
hmmer-3.3/bin/hmmbuild dup/domain/all_domain.hmm dup/domain/all_domain/all_domain.aln 
parallel hmmer-3.3/bin/hmmsearch --domtblout dup/domain/domian_{}_hmm_out.txt dup/domain/all_domain.hmm dup/ref/protein_db/{}.protein.fa ::: Alyrata_384_v2.1 Athaliana_447_Araport11 BrapaFPsc_277_v1.3 Crubella_474_v1.1 
```
##### d. extract the domain sequence from the bacteria *Bidobacterium longicum* using the perl script domain_xulie.pl, make the alignment among all the extracted domain seuqneces, use iqtree to construct the Phylogenetic tree based on the domain alignment
```
hmmer-3.3/bin/hmmsearch --domtblout dup/domain/glyco_hyba1_hmm_out.txt --cut_tc dup/data/Glyco_hydro_127.hmm dup/domain/HYBA1_BIFL2_prot.fa
perl scripts/domain_xulie.pl dup/domain/glyco_hyba1_hmm_out.txt dup/domain/HYBA1_BIFL2_prot.fa dup/domain/HYBA1_domain.fa 1
ml load mafft
export PATH=$PATH:$HOME/iqtree-1.6.12-Linux/bin/
cd dup/domain
mafft --localpair --maxiterate 1000 dup/domain/all_domain_bif.fa > dup/domain/all_domain_bif.aln
iqtree -s all_domain_bif.aln -nt 8 -m MFP -mtree -b 1000
```
Use this alignment to check the catalytically critical residues are conserved in all the identified Brassicaceae protein domains.
 
##### e. verify the function domain on three website 
Using the extracted protein sequences from selected 4 Brassicales species, check the results from SMART, NCBI CDD, Pfam. Further confirmed that the Glyco_hydro_127.hmm is the only functional motif contained in all these homologs. 
SMART <http://smart.embl.de/>  
CDD <https://www.ncbi.nlm.nih.gov/cdd/>   
Pfam <https://pfam.xfam.org/>  

##### f. use [signalP v6.0](https://services.healthtech.dtu.dk/service.php?SignalP-6.0) and [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM) to check whether the iddentified Brassicaceae proteins had the signal and/or transmemberane peptides. 

##### g. use [SUBA4](https://suba.live) to predict AT5G12950 and AT5G12960 subcelluar localization. 

##### h. use [MEME v5.4.1](https://meme-suite.org/meme/tools/meme) to predict the top 10 conserved motifs in the GH127 domain of each identified Brassicaceae putative beta-L-arabinofuranosidase. 

##### i. download the protein structure of AT5G12950 and AT5G12960, and analyzed them by PyMOL. 
