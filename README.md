# Neofunctionalization of tandem duplicate genes encoded novel putative beta-L-arabinofuranosidases in Arabidopsis thaliana analysis
Workflow and scripts used for this project bioinformatics analysis


![](https://github.com/tobytaogla/Phragmites-australis-transcriptome-optimal-assembly/blob/main/Flowchart.png)


## 1. Introduction
Tandem duplication, one of the major types of duplication, provides the raw material for the evolution of novel functions. In this study, we identified one pair of tandem duplicate genes *AT5G12950* and *AT5G12960* that originated within the last 16 million years after the split of *Arabidopsis* from the Capsella-Boechera ancestor. We systematically used bioinformatic tools to define their putative protein function as the only two beta-L-arabinofuranosidases that release L-Arabinose from the beta-L-Araf-containing molecules in Arabidopsis thaliana. Comprehensive transcriptomic and proteomic analyses using different datasets revealed distinct expression patterns between the two duplicate genes, which resulted from the regulatory neofunctionalization. Both functionality tests and population genetic analysis showed that the protein-coding regions of the two duplicate genes are under functional constraints. We further collected data from two types of phenotyping measurements to indicate that *AT5G12950* and *AT5G12960* have roles that contribute to different phenotypic effects. Overall, the A. thaliana tandem duplicate genes *AT5G12950* and *AT5G12960* represent the first putative beta-L-arabinofuranosidases coding genes characterized in *A. thaliana*. After duplication, the regulatory divergence drove the new duplicate gene to develop neofunction and contribute to a different phenotypic evolution in *Arabidopsis*. 


## Data 
Sequencing data information could be also found in supplementary Data 7.

### 1. RNA_seq raw data:
Raw illumina sequencing data 
***Arabidopsis thaliana*** used tissues: pollen * 3, stigma * 3, leaf, floral bud, root, stem.  	
***Arabidopsis lyrata*** used tissues: pollen * 3, stigma * 3, leaf, floral bud, root, stem.  
***Capsella rubella*** used tissues: pollen * 3, stigma * 3, leaf, floral bud, root, stem.  
***Brassica rapa*** used tissues: pollen * 3, stigma * 3
	
### 2. Phytozome data for RNA\_seq: 
Download from [Phytozome 13] (https://phytozome-next.jgi.doe.gov/)  
***Arabidopsis lyrata***: Alyrata\_384\_v1.fa, Alyrata_384\_v2.1.gene.gff3,  Alyrata\_384\_v2.1.transcript.fa.    
***Arabidopsis thaliana***: Athaliana\_447\__TAIR10.fa, Athaliana\_447\_Araport11.gene.gff3,  Athaliana\_447\_Araport11.transcript.fa.  
***Brassica rapa***: BrapaFPsc\_277\_v1.3.gene.gff3,  BrapaFPsc\_277\_v1.3.transcript.fa,  BrapaFPsc\_277\_v1.fa.   
***Capsella rubella***: Crubella\_474\_v1.1.gene.gff3,  Crubella\_474\_v1.1.transcript.fa,  Crubella\_474\_v1.fa.

### 3. Phytozome data to identify orthologs
Download from [Phytozome 13] (https://phytozome-next.jgi.doe.gov/)  
**Protein fasta file**: Alyrata\_384\_v2.1, Athaliana\_167\_TAIR10, Athaliana\_447\_Araport11, BrapaFPsc\_277\_v1.3, Cpapaya\_113\_ASGPBv0.4, Crubella\_474\_v1.1, Esalsugineum\_173\_v1.0, (Osativa\_323\_v7.0*), Bstricta\_278\_v1.2, Cgrandiflora\_266\_v1.1, Graimondii\_221\_v2.1, Tcacao\_523\_v2.1.  
**CDS fasta file**: Alyrata\_384\_v2.1, Athaliana\_447\_Araport11, BrapaFPsc\_277\_v1.3, Crubella\_474\_v1.1, Esalsugineum\_173\_v1.0, Bstricta\_278\_v1.2, Cgrandiflora\_266\_v1.1, Graimondii\_221\_v2.1, Tcacao\_523\_v2.1.   
**Genome fasta file**: Alyrata\_384\_v2.1, Athaliana\_167\_TAIR10, Athaliana\_447\_Araport11, BrapaFPsc\_277\_v1.3, Crubella\_474\_v1.1, Esalsugineum\_173\_v1.0, Bstricta\_278\_v1.2, Cgrandiflora\_266\_v1.1.  
**gff3 file**: Alyrata\_384\_v2.1, Athaliana\_167\_TAIR10, Athaliana\_447\_Araport11, BrapaFPsc\_277\_v1.3, Crubella\_474\_v1.1, Esalsugineum\_173\_v1.0, Bstricta\_278\_v1.2, Cgrandiflora\_266\_v1.1.   
*Osativa was not used at the end

### 4. Population seuqence data across 1,135 accessions
The sequence of the two duplicated *A. thaliana* genes (*AT5G12950* and *AT5G12960*) across 1,135 accessions was downloaded from the [*Arabidopsis* 1001 Genomes Project] (https://1001genomes.org/tools.html)  
Choose Pseudogenomes, select all the accessions, search for the desired gene, and download the data.  

### 5. Epigentics sequencing data
Raw data downloaded from the [NCBI Sequence Read Archive (SRA)] (https://www.ncbi.nlm.nih.gov/sra)  

Chip-Seq data:
H3K27ac  (leaf SRR1509474, floral bud SRR1509479)
H3K4me1  (leaf SRR1509477, floral bud SRR1509476)
H3K27me3 (leaf SRR1509472, floral bud SRR1509478) 

MNase-Seq data: (leaf SRR1536110, floral bud SRR1536143)

### 6. Other data
##### [Glyco_hydro_127 (PF07944)] (https://pfam.xfam.org/family/Glyco_hydro_127#tabview=tab6)  
Beta-L-arabinofuranosidase, GH127.  
One member of this family, from *Bidobacterium longicum*, UniProtKB:E8MGH8 has been characterised as an unusual beta-L-arabinofuranosidase enzyme, EC:3.2.1.185. It rleases l-arabinose from the l-arabinofuranose (Araf)-beta1,2-Araf disaccharide and also transglycosylates 1-alkanols with retention of the anomeric configuration. Terminal beta-l-arabinofuranosyl residues have been found in arabinogalactan proteins from a mumber of different plant species. beta-l-Arabinofuranosyl linkages with 1-4 arabinofuranosides are also found in the sugar chains of extensin and solanaceous lectins, hydroxyproline (Hyp)2-rich glycoproteins that are widely observed in plant cell wall fractions. The critical residue for HypBA1 catalytically important residues (Glu322, Glu338, Cys340, Cys417, and Cys418).

##### Protein predicted structure
Download from the [AlphaFold DB] (https://alphafold.ebi.ac.uk/) 
[AT5G12950] (https://alphafold.ebi.ac.uk/entry/Q9LXU4)
[AT5G12960] (https://alphafold.ebi.ac.uk/entry/Q84W43)

## 2. Software
##### Linux software:
Trimmomatic-0.35, star/2.6.1d, diamond/2.0.6, iqtree-1.6.12, seqkitv0.16.1, fastqc v0.11.8, subread-2.0.3, MAFFT v7.487, FastQCv0.11.8, TMHMM-2.0, ncbi-blast-2.8.1+, PAML4.9j, PAL2NAL v11, KaKs_Caculator2.0, HMMER-3.3, samtools-1.9, bedtools 2.25.0, bedGraphToBigWig v4, pyfaidx, deeptools 3.5.1, Sambamba-0.8.2, macs2 2.2.7.1, fastp/0.19.4, bowtie2 2.3.4.3, 

##### R package: 
ComplexHeatmap, PopGenome 2.7.5, FactoMineR 2.4

##### Mac software: 
Genevestigator 9.1.0, Jalview 2.11.2.2, TBtools v1.09852,  IGV_2.9.4, PyMOL 2.5.2

##### Windows:
GeneDoc

##### Bioinformatics online tools
iTOL <https://itol.embl.de/>.   
DeepTMHMM <https://dtu.biolib.com/DeepTMHMM>   
meme v5.4.1 <https://meme-suite.org/meme/tools/meme>  
signalP v5.0 <https://services.healthtech.dtu.dk/service.php?SignalP-5.0>  
CDD <https://www.ncbi.nlm.nih.gov/cdd/>  
Pfam <http://pfam.xfam.org/>    
SMART <http://smart.embl.de/>   
ExPASy <https://web.expasy.org/protparam/>   
GSDS <http://gsds.gao-lab.org/>    
Plant CARE <http://bioinformatics.psb.ugent.be/webtools/plantcare/html/>   
TSSFinder <http://sucest-fun.org/wsapp/tssfinder/>  
WoLF PSORT <https://wolfpsort.hgc.jp/>  
CELLO <http://cello.life.nctu.edu.tw/>  
AlphaFold DB <https://alphafold.ebi.ac.uk/>  
SUBA4 <https://suba.live>   
ATHENA <http://athena.proteomics.wzw.tum.de/>
