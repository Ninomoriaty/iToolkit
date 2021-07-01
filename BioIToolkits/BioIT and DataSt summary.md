Guideline Level

<table><tr><td  style="background: #3ADDBA !important;"><center><font style="font-size:45px; color:black !important;"><b>
    Bioinformatics Toolkits and Data Structure in Bioinformatic
    </b></font></center></td></tr></table>
# Ref

Tips: Needed to be assigned to corresponding  Section

1. (NCBI) BLAST User guide: https://www.ncbi.nlm.nih.gov/books/NBK279690/
2. Another BLAST Usage Notebook: https://www.ccg.unam.mx/~vinuesa/tlem/pdfs/Bioinformatics_explained_BLAST.pdf



# Lk

- Lk: common workflow will be stored in this folder and detailed records and manual can be found in the Wiznote
- E:\0 Science\Biology\QB Quantitative Biology-计量生物学\BioI-生信\BioIToolkit\NCBI\BLAST



# Pcp-Main

## Bioinformatic Toolkit Section

- This note book summarise all bioinformatic tools for different usage and provide proper Lk for further manual or usage.
  - Thus this note book focus on complete collection and reference guideline rather than detailed usage.
  - More information should be found through Lk.
- To simplify the and make the search for particular tool more easily, this notebook only contain: 
  - Ref: The website for download, documents and further information
  - Lk: The Wiznote note for further usage
  - F: The main functions of the tool
  - Mt-Install: The basic installation method
  - Tips: Reminder for some preference settings or prerequisites.
- The pipelines for different stages can be found in Lk notes and you should check them in corresponding Wiznote note.



## Database in BioI Section

- This note book also contain the collection of bioinformatic database for different function.
-  



## Data Structure in BioI Section

- Talking of database, the bioinformatic-specific data structures will also be introduced in this note book and you can find the 



## High-Performance Computation Section

- This part introduce some usage of the high performance machine during my learning experience.





<table><tr><td  style="background: #3ADDBA !important;"><center><font style="font-size:45px; color:black !important;"><b>
    Bioinformatics Toolkits
    </b></font></center></td></tr></table>

# 1 Bioinformatics Toolkits

# Knowledge Mapping for papers

### Citespace



# Synthetic Circuit Design

### iGEM cello



# Wet Experiment Manager

## Integrated

### Vector NTI Express Designer



## Primer/Adaptor Design

### Primer Premier



# Data Download and Data Conversion

### SRAToolkit

**Ref**

- GitHub Documentation:https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit
- 

**Obj1+2**

- Datasets from NCBI database
- 

**Mt-faster-dump**



**F**

- Standard method to download data from NCBI
- sequence data conversion
- 

### ISPEMRA

**Func**

- Alternative/Integrated download method for NCBI and EBI data.

### seqtk

**Func**

- Convert fastq to fasta
- parse and filter sequence data with hard conditions
- 



# Sequencing and Assembly

## k-mers and k-mer spectra

### Jellyfish



### kmc



### KAT





## RNA-Seq

### Rascaf

Ref: https://github.com/mourisl/Rascaf

### L-RNA_SCAFFOLDER

Ref: https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-604

### SCUBAT

Ref: https://github.com/GDKO/SCUBAT2

### RNAPATH in ERANGE

Ref: https://genome.cshlp.org/content/20/12/1740.full.html



## Sequence Alignment

### BLAST

**Ref**:

1. Web-based toolkit: https://blast.ncbi.nlm.nih.gov/Blast.cgi

**Soc**

1. NCBI

**Def**: identify  putatively homologous sequences in a database

**Obj2**:

1. DNA sequences
2. Protein sequences
3. Phylogenetic tree groups

**Ab**:

1. Data Ab:
   1. DNA
      1. nucleotide sequences
   2. Protein
      1. Amino acid sequences
      2. Distance in genetic code and components from different genes in the genome will make the protein sequences hard to match
      3. The chemical similarity in some amino acid and protein structures make the protein alignment more challenging.
      4. Existing alignments and structures may not suit the new structures of the sequencing proteins.

**Mt-Install**

```bash
sudo apt-get install ncbi-blast+
```

**Mt-BLAST process**:

1. <u>Build an index</u> of short words for all sequences in a database
2. <u>Search matching short words</u> in the query sequences  
3. Extend each match to the left or the right to <u>make the longest match</u> (High scoring Segment Pair, HSP)
4. Stop the extending when the match score does go up
5. Try to merge HSPs to get even longer matches

**Mt-Use BLAST in command line**

​    Ref: https://www.ncbi.nlm.nih.gov/books/NBK25501/

1. build the BLAST database

   1. Select proper database

      1. NCBI database list: https://www.ncbi.nlm.nih.gov/guide/all/
         1. nr : non-redundant protein database
         2. nt : non-redundant nucleotide database
         3. swiss-prot : manually annotated protein database
         4. pat : protein sequences from patents
         5. month : recently added sequences
         6. env_nr : environmental sequences
         7. est_human : only human Expressed Sequence Tag sequences
      2. Access the database remotely when the database is relatively large

   2. Get the datasets from the database

      ```bash
      # check the workplace first
      pwd
      
      # download the datasest
      wget http://nematodes.org/teaching/bioinformatics/nem.fasta
      
      # get particular database txt from the html
      wget -O testsequence.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_015119.1&strand=1&seq_start=1426&seq_stop=2962&rettype=fasta&retmode=text"
      # Explanation
      -O testsequence.fasta (the Output file)
      efetch.fcgi? (there are many NCBI utilities, described in more detail at https://www.ncbi.nlm.nih.gov/books/NBK25501/)
      db=nuccore (which database)
      id=NC_015119.1 (sequence AccessionID)
      strand=1 (top strand)
      seq_start=1426 (start from base 1426 of the sequence)
      seq_stop=2962 (stop at base 2962 of the sequence)
      rettype=fasta (return in fasta format)
      retmode=text (return in plain text format please!)
      
      ```

   3. Inspect the data and see if the data is correct

      ```bash
      cat nem.fasta | awk '{
      	# Is this line a fasta header
      	if(substr($1,1,1)==">")
      	{
      	 print "this is a header line: " $0; 
      	 hline=FNR ;
      	}
      	# Identify the first line of the sequence
      	# ...and then the first character of that line
      	if(FNR==hline+1)
      	{
      	 print "First line is:" $0;
           first_seq_character=substr($0,1,1)
           print "First character is: " first_seq_character
           print first_seq_character > "first_seq_character.txt" 
          }
      }'
      ```

      Normally, FASTA data will have the structure described in the Data Structure in BioI part in this notebook.

   4. Build the database structure that suit the usage of blast analysis

      ```bash
      # help manual
      makeblastdb --help
      
      # Common usage
      makeblastdb -in nem.fasta -dbtype prot -out nem
      # -in the input of dataset
      # -dbtype the type of database
      # -out the output NAME of the database for blast
      
      # Input
      nem.fasta: database input 
      
      # Output
      nem.phr: protein sequences headers(Names of the protein and basic information)
      nem.pin: protein index(The locations of the short words)
      nem.psq: compressed protein sequences
      ```

   5. Get the query sequences

      From a summary accession list provided by NCBI:

      ```bash
      while read acc
      do
      echo -e "Downloading $acc"
      wget -q -O $acc.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=$acc&strand=1&rettype=fasta&retmode=text"
      done < 590accessions.txt 
      ```

      Actually, `efetch` will be a better choice

      ```bash
      
      ```

   6. 

2. Search and match targeted sequences

   ```bash
   # DNA - DNA: 
   blastn
   megablast
   tblastx
   
   # DNA - Protein
   blastx -db nem -query testsequence.fasta
   blastx -db nem -query testsequence.fasta -outfmt 7 > blastoutput1.out
   # -db: database
   # -query: sequences needed to match to the database
   # -outfmt 7: the output format tsv format
   
   # Protein - DNA
   tblastn
   
   # Protein -Protein
   blastp
   quickblastp
   psi-blast
   
   ```

   **Mt-Choose proper methods for pairs of query and database**

   Tips: The first is the query and the second is the database

   1. DNA - DNA: blastn, megablast, tblastx, [Magic BLAST](https://ncbi.github.io/magicblast/)
   2. DNA - Protein: blastx
   3. Protein - DNA: tblastn
   4. Protein -Protein: blastp, quickblastp, psi-blast

3. Inspect and interpret the results

   1. **Comment** lines start with #
      1. **Fields** in the comments section are separated by a space
   2. **Data Section**
      1. each data section is a HSP
      2. separated by tab or the output format
   3. 

4. 

5. 

**Mt-Bitscore Calculation**

1. Ab:
   1. The scoring matrix
      1. [Karlin-Altschul statistics](http://www.drive5.com/usearch/manual/karlin_altschul.html) are used to convert the HSP bitscore to an expectation value (E-value) 
      2. bigger HSP bitscore are "better" and smaller E-values are "better" 
      3. not all HSP are equally comparable
   2. The length of query
   3. the word size
   4. the size of the database
   5. 

**Mt-Scoring matrix**

Ref: [Selecting the Right Similarity-Scoring Matrix](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848038/)

DNA

1. The simplest example: just count the match number and give it 1 value and the unmatched nucleotide get 0 value.
2. 

Protein

1. **p**oint **a**ccepted **m**utation (**PAM**): the replacement of a single amino acid in the primary structure of a protein with another single amino acid, which is accepted by the processes of natural selection based on global (end to end) alignments of highly similar proteins.

   ![Screenshot from 2020-11-13 19-55-40.png](https://i.loli.net/2020/11/14/PNIgrQAHKZYtXaV.png)

   1. Cat:
      1. PAM N corresponds to <u>N mutations in DNA sequence</u> per 100 amino acids. N can be > 100
      2. PAM250 is one of the more common used scoring matrices.
   2. Ab/Pcp-PAM model naming convention: 
      1. high number for distantly related proteins, low number for closely related protein 

2. **BLOSUM matrices**: BLOSUM (**BLO**cks **SU**bstitution **M**atrix) matrix is a substitution matrix used for scoring sequence alignment of <u>evolutionary divergent</u> protein sequences based on local alignments.

   ![Screenshot from 2020-11-13 19-55-46.png](https://i.loli.net/2020/11/14/DeAm4h6VLpMzgNq.png)

   1. Ref: https://www.pnas.org/content/89/22/10915.long
   2. Ab
      1. All BLOSUM matrices are based on observed alignments; they are not extrapolated from comparisons of closely related proteins like the PAM Matrices.
   3. Mt
      1. scanned the BLOCKS database for very conserved regions of protein families(UNGAPPED local alignments)
      2. counted the relative frequencies of amino acids and their substitution probabilities.
      3. 
   4. Cat
      1. BLOSUM62 is the matrix built using sequences with <u>no more than 62% similarity</u>; BLOSUM62 is the <u>default matrix used for protein BLAST</u>.
      2. BLOSUM80 is for less divergent alignments, whereas BLOSUM 50 is for more divergent alignments.

**F**:

1. Predict statistical evolution relationship but not the certain relationship
2. Pairwise alignment
3. <u>Functional identification</u> of novel sequences
4. <u>annotation</u> of genomes
5. <u>localisation</u> of particular sequences in <u>genomes</u>
6. screening of <u>environments</u> using set of sequences
7. <u>clustering</u> of sequence data to identify sequence families
8. To be continued...



### Clustalo

**Ref**: http://www.ebi.ac.uk/Tools/msa/clustalo/

**Source**: EBI

**F**:

1. You can use the `Clustalo Omega` website toolkit to align the sequences online. The result and summary will be valid for 7 days.
2. 

**Cat**:

1. Clustalo
2. Clustalw

**Mt-install**

```bash
sudo apt install clustalo
```

### bowtie2

### bwa

**F**

- Read mapper

### hisat2

### picard

### IGV

**Def**

- integrated genome viewer

**F**

- 



## Sequencing Quality Control(QC)

<font color=orange size=6>Quality Assessment</font>

### FastQC

**Ref**

- Website https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

**Lk**

- Further usage in sequencing assessment: Wiznote > 0 Science > QB-计量生物学 > BioI-生信 > BioIToolkit > Sequencing Analysis Research

**Mt-Install from zip file**

- Download the zip file from the website

- Extract the folder from the zip file

- ```bash
  # Access the fastqc program
  chmod 755 ./fastqc
  # Settle the PATH for fastqc in .bashrc
  export PATH=~/[Your dir]/FastQC:$PATH
  # Try
  fastqc
  ```

**Mt-Install from the pkg store**

- ```bash
  # From apt
  sudo apt install fastqc
  # Frome conda
  conda install
  ```

### NanoPlot

**Ref**

**Obj1**

- 

**F**

- Visualise the Quality Assessment of Nanopore genome data



------

<font color=orange size=6>Quality Improvement</font>

### Cutadapt

**Ref**

https://cutadapt.readthedocs.io/en/stable/

**F**

- 







## Genome Assembly

<font color=orange size=6>Short Read</font>

### Quiver

Ref: https://github.com/PacificBiosciences/GenomicConsensus

F: Uses shorter PB reads to correct longer ones

RSI, RSII early 

### Arrow

Ref: https://github.com/PacificBiosciences/GenomicConsensus

RSII and SEQUEL

### Nanopolish

Ref: https://github.com/jts/nanopolish

F: Nanopore only; uses original voltage data to correct errors

------

<font color=orange size=6>Long Read</font>

### PBJelly

Cat

### Juicebox

Cat+F: HiC Contact Map

### Hierarchical Genome Assembly Process (HGAP)

Ref: https://github.com/PacificBiosciences/Bioinformatics-Training/wiki/HGAP

Pcp

- Use the long reads to construct pre-assembled reads (longest) as corrected seeds
- Then use the pre-assembled reads to assemble the finished genome

Cat+F: 

- Long Read Assembly
- reduce the influence of relative short reads by using pre-assembly method.

### Canu 

Ref: https://github.com/marbl/canu

Obj1: PacBio or Nanopore

Pcp

- OLC

Pc

1. Correct
2. Trim
3. OLC Assembly

Cat+F: 

- Long Read Assembly
- 

### FALCON / FALCON-unzip

Ref: https://github.com/PacificBiosciences/FALCON

Obj1: PacBio

Cat+F: 

- Long Read Assembly
- Show the haplotype-resolved assembly graph and the SNPs and SV situation for the assembly. 

### Miniasm 

Ref: https://github.com/lh3/miniasm

Obj1: PacBio and Nanopore

Pcp

- OLC
- Use k-mers minimisers to allow fast matching and allow low Q sequence

Ab

- Able to deal with noisy long read data
- No error correction
- Fast
- 

Cat+F: Long Read Assembly

### wtdbg2/RedBean

Ref: https://github.com/ruanjue/wtdbg

Pcp

- De Bruijn Graph(DBG) approach
- Use FBG of K-bins(a K-bin is 256 bases)
- Polish after assembly

Ab

- High quality

Pc

1. Binning the reads
2. Match the overlaps of each reads' k-bin and construct the DBG
3. Clean the DBG
4. Consensus

Cat+F: Long Read Assembly



## Genome Assembly Quality Control

<font color=orange size=6>Quality Assessment</font>

### BUSCO

**Ref**

- GitLab: https://gitlab.com/ezlab/busco
- Documentation: https://busco.ezlab.org/busco_userguide.html



### CEGMA

**Tips**

- This benchmark was discontinued from 2015 and you should choose BUSCO to do the analysis

**Ref**

- 

### assembly-stats

**Ref**

- GitHub: https://github.com/rjchallis/assembly-stats
- Website and Documentation: https://assembly-stats.readme.io/docs

**Tips**

- This project is not yet completed and the documentation can be the worst manual for me.

**Func**

- Visualise the statistics of contigs in genome assembly

**Mt-Install**

- For GitHub
  - Download repo from the GitHub
  - Set up the PATH for pl folder of the repo in .bashrc

**Mt-Usage**

- First, use the Perl Script in the pl folder of the repo to generate json object

  ```bash
  # Remember the json object should have the .assembly-stats.json for HTML usage
  perl asm2stats.pl genome_assembly.fa > output.assembly-stats.json
  ```

  

- Then insert the json object into the HTML file of the repo (Tips: The json)

  ```html
  <div id="assembly_stats">
  <script>
    d3.json("Danaus_plexippus_v3.assembly-stats.json", function(error, json) {
      if (error) return console.warn(error);
      asm = new Assembly (json);
      asm.drawPlot('assembly_stats');
    })
  </script>
  ```

- Use WebStorm open the porter and visualise the HTML page with &assembly=[name_before_.assembly-stats.json ] modification

Mt-Example

http://localhost:63342/ninomoriaty/assembly-stats/assembly-stats.html?path=json/&assembly=Danaus_plexippus_v3&view=circle&altView=cumulative&altView=table

Mt-Practice

http://localhost:63342/ninomoriaty/assembly-stats/assembly-stats.html?path=json/&assembly=output2&view=circle&altView=cumulative&altView=table

http://localhost:63342/ninomoriaty/assembly-stats/assembly-stats.html?path=json/&assembly=Assembly_SRR9198433&view=circle&altView=cumulative&altView=table

### assembly_stats

**Ref**

- 

**Func**

- Quick Check 'Contig Stats' and 'Scaffold Stats' for a genome assembly

**Mt-Install**

```bash
# Tips: assembly_stats not assembly-stats
pip install assembly_stats
```

### GenomeQC

Ref

- GitHub: https://github.com/HuffordLab/GenomeQC

### SQUAT

**Ref**

- Article: https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-019-5445-3
- 



## Genome Annotation

<font color=orange size=6>Ref</font>

- Yandell, M. and Ence, D. (2012) A beginner’s guide to eukaryotic genome annotation. Nature Reviews Genetics 13: 329-342.
- Stanke, M., Diekhans, M., Baertsch, R. and Haussler, D. (2008) Using native and syntenically mapped cDNA alignments to improve de novo gene finding. Bioinformatics 24: 637-644.
- Guigó, R., Flicek, P., Abril, J.F., Reymond, A., Lagarde, J., Denoeud, F., Antonarakis, S., Ashburner, M., Bajic, V.B., Birney, E., Castelo, R., Eyras, E., Ucla, C., Gingeras, T.R., Harrow, J., Hubbard, T., Lewis, S.E. and Reese, M.G. (2006) EGASP: the human ENCODE Genome Annotation Assessment Project. Genome Biology 7(Suppl 1): S2.

------

<font color=orange size=6>Pcp</font>

- The optimal parameters for annotation are species specific.
- 

------

<font color=orange size=6> Cat-General </font>

Tips: General method usually combine ab initio and homology-based annotation

### NCBI Prokaryotic Genome Annotation Process

Ref: https://www.ncbi.nlm.nih.gov/genome/annotation_prok/process

### EnsEMBL genome annotation project

Ref: http://ensemblgenomes.org/info/data/annotation

#### Mt-BRAKER pipeline

**Ref**: Hoff et al. (2019) Whole-Genome Annotation with BRAKER. doi: 10.1007/978-1-4939-9173-0_5

**Cc**

- Training: an evaluation – perhaps one of a series – to help us set parameter values.
- Testing: a final evaluation to see how well our parameter values work
- Sensitivity measures the proportion of positives that are correctly identified.
- Specificity measures the proportion of negatives that are correctly identified.
- 

**Pcp**

- Use GeneMark-ET to produce initial set of protein coding genes for training
- Use AUGUSTUS to final gene prediction (machine learning)

**Pc**

1. GeneMark-ES/ET 
   1. Input: genome.fasta data and extrinsic data (like RNA-seq)
   2. output: genemark.gtf
2. Filter gene
3. AUGUSTUS Training
4. AUGUSTUS Prediction
5. Evaluation with pre-defined known gene
6. Evaluation with benchmark
   1. **BUSCO**: Benchmarking Universal Single-Copy Orthologs
      1. BUSCO genes are expected to be found in 1 copy per genome.
   2. **CEGMA**: Core Eukaryotic Genes Mapping Approach
      1. Each eukaryote is expected to have a set of core genes.
      2. Did our annotation correctly predict the CEGMA core?
7. Final output: GFF3 format for the annotation
8. Prove gene function
   1. Test of molecular function: like knockdowns etc. molecular experiments
   2. Test of adaptive function: field selection experiments

**Mt-Improve**

- Use a high quality and high contiguity genome
- Mask for repeats
- Consider organism of interest and clade specific features



### ANNOVA

Cat: General Gene Annotation and other function

------

<font color=orange size=6>Cat- Function annotation</font>

### Gene Ontology

Cat: 

### InterProScan

Cat:

------

<font color=orange size=6>Cat-de novo genome annotation</font>

### GeneMark-ET

Cat: de novo genome annotation

### SNAP

Cat: de novo genome annotation

### GenScan software

Cat: de novo genome annotation

------

<font color=orange size=6>Cat-Homology-based genome annotation</font>

### TBLASTN

Cat: Homology-based genome annotation

### GeneWise

Cat: Homology-based genome annotation

### AUGUSTUS software

Ref

- http://bioinf.uni-greifswald.de/augustus/

Cat: Homology-based genome annotation

### HMMER

Ref

- 

Func

- sequence homologs
- making sequence alignments



## Conservation Analysis

### For taxonomic group

### For genes





## EMBOSS

Tips: European Molecular Biology Open Software Suite

**Ref**:http://emboss.sourceforge.net/

**F**:

1. provide small and single-purpose bioinformatics tools
2. it provide online [GUI version](http://emboss.bioinformatics.nl/) and command line version

**Mt-Install**

```bash
sudo apt install emboss
```

**Mt-Plotcon**:

1. plot conservation of a sequence  alignment



## Entrez Programming Utilities

**Ref**: 

1. https://www.ncbi.nlm.nih.gov/books/NBK25501/
2. https://www.ncbi.nlm.nih.gov/books/NBK179288/

**Mt-install**

1. Current method can be found in https://www.ncbi.nlm.nih.gov/books/NBK179288/

2. But actually with Ubuntu, you can just install through `apt`

   ```bash
   sudo apt install ncbi-entrez-direct
   ```

   

3. 

**Mt-How to get all sequences for taxonomic groups with `wget`**

1. The `usehistory=y` parameter will generate the Web environment (`&WebEnv`) and query key (`&query_key`) parameters that will specify the location of the retrieved GI numbers (NOTE!) on the Entrez history server. 

2. Let's specify the retrieval type `rettype` as "fasta" and return mode `retmode` as "text". 

3. get the WebEnv and query_key

   ```bash
   # Set the taxonomic group and database
   beast="Cosmoscarta"
   db="nucleotide"
   
   # wget the WebEnv and query_key
   wget -qO- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=$db&term=$beast&usehistory=y" | awk '{
   ## Is this the third line
   if(FNR==3)
   {
   split($0,xml_line_split,"<|>");
   query_key=xml_line_split[17];
   WebEnv=xml_line_split[21];
   print query_key > "query_key";
   print WebEnv > "WebEnv";
   exit;}}'
   
   # get the query_key and WebEnv for the sequences download.
   query_key=$(cat query_key)
   echo $query_key
   WebEnv=$(cat WebEnv)
   echo $WebEnv
   
   # do the esearch component of the search and get the fasta format sequences
   wget -qO Seqs_I_wanted.fasta "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$db&query_key=$query_key&WebEnv=$WebEnv&rettype=fasta&retmode=text"
   ```

4. 

**Mt-Use the `edirect` method to get all sequences for taxonomic groups **

1. `esearch` and `efetch`

   ```bash
   # get the 
   esearch -db nucleotide -query "Cosmoscarta[organism]"
   
   # get the UID values
   esearch -db nucleotide -query "Cosmoscarta[organism]" | efetch -db nucleotide -format uid > Cosmoscarta.nuc.gis
   
   # get the accession values
   esearch -db nucleotide -query "Cosmoscarta[organism]" | efetch -db nucleotide -format acc > Cosmoscarta.nuc.acc
   
   # get the fasta sequences
   esearch -db nucleotide -query "Cosmoscarta[organism]" | efetch -db nucleotide -format fasta > Cosmoscarta.nuc.fa
   
   # get the record in GenBank format
   esearch -db nucleotide -query 1802791303 | efetch -db nucleotide -format gb > Cosmoscarta.first.gb
   
   ```

2. 

3. 



## Pullseq



## SAMtools

### SAMtools



Document: 

Mt-Common Usage

```bash

# Simple stat summary
samtools flagstat [file] # show the SIMPLE smmary statistics on the file

# Convert files Types
samtools view -b [file.sam] # SAM > BAM
samtools view -h [file.bam] # BAM > SAM
samtools fastq [file.bam] # BAM > FASTQ
samtools fasta [file.bam] # BAM > FASTA
```

Mt-Installation

1. http://www.htslib.org/download/

### BEDtools

Document:https://bedtools.readthedocs.io/en/latest/index.html

1. jkl

   ```bash
   ~/seqtk-master/seqtk trimfq SRR026762.fastq.gz | ~/seqtk-master/seqtk seq -a - > SRR026762.trimmed.fasta
   ```

2. 

### BCFtools

### HTSlib



### STAR



# Pipeline Programming Language

Tips: This section is not for the detailed pipeline constructions and it provides a index to link those pipelines constructors.

### Nextflow

**Lk**: https://www.nextflow.io/

**Mt-install**

```bash
# Remote usage
java --version # check if the version is sutiable for the usage of nextflow
curl -s https://get.nextflow.io | bash # download the makefile and use bash to execute the shell script to install and config
./nextflow run hello # test if the installation is complete
export PATH=/home/ninomoriaty/Nextflow:$PATH # add to .bashrc

# Download by Bioconda
conda install -c bioconda nextflow
```



### Snakemake

**Lk**: 

**Mt-install**

**Add-BBMap** 



# Physical Biology and Structure

### PyMOL

**Ref**: 

**Mt-Installation**

**Mt-Usage**



# Phylogenetic/Evolution Relationship

### X Mega / MEGA-X / MEGA 11

**Ref**

- Official website: https://www.megasoftware.net/

**Mt-Installation**

- Windows
  1. Find the 
- Ubuntu
  - Download from the web browser
    - Actually it is similar with the process in Windows
    - But the software across all platforms is only the MEGA X.
  - Download from the conda







### BEAST

Tips: Bayesian Evolutionary Analysis Sampling Trees

**Ref** 

- Main Page: 
- Community and Tutorial:https://beast.community/first_tutorial
- 

**Mt-Installation**

- Ref: https://beast.community/install_on_unix
- 



# Cluster for Community

### Netlogo

- 





<table><tr><td  style="background: #3ADDBA !important;"><center><font style="font-size:45px; color:black !important;"><b>
    Database in BioI
    </b></font></center></td></tr></table>

# 2 Database in BioI

# Integrated Database

### NCBI

Ref: http://www.ncbi.nlm.nih.gov/

Tutorial: https://www.ncbi.nlm.nih.gov/books/NBK1058/

#### NCBI taxonomy

1. **Ref**: http://www.ncbi.nlm.nih.gov/taxonomy
2. **F:** 
   1. Search database in the taxonomic group with the taxon name
   2. Get the taxonomy ID and it could be used in other database in NCBI
   3. 
3. 

#### GenBank

1. **Pcp-How to use the search syntax in NCBI**
   1. [] could limit the search range
      1. [Organism:exp] protein sequence
      2. [Gene Name] gene name
      3. NOT Partial will not include the partial gene/sequence(in another word, it will be complete)
   2. () could prioritise the statement
   3. logic operators
      1. AND
      2. OR
      3. NOT
   4. 
2. **Cpn-NCBI data information page**
   1. **Metadata**: The information about the experiment of the dataset and the sample information of the result
   2. **Annotation**: The ontological and functional annotation and index of the sequence
   3. **Sequences**: The sequences of the sample.
   4. 
3. **F:**
   1. Use the taxonomy ID to find query sequences



### European Institute of Bioinformatics(EMBL-EBI)

**Ref**

- https://www.ebi.ac.uk/





# Nucleotide Database

### European Nucleotide Archive(ENA)

**Ref**

- https://www.ebi.ac.uk/ena/browser/home

**Obj2**

- Collect free+public+annotated DNA and RNA sequences
- 



### Ensemble BioMart



# Protein Database

### PROSITE



### UniPort/Swiss Port

Cc



Ab-Total Database

- Protein Resource Database
- Data is manually curated by domain experts
- 



Ab-Every part of Annotation Data

- Protein first/second/quaternary structure
- Protein domain, sites
- protein function
- protein post translation modifications
- Similarity (Homologene)
- Potential Sequence variants 
- Disease
- Gene Ontology Information



### PDB

**Ref**:https://www.rcsb.org/



# Functional omics Database

### GEO



### KEGG Pathway



# Drug Discovery Database

### Drug Bank

**Ref**: https://go.drugbank.com/



### ChEMBL

**Ref**: https://www.ebi.ac.uk/chembl/



### ZINC

**Ref**: http://zinc15.docking.org





<table><tr><td  style="background: #3ADDBA !important;"><center><font style="font-size:45px; color:black !important;"><b>
    Data Structure in BioI
    </b></font></center></td></tr></table>

# 3 Data Structure in BioI

# Pcp-codes and symbols

## **IUPAC Ambiguity Codes**:

Tips: If needed, you can remove the nucleotide codes which don't list in this sheet.

| IUPAC Code | Meaning          | Complement |
| ---------- | ---------------- | ---------- |
| A          | A                | T          |
| C          | C                | G          |
| G          | G                | C          |
| T/U        | T                | A          |
| M          | A or C           | K          |
| R          | A or G           | Y          |
| W          | A or T           | W          |
| S          | C or G           | S          |
| Y          | C or T           | R          |
| K          | G or T           | M          |
| V          | A or C or G      | B          |
| H          | A or C or T      | D          |
| D          | A or G or T      | H          |
| B          | C or G or T      | V          |
| N          | G or A or T or C | N          |







# Fomats

## GenPept

**Cpn**

1. Metadata
2. annotation
3. sequence

**Ab**

1. A complex and strictly defined structure



## GenBank(gb)



**Cpn**

1. CDS: the range of coding segment
2. 



**e.g.**

```
LOCUS       AJ223353                 808 bp    mRNA    linear   PRI 07-OCT-2008
DEFINITION  Homo sapiens mRNA for histone H2B, clone pJG4-5-15.
ACCESSION   AJ223353
VERSION     AJ223353.1  GI:3255998
KEYWORDS    histone H2B.
SOURCE      Homo sapiens (human)
  ORGANISM  Homo sapiens
            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;
            Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;
            Catarrhini; Hominidae; Homo.
REFERENCE   1
  AUTHORS   Lorain,S., Quivy,J.P., Monier-Gavelle,F., Scamps,C., Lecluse,Y.,
            Almouzni,G. and Lipinski,M.
  TITLE     Core histones and HIRIP3, a novel histone-binding protein, directly
            interact with WD repeat protein HIRA
  JOURNAL   Mol. Cell. Biol. 18 (9), 5546-5556 (1998)
   PUBMED   9710638
REFERENCE   2  (bases 1 to 808)
  AUTHORS   Lipinski,M.
  TITLE     Direct Submission
  JOURNAL   Submitted (08-JAN-1998) Lipinski M., CNRS URA 1156, Institut
            Gustave Roussy, Rue Camille Desmoulins, 94805, FRANCE
FEATURES             Location/Qualifiers
     source          1..808
                     /organism="Homo sapiens"
                     /mol_type="mRNA"
                     /db_xref="taxon:9606"
                     /clone="pjG4-5-15"
                     /cell_line="HeLa"
                     /tissue_type="cervix carcinoma"
                     /dev_stage="adult"
     gene            1..808
                     /gene="HIRIP2"
     CDS             29..409
                     /gene="HIRIP2"
                     /codon_start=1
                     /product="Histone H2B"
                     /protein_id="CAA11277.1"
                     /db_xref="GI:3255999"
                     /db_xref="GOA:P58876"
                     /db_xref="HGNC:HGNC:4747"
                     /db_xref="InterPro:IPR000558"
                     /db_xref="InterPro:IPR007125"
                     /db_xref="InterPro:IPR009072"
                     /db_xref="UniProtKB/Swiss-Prot:P58876"
                     /translation="MPEPTKSAPAPKKGSKKAVTKAQKKDGKKRKRSRKESYSVYVYK
                     VLKQVHPDTGISSKAMGIMNSFVNDIFERIAGEASRLAHYNKRSTITSREIQTAVRLL
                     LPGELAKHAVSEGTKAVTKYTSSK"
     regulatory      769..774
                     /regulatory_class="polyA_signal_sequence"
                     /gene="HIRIP2"
ORIGIN      
        1 ggggagtgtg ctaactatta acgctacgat gcctgaacct accaagtctg ctcctgcccc
       61 aaagaagggc tccaagaagg cggtgactaa ggctcagaag aaggacggga agaagcgcaa
      121 gcgcagccgc aaggagagct attcagtgta tgtgtacaag gtgctgaagc aggtccatcc
      181 cgacaccggc atctcttcca aggcaatggg gatcatgaat tccttcgtca acgacatctt
      241 cgagcgcatc gcaggcgagg cttcccgcct ggcgcattac aacaagcgct cgaccatcac
      301 ctccagggag atccagacgg ccgtgcgcct gctgcttccg ggggagctgg ccaagcacgc
      361 cgtgtcggag ggcaccaagg ccgtcaccaa gtacaccagt tccaagtaac tttgccaagg
      421 gagagacatg aagacagagg agaaatgaat gcataaaata actgataata tgaatctata
      481 catagaactt aggaagtctc atctgcctga aaatgactgt gtggatccca cccaaatcca
      541 actcatcctg gtttgctgca cactggttca tcaaaagaag gttaccgagg ggaaggaact
      601 aaaggtgttt gcacttcatg ttactttttg agtttataaa cataaaaaca gaatttactt
      661 ctgttacaga cctagttact gggaattcat tacttgccat ggactacctt tgctaagaaa
      721 agtctgaatg agaagatggc aggacgtctg aaaaaaaaag ttataattaa taaaatctgc
      781 ggagaattgt aaaaaaaaaa aaaaaaaa
//

```





## FASTA format

**Cpn**

1. Brief description/header(with the product information and organism details)
2. Sequnces

**Ab**

1. Brief and standard for many downstream analysis

**e.g.**

```
>gi|3255998|emb|AJ223353.1| Homo sapiens mRNA for histone H2B, clone pJG4-5-15
GGGGAGTGTGCTAACTATTAACGCTACGATGCCTGAACCTACCAAGTCTGCTCCTGCCCCAAAGAAGGGC
TCCAAGAAGGCGGTGACTAAGGCTCAGAAGAAGGACGGGAAGAAGCGCAAGCGCAGCCGCAAGGAGAGCT
ATTCAGTGTATGTGTACAAGGTGCTGAAGCAGGTCCATCCCGACACCGGCATCTCTTCCAAGGCAATGGG
GATCATGAATTCCTTCGTCAACGACATCTTCGAGCGCATCGCAGGCGAGGCTTCCCGCCTGGCGCATTAC
AACAAGCGCTCGACCATCACCTCCAGGGAGATCCAGACGGCCGTGCGCCTGCTGCTTCCGGGGGAGCTGG
CCAAGCACGCCGTGTCGGAGGGCACCAAGGCCGTCACCAAGTACACCAGTTCCAAGTAACTTTGCCAAGG
GAGAGACATGAAGACAGAGGAGAAATGAATGCATAAAATAACTGATAATATGAATCTATACATAGAACTT
AGGAAGTCTCATCTGCCTGAAAATGACTGTGTGGATCCCACCCAAATCCAACTCATCCTGGTTTGCTGCA
CACTGGTTCATCAAAAGAAGGTTACCGAGGGGAAGGAACTAAAGGTGTTTGCACTTCATGTTACTTTTTG
AGTTTATAAACATAAAAACAGAATTTACTTCTGTTACAGACCTAGTTACTGGGAATTCATTACTTGCCAT
GGACTACCTTTGCTAAGAAAAGTCTGAATGAGAAGATGGCAGGACGTCTGAAAAAAAAAGTTATAATTAA
TAAAATCTGCGGAGAATTGTAAAAAAAAAAAAAAAAAA
```









## FASTQ format

**Soc**:[Illumina's fastq-files-explained.html](https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html)

**e.g.**

```
@K00114:84:H2N7KBBXX:8:1101:6105:1457 1:N:0
ATTTCCTTTTCGTGGATACGACANAGGATGATGCAAGCAGCGAGTAAACGATTTGATGAGGCTTAGCGCATCATGAATAGCACAGATGATGTGGAGAAAG
+
AAFFFKFFKKKKKKKKKKKFKKK#KFKFFFAKKKKF7,AFAKAAKKKKKK,FAFFFKA7FFKKFK,,7AFFKAKKKKKF,7FKK7AFFF,A,,,7AFKF7
@K00114:84:H2N7KBBXX:8:1101:6856:1457 1:N:0
TGAAGAAATAAAACAAGTAACAGNGACACCAACTCAGCTATCAGCCGCTTCCGTGATTAGACGCAGCCGTTCAACGGCCATTGCGAGCAAGGCAAATGTA
+
AAFFFKKKKKKKKKKKKKKKKKK#KKKKKKKKKKKKKKKKKFKFFAFFKKKKKKKKKKKFKKFFAF,(7AFAFFKFFKFFKKKKKK7FKKKKKKKKKA,A
...
```



## BLAST output

**e.g.**

```
# BLASTX 2.2.28 + 
# Query: gi|322830704:1426-2962 Boreus elegans mitochondrion, complete genome
# Database: nem
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 405 hits found
gi|322830704:1426-2962  gi|225622197|ref|YP_002725710.1|        61.71   444     169     1       4       1335    12      454     4e-145   432
gi|322830704:1426-2962  gi|288903410|ref|YP_003434131.1|        61.99   442     167     1       10      1335    14      454     2e-144   430
gi|322830704:1426-2962  gi|326536058|ref|YP_004300493.1|        62.44   442     165     1       10      1335    14      454     2e-144   430
gi|322830704:1426-2962  gi|188011119|ref|YP_001905892.1|        61.31   442     170     1       10      1335    14      454     4e-144   429
gi|322830704:1426-2962  gi|225622184|ref|YP_002725698.1|        61.49   444     170     1       4       1335    12      454     6e-144   429
gi|322830704:1426-2962  gi|171260186|ref|YP_001795390.1|        61.54   442     169     1       10      1335    15      455     6e-144   429
gi|322830704:1426-2962  gi|288903312|ref|YP_003434040.1|        61.99   442     167     1       10      1335    14      454     6e-144   429
gi|322830704:1426-2962  gi|49146527|ref|YP_026087.1|    61.54   442     169     1       10      1335    14      454     9e-144   428
gi|322830704:1426-2962  gi|296315379|emb|CAQ56159.1|    61.54   442     169     1       10      1335    14      454     1e-143   428
gi|322830704:1426-2962  gi|160694986|ref|YP_001552232.1|COX1_21467      61.99   442     167     1       10      1335    14      454     1e-143   428
...

```



## BED format





## SAM format

Soc:

1. SAMtools
2.  https://samtools.github.io/hts-specs/SAMv1.pdf

Soc-Produce SAM file

1. bowtie, bwa, STAR, etc. will output in SAM
2. 

Soc-Need SAM file

Cpn

1. header section (@): contigs in the reference genome, mapped software
2. data section: sequences listed in lines

F:

1. capture sequences and quality information, as found in FASTQ files
2. capture how the reads align to a reference genome

e.g.

```
@HD VN:1.5 SO:coordinate
@SQ SN:ref LN:45
r001   99 ref  7 30 8M2I4M1D3M = 37  39 TTAGATAAAGGATACTG *
r002    0 ref  9 30 3S6M1P1I4M *  0   0 AAAAGATAAGGATA    *
r003    0 ref  9 30 5S6M       *  0   0 GCCTAAGCTAA       * SA:Z:ref,29,-,6H5M,17,0;
r004    0 ref 16 30 6M14N5M    *  0   0 ATAGCTTCAGC       *
r003 2064 ref 29 17 6H5M       *  0   0 TAGGC             * SA:Z:ref,9, + ,5S6M,30,1;
r001  147 ref 37 30 9M         =  7 -39 CAGCGGCAT         * NM:i:1
...
```



## BAM format

Def: Native compressed format of SAM

F

1. perform more efficient operations on the BAM files than SAM files



## VCF format

**Def**: Variant Call Format is a text file format (often stored in a compressed manner as .vcf.gz) to store variants

**Soc**: https://github.com/samtools/hts-specs

**Cpn**:

1. meta-information lines which introduce the sample information

2. 1. Ab: preceded with ##

3. a header line + 3 = one data section

4. 1. Cpn:

   2. 1. CHROM: chromosome
      2. POS: position
      3. ID: variant identifier
      4. REF: reference base(s)
      5. ALT: alternate base(s)
      6. QUAL: quality (not used in our file)
      7. FILTER: filter status (not us ed in our file)
      8. INFO: additional information – in our file this field contain the variant type (SNV, insertion, or deletion)

   3. Pcp: 

   4. 1. insertion and deletion's REF and ALT will start with the base that is present before the insertion or deletion event

5. data lines + 2 = one data section 

6. genotype information


**e.g.**

```
##fileformat=VCFv4.0
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=1000GenomesPilot-NCBI36
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
20     1234567 microsat1 GTCT   G,GTACT 50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3
...

```





## CCF format

## MAF format

 



<table><tr><td  style="background: #3ADDBA !important;"><center><font style="font-size:45px; color:black !important;"><b>
    High-Performance Computation
    </b></font></center></td></tr></table>

# 3 High-Performance Computation

# Super Computer

## Tianhe2

Source: Guangzhou SYSU

## Cirsus



# Cluster Server

## Eddie

**Ref**: https://www.wiki.ed.ac.uk/display/ResearchServices/Eddie

**Source**: the University of Edinburgh

**Def**: Eddie Mark 3 is the third iteration of the University's <u>compute cluster</u> and is available to all University of Edinburgh researchers. It consists of over 7000 cores with up to 3 TB of memory.

**Cpn**:

1. **node**:
   1. Cat-login nodes and clusters
      1. login nodes: The nodes for interact with the Eddie server and access to the cluster by queuing system
      2. Clusters/Working nodes: The nodes for executing the tasks submitted by users
      3. Grid Engine/Queuing System: The system for arrange the working sequences for different tasks.
   2. Cat-level of cores and memory
      1. Soc-Service: https://www.wiki.ed.ac.uk/display/ResearchServices/Eddie
      2. Soc-Storage: https://www.wiki.ed.ac.uk/display/ResearchServices/Storage
      3. Lk: you can find the "The list of all nodes available as well as their specifications" in the Eddie Wiznote in BPSM folder
      4. standard
         1. 16 core per nodes
         2. 64 GB
      5. There will be a list of the nodes that could be used in the official website 
2. **Cores**: The processing units(Threads) of each of the computers in the cluster
   1. Available memory for each core is the average of the whole memory
3. **Memory:** RAM/memory

**Ab**:

1. Parallelisation: divide the tasks into a number of cores.
2. Different research groups have different Eddie directory and different DataStore space.

**Mt-The workflow of Eddie**:

P:

1. You should be within the university network. VPN and other methods are also acceptable.
2. 

Pc:

1. **Login in the Eddie**

   ```bash
   ssh -X s2059232@eddie.ecdf.ed.ac.uk
   ssh -X s2059232@129.215.120.98
   ```

   F:

   1. The login space is the home directory for your account/UUN and the documents here will be backed up (up to 2GB space)
   2. The scratch space use up to 20-50TB  and it is not backed up, which will be clean after one month.
   3. 

2. **Copy the dataset to the scratch space/Move files from MSc server**

   1. Lk: DataStore: university's file system for storing research datahttps://www.wiki.ed.ac.uk/display/ResearchServices/Data+Staging

      ```bash
      scp [options] [dir/file] [destination_dir]
      ## e.g.
      # copy to the home dir
      scp ./environment.yml s2059232@eddie.ecdf.ed.ac.uk:/home/s2059232
      # copy to the scratch space
      scp ./environment.yml s2059232@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s2059232/
      # copy to the home dir
      scp s2059232@eddie.ecdf.ed.ac.uk:/home/s2059232/environment.yml ./
      # copy from the scratch space
      scp s2059232@eddie.ecdf.ed.ac.uk:/exports/eddie/scratch/s2059232/file ./
      # move the folder and files in a dir to Eddie
      scp -r ./environment.yml s2059232@eddie.ecdf.ed.ac.uk:/home/s2059232/
      
      # options
      -r # recursive
      
      # Tips: After enter the interactive node as ii part, you can use the Unix command to move files from Eddie to DataStore.
      ```

   2. Att: BEFORE you log in to eddie (or from another terminal), you will need to scp a file from the bioinformatics server to eddie. 

   3. Att: Other Edinburgh University storage solutions can be used to store files and move the files to Eddie with a staging node.

3. **Start an interactive session and login in the login nodes:**

   Tips: Lk about how to settle the interactive sesssion:https://www.wiki.ed.ac.uk/display/ResearchServices/Interactive+Sessions 

   ```bash
   # interactive sessions within the Eddie
   qlogin
   # indicate the core and memory th
   qlogin -pe sharedmem < NUMBER OF CORES > -l h_vmem=< MEMORY PER CORE >
   # Log into a staging node
   qlogin -q staging
   ```

   1. Load modules

      ```bash
      # usage
      module load <modules needed to use in Eddie>
      module available # check available modules
      ## e.g.
      module load htop
      module load python
      
      ```

      F: However, this function could not be tested in my computer thus you may need to try it when the section is needed.

4. **Write and configure a shell script**

   1. Lk: https://www.wiki.ed.ac.uk/display/ResearchServices/Job+Submission

   2. Mt: nano or nedit

   3. Mt-Add task description for queuing system

      Tips: This description must be added to the top of the script and start with `#$`. But **after the shebang line.**

      ```latex
      #!/bin/bash %shebang line should be in the first line as well
      
      % The Grid Engine Instrcution Part
      #$ -cwd % Tell Eddie to run the code in the current working directroy
      #$ -t 1-100 % Grid Engine parameter -t which will give each of our jobs a number (1-100), called the SGE_TASK_ID (for files with numeric suffixes or prefixes and the $SGE_TASK_ID can be called in the shell script. Also ues to parallelise these tasks and consider them as a collection of tasks.)
      #$ -N Hello % Give a name for the task and Eddie will also give the number code for the task name
      #$ -l h_rt=01:01:10 % Specify the time for the job. h_rt means the hoour time: 1 hr and 1 min and 1 sec.
      #$ -l h_vmem=1G % Specify the total memory(seems that it is not the memory needed for each core?)
      
      % The main program
      ####### your commands or original shell script ######
      
      % This output part should be placed after the main project.
      #$ -o hello.o % Specify the output directory and output file
      #$ -e hello.e % Specify the error directory and error report file
      
      ```

      Tips: To make the command after `#$` visible, we use latex as the highlight marker of code section. However, these code sections should be considered as part of the shell script.

   4. Tips: but actually you can use the method in the i part to move scripts into Eddie and execute them as Unix command line. 

   1. Add the task into the Eddie clusters and queuing system

5. **Submit and complete the analysis**

   ```bash
   # Check the manual of the qsub
   man qsub
   # Submit the shell script and wait for the execution
   qsub <your_script_or_other_execuatble_files>
   # Display the status of the Eddie
   qstat # The state column tells us that there are six running (r state) and the rest are qwstate (queued and waiting). If no output, it means that there is no works for Eddie to do in this accout.
   # Kill a job you subimitted to Eddie by the task number or the job name
   qdel <job_number_or_name>
   # After the job has been finised, get useful info from the Eddie report
   qacct -j <job_number_or_name>
   ```

   

   1. 

6. Copy the results elsewhere

7. Ask to work inside a node interactively

8. Check the processes

9. 































































































































































