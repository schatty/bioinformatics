# Antibiotic resistance factors analysis with alignment to reference and variant calling

__Daria Tsyba, Igor Kuznetsov__

## Abstract

Microbial resistance to antibiotics is a major problem of bacterial infections. The ability to identify proteins that affects bacteria resistance is crucial to develop new drugs. Presented study proposes the computational approach to analyse mutated genome of antibiotic resistant bacteria. We demonstrate results on the strain of E. coli resistant to the ampicillin. We show that use of such bioinformatics techniques as alignment to reference and variant calling is efficient to locate mutated genes that may cause antibiotic resistance. Proposed sequence alignment pipeline has a potential to be used for identifying mutated genes in similar antibiotic-resistance analysis problems.

## Introduction

To fight against microbes and infections antibiotics are the vital component of the modern medicine that can save millions of lives [1]. However, due to the enormous use of the antibiotics, the advent of many resistant strains can be witnessed. As a result much of the modern approach in fighting against the diseases is focused on modification of existing antibiotics. To facilitate of the development of new drugs it is important to understand what bacteria genes are responsible for certain antibiotic resistance. 

The computational approach of bioinformatics can shed light on difference between original and mutated genomes of certain bacterias to infer proteins responsible to the drug resistance [2]. In present study we demonstrate how mutated genes can be located using such bioinformatics techniques as sequence alignment and variant calling. Proposed pipeline is demonstrated on E.coli resistant to the antibiotic ampicillin and facilitate to locate genes and associated proteins that causes resistance to the antibiotic. 

The rest of the paper organized as follows. Section Materials and methods describes tools and techniques used to identify mutated genes. Results section exposes numerical results of sequence preprocessing and variant calling steps as well as disclosure of the role of mutated proteins for the drug resistance. Discussion section affects results interpretation.

## Materials and Methods

Given the reference parental sequence of the bacteria and sequencing reads of the same bacteria it is possible to locate genes that mutated and therefore has potential to be important for the antibiotic resistance. 

The parental sequence of the E.coli was obtained from NCBI database [3]. Forward and reverse strands for the same bacteria was obtained from Theodosius Dobzhansky Center for Genome Bioinformatics database [4].

We use fastqc tool [5] for evaluating basic fastq statistics of the given raw Illumina sequencing reads from shotgun sequencing routine. The initial analysis of the reads showed that quality scores across 96-97 bases are lover than 20-score for both of the forward and reverse strands which can be considered as anomaly. Also anomaly for the per tile sequence quality was encountered for the forward strand. We use Trimmomatic [6] program to improve the overall quality of the sequencing reads before proceeding with following analysis. The method is to remove low quality base calls in the given reads. The trimming procedure of Trimmomatic was set to cut bases of the start and the end of a read if quality was below 20. Trimming was performed with sliding window approach, with window size of 10 and average quality within the window of size 20. Also all reads with a length below 20 were dropped. The size of 20 was chosen as additional experiments showed that bigger trimming parameter cause too much reads loss (18.6% of the reads were dropped with trimming parameter set to 30, whereas only 3.5% with trimming parameter of 20).

To match sequence reads to the original genome we use the alignment procedure based on Burrows-Wheeler Transform algorithm [7]. We used BWA tool [8] as it is optimized for ‘long’ next generation sequencing reads starting from 100bp. This reads can contain several mutations, insertions or deletions which can be hard to process with earlier Needleman-Wunsh based alignment tools [9]. The reference file was indexed with bwtsw algorithm and subsequently aligned with bwa-mem tool.

To estimate the percentage of the mapped reads we use module flagstat of the samtools [10] program. Before following steps obtained, bam file was sorted and indexed with default parameters of the samtools. The next step of the analysis concerns distinguishing actual mutations from the sequencing errors. In order to accomplish this computational procedure we made an intermediate file called mpileup with samtools program.

For the final step of calling actual variants we use program VarScan [11] given mpileup file as input. We use 50% threshold as minimum percentage of on-reference bases as a position required to call it mutation in the sample. This threshold is chosen as sequences are from isogenic population that resistant to ampicillin, therefore original SNP are expected in majority of genes. 

## Results

In the step of sequence reads preprocessing we use trimming technique to improve the overall quality of our reads. The figure 1 shows the difference between forward and reverse reads before and after trimming procedure regarding per base quality. Per tile quality depicted on figure 2.

![title](https://github.com/schatty/bioinformatics/blob/master/mini-projects/antibiotic-resistance-factors/init_per_base_seq_quality_1.png)

Evaluation of the tripping procedure improve quality of the reads with respect to the per base sequence quality. However, per tile quality for the forward strand remained to be considered as anomaly. The use of the tripping procedure caused little reduction in number of reads. The threshold was set in a way to cause no more than 4% reads reduction. The exact numberers of reads remaining after processing are shown in the table 1.

| Type                   | Number                 |
| ---------------------- | ---------------------- |
| Raw                    | 455876                 |
| Remain after trimming  | 439767                 |
| Aligned reads          | 439195                 |

The genes that were mutated and located with variant calling are ftsI, acrB, rybA (mntS), envZ, rsgA. Mutation of rsgA is synonymous and other mutations are missense: G544A in ftsI, G569L in acrB, V241G in envZ. The mutation in rybA (mntS) affects non-coding part of gene.

## Discussion

The most evident gene involved in mechanism of antibiotics resistance in our case is ftsI. This gene encodes peptidoglycan D,D-transpeptidase (penicillin-binding protein-3, PBP-3) - essential enzyme for cross-linking of the peptidoglycan in cell wall at the division septum.  Beta-lactam antibiotics like ampicillin inhibit the activity of this transpeptidase by binding to the catalytic serine [12]. So changes in catalytic site structure can damage its interaction with beta-lactams and be the source of antibiotic resistance due to target site alteration. During this study we find one missense mutation in ftsI gene leading to G544A replacement in its product. This replacement affects periplasmic domain - the crucial part of PBP-3. The role of ftsI mutations in development of bacterial resistance to beta-lactam antibiotics has been shown earlier [13]. 

The second potential candidate for gene of resistance is acrB - its product is a part of a drug efflux protein complex AcrA-AcrB-AcrZ-TolC with broad substrate specificity that uses the proton motive force to export substrates [14]. Previously the association between this efflux system and beta-lactam resistance has been reported but with more effectiveness of AcrD subunit [15]. We suppose replacement G569L in acrB product may increase its selectivity to beta-lactams and take part in resistance formation. 

EnvZ is a membrane-associated protein kinase/phosphatase in Escherichia coli involved in osmoregulation due to regulating the phosphorylation state of the transcription factor OmpR [16]. Despite of non evident relation between osmoregulation and antibiotic resistance combination of mutations in envZ, ftsI, acrB and some other genes has been reported to be the potential cause of high-level resistance to other beta-lactam antibiotic - carbapenem [17]. 

There is no evidence of rybA or rsgA mutations influence on ampicillin resistance in our research.  Consequently, the most probable causes of ampicillin resistance in experimental E.coli are target site modification and drug efflux.

Because of common usage of ampicillin clinical guidelines in this case depend on site of E.coli infection, presence of other pathogens (mixt-infection), patient age, renal and hepatic function and many-many other clinical factors. But the main recommendation is to change antibacterial therapy to other class of antibiotic drugs - non-beta-lactam. Switching within beta-lactam class isn’t effective because of changed target site in PBP-3 and increased drug efflux. Administration of fluoroquinolones, aminoglycosides or macrolides for empiric therapy should be considered. Of course, before the start of empiric therapy samples should be obtained and antimicrobial susceptibility testing with other available antibiotics should be performed for effective definitive antimicrobial therapy later. The common recommendations for clinicians include rational usage of antibiotics, prescription of an appropriate dose, observation over complement course of antibiotic treatment and avoiding reserve antibiotics without invariable indications [18].

## Resources

* [1] Levy S.B. From Tragedy the Antibiotic Age is Born. In: The Antibiotic Paradox. Springer, Boston, MA. (1992) pp. 1–12
* [2] Su M., Satola S.W., Read T.D. Genome-based prediction of bacterial antibiotic resistance. J Clin Microbiol. (2019) 57:e01405–18. 10.1128/JCM.01405-18
* [3] https://www.ncbi.nlm.nih.gov/genome/167
* [4] http://public.dobzhanskycenter.ru/mrayko/amp_res_1.fastq.zip, http://public.dobzhanskycenter.ru/mrayko/amp_res_2.fastq.zip 
* [5] Andrews S. FastQC: a quality control tool for high throughput sequence data. (2010) Available online at:http://www.bioinformatics.babraham.ac.uk/projects/fastqc
* [6] Bolger A.M., Lohse M., Usadel B. Trimmomatic: A flexible trimmer for Illumina Sequence Data. (2014) Bioinformatics, btu170.
* [7] The Burrows-Wheeler Transform: Data Compression, Suffix Arrays, and Pattern Matching (1 ed.). Springer Publishing Company, Incorporated.
* [8] Li H., Durbin R. Fast and accurate long-read alignment with Burrows-Wheeler Transform. Bioinformatics, Epub.  (2010) http://bio-bwa.sourceforge.net
* [9] Needleman S.B., Wunsch C.D. A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol (1970) 48(3):443-453
* [7] Li H., Handsaker B., Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics. (2009) 25, 2078-9. http://www.htslib.org/ 
* [11] Koboldt D., Zhang Q., Larson D., Shen D., McLellan M., Lin L., Miller C., Mardis E., Ding L., Wilson R. VarScan 2: Somatic mutation and copy number alteration discovery in cancer by exome sequencing Genome Research. (2012) URL: http://varscan.sourceforge.net 
* [12] Nguyen-Disteche M., Fraipont C., Buddelmeijer N., Nanninga N. The structure and function of Escherichia coli penicillin-binding protein 3. Cell. Mol. Life Sci. (1998) 54:309-316
* [13] Hedge P.J., Spratt B.G. Resistance to beta-lactam antibiotics by re-modelling the active site of an E. coli penicillin-binding protein. Nature. 1985 Dec 5-11;318(6045):478-80.
* [14] Hobbs E.C. et al. Conserved small protein associates with the multidrug efflux pump AcrB and differentially affects antibiotic resistance. Proc Natl Acad Sci U S A. 2012 Oct 9;109(41):16696-701.
* [15] Kobayashi N. et al. β-Lactam selectivity of multidrug transporters AcrB and AcrD resides in the proximal binding pocket. J Biol Chem. 2014 Apr 11;289(15):10680-90.
* [16] Cai S.J., Inouye M. EnvZ-OmpR interaction and osmoregulation in Escherichia coli. J Biol Chem. 2002 Jul 5;277(27):24155-61.
* [17] Adler M. et al. Combinations of mutations in envZ, ftsI, mrdA, acrB and acrR can cause high-level carbapenem resistance in Escherichia coli. J Antimicrob Chemother. 2016 May;71(5):1188-98.
* [18] Leekha S., Terrell C.L., Edson R.S. General Principles of Antimicrobial Therapy. Mayo Clin Proc. 2011 Feb; 86(2): 156–167.
