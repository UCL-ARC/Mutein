# Summary of Variant Calling Methods use in the Dataset Publications

What steps used in their variant calling pipelines, which tools, did they has matched normal (eg blood) samples? Technical notes on the datasets, such as which panel of genes did they use for target capture.

### deepSNV R package
https://github.com/gerstung-lab/deepSNV seems to be the current official home of the package.
https://github.com/gerstung-lab/deepSNV/tree/master/vignettes details three different algorithm: deepSNV, Shearwater and ShearwaterML. Martincorena2015 used ShearwaterML.

### Martincorena2015
Eyelid surgery samples from 4 adults 55-73 years old. Indiv 1: 24 samples. Indiv 2: 92 samples. Indiv 3:91 samples. Indiv 4: 28 samples.
Reference GRCh37. A supplement Excel sheet lists all the somatic mutations found. "The R code used for variant calling (ShearwaterML) is available in
Bioconductor (deepSNV package)".

From the supplementary material pdf:
```
A panel of 74 genes (see list below) was chosen to perform ultra-deep targeted
sequencing. These genes were chosen based on the following criteria: (a) genes often
involved in cutaneous squamous cell carcinomas, basal cell carcinomas or melanomas,
(b) genes found to be involved in a wide range of cancers, or (c) genes frequently
mutated in skin samples from the COSMIC database (51), even if not known to play a
driver role in cancer:
ADAM29, ADAMTS18, AJUBA, AKT1, AKT2, APOB, ARID1A, ARID2, AURKA,
BAI3, BRAF, CASP8, CCND1, CDH1, CDKN2A, CR2, CREBBP, CUL3, DICER1,
EGFR, EPHA2, ERBB2, ERBB3, ERBB4, EZH2, FAT1, FAT4, FBXW7, FGFR1, FGFR2,
FGFR3, FLG2, GRIN2A, GRM3, HRAS, IRF6, KCNH5, KEAP1, KRAS, MET, MLL,
MLL2, MLL3, MUC17, NF1, NFE2L2, NOTCH1, NOTCH2, NOTCH3, NOTCH4, NRAS,
NSD1, PCED1B, PIK3CA, PLCB1, PPP1R3A, PREX2, PTCH1, PTEN, PTPRT, RB1,
RBM10, SALL1, SCN11A, SCN1A, SETD2, SMAD4, SMO, SOX2, SPHKAP, SUFU,
TP53, TP63 and TRIOBP
A custom bait capture was designed using Agilent SureDesign, targeting the exonic
sequences of these 74 genes. In addition, probes were designed against 610 SNPs
regularly scattered throughout the genome and against 1,124 SNPs within or around the
74 selected genes, to perform copy number analysis (see supplementary methods S1.6).
Repetitive or low complexity regions were carefully removed to maximise on-target
coverage. The total size of the targeted regions was 0.67 Mb.
Sequencing of paired-end 75bp reads was performed on Illumina HiSeq 2000 or
2500 machines. After removing reads for off-target capture (mean 32.6%) and PCR
duplicates, the average on-target coverage across samples was 500.2x (ranging from
107.9x to 1163.0x; first quartile 373.7x and third quartile 598.1x).
```

Data processing pipeline:

- BWA align paired end reads
- Picard MarkDuplicates
- ShearWaterML

Is that it? Not sure if the actual scripts are available or not.

```
Typically, somatic mutations are called by detecting mismatches present in a tumor
sample that are absent in a matched normal sample. The matched normal sample is often
a sample of normal tissue sequenced to moderate coverage (typically 20-40X), often
blood from convenience but not necessarily so. Here, instead of using a single matched
normal sample as a reference, mutations in each sample were called against all other
normal samples from the same patient.
```

Perhaps they combined all the other samples, so that they only get one set of calls per samples,
rather than doing all-versus-all like we are doing!


### Buscarlet2017

### Martincorena2018
esophagus harvested from nine deceased organ donors aged 20-75 while other organs were being harvested.

From the supplement:
```
we used an Agilent SureSelect custom bait capture design covering 74 cancer
genes. In addition, we targeted 610 SNPs regularly scattered across the genome for copy number
analysis and 1,124 SNPs within or around the 74 target genes for targeted copy number analysis.
This was the same design used in our previous study on sun-exposed skin.

The list of genes selected for ultra-deep targeted sequencing is shown below:
ADAM29, ADAMTS18, AJUBA, AKT1, AKT2, APOB, ARID1A, ARID2, AURKA, BAI3, BRAF,
CASP8, CCND1, CDH1, CDKN2A, CR2, CREBBP, CUL3, DICER1, EGFR, EPHA2, ERBB2,
ERBB3, ERBB4, EZH2, FAT1, FAT4, FBXW7, FGFR1, FGFR2, FGFR3, FLG2, GRIN2A,
GRM3, HRAS, IRF6, KCNH5, KEAP1, KMT2A, KMT2C, KMT2D, KRAS, MET, MUC17, NF1,
NFE2L2, NOTCH1, NOTCH2, NOTCH3, NOTCH4, NRAS, NSD1, PCED1B, PIK3CA, PLCB1,
PPP1R3A, PREX2, PTCH1, PTEN, PTPRT, RB1, RBM10, SALL1, SCN11A, SCN1A, SETD2,
SMAD4, SMO, SOX2, SPHKAP, SUFU, TP53, TP63 and TRIOBP
```

```
Paired-end reads were aligned with BWA (41) and PCR duplicates were marked
using Pircard (http://broadinstitute.github.io/picard/). We then performed indel realignment on
the resulting bam files using IndelRealigner from GATK
```

```
ShearwaterML was used in the way described in our previous study (7), with three
modifications. First, we extended the original model to detect insertions as well as deletions.
Second, we allowed the overdispersion parameter to vary per site. Third, as described below,
instead of using all other samples from a donor as a matched normal panel, we used samples
from different donors as controls, filtering germline mutations after variant calling
```

### Keogh2018
### Leesix2019
### Zhang2019
### Yokoyama2019
### Brunner2019
### Lawson2020
### Fowler2021
### Yamaguchi2022

