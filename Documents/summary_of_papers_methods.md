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

Brain samples from Parkinson's disease (PD) or Lewy body (LB) (n=20), Alzheimer's disease (AD, n=20), age matched healthy controls (n=14).
Brain regions were: cerebellum (CB, n=54), Entorhinal cortex (EC, n=53), frontalcortex (FC, n=32), medulla (Med, n=24), Cingulate (Cin, n=10). Blood (n=6). See Sup. Tab. 4.

We sequenced all exons of 56 genes known to cause, or predispose to, common neurodegenerative disorders (132,617 base pairs (bp)) (Supple-
mentary Table 1, left), and 46 control genes expressed at low levels in the brain which are typically associated with cancer
(152,519 bp) (Supplementary Table 1, right, subsequently referred to as ‘cancer’ genes)

Bioinformatics: trim galore, bwa mem, mark duplicates, GATK indel realigner, GATK base recalibrator.
Germline: GATK haplotypecaller, varscan2 single sample VAR>20%
Somatic: mutect2 and varscan2, deepSNV

### Leesix2019
### Zhang2019
### Yokoyama2019
### Brunner2019

### Lawson2020
2097 bladder microbiopsies from 20 individuals using targeted (n = 1914 microbiopsies), whole-exome (n = 655), and whole-genome (n = 88) sequencing.

we studied 1647 microbiopsies from 15 deceased transplant organ donors and 450 microbiopsies from five patients with bladder cancer
(table S1)

targeted sequencing of 321 cancer-associated genes for 1914 microbiopsies (median coverage of 89×)

we performed whole-exome sequencing of 655 microbiopsies (median coverage of 72×)

and whole-genome resequencing of 88 microbiopsies dominated by large clones (median cover-age of 33×)

### Fowler2021

```
Normal skin from chronically and intermittently sunexposed sites was collected from 35 Caucasian donors whose
ages ranged from 26 to 79 with a balance of males and
females (Fig.  1B; Supplementary Table S1)
```

```
Normal skin samples were collected from patients undergoing wide
local excision after initial melanoma excision, patients undergoing
browplexy, or deceased organ donors from whom organs were being
retrieved for transplantation
```

```
A total of 1,261 2-mm2
pieces were sequenced at an average coverage of 690× using a
bait set of 74 cancer-associated genes (Fig. 1C; Supplementary
Tables S2 and S3; ref. 1). Mutations were called using the
ShearwaterML algorithm, which detects mutations present in
1% or less of nucleated cells in the sample (5, 6). In a total area
of 25.2 cm2 of skin sampled across all donors, we identified
47,977 single-base substitutions (SBS), 3,824 double-base substitutions (DBS),
and 2,090 small (<200 bp) insertion or deletion events (indels) after merging
mutations shared between
adjacent samples (Fig. 1D; Supplementary Table S4; Methods)
```

p18 main pdf
```
The two sequence capture bait sets used in this study have been
described previously ( 5, 9 ). The “grid” bait set contains a set of 74
genes recurrently mutated in SCC and BCC as well as genes commonly mutated in
other epithelial cancers. The “punches + follicles”
bait set is a broader range of genes frequently mutated in a range of cancers
based on the COSMIC cancer gene census ( https://cancer.sanger.
ac.uk/census ). Samples were sequenced with each bait set as detailed
below using fat/dermis from the same patient as a germline control. A
list of all genes covered by the bait sets can be found in Supplementary
Table S3, and metrics of the bait sets are summarized below.
```

```
BAM files were mapped to the GRCh37d5 reference genome
using BWA-mem (version 0.7.17; ref. 21 ) and targeted sequencing
was aligned using the GATK tool IndelRealigner (version 3.6.0;
ref. 22 ). Duplicate reads were marked using Biobambam2 (Biobambam2 version
2.0.86. https://gitlab.com/german.tischler/biobambam2 ,
https://www.sanger.ac.uk/science/tools/biobambam ). Depth
of coverage was calculated using Samtools (version 0.1.18) to exclude
reads which were unmapped, not in the primary alignment, failing
platform/vendor quality checks, or were PCR/Optical duplicates.
BEDTools (version 2.23.0) coverage program was then used to calculate
the depth of coverage per base across samples (Supplementary
Table S2).
 To determine which samples were suitable for WGS, we used the
VAF from targeted sequencing to determine which samples were
clonal. For the 0.25-mm samples, WGS was performed on all clonal
samples in the individuals studied. For the 2-mm 2 samples, only eight
out of 1,261 samples were clonal, and all of these were sequenced.

For targeted sequencing data, subclonal mutation variant calling
was made using the deepSNV R package (also commonly referred to
as ShearwaterML), version 1.21.3, available at https://github.com/
gerstung-lab/deepSNV , used in conjunction with R version 3.3.0
(2016-05-03; ref. 5 ).
 deepSNV makes use of statistical testing to differentiate sequencing
 errors from true low-frequency mutations and has been shown
to be reliable down to a detection limit of 1/10,000 alleles ( 23 ). The
statistical tests compare by position and strand between skin samples
and a panel of control samples to estimate how likely an observed
nucleotide is a sequencing error or a true variant. Combining the
information for each strand generates a single value used for fi ltering
false-positive variants. It was noted in development that the performance
of deepSNV was not strongly dependent upon P, q -values,
or PCR amplifi cations, and its sensitivity can be increased through
higher sequencing depths. A q -value of 0.01 was used to fi lter the
variant calls.
 Fat and dermis was used as the germline sample for each donor,
and sequenced as outlined previously. Aligned germline BAM fi les
for each corresponding sample type (grids, punches, and follicles)
were provided to deepSNV, excluding the germline for the sample
being analyzed, to form a normal sample panel used for statistical
testing and false-positive variant identifi cation. Variants called from
a donor’s germline sample were subtracted from the list of variants
called from nongermline samples belonging to the donor.
 Mutations called by ShearwaterML were fi ltered using the following criteria:
● Positions of called SNVs must have a coverage of at least 100 reads
(10 reads for 0.25 mm diameter punch biopsy samples).
● Germline variants called from the same individual were removed
from the list of called variants.
● The P values of the putative mutations were adjusted with FDR
and fi ltered with a q -value threshold of 0.01.
● Mutations not present in at least one read from both strands were
removed.
● Pairs of SNVs on adjacent nucleotides within the same sample are
merged into a dinucleotide variant if at least 90% of the mapped
DNA reads containing at least one of the SNV pair contained both
SNVs.
● Identical mutations found in multiple contiguous tissue biopsies
are merged and considered as a single clone in order to prevent
duplicate clone counting. 

For the hair follicle samples, due to the very low input quantity
of DNA, particularly toward the middle and base of the follicle,
false-positive variant calls attributed to sequencing artifacts are
more likely. We therefore applied a strict method of variant calling
with ShearwaterML in order to be as conservative as possible. This
involved, for all samples within each follicle, calling variants against a

custom normal panel consisting of all other follicle samples from all
donors, plus the dermis/fat samples from all other donors. The corresponding dermis/fat sample for the donor of that follicle was then
used to remove germline variants. For follicle-spanning mutations,
only follicles where adjacent segments were successfully sequenced
or where the same mutation is present within the same follicle, i.e.,
base and top, are shown.
Shearwater was run with a normal panel of >24,000, >31,000, and
>12,000× mean coverage depth for the 2-mm2 grid samples, hair follicles, and punch biopsy samples, respectively.
For WGS, data variants were called using the CaVEMan and
Pindel algorithms (24, 25). For SNVs, CaVEMan was run with the
major copy number set to 10 and the minor copy number set to 2.
Only SNVs that passed all CaVEMan filters were kept. Additional
filtering to remove mapping artifacts associated with BWA-MEM
were: the median alignment score of reads supporting a variant had
to be at least 140 and the number of clipped reads equal to zero. In
addition, the proportion of mutant reads present in the matched
sample also had to be zero. Variants with at least one mutant read
present in the matched sample were also removed. Two SNVs called
at adjacent positions within the same sample were merged to form a
doublet-base substitution if at least 90% of the mapped DNA reads
containing at least one of the SNV pair contained both SNVs. Small
(<200 bp) insertions and deletions were called using Pindel. Only
indels that passed all Pindel filters were kept. For the punch samples
only, variants were filtered to remove a large excess of single base
pair insertions at homopolymers of length five or more, an artifact
likely caused by PCR amplification of low-input DNA concentrations during WGS. Indels were then classed as clonal if the VAF was
at least 0.3.
Variants were annotated using VAGrENT (26). Full lists of called
variants from 2-mm2 gridded samples, 0.25 mm diameter punch
samples, and hair follicles are shown in Supplementary Tables S4,
S7, and S8, respectively
```

```
1. Martincorena I, Roshan A, Gerstung M, Ellis P, Van Loo P, McLaren S,
et al. Tumor evolution. High burden and pervasive positive selection
of somatic mutations in normal human skin. Science 2015;348:880–6
5. Martincorena I, Fowler JC, Wabik A, Lawson ARJ, Abascal F, Hall
MWJ, et  al. Somatic mutant clones colonize the human esophagus
with age. Science 2018;362:911–7.
6. Gerstung M, Papaemmanuil E, Campbell PJ. Subclonal variant calling
with multiple samples and prior knowledge. Bioinformatics 2014;30:
1198–204
```

### Yamaguchi2022
Note: no access to this data without doing another access request.
```
In this work, we perform target-gene sequencing, whole-exome
sequencing (WES), and whole-genome sequencing (WGS) for 1311
endometrial glands from 37 women across a wide range of ages.
```

```
 We conducted target-gene sequencing of 112 genes in 891
normal endometrial glands from the uteri of 32 subjects ranging
in age from 21 to 53 years old (Fig. 1a, Supplementary Table 1
and Supplementary Data)
```

```
Code for statistical analyses on mutation burden, mutational signatures, and timing of
genome events including clonal expansions and copy neutral loss-of-heterozygosity is
deposited on GitHub at https://github.com/HirofumiNakaoka/endometrium_natcommun_2021.
```

```
The target-gene sequencing of single endometrial glands
for 112 genes was performed as described in our previous studies with some
modifications15,59,62,63. Briefly, 112 genes were selected (Supplementary Data)
based on WES data for ovarian endometriosis and normal uterine endometrium15,
the mutation profiles in endometriosis-related ovarian cancer64 and in endometrial
cancer48, and genes involved in DNA repair pathways65,66.

```
