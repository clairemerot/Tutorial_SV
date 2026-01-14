# Tutorial_SV
This is half-a-day training about SV detection using diverse sets of data. 
We will use a dataset provided by my collaborators and myself (C. Mérot) including short-reads, long-reads and assemblies of the seaweed fly _Coelopa frigida_. This species is interesting for its large polymorphic inversions.
To ensure efficient computational time, we will work on a dummy genome. I also made dummy labels for the samples to simplify the dataset. 

## 01 - Using SNPs & local PCA to detect haploblocks putatively representing large rearrangements
# Data and pre-processing
In preparation for this analysis I have selected 22 samples (called A to V) sequenced with illumina short-reads 150 paired-end. 
On those samples we would like to call genetic variants, starting with SNPs. The distribution of SNP variation and linkage disequilibrium can be very informative to detect putative non-recombining haploblocks, which may be large structural rearrangements like inversions.

To save energy (running it only once) and time), I pre-processed the data for you. I used (fastp https://github.com/OpenGene/fastp) to trim for adaptor and assess quality. Then I used bwa-mem (https://arxiv.org/abs/1303.3997) to align the sequences on the genome, and samtools to sort and index the resulting bam files. I also ran a deduplication module from Picardtools to remove PCR-duplicates (https://github.com/broadinstitute/picard/tree/master/src/main/java/picard/sam/markduplicates). Then I used bcftools mpile up to call SNPs into a vcf file (https://samtools.github.io/bcftools/howtos/variant-calling.html). 

#Local PCAs with lostruct
First, we will use a R package called lostruct (https://github.com/petrelharp/local_pca) which perform local PCA along the genome and groups together (using a MDS - multidimensional scaling analysis) the windows that look alike each other. To fully understand the idea of the method you may want to read the paper.
Local PCA Shows How the Effect of Population Structure Differs Along the Genome, Han Li and Peter Ralph, Genetics January 1, 2019 vol. 211 no. 1 289-304.

If you want to run yourself all the initial steps to go from the vcf of SNPs to the matrix of PCA-by-windows, you can follow the detailled tutorial here [insert link]. If you are short in time, I recommend that you use the pcamatrix that I made for you here: [insert path]

Briefly, we do several steps to convert into a bcf. Then we use a function in lostruct to make windows of your chosen size. We suggest to use window of 1000 snps. Typically with whole genome you may first run by windows of 5000 snps (or more) for a first look, and then refine with smaller windows. The analysis can be run chromosome by chromosome (as in the paper) or on the entire genome. Here, we have just one scaffold. Then, lostruct run the PCA on all windows. Here we choose to consider k=npc=2 because they usually capture most variance for each local PCA

It outputs a matrix pcs in which each row give the first k eigenvalues and k eigenvectors for each window. This gives you a matrix with 47 columns (3 columns of info, 22 columns with PC1 score for each individual, and 22 column with PC2 score for each individual). It has as many rows as windows. I added 3 columns of information about the window position. 

In a terminal, you can have a look at the matrix with
```
less
```
(use "q" to exit less visualisation in a terminal)

Back to work, we will run the end of lostruct procedure. you can do it either on the terminal or in Rstudio on your computer
```
library(lostruct)
#load matrix
pca_matrix_noNA<-read.table("pca_matrix.txt", sep="\t", header=T, stringsAsFactors=FALSE)
head(pca_matrix_noNA)
#split columns with positions information and PC
window_pos_noNA<-pca_matrix_noNA[,1:3]
pcs_noNA<-as.matrix(pca_matrix_noNA[,4:dim(pca_matrix_noNA)[2]])
```
The lostruct procedure proposes to compute pairwise distances between those windows and visualise it with a MDS (multidimensional scaling). Our goal is to identify groups of windows which display similar PCA pattern.This is done with the following functions (we uses 2 PC per window as above, and will look at the 1st 10 axes of the MDS)

```
pcdist <- pc_dist(pcs_noNA,npc=2)
mds_axe<-cmdscale(pcdist, k=10)
head(mds_axe)

#again the mds file is missing position information so:
mds_matrix<-cbind(window_pos_noNA, mds_axe)
write.table(mds_matrix, "mds_matrix.txt", sep="\t", row.names=FALSE, quote=FALSE)
```
What do you think? Do we have a region where we could suspect non-recombining blocks revealing a putative inversion?
How could we do to genotype our individuals for this putative variant?


A recent tool was published which offers a better visualisation of PC scores along the genome. This is called Winpca - https://github.com/MoritzBlumer/winpca.
We will play a little with this package to explore our dataset.

In a terminal:
```
#winpca code
```


Using either information from lostruct or winpca, may you identify approximative coordinates for the suspected rearrangement?

A step further... 
It may be interesting to run population genetics classic analysis to better capture the variation within and across the putative inversion.
For example, you may be interested in calculating LD across distant windows. You may also want to calculate FST statistics between homokaryotypic groups or pi diversity within each group. If you have time you can return to this point and re-apply the methods you learnt over the last few days to the groups determined by the inversion.
Some examples of code here: [insert link]





## 02 - Comparing assemblies for large rearrangements
The ideal dataset to confirmlarge rearrangements would be fully-assembled assemblies, although those will rarely be available for non-model species.
Here, let's say that we are lucky and two different individual have been sequenced with long-reads and assembled as chromosome-scale assembly. One individual has been picked in the A/A group (ref.fasta) and one in the B/B group (alt.fasta).

We will work in the second folder 02_assemblies.

To compare those assemblies we need to align them to each other. For this we will use Minimap 2. https://github.com/lh3/minimap2 . You can check on the Minimap manual what are the different options. This will produce a file listing syntenic and aligned regions in a paf format. 
This file can be opened with a regular text editor to see what it looks like.  

```
minimap2 -x asm20 -c --eqx ~/workshop_materials/structural_variants/assemblies/ref.fasta ~/workshop_materials/structural_variants/assemblies/alt.fasta > 02_assemblies/alt-to-ref20.paf
```

Let's quickly visualise this alignment. D-genies offers a quick dotplot between two genomes. You can either upload the two fasta files or directly the paf. (Please favour the paf file for today to avoid redudant demands on the server)

What do you think? Do you visualise some rearrangements? Does it fit what you uspected from the PCA analysis?

Now we may want to do nicer plots. A recent R package will help you do nice ribbon plot. Let's try SVbyEye. https://github.com/daewoooo/SVbyEye
Here is a very quick code to get a first figure. You can try it in Rstudio. If you have time you may explore more options from this package.

```
library(SVbyEye)
paf.table <- readPaf("02_assemblies/alt-to-ref20.paf",
  include.paf.tags = TRUE, restrict.paf.tags = "cg")
paf.table
plotMiro(paf.table = paf.table, color.by = "direction")
```

Using the paf file and/or the R package, can you extract the breakpoints of the suspected rearrangement? How does it compare to the putative breakpoints identified by indirect PCA methods based on recombination suppression?

## 03 - Assessing breakpoints and detecting other SVs with long-reads
We have an additionnal sample (fri59) which has been sequenced with Oxford nanopore technology.
I did basic filtering steps and aligned it, using minimap2, on our dummy genome. The resulting bamfile is available here ~/workshop_materials/structural_variants/ONT/indA59.bam
Let's use it to dig into the breakpoints of our large rearrangement, but also to look for small structural variants all along the genome.

To do so, we will rely on a widely-used software, Sniffles 2 https://github.com/fritzsedlazeck/Sniffles . 
Be aware that several SV callers are available and, depending on your data, you may favour one over the others or decide on using several tools.
```
sniffles -i ~/workshop_materials/structural_variants/ONT/indA59.bam --reference ~/workshop_materials/structural_variants/assemblies/ref.fasta  --vcf 03_LR/fri59.vcf –minsupport 10
```

Next, we will try to visualize our SVs to assess how confident we can be in the detection. We suggest to use IGV which may be heavy to use on a large chromosome but will be ok on our dummy genome.
Please go on the desktop and start IGV.
You can upload the genome, then you can upload your vcf, and then you bam file. You can try playing with the zoom to visualise the SVs and the read support.

As you may notice in the vcf track we have a lot of long SVs which are most likely wrong. Some have the label "PASS" in the vcf while others have the label "GT" and a genotype 0/0 meaning they have much less upport and are likely wrong.
you can hide the ones that do not pass filters (in light blue) in the interface of IGV. Alternatvely you can filter them on the server with bcftools.

```
bcftools view --apply-filters PASS fri59.vcf > fri59_pass.vcf
```

Please take some time to explore the visualisation. For example, you can try finding a deletion with good support and a deletion which you would be less confident in (same for DUP or INS). You can click on each variant in the vcf panel to get more information. In the lower panel you can visualize coverage and then each raw read. 
When you clik on each read you can get the information about where it maps, whether it is split and map elsewhere, etc. 

You can also zoom on the breakpoints previously identified with the assembly comparison for the very large rearrangement. What do they look like here? Are they supported by long-reads? Which genotype is this individual fri59? 



## 04 - Using population dataset of short-reads to detect SVs

## 05 - Towards graph-based analysis
