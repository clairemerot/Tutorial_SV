## 01 - Using SNPs & local PCA to detect haploblocks putatively representing large rearrangements

### Data and pre-processing
In preparation for this analysis I have selected 22 samples (called A to V) sequenced with illumina short-reads 150 paired-end. 
On those samples we would like to call genetic variants, starting with SNPs. The distribution of SNP variation and linkage disequilibrium can be very informative to detect putative non-recombining haploblocks, which may be large structural rearrangements like inversions.

To save energy (running it only once) and time), I pre-processed the data for you. I used (fastp https://github.com/OpenGene/fastp) to trim for adaptor and assess quality. Then I used bwa-mem (https://arxiv.org/abs/1303.3997) to align the sequences on the genome, and samtools to sort and index the resulting bam files. I also ran a deduplication module from Picardtools to remove PCR-duplicates (https://github.com/broadinstitute/picard/tree/master/src/main/java/picard/sam/markduplicates). Then I used bcftools mpile up to call SNPs into a vcf file (https://samtools.github.io/bcftools/howtos/variant-calling.html). 

### Local PCAs with lostruct
#### pre-process
First, we will use a R package called lostruct (https://github.com/petrelharp/local_pca) which perform local PCA along the genome and groups together (using a MDS - multidimensional scaling analysis) the windows that look alike each other. To fully understand the idea of the method you may want to read the paper.

Local PCA Shows How the Effect of Population Structure Differs Along the Genome, Han Li and Peter Ralph, Genetics January 1, 2019 vol. 211 no. 1 289-304.https://doi.org/10.1534/genetics.118.301747

If you want to run yourself all the initial steps to go from the vcf of SNPs to the matrix of PCA-by-windows, you can follow the detailled tutorial [here](https://github.com/clairemerot/Tutorial_SV/blob/main/01_pca_haploblocks/preprocess.md). If you are short in time, I recommend that you use the pcamatrix that I made for you. In a terminal, you can have a look at the matrix with
```
less ~/workshop_materials/structural_variants/SNPs/pca_matrix.txt
```
(use "q" to exit less visualisation in a terminal)

Briefly, we do several steps to convert into a bcf. Then we use a function in lostruct to make windows of your chosen size. We suggest to use window of 1000 snps. Typically with whole genome you may first run by windows of 5000 snps (or more) for a first look, and then refine with smaller windows. The analysis can be run chromosome by chromosome (as in the paper) or on the entire genome. Here, we have just one scaffold. Then, lostruct run the PCA on all windows. Here we choose to consider k=npc=2 because they usually capture most variance for each local PCA

It outputs a matrix pcs in which each row give the first k eigenvalues and k eigenvectors for each window. This gives you a matrix with 47 columns (3 columns of info, 22 columns with PC1 score for each individual, and 22 column with PC2 score for each individual). It has as many rows as windows. I added 3 columns of information about the window position. 

#### analysis
You will run the end of lostruct procedure. you can do it either on the terminal or in Rstudio on your computer

If you have not run the preparatory steps, please copy the pca_matrix already prepared to your own directory.

```
cp ~/workshop_materials/structural_variants/SNPs/pca_matrix.txt 01_pca_haploblocks/pca_matrix.txt
```

```
library(lostruct)
#load matrix
pca_matrix<-read.table("01_pca_haploblocks/pca_matrix.txt", sep="\t", header=T, stringsAsFactors=FALSE)
head(pca_matrix)
#split columns with positions information and PC
window_pos<-pca_matrix[,1:3]
pcs<-as.matrix(pca_matrix[,4:dim(pca_matrix)[2]])
```
The lostruct procedure proposes to compute pairwise distances between those windows and visualise it with a MDS (multidimensional scaling). Our goal is to identify groups of windows which display similar PCA pattern.This is done with the following functions (we uses 2 PC per window as above, and will look at the 1st 10 axes of the MDS)

```
pcdist <- pc_dist(pcs,npc=2)
mds_axe<-cmdscale(pcdist, k=10)
head(mds_axe)

#again the mds file is missing position information so:
mds_matrix<-cbind(window_pos, mds_axe)
write.table(mds_matrix, "01_pca_haploblocks/mds_matrix.txt", sep="\t", row.names=FALSE, quote=FALSE)
```

####  Using the MDS
We are looking for regions with similar clustering of individuals that may reveal population structure, chromosomal rearragements, sex, non recombining haploblocks, etc. Exploring the MDS is a way to detect such heterogeneity in the genome.

Let's load the mds and plot the first axes
```
library(ggplot2)
head(mds_matrix)
#for ggplot you need to rename the columns
colnames(mds_matrix)<-c("chrom","start","end","mds1","mds2","mds3","mds4","mds5","mds6","mds7","mds8","mds9","mds10")
#we will also add a column with the mid position of each window
mds_matrix$midpos<-(mds_matrix$start+mds_matrix$end)/2

ggplot(mds_matrix, aes(x=mds1, y=mds2, colour=chrom))+
  geom_point()+
  theme_classic()
```
A group of windows seems to cluster at high values on mds1. Let's see whether those windows are also in the same area of our genome by plotting the mds1 scores along the chromosome:

```
ggplot(mds_matrix, aes(x=midpos, y=mds1, colour=chrom))+
  geom_point()+
  theme_classic()
```
What do you think? Do we have a region where we could suspect non-recombining blocks revealing a putative inversion?
How could we do to genotype our individuals for this putative variant?

NB: Here we live in a simplified world... Most of the time, if you run the lostruct analysis on all the genome or on each chromosome it may be worth exploring mds2, mds3, etc. It may also be relevant to decide a threshold on mds score (for exemple 3 times the standard deviation) beyond which windows will be considered outliers and will then receive further attention.

#### Exploring PCAs themselves
We may want to simply plot the pca for some windows. This is not the most easy because remember the format is a bit tricky
Look at the format. We have 3 columns for position, total eigen values, eigvalue of PC1, of PC2 and then 22 values for PC1 scores of all our samples, and 22 values for PC2 scores of all samples
```
  chrom  start    end     total     lam_1     lam_2      PC_1_A      PC_1_B      PC_1_C
1  Chr1     35  19975 0.7825983 0.3847001 0.3071047  0.27546442  0.23256006  0.01978089
2  Chr1  20028  39517 0.8102097 0.3601706 0.3408839 -0.04928240  0.01148246 -0.07928459
3  Chr1  39537  61771 0.9387336 0.4584682 0.4118569  0.11736869 -0.17638771  0.03902070
4  Chr1  61772  86432 0.8102367 0.5066009 0.2716244 -0.09503886 -0.14427806 -0.03715027
5  Chr1  86452 108171 0.7496138 0.4306277 0.2852556  0.07593774 -0.20193475 -0.11323903
6  Chr1 108175 122831 0.7926275 0.3967518 0.3167858 -0.26975570 -0.43745346 -0.14152673
```

So to get information for a window we can do something like:
```
Nind<-22
i=400 #for the 400th window

pc1_i<-t(pca_matrix[i, 7:(Nind+6)]) #scores along PC1
pc2_i<-t(pca_matrix[i, (Nind+7):(2*Nind+6)]) #scores along PC2
var1<-round(pca_matrix[i, 5]/pca_matrix[i, 4],2)*100 # % of variance explained by PC1
var2<-round(pca_matrix[i, 6]/pca_matrix[i, 4],2)*100 # % of variance explained by PC2
midpos_i<-(pca_matrix[i, 2]+pca_matrix[i, 3])/2 #average position of the window
window_i<-paste(pca_matrix[i, 1], midpos_i , sep="_") #paste the name of CHR and the midposition

plot(pc1_i, pc2_i, pch=20, xlab=paste("PC1", var1 , "%"), ylab=paste("PC2", var2, "%"), main=window_i)
```
You may try to visualise a PCA within and outside your region of interest. What do you think?

In the folder 00_ressources you will find a file classifying the samples as AA, AB and BB (info_ind.txt). Those could be imagined as populations or groups with different versions of a rearrangements. I extracted this information from a PCA in our region of interest and confirmed it with a PCR marker. Does it correspond to the groups you observe?

### WinPCA - visualising PC1 along the genome.

A recent tool was published which offers a better visualisation of PC scores along the genome. This is called Winpca - https://github.com/MoritzBlumer/winpca.
We will play a little with this package to explore our dataset. The principle is to perform pca on windows along the genome (here defined by a fixed size (while before we defined them based on the number of SNPs). Then, the 1st of the 2nd PCs are plotted along the chromosome. Each line is a sample and this line connects its coordinate on PC1 of all the windows. Of course, PCs can also be flipped to improve visualisation.

Here we pick windows of 10000bp incremented by 10000bp. We will use genotypes (GT) for biallelic SNPs encoded in the vcf. Note that if you are processing low-depth data with angsd winpca can also use GL and pcangsd.

To run winpca in the terminal
```
#filter for biallelic SNPs - otherwise we run into an error:
bcftools view -m2 -M2 -v snps ~/workshop_materials/structural_variants/SNPs/SNPs.vcf > 01_pca_haploblocks/SNPs_biallelic.vcf

#lets load an adequate conda env
conda activate winpca

#Then let's run winpca code

#make the PCAs
#-w for window size -i for increment size --np to remove filters creating an error -v GT to precise the type of data.
#then there are three argument "$PREFIX" "$VCF" "$REGION"
winpca pca -w 10000 -i 10000 --np -v GT 01_pca_haploblocks/winpca_out 01_pca_haploblocks/SNPs_biallelic.vcf Chr1:1-9999999

#polarize them
winpca polarize 01_pca_haploblocks/winpca_out

#plot PCs
winpca chromplot 01_pca_haploblocks/winpca_out Chr1:1-9999999

#get out of the env
conda deactivate
```
Winpca outputs a .html file in which you can directly visualise PC1 along the chromosome (below) and the proportion of variance capture by pc1.

You may want to do your own plots. To do so, we can unzip the .tsv.gz files formed using the command 
```gunzip 01_pca_haploblocks/*.tsv.gz```

Then we can go in Rstudio for some reformatting and plotting. We will use libraries which help making "tidy data" such as tidyr and dplyr, as well as ggplot2 as before for plotting.
We will re-color the plot based on haplogroups determined in our previous PCAs.

```
library(ggplot2)
library(tidyr)
library(dplyr)
#read pc1 scores
winpca_scores<- read.delim("01_pca_haploblocks/winpca_out.pc_1.tsv")
head(winpca_scores)

#read stats
stat_pca<-read.delim("01_pca_haploblocks/winpca_out.stat.tsv")
head(stat_pca)

#read information about the samples 
info<-read.delim("00_ressources/info_ind.txt")
head(info)

#let's reformat the data as tidy data and add some information about the samples
winpca_scores_long<-pivot_longer(winpca_scores, cols=A:V)
head(winpca_scores_long)
winpca_scores_long_info<-left_join(winpca_scores_long,info)

head(winpca_scores_long_info)
ggplot(data=winpca_scores_long_info, mapping=aes(x=pos, y=value, group=name, col=factor(haplogroup)))+ geom_line()+ theme_classic()

```

Using either information from lostruct or winpca, may you identify approximative coordinates for the suspected rearrangement?


### A step further... 
It may be interesting to run population genetics classic analysis to better capture the variation within and across the putative inversion.

It may be relevant to calculate LD across distant windows. You may want to calculate FST statistics between homokaryotypic groups or pi diversity within each group. If you have time you can return to this point and re-apply the methods you learnt over the last few days using the vcf and treating as populations the haplogroups determined by the inversion (either by yourself based on the pc scores or kmeans approaches, or using the info_ind.txt file provided in 00_resssources).

You can find this approach of localPCAs used in several recent papers. Here are a few that I know:
Huang K, Andrew RL, Owens GL, Ostevik KL, Rieseberg LH. Multiple chromosomal inversions contribute to adaptive divergence of a dune sunflower ecotype. Mol Ecol. 2020; 29: 2535–2549. https://doi.org/10.1111/mec.15428
Mérot, C., Berdan, E. L., Cayuela, H., Djambazian, H., Ferchaud, A. L., Laporte, M., ... & Bernatchez, L. (2021). Locally adaptive inversions modulate genetic variation at different geographic scales in a seaweed fly. Molecular biology and evolution, 38(9), 3953-3971.

To get an overview of the kind of approach you may also visualise this summary [here](https://github.com/clairemerot/Tutorial_SV/blob/main/00_ressources/images/merot_huang.jpg)


