# Tutorial_SV
This is half-a-day training about SV detection using diverse sets of data. 
We will use a dataset provided by my collaborators and myself (C. Mérot) including short-reads, long-reads and assemblies of the seaweed fly _Coelopa frigida_. This species is interesting for its large polymorphic inversions.
To ensure efficient computational time, we will work on a dummy genome  of 10Mb. I also made dummy labels for the samples to simplify the dataset. 

To make our lives easier, I suggest that you clone this whole repository in your working directory so taht we all use the same file/folder architecture. 
Very important ** all commands are written to be run from the Tutorial_SV folder (and not from each subfolder)**. 

```git clone https://github.com/clairemerot/Tutorial_SV ```

## 01 - Using SNPs & local PCA to detect haploblocks putatively representing large rearrangements
Using 22 flies sequenced with short-reads, I called SNPs (stored in a vcf file). Using this dataset we will perform local PCAs along the genome (cutting the genome in small window on which we decompose genetic variation into principal compenents). We are looking for regions with unexpected population structure. Typically, non-recombining haploblocks, like chromosomal inversions or other regions, will display three clusters (homozygote for one arrangement/one haplotypic block, heterozygotes and homozygote for the other arrangement). If we find such a region (spoiler alert... yes!), we will explore it more deeply with this dataset and others. 

Please follow the tutorial [here](https://github.com/clairemerot/Tutorial_SV/tree/main/01_pca_haploblocks/README.md).

Methods comes from those papers:

Local PCA Shows How the Effect of Population Structure Differs Along the Genome, Han Li and Peter Ralph, Genetics January 1, 2019 vol. 211 no. 1 289-304. https://doi.org/10.1534/genetics.118.301747

L Moritz Blumer, Jeffrey M Good, Richard Durbin, WinPCA: a package for windowed principal component analysis, Bioinformatics, Volume 41, Issue 10, October 2025, btaf529, https://doi.org/10.1093/bioinformatics/btaf529

## 02 - Comparing assemblies for large rearrangements
Next, we have a second fully-assembled genome. Obviously, this makes a perfect dataset to detect all types of variants, including very long ones (but beware of assembly errors!). We will use this comparison to confirm the putative rearrangement that we suspected during the localPCA analysis.
We will do easy plotting of genome synteny with dotplots using a stan web application called D-Genies and R package making ribbon plots (SVbyEye).

Here is the tutorial [here](https://github.com/clairemerot/Tutorial_SV/blob/main/02_assemblies/README.md)

Methods comes from those papers:

Cabanettes F, Klopp C. 2018. D-GENIES: dot plot large genomes in an interactive, efficient and simple way. PeerJ 6:e4958 https://doi.org/10.7717/peerj.4958

David Porubsky, Xavi Guitart, DongAhn Yoo, Philip C Dishuck, William T Harvey, Evan E Eichler, SVbyEye: a visual tool to characterize structural variation among whole-genome assemblies, Bioinformatics, Volume 41, Issue 6, June 2025, btaf332, https://doi.org/10.1093/bioinformatics/btaf332

## 03 - Assessing breakpoints and detecting other SVs with long-reads

tutorial [here](https://github.com/clairemerot/Tutorial_SV/blob/main/03_LR/README.md)
## 04 - Using population dataset of short-reads to detect SVs

tutorial [here](https://github.com/clairemerot/Tutorial_SV/blob/main/04_SR/README.md)
## 05 - Towards graph-based analysis

tutorial [here](https://github.com/clairemerot/Tutorial_SV/blob/main/05_graphs/README.md)
