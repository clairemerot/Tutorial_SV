# Tutorial_SV
This is half-a-day training about SV detection using diverse sets of data. 
We will use a dataset provided by my collaborators and myself (C. Mérot) including short-reads, long-reads and assemblies of the seaweed fly _Coelopa frigida_. This species is interesting for its large polymorphic inversions.
To ensure efficient computational time, we will work on a dummy genome  of 10Mb. I also made dummy labels for the samples to simplify the dataset. 

To make our lives easier, I suggest that you clone this whole repository in your working directory:

```git clone https://github.com/clairemerot/Tutorial_SV ```

We will all use the same file/folder architecture. Please open this new folder:

```cd Tutorial SV```

Now: VERY IMPORTANT ** all commands are written to be run from the Tutorial_SV folder (and not from each subfolder)**. 


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
Long-read sequencing is now provided very relevant data to detect all kinds of SVs because sequence length (10-100kb) can fully cover most of the structural variants, or pinpoint the breakpoints of larger rearrangements. For this section of the tutorial we have the bam file (reads aligned) from individual fri59 which I already aligned on the genome. We will call SVs using Sniffles and visualise reads aligned on the genome & the SVs with IGV (https://igv.org/).

Sedlazeck, F.J., Rescheneder, P., Smolka, M. et al. Accurate detection of complex structural variations using single-molecule sequencing. Nat Methods 15, 461–468 (2018). https://doi.org/10.1038/s41592-018-0001-7

The tutorial is available [here](https://github.com/clairemerot/Tutorial_SV/blob/main/03_LR/README.md)

## 04 - Using population dataset of short-reads to detect SVs
Long-reads are amazing but hardly applicable at population-scale (for now). Let's try to valorize our short-read sequences to also call SVs on a higher number of samples. In the interest of time, I have selected 4 samples (A,B,C,D) and provided you the alignement on our dummy genome (.bam). We will use Delly, a SV caller, to call and genotype SVs from those 4 samples. Note that, although ths software was originally published in 2012, it is maintained, updated and now also include long-read sequences.

Tobias Rausch, Thomas Zichner, Andreas Schlattl, Adrian M. Stütz, Vladimir Benes, Jan O. Korbel, DELLY: structural variant discovery by integrated paired-end and split-read analysis, Bioinformatics, Volume 28, Issue 18, September 2012, Pages i333–i339, https://doi.org/10.1093/bioinformatics/bts378

The tutorial is [here](https://github.com/clairemerot/Tutorial_SV/blob/main/04_SR/README.md)

## 05 - Towards graph-based analysis
Graphs are another way of representing genomes, which aim to better integrate their inherent variability. When genomes are similar, there is just one path and all variants represent an alternative path. Full pangenome graphs can be constructed from several assemblies. With a list of variants, one can also build a "variation graph" based on a reference paths and many bubbles corresponding to variants.

On those graphs, it is possible to map short-reads and ask on which path they map to extract their genotype at our variants of interest! We will do that to combine our two sets of data: we build a graph including all the variants detected by sniffles based on long reads, and then we genotype a sample (eg. A) sequenced with short-reads.

The tutorial is  [here](https://github.com/clairemerot/Tutorial_SV/blob/main/05_graphs/README.md)
