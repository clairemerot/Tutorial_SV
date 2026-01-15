# Tutorial_SV
This is half-a-day training about SV detection using diverse sets of data. 
We will use a dataset provided by my collaborators and myself (C. Mérot) including short-reads, long-reads and assemblies of the seaweed fly _Coelopa frigida_. This species is interesting for its large polymorphic inversions.
To ensure efficient computational time, we will work on a dummy genome. I also made dummy labels for the samples to simplify the dataset. 

## 01 - Using SNPs & local PCA to detect haploblocks putatively representing large rearrangements
Using 22 flies sequenced with short-reads, I called SNPs (stored in a vcf file). Using this dataset we will perform local PCAs along the genome (cutting the genome in small window on which we decompose genetic variation into principal compenents). We are looking for regions with unexpected population structure. Typically, non-recombining haploblocks, like chromosomal inversions or other regions, will display three clusters (homozygote for one arrangement/one haplotypic block, heterozygotes and homozygote for the other arrangement). If we find such a region (spoiler alert... yes!), we will explore it more deeply with this dataset and others. 

Methods comes from those papers:
Local PCA Shows How the Effect of Population Structure Differs Along the Genome, Han Li and Peter Ralph, Genetics January 1, 2019 vol. 211 no. 1 289-304. https://doi.org/10.1534/genetics.118.301747
L Moritz Blumer, Jeffrey M Good, Richard Durbin, WinPCA: a package for windowed principal component analysis, Bioinformatics, Volume 41, Issue 10, October 2025, btaf529, https://doi.org/10.1093/bioinformatics/btaf529

They have been used in empircial studies such as:
Huang K, Andrew RL, Owens GL, Ostevik KL, Rieseberg LH. Multiple chromosomal inversions contribute to adaptive divergence of a dune sunflower ecotype. Mol Ecol. 2020; 29: 2535–2549. https://doi.org/10.1111/mec.15428
Mérot, C., Berdan, E. L., Cayuela, H., Djambazian, H., Ferchaud, A. L., Laporte, M., ... & Bernatchez, L. (2021). Locally adaptive inversions modulate genetic variation at different geographic scales in a seaweed fly. Molecular biology and evolution, 38(9), 3953-3971.


## 02 - Comparing assemblies for large rearrangements

## 03 - Assessing breakpoints and detecting other SVs with long-reads

## 04 - Using population dataset of short-reads to detect SVs

## 05 - Towards graph-based analysis
