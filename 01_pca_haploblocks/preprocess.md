## Step 1 Sliding PCA : skipped step of preparing files

####  prepare the vcf (in the terminal) [do not run]
First you need to  convert the vcf into a bcf, sort it and index it.
The package lostruct can also take a vcf but it is heavy to load the whole file. Instead they made function that is able to cut the vcf window by window to avoid overloading the memory (but this requires a sorted, indexed bcf)
```
#we need to sort the vcf
(grep ^"#"  ~/workshop_materials/structural_variants/SNPs/SNPs.vcf; grep -v ^"#" SNPs.vcf | sort -k1,1 -k2,2n) > 01_pca_haploblocks/SNPs_sorted.vcf

#then compress it and index it
bgzip -c 01_pca_haploblocks/SNPs_sorted.vcf > 01_pca_haploblocks/SNPs_sorted.vcf.gz
tabix -fp vcf 01_pca_haploblocks/SNPs_sorted.vcf.gz

#then convert it to a bcf and sort it
bcftools convert -O b 01_pca_haploblocks/SNPs_sorted.vcf.gz > 01_pca_haploblocks/SNPs_sorted.bcf
bcftools index 01_pca_haploblocks/SNPs_sorted.bcf
```
To check whether this has worked and produced files, you can do a quick ```ls -lh 01_pca_haploblocks``` which will show you the size of the files in human-readable format

Now we are good to work in R with the library. This requires some computational power and memory, so we suggest to make the initial steps in R command lines on the server and then copy the output files to visualise on your local Rstudio

####  run lostruct (in the terminal)
To start R in command line, just type "R". Now you have a R console and we wil run the lostruct procedure
```
#open library
library(lostruct)
snps <- vcf_windower("01_pca_haploblocks/SNPs_sorted.bcf",size=1000,type='snp', sites= vcf_positions("01_pca_haploblocks/SNPs_sorted.bcf"))
```
This function makes windows out of the given data file of your chosen size. You can choose the size of the window with "size" and on which variable you want to split ('snp' or 'bp'). We suggest to use window of 100 snp since we are not very dense (RAD data) and we don't have a lot of snps.Typically with whole genome you may first run by windows of 1000 or 5000 snps for a first look, and then refine with smaller windows. The analysis can be run chromosome by chromosome (as in the paper) or on the entire genome. Here, we are going for the entire genome.
You can display for instance the 5th window and know its location by doing
```
snps(5)
region(snps) (5)
```
We can now run the PCA on all windows. Here we choose to consider k=npc=2 because they usually capture most variance for each local PCA
```
pcs <- eigen_windows(snps,k=2)
dim(pcs) #check dimension
head (pcs[,1:10]) #look at the first 10 columns
```
In the matrix pcs, each rows give the first k eigenvalues and k eigenvectors for each window. This gives you a matrix with 47 columns (3 columns of info, 22 columns with PC1 score for each individual, and 22 column with PC2 score for each individual). It has as many rows as windows (463 with windows of 1000 SNPs)

As you see we don't know the position of each window, we will get it with the function regions, remove the NA windows and expand the pca matrix to include the position information we retrieve before and export the file
```
#retrieve positions
window_pos<-region(snps)()
head(window_pos)
#merge
pca_matrix<-cbind(window_pos, pcs)
head (pca_matrix[,1:10])
#save the file
write.table(pca_matrix, "01_pca_haploblocks/pca_matrix.txt", sep="\t", row.names=FALSE, quote=FALSE)

#keep windows without NA: because of msising data in the vcf, some windows may not computed by pca (particularly true when doing small window with 100SNPs).
In our example we don't need it but here is how to remove them. 
pcs_noNA<-pcs[-which(is.na(pcs[,1])),]
dim(pcs_noNA) # if the dimension is > 0, one should remove the NA -lines # in our case we don't need
window_pos_noNA<- window_pos[-which(is.na(pcs[,1])),]
pca_matrix_noNA<-cbind(window_pos_noNA, pcs_noNA)
head (pca_matrix_noNA[,1:10])
#save the file
write.table(pca_matrix_noNA, "01_pca_haploblocks/pca_matrix_noNA.txt", sep="\t", row.names=FALSE, quote=FALSE)
```
