## 04 - Using population dataset of short-reads to detect SVs
Although short-reads are short (usually around 150bp), they come in pairs (here they are paired end) and this provides some information to detect structural variants. 

For example, if the pairs are further away than expected given the insert size, one may suspect an insertion in this area. If the paired-end map in reverse direction, they may signals an inversion. Finally the coverage provides a lot of data about copy-number-variants such as duplications, deletions, etc... 

Based on this idea, several SV callers have emerged a few years ago. However, athough those programms do their best, they are also confounded by repeats, mapping issues of short-reads, complex areas of the genome, etc leading sometimes to high rates of false positives. Detecting SVs from short-reads is nevertheless valuable given the most efficent cost of short-reads (relatively to long-reads and assemblies), which allows to really get variants at population level (on many samples!).

Here we will use just one possible caller but it is recommended to use several callers and carefully consider the outcomes. I picked Delly as it is still maintained, can also be run on long-reads and is easy to use.

Like for section 1, I have trimmed the fastq files, mapped fastq1/fastq2 in pairs on the genome, removed PCR duplicates, and sorted and indexed the resulting bam files. Those files are available for 4 samples here: ```ls ~/workshop_materials/structural_variants/SR ```. If you are curious about what bam files look like you can display them with samtools. Here is a line of code to look at the first 100 lines with the header (option -h).

```samtools view -h ~/workshop_materials/structural_variants/SR/A.bam | head -n 100 ```

### Calling SV with Delly
Let's run Delly on one sample first (sample B). We will use the following options -t to choose the type of SVs (here all)
-o the output file -g the reference genome -q 20 to keep alignement above a minimum quality of 20. This may take up to 10 min. Feel free to take a break :-) 

```
output=04_SR/B_SVSR.bcf
ref=~/workshop_materials/structural_variants/assemblies/ref.fasta
bam=~/workshop_materials/structural_variants/SR/B.bam

delly call -t ALL \
-o $output \
-g $ref \
-q 20 $bam 

#next lets convert the bcf to vcf
bcftools view 04_SR/B_SVSR.bcf > 04_SR/B_SVSR.vcf

echo "nb of SV detected by Delly"
grep -v ^\#\# 04_SR/B_SVSR.vcf | wc -l
```

### Calling on more samples | Visualising SVs
It is also possible to do all the samples at once - but slightly longer. Below is the code for several bam files. Since it is long, you may want to have a coffee break or play with the vcf of sample B in the mean time (see next section).

### warning - slightly long to run
```
output=04_SR/all_SVSR.bcf
ref=~/workshop_materials/structural_variants/assemblies/ref.fasta

delly call -t ALL \
-o $output \
-g $ref \
-q 20 \
~/workshop_materials/structural_variants/SR/A.bam \
~/workshop_materials/structural_variants/SR/B.bam \
~/workshop_materials/structural_variants/SR/C.bam \
~/workshop_materials/structural_variants/SR/D.bam 

#-t to choose the type of SVs - here all
# -o the output file
# -g the reference genome
# -q 20 to keep alignement above a minimum quality of 20

#next lets convert the bcf to vcf
bcftools view 04_SR/all_SVSR.bcf > 04_SR/all_SVSR.vcf

echo "nb of SV detected by Delly"
grep -v ^\#\# 04_SR/all_SVSR.vcf | wc -l
```

### Visualising and validating
Open the vcf (for exemple with ```less 04_SR/B_SVSR.vcf```) in a different terminal.

-> What do you think about the information about SVs? What differs and does not differ with the vcf produced by Sniffles?

Like we did for the SVs detected from long-reads we can also open IGV in the desktop and upload the new vcf (04_SR/B_SVSR.vcf) and the B.bam (or the 4 bams). You can walk along the genome to evaluate visually how the reads seem to support (or not) some SVs.
Reference is here ~/workshop_materials/structural_variants/assemblies/ref.fasta
Bam is here ~/workshop_materials/structural_variants/SR/B.bam

-> What do you think about our large inversion? Was it detected by Delly? 

-> What do you think of the SVs? Is it readable? 

Let's filter a little our vcf. For example, we may want to focus only on deletions which pass the built-in filters. Then, we come back to the visualisation.

```
#filter for only PASS deletions
bcftools filter -i 'FILTER="PASS" && INFO/SVTYPE="DEL"' -o 04_SR/B_SVSR_DEL.vcf -Ov 04_SR/B_SVSR.vcf

#filter for length - more tricky because it is not directly encoded in the INFO field. let's play with different tools to add it

#we need to add a field for SVLEN
##step 1 export position 
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/END\n' 04_SR/B_SVSR_DEL.vcf > 04_SR/B_SVSR.vcf.info

#step 2 calculate length
#I provide here a R script but feel free to write yours
Rscript ~/workshop_materials/structural_variants/scripts/add_info_bcf.r 04_SR/B_SVSR.vcf.info
bgzip 04_SR/B_SVSR.vcf.info.annot
tabix -s1 -b2 -e2 04_SR/B_SVSR.vcf.info.annot.gz

##step3 prepare the header
echo -e '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">' > 04_SR/B_SVSR.vcf.info.annot.hdr

##step4 run bcftools annotate
#-a is the annotation file (tabix and bgzip, it needs at least CHROM and POS, -h are the header lines to add, -c are the meaning of the column in the annotation file
bcftools annotate -a 04_SR/B_SVSR.vcf.info.annot.gz -h 04_SR/B_SVSR.vcf.info.annot.hdr -c CHROM,POS,ID,INFO/END,INFO/SVLEN 04_SR/B_SVSR.vcf > 04_SR/B_SVSR_bis.vcf

#filter vcf -i (include, -O vcf format -o
bcftools filter -i'INFO/SVLEN<=1000 && INFO/SVLEN>=-1000' -o 04_SR/B_SVSR_bis_1k.vcf -Ov 04_SR/B_SVSR_bis.vcf
echo "total number of SVs < 1kb"
grep -v ^\#\# 04_SR/B_SVSR_bis_1k.vcf | wc -l

```

-> Can you find some deletions well supported by reads ? 

When the vcf for the 4 samples is ready, open also the all_SVSR.vcf. You can now notice 4 columns corresponding to genotypes and other informations about each sample. 

We can also be interested in assessing the fraction of missing data and the frequency of alternative alleles. To do this we may use a bcftools plug-in which adds the fraction of missing data. Let's also add the frequency of the alternative allele (MAF) and export it.

```
#using the plug in to add infor in the vcf
bcftools +fill-tags 04_SR/all_SVSR.vcf -- -t 'INFO/F_MISSING','INFO/MAF' > 04_SR/all_SVSR_missingMAF.vcf

#exporting
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/END\t%INFO/F_MISSING\t%INFO/MAF\n' 04_SR/all_SVSR_missingMAF.vcf > 04_SR/all_SVSR_missingMAF.info

```

Open the info file to more easily visualize our new stats.

-> What do you think of allelic frequency? And missing data?

If you have time you can try to implement new filters to better restrict your set of SVs.

### Re-genotyping
Delly can also be used to re-genotype all the samples wihht the objective to obtain the most accurate genotypes for all samples and for all SVs. 

The command is nearly identical, except that we give the vcf, with -v.

```
ref=~/workshop_materials/structural_variants/assemblies/ref.fasta

delly call -t ALL \
-o 04_SR/all_SVSR_regenotyped.bcf \
-g  $ref \
-q 20 -v 04_SR/all_SVSR.vcf \
~/workshop_materials/structural_variants/SR/A.bam \
~/workshop_materials/structural_variants/SR/B.bam \
~/workshop_materials/structural_variants/SR/C.bam \
~/workshop_materials/structural_variants/SR/D.bam 

#next lets convert the bcf to vcf
bcftools view 04_SR/all_SVSR_regenotyped.bcf > 04_SR/all_SVSR_regenotyped.vcf

echo "nb of SV detected by Delly"
grep -v ^\#\# 04_SR/all_SVSR_regenotyped.vcf | wc -l
```

-> How many SVs are you analysing in total? Same?

-> Do you notice differences with the earlier vcf?

We could expect to have less missing data and more accurate genotype. For example, one can imagine that a given SV was only called in two samples with high-confidence (maybe because they are alternate homozygotes)... Now by regenotyping, we can ask what is the genotype in the other two samples for which it wasn't called (for exemple due to a lower coverage if it was heterozygote or less deeply sequenced). 

The re-genotyping step also refines the breakpoints, and unfortunately re-name the ID of each SVs. While the refinment is an interesting step, the re-ordering of SVs makes it difficult to easily track what happen to a specific SV.





