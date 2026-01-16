### Calling SV with Delly
Although short-reads are short (usually around 150bp), they come in pairs (here they are paired end) and this provides some information to detect structural variants. For example, if the pairs are further away than expected given the insert size, one may suspect an insertion in this area. If the paired-end map in reverse direction, they may signals an inversion. Finally the coverage provides a lot of data about copy-number-variants such as duplications, deletions, etc... Based on this idea, several SV callers have emerged a few years ago. However, athough those programms do their best, they are also confounded by repeats, mapping issues of short-reads, complex areas of the genome, etc leading sometimes to high rates of false positives. Detecting SVs from short-reads is nevertheless valuable given the most efficent cost of short-reads (relatively to long-reads and assemblies), which allows to really get variants at population level (on many samples!).

Here we will use just one possible caller but it is recommended to use several callers and carefully consider the outcomes. I picked Delly as it is still maintained, can also be run on long-reads and is easy to use.

The caller 

```
delly call -t ALL \
-o 04_SR/B_SVSR.bcf \
-g  ~/workshop_materials/structural_variants/assemblies/ref.fasta \
-q 20 \
~/workshop_materials/structural_variants/SR/B.bam 

#-t to choose the type of SVs - here all
# -o the output file
# -g the reference genome
# -q 20 to keep alignement above a minimum quality of 20


#next lets convert the bcf to vcf
bcftools view 04_SR/B_SVSR.bcf > 04_SR/B_SVSR.vcf

echo "nb of SV detected by Delly"
grep -v ^\#\# 04_SR/B_SVSR.vcf | wc -l
```

It is also possible to do all the samples at once - but slightly longer. Please let it run while we will go visualising the call from the 1st run.

```
delly call -t ALL \
-o 04_SR/all_SVSR.bcf \
-g  ~/workshop_materials/structural_variants/assemblies/ref.fasta \
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
