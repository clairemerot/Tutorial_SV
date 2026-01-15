


```
delly call -t ALL \
-o 04_SR/allSV_SR.bcf \
-g  ~/workshop_materials/structural_variants/assemblies/ref.fasta \
-q 20 \
~/workshop_materials/structural_variants/SR/B.bam 

#-t to choose the type of SVs - here all
# -o the output file
# -g the reference genome
# -q 20 to keep alignement above a minimum quality of 20
# -t using 4 threads
#-s
```

It is also possible to do all the samples at once - but slightly longer. You may want to let it run while visualising the call from the 1st run.

```
delly call -t ALL \
-o 04_SR/allSV_SR.bcf \
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
# -t using 4 threads
#-s
```
