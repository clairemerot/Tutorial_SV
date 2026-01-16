In this folder we will explore graphs for SV


### Build a graph including variants

```
ref= ~/workshop_materials/structural_variants/assemblies/ref.fasta
vcf_file=~/workshop_materials/structural_variants/ONT/fri59_pass.vcf
prefix=05_graphs/variantgraph

vg autoindex --workflow giraffe -r $ref -v $vcf_file -p $prefix

```

### Genotype other samples on the graph

```
fq1=~/workshop_materials/structural_variants/SR/A_R1.fastq
fq2=~/workshop_materials/structural_variants/SR/A_R2.fastq

vg giraffe -Z 05_graphs/variantgraph.giraffe.gbz -m 05_graphs/variantgraph.shortread.withzip.min -d 05_graphs/variantgraph.dist -f $fq1 -f $fq2 > 05_graphs/A_mapped.gam
```
