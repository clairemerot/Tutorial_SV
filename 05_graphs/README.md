## 05 - Towards graph-based analysis
Now we will explore genome graphs and see how they can help us better manipulate SVs.

*** To install?***
See the desktop app maybe?
https://github.com/TF-Chan-Lab/panGraphViewer
https://github.com/vgteam/sequenceTubeMap
https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md

### Build a graph including variants
Using the programm vg (https://github.com/vgteam/vg), we will construct a simple graph which includes the reference path and variants detected earlier using sniffles and our sample fri59. To limit computational time, we will restrict to a "mini" subset of our genome and vcf (<1Mb long).

```
ref=~/workshop_materials/structural_variants/assemblies/ref_mini.fasta
vcf_file=~/workshop_materials/structural_variants/ONT/fri59_filtered_mini.vcf
prefix=05_graphs/variantgraph

vg autoindex --workflow giraffe -r $ref -v $vcf_file -p $prefix

```
We get some warnings. Why? 

### Genotype other samples on the graph
Let's focus on sample A, we will come back to the original paired-end sequences in the fastq files. Those will be mapped directly on the graph using vg giraffe.

```
fq1=~/workshop_materials/structural_variants/SR/A_R1.fastq
fq2=~/workshop_materials/structural_variants/SR/A_R2.fastq
output=05_graphs/A_mapped.gam

vg giraffe -Z 05_graphs/variantgraph.giraffe.gbz -m 05_graphs/variantgraph.shortread.withzip.min -d 05_graphs/variantgraph.dist -f $fq1 -f $fq2 > $output

vg pack -Q 5 -x 05_graphs/variantgraph.giraffe.gbz -g 05_graphs/A_mapped.gam -o 05_graphs/A_mapped.gam.pack

vg snarls 05_graphs/variantgraph.giraffe.gbz > 05_graphs/variantgraph.snarls
vg call -a -k 05_graphs/A_mapped.gam.pack -r 05_graphs/variantgraph.snarls -f $ref 05_graphs/variantgraph.giraffe.gbz >05_graphs/A.vcf
```

-> open the vcf: how many variants could you genotype? What does the format look like? 
-> How long was it to run?

### Visualisation

### Pangenome graph?
