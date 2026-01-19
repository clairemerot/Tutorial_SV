## 05 - Towards graph-based analysis
Now we will explore genome graphs and see how they can help us match data to study our SVs.

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
Let's focus on sample A, we will come back to the original paired-end sequences in the fastq files. Those will be mapped directly on the graph using vg giraffe and formatted/packed.

```
fq1=~/workshop_materials/structural_variants/SR/A_R1.fastq
fq2=~/workshop_materials/structural_variants/SR/A_R2.fastq
output=05_graphs/A_mapped.gam

vg giraffe -Z 05_graphs/variantgraph.giraffe.gbz -m 05_graphs/variantgraph.shortread.withzip.min -d 05_graphs/variantgraph.dist -f $fq1 -f $fq2 > $output

vg pack -Q 5 -x 05_graphs/variantgraph.giraffe.gbz -g 05_graphs/A_mapped.gam -o 05_graphs/A_mapped.gam.pack
```

If we had a lot of time/computation power, we could do it for many samples sequenced with short-reads. Let's stick to ind A for now and try genotyping the SVs that exists in the graph.

Our first step is to extract the snarls (the bubbles representing the SVs). Then we will call genotype for individual A for all those variants. 

```
vg snarls 05_graphs/variantgraph.giraffe.gbz > 05_graphs/variantgraph.snarls
vg call -a -k 05_graphs/A_mapped.gam.pack -r 05_graphs/variantgraph.snarls -f $ref 05_graphs/variantgraph.giraffe.gbz >05_graphs/A.vcf
```

-> open the vcf: how many variants could you genotype? What does the format look like? 

-> How long was it to run?

VG graphs can also be converted to gfa. You may used this gfa to visualise your graph in Bandage - a desktop app. You cna play with the visualisation and try finding a bubble.

```
 vg convert -f 05_graphs/variantgraph.giraffe.gbz > 05_graphs/variantgraph.giraffe.gfa
```

### Pangenome graph
Before we have built a kind of "variant-aware graph" which is still deeply anchored in a reference. The idea of pangenomes is to get reference-free and to build a pangenome graph which include several assemblies. 

To do so we will use pggb. https://github.com/pangenome/pggb 

I prepared a fasta file which includes our two chromosomes ref/alt as used in section 2 but formatted so that pggb accepts it (format is a little picky). The first thing you will need to do is indexing it:

```
samtools faidx  ~/workshop_materials/structural_variants/assemblies/ref_alt.fasta
```

Next, we will run pggb on this file and build a pangenome graph in gfa format.
```
#we do a folder for the results
mkdir 05_graphs/pggb_graph

pggb -i ~/workshop_materials/structural_variants/assemblies/ref_alt.fasta -o 05_graphs/pggb_graph/
```

This run produces a lot of files. Those includes files for visualisation (*.png) which you may open on the desktop.

-> What do you observe on the 2D representation?
