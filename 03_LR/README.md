## 03 - Assessing breakpoints and detecting other SVs with long-reads
We have an additionnal sample (fri59) which has been sequenced with Oxford nanopore technology providing long-reads (from a few kbs to ~70 kb, with a N50 of 20kb)

I did basic filtering steps and aligned it, using minimap2, on our dummy genome. The resulting bamfile is available here ~/workshop_materials/structural_variants/ONT/indA59.bam

Let's use it to dig into the breakpoints of our large rearrangement. We can wonder "If we had not suspected this inversion from population data, could we have found it with long-reads?". We also have the perfect data to look for all type and length of structural variants along the genome. 

### Call SVs with long-read data
To detect structural variants, we will rely on a widely-used software, Sniffles 2 https://github.com/fritzsedlazeck/Sniffles . 

Be aware that several SV callers are available and, depending on your data, you may favour one over the others or decide on using several tools.

```
sniffles -i ~/workshop_materials/structural_variants/ONT/indA59.bam --reference ~/workshop_materials/structural_variants/assemblies/ref.fasta --minsupport 10 --vcf 03_LR/fri59.vcf

echo "nb of SV detected by Sniffles"
grep -v ^\#\# 03_LR/fri59.vcf | wc -l
```

You can open the vcf and explore it.

-> How many SVs were detected? 

-> What are the different pieces of information available for each variant? 

-> How does it differ from a vcf of SNPs? 

-> What kind of filters or classifications may you want to implement for subsequent analysis?

### Visualise SVs and validate breakpoints 

Next, we will try to visualize our SVs to assess how confident we can be in the detection. 

I suggest to use IGV which may be heavy to use on a large chromosome but will be ok on our dummy genome.

Please go on the desktop and start IGV.

You can upload the genome (ref.fasta), then you can upload your vcf produced with sniffles, and then your bam file. You can try playing with the zoom to visualise the SVs and the read support.

As you may notice in the vcf track we have a lot of long SVs which are most likely wrong. Some have the label "PASS" in the vcf while others have the label "GT" and a genotype 0/0 meaning they have much less upport and are likely wrong.

In IGV, we can hide the ones that do not pass filters (in light blue) directly in the interface. Alternatively, you can filter them on the server with bcftools.

```
bcftools view --apply-filters PASS 03_LR/fri59.vcf > 03_LR/fri59_pass.vcf
```

Please take some time to explore the visualisation. For example, you can try finding a deletion with good support and a deletion which you would be less confident in (same for DUP or INS). 

You can click on each variant in the vcf panel to get more information. In the lower panel you can visualize the coverage in this area and get information about each raw read. 

When you click on each read you can get the information about where it maps, whether it is split and map elsewhere, etc. 

You can also zoom on the breakpoints previously identified with the assembly comparison for the very large rearrangement. 

-> What do they look like here? 

-> Are they supported by long-reads? 

-> Which genotype is this individual "fri59"? 

