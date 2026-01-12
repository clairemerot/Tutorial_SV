# Tutorial_SV
This is half-a-day training about SV detection using diverse sets of data. 

## 01 - Using SNPs & local PCA to detect haploblocks putatively representing large rearrangements

## 02 - Comparing assemblies for large rearrangements
The ideal dataset to confirmlarge rearrangements would be fully-assembled assemblies, although those will rarely be available for non-model species.
Here, let's say that we are lucky and two different individual have been sequenced with long-reads and assembled as chromosome-scale assembly. One individual has been picked in the A/A group (ref.fasta) and one in the B/B group (alt.fasta).

We will work in the second folder 02_assemblies.

To compare those assemblies we need to align them to each other. For this we will use Minimap 2. https://github.com/lh3/minimap2 . You can check on the Minimap manual what are the different options. This will produce a file listing syntenic and aligned regions in a paf format. 
This file can be opened with a regular text editor to see what it looks like.  

```
minimap2 -x asm20 -c --eqx ~/workshop_materials/structural_variants/assemblies/ref.fasta ~/workshop_materials/structural_variants/assemblies/alt.fasta > 02_assemblies/alt-to-ref20.paf
```

Let's quickly visualise this alignment. D-genies offers a quick dotplot between two genomes. You can either upload the two fasta files or directly the paf. (Please favour the paf file for today to avoid redudant demands on the server)

What do you think? Do you visualise some rearrangements? Does it fit what you uspected from the PCA analysis?

Now we may want to do nicer plots. A recent R package will help you do nice ribbon plot. Let's try SVbyEye. https://github.com/daewoooo/SVbyEye
Here is a very quick code to get a first figure. You can try it in Rstudio. If you have time you may explore more options from this package.

```
library(SVbyEye)
paf.table <- readPaf("02_assemblies/alt-to-ref20.paf",
  include.paf.tags = TRUE, restrict.paf.tags = "cg")
paf.table
plotMiro(paf.table = paf.table, color.by = "direction")
```


## 03 - Assessing breakpoints and detecting other SVs with long-reads

## 04 - Using population dataset of short-reads to detect SVs

## 05 - Towards graph-based analysis
