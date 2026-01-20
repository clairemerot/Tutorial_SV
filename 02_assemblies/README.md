## 02 - Comparing assemblies for large rearrangements
The ideal dataset to confirm large rearrangements would be fully-assembled assemblies, although those will rarely be available for non-model species.

Here, let's say that we are lucky and two different individual have been sequenced with long-reads and assembled as chromosome-scale assembly. One individual has been picked in the AA group (ref.fasta) and one in the BB group (alt.fasta).

We will work in the second folder 02_assemblies.

### Assemblies alignement
To compare those assemblies we need to align them to each other. For this we will use Minimap 2. https://github.com/lh3/minimap2 . You can check on the Minimap manual what are the different options. This will produce a file listing syntenic and aligned regions in a paf format. 
This file can be opened with a regular text editor to see what it looks like.  

```
ref=~/workshop_materials/structural_variants/assemblies/ref.fasta
alt=~/workshop_materials/structural_variants/assemblies/alt.fasta
paf_output=02_assemblies/alt-to-ref20.paf
minimap2 -x asm20 -c --eqx $ref $alt > $paf_output
```

Let's quickly visualise this alignment. D-genies offers a quick dotplot between two genomes. You can either upload the two fasta files or directly the paf. 

Please favour the paf file for today to avoid redudant demands on the server. You need to go in the second tab "plot alignment", then upload your paf file in "Alignement file", the ref.fasta in the "Target file" and the alt.fasta in the "Query file". The paf file is inside your Tutorial_SV/02_assemblies while the fasta are in our data dungeon "~/workshop_materials/structural_variants/assemblies/"

What do you think? Do you visualise some rearrangements? Does it fit what you suspected from the PCA analysis?

You can further explore the interactive view of D-Genies. 

### Ribbon plots and exploration of the assemblies
Now we may want to do nicer plots. A recent R package will help you do nice ribbon plot. Let's try SVbyEye. https://github.com/daewoooo/SVbyEye
Here is a very quick code to get a first figure. You can try it in Rstudio. If you have time you may explore more options from this package.

```
library(SVbyEye)
paf.table <- readPaf("02_assemblies/alt-to-ref20.paf",
  include.paf.tags = TRUE, restrict.paf.tags = "cg")
paf.table
plotMiro(paf.table = paf.table, color.by = "direction")
```

Using the paf file and/or the R package, can you extract the breakpoints of the suspected rearrangement? How does it compare to the putative breakpoints identified by indirect PCA methods based on recombination suppression?
