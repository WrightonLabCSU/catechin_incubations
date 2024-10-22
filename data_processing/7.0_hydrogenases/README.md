## step 1: search Hydrogenase HMMs against MAG database genes
> note, using CAMPER to do the HMM search
```
camper_annotate -i genes.faa -o hyd_annotations --camper_hmm_loc 7.1_hydrogenases.hmm --camper_hmm_cutoffs_loc 7.2_hydrogenase_HMM_scores.txt --threads 15
```

## step 2: pull amino acid seqs of genes that get annotated in step 1
```
#!/bin/bash
#$1 list of genes
for i in $(<$1)
do
grep -A 1 "$i" genes.faa >> fefe_genes.faa
done
echo
```

```
#!/bin/bash
#$1 list of genes
for i in $(<$1)
do
grep -A 1 "$i" genes.faa >> nife_genes.faa
done
echo
```

## step 3: for both fefe and nife hydrogenases, combine with reference sequences, align, and build tree
```
cat fefe_genes.faa fefe_refs.faa > combined_ref_pulled_fefe.faa

#Align your sequences, allow mafft to select the approach:
mafft --auto --anysymbol combined_ref_pulled_fefe.faa > combined_ref_pulled_fefe.mafftauto.faa

#Trim the alignment
trimal -gappyout -in combined_ref_pulled_fefe.mafftauto.faa -out combined_ref_pulled_fefe.mafftauto.gappyout.faa

#Run iqtree with extended model selection (ModelFinderPlus) followed by tree inference, 1000 fast bootstraps using
# UFBoot and SH-aLRT test  (see iqtree FAQ for explanation)
#See http://www.iqtree.org/doc/ for documentation
iqtree -s  combined_ref_pulled_fefe.mafftauto.gappyout.faa -alrt 1000 -bb 1000 -m MFP -nt AUTO -ntmax 10
```

```
cat nife_genes.faa nife_refs.faa > combined_ref_pulled_nife.faa

#Align your sequences, allow mafft to select the approach:
mafft --auto --anysymbol combined_ref_pulled_nife.faa > combined_ref_pulled_nife.mafftauto.faa

#Trim the alignment
trimal -gappyout -in combined_ref_pulled_nife.mafftauto.faa -out combined_ref_pulled_nife.mafftauto.gappyout.faa

#Run iqtree with extended model selection (ModelFinderPlus) followed by tree inference, 1000 fast bootstraps using
# UFBoot and SH-aLRT test  (see iqtree FAQ for explanation)
#See http://www.iqtree.org/doc/ for documentation
iqtree -s  combined_ref_pulled_nife.mafftauto.gappyout.faa -alrt 1000 -bb 1000 -m MFP -nt AUTO -ntmax 10
```
## step 4: visualize in iTOL and pull placement relative to references
Use `7.3_combined_ref_pulled_fefe.mafftauto.gappyout.faa.treefile` and `7.4_combined_ref_pulled_nife.mafftauto.gappyout.faa.treefile`

## step 5: make nicer looking tree in ggtree
Use `7.5_hydrogenase_trees.R`




