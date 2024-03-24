# Sex- and Cell type-specific impact of genetic variation associated with neurodevelopmental disorders on the striatum

## Highlights

-	Spiny projection neurons (SPNs) exhibit sex- and cell type-specific transcriptomic changes in 16p11.2 hemideletion (16p11.2 del/+) mice.
-	Sex-specific transcriptomic changes are influenced by wildtype sex differences.
-	In males, but not in females, D2-SPNs specific 16p11.2 del/+ mediates hyperactive behavior.
-	D2-SPNs specific 16p11.2 del/+ in dorsal lateral striatum mediates the male hyperactive behavior.
 
## Summary

Sex is a major factor shaping the manifestation and progression of neurodevelopmental disorders (NDDs), however, it remains poorly understood. Hemideletion of the 16p11.2 region (16p11.2 del/+) is associated with NDDs, displaying sex-specific striatum-related phenotypes resembling NDDs. Striatal circuits, crucial for locomotor control, consist of two distinct pathways: the direct (D1) and indirect (D2) projecting spiny projection neurons (SPNs).  In this study, we establish the causal relationship between the 16p11.2 locus and specific cell types in the striatum in male and female mice. Using snRNA-seq, we identified sex-specific and cell type-specific transcriptomic changes in the striatum of 16p11.2 del/+ mice. Interestingly, these sex-specific transcriptomic changes are influenced by wildtype sex differences. Further pathway analysis unveiled differential gene expression changes linked to synaptic plasticity in D1- and D2-SPNs, as well as GABA signaling pathway changes in D1-SPNs. Besides, 16p11.2 del/+ mediates distinct synaptic function changes in miniature inhibitory postsynaptic current (mIPSC) between D1- and D2-SPNs. Behaviorally, we observe that conditional 16p11.2 del/+ within D2-SPNs phenocopies hyperactivity in male mice, but not within D1-SPNs. Moreover, we found 16p11.2 del/+ within D2-SPNs in the dorsal lateral striatum resulted in a profound hyper locomotor phenotype. Our work reveals that loci linked to NDDs act distinctly in different striatal circuits, impacting behavior in a sex- and cell type-specific manner. This proposes a novel perspective on a potential mechanism of male vulnerability or female resilience to NDDs.

## R Script

The R script "snRNAseq males females WT 16p integrated analysis R script.R" contains the code necessary for the integrated analysis of the single nuclei RNA-seq from males and females, wildtype and 16p11.2 del/+. Before running this script, it is assumed that the CellRanger 10x Genomics pipeline was used on each individual sample.
To read the CellRanger outputs, please replace my absolute path by your own.

## Quadrant plot

The R script "quadrant plot.R" contains the code necessary to reproduce a quadrant plot.
Prior to run this script...

```
join -o 1.1,1.2,2.2 -1 1 -2 1 -e "NA" <(sort -k1,1 SigndSPNs.txt) <(sort -k1,1 dSPNsFemalemale.txt) > dSPNsDelSexComparison.txt
```
