# A chromosome region linked to neurodevelopmental disorders acts in distinct neuronal circuits in males and females to control locomotor behavior

## Highlights

-	16p11.2 hemideletion (16p11.2 del/+) induces sex- and cell type-specific transcriptomic signatures in spiny projection neurons (SPNs).
- Transcriptomic changes in GABA signaling in D1-SPNs align with changes in inhibitory synapse function.
- 16p11.2 del/+ in D2-SPNs causes hyperactivity in males but not females.
- 16p11.2 del/+ in D2-SPNs in the dorsal lateral striatum drives hyperactivity in males.
- 16p11.2 del/+ in cortex drives hyperactivity in both sexes.
 
## Summary

Biological sex shapes the manifestation and progression of neurodevelopmental disorders (NDDs). These disorders often demonstrate male-specific vulnerabilities; however, the identification of underlying mechanisms remains a significant challenge in the field. Hemideletion of the 16p11.2 region (16p11.2 del/+) is associated with NDDs, and mice modeling 16p11.2 del/+ exhibit sex-specific striatum-related phenotypes
relevant to NDDs. Striatal circuits, crucial for locomotor control, consist of two distinct pathways: the direct and indirect pathways originating from D1 dopamine receptor (D1R) and D2 dopamine receptor (D2R) expressing spiny projection neurons (SPNs), respectively. In this study, we define the impact of 16p11.2 del/+ on striatal circuits in male and female mice. Using snRNA-seq, we identify sex- and cell type-specific
transcriptomic changes in the D1- and D2-SPNs of 16p11.2 del/+ mice, indicating distinct transcriptomic signatures in D1-SPNs and D2-SPNs in males and females, with a ~5-fold greater impact in males. Further pathway analysis reveals differential gene expression changes in 16p11.2 del/+ male mice linked to synaptic plasticity in D1- and D2-SPNs and GABA signaling pathway changes in D1-SPNs. Consistent with our
snRNA-seq study revealing changes in GABA signaling pathways, we observe distinct changes in miniature inhibitory postsynaptic currents (mIPSCs) in D1- and D2-SPNs from 16p11.2 del/+ male mice. Behaviorally, we utilize conditional genetic approaches to introduce the hemideletion selectively in either D1- or D2-SPNs and find that conditional hemideletion of genes in the 16p11.2 region in D2-SPNs causes
hyperactivity in male mice, but hemideletion in D1-SPNs does not. Within the striatum, hemideletion of genes in D2-SPNs in the dorsal lateral striatum leads to hyperactivity in males, demonstrating the importance of this striatal region. Interestingly, conditional 16p11.2 del/+ within the cortex drives hyperactivity in both sexes. Our work reveals that a locus linked to NDDs acts in different striatal circuits, selectively impacting behavior in a sex- and cell type-specific manner, providing new insight into male vulnerability for NDDs.

## Pre-processing raw fastq files with CellRanger 10x Genomics pipeline

This is an example of bash script to submit the CellRanger 10x Genomics pipeline on University of Iowa Argon cluster node for one sample:

```
#!/bin/bash
#$ -m e

cellranger count --id=Sample796 --transcriptome=/Shared/NEURO/AbelLab/refdata-cellranger-mm10-3.0.0 --fastqs=/Shared/NEURO/AbelLab/Yann/scRNA-seq_16p11/raw_reads/Sample_796/ --sample=796_Lane7_20190306000,796_Lane7_20190306001,796_Lane7_20190306002,796_Lane7_20190306003,796_Lane7_20190313000,796_Lane7_20190313001,796_Lane7_20190313002,796_Lane7_20190313003,796_Lane8_20190306000,796_Lane8_20190306001,796_Lane8_20190306002,796_Lane8_20190306003,796_Lane8_20190313000,796_Lane8_20190313001,796_Lane8_20190313002,796_Lane8_20190313003
```

## R Script

The R script "snRNAseq males females WT 16p integrated analysis R script.R" contains the code necessary for the integrated analysis of the single nuclei RNA-seq from males and females, wildtype and 16p11.2 del/+. after the CellRanger 10x Genomics pipeline was used on each individual sample.
To read the CellRanger outputs, please replace my absolute path by your own such as:

```
hetmale1.data <- Read10X(data.dir = "/your_own_path/Sample796/outs/filtered_feature_bc_matrix")
```

## Quadrant plot

The R script "quadrant plot.R" contains the code necessary to reproduce a quadrant plot.
Prior to run this script, the input file is generated as followed:

```
join -o 1.1,1.2,2.2 -1 1 -2 1 -e "NA" <(sort -k1,1 SigndSPNs.txt) <(sort -k1,1 dSPNsFemalemale.txt) > dSPNsDelSexComparison.txt
```

where SigndSPNs.txt is a tab delimited text file with one column of gene names and one column of fold change from the differential gene expression table of dSPNs between WT males and 16p11.2 del males.
dSPNsFemalemale.txt is a tab delimited text file with one column of gene names and one column of fold change from the differential gene expression table of dSPNs between WT males and WT females.
