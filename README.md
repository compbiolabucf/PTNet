# PTNet
We provide the source code of paper *In Silico Model for miRNA-medicated Regulatory Network in Cancer* for Nucleic Acids Research submission.

## Goal
 - PTNet is a graph-based learning model which simulates the miRNAs (microRNAs) that regulate gene expression post-transcriptionally *in silico*.
 - Our model estimates the protein levels by considering the miRNA-mRNA interaction network, the mRNA expression and the miRNA expression.

## Framework
<p align="center">
  <img src="PTNet.png" width="400">
  <figcaption>An illustration of the proposed graph-based learning model on miRNA-mRNA bipartite graph to estimate the protein expression levels.</figcaption>
</p>

## Required Python packages
 - Numpy
 - Pandas
 - Scipy

## Data
 - bipartite_targetscan.xlsx: miRNA-mRNA interaction network
 - gencode23.csv: mRNA to gene symbol conversion (annotation)
 - Breast cancer: (1) [mRNA_value.csv]: mRNA expression data, (2) mRNA_sample.csv: sample names for mRNA expression (mRNA_value.csv), (3) mRNA_feature.csv: feature names for mRNA expression, (4) miRNA.csv: miRNA expression data, (5) spectral_count.tsv: ground truth protein expression.

## Code
 - bipartite_breast.py: run code to generate predicted protein expression from mRNA expression, miRNA expression and their interaction network.
 - bipartite_ovarian.py: run code to generate predicted protein expression from mRNA expression, miRNA expression and their interaction network.

[mRNA_value.csv]: <https://drive.google.com/file/d/1Ez5A7Znbrms4O2MJHr0h_02aiCaIY4Cg/view?usp=sharing>
