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
 - 

## Code
 - bipartite_breast.py: run code to generate predicted protein expression from mRNA expression, miRNA expression and their interaction network.
 - bipartite_ovarian.py: run code to generate predicted protein expression from mRNA expression, miRNA expression and their interaction network.
