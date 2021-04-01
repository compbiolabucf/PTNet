# PTNet
We provide the source code of paper *In Silico Model for miRNA-mediated Regulatory Network in Cancer* for Briefings in Bioinformatics submission.

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

## Sample data
 - **bipartite_targetscan.csv**: miRNA-mRNA interaction network, please note that the interaction network is an mRNA by miRNA matrix in .csv format. Each row represents one mRNA isoform, each column represents one miRNA, value equals to 1 means that there is an interaction between the corresponding mRNA (row index) and miRNA (column index), 0 means no interaction.
 - **Ovarian cancer RNA-Seq data**: 
   * **[mRNA.csv]**: mRNA expression data of TCGA Ovarian cancer patient samples, please note that the mRNA expression data is a feature by sample matrix in .csv format. Each row represents one mRNA isoform, and each column represents one patient sample. 
   * **miRNA.csv**: miRNA expression data of TCGA Ovarian cancer patient samples, please note that the miRNA expression data is a feature by sample matrix in .csv format. Each row represents one miRNA, and each column represents one patient sample.

## Code
 - **PTNet.py**: run code to generate predicted protein expression from mRNA expression, miRNA expression and their interaction network.
 - **Command to run the code**: ``` $ python3 PTNet.py mRNA_expression miRNA_expression miRNA-mRNA_interaction_network alpha ```. **alpha** is the coefficient in bipartite network propagation ranging from 0 to 1. Higher value of alpha puts more emphasis on the network, lower value puts more emphasis on the initial mRNA values.

## Run PTNet with TCGA Ovarian cancer patient samples 
```sh
$ python3 PTNet.py mRNA.csv miRNA.csv bipartite_targetscan.csv 0.6 
```
All input data and output should be in the same folder as the code. **predicted_protein.csv** will be generated as output which contains the protein expressions estimated by PTNet.

[mRNA.csv]: <https://drive.google.com/file/d/18WrnFyqQcp7GjZc9YdvTtt6acJTHYkLU/view?usp=sharing>
