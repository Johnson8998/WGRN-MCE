This is the code associated with the article “Weighted Gene Networks Derived from Multi-Omics Reveal Core Cancer Genes in Lung Cancer”.

## Data Sources
The BioGrid data can be found at https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.226/. TCGA-LUAD and TCGA-LUSC can be found at The Cancer Genome Atlas (TCGA) Research Network: http://cancergenome.nih.gov.


This repository contains the code associated with the article “Weighted Gene Networks Derived from Multi-Omics Reveal Core Cancer Genes in Lung Cancer”.

## Code Files and Their Functions

### 1. `MCE.m`
- **Function**: Calculating Markov chain entropy.
- **Input**:
  - `data`: Gene expression dataset with dimensions `p*n`, where `p` represents the number of genes and `n` represents the number of samples.
  - `net`: A `p*p` adjacency matrix of the gene network (usually a Protein-Protein Interaction (PPI) network). Note that the diagonal elements of this matrix must be 1, not 0.
- **Output**:
  - `ge`: A vector of length `n`, where `ge(i)` is the Markov chain entropy (MCE) of sample `i`.
- **Note**: In the code, `p_ij` represents the transfer probability.

### 2. `MFE`
- **Function**: Calculating Flow chain entropy.

### 3. `wpr_cluster.py`
- **Function**: Calculate weighted PageRank centrality and perform K-means clustering.
- **Input**: A side table in the following format:
```
Gene1  Gene2  Weight
```

