This is the code associated with the article “Weighted Gene Networks Derived from Multi-Omics Reveal Core Cancer Genes in Lung Cancer”.

MCE.m : Calculating Markov chain entropy
Input: 
data: gene expression dataset, p*n where p is the number of genes and n is the number of samples.
net: p*p adjacency matrix of gene network (usually PPI). The diagnal is 1 NOT 0!!!
Output:
ge: a length n vector, ge(i) is the MCE of sample i
In the code p_ij is the transfer probability


MFE : Calculating Flow chain entropy

wpr_cluster.py : Calculate weighted pagerank centrality and use Kmeans clustering. The input is a side table in the following format:
Gene1  Gene2  Weight
