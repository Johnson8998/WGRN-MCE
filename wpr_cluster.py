import os
import pandas as pd
import networkx as nx
from sklearn.cluster import KMeans

# Folder path
folder_path = 'D:\\hqcidea1\\network\\MCE1_LUSC_tumor'

# Output file paths
all_genes_output_file = 'MCE1_LUSC_tumor_clustering_results(all).csv'
driver_genes_output_file = 'MCE1_LUSC_tumor_clustering_results(driver).csv'

# Read the driver gene file
driver_genes_df = pd.read_csv('D:\\hqcidea1\\Census_LUNG.csv')
driver_genes = driver_genes_df['Gene Symbol'].values

# Initialize dictionaries to save the results for all genes and driver genes
all_genes_matrix = {}
driver_genes_matrix = {}

# Read files numbered from 1 to 370 in sequence
for i in range(1, 371):
    # Construct the file path
    file_path = os.path.join(folder_path, f'MCE1_LUSC_tumor_{i}.csv')
    
    # Read the data
    data = pd.read_csv(file_path)
    
    # Create a NetworkX graph and add weighted edges
    G = nx.Graph()
    for index, row in data.iterrows():
        G.add_edge(row['Gene1'], row['Gene2'], weight=row['Weight'])
    
    # Use the PageRank algorithm to calculate the weighted PageRank score for each gene
    pagerank_scores = nx.pagerank(G, weight='weight')
    
    # Convert the PageRank scores to a DataFrame format
    pagerank_df = pd.DataFrame(list(pagerank_scores.items()), columns=['Gene', 'PageRank'])
    
    # Use the PageRank scores as features for KMeans clustering
    n_clusters = 10  # You can adjust the number of clusters
    kmeans = KMeans(n_clusters=n_clusters)
    pagerank_df['Cluster'] = kmeans.fit_predict(pagerank_df[['PageRank']])
    
    # Use the file number as the column name
    file_column = f'MCE1_LUSC_tumor_{i}'
    
    # Save the results for all genes to the dictionary
    for gene, cluster in zip(pagerank_df['Gene'], pagerank_df['Cluster']):
        if gene not in all_genes_matrix:
            all_genes_matrix[gene] = {}
        all_genes_matrix[gene][file_column] = cluster
    
    # Filter out the driver genes and save the results to the dictionary
    driver_gene_clusters = pagerank_df[pagerank_df['Gene'].isin(driver_genes)]
    for gene, cluster in zip(driver_gene_clusters['Gene'], driver_gene_clusters['Cluster']):
        if gene not in driver_genes_matrix:
            driver_genes_matrix[gene] = {}
        driver_genes_matrix[gene][file_column] = cluster

# Convert the dictionaries to DataFrames and fill in the missing values
all_genes_df = pd.DataFrame(all_genes_matrix).T.fillna(-1)  # Fill missing values with -1
driver_genes_df = pd.DataFrame(driver_genes_matrix).T.fillna(-1)  # Fill with -1 to indicate the gene is not in the file

# Add the file names as the first row and the gene symbols as the first column
all_genes_df.index.name = 'Gene Symbol'
driver_genes_df.index.name = 'Gene Symbol'

# Save as CSV files
all_genes_df.to_csv(all_genes_output_file)
driver_genes_df.to_csv(driver_genes_output_file)

# Print the save paths
print(f"The clustering results for all genes have been saved in matrix form: {all_genes_output_file}")
print(f"The clustering results for driver genes have been saved in matrix form: {driver_genes_output_file}")