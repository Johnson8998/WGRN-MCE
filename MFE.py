import numpy as np

def compute_mfe(pi_k, p_k):
    """
    Compute the Markov Flow Entropy (MFE) for each gene in sample k.
    
    Parameters:
    - pi_k: Normalized expression of genes (1D array-like)
    - p_k: Transfer probability matrix (2D array-like)
    
    Returns:
    - MFE values for each gene (1D array)
    """
    # Calculate MFE for each gene
    mfe_values = -np.sum(pi_k[:, None] * p_k * np.log(pi_k[:, None] * p_k + 1e-10), axis=1)
    return mfe_values

# Example Input
# Normalized expression for all genes
pi_k = np.array([0.1, 0.2, 0.3, 0.4])

# Transfer probability matrix indicating the transfer probabilities to genes
p_k = np.array([
    [0.1, 0.6, 0.2, 0.1],  # Gene 1
    [0.3, 0.2, 0.4, 0.1],  # Gene 2
    [0.2, 0.5, 0.2, 0.1],  # Gene 3
    [0.4, 0.1, 0.3, 0.2]   # Gene 4
])

# Calculate MFE for all genes
mfe_values = compute_mfe(pi_k, p_k)
print("Markov Flow Entropy (MFE) for each gene:", mfe_values)