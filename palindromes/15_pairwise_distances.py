from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Read your MFA file
alignment = AlignIO.read("/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/left_fastas/Q6_aligned.mfa_bkp", "fasta")

# Calculate pairwise distances using identity
calculator = DistanceCalculator('identity')
distance_matrix = calculator.get_distance(alignment)

# Convert distances to similarities (identity distances are 1 - similarity)
similarities = np.array([[1 - distance_matrix[i, j] for j in range(len(alignment))] 
                        for i in range(len(alignment))])

# Convert to percentages
similarities_percent = similarities * 100

# Create heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(similarities_percent, annot=True, cmap='RdYlGn', 
            xticklabels=[seq.id for seq in alignment],
            yticklabels=[seq.id for seq in alignment],
            fmt='.1f', vmin=0, vmax=100)
plt.title('Pairwise Sequence Similarity (%)')
plt.savefig("pairwise.png")