from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os

def calculate_msa_similarity_stats(mfa_file):
    """Calculate similarity statistics for a single MSA file"""
    try:
        # Read the alignment
        alignment = AlignIO.read(mfa_file, "fasta")
        
        # Calculate pairwise distances using identity
        calculator = DistanceCalculator('identity')
        distance_matrix = calculator.get_distance(alignment)
        
        # Convert distances to similarities
        n = len(alignment)
        similarities = np.array([[1 - distance_matrix[i, j] for j in range(n)] 
                                for i in range(n)])
        
        # Convert to percentages
        similarities_percent = similarities * 100
        
        # Get upper triangle (excluding diagonal) for summary stats
        upper_triangle = similarities_percent[np.triu_indices_from(similarities_percent, k=1)]
        
        # Calculate statistics
        stats = {
            'mean': np.mean(upper_triangle),
            'median': np.median(upper_triangle),
            'min': np.min(upper_triangle),
            'max': np.max(upper_triangle),
            'std': np.std(upper_triangle),
            'n_sequences': len(alignment),
            'n_comparisons': len(upper_triangle),
            'similarities': upper_triangle
        }
        
        return stats
        
    except Exception as e:
        print(f"Error processing {mfa_file}: {e}")
        return None

# List of your MFA files
mfa_files = [
    "/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/left_fastas/Q1_aligned.mfa_bkp",
    "/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/left_fastas/Q2_aligned.mfa_bkp", 
    "/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/left_fastas/Q3_aligned.mfa_bkp",
    "/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/left_fastas/Q4_aligned.mfa_bkp",
    "/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/left_fastas/Q5_aligned.mfa_bkp",
    "/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/left_fastas/Q6_aligned.mfa_bkp",
    "/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/left_fastas/Q7.7_aligned.mfa_bkp",
    "/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/left_fastas/Q10_aligned.mfa"

]

# Calculate statistics for each MSA
results = []
violin_data = []  # For violin plot

if os.path.exists("pairwise.csv"):
    print("CSV exists")
    df = pd.read_csv("pairwise.csv")
    # Recalculate for violin plot data
    for i, mfa_file in enumerate(mfa_files):
        stats = calculate_msa_similarity_stats(mfa_file)
        if stats:
            # msa_name = f"MSA {i+1}"
            msa_name = mfa_file.split('/')[-1].split("_")[0]

            for similarity in stats['similarities']:
                violin_data.append({'msa_name': msa_name, 'similarity': similarity})
else:
    for i, mfa_file in enumerate(mfa_files):
        stats = calculate_msa_similarity_stats(mfa_file)
        if stats:
            msa_name = f"MSA {i+1}"
            stats['msa_name'] = msa_name
            stats['file'] = mfa_file
            results.append(stats)
            
            # Collect data for violin plot
            for similarity in stats['similarities']:
                violin_data.append({'msa_name': msa_name, 'similarity': similarity})

    # Convert to DataFrame for easier plotting
    df = pd.DataFrame(results)
    df.to_csv("pairwise.csv")

violin_df = pd.DataFrame(violin_data)
# Create bar plot
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
fig.suptitle('Sequence Similarity Statistics Across MSAs', fontsize=16, fontweight='bold')

# Mean similarity
axes[0, 0].bar(df['msa_name'], df['mean'], color='skyblue', alpha=0.7)
axes[0, 0].set_title('Mean Similarity (%)')
axes[0, 0].set_ylabel('Similarity (%)')
axes[0, 0].tick_params(axis='x', rotation=45)

# Standard deviation (variability)
axes[0, 1].bar(df['msa_name'], df['std'], color='lightcoral', alpha=0.7)
axes[0, 1].set_title('Standard Deviation of Similarity (%)')
axes[0, 1].set_ylabel('Standard Deviation (%)')
axes[0, 1].tick_params(axis='x', rotation=45)

# Min and Max similarity
x_pos = np.arange(len(df))
width = 0.35
axes[1, 0].bar(x_pos - width/2, df['min'], width, label='Min', color='orange', alpha=0.7)
axes[1, 0].bar(x_pos + width/2, df['max'], width, label='Max', color='green', alpha=0.7)
axes[1, 0].set_title('Min and Max Similarity (%)')
axes[1, 0].set_ylabel('Similarity (%)')
axes[1, 0].set_xticks(x_pos)
axes[1, 0].set_xticklabels(df['msa_name'], rotation=45)
axes[1, 0].legend()

# Number of sequences
axes[1, 1].bar(df['msa_name'], df['n_sequences'], color='mediumpurple', alpha=0.7)
axes[1, 1].set_title('Number of Sequences')
axes[1, 1].set_ylabel('Count')
axes[1, 1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.show()
plt.savefig("parwise_min_max_labels.png")

# Print summary table
print("\nSummary Statistics for All MSAs:")
print("=" * 80)
print(f"{'MSA':<8} {'Mean':<8} {'Median':<8} {'Min':<8} {'Max':<8} {'Std':<8} {'N_Seqs':<8}")
print("-" * 80)
for _, row in df.iterrows():
    print(f"{row['msa_name']:<8} {row['mean']:<8.2f} {row['median']:<8.2f} "
          f"{row['min']:<8.2f} {row['max']:<8.2f} {row['std']:<8.2f} {row['n_sequences']:<8}")

# Alternative: Simple single bar plot showing mean similarity
plt.figure(figsize=(10, 6))
bars = plt.bar(df['msa_name'], df['mean'], color='steelblue', alpha=0.7)
plt.title('Mean Pairwise Similarity Across MSAs', fontsize=14, fontweight='bold')
plt.ylabel('Mean Similarity (%)', fontsize=12)
plt.xlabel('Multiple Sequence Alignment', fontsize=12)

# Add value labels on bars
for bar, value in zip(bars, df['mean']):
    plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
             f'{value:.1f}%', ha='center', va='bottom', fontweight='bold')

plt.xticks(rotation=45)
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.show()
plt.savefig("parwise_mean_labels.png")

# Optional: Create a comparison plot showing variability
plt.figure(figsize=(12, 6))
plt.errorbar(df['msa_name'], df['mean'], yerr=df['std'], 
             fmt='o', capsize=5, capthick=2, ecolor='red', 
             markersize=8, color='steelblue', markerfacecolor='lightblue')
plt.title('Mean Similarity with Standard Deviation', fontsize=14, fontweight='bold')
plt.ylabel('Similarity (%)', fontsize=12)
plt.xlabel('Multiple Sequence Alignment', fontsize=12)
plt.grid(axis='y', alpha=0.3)
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig("parwise_var.png")


# New: Violin plot showing distribution of similarities
plt.figure(figsize=(12, 8))
sns.boxplot(data=violin_df, x='msa_name', y='similarity')
plt.title('Distribution of Pairwise Similarities Across MSAs', fontsize=14, fontweight='bold')
plt.ylabel('Similarity (%)', fontsize=12)
plt.xlabel('Multiple Sequence Alignment', fontsize=12)
plt.xticks(rotation=45)
plt.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig("pairwise_violin.png")
plt.show()