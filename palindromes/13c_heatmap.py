import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser(description='Process components directory')
parser.add_argument('count_table', type=str, 
                        default='output_count_table.csv',
                        help='Path to counts table')

args = parser.parse_args()

# Read the TSV file
df = pd.read_csv(f"{args.count_table}", index_col=0, sep=',')

Q_TO_P_MAPPING = {
    "Q6": "P9",
    "Q1": "P8",
    "Q2": "P7",
    "Q3": "P6",
    "Q4": "P5",
    "Q5": "P4",
    "Q7.7": "P3",
    "Q8": "P2 inversion",
    "Q9.2": "DAZ",
    "Q11.1": "Red",
    "HG02071_chrY.Q9.2": "P1 truncated",
    "HG00621_chrY.Q9.1": "P1",
}
df.columns = [Q_TO_P_MAPPING.get(c, c) for c in df.columns]
df = df.clip(lower=0)

# Remove Q11.1 column if it exists
# if 'Q11.1' in df.columns:
#     df = df.drop('Q11.1', axis=1)

# print(df.index)
# print("HG02984" in df.index)
# print("HG03456" in df.index)
print((df['DAZ'] == 1).sum() )
print((df['DAZ'] == 2).sum() )
print((df['DAZ'] == 3).sum() )
col = df.pop('P9')
df.insert(0, 'P9', col)
col = df.pop("Red")
df.insert(8, "Red", col)
# put p9 as the first begininng
# print(df.columns)
# assert False
# Create a custom colormap - white for 0, increasing shades of blue for 1-3
colors = ['white', '#90CAF9', '#42A5F5', '#1976D2']

artifacts={"P9":[],
           "P8":[],
           "P7":[],
           "P6":[],
           "P5":["HG01243", "HG02602", "HG02083", "NA20805", "HG00140", "NA18612", "HG03130", "HG02258"],
           "P4":["HG02027", "HG02074", "HG00140", "NA21093", "HG02145"],
           "P3":["HG02647", "HG02668", "NA19384", "HG02258", "HG04157", "HG01433", "HG00280", "HG03492"],
           "P1":["HG02145", "HG03098", "NA19347", "HG02145", "HG04187", "HG00140", "HG00321", "NA20809", "HG01433", "HG02132" ]} #both for P1 and P2


from matplotlib.colors import ListedColormap
cmap = ListedColormap(colors)

# Create the plot
fig, ax = plt.subplots(figsize=(12, 18))

# Create the heatmap using pcolormesh for better grid alignment
X, Y = np.meshgrid(np.arange(len(df.columns)+1), np.arange(len(df.index)+1))
im = ax.pcolormesh(X, Y, df.values, cmap=cmap, vmin=0, vmax=3, 
                   edgecolors=(0.5, 0.5, 0.5, 0.3), linewidth=0.5)

# Overlay pink color for artifact cells
for palindrome, samples in artifacts.items():
    # Special handling for P1 - apply to multiple columns
    if palindrome == "P1":
        target_columns = ["P1", "P1 truncated", "P2 inversion"]
    else:
        target_columns = [palindrome]
    
    for col_name in target_columns:
        if col_name not in df.columns:
            continue
        col_idx = df.columns.get_loc(col_name)
        
        for sample in samples:
            if sample not in df.index:
                continue
            row_idx = df.index.get_loc(sample)
            value = df.loc[sample, col_name]
            
            # Check if cell is white (value 0)
            if value != 0:
                raise ValueError(f"Artifact cell [{sample}, {col_name}] has non-zero value: {value}. Cannot overlay pink on non-white cell.")
            
            # Draw pink rectangle that fills the cell exactly
            rect = plt.Rectangle((col_idx, row_idx), 1, 1, 
                                facecolor='pink', edgecolor=(0.5, 0.5, 0.5, 0.3), linewidth=0.5, zorder=3)
            ax.add_patch(rect)

# Set ticks and labels
ax.set_xticks(np.arange(len(df.columns)) + 0.5)
ax.set_yticks(np.arange(len(df.index)) + 0.5)
ax.set_xticklabels(df.columns, rotation=45, ha='right', fontsize=10)
ax.set_yticklabels(df.index, fontsize=8)

# Set axis limits
ax.set_xlim(0, len(df.columns))
ax.set_ylim(0, len(df.index))
ax.invert_yaxis()  # Invert to match standard heatmap orientation

# Create legend with all colors including pink artifact
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='white', edgecolor='gray', label='0'),
    Patch(facecolor='#90CAF9', edgecolor='gray', label='1'),
    Patch(facecolor='#42A5F5', edgecolor='gray', label='2'),
    Patch(facecolor='#1976D2', edgecolor='gray', label='3'),
    Patch(facecolor='pink', edgecolor='gray', label='assembly artefact')
]
ax.legend(handles=legend_elements, title='Count', loc='upper center', 
          bbox_to_anchor=(0.5, 1.15), ncol=5, fontsize=10, framealpha=0.9)

# Set title
plt.title('Palindrome Data Heatmap', fontsize=16, fontweight='bold', pad=20)

# Adjust layout to prevent label cutoff
plt.tight_layout()

# Save the plot
plt.savefig('palindrome_heatmap.png', dpi=300, bbox_inches='tight')


# Show the plot
# plt.show()
plt.savefig('palindrome_heatmap.pdf', bbox_inches='tight')

print("Heatmap saved as 'palindrome_heatmap.png' and 'palindrome_heatmap.pdf'")
