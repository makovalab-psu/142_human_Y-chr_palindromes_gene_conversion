#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

# Columns from output_elements_table.csv to include (Q6 has no P mapping)
REPORTED_PALINDROMES = ["Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7.7", "Q8", "Q9.2", "HG02071_chrY.Q9.2", "HG00621_chrY.Q9.1"]
 
  

# Mapping from Q names to P names
Q_TO_P_MAPPING = {
    "Q6": "P9",
    "Q1": "P8",
    "Q2": "P7",
    "Q3": "P6",
    "Q4": "P5",
    "Q5": "P4",
    "Q7.7": "P3",
    "Q8": "P2_inversion",
    "Q9.2": "DAZ",
    "HG02071_chrY.Q9.2": "P1",
    "HG00621_chrY.Q9.1" : "P1"
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="output_elements_table.csv", required=True)
    parser.add_argument("-p", "--path", help="path to lastz self alignments", required=True)
    parser.add_argument("-o", "--out", help="output plot", required=True)
    args = parser.parse_args()

    data_records = []

    df_csv = pd.read_csv(args.input, index_col=0)

    for q_name in REPORTED_PALINDROMES:
        if q_name not in df_csv.columns:
            continue
        for sample_id, cell in df_csv[q_name].items():
            if pd.isna(cell) or str(cell).strip() == "":
                continue
            palindrome_ids = list(set(p.strip() for p in str(cell).split(",") if p.strip()))
            for sample in palindrome_ids:
                contig = sample.split(".", 1)[0]
                palindrome_file = f"{args.path}/{contig}.pal"

                with open(palindrome_file, 'r') as file:
                    for palindrome_line in file:
                        palindrome_fields = palindrome_line.strip().split("\t")
                        if len(palindrome_fields) > 11 and palindrome_fields[9] == sample:
                            arm_length = int(palindrome_fields[10])
                            spacer_length = int(palindrome_fields[11])
                            identity_pct = float(palindrome_fields[7].replace("%", ""))
                            repeat_content_pct = float(palindrome_fields[8].replace("%", ""))
                            data_records.append({
                                'sample_id': sample_id,
                                'location': q_name,
                                'arm_length': arm_length,
                                'spacer_length': spacer_length,
                                'identity_pct': identity_pct,
                                'repeat_content_pct': repeat_content_pct
                            })

    df = pd.DataFrame(data_records)
    
    # Map Q names to P names
    df['location'] = df['location'].map(Q_TO_P_MAPPING)
    
    # Order: P1, P2, P3, P4, P5, P6, P7, P8
    p_order = ["P1", "P2_inversion", "DAZ", "P3", "P4", "P5", "P6", "P7", "P8", "P9"]
    df['location'] = pd.Categorical(df['location'], categories=p_order, ordered=True)
    
    # Calculate counts for each location
    location_counts = df['location'].value_counts().to_dict()
    df['location_with_count'] = df['location'].apply(lambda x: f"{x} (n={location_counts[x]})")
    
    # Create categorical with counts in order
    # labels_with_counts = [f"{p} (n={location_counts.get(p, 0)})" for p in p_order]
    # df['location_with_count'] = pd.Categorical(df['location_with_count'], categories=labels_with_counts, ordered=True)

    print(df)
    df.to_csv(f'{args.out}.csv')

    print("\nSamples per location:")
    for loc in p_order:
        print(f"  {loc}: {location_counts.get(loc, 0)}")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    df['spacer_length_kb'] = df['spacer_length'] / 1000

    # features = ['arm_length', 'spacer_length', 'identity_pct', 'repeat_content_pct']
    features = ['arm_length', 'spacer_length_kb', 'identity_pct', 'repeat_content_pct']
    titles = ['Arm length', 'Spacer length', 'Identity between arms', 'Repeat content']

    label_configs = {
        'arm_length':         {'fmt': lambda v: f'{v/1000:.1f}K', 'offset_px': 6},
        'spacer_length_kb':      {'fmt': lambda v: f'{v:.1f}K', 'offset_px': 6},
        'identity_pct':       {'fmt': lambda v: f'{v:.3f}',       'offset_px': -8},
        'repeat_content_pct': {'fmt': lambda v: f'{v:.2f}',       'offset_px': 6},
    }

    panel_labels = ['a', 'b', 'c', 'd']
    y_labels = ['Length in megabase pairs', 'Lenght in kilobase pairs', 'Percentage', 'Percentage']

    for i, (feature, title) in enumerate(zip(features, titles)):
        ax = axes[i//2, i%2]
        sns.boxplot(data=df, x='location', y=feature, hue='location', ax=ax, legend=False)
        ax.set_title(title)
        ax.set_xlabel('')
        ax.set_ylabel(y_labels[i])
        ax.tick_params(axis='x', rotation=45)
        ax.yaxis.offsetText.set_visible(False)
        ax.text(-0.12, 1.05, panel_labels[i], transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top')
        cfg = label_configs[feature]
        medians = df.groupby('location', observed=True)[feature].median()
        for j, loc in enumerate(p_order):
            if loc in medians:
                ax.annotate(
                    cfg['fmt'](medians[loc]),
                    xy=(j, medians[loc]),
                    xytext=(0, cfg['offset_px']),
                    textcoords='offset pixels',
                    ha='center',
                    va='bottom' if cfg['offset_px'] >= 0 else 'top',
                    fontsize=8, fontweight='bold'
                )

    plt.tight_layout()
    plt.savefig(args.out)


if __name__ == "__main__":
    main()
