#!/usr/bin/env python3

import argparse
import sys
import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt

# REPORTED_PALINDROMES = ["Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7.7", "Q10"]
REPORTED_PALINDROMES = ["Q1", "Q2", "Q3", "Q4", "Q5", "Q7.9", "Q8.1", "Q9.2"] # Q's for HG02071

# Mapping from Q names to P names
Q_TO_P_MAPPING = {
    "Q1": "P8",
    "Q2": "P7",
    "Q3": "P6",
    "Q4": "P5",
    "Q5": "P4",
    "Q7.9": "P3",
    "Q8.1": "P2",
    "Q9.2": "P1"
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="input table produced by step 08", required=True)
    parser.add_argument("-p", "--path", help="path to lastz self alignments", required=True)
    parser.add_argument("-o", "--out", help="output plot", required=True)
    args = parser.parse_args()

    data_records = []

    table_file = args.input
    with open(table_file,'r') as f:
        for line in f:
            line = line.strip()
            if ("===") in line or line == "":
                if ("DON'T") in line:
                    break
                continue
            fields = line.split("\t")

            if ("name") in fields[0]:
                continue
            
            q_name = ".".join(fields[0].split(".")[1:])
            if q_name in REPORTED_PALINDROMES:
                for field in fields[1:]:
                    samples = field.split(";")
                    if (';' in field):
                        print(field)
                        # if q_name != "Q10":
                        #     continue
                        # else 
                    samples = list(set(samples))
                    for sample in samples:

                        if sample == "":
                            continue

                        contig, palindrome_id = sample.split(".",1)
                        sample_id = contig.split("_")[0]
                        palindrome_file = f"{args.path}/{contig}.pal"

                        with open(palindrome_file, 'r') as file:
                            for palindrome_line in file:
                                palindrome_line = palindrome_line.strip()
                                palindrome_fields = palindrome_line.split("\t")

                                if palindrome_fields[9] == sample:
                                    arm_length = int(palindrome_fields[10])
                                    spacer_length = int(palindrome_fields[11])
                                    identity_pct = float(palindrome_fields[7].replace("%", ""))
                                    repeat_content_pct = float(palindrome_fields[8].replace("%", ""))
                                

                                    record = {
                                        'sample_id': sample_id,
                                        'location': q_name,
                                        'arm_length': arm_length,  # Use np.nan for missing
                                        'spacer_length': spacer_length,
                                        'identity_pct': identity_pct,
                                        'repeat_content_pct':repeat_content_pct
                                    }
                                    data_records.append(record)

    df = pd.DataFrame(data_records)
    
    # Map Q names to P names
    df['location'] = df['location'].map(Q_TO_P_MAPPING)
    
    # Order: P1, P2, P3, P4, P5, P6, P7, P8
    p_order = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"]
    df['location'] = pd.Categorical(df['location'], categories=p_order, ordered=True)
    
    # Calculate counts for each location
    location_counts = df['location'].value_counts().to_dict()
    df['location_with_count'] = df['location'].apply(lambda x: f"{x} (n={location_counts[x]})")
    
    # Create categorical with counts in order
    labels_with_counts = [f"{p} (n={location_counts.get(p, 0)})" for p in p_order]
    df['location_with_count'] = pd.Categorical(df['location_with_count'], categories=labels_with_counts, ordered=True)

    print(df)


    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    features = ['arm_length', 'spacer_length', 'identity_pct', 'repeat_content_pct']
    titles = ['Arm Length', 'Spacer Length', 'Identity between arms (%)', 'Repeat Content (%)']

    for i, (feature, title) in enumerate(zip(features, titles)):
        ax = axes[i//2, i%2]
        sns.boxplot(data=df, x='location_with_count', y=feature, hue='location_with_count', ax=ax, legend=False)
        ax.set_title(title)
        ax.set_xlabel('Location')
        ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(args.out)


if __name__ == "__main__":
    main()