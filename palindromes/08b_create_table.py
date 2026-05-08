#!/usr/bin/env python3

import os
import sys
import glob
from collections import defaultdict

def parse_lastz_file(filepath):
    """
    Parse a LASTZ file and extract high-quality matches and all sample sequences.
    Returns a tuple: (matches, all_sample_sequences)
    - matches: list of tuples (name1, name2) where both cov1 and cov2 > 0.95
    - all_sample_sequences: set of all name1 values in the file
    """
    matches = []
    all_sample_sequences = set()
    
    with open(filepath, 'r') as f:
        # Read header
        header_line = f.readline().strip()
        if not header_line.startswith('#'):
            print(f"Warning: {filepath} doesn't have proper header", file=sys.stderr)
            return matches, all_sample_sequences
        
        # Parse header to find column positions
        header_cols = header_line[1:].split('\t')
        try:
            name1_idx = header_cols.index('name1')
            name2_idx = header_cols.index('name2')
            cov1_idx = header_cols.index('cov1')
            cov2_idx = header_cols.index('cov2')
        except ValueError as e:
            print(f"Error: Required column not found in {filepath} - {e}", file=sys.stderr)
            return matches, all_sample_sequences
        
        # Process data rows
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            cols = line.split('\t')
            if len(cols) <= max(name1_idx, name2_idx, cov1_idx, cov2_idx):
                continue
            
            try:
                name1 = cols[name1_idx]
                name2 = cols[name2_idx]
                cov1 = float(cols[cov1_idx])
                cov2 = float(cols[cov2_idx])
                
                # Collect all sample sequences
                all_sample_sequences.add(name1)
                
                # Check if both coverages are > 95%
                if cov1 > 0.95 and cov2 > 0.95:
                    matches.append((name1, name2))
            
            except (ValueError, IndexError):
                continue
    
    return matches, all_sample_sequences

def analyze_all_lastz_files(directory=".", ref = "HG002"):
    """
    Analyze all .lastz files in the directory and create matching table.
    """
    # Find all .lastz files
    lastz_files = glob.glob(os.path.join(directory, f"*{ref}.lastz"))
    
    if not lastz_files:
        print("No .lastz files found in directory", file=sys.stderr)
        return
    
    # Dictionary to store matches: {hg002_name: {sample: [sample_names]}}
    match_table = defaultdict(lambda: defaultdict(list))
    # Dictionary to store all sequences per sample: {sample: set(all_sequences)}
    all_sequences_per_sample = defaultdict(set)
    # Dictionary to store matched sequences per sample: {sample: set(matched_sequences)}
    matched_sequences_per_sample = defaultdict(set)
    
    all_samples = set()
    all_ref_names = set()
    
    # Process each file
    for filepath in sorted(lastz_files):
        filename = os.path.basename(filepath)
        # Extract sample name (everything before first underscore)
        sample_name = filename.split('_')[0]
        all_samples.add(sample_name)
        
        print(f"Processing {filename}...", file=sys.stderr)
        
        matches, all_sample_sequences = parse_lastz_file(filepath)
        
        # Store all sequences for this sample
        all_sequences_per_sample[sample_name] = all_sample_sequences
        
        for name1, name2 in matches:
            match_table[name2][sample_name].append(name1)
            matched_sequences_per_sample[sample_name].add(name1)
            all_ref_names.add(name2)
    
    # Create output table for matched sequences
    all_samples = sorted(all_samples)
    all_ref_names = sorted(all_ref_names)
    
    print(f"=== SEQUENCES THAT MATCH {ref} ===")
    print(f"{ref}_name\t" + "\t".join(all_samples))
    
    # Print data rows for matched sequences
    for ref_name in all_ref_names:
        row = [ref_name]
        for sample in all_samples:
            sample_matches = match_table[ref_name][sample]
            if sample_matches:
                # Join multiple matches with semicolon
                row.append(";".join(sample_matches))
            else:
                row.append("")
        print("\t".join(row))
    
    print(f"\n=== SEQUENCES THAT DON'T MATCH ANY {ref} LOCATION ===")
    print("Sample\tUnmatched_sequences")
    
    # Print unmatched sequences for each sample
    for sample in all_samples:
        all_seqs = all_sequences_per_sample[sample]
        matched_seqs = matched_sequences_per_sample[sample]
        unmatched_seqs = all_seqs - matched_seqs
        
        if unmatched_seqs:
            unmatched_list = sorted(list(unmatched_seqs))
            print(f"{sample}\t{';'.join(unmatched_list)}")
        else:
            print(f"{sample}\t")

def main():
    ref = "HG002"
    directory = "."
    if len(sys.argv) > 1:
        directory = sys.argv[1]
        if len(sys.argv) > 2:
            ref = sys.argv[2]
    
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a directory", file=sys.stderr)
        sys.exit(1)
    
    analyze_all_lastz_files(directory, ref)

if __name__ == "__main__":
    main()