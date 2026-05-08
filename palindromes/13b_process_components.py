import argparse
import os
import pandas as pd

from datetime import datetime




def main():
    parser = argparse.ArgumentParser(description='Process components directory')
    parser.add_argument('components_dir', type=str, nargs='?', 
                        default='data/palindromes_A8K_S500K/all_v_all/individual_alignments/subtables',
                        help='Path to the components directory')
    parser.add_argument('sample_list', type=str, nargs='?',
                        default='data/phylotree/list_of_IDs.txt',
                        help='Path to the sample list file')
    parser.add_argument('palindrome_dir', type=str, nargs='?',
                        default='data/palindromes_A8K_S500K',
                        help='Path to the palindrome directory')
    args = parser.parse_args()
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    components_dir = args.components_dir
    sample_list = args.sample_list
    palindrome_dir = args.palindrome_dir
    print(f"Processing components directory: {components_dir}")
    print(f"Sample list: {sample_list}")
    print(f"Palindrome directory: {palindrome_dir}")

    # Dictionary to store Q id -> list of samples
    palindrome_dict = {}
    
    # Dictionary for specific elements
    specific_elements = {
        'HG02071_chrY.Q9.2': None,
        'HG00621_chrY.Q9.1': None,
        
    }

    #construct map of homologous palindromes
    for filename in os.listdir(components_dir):
        filepath = os.path.join(components_dir, filename)
        if os.path.isfile(filepath):
            with open(filepath, 'r') as f:
                first_line = f.readline().strip()
                elements = first_line.split(',')

                # Check for specific elements
                for specific_key in specific_elements.keys():
                    if any(specific_key in elem for elem in elements):
                        specific_elements[specific_key] = elements
                        print(f"File: {filename}, Found {specific_key}, Total samples: {len(elements)}")
                
                # Check if HG002 is in any element
                hg002_element = None
                for element in elements:
                    if 'HG002_' in element:
                        hg002_element = element
                        break
                
                # If HG002 found, extract Q id and store all elements
                if hg002_element:
                    if 'Q' in hg002_element:
                        palindrome_name = hg002_element[hg002_element.index('Q'):]
                        
                        # Add all elements to dictionary
                        if palindrome_name not in palindrome_dict:
                            palindrome_dict[palindrome_name] = []
                        palindrome_dict[palindrome_name].extend(elements)
                        
                        print(f"File: {filename}, HG002: {hg002_element}, Palindrome: {palindrome_name}, Total samples: {len(elements)}")
                    else:
                        print(f"File: {filename}, HG002 element has no Q: {hg002_element}")
    
    # Print summary
    print(f"\nFound {len(palindrome_dict)} unique palindromes:")
    for q_id, samples in palindrome_dict.items():
        print(f"{q_id}: {len(samples)} samples")
    
    # Load sample list
    print(f"\nLoading sample list from: {sample_list}")
    with open(sample_list, 'r') as f:
        sample_names = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(sample_names)} samples")
    
    # Define specific Q ids to report in order
    desired_q_ids = ['Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7.7', 'Q9.2', 'Q8', 'Q11.1']
    
    # Filter to only include Q ids that exist in palindrome_dict
    q_ids = [q_id for q_id in desired_q_ids if q_id in palindrome_dict]
    
    # Add specific element columns
    all_columns = q_ids + ['HG02071_chrY.Q9.2', 'HG00621_chrY.Q9.1']
    
    # Create first table: elements containing each sample
    elements_table = pd.DataFrame(index=sample_names, columns=all_columns)
    
    # Create second table: count of elements
    count_table = pd.DataFrame(index=sample_names, columns=all_columns)
    
    # Fill the tables for Q ids
    for q_id in q_ids:
        elements_list = palindrome_dict[q_id]
        for sample in sample_names:
            # Find all elements starting with this sample name followed by underscore
            matching_elements = [elem for elem in elements_list if elem.startswith(sample + '_')]
            
            # Store in tables
            elements_table.at[sample, q_id] = ','.join(matching_elements) if matching_elements else ''
            count_table.at[sample, q_id] = len(matching_elements)
    
    # Fill the tables for specific elements
    for specific_key, elements_list in specific_elements.items():
        if elements_list is not None:
            for sample in sample_names:
                # Find all elements starting with this sample name followed by underscore
                matching_elements = [elem for elem in elements_list if elem.startswith(sample + '_')]
                
                # Store in tables
                elements_table.at[sample, specific_key] = ','.join(matching_elements) if matching_elements else ''
                count_table.at[sample, specific_key] = len(matching_elements)
    
    # Save tables
    elements_output = f"output_elements_table{timestamp}.csv"
    count_output = f"output_count_table_{timestamp}.csv"
    
    elements_table.to_csv(elements_output)
    count_table.to_csv(count_output)
    
    print(f"\nTables saved:")
    print(f"  Elements table: {elements_output}")
    print(f"  Count table: {count_output}")
    
    # Create BED file from all elements in the table
    print(f"\nCreating BED file from elements...")
    bed_entries = []
    
    # Iterate through all cells in the elements table
    for sample in sample_names:
        for col_header in all_columns:
            cell_value = elements_table.at[sample, col_header]
            if cell_value:  # If not empty
                # Split by comma in case there are multiple elements
                elements_in_cell = cell_value.split(',')
                
                for element in elements_in_cell:
                    element = element.strip()
                    if not element:
                        continue
                    
                    # Extract sample name and Q id from element (e.g., HG03098_chrY.Q2)
                    if '_' not in element or '.' not in element:
                        continue
                    
                    # Parse element to get filename
                    parts = element.split('.')
                    if len(parts) < 2:
                        continue
                    
                    base_name = parts[0]  # e.g., HG03098_chrY
                    q_name = '.'.join(parts[1:])  # e.g., Q2
                    
                    # Construct .pal file path
                    pal_file = os.path.join(palindrome_dir, f"{base_name}.pal")
                    
                    if not os.path.isfile(pal_file):
                        continue
                    
                    # Read the .pal file and find the matching line
                    try:
                        with open(pal_file, 'r') as f:
                            for line in f:
                                if line.startswith('#'):
                                    continue
                                
                                fields = line.strip().split('\t')
                                if len(fields) < 10:
                                    continue
                                
                                # Check if column 9 (index 9) matches the element
                                if fields[9] == element:
                                    # Create two BED entries for arm A and arm B
                                    chrom = fields[0]
                                    armA_start = fields[1]
                                    armA_end = fields[2]
                                    armB_start = fields[5]
                                    armB_end = fields[6]
                                    pct_id = fields[7].replace('%', 'pct')
                                    
                                    # Determine which Q column this came from
                                    if col_header.startswith('Q'):
                                        hg002_q = col_header
                                    else:
                                        # For special columns like HG02071_chrY.Q9.2, simplify to HG02071Q9.2
                                        hg002_q = col_header.replace('_chrY.', '')
                                    
                                    # Map Q designations to P designations
                                    q_to_p = {
                                        'Q1': 'P8',
                                        'Q2': 'P7',
                                        'Q3': 'P6',
                                        'Q4': 'P5',
                                        'Q5': 'P4',
                                        'Q7.7': 'P3',
                                        'Q9.2': 'DAZ',
                                        'Q11.1': 'RED',
                                        'Q8': 'P2Long',
                                        'HG00621Q9.1': 'P1Long',
                                        'HG02071Q9.2': 'P1Long'
                                    }
                                    
                                    p_designation = q_to_p.get(hg002_q, '')
                                    
                                    # Build name: element_arm_pct_id_hg002_q_p_designation
                                    name_armA = f"{element}A|{pct_id}|{hg002_q}|{p_designation}" if p_designation else f"{element}A|{pct_id}|{hg002_q}"
                                    name_armB = f"{element}B|{pct_id}|{hg002_q}|{p_designation}" if p_designation else f"{element}B|{pct_id}|{hg002_q}"
                                    
                                    # Add arm A
                                    bed_entries.append(f"{chrom}\t{armA_start}\t{armA_end}\t{name_armA}\n")
                                    
                                    # Add arm B
                                    bed_entries.append(f"{chrom}\t{armB_start}\t{armB_end}\t{name_armB}\n")
                                    
                                    break
                    except Exception as e:
                        print(f"Error processing {pal_file}: {e}")
    
    # Save BED file
    bed_output = f"output_palindromes{timestamp}.bed"
    with open(bed_output, 'w') as f:
        f.writelines(bed_entries)
    
    print(f"BED file saved: {bed_output} ({len(bed_entries)} entries)")







if __name__ == "__main__":
    main()