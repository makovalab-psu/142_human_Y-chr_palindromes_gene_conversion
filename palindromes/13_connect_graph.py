#!/usr/bin/env python3

import glob
import networkx as nx
import numpy as np
import pandas as pd
import os
import json

INPUT_DIR="/storage/home/kxp5629/proj/10_HPRC_R2_Y/data/palindromes_A8K_S500K/all_v_all/individual_alignments"

def parse_line(line):
    line = line.strip()
    if line.startswith("#"):
        return None
    fields = line.split("\t")

    palindrome_A = fields[0]
    sample_A = palindrome_A.split("_")[0]
    palindrome_B = fields[5]
    sample_B = palindrome_B.split("_")[0]

    # By ignoring self alignments I lose track of singletons that align only to themselves e.g. HG03209 Q1
    # if sample_A == sample_B:
    #     return None

    cov1 = int(fields[2])/int(fields[1])
    cov2 = int(fields[8])/int(fields[7])
    ident = fields[11]

    return palindrome_A, palindrome_B, cov1, cov2, ident

def main():

    G = nx.Graph()

    pattern = f"{INPUT_DIR}/*.lastz"
    files = glob.glob(pattern)

    best_matches = {}
    all_nodes = set()
    
    for file_path in files:
    
        file_name = file_path.split("/")[-1]

        if os.stat(file_path).st_size != 0:
            with open(file_path, 'r') as f:
                for line in f:
                    result = parse_line(line)
                    if result == None:
                        continue
                    sample_1, sample_2, cov1, cov2, ident = result

                    all_nodes.add(sample_1)
                    all_nodes.add(sample_2)

                    if sample_1 == sample_2:
                        continue

                    if sample_1 not in best_matches or cov1 > best_matches[sample_1][1]:
                        best_matches[sample_1] = (sample_2, cov1, cov2, ident)
                    if sample_2 not in best_matches or cov2 > best_matches[sample_2][1]:
                        best_matches[sample_2] = (sample_1, cov2, cov1, ident)



                    if cov1 > 0.8 and cov2 > 0.8:
                        G.add_edge(sample_1,sample_2)

    G.add_nodes_from(all_nodes)

    # nx.write_graphml(G, f"{INPUT_DIR}/palindrome_graph.graphml")
    singletons = [node for node in G.nodes() if G.degree(node) == 0]
    with open(f"{INPUT_DIR}/best_matches.txt", "w") as text_file:
        print("sample\tpalindrome\tbest_match\tcov1\tcov2\tidentity", file = text_file)
        for entry in singletons:
            sample = entry.split("_")[0]
            tab_split = "NA\tNA\tNA\tNA"
            try:
                tab_split = "\t".join([str(x) for x in best_matches[entry]])
            except:
                pass
                

            print(f"{sample}\t{entry}\t{tab_split}", file = text_file)

    #print individual components
    components = list(nx.connected_components(G))

    subtables_dir = f"{INPUT_DIR}/subtables"

    os.makedirs(subtables_dir, exist_ok=True)

    for i,component in enumerate(components):
        nodes = list(component)

        matrix = np.zeros((len(nodes), len(nodes)))
        for j, node1 in enumerate(nodes):
                for k, node2 in enumerate(nodes):
                    if G.has_edge(node1, node2):
                        matrix[j, k] = 1

        df = pd.DataFrame(matrix, index=nodes, columns=nodes)

        print(f"Component {i+1} ({len(nodes)} nodes):")
        print(df)
        print("\n" + "="*50 + "\n")

        # Save to CSV if desired
        df.to_csv(f"{subtables_dir}/component_{i+1}.csv")


if __name__ == "__main__":
    main()