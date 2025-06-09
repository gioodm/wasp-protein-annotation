#!/usr/bin/env python3.10
# @author: Giorgia Del Missier

import os
import requests
import multiprocessing
from functools import partial
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(description="Download structures from AlphaFold DB using UniProt IDs")
    parser.add_argument('--input_file', required=True, help="Input file with UniProt IDs or a FASTA file")
    parser.add_argument('--output_dir', required=True, help="Directory where the structures will be saved")
    parser.add_argument('--num_cores', type=int, default=16, help="Number of CPU cores to use for downloading")
    
    return parser.parse_args()

def download_structure(uniprot_id, output_dir):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(url)
    
    if response.status_code == 200:
        file_path = os.path.join(output_dir, f"{uniprot_id}.pdb")
        with open(file_path, 'wb') as f:
            f.write(response.content)
    else:
        print(f"Failed to download {uniprot_id} (status code {response.status_code})")

def download_all_structures(uniprot_ids, output_dir, num_cores):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Use multiprocessing to download structures in parallel
    pool = multiprocessing.Pool(processes=num_cores)
    download_func = partial(download_structure, output_dir=output_dir)
    
    pool.map(download_func, uniprot_ids)
    pool.close()
    pool.join()

def read_uniprot_ids(input_file):
    uniprot_ids = []

    if input_file.endswith(".fasta") or input_file.endswith(".fa"):
        for record in SeqIO.parse(input_file, "fasta"):
            uniprot_id = record.id.split("|")[1] if '|' in record.id else record.id # Extracting UniProt ID from the FASTA header
            uniprot_ids.append(uniprot_id)
    else:
        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line:  # Only add non-empty lines
                    uniprot_ids.append(line)

    return uniprot_ids

def main():
    args = parse_args()

    # Read UniProt IDs either from a list of UniProt IDs or a FASTA file
    uniprot_ids = read_uniprot_ids(args.input_file)
    download_all_structures(uniprot_ids, args.output_dir, args.num_cores)

if __name__ == '__main__':
    main()
