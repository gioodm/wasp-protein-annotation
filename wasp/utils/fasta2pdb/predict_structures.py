import os
import subprocess
import multiprocessing
from functools import partial
import argparse
import requests

def download_fasta(uniprot_id, fasta_dir):
    """Download FASTA file from UniProt for a given ID."""
    fasta_path = os.path.join(fasta_dir, f"{uniprot_id}.fasta")
    if os.path.exists(fasta_path):
        return fasta_path
    
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(fasta_path, 'w') as f:
            f.write(response.text)
        return fasta_path
    else:
        raise Exception(f"Failed to download FASTA for {uniprot_id}")

def predict_structure(uniprot_id, output_dir, fasta_dir):
    """Predict structure for a single UniProt ID."""
    try:
        # Create output directory
        uniprot_output_dir = os.path.join(output_dir, uniprot_id)
        os.makedirs(uniprot_output_dir, exist_ok=True)
        
        # Download FASTA file
        fasta_path = download_fasta(uniprot_id, fasta_dir)
        
        # Run ColabFold
        command = [
            "colabfold_batch",
            "--use-gpu",
            "--num-recycle", "3",
            "--output-dir", uniprot_output_dir,
            fasta_path
        ]
        
        subprocess.run(command, check=True)
        print(f"Prediction completed for {uniprot_id}")
    except Exception as e:
        print(f"Failed to predict structure for {uniprot_id}: {str(e)}")

def predict_all_structures(uniprot_ids, output_dir, fasta_dir, num_cores):
    """Predict structures for all UniProt IDs in parallel."""
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(fasta_dir, exist_ok=True)
    
    pool = multiprocessing.Pool(processes=num_cores)
    predict_func = partial(predict_structure, 
                         output_dir=output_dir,
                         fasta_dir=fasta_dir)
    
    pool.map(predict_func, uniprot_ids)
    pool.close()
    pool.join()

def parse_args():
    parser = argparse.ArgumentParser(description="Predict structures using ColabFold from UniProt IDs.")
    parser.add_argument('-i', '--input_file', required=True, 
                       help="Input file with UniProt IDs, one per line.")
    parser.add_argument('-o', '--output_dir', required=True, 
                       help="Directory where predicted structures will be saved.")
    parser.add_argument('-f', '--fasta_dir', default="fasta_files",
                       help="Directory to store downloaded FASTA files.")
    parser.add_argument('-c', '--num_cores', type=int, default=4,
                       help="Number of CPU cores to use (default: 4).")
    return parser.parse_args()

def read_uniprot_ids(input_file):
    with open(input_file, 'r') as f:
        return [line.strip() for line in f.readlines() if line.strip()]

if __name__ == "__main__":
    args = parse_args()
    uniprot_ids = read_uniprot_ids(args.input_file)
    predict_all_structures(uniprot_ids, args.output_dir, args.fasta_dir, args.num_cores)
