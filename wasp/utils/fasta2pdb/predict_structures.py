import os
import subprocess
import multiprocessing
from functools import partial
import argparse

def predict_structure(uniprot_id, output_dir):
    # Create an output folder for each UniProt ID
    uniprot_output_dir = os.path.join(output_dir, uniprot_id)
    os.makedirs(uniprot_output_dir, exist_ok=True)

    # Run ColabFold for each UniProt ID
    command = f"colabfold_batch --use-gpu --output-dir {uniprot_output_dir} {uniprot_id}.fasta {uniprot_output_dir}"
    
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Prediction completed for {uniprot_id}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to predict structure for {uniprot_id}: {e}")

def predict_all_structures(uniprot_ids, output_dir, num_cores):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Use multiprocessing to predict structures in parallel
    pool = multiprocessing.Pool(processes=num_cores)
    predict_func = partial(predict_structure, output_dir=output_dir)
    
    pool.map(predict_func, uniprot_ids)
    pool.close()
    pool.join()

def parse_args():
    parser = argparse.ArgumentParser(description="Predict structures using ColabFold from UniProt IDs.")
    parser.add_argument('-i', '--input_file', required=True, help="Input file with UniProt IDs, one per line.")
    parser.add_argument('-o', '--output_dir', required=True, help="Directory where the predicted structures will be saved.")
    parser.add_argument('-c', '--num_cores', type=int, default=4, help="Number of CPU cores to use for predictions (default: 4).")
    
    return parser.parse_args()

def read_uniprot_ids(input_file):
    with open(input_file, 'r') as f:
        return [line.strip() for line in f.readlines() if line.strip()]

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_args()

    # Read UniProt IDs from input file
    uniprot_ids = read_uniprot_ids(args.input_file)

    # Predict structures using specified number of cores
    predict_all_structures(uniprot_ids, args.output_dir, args.num_cores)
