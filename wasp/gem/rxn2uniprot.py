#!/usr/bin/env python3.10
# @author: Giorgia Del Missier

import requests
import time

def get_UniProt(protein, max_retries=3, delay=2):
    """
    Retrieve UniProt IDs associated with a given Rhea ID or EC number using the UniProt REST API.
    """
    UNIPROT_API = "https://rest.uniprot.org/uniprotkb/stream"

    # Determine query string
    if protein.startswith("RHEA:"):
        rhea_id = protein.split(":")[1]
        query = f'(cc_catalytic_activity_exp:"rhea:{rhea_id}") AND (database:alphafolddb)'
    elif protein.startswith("EC:"):
        ec_id = protein.split(":")[1]
        query = f'(cc_catalytic_activity_exp:"EC:{ec_id}") AND (database:alphafolddb)'
    else:
        print(f"Unrecognized protein format: {protein}")
        return []

    params = {"format": "list", "query": query}

    # Attempt with retries
    for attempt in range(1, max_retries + 1):
        try:
            response = requests.get(UNIPROT_API, params=params, timeout=10)

            if response.status_code == 200:
                uniprot_ids = response.text.strip().split()
                if uniprot_ids:
                    return uniprot_ids
                else:
                    print(f"[{protein}] No UniProt IDs found (empty response).")
                    return []
            else:
                print(f"[{protein}] Attempt {attempt}: HTTP {response.status_code} - {response.text.strip()[:200]}")

        except requests.exceptions.RequestException as e:
            print(f"[{protein}] Attempt {attempt}: Request failed - {e}")

        if attempt < max_retries:
            time.sleep(delay)

    print(f"[{protein}] All attempts failed.")
    return []

def map_orphan_rxns(input_file, output_file, output_ids_file):
    """
    Map orphan reactions to UniProt IDs using Rhea or EC numbers.
    """
    all_uniprot_ids = []

    with open(output_ids_file, "w") as f_ids, open(output_file, "w") as f_out:
        cols = ["#rxn ID", "rxn (extended name)", "rxn codes (Rhea ID/EC number)", "UniProt IDs"]
        f_out.write("\t".join(cols) + "\n")
        
        # Dictionary to cache UniProt IDs for Rhea/EC numbers
        rhea2uniprot = {}

        # Read the input file line by line
        with open(input_file) as f_orphans2rhea:
            for line in f_orphans2rhea:
                if line.startswith("#"):
                    continue

                line = line.strip().split("\t")
                orphan = line[0]
                rxn_name = line[1] if len(line) > 1 else ''
                rxn_codes = line[2:] if len(line) > 2 else []

                orphan_uniprot_ids = []
                uniprot_ids = []

                for rxn in rxn_codes:
                    rxn = rxn.strip(" , ").split(", ")
                    for r in rxn:
                        if r in rhea2uniprot:
                            uniprot_ids = rhea2uniprot[r]
                        else:
                            if r != "":
                                uniprot_ids = get_UniProt(r)
                                rhea2uniprot[r] = uniprot_ids
                        
                        orphan_uniprot_ids.extend(uniprot_ids)

                orphan_uniprot_ids = list(set(orphan_uniprot_ids))
                all_uniprot_ids.extend(orphan_uniprot_ids)

                fields = [orphan, rxn_name, ", ".join(rxn_codes), ", ".join(orphan_uniprot_ids)]
                f_out.write("\t".join(fields) + "\n")

        all_uniprot_ids = list(set(all_uniprot_ids))

        for uniprot_id in all_uniprot_ids:
            f_ids.write(f"AF-{uniprot_id}-F1-model_v4\n")