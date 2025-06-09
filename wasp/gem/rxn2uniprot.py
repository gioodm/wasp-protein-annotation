#!/usr/bin/env python3.10
# @author: Giorgia Del Missier

import argparse
import requests

def get_UniProt(protein):
    """
    Retrieve UniProt IDs associated with a given Rhea ID or EC number using the UniProt REST API.
    """
    UNIPROT_API = "https://rest.uniprot.org/uniprotkb/"
    response = None

    try:

        if protein.startswith("RHEA"):
            rhea_id = protein.split(":")[1]
            # Query UniProt API for UniProt IDs associated with the Rhea ID
            response = requests.get(f"{UNIPROT_API}stream?format=list&query=(cc_catalytic_activity_exp:\"rhea:{rhea_id}\") AND (database:alphafolddb)")

        elif protein.startswith("EC"):
            ec_id = protein.split(":")[1]
            # Query UniProt API for UniProt IDs associated with the EC number
            response = requests.get(f"{UNIPROT_API}stream?format=list&query=(cc_catalytic_activity_exp:\"EC:{ec_id}\") AND (database:alphafolddb)")

        # Extract UniProt IDs from the response
        uniprot_ids = response.text.strip().split()
        return uniprot_ids
    
    except Exception as e:
        print(f"Error retrieving UniProt IDs for {protein}: {e}")
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