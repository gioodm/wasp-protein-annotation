#!/usr/bin/env python3.10
# @author: Giorgia Del Missier

import argparse
import cobra
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description='Process GEM model to find orphan reactions and their annotations.')
    parser.add_argument('--input', required=True, help='Input GEM model file.')
    parser.add_argument('--extension', required=True, choices=['xml', 'sbml', 'mat', 'json'], help='GEM model file extension.')
    parser.add_argument('--output', required=True, help='Output file name (.tsv) containing orphan reactions and annotations.')
    return parser.parse_args()

def load_model(input_file, extension):
    """
    Load the GEM model from the specified file and format.
    """
    if extension in ["xml", "sbml"]:
        return cobra.io.read_sbml_model(input_file)
    elif extension == "json":
        return cobra.io.load_json_model(input_file)
    elif extension == "mat":
        return cobra.io.load_matlab_model(input_file)
    else:
        raise ValueError("File extension not accepted")

def find_orphan_rxns(model):
    """
    Find orphan reactions in the model. These are reactions with no associated genes.
    """
    return [(reaction.id, idx+1) for idx, reaction in enumerate(model.reactions) if len(reaction.genes) == 0]

def find_exchange_rxns(model, inclObjFlag=False, irrevFlag=False):
    """
    Find exchange reactions in the model. Exchange reactions typically involve only one metabolite.
    """
    S = cobra.util.create_stoichiometric_matrix(model)
    if not irrevFlag:
        # Identify exchange reactions in reversible models
        selExc = (S.sum(axis=0) == -1) | (S.sum(axis=0) == 1)
        selExc &= (S != 0).sum(axis=0) == 1
        if hasattr(model, 'objective_coefficients'):
            selExc[model.objective_coefficients != 0] = inclObjFlag
        selUpt = (model.lower_bounds < 0) & selExc if hasattr(model, 'lower_bounds') else []
    else:
        # Identify exchange reactions in irreversible models
        selExc = (abs(S).sum(axis=0) == 1) & ((S != 0).sum(axis=0) == 1)
        if hasattr(model, 'objective_coefficients'):
            selExc[model.objective_coefficients != 0] = inclObjFlag
        selUpt = (S.sum(axis=0) == 1) & ((S != 0).sum(axis=0) == 1)
    return selExc, selUpt

def get_rxns_annotation(model, orphans):
    """
    Retrieve annotations for orphan reactions, excluding non-enzymatic ones.
    """
    exc_words = ['diffusion', 'spontaneous', 'slime', 'biomass', 'pseudoreaction', 'leak', 
                 'non enzymatic', 'non enz', 'artificial', 'pseudo', 'absorption', 'non-enzymatic']
    data = []
    for rid in orphans:
        rxn = model.reactions.get_by_id(rid)
        if not any(word in rxn.name.lower() for word in exc_words):
            data.append({'rxn_id': rid, 'rxn_name': rxn.name, **rxn.annotation})
    return pd.DataFrame(data)

def main():

    args = parse_arguments()
    model = load_model(args.input, args.extension)
    
    orphans = find_orphan_rxns(model)
    indexes = [i[1] for i in orphans]
    
    selExc, _ = find_exchange_rxns(model, inclObjFlag=True, irrevFlag=False)
    non_exc_rxns = [idx for idx in indexes if selExc[idx-1]]
    orphans_ids = [r[0] for r in orphans if r[1] not in non_exc_rxns]

    print("Number of orphan reactions in the model:", len(orphans_ids))
    
    df = get_rxns_annotation(model, orphans_ids)
    df.to_csv(args.output, sep='\t', index=False)
    print(f"Orphan reactions and annotations written to {args.output}")

if __name__ == '__main__':
    main()
