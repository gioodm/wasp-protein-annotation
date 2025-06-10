#!/usr/bin/env python3.10
# @author Giorgia Del Missier

#!/usr/bin/env python3.10
# Author: Giorgia Del Missier

import os
import sys
import multiprocessing
import pandas as pd
import numpy as np
import altair as alt
from altair_saver import save
from pathlib import Path

# Get the directory where safe.py is located
safepy_dir = str(Path(__file__).parent.resolve())
sys.path.append(safepy_dir)

from wasp.main.safepy import safe

# Set random seed for reproducibility
np.random.seed(0)

def process_identifier(args_list):
    """Process protein identifiers for SAFE analysis.
    
    Args:
        args_list (tuple): A tuple containing the identifier and the dataframe.
    """
    i, df = args_list
    subdf = df[['#Cluster', 'UniProt ID', i]]
    flattened_ids = subdf[i].dropna().astype(str).str.split(';').explode()
    final_matrix = pd.DataFrame()
    
    if not flattened_ids.empty:
        unique_ids = flattened_ids.dropna().unique()
        final_matrix = pd.DataFrame(columns=unique_ids)

        for index, row in subdf.iterrows():
            protein = row['UniProt ID']
            ids = row[i] if not pd.isna(row[i]) else []
            final_matrix.loc[protein] = [1 if j in ids else 0 for j in unique_ids]

    output_file = f"{wd}/SAFE/{i}_matrix.txt"
    final_matrix.to_csv(output_file, sep='\t')

def perform_safe_analysis(i, wd, radius, sf, enrichment_thr=0.05):
    """Perform SAFE analysis on the given identifier.
    
    Args:
        i (str): Identifier for the analysis.
        wd (str): Working directory for output files.
        radius (float): Neighborhood radius for SAFE analysis.
        sf (safe.Safe): SAFE analysis instance.
        enrichment_thr (float): Threshold for enrichment significance.
        
    Returns:
        tuple: The identifier and its significant results.
    """
    sf.load_attributes(attribute_file=f'SAFE/{i}_matrix.txt')
    sf.define_neighborhoods(neighborhood_radius=float(radius))
    sf.compute_pvalues(background='network', num_permutations=1000, processes=64)
    sf.print_output_files(output_dir=f'{wd}/SAFE/{i}_')

    df_tmp = pd.read_csv(f"{wd}/SAFE/{i}_node_properties_annotation.txt", sep="\t")
    
    significant = {}
    for index, row in df_tmp.iterrows():
        key = row['key']
        significant[key] = [
            (col, round(row[col], 3)) 
            for col in df_tmp.columns[3:] 
            if row[col] > -np.log10(enrichment_thr) and not col.startswith("Unnamed")
        ]
    
    return i, significant

def prepare_df(d):
    """Prepare a DataFrame from the provided data dictionary.
    
    Args:
        d (dict): A dictionary where keys are identifiers and values are dataframes.
        
    Returns:
        pd.DataFrame: A stacked DataFrame suitable for plotting.
    """
    df = pd.DataFrame()

    for i in d:
        tmp_df = pd.DataFrame(d[i]).transpose()
        tmp_df = tmp_df.stack().reset_index()
        tmp_df.columns = ['id', 'group', 'count']
        tmp_df['category'] = i

        df = pd.concat([df, tmp_df])

    return df

def make_barcharts(txid, coll, outfig, iteration):

    df_taxid = prepare_df(txid)
    df_collateral = prepare_df(coll)

    chart_taxid = alt.Chart(df_taxid).mark_bar().encode(

        # tell Altair which field to group columns on
        x=alt.X('group:N', title=None, sort=None),  # Add sort argument

        # tell Altair which field to use as Y values and how to calculate
        y=alt.Y('sum(count):Q',
            axis=alt.Axis(
                grid=False,
                title=None)),

        # tell Altair which field to use to use as the set of columns to be represented in each group
        column=alt.Column('id:N', sort=alt.SortField(field='id', order='ascending'), title=None),  # Add sort argument

        # tell Altair which field to use for color segmentation 
        color=alt.Color('category:N', sort=None,
            scale=alt.Scale(
                # make it look pretty with an enjoyable color pallet
                range=['#001219','#005F73','#94D2BD','#EE9B00','#BB3E03','#9B2226'],
            )),

        order=alt.Order(
        # Sort the segments of the bars by this field
        'category',
        sort='ascending')
            )\
            .configure_view(
            # remove grid lines around column clusters
            strokeOpacity=0    
    ).properties(
    width=150)

    chart_taxid = chart_taxid.configure_legend(labelFontSize=8)
    chart_taxid.save(f"{outfig}_taxid_iter{iteration}.png", ppi=200)

    chart_collateral = alt.Chart(df_collateral).mark_bar().encode(

        # tell Altair which field to group columns on
        x=alt.X('group:N', title=None, sort=None),  # Add sort argument

        # tell Altair which field to use as Y values and how to calculate
        y=alt.Y('sum(count):Q',
            axis=alt.Axis(
                grid=False,
                title=None)),

        # tell Altair which field to use to use as the set of columns to be represented in each group
        column=alt.Column('id:N', sort=alt.SortField(field='id', order='ascending'), title=None),  # Add sort argument

        # tell Altair which field to use for color segmentation 
        color=alt.Color('category:N', sort=None,
            scale=alt.Scale(
                # make it look pretty with an enjoyable color pallet
                range=['#001219','#005F73','#94D2BD','#EE9B00','#BB3E03','#9B2226'],
            )),

        order=alt.Order(
        # Sort the segments of the bars by this field
        'category',
        sort='ascending')
            )\
            .configure_view(
            # remove grid lines around column clusters
            strokeOpacity=0    
    ).properties(
    width=150)

    chart_collateral = chart_collateral.configure_legend(labelFontSize=8)
    chart_collateral.save(f"{outfig}_collateral_iter{iteration}.png", ppi=200)

def analyze_safe_results(df_anno, df_safe, ids, allp, output, wd, iteration, nan):

    taxid, collateral = {}, {}

    for key in ["total", "previously annotated", "not annotated", "nan2nan", "new annotation by WASP", "nan2nan in dark clusters"]:
        taxid[key], collateral[key] = {}, {}

        for identifier in ids:
            taxid[key][f'{ids.index(identifier) + 1}.{identifier}'] = {'total': 0, 'annotation': 0, 'nan2nan / new': 0, 'nan2nan (dark)': 0}
            collateral[key][f'{ids.index(identifier) + 1}.{identifier}'] = {'total': 0, 'annotation': 0, 'nan2nan / new': 0, 'nan2nan (dark)': 0}

    list_nan2nan = list()

    cols = ["UniProt ID"] + ids + ["Organism"]
    taxid_outdf, collateral_outdf = pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

    nclusters = df_anno["#Cluster"].unique()
    for c in nclusters:
        
        subdf_anno = df_anno[df_anno["#Cluster"] == c]
        subdf_safe = df_safe[df_safe["UniProt ID"].isin(subdf_anno["UniProt ID"])]

        dark_anno = subdf_anno[ids].apply(lambda col: col.isna().all()).to_dict()
        dark_safe = subdf_safe[ids].apply(lambda col: col.isna().all()).to_dict()
        
        tmp_taxid, tmp_collateral = pd.DataFrame(columns=cols), pd.DataFrame(columns=cols)

        for i in ids:

            taxid['total'][f'{ids.index(i) + 1}.{i}']['total'] += len(subdf_anno[subdf_anno['UniProt ID'].isin(allp['UniProt ID'])])
            collateral['total'][f'{ids.index(i) + 1}.{i}']['total'] += len(subdf_anno[~subdf_anno['UniProt ID'].isin(allp['UniProt ID'])])

            if dark_anno[i] and dark_safe[i]:
                proteins = subdf_anno[["UniProt ID", "Organism"]]
                proteins_taxid = proteins[proteins['UniProt ID'].isin(allp['UniProt ID'])]
                list_nan2nan.extend(proteins_taxid["UniProt ID"].tolist())

                taxid['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(proteins_taxid)
                taxid['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += len(proteins_taxid)
                taxid['nan2nan in dark clusters'][f'{ids.index(i) + 1}.{i}']['nan2nan (dark)'] += len(proteins_taxid)

                collateral['nan2nan'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += len(proteins[~proteins['UniProt ID'].isin(allp['UniProt ID'])])
                collateral['not annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += len(proteins[~proteins['UniProt ID'].isin(allp['UniProt ID'])])
                collateral['nan2nan in dark clusters'][f'{ids.index(i) + 1}.{i}']['nan2nan (dark)'] += len(proteins[~proteins['UniProt ID'].isin(allp['UniProt ID'])])
            else:
                for index, row in subdf_safe.iterrows():

                    if row[i] == row[i]: 
                        if pd.isna(subdf_anno[subdf_anno["UniProt ID"] == row["UniProt ID"]][i].values[0]):
                            if row['UniProt ID'] in allp['UniProt ID'].tolist():
                                taxid['new annotation by WASP'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += 1
                            else:
                                collateral['new annotation by WASP'][f'{ids.index(i) + 1}.{i}']['nan2nan / new'] += 1
                        else:
                            if row['UniProt ID'] in allp['UniProt ID'].tolist():
                                taxid['previously annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += 1
                            else:
                                collateral['previously annotated'][f'{ids.index(i) + 1}.{i}']['annotation'] += 1

                    if row['UniProt ID'] in allp['UniProt ID'].tolist():
                        tmp_taxid = pd.concat([tmp_taxid, subdf_safe[subdf_safe["UniProt ID"] == row["UniProt ID"]]])
                    else:
                        tmp_collateral = pd.concat([tmp_collateral, subdf_safe[subdf_safe["UniProt ID"] == row["UniProt ID"]]])

        taxid_outdf = pd.concat([taxid_outdf, tmp_taxid])
        collateral_outdf = pd.concat([collateral_outdf, tmp_collateral])

    taxid_outdf.to_csv(f"{wd}/taxid_{output}_iter{iteration}.txt", sep="\t", index=False)
    collateral_outdf.to_csv(f"{wd}/collateral_{output}_iter{iteration}.txt", sep="\t", index=False)
    
    if len(list_nan2nan) > 0:
        nan_list_df = pd.DataFrame(list_nan2nan, columns=["nan2nan_list"])
        nan_list_df.to_csv(f"{wd}/{nan}_iter{iteration}.txt", sep="\t", index=False)

    return taxid, collateral

def run_safe_analysis(input_file, wd, edgelist, iteration, identifiers, output, outfig, nan_file):
    """Run the SAFE analysis pipeline with the specified arguments.
    """

    if not os.path.exists(f"{wd}/SAFE"):
        os.makedirs(f"{wd}/SAFE")

    df_anno = pd.read_csv(input_file, sep="\t", header=0)
    df_ids = pd.read_csv(identifiers, sep="\t", header=0)
    all_proteins = df_ids[['UniProt ID']]

    safe_file = safe.Safe(
        input_network=edgelist, 
        input_node_attribute=None, 
        neighborhood_radius=float(0.5), 
        output_dir=wd, 
        processes=64, 
        output_prefix=None
    )

    ids = df_anno.columns[2:]

    pool = multiprocessing.Pool(processes=64)
    pool.map(process_identifier, [(i, df_anno) for i in ids])

    manager = multiprocessing.Manager()
    results = manager.dict()

    results_list = pool.starmap(perform_safe_analysis, [(i, wd, 0.5, safe_file) for i in ids])
    
    for res in results_list:
        key, value = res
        results[key] = value

    taxid, collateral = analyze_safe_results(df_anno, df_anno, ids, all_proteins, output, wd, iteration, nan_file)
    make_barcharts(taxid, collateral, outfig, iteration)

