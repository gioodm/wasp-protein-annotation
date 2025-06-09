#! /usr/bin/env python
"""This file contains the code for the SAFE class and command-line access."""
import re
import os
from pathlib import Path
import logging

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scipy.io as spio
import pandas as pd
import zipfile
import random
import shutil
import ast

from os.path import expanduser
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist, squareform
from scipy.stats import gaussian_kde
from scipy.optimize import fmin
from collections import Counter
from xml.dom import minidom


def load_network_from_txt(filename, layout='spring_embedded', node_key_attribute='key', verbose=True):
    """
    Loads network from tab-delimited text file and applies a network layout.

    :param str filename: Path to the network file
    :param str layout: Name of the network layout to be applied
    :param str node_key_attribute: Name of the node attribute to be used as key for matching with additional attributes
    :param bool verbose:
    :return:

    Notes:
        1. Header (column names) should not be present in the file with '.txt' extension.
    """

    filename = re.sub('~', expanduser('~'), filename)
    
    # Get the number of columns first
    if not Path(filename).suffix=='.gz':
        with open(filename, 'r') as f:
            first_line = f.readline()
            num_cols = len(first_line.split('\t'))
    else:
        ## extract the file
        import gzip
        with gzip.open(filename,'rt') as f:
            first_line = f.readline()
            num_cols = len(first_line.split('\t'))
            
    ## parameters provided to `pd.read_table`
    if Path(filename).suffixes[0]=='.txt':
        kws_read_table=dict(header=None,
                            # names=None,
                           )
    elif Path(filename).suffixes[0]=='.tsv':
        kws_read_table=dict(
                            header=0, # column names in the first line
                            names=range(num_cols),
                           )    
    else:
        raise ValueError(f'extension {Path(filename).suffixes[0]} not supported')
            
    if num_cols == 3:

        data = pd.read_table(filename, sep='\t', 
                             dtype={0: str, 1: str, 2: float},
                             **kws_read_table,
                            )        
        data = data.rename(columns={0: 'node_key1', 1: 'node_key2', 2: 'edge_weight'})
        data['node_label1'] = data['node_key1']
        data['node_label2'] = data['node_key2']

    elif num_cols == 5:

        data = pd.read_table(filename, sep='\t', 
                             # dtype={0: str, 1: str, 2: str, 3: str, 4: float},
                             **kws_read_table,
                            )
        
        data = data.rename(
            columns={0: 'node_label1', 1: 'node_key1', 2: 'node_label2', 3: 'node_key2', 4: 'edge_weight'})

    else:

        raise ValueError('Unknown network file format. 3 or 5 columns are expected.')

    # Merge nodes1 and nodes2 and drop duplicates
    t1 = data[['node_label1', 'node_key1']]
    t2 = data[['node_label2', 'node_key2']].rename(columns={'node_label2': 'node_label1', 'node_key2': 'node_key1'})
    nodes = pd.concat([t1, t2], ignore_index=True).drop_duplicates()

    # Re-number the node index
    nodes = nodes.reset_index(drop=True)

    # Add the node index to network data
    nodes = nodes.reset_index().set_index('node_label1')
    data['node_index1'] = nodes.loc[data['node_label1'], 'index'].values
    data['node_index2'] = nodes.loc[data['node_label2'], 'index'].values

    # Create the graph
    G = nx.Graph()

    # Add the nodes & their attributes
    nodes = nodes.reset_index().set_index('index')
    G.add_nodes_from(nodes.index.values)

    for n in G:
        G.nodes[n]['label'] = nodes.loc[n, 'node_label1']
        G.nodes[n][node_key_attribute] = nodes.loc[n, 'node_key1']

    # Add the edges between the nodes
    edges = [tuple(x) for x in data[['node_index1', 'node_index2', 'edge_weight']].values]
    for t in edges:
        n = (int(t[0]), int(t[1]))
        G.add_edge(n[0], n[1], weight=t[2])

    G = apply_network_layout(G, layout=layout, weight="weight")
    G = calculate_edge_lengths(G, verbose=verbose)

    return G


def apply_network_layout(G, layout='kamada_kawai', verbose=True, weight="weight"):

    if layout == 'kamada_kawai':

        if verbose:
            logging.info('Applying the Kamada-Kawai network layout... (may take several minutes)')

        pos = nx.kamada_kawai_layout(G, weight="weight")

    elif layout == 'spring_embedded':

        if verbose:
            logging.info('Applying the spring-embedded network layout... (may take several minutes)')

        pos = nx.spring_layout(G, k=0.2, iterations=100, weight="weight", seed=0)

    for n in G:
        G.nodes[n]['x'] = pos[n][0]
        G.nodes[n]['y'] = pos[n][1]

    return G


def calculate_edge_lengths(G, verbose=True):

    # Calculate the lengths of the edges

    if verbose:
        logging.info('Calculating edge lengths...')

    x = np.matrix(G.nodes.data('x'))[:, 1]
    y = np.matrix(G.nodes.data('y'))[:, 1]

    node_coordinates = np.concatenate([x, y], axis=1)
    node_distances = squareform(pdist(node_coordinates, 'euclidean'))

    adjacency_matrix = np.array(nx.adjacency_matrix(G).todense())
    adjacency_matrix = adjacency_matrix.astype('float')
    adjacency_matrix[adjacency_matrix == 0] = np.nan

    edge_lengths = np.multiply(node_distances, adjacency_matrix)

    edge_attr_dict = {index: v for index, v in np.ndenumerate(edge_lengths) if ~np.isnan(v)}
    nx.set_edge_attributes(G, edge_attr_dict, 'length')

    return G


def read_attributes(attribute_file='', node_label_order=None, mask_duplicates=False, fill_value=np.nan, verbose=True):

    node2attribute = pd.DataFrame()
    attributes = pd.DataFrame()

    if isinstance(attribute_file, str):

        file_name = re.sub('~', expanduser('~'), attribute_file)
        [_, file_extension] = os.path.splitext(file_name)

        if file_extension == '.mat':
            mat = load_mat(file_name)

            node2attribute = pd.DataFrame(data=np.transpose(mat['go']['term2orf']),
                                          index=mat['go']['orfs'],
                                          columns=mat['go']['term_ids'])
            node2attribute = node2attribute.apply(pd.to_numeric, downcast='unsigned')

            data = {'id': mat['go']['term_ids'], 'name': mat['go']['term_names']}
            attributes = pd.DataFrame(data=data)

        elif (file_extension == '.txt') or (file_extension == '.gz'):

            node2attribute = pd.read_csv(file_name, sep='\t', dtype={0: str})
            node2attribute.set_index(node2attribute.columns[0], drop=True, inplace=True)
            node2attribute = node2attribute.apply(pd.to_numeric, downcast='float', errors='coerce')

            data = {'id': np.arange(len(node2attribute.columns)), 'name': node2attribute.columns}
            attributes = pd.DataFrame(data=data)

            node2attribute.columns = np.arange(len(node2attribute.columns))

    elif isinstance(attribute_file, pd.DataFrame):

        node2attribute = attribute_file
        data = {'id': np.arange(len(node2attribute.columns)), 'name': node2attribute.columns}
        attributes = pd.DataFrame(data=data)

    # Force all values to numeric
    node2attribute = node2attribute.apply(pd.to_numeric, errors='coerce')

    # Force attribute names to be strings
    attributes['name'] = attributes['name'].astype(str)

    # Averaging duplicate rows (with notification)
    if not node2attribute.index.is_unique:
        logging.info('\nThe attribute file contains multiple values for the same labels. Their values will be averaged.')
        node2attribute = node2attribute.groupby(node2attribute.index, axis=0).mean()
        
    if not node_label_order:
        node_label_order = node2attribute.index.values

    node_label_in_file = node2attribute.index.values
    node_label_not_mapped = [x for x in node_label_in_file if x not in node_label_order]

    node2attribute = node2attribute.reindex(index=node_label_order, fill_value=fill_value)

    # If necessary, de-duplicate the network nodes (leave one node per gene)
    if mask_duplicates:

        # Keep a random node every time
        idx = np.random.permutation(np.arange(len(node2attribute)))
        mask_dups = node2attribute.iloc[idx].index.duplicated(keep='first')

        num_dups = mask_dups.sum()
        logging.info('\nThe network contains %d nodes with duplicate labels. '
              'Only one random node per label will be considered. '
              'The attribute values of all other nodes will be set to NaN.' % num_dups)
        node2attribute.iloc[idx[mask_dups], :] = np.nan

    node2attribute = node2attribute.values

    if verbose:
        logging.info('\nAttribute data provided: %d labels x %d attributes' % (len(node_label_in_file), attributes.shape[0]))

        # Notification about labels **not** mapped onto the network
        n = np.min([len(node_label_not_mapped), 3])
        m = len(node_label_not_mapped) - n
        if n > 0:
            msg1 = ', '.join(node_label_not_mapped[:n])
            msg2 = format(' and %d other labels in the attribute file were not found in the network.' % m)
            logging.info(msg1 + msg2)

        n_nlm = len(node_label_in_file) - len(node_label_not_mapped)
        logging.info('\nAttribute data mapped onto the network: %d labels x %d attributes' % (n_nlm, attributes.shape[0]))
        logging.info('Values: %d NaNs' % np.sum(np.isnan(node2attribute)))
        logging.info('Values: %d zeros' % np.sum(node2attribute[~np.isnan(node2attribute)] == 0))
        logging.info('Values: %d positives' % np.sum(node2attribute[~np.isnan(node2attribute)] > 0))
        logging.info('Values: %d negatives' % np.sum(node2attribute[~np.isnan(node2attribute)] < 0))

    return attributes, node_label_order, node2attribute


def plot_network(
    G,
    ax=None,
    foreground_color = '#ffffff',
    background_color='#000000',
    random_sampling_edges_min=1000000,
    title='Network',
    ):
    """
    Plot/draw a network.
    
    Note: 
        The default attribute names 
        gene ids: label_orf
        gene symbols: label
    """
    if background_color == '#ffffff':
        foreground_color = '#000000'

    node_xy = get_node_coordinates(G)

    if ax is None:
        fig, ax = plt.subplots(figsize=(20, 10), facecolor=background_color, edgecolor=foreground_color)
        fig.set_facecolor(background_color)

    # Randomly sample a fraction of the edges (when network is too big)
    edges = tuple(G.edges())
    if len(edges) >= random_sampling_edges_min:
        logging.warning(f"Edges are randomly sampled because the network (edges={len(edges)}) is too big (random_sampling_edges_min={random_sampling_edges_min}).")
        edges = random.sample(edges, int(len(edges)*0.1))

    nx.draw(G, ax=ax, pos=node_xy, edgelist=edges,
            node_color=foreground_color, edge_color=foreground_color, node_size=10, width=1, alpha=0.2)

    ax.set_aspect('equal')
    ax.set_facecolor(background_color)

    ax.grid(False)
    ax.invert_yaxis()
    ax.margins(0.1, 0.1)

    ax.set_title(title, color=foreground_color)

    plt.axis('off')

    try:
        fig.set_facecolor(background_color)
    except NameError:
        pass

    return ax


def mark_nodes(
    x,
    y,
    kind: list,
    ax=None,
    foreground_color='#ffffff',
    background_color='#000000',
    labels=None, # subset the nodes by labels
    label_va='center',
    legend_label: str=None,
    test=False,
    **kws,
    ):
    """
    Show nodes.

    Parameters:
        s (str): legend name (defaults to '').
        kind (str): 'mark' if the nodes should be marked, 'label' if nodes should be marked and labeled.
    """
    if ax is None:
        ax = plt.gca() # get current axes i.e. subplot
    if isinstance(kind, str):
        kind = [kind]
        
    if 'mark' in kind:
        ## mark the selected nodes with the marker +
        sn1 = ax.scatter(x, y, **kws)

    if 'label' in kind:
        ## show labels e.g. gene names
        if test:print(x,y,labels)
        assert len(x)==len(labels), f"len(x)!=len(labels): {len(x)}!={len(labels)}"
        if test:ax.plot(x, y, 'r*')
        for i,label in enumerate(labels):
            ax.text(x[i], y[i], label,
                    fontdict={'color': 'white' if background_color== '#000000' else 'k', 'size': 14, 'weight': 'bold'},
                    # bbox={'facecolor': 'black', 'alpha': 0.5, 'pad': 3},
                    ha='center',
                    va=label_va,
                   )

    if not legend_label is None:
        # Legend
        leg = ax.legend([sn1], [legend_label], loc='upper left', bbox_to_anchor=(0, 1),
                        title='Significance', scatterpoints=1, fancybox=False,
                        facecolor=background_color, edgecolor=background_color)

        for leg_txt in leg.get_texts():
            leg_txt.set_color(foreground_color)

        leg_title = leg.get_title()
        leg_title.set_color(foreground_color)
        
    return ax


def get_node_coordinates(graph,labels=[]):

    x = dict(graph.nodes.data('x'))
    y = dict(graph.nodes.data('y'))    

    ds = [x, y]
    pos = {}
    for k in x:
        pos[k] = np.array([d[k] for d in ds])

    node_xy_list=list(pos.values())

    if len(labels)==0:
        return  np.vstack(node_xy_list)    
    else:

        # Get the co-ordinates of the nodes
        node_labels = nx.get_node_attributes(graph, 'label')
        node_labels_dict = {k: v for v, k in node_labels.items()}
        
        # TODOs: avoid determining the x and y again.
        x = list(dict(graph.nodes.data('x')).values())
        y = list(dict(graph.nodes.data('y')).values())

        # x_offset = (np.nanmax(x) - np.nanmin(x))*0.01

        idx = [node_labels_dict[x] for x in labels if x in node_labels_dict.keys()]

        # Labels found in the data
        labels_found = [x for x in labels if x in node_labels_dict.keys()]
        x_idx = [x[i] for i in idx]
        y_idx = [y[i] for i in idx]

        # Print out labels not found
        labels_missing = [x for x in labels if x not in node_labels_dict.keys()]
        if labels_missing:
            labels_missing_str = ', '.join(labels_missing)
            logging.warning('These labels are missing from the network (case sensitive): %s' % labels_missing_str)    
        
        node_xy_list=[x_idx,y_idx]
        
        return np.vstack(node_xy_list).T, labels_found

