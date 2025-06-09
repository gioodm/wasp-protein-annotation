#!/usr/bin/env python3.10
# @author: Giorgia Del Missier

import requests
import re
from requests.adapters import HTTPAdapter, Retry

# Regular expression to find the next link in the headers
re_next_link = re.compile(r'<(.+)>; rel="next"')

# Setting up retries for HTTP requests
retries = Retry(total=25, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    """
    Extracts the next link from the headers if available.

    Parameters:
    - headers: HTTP response headers

    Returns:
    - Next link URL if available, otherwise None
    """
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_batch(batch_url):
    """
    Generator to handle paginated UniProt API responses.

    Parameters:
    - batch_url: Initial URL to fetch the batch

    Yields:
    - Response and total results count
    """
    while batch_url:
        response = session.get(batch_url, timeout=25)
        response.raise_for_status()
        total = response.headers.get("x-total-results", "Unknown")
        yield response, total
        batch_url = get_next_link(response.headers)

def get_UniProt(chunk):
    """
    Fetches UniProt annotations for a chunk of protein IDs.

    Parameters:
    - chunk: List of UniProt IDs

    Returns:
    - Dictionary of UniProt IDs and their annotations
    """
    UNIPROT_API = "https://rest.uniprot.org/"
    pnames = "%20OR%20".join(chunk)

    annotation = {}
    url = f"{UNIPROT_API}/uniprotkb/search?query=accession={pnames}&fields=accession,xref_pfam,xref_panther,xref_gene3d,ec,rhea,go,length,annotation_score,organism_name&format=tsv&size=500"
    
    for batch, total in get_batch(url):
        for r in batch.text.splitlines()[1:]:
            r = r.split('\t')
            r = [item.strip(";") for item in r]
            r[4] = r[4].replace("; ", ";")
            r[5] = r[5].replace(" ", ";")
            r[6] = ';'.join(['GO:' + g for g in re.findall(r'GO:(\d+)', r[6])])
            if ':' in r[2]:
                ids = r[2].split(";")
                r[2] = ';'.join(set([pthr for pthr in ids if ':' in pthr]))
            annotation[r[0]] = r[1:]

    return annotation

def fetch_annotations(input_file, output_file):
    """
    Fetch UniProt annotations based on the provided cluster information.

    Parameters:
    - input_file: Path to the input file containing cluster data
    - output_file: Path to the output file for saving annotations
    """
    cols = ["#Cluster", "UniProt ID", "Pfam", "PANTHER", "CATH", "EC number", "Rhea ID", "GO terms", "Length", "Annotation score", "Organism"]

    with open(input_file) as fin, open(output_file, "w") as fout:
        fout.write("\t".join(cols) + "\n")
        lines = fin.readlines()

        for line in lines[1:]:
            line = line.strip().split("\t")
            cluster_id = line[0]
            nodes = line[-1].split(",")
            nodes = [node.strip(" ' ' ") for node in nodes]

            # Filtering valid UniProt IDs
            up_nodes = [n for n in nodes if re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$', n) and len(n) in (6, 10)]
            chunksize = 150
            chunks = [up_nodes[x:x+chunksize] for x in range(0, len(up_nodes), chunksize)]

            # Fetch annotations for each chunk
            results = {}
            for chunk in chunks:
                results.update(get_UniProt(chunk))

            # Write annotations to the output file
            for node in nodes:
                if node in results:
                    fields = [cluster_id, node] + results[node]
                    fout.write("\t".join(str(item) for item in fields) + "\n")
                else:
                    fields = [cluster_id, node] + [""] * 9
                    fout.write("\t".join(str(item) for item in fields) + "\n")