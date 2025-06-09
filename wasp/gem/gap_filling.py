#!/usr/bin/env python3.10# @author: Giorgia Del Missier
#!/usr/bin/env python3.10
# @author: Giorgia Del Missier

import numpy as np
import networkx as nx
import itertools

np.random.seed(0)

def parse_m8(infile, eval_thr, bits_thr, tm_thr):
    """
    Parses the Foldseek output files.
    """
    all_queries = dict()

    with open(infile) as fs_out:
        for line in fs_out:
            line = line.strip().split()
            query, target, evalue, bitscore, tmscore = line[0], line[1], float(line[-2]), int(line[-1]), float(line[-3])

            # Extract query and target IDs from headers
            query = query.split("-")[1] if "-" in query else query
            target = target.split("-")[1] if "-" in target else target

            if query != target:
                if bitscore > bits_thr and evalue < eval_thr and tmscore > tm_thr:
                    all_queries.setdefault(query, []).append((target, tmscore))

    return all_queries

def fill_gaps(allq, rxn2up, tms, fout):
    """
    Fills gaps in GEM models using structural alignment results.
    """
    cols = ["#rxn ID", "rxn (extended name)", "rxn codes (Rhea ID/EC number)", "UniProt IDs (other organisms)", "top hit (target organism)", "other hits (target organism)"]

    with open(fout, 'w') as f_out:
        f_out.write("\t".join(cols) + "\n")

        for key, values in rxn2up.items():
            rxn_codes = ', '.join([str(x) for x in values[0] if x])

            if not values[1]:
                fields = [key[0], key[1], '', '', '', '']
                f_out.write("\t".join(map(str, fields)) + "\n")
                continue

            up_ids = values[1]
            if len(up_ids) == 1:
                up_id = up_ids[0]
                if up_id in allq:
                    hits = sorted(allq[up_id], key=lambda x: x[1], reverse=True)
                    fields = [key[0], key[1], rxn_codes, up_id, hits[0], str(hits[1:]).strip('[]')]
                else:
                    fields = [key[0], key[1], rxn_codes, up_id, '', '']
                f_out.write("\t".join(map(str, fields)) + "\n")
                continue

            G = nx.Graph()
            for u, v in itertools.combinations(up_ids, 2):
                if (u, v) in tms:
                    G.add_edge(u, v)

            Gcc = sorted(nx.connected_components(G), key=len, reverse=True)
            if not Gcc:
                fields = [key[0], key[1], rxn_codes, ', '.join(up_ids), '', '']
                f_out.write("\t".join(map(str, fields)) + "\n")
                continue

            for component in Gcc:
                hits = [allq[j] for j in component if j in allq]
                if not hits:
                    fields = [key[0], key[1], rxn_codes, ', '.join(component), '', '']
                    f_out.write("\t".join(map(str, fields)) + "\n")
                    continue

                id_dict = {}
                for hit_list in hits:
                    for target, score in hit_list:
                        id_dict.setdefault(target, []).append(score)

                output = sorted([(k, round(sum(v) / len(v), 3)) for k, v in id_dict.items() if len(v) == len(hits)], key=lambda x: x[1], reverse=True)
                if output:
                    fields = [key[0], key[1], rxn_codes, ', '.join(component), output[0], str(output[1:]).strip('[]')]
                else:
                    fields = [key[0], key[1], rxn_codes, ', '.join(component), '', '']
                f_out.write("\t".join(map(str, fields)) + "\n")

def run_gap_filling(input_file, input_db, input_rxn2up, output_file, eval_thr, bits_thr, tm_thr):
    all_queries = parse_m8(input_file, eval_thr, bits_thr, tm_thr)

    rxn2up_dict = {}
    with open(input_rxn2up) as f_rxn2up:
        for line in f_rxn2up:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            if len(line) == 4:
                rxn_codes = [item for item in line[2].split(', ')]
                rxn_up = [item for item in line[3].split(', ')]
                rxn2up_dict[(line[0], line[1])] = (rxn_codes, rxn_up)

    tms_dict = {}
    with open(input_db) as f_db:
        for line in f_db:
            elements = line.strip().split("\t")
            query_id = elements[0].split("-")[1]
            target_id = elements[1].split("-")[1]
            tmscore = float(elements[-3])
            if tmscore > tm_thr:
                tms_dict[(query_id, target_id)] = tmscore

    fill_gaps(all_queries, rxn2up_dict, tms_dict, output_file)
