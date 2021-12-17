import sys
import dendropy
from dendropy.calculate.treecompare \
    import false_positives_and_negatives
import numpy as np
import pandas as pd
from os import listdir
from copy import deepcopy
from multiprocessing import Pool

def Robinson_Foulds(ip):
    tr1 = ip[0]
    tr2 = ip[1]
    g1 = ip[2]
    g2 = ip[3]

    lb1 = set([l.taxon.label for l in tr1.leaf_nodes()])
    lb2 = set([l.taxon.label for l in tr2.leaf_nodes()])
    com = lb1.intersection(lb2)
    tns = dendropy.TaxonNamespace(com)
    tr1.retain_taxa_with_labels(com)
    tr1.migrate_taxon_namespace(tns)
    tr2.retain_taxa_with_labels(com)
    tr2.migrate_taxon_namespace(tns)
    tr1.update_bipartitions()
    tr2.update_bipartitions()
    [fp, fn] = false_positives_and_negatives(tr1, tr2)
    nl = len(com)
    rf = (fp+fn)/(2.0*nl-6.0)
    return {'Taxa1':g1, 'Taxa2':g2, 'RF':rf}

if __name__ == '__main__':
    filepath = sys.argv[1]
    outpath = sys.argv[2]

    files = listdir(filepath)
    files[:] = [f for f in files if 'best' in f]
    files.sort()

    Trees = {}
    for f in files:
        t = dendropy.Tree.get(path = filepath+f, schema = "newick", rooting='force-unrooted')
        Trees[f] = t
    
    inputs = []
    for t1 in files:
        for t2 in files:
            tr1 = deepcopy(Trees[t1])
            tr2 = deepcopy(Trees[t2])
            i = [tr1, tr2, t1, t2]
            inputs.append(i)
    
    with Pool(64) as P:
        out = P.map(Robinson_Foulds, inputs)

    df_op = pd.DataFrame(out)
    df_op.to_csv(outpath, sep = '\t')