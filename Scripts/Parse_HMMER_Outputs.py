#!/usr/bin/env python
# coding: utf-8

# In[22]:


import pandas as pd
from sys import exit
from os.path import isdir, isfile
from os import listdir, mkdir
import Bio.SeqIO as IO
import numpy as np
import sys

# In[23]:


def Parse_HMM_Alignments(align):
    OP = []
    counts = 0
    for a in align.split("//\n"):
        o = a.split("\n")
        align_flag = False
        ctr = 0
        while ctr < len(o):
            if o[ctr].lstrip().startswith("Query:"):
                query = o[ctr].lstrip().replace("Query:","").lstrip().split(' ')[0]
                ctr += 1
            elif o[ctr].lstrip().startswith('>>'):
                print(o[ctr])
                subject = o[ctr].lstrip().replace(">>","").rstrip()
                aligns = []
                if "[No individual domains that satisfy reporting thresholds (although complete target did)]" in o[ctr+1]:
                    ctr += 2
                    continue
                ctr += 2
                while ("Alignments for each domain" not in o[ctr]) and ("Internal pipeline statistics" not in o[ctr]):
                    if o[ctr].lstrip().startswith("--"):
                        ctr += 1
                        continue
                    values = o[ctr].lstrip().replace("[.","").replace(".]","").replace("..","").replace("[]","").split(' ')
                    values = [t for t in values if t != '']
                    ctr += 1
                    header = ['No', 'Score','bias','c-Evalue', 'i-Evalue','hmmfrom','hmmto',
                              'alifrom','alito','envfrom', 'envto', 'acc']
                    print(values)
                    values = [float(v) for v in values if v != '!' and v != '?']
                    d = dict(zip(header, values))
                    d['Query'] = query
                    d['Subject'] = subject
                    d['Subject_Alignment'] = ''
                    d['Alignment_Description'] = ''
                    d['Query_Alignment'] = ''
                    aligns.append(d)
                
                if o[ctr].lstrip().startswith('Alignments for each domain') :
                    cur = 0
                    align_flag = True
                    ctr += 2
            
            if(align_flag):
                if (o[ctr].lstrip().startswith("==")):
                    cur += 1
                    ctr += 1
                tgt = o[ctr].lstrip().split(' ')
                tgt = [t for t in tgt if t != '']
                mat = o[ctr+1].lstrip().replace(" ","S")
                qry = o[ctr+2].lstrip().split(' ')
                qry = [t for t in qry if t != '']
                aligns[cur]['Subject_Alignment'] += tgt[2]
                aligns[cur]['Query_Alignment']  += qry[2]
                aligns[cur]['Alignment_Description']  += 'S'*(len(mat)-len(tgt[2]))+mat
                ctr += 4

                if o[ctr].startswith('Internal') or o[ctr].lstrip().startswith('>>') :
                    counts += 1
                    OP += aligns
                    align_flag = False
                if o[ctr].startswith('Internal'):
                    break
            else:
                ctr += 1
            
            
    df_OP = pd.DataFrame(OP)
    return df_OP
    
def Calculate_Untagged_Gene_Lengths(filepath):
    d = {}
    record_dict = IO.to_dict(IO.parse(filepath, "fasta"))
    for key in record_dict.items():
        d[key[0]] = len(key[1].seq)
    df_lengths = pd.DataFrame(data = {'Query':list(d.keys()), 'Length':list(d.values())})
    return df_lengths

def Load_HMMER_Outputs(filepath, df_lengths):
    df = pd.read_csv(filepath)
    del df['Unnamed: 0']
   
    return df

def Return_Breadth_of_Coverage(grp):
    grp = grp.sort_values(by = 'alifrom')
    starts = grp['alifrom'].tolist()
    ends = grp['alito'].tolist()
    breadth = 0
    for i in range(1, len(starts)):
        prev_end = ends[i-1]
        curr_starts = starts[i]
        dist = curr_starts - prev_end
        if dist > 0:
            breadth += dist
    breadth = (starts[0] - 0) + breadth + (grp.iloc[0]['Length'] - ends[-1])
    breadth = breadth/grp.iloc[0]['Length'] * 100.0
    norm = np.sum(grp['alito'] - grp['alifrom'])
    score = np.sum((grp['alito'] - grp['alifrom'])*grp['Score']/(norm))
    cevlaue = np.sum((grp['alito'] - grp['alifrom'])*grp['c-Evalue']/norm)
    ievlaue = np.sum((grp['alito'] - grp['alifrom'])*grp['i-Evalue']/norm)
    
    return pd.Series({'breadth':breadth, 'Counts':len(grp), 'Avg. Score':score, 
                      'Avg. c-Evalue':cevlaue, 'Avg. i-Evalue':ievlaue})


#filepath = '/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Untagged_Accessory_Genes_Proteins_HMM_Alignment/Out/'
#tophit_dir = '/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Untagged_Accessory_Genes_Proteins_HMM_Alignments_Tophits/'
#untagged_gene_seq = '/fs/cbcb-scratch/hsmurali/CMSC829A/Data/Untagged_Accessory_Genes_Proteins/'

filepath = sys.argv[1]
tophit_dir = sys.argv[2]
untagged_gene_seq = sys.argv[3]
f = sys.argv[4]

if not isdir(tophit_dir):
    mkdir(tophit_dir)

try:
    if f.startswith("GCF") and f.endswith(".out"):
        lines = open(filepath+f).readlines()
        out = ''
        for line in lines:
            if line.startswith('#') or line.startswith('--') or line == "\n":
                continue
            else:
                out += line
        df_HMMER = Parse_HMM_Alignments(out)
        df_lengths = Calculate_Untagged_Gene_Lengths(untagged_gene_seq+f.replace(".out",".faa"))
        df_HMMER = pd.merge(df_HMMER, df_lengths, on = ['Query'], how = 'left')
        df_HMMER = df_HMMER[df_HMMER['Score'] >= 0]
        df_HMMER_grp = df_HMMER.groupby(['Query','Subject']).apply(Return_Breadth_of_Coverage)
        df_seqs = df_HMMER[['Query','Subject','Length','Subject_Alignment','Alignment_Description','Query_Alignment']]
        df_seqs = df_seqs.set_index(['Query', 'Subject'])
        df_HMMER_grp = df_HMMER_grp.join(df_seqs, how = 'left')
        df_HMMER_grp = df_HMMER_grp.reset_index()
        min_breadth = 60
        df_HMMER_filter = df_HMMER_grp[df_HMMER_grp['breadth'] >= min_breadth]
        df_HMMER_filter = df_HMMER_filter.loc[df_HMMER_filter.groupby('Query')['Avg. Score'].idxmax()]
        df_HMMER_filter.to_csv(tophit_dir+f, sep = '\t')
        print(f)
except IndexError:
    print('---->', f)


# In[ ]:




