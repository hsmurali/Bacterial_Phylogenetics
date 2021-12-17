import sys
import numpy as np
from Bio import SeqIO as IO
from pymsa import MSA, PercentageOfNonGaps, PercentageOfTotallyConservedColumns 

def Load_MSA(filepath):
    record_dict = IO.to_dict(IO.parse(filepath, "fasta"))
    d = []
    for key in record_dict.items():
        d.append((key[0], key[1].seq))
    return d

def Hamming_Distance(s1, s2):
    result = 0
    for x,(p,q) in enumerate(zip(s1,s2)):
        if p!=q:
            result+=1
    result = result/len(s1)
    return result

def Return_Pairwise_Hamming_Distance(aligmments):
    distances = []
    for i in range(0, len(alignments)):
        for j in range(i+1, len(alignments)):
            distances.append(Hamming_Distance(alignments[i][1], alignments[j][1]))
    return distances

def MSA_Stats(alignments):
    n = len(alignments)
    distances = Return_Pairwise_Hamming_Distance(alignments)
    avg_p_dist = np.mean(distances)
    median_p_dist = np.median(distances)
    std_p_dist = np.std(distances)
    min_p_dist = np.min(distances)
    max_p_dist = np.max(distances)

    aligned_sequences = list(pair[1] for pair in alignments)
    sequences_id = list(pair[0] for pair in alignments)
    msa = MSA(aligned_sequences, sequences_id)
    
    non_gaps = PercentageOfNonGaps(msa)
    totally_conserved_columns = PercentageOfTotallyConservedColumns(msa)
    frac_non_gaps = non_gaps.compute() 
    frac_conserved = totally_conserved_columns.compute()

    d = {'Num_Seqs':n, 'Avg_P_Dist':avg_p_dist, 'Median_P_Dist':median_p_dist,
         'Min_P_Dist':min_p_dist, 'Max_P_Dist':max_p_dist,
         'Dev_P_Dist': std_p_dist, 'Frac_Non_Gaps':frac_non_gaps, 'Frac_Conserved':frac_conserved}
    return d

alignment_file = sys.argv[1]
outfile = sys.argv[2]
alignments = Load_MSA(alignment_file)
d = MSA_Stats(alignments)

f = open(outfile.replace(".faa",".txt"),'w')
f.write(str(d))
f.close()
