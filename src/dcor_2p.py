# distance correlation - flu virus - sequences
# uday yallapragada

from numpy import array,transpose, var
from scipy.stats import mode
from Bio import AlignIO
from math import isnan
import sys,copy
from minepy import *
import numpy as np
from scipy.spatial.distance import pdist, squareform

def generate_binary_sequence(fasta_list):
    """ Generates a binary sequence out of an MSA.

    Arguments:
    fasta_list -- A list of protein sequences.
    """
    mod = [mode(x)[0][0] for x in transpose(array([list(z) for z in fasta_list]))]
    return array([[1 if x==mod[i] else 0 for i,x in enumerate(y)] for y in fasta_list])

#we need this to concatenate sequences with same strain name from two proteins
def generate_list_with_strain_name(binary_sequence_list, align):
    x = copy.copy(binary_sequence_list)
    for counter, record in enumerate(x):
        strain_full_str = align[counter].description.split('|')[3]
        seq_id = align[counter].description.split('|')[0]
        if (strain_full_str.startswith('Organism')):
            strain = strain_full_str[29:]
            record.insert(0, strain)
        else:
            strain_full_str = align[counter].id.split('|')[2]
            strain = strain_full_str[29:]
            record.insert(0, strain)
        record.append(9)        
    return x

def find(list_of_tuples, value):
    try:
        return next(x for x in list_of_tuples if value in x)
    except:
        return None

#concatenate two sequence lists based on strain name of each sequence
def concatTwoLists(bin_seq_list1, bin_seq_list2):
    final_list = list()
    for record1 in bin_seq_list1:
        print(record1)
        strain = record1[0]
        record2 = find(bin_seq_list2,strain)
        print(record2)
        if (record2):
            final_list.append(record1[1:]+record2[1:])
    print(final_list)
    final_listT = transpose(array([list(z) for z in final_list])).tolist()
    return final_listT

def distcorr(X, Y):

    X = np.atleast_1d(X)
    Y = np.atleast_1d(Y)
    if np.prod(X.shape) == len(X):
        X = X[:, None]
    if np.prod(Y.shape) == len(Y):
        Y = Y[:, None]
    X = np.atleast_2d(X)
    Y = np.atleast_2d(Y)
    n = X.shape[0]
    if Y.shape[0] != X.shape[0]:
        raise ValueError('Number of samples must match')
    a = squareform(pdist(X))
    b = squareform(pdist(Y))
    A = a - a.mean(axis=0)[None, :] - a.mean(axis=1)[:, None] + a.mean()
    B = b - b.mean(axis=0)[None, :] - b.mean(axis=1)[:, None] + b.mean()
    
    dcov2_xy = (A * B).sum()/float(n * n)
    dcov2_xx = (A * A).sum()/float(n * n)
    dcov2_yy = (B * B).sum()/float(n * n)
    dcor = np.sqrt(dcov2_xy)/np.sqrt(np.sqrt(dcov2_xx) * np.sqrt(dcov2_yy))
    return dcor

def performDCor(combined_list, l1, l2):
    print (l1, l2, len(combined_list))
    outfile2 = open(prots+"_dcor"+".csv", 'w')
    outfile2.writelines('first,' + 'second,' + 'dCor,' + '\n')
 
    for counter1 in range(0, l1-1):
        for counter2 in range(l1, len(combined_list)-1):
            dCor = distcorr(combined_list[counter1],combined_list[counter2])
            if not(isnan(dCor)):
                outfile2.writelines(str(counter1)+','+str(counter2-l1)+','+ "{:.3f}".format(dCor))
            else:
                outfile2.writelines(str(counter1)+','+str(counter2-l1)+','+str(0.0))          
            outfile2.write('\n')
            outfile2.flush()
            
    outfile2.close()

if __name__ == '__main__':
    align1 = AlignIO.read(sys.argv[1], 'fasta')
    align2 = AlignIO.read(sys.argv[2], 'fasta')
    prots = sys.argv[3]
    
    seq1 = [x.seq for x in align1]
    seq2 = [x.seq for x in align2]
    
    bin_seq2 = generate_binary_sequence(seq2).tolist()
    bin_seq1 = generate_binary_sequence(seq1).tolist()
    
    bin_seq_list1 = generate_list_with_strain_name(bin_seq1, align1)
    bin_seq_list2 = generate_list_with_strain_name(bin_seq2, align2)
    
    l1=len(bin_seq_list1[0])
    l2=len(bin_seq_list2[0])
    
    print (l1, " ", l2)
    final_listT = concatTwoLists(bin_seq_list1, bin_seq_list2)
    performDCor(final_listT, l1, l2)
    