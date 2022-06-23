# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:44:56 2022

@author: wyjsu
"""
from Bio import SeqIO
import numpy as np


def comparison(reference, strings):
## finds the similarity between two strings of the same length
# reference: fasta file containing the reference sequence(s)
# strings: fasta file containing the string(s) that need to be compared to the reference(s)
#returns:
    #ref_ids: all the id's of reference sequences
    #str_ids: all the id's of string sequences
    #similarity: similarity between references and strings, with strings along columns and references over the rows.
    
    #read from files
    ref_ids = []    #list of strings that have the id names of the references
    ref_seqs = []   #list of strings with the sequences of the references
    ref_sequences = SeqIO.parse(open(reference),'fasta')
    for ref in ref_sequences:
        ref_ids.append(ref.id)
        ref_seqs.append(str(ref.seq))
        
    str_ids = []    #list of strings that have the id names of the comparison strings
    str_seqs = []   #list of strings with the sequences of the comparison strings
    str_sequences = SeqIO.parse(open(strings),'fasta')
    for string in str_sequences:
        str_ids.append(string.id)
        str_seqs.append(str(string.seq))
        

    #comparison    
    similarity = np.zeros((len(ref_ids), len(str_ids))) #output array with similarities
    #loop through every character of every string compared to every reference
    for ii, ref in enumerate(ref_seqs):
        for jj, string in enumerate(str_seqs):
            for kk, char in enumerate(string):
                if char == ref[kk]:
                    similarity[ii,jj] += 1
    #normalize by string length
    similarity = similarity/len(ref[0])
    
    return ref_ids, str_ids, similarity

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description = 'compares two sequences and calculates how alike they are')
    ap.add_argument('reference', metavar = 'ref',  type=str, \
                    help="Input fasta file")
    ap.add_argument('sequence', metavar = 'seq',  type=str, \
                    help="Input fasta file")

    args = args = vars(ap.parse_args())
    
    ref = args['reference']
    seq = args['sequence']
    
    print(comparison(ref,seq)[2])
    