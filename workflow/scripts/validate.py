# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 22:06:15 2022

@author: wyjsu
"""

import pandas as pd
from Bio import SeqIO
import argparse
import os

def validate(comparison_vcfs, caller,  reference, sample_index, out_dir):
## validate validates the called vcfs compared to the vcfs made from the true sequences.
    #comparison_vcfs: vcfs with all snps from the called and true sequences compared
    #ref column is called, comparison is true vcf
    
    #create variables to store relevant information in
    snps = [] #list of chromsome positions that have been checked
    fp = [] #list of chromosome positions that have been marked false positive (for the moment)
    miscalls = [] #snp called different snp from true snp
    tp = [] #list of true positives
    TP = 0
    FN = 0
    
    #print(len(comparison_vcfs))
    #read files
    for file in comparison_vcfs:
        vcf = pd.read_csv(file, sep='\t')
        chr_pos = list(vcf["POS1"])
        call = list(vcf["ALT1"])
        true = list(vcf["ALT2"])
        
        for ii in range(len(chr_pos)):
            if str(chr_pos[ii]) not in snps :
                if call[ii] == true[ii]:
                    TP += 1 #true negatives are not included in comparison vcf
                    tp.append(chr_pos[ii])
                else:
                    if call[ii] == '.': #if no snp called, but there is a snp
                        FN += 1 #if not called, it won't be called in other vcf -> false negative
                    elif true[ii] == '.':
                        fp.append(chr_pos[ii]) #snp that is called could be from other sequence
                    else: #both don't match reference, but also don't match each other
                        miscalls.append(chr_pos[ii])
                if chr_pos[ii] != '.':
                    snps.append(str(chr_pos[ii]))

            elif str(chr_pos[ii]) in fp: #has been marked as fp for now
                if call[ii] == true[ii]:
                    TP += 1 #snp is from other sequence
                    fp.remove(str(chr_pos[ii]))
                    tp.append(chr_pos[ii])
                    #print("false positive re-evaluated")      
 
            elif str(chr_pos[ii]) in miscalls: #has been marked as miscall for now
                if call[ii] == true[ii]:
                    TP += 1 #snp is from other sequence
                    miscalls.remove(str(chr_pos[ii]))
        #print(" ".join(str(loc) for loc in snps))         
    
    #if there are still miscalls this will be outputted and the samples with miscalls thrown away
        
    
    ##calculate remaining values
    for ref in SeqIO.parse(open(reference),'fasta'):
      	ref_seq = str(ref.seq)
            
    FP = len(fp)
    n_miscalls = len(miscalls)
    TN = len(ref_seq) - TP - FN - FP - n_miscalls
      
    ##store results for output
    dict = {}
    dict['TP'] = TP
    dict['FP'] = FP
    dict['TN'] = TN
    dict['FN'] = FN
    dict['miscalls'] = n_miscalls
    dict['called_pos'] = [tp]
        
    df = pd.DataFrame(dict, index = [str(sample_index)])
    #df['called_pos'] = df['called_pos'].astype('object')
    #df.loc[sample_index,'called_pos'] = [fp]
        
    df.to_csv(open(os.path.join("{}/{}".format(os.getcwd(),out_dir),'analysis_{}.csv'.format(caller)),"w"))
    print(df)        

    return df
    
if __name__ == '__main__':
    
    ap = argparse.ArgumentParser(description = 'validate validates the called vcfs compared to the vcfs made from the true sequences.')
    ap.add_argument("comparison_vcfs", metavar = 'vcfs',  nargs = '*', type=str,  \
                    help="vcfs with all snps from the called and true sequences compared")
    ap.add_argument("-r","--reference", metavar='ref',  type=str,    \
                    help='fasta file of reference sequence')
    ap.add_argument("-c","--caller", metavar='caller',  type=str,    \
                    help='name of variant caller used')
    ap.add_argument("-i","--index", metavar='ref',  type=int,    \
                    help='index of sample')
    ap.add_argument("-d","--directory", metavar='dir', default='', type=str,    \
                    help='relative path to output directory')
    
    try:
        vcfs = snakemake.input['comparison_vcfs']
        ref = snakemake.input['reference']
        idx = snakemake.params['sample_index']
        caller = snakemake.params['caller']
        out_dir = snakemake.params['out_dir']
        
    except:
        args = vars(ap.parse_args())
        vcfs = args['comparison_vcfs']
        ref = args['reference']
        idx = args['index']
        caller = args['caller']
        out_dir = args['directory']
    
    validate(vcfs, caller, ref,idx,out_dir)
    

            
    
