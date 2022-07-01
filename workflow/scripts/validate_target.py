# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 14:04:25 2022

@author: wyjsu
"""

import pandas as pd
from Bio import SeqIO
import argparse

def validate_target(target_vcf, target_index, comparison_vcfs, caller, reference, sample_index, outdir):
## validate_seq validates the called vcfs compared to the vcf made for one specific target sequence.
## The other vcfs are taking into consideration to remove false positives
    #target_vcf: vcf with called snps compared to target snps
    #comparison_vcfs: vcfs with all snps from the called and true sequences
    #ref column is called, comparison is true vcf
    
    #create variables to store relevant information in
    snps = [] #list of chromsome positions that have been checked
    fp = [] #list of chromosome positions that have been marked false positive (for the moment)
    miscalls = [] #snp called different snp from true snp
    TP = 0
    FN = 0
    ambig = 0 #ambigious calls: SNP in target but also other mutated sequence, other mutated sequence is called
    
    target = pd.read_csv(target_vcf, sep='\t')
    chr_pos = list(target["POS1"])
    call = list(target["ALT1"])
    true = list(target["ALT2"])
    
    for ii in range(len(chr_pos)):
        if call[ii] == true[ii]:
            TP += 1 #true negatives are not included in comparison vcf
        else:
            if call[ii] == '.': #if no snp called, but there is a snp
                FN += 1 #if not called, it won't be called in other vcf -> false negative
            elif true[ii] == '.':
                fp.append(chr_pos[ii]) #snp that is called could be from other sequence
            else: #both don't match reference, but also don't match each other
                miscalls.append(chr_pos[ii])
        if chr_pos[ii] != '.':
            snps.append(str(chr_pos[ii]))
    
    #check if false positives and miscalls are not from other sequences
    for file in comparison_vcfs:
        #read files
        vcf = pd.read_csv(file, sep='\t')
        chr_pos = list(vcf["POS1"])
        call = list(vcf["ALT1"])
        true = list(vcf["ALT2"])
        
        
        for ii in range(len(chr_pos)):
            #no false negatives from other files, as we are looking at a specific target
        
            if str(chr_pos[ii]) in fp: #has been marked as fp for now
                if call[ii] == true[ii]:
                    #snp is from other sequence, mark as true negative (true negatives calculated at end)
                    fp.remove(str(chr_pos[ii]))
            
            elif str(chr_pos[ii]) in miscalls: #has been marked as miscall for now
                if call[ii] == true[ii]:
                    ambig += 1 #snp is from other sequence, can't determine whether was called correctly
                    miscalls.remove(str(chr_pos[ii]))
    
    #if there are still miscalls this will be outputted and the samples with miscalls thrown away
    
    for ref in SeqIO.parse(open(reference),'fasta'):
      	ref_seq = str(ref.seq)
            
    FP = len(fp)
    n_miscalls = len(miscalls)
    TN = len(ref_seq) - TP - FN - FP - n_miscalls - ambig
    
    
    dict = {}
    dict['TP'] = TP
    dict['FP'] = FP
    dict['TN'] = TN
    dict['FN'] = FN
    dict['miscalls'] = n_miscalls
    dict['ambigious'] = ambig

    df = pd.DataFrame(dict, index = [str(sample_index)])
    
    df.to_csv(open(outdir+"/","analysis_target_{}_s{}.csv".format(caller,target_index),"w"))
    
    return df
    
if __name__ == '__main__':
    
    ap = argparse.ArgumentParser(description = 'validate validates the called vcfs compared to the vcfs made from the true sequences.')
    ap.add_argument("-c","--comparison_vcfs", metavar = 'vcfs',  nargs = '*', type=str,  \
                    help="vcfs with all snps from the called and true sequences compared")
    ap.add_argument("-r","--reference", metavar='ref',  type=str,    \
                    help='fasta file of reference sequence')
    ap.add_argument("-i","--sample-index", metavar='sidx',  type=int,    \
                    help='index of sample')
    ap.add_argument("-c","--caller", metavar='caller',  type=str,    \
                    help='name of variant caller used')
    ap.add_argument("-d","--directory", metavar='dir', default='', type=str,    \
                    help='relative path to output directory')
    

    try:
        vcfs = snakemake.input['comparison_vcfs']
        ref = snakemake.input['reference']
        sample_idx = snakemake.params['sample_index']
        caller = snakemake.params['caller']
        out_dir = snakemake.params['out_dir']
    except:
        args = vars(ap.parse_args())
        vcfs = args['comparison_vcfs']
        ref = args['reference']
        sample_idx = args['sample-index']
        caller = args['caller']
        out_dir = args['directory']
    
    for ii,file in enumerate(vcfs.sort()):
        validate_target(file, ii, vcfs, caller, ref, sample_idx, out_dir)