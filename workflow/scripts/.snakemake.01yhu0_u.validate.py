
######## Snakemake header ########
import sys; sys.path.insert(0, "/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/.conda/envs/art/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X_\x00\x00\x00ploidy_3/20fcov/0.33_0.33_0.33props/sample12/vcf/comparison/freebayes_seq_0.diff.sites_in_filesq\x06X_\x00\x00\x00ploidy_3/20fcov/0.33_0.33_0.33props/sample12/vcf/comparison/freebayes_seq_1.diff.sites_in_filesq\x07X_\x00\x00\x00ploidy_3/20fcov/0.33_0.33_0.33props/sample12/vcf/comparison/freebayes_seq_2.diff.sites_in_filesq\x08Xr\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/IMBvariant/reference/chromosomes/CEN.PK113-7D_chr01_pilon.faq\te}q\n(X\x06\x00\x00\x00_namesq\x0b}q\x0c(X\x0f\x00\x00\x00comparison_vcfsq\rK\x00K\x03\x86q\x0eX\t\x00\x00\x00referenceq\x0fK\x03N\x86q\x10uh\rcsnakemake.io\nNamedlist\nq\x11)\x81q\x12(h\x06h\x07h\x08e}q\x13h\x0b}q\x14sbh\x0fh\tubX\x06\x00\x00\x00outputq\x15csnakemake.io\nOutputFiles\nq\x16)\x81q\x17XL\x00\x00\x00ploidy_3/20fcov/0.33_0.33_0.33props/sample12/analysis/analysis_freebayes.csvq\x18a}q\x19h\x0b}q\x1asbX\x06\x00\x00\x00paramsq\x1bcsnakemake.io\nParams\nq\x1c)\x81q\x1d(X\x02\x00\x00\x0012q\x1eX\t\x00\x00\x00freebayesq\x1fX5\x00\x00\x00ploidy_3/20fcov/0.33_0.33_0.33props/sample12/analysisq e}q!(h\x0b}q"(X\x0c\x00\x00\x00sample_indexq#K\x00N\x86q$X\x06\x00\x00\x00callerq%K\x01N\x86q&X\x07\x00\x00\x00out_dirq\'K\x02N\x86q(uh#h\x1eh%h\x1fh\'h ubX\t\x00\x00\x00wildcardsq)csnakemake.io\nWildcards\nq*)\x81q+(X#\x00\x00\x00ploidy_3/20fcov/0.33_0.33_0.33propsq,h\x1eh\x1fe}q-(h\x0b}q.(X\x06\x00\x00\x00outdirq/K\x00N\x86q0X\x02\x00\x00\x00iiq1K\x01N\x86q2X\x06\x00\x00\x00callerq3K\x02N\x86q4uX\x06\x00\x00\x00outdirq5h,X\x02\x00\x00\x00iiq6h\x1eh%h\x1fubX\x07\x00\x00\x00threadsq7K\x01X\t\x00\x00\x00resourcesq8csnakemake.io\nResources\nq9)\x81q:(K\x01K\x01e}q;(h\x0b}q<(X\x06\x00\x00\x00_coresq=K\x00N\x86q>X\x06\x00\x00\x00_nodesq?K\x01N\x86q@uh=K\x01h?K\x01ubX\x03\x00\x00\x00logqAcsnakemake.io\nLog\nqB)\x81qC}qDh\x0b}qEsbX\x06\x00\x00\x00configqF}qG(X\t\x00\x00\x00referenceqHXr\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/IMBvariant/reference/chromosomes/CEN.PK113-7D_chr01_pilon.faqIX\x04\x00\x00\x00nameqJX\x04\x00\x00\x00chr1qKX\t\x00\x00\x00mutationsqLK\xc8X\x07\x00\x00\x00repeatsqMK\x1eX\x05\x00\x00\x00propsqN]qO(]qP(G?\xec\xcc\xcc\xcc\xcc\xcc\xcdG?\xa9\x99\x99\x99\x99\x99\x9aG?\xa9\x99\x99\x99\x99\x99\x9ae]qQ(G?\xeb333333G?\xb9\x99\x99\x99\x99\x99\x9aG?\xa9\x99\x99\x99\x99\x99\x9ae]qR(G?\xd5\x1e\xb8Q\xeb\x85\x1fG?\xd5\x1e\xb8Q\xeb\x85\x1fG?\xd5\x1e\xb8Q\xeb\x85\x1fe]qS(G?\xdc\xcc\xcc\xcc\xcc\xcc\xcdG?\xdc\xcc\xcc\xcc\xcc\xcc\xcdG?\xb9\x99\x99\x99\x99\x99\x9aeeX\x05\x00\x00\x00f_covqT]qU(K\x14K2KdK\xc8eX\x0b\x00\x00\x00working_dirqVX<\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/outputqWX\x0b\x00\x00\x00rename_fileqXXX\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/rename_chr.txtqYX\x0e\x00\x00\x00varscan_filterqZX_\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/workflow/scripts/fpfilter.plq[X\x07\x00\x00\x00callersq\\]q](X\x04\x00\x00\x00gatkq^X\n\x00\x00\x00gatk_upperq_X\x0f\x00\x00\x00gatk_unfilteredq`X\x07\x00\x00\x00varscanqaX\t\x00\x00\x00freebayesqbX\x14\x00\x00\x00freebayes_unfilteredqceX\x0f\x00\x00\x00experiment_nameqdX\x08\x00\x00\x00ploidy_3qeuX\x04\x00\x00\x00ruleqfX\x08\x00\x00\x00validateqgub.')
######## Original script #########
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
    

            
    
