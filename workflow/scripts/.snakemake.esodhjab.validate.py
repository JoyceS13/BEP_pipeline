
######## Snakemake header ########
import sys; sys.path.insert(0, "/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/.conda/envs/art/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X[\x00\x00\x00proportions2/50fcov/0.7_0.3props/sample7/vcf/comparison/freebayes_seq_0.diff.sites_in_filesq\x06X[\x00\x00\x00proportions2/50fcov/0.7_0.3props/sample7/vcf/comparison/freebayes_seq_1.diff.sites_in_filesq\x07Xr\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/IMBvariant/reference/chromosomes/CEN.PK113-7D_chr01_pilon.faq\x08e}q\t(X\x06\x00\x00\x00_namesq\n}q\x0b(X\x0f\x00\x00\x00comparison_vcfsq\x0cK\x00K\x02\x86q\rX\t\x00\x00\x00referenceq\x0eK\x02N\x86q\x0fuh\x0ccsnakemake.io\nNamedlist\nq\x10)\x81q\x11(h\x06h\x07e}q\x12h\n}q\x13sbh\x0eh\x08ubX\x06\x00\x00\x00outputq\x14csnakemake.io\nOutputFiles\nq\x15)\x81q\x16XH\x00\x00\x00proportions2/50fcov/0.7_0.3props/sample7/analysis/analysis_freebayes.csvq\x17a}q\x18h\n}q\x19sbX\x06\x00\x00\x00paramsq\x1acsnakemake.io\nParams\nq\x1b)\x81q\x1c(X\x01\x00\x00\x007q\x1dX\t\x00\x00\x00freebayesq\x1eX1\x00\x00\x00proportions2/50fcov/0.7_0.3props/sample7/analysisq\x1fe}q (h\n}q!(X\x0c\x00\x00\x00sample_indexq"K\x00N\x86q#X\x06\x00\x00\x00callerq$K\x01N\x86q%X\x07\x00\x00\x00out_dirq&K\x02N\x86q\'uh"h\x1dh$h\x1eh&h\x1fubX\t\x00\x00\x00wildcardsq(csnakemake.io\nWildcards\nq))\x81q*(X \x00\x00\x00proportions2/50fcov/0.7_0.3propsq+h\x1dh\x1ee}q,(h\n}q-(X\x06\x00\x00\x00outdirq.K\x00N\x86q/X\x02\x00\x00\x00iiq0K\x01N\x86q1X\x06\x00\x00\x00callerq2K\x02N\x86q3uX\x06\x00\x00\x00outdirq4h+X\x02\x00\x00\x00iiq5h\x1dh$h\x1eubX\x07\x00\x00\x00threadsq6K\x01X\t\x00\x00\x00resourcesq7csnakemake.io\nResources\nq8)\x81q9(K\x01K\x01e}q:(h\n}q;(X\x06\x00\x00\x00_coresq<K\x00N\x86q=X\x06\x00\x00\x00_nodesq>K\x01N\x86q?uh<K\x01h>K\x01ubX\x03\x00\x00\x00logq@csnakemake.io\nLog\nqA)\x81qB}qCh\n}qDsbX\x06\x00\x00\x00configqE}qF(X\t\x00\x00\x00referenceqGXr\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/IMBvariant/reference/chromosomes/CEN.PK113-7D_chr01_pilon.faqHX\x04\x00\x00\x00nameqIX\x04\x00\x00\x00chr1qJX\t\x00\x00\x00mutationsqKK\xc8X\x07\x00\x00\x00repeatsqLK\x1eX\x05\x00\x00\x00propsqM]qN(]qO(G?\xe0\x00\x00\x00\x00\x00\x00G?\xe0\x00\x00\x00\x00\x00\x00e]qP(G?\xe3333333G?\xd9\x99\x99\x99\x99\x99\x9ae]qQ(G?\xe9\x99\x99\x99\x99\x99\x9aG?\xc9\x99\x99\x99\x99\x99\x9ae]qR(G?\xe6ffffffG?\xd3333333e]qS(G?\xec\xcc\xcc\xcc\xcc\xcc\xcdG?\xb9\x99\x99\x99\x99\x99\x9ae]qT(G?\xeeffffffG?\xa9\x99\x99\x99\x99\x99\x9ae]qU(G?\xef\\(\xf5\xc2\x8f\\G?\x94z\xe1G\xae\x14{e]qV(G?\xef\xae\x14z\xe1G\xaeG?\x84z\xe1G\xae\x14{eeX\x05\x00\x00\x00f_covqW]qX(K\nK\x14K2KdK\xc8eX\x0b\x00\x00\x00working_dirqYX<\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/outputqZX\x0b\x00\x00\x00rename_fileq[XX\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/rename_chr.txtq\\X\x0e\x00\x00\x00varscan_filterq]X_\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/workflow/scripts/fpfilter.plq^X\x07\x00\x00\x00callersq_]q`(X\n\x00\x00\x00gatk_upperqaX\x0f\x00\x00\x00gatk_unfilteredqbX\x07\x00\x00\x00varscanqcX\x12\x00\x00\x00varscan_unfilteredqdX\t\x00\x00\x00freebayesqeX\x14\x00\x00\x00freebayes_unfilteredqfeX\x0f\x00\x00\x00experiment_nameqgX\x0c\x00\x00\x00proportions2qhuX\x04\x00\x00\x00ruleqiX\x08\x00\x00\x00validateqjub.')
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
    

            
    
