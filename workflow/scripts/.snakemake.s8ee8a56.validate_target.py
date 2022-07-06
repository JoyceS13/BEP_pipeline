
######## Snakemake header ########
import sys; sys.path.insert(0, "/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/.conda/envs/art/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(X]\x00\x00\x00proportions2/20fcov/0.5_0.5props/sample17/vcf/comparison/gatk_upper_seq_0.diff.sites_in_filesq\x06X]\x00\x00\x00proportions2/20fcov/0.5_0.5props/sample17/vcf/comparison/gatk_upper_seq_1.diff.sites_in_filesq\x07Xr\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/IMBvariant/reference/chromosomes/CEN.PK113-7D_chr01_pilon.faq\x08e}q\t(X\x06\x00\x00\x00_namesq\n}q\x0b(X\x0f\x00\x00\x00comparison_vcfsq\x0cK\x00K\x02\x86q\rX\t\x00\x00\x00referenceq\x0eK\x02N\x86q\x0fuh\x0ccsnakemake.io\nNamedlist\nq\x10)\x81q\x11(h\x06h\x07e}q\x12h\n}q\x13sbh\x0eh\x08ubX\x06\x00\x00\x00outputq\x14csnakemake.io\nOutputFiles\nq\x15)\x81q\x16(X[\x00\x00\x00proportions2/20fcov/0.5_0.5props/sample17/analysis_target/analysis_target_gatk_upper_s0.csvq\x17X[\x00\x00\x00proportions2/20fcov/0.5_0.5props/sample17/analysis_target/analysis_target_gatk_upper_s1.csvq\x18e}q\x19h\n}q\x1asbX\x06\x00\x00\x00paramsq\x1bcsnakemake.io\nParams\nq\x1c)\x81q\x1d(X\x02\x00\x00\x0017q\x1eX\n\x00\x00\x00gatk_upperq\x1fX9\x00\x00\x00proportions2/20fcov/0.5_0.5props/sample17/analysis_targetq e}q!(h\n}q"(X\x0c\x00\x00\x00sample_indexq#K\x00N\x86q$X\x06\x00\x00\x00callerq%K\x01N\x86q&X\x07\x00\x00\x00out_dirq\'K\x02N\x86q(uh#h\x1eh%h\x1fh\'h ubX\t\x00\x00\x00wildcardsq)csnakemake.io\nWildcards\nq*)\x81q+(X \x00\x00\x00proportions2/20fcov/0.5_0.5propsq,h\x1eh\x1fe}q-(h\n}q.(X\x06\x00\x00\x00outdirq/K\x00N\x86q0X\x02\x00\x00\x00iiq1K\x01N\x86q2X\x06\x00\x00\x00callerq3K\x02N\x86q4uX\x06\x00\x00\x00outdirq5h,X\x02\x00\x00\x00iiq6h\x1eh%h\x1fubX\x07\x00\x00\x00threadsq7K\x01X\t\x00\x00\x00resourcesq8csnakemake.io\nResources\nq9)\x81q:(K\x01K\x01e}q;(h\n}q<(X\x06\x00\x00\x00_coresq=K\x00N\x86q>X\x06\x00\x00\x00_nodesq?K\x01N\x86q@uh=K\x01h?K\x01ubX\x03\x00\x00\x00logqAcsnakemake.io\nLog\nqB)\x81qC}qDh\n}qEsbX\x06\x00\x00\x00configqF}qG(X\t\x00\x00\x00referenceqHXr\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/IMBvariant/reference/chromosomes/CEN.PK113-7D_chr01_pilon.faqIX\x04\x00\x00\x00nameqJX\x04\x00\x00\x00chr1qKX\t\x00\x00\x00mutationsqLK\xc8X\x07\x00\x00\x00repeatsqMK\x1eX\x05\x00\x00\x00propsqN]qO(]qP(G?\xe0\x00\x00\x00\x00\x00\x00G?\xe0\x00\x00\x00\x00\x00\x00e]qQ(G?\xe3333333G?\xd9\x99\x99\x99\x99\x99\x9ae]qR(G?\xe9\x99\x99\x99\x99\x99\x9aG?\xc9\x99\x99\x99\x99\x99\x9ae]qS(G?\xe6ffffffG?\xd3333333e]qT(G?\xec\xcc\xcc\xcc\xcc\xcc\xcdG?\xb9\x99\x99\x99\x99\x99\x9ae]qU(G?\xeeffffffG?\xa9\x99\x99\x99\x99\x99\x9ae]qV(G?\xef\\(\xf5\xc2\x8f\\G?\x94z\xe1G\xae\x14{e]qW(G?\xef\xae\x14z\xe1G\xaeG?\x84z\xe1G\xae\x14{eeX\x05\x00\x00\x00f_covqX]qY(K\nK\x14K2KdK\xc8eX\x0b\x00\x00\x00working_dirqZX<\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/outputq[X\x0b\x00\x00\x00rename_fileq\\XX\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/rename_chr.txtq]X\x0e\x00\x00\x00varscan_filterq^X_\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/workflow/scripts/fpfilter.plq_X\x07\x00\x00\x00callersq`]qa(X\n\x00\x00\x00gatk_upperqbX\x0f\x00\x00\x00gatk_unfilteredqcX\x07\x00\x00\x00varscanqdX\x12\x00\x00\x00varscan_unfilteredqeX\t\x00\x00\x00freebayesqfX\x14\x00\x00\x00freebayes_unfilteredqgeX\x0f\x00\x00\x00experiment_nameqhX\x0c\x00\x00\x00proportions2qiuX\x04\x00\x00\x00ruleqjX\x0f\x00\x00\x00validate_targetqkub.')
######## Original script #########
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
    
    df.to_csv(open("{}/analysis_target_{}_s{}.csv".format(outdir,caller,target_index),"w"))
    
    return df
    
if __name__ == '__main__':
    
    ap = argparse.ArgumentParser(description = 'validate validates the called vcfs compared to the vcfs made from the true sequences.')
    ap.add_argument("comparison_vcfs", metavar = 'vcfs',  nargs = '*', type=str,  \
                    help="vcfs with all snps from the called and true sequences compared")
    ap.add_argument("-r","--reference", metavar='ref',  type=str,    \
                    help='fasta file of reference sequence')
    ap.add_argument("-i","--sampleindex", metavar='sidx',  type=int,    \
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
        sample_idx = args['sampleindex']
        caller = args['caller']
        out_dir = args['directory']
    
    vcfs.sort()
    for ii,file in enumerate(vcfs):
        validate_target(file, ii, vcfs, caller, ref, sample_idx, out_dir)
