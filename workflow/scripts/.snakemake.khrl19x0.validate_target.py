
######## Snakemake header ########
import sys; sys.path.insert(0, "/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/.conda/envs/art/lib/python3.6/site-packages"); import pickle; snakemake = pickle.loads(b'\x80\x03csnakemake.script\nSnakemake\nq\x00)\x81q\x01}q\x02(X\x05\x00\x00\x00inputq\x03csnakemake.io\nInputFiles\nq\x04)\x81q\x05(Xl\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/vcf/comparison/gatk_unfiltered_seq_0.diff.sites_in_filesq\x06Xl\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/vcf/comparison/gatk_unfiltered_seq_1.diff.sites_in_filesq\x07Xl\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/vcf/comparison/gatk_unfiltered_seq_2.diff.sites_in_filesq\x08Xl\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/vcf/comparison/gatk_unfiltered_seq_3.diff.sites_in_filesq\tXl\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/vcf/comparison/gatk_unfiltered_seq_4.diff.sites_in_filesq\nXr\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/IMBvariant/reference/chromosomes/CEN.PK113-7D_chr01_pilon.faq\x0be}q\x0c(X\x06\x00\x00\x00_namesq\r}q\x0e(X\x0f\x00\x00\x00comparison_vcfsq\x0fK\x00K\x05\x86q\x10X\t\x00\x00\x00referenceq\x11K\x05N\x86q\x12uh\x0fcsnakemake.io\nNamedlist\nq\x13)\x81q\x14(h\x06h\x07h\x08h\th\ne}q\x15h\r}q\x16sbh\x11h\x0bubX\x06\x00\x00\x00outputq\x17csnakemake.io\nOutputFiles\nq\x18)\x81q\x19(Xj\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/analysis_target/analysis_target_gatk_unfiltered_s0.csvq\x1aXj\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/analysis_target/analysis_target_gatk_unfiltered_s1.csvq\x1bXj\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/analysis_target/analysis_target_gatk_unfiltered_s2.csvq\x1cXj\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/analysis_target/analysis_target_gatk_unfiltered_s3.csvq\x1dXj\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/analysis_target/analysis_target_gatk_unfiltered_s4.csvq\x1ee}q\x1fh\r}q sbX\x06\x00\x00\x00paramsq!csnakemake.io\nParams\nq")\x81q#(X\x02\x00\x00\x0025q$X\x0f\x00\x00\x00gatk_unfilteredq%XC\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05props/sample25/analysis_targetq&e}q\'(h\r}q((X\x0c\x00\x00\x00sample_indexq)K\x00N\x86q*X\x06\x00\x00\x00callerq+K\x01N\x86q,X\x07\x00\x00\x00out_dirq-K\x02N\x86q.uh)h$h+h%h-h&ubX\t\x00\x00\x00wildcardsq/csnakemake.io\nWildcards\nq0)\x81q1(X*\x00\x00\x00ploidy5/100fcov/0.5_0.2_0.15_0.1_0.05propsq2h$h%e}q3(h\r}q4(X\x06\x00\x00\x00outdirq5K\x00N\x86q6X\x02\x00\x00\x00iiq7K\x01N\x86q8X\x06\x00\x00\x00callerq9K\x02N\x86q:uX\x06\x00\x00\x00outdirq;h2X\x02\x00\x00\x00iiq<h$h+h%ubX\x07\x00\x00\x00threadsq=K\x01X\t\x00\x00\x00resourcesq>csnakemake.io\nResources\nq?)\x81q@(K\x01K\x01e}qA(h\r}qB(X\x06\x00\x00\x00_coresqCK\x00N\x86qDX\x06\x00\x00\x00_nodesqEK\x01N\x86qFuhCK\x01hEK\x01ubX\x03\x00\x00\x00logqGcsnakemake.io\nLog\nqH)\x81qI}qJh\r}qKsbX\x06\x00\x00\x00configqL}qM(X\t\x00\x00\x00referenceqNXr\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/IMBvariant/reference/chromosomes/CEN.PK113-7D_chr01_pilon.faqOX\x04\x00\x00\x00nameqPX\x04\x00\x00\x00chr1qQX\t\x00\x00\x00mutationsqRK\xc8X\x07\x00\x00\x00repeatsqSK2X\x05\x00\x00\x00propsqT]qU(]qV(G?\xe9\x99\x99\x99\x99\x99\x9aG?\xa9\x99\x99\x99\x99\x99\x9aG?\xa9\x99\x99\x99\x99\x99\x9aG?\xa9\x99\x99\x99\x99\x99\x9aG?\xa9\x99\x99\x99\x99\x99\x9ae]qW(G?\xe0\x00\x00\x00\x00\x00\x00G?\xc9\x99\x99\x99\x99\x99\x9aG?\xc3333333G?\xb9\x99\x99\x99\x99\x99\x9aG?\xa9\x99\x99\x99\x99\x99\x9ae]qX(G?\xc9\x99\x99\x99\x99\x99\x9aG?\xc9\x99\x99\x99\x99\x99\x9aG?\xc9\x99\x99\x99\x99\x99\x9aG?\xc9\x99\x99\x99\x99\x99\x9aG?\xc9\x99\x99\x99\x99\x99\x9aeeX\x05\x00\x00\x00f_covqY]qZ(K2KdK\xc8M\xf4\x01eX\x0b\x00\x00\x00working_dirq[X<\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/outputq\\X\x0b\x00\x00\x00rename_fileq]XX\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/config/rename_chr.txtq^X\x0e\x00\x00\x00varscan_filterq_X_\x00\x00\x00/tudelft.net/staff-umbrella/YeastVariantCalling/wsung/BEP_pipeline/workflow/scripts/fpfilter.plq`X\x07\x00\x00\x00callersqa]qb(X\x04\x00\x00\x00gatkqcX\t\x00\x00\x00freebayesqdX\n\x00\x00\x00gatk_upperqeX\n\x00\x00\x00gatk_lowerqfX\x0f\x00\x00\x00gatk_unfilteredqgX\x07\x00\x00\x00varscanqhX\x14\x00\x00\x00freebayes_unfilteredqieX\x0f\x00\x00\x00experiment_nameqjX\x07\x00\x00\x00ploidy5qkuX\x04\x00\x00\x00ruleqlX\x0f\x00\x00\x00validate_targetqmub.')
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
