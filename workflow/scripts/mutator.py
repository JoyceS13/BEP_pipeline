# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 11:21:04 2022

@author: wyjsu
"""
from Bio import SeqIO
import random as rand
import argparse
import os

def mutator(reference, mutations, clones, file_name = '' , out_dir):
## mutator introduces a {mutations} number of mutations to a {clones} number of clones from a reference sequence
    #reference: Fasta file containing reference sequence
    #mutations: Integer with the number of mutations that need to be introduced
    #clones: number of clones to make
    #file_name: optionally give file_name
#makes an output file with 
    
    #read file
    for ref in SeqIO.parse(open(reference),'fasta'):
    	ref_id = ref.id
    	ref_seq = str(ref.seq)
    clones = int(clones)
    mutations = int(mutations)

    #prepare output file
    if file_name == '':
        file_name = '{}_{}mut'.format(ref_id,mutations)
    
    #mutates and adds clone to file
    for ii in range(clones):
        clone = list(ref_seq)
        random = []
        for jj in range(mutations):
            #potential nucleotides
            nucl = ['A','T','C','G']
            #chooses random index to mutate
            idx = rand.randint(0,len(clone)-1)
            #if the site has already been mutated a new site is chosen
            while idx in random:
                idx = rand.randint(0,len(clone)-1)
            nucl.remove(clone[idx]) #removes read nucleotide from potential nucleotide list
            clone[idx] = rand.choice(nucl) #assigns a new nucleotide that is different from the old nucleotide
        clone = "".join(clone)   
        f = open(os.path.join("{}/{}".format(os.getcwd(),out_dir),file_name+'_'+str(ii)+'.fasta','w'))
        f.write('>{}_{}mut_clone_{}\n{}\n'.format(ref_id,mutations,ii,clone))
        f.close()
              
    similarity = 1 - mutations/len(ref_seq)
    print('The similarity of the clones to the reference is ',similarity)
    
if __name__ == '__main__':
    
    ap = argparse.ArgumentParser(description = 'mutator introduces a {mutations} number of mutations to a {clones} number of clones from a reference sequence')
    ap.add_argument("in_file", metavar = 'in_file',  type=str,  \
                    help="Input fasta file")
    ap.add_argument("mutations", metavar='mut',  type=int,    \
                    help='number of mutations to be introduced')
    ap.add_argument("clones", metavar = 'N',  type = int,   \
                    help='number of clones to be produced')
    ap.add_argument('-o','--outfile', nargs='?', required=False, default = '', type = str, \
                        help = 'prefix/name of out file' )
    ap.add_argument('-d','--directory', nargs='?', required=False, default = '', type = str, \
                        help = 'path to output directory' )
    #if there are snakemake variables, use those, otherwise get arguments from command line
    try:
        in_file = snakemake.input[0]
        mutations = snakemake.params["mutations"]
        clones = snakemake.params["pop_size"]
        outfile = snakemake.params["out_prefix"]
        out_dir = snakemake.params["out_dir"]
    except:
        args = vars(ap.parse_args())
        in_file = args['in_file']
        mutations = args['mutations']
        clones = args['clones']
        outfile = args['outfile']
        out_dir = args['directory']
    
    mutator(in_file, mutations, clones, outfile,out_dir)
    
    
            
    
