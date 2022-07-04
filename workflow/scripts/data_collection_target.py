# -*- coding: utf-8 -*-
"""
Created on Fri Jul  1 00:51:00 2022

@author: wyjsu
"""

import pandas as pd
import argparse
import os

def data_collection(csvs, file_name, out_dir):
    df = pd.DataFrame()
    for file in csvs:
        if df.empty:
            df =  pd.read_csv(file, index_col = 0,header=0)
            #print(df)
        else: 
            temp = pd.read_csv(file, index_col=0,header=0)
            print(temp)
            df = df.append(temp)
            #print(df)
        
    df.to_csv(open(out_dir+"/"+file_name,"w"))
    
    
if __name__ == '__main__':
    
    ap = argparse.ArgumentParser(description = 'validate validates the called vcfs compared to the vcfs made from the true sequences.')
    ap.add_argument("csvs", metavar = 'csvs',  nargs = '*', type=str,  \
                    help="csvs with all individual sample results")
    ap.add_argument("-f","--filename", metavar='caller',  type=str,    \
                    help='name of output file (incl. csv)')
    ap.add_argument("-o","--output", metavar='dir',  default = '',type=str,    \
                    help='relative path to output directory')

    try:
        csvs = snakemake.input
        file_name = snakemake.params['file_name']
        out_dir = snakemake.params['out_dir']
        
    except:
        args = vars(ap.parse_args())
        csvs = args['csvs']
        file_name = args['filename']
        out_dir = args['directory']
    
    data_collection(csvs, file_name, out_dir)
