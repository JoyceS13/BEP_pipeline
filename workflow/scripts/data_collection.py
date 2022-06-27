# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 01:27:23 2022

@author: wyjsu
"""

import pandas as pd
import argparse
import os

def data_collection(csvs, caller,out_dir):
    df = pd.DataFrame()
    for file in csvs:
        df.append( pd.read_csv(file, sep='\t'))
        
    df.to_csv(open(os.path.join("{}/{}".format(os.getcwd(),out_dir),"results_{}.csv".format(caller)),"w"),index=False)
    
    
if __name__ == '__main__':
    
    ap = argparse.ArgumentParser(description = 'validate validates the called vcfs compared to the vcfs made from the true sequences.')
    ap.add_argument("csvs", metavar = 'csvs',  nargs = '*', type=str,  \
                    help="csvs with all individual sample results")
    ap.add_argument("-c","--caller", metavar='caller',  type=str,    \
                    help='name of variant caller used')
    ap.add_argument("-d","--directory", metavar='dir',  default = '',type=str,    \
                    help='relative path to output directory')

    try:
        csvs = snakemake.input
        caller = snakemake.params['caller']
        out_dir = snakemake.params['out_dir']
        
    except:
        args = vars(ap.parse_args())
        csvs = args['csvs']
        caller = args['caller']
        out_dir = args['directory']
    
    data_collection(csvs, caller,out_dir)
    
