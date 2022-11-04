# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 18:23:01 2022

@author: yl183
"""
import pandas as pd 
import argparse
import numpy as np
from collections import defaultdict
import bisect

def dirichlet_sample(alphas):
    y_vals = []
    x_vals = []
    sum_gamma = 0
    for i in range(0, len(alphas)):
        y_vals.append(np.random.gamma(alphas[i],1.0))
        sum_gamma += y_vals[i]
    for j in range(0, len(alphas)):
        x_vals.append(y_vals[j]/sum_gamma)
    return x_vals

def sample_reads(nreads,d):
    count = defaultdict(int)
    cs = np.cumsum(d)
    L = []
    for i in range(0, nreads):
        s = bisect.bisect(cs, np.random.random())
        count[s] += 1
    for i in range(0, len(d)):
        L.append(count[i])
    return L
        


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('-in','--path to cnp matrix', required = True, help = "Path to cnp matrix ")
    ap.add_argument("-nreads","--number of reads",required=True, help="Number of reads to be sampled for each cell")
    ap.add_argument("-out","--path to output rc matrix",required=True, help="Path to generated read count matrix ")
    args = vars(ap.parse_args())
    if args['path to cnp matrix']!=None:
        path = args['path to cnp matrix']
    if args['number of reads']!=None:
        nreads = args['number of reads']
    if args['path to output rc matrix']!=None:
        out = args['path to output rc matrix']
    data = pd.read_csv(path, sep = '\t')
    data_prefix = data[['CHR','START','END ']]
    data = data.drop(['CHR','START','END '], axis = 1)
    over_dispersed = True
    nu = 4
    nreads = int(nreads)
    if not over_dispersed:
         data = data.div(data.sum(axis=0), axis=1)
    else:
        data = data.mul(nu)
        for column in data:
            alphas = data[column]
            data[column] = dirichlet_sample(alphas)
    for column in data:
        d = data[column]
        data[column] = sample_reads(nreads, d)
    data = data_prefix.join(data)
    data.to_csv(out, index = False)