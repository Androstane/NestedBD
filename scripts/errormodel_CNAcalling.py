# -*- coding: utf-8 -*-
import scipy
from scipy.stats import poisson
from scipy.stats import nbinom
from scipy.stats import norm
from scipy.special import binom
from scipy.special import loggamma
import math
from collections import defaultdict
from ete3 import Tree
import pandas as pd
import argparse
import sys

class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)


def bd_prob(ancestor,child,t):
    p = 0
    #coef = binom(i,j)
    if t <= 1.0:
        for m in range(1, min(child, ancestor)+1):
            p_m = binom(ancestor, m) * binom(child - 1, m - 1) * pow(t, child + ancestor -2*m)
            p += p_m
        p = p * pow(1/(t+1), child + ancestor)
    else:
        for m in range(1, min(child, ancestor)+1):
            p_m = binom(ancestor, m) * binom(child - 1, m - 1) * pow(t, -2*m)
            p += p_m
        p = p * pow(t/ (t+1), child + ancestor)
    return p

#i: child
#j: ancestor
#compute P(i|j, t) under the birth-death model
def transition_prob(t):
    d = defaultdict(int)
    for i in range(0, 10):
        for j in range(0, 10):
            index = i * 10 + j
            if i != 0 and j == 0:
                p = 0
            elif i == 0 and j == 0:
                p = 1
            elif i == 0:
                p =  pow(t, j - 1) / pow((1 + t), j + 1)
            elif j == 1:
                p = pow(t, i - 1)/pow((1 + t), i + 1)
            else:
                p = bd_prob(j, i, t)
            d[index] = p
    return d

#i: child
#j: ancestor
#compute P(i|j, t) under the birth-death model
def transition_prob_single(child, ancestor, t):
    i, j = child, ancestor
    if i != 0 and j == 0:
        p = 0
    elif i == 0 and j == 0:
        p = 1
    elif i == 0:
        p =  pow(t, j - 1) / pow((1 + t), j + 1)
    elif j == 1:
        p = pow(t, i - 1)/pow((1 + t), i + 1)
    else:
        p = bd_prob(j, i, t)
    return p

def poisson_prob(rc, expected_rc, ploidy):
    d = defaultdict(int)
    for i in range(0, 10):
        mu = max(5, expected_rc * i/ploidy)
        #print(mu)
        x = poisson.pmf(rc, mu)
        d[i] = x
    d = list(d.values())
    return [float(i)/sum(d) for i in d]

def nb_prob(rc, expected_rc, ploidy, theta):
    d = defaultdict(int)
    for i in range(0, 10):
        mu = max(5, expected_rc * i/ploidy)
        #print(mu)
        #x = nbinom.pmf(rc, n=theta , p = mu/(mu + theta))
        val1 = loggamma(rc + theta) - loggamma(rc+1) - loggamma(theta)
        val2 = theta * math.log(theta/(theta+mu))
        val3 = rc*math.log(mu/(theta+mu))
        x = math.exp(val1 + val2 + val3)
        d[i] = x
    d = list(d.values())
    #print(d)
    return [float(i)/sum(d) for i in d]


def normal_prob(n, var):
    d = defaultdict(int)
    for i in range(0, 10):
        d[i] = norm.pdf(n, loc = i, scale = math.sqrt(var))
    d = list(d.values())
    return [float(i)/sum(d) for i in d]

def calc_partials(L1, L2, t):
    L_node, C_node = [], []
    prob_mat = transition_prob(t)
    #i:ancestor
    #j:child
    for i in range(0, 10):
        index = 0
        max_ = 0
        for j in range(0, 10):
            val =  prob_mat[j*10 + i] * L1[j] * L2[j]
            if val >= max_:
                index = j
                max_ = val
        L_node.append(max_)
        C_node.append(index)
    return L_node, C_node


def calc_leafpartials(error_model,t, rc, expected_rc):
    prob_mat = transition_prob(t)
    L_leaf, C_leaf = [], []
    if error_model == 'p':
        error_partial = poisson_prob(rc, expected_rc, ploidy)
    elif error_model == 'nb':
        error_partial = nb_prob(rc, expected_rc, ploidy, theta)
    #print(error_partial)
    elif error_model == 'normal':
        error_partial = normal_prob(rc, var)
    #print(error_partial)
    for i in range(0, 10):
        index, max_ = 0,0
        for j in range(0, 10):
            val = prob_mat[j*10 + i] * error_partial[j]
            if val >= max_:
                index = j
                max_ = val
        L_leaf.append(max_)
        C_leaf.append(index)
    return L_leaf, C_leaf


#nbin: total number of bin
#n: current bin computed by traverse_up
def traverse_up(tree, nbin, n, error_model):
    #print(tree.name)
    t = tree.dist
    if tree.is_leaf():
        #print(n)
        rc = rc_matrix[tree.name][n]
        #print(rc)
        expected_rc = sum_dict[tree.name]/nbin
        L[tree.name], C[tree.name] = calc_leafpartials(error_model,t, rc, expected_rc)
    #both children were visited
    #compute
    else:
        if tree.children[0].name not in L:
            traverse_up(tree.children[0],nbin, n, error_model)
        if tree.children[1].name not in L:
            traverse_up(tree.children[1], nbin, n, error_model)
        L[tree.name], C[tree.name] = calc_partials(L[tree.children[0].name], L[tree.children[1].name],t)
#L, C should contain all internal node as key and an array of length 10 with value of Li and Ci as value.

def assign(tree,orig_dist, R):
    #print(tree.name)
    #compute root partials
    root_partial = []
    if tree.is_root():
        L_node = L[tree.name]
        left = L[tree.children[0].name]
        right = L[tree.children[1].name]
        for i in range(0, 9):
            if i == 0:
                val =  pow(orig_dist, 1) / pow((1 + orig_dist), 2);
            else:
                #val = transition_prob(i, 2, orig_dist)
                val = transition_prob_single(i, 2, orig_dist)
            root_partial.append(val)
        max_ = 0
        max_index = 0
        #print(root_partial, left, right)
        for j in range(0, 9):
            temp = root_partial[j] * left[j] * right[j]
            if temp > max_:
                max_ = temp
                max_index = j
        R[tree.name] = max_index
        if tree.children[0] != None:
            assign(tree.children[0], orig_dist, R)
        if tree.children[1] != None:
            assign(tree.children[1], orig_dist, R)
    else:
        name = tree.name
        val = R[tree.up.name]
        R[name] = C[name][val]
        if tree.children:
            #if tree.children[0] != None:
            assign(tree.children[0], orig_dist, R)
            #if tree.children[1] != None:
            assign(tree.children[1], orig_dist, R)

if __name__ == "__main__":
        #orig_stdout = sys.stdout
        sys.stdout = Unbuffered(sys.stdout)
        #print('here')
        ap = argparse.ArgumentParser()
        ap.add_argument("-in","--path to input",required=True, help="Path to the input folder containing the copy number profiles ")
        ap.add_argument("-od","--origin_dist of tree",required=True, help="origin_dist of the NestedBD inferred tree ")
        ap.add_argument("var","--variance of the normal error model",required=True, help="variance of the normal error model")
        args = vars(ap.parse_args())
        if args['path to input']!=None:
                path = args['path to input']
        if args['origin_dist of tree']!=None:
                path = args['origin_dist of tree']
        if args['variance of the normal error model']!=None:
                var = args['variance of the normal error model']
        #f = open(path + '/infer_prof_out.txt', 'w+')
        #sys.stdout = f
        #out = 'BD_output/simulated/100celltree/'+ ext[ind]
        #path = 'BD_data/100celltree/' + ext[ind]
        tree_file = path + '/MCCtree/mcc.nw'
        #print(tree_file)
        rc_matrix = pd.read_csv(path + '/Ginkgo.cnp', sep = '\t')
        rc_matrix.drop(columns = ['CHR', 'START','END'], inplace = True)
        rc_matrix = rc_matrix.iloc[: , 1:]
        #tree = Tree(tree_file)
        tree = Tree(tree_file)
        i = 0
        for node in tree.traverse():
            if node.name == '':
                node.name = 'node' + str(i)
                #print(node.name)
                i += 1
        total_bin = rc_matrix.shape[0]
        orig_dist = orig_dists[ind]
        error_model = 'normal'
        output = pd.DataFrame(columns = rc_matrix.columns)
        sum_dict = rc_matrix.sum().to_dict()
        for n in range(0, total_bin):
            L = defaultdict(list)
            C = defaultdict(list)
            R = {}
            traverse_up(tree,total_bin, n, error_model)
            assign(tree, orig_dist, R)
            leaf_R = {k:v for k, v in R.items() if k.startswith('leaf')}
            #output.loc[0, leaf_R.keys()] = leaf_R.values()
            output = output.append(leaf_R, ignore_index = True)
        output.to_csv(path + '/cnp.csv', index = False)
        #sys.stdout = orig_stdout
        #f.close()
