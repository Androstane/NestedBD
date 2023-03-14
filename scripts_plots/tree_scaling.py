import dendropy
import numpy as np
import random
#tree= dendropy.Tree.get(path="test.MCC", schema="nexus")
#tree.write(path = 'tt.nw', schema='newick')
import dendropy
import numpy as np
import random
#random.uniform(0,1)


def leaves_under(node):
    L = []
    for leaf in node.leaf_iter():
        L.append(leaf.taxon.label)
    return L

def scale_branch(tree, out_file,dict_file):
    internal_node_count = 0
    internal_node_dict = {}
    tree= dendropy.Tree.get(path=tree, schema="nexus")
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
                #print('root')
                node.scaled_branch = float(node.annotations['length'].value) * float(node.annotations['rate'].value)
                node.rate = 1
        else:
            node.scaled_branch = float(node.annotations['length'].value) * float(node.annotations['rate'].value)
        #print(node.scaled_branch)
        node.edge.length = node.scaled_branch
    for node in tree.preorder_node_iter():
        if node.taxon != None:
            #node is leaf
            pass
        else:
            #node is internal node or root
            L = []
            node.label = 'node' + str(internal_node_count)
            internal_node_dict[node.label] = leaves_under(node)
            internal_node_count += 1
    tree.write(path = out_file, schema='newick')
    with open(dict_file,'w') as out_dict:
          out_dict.write(str(internal_node_dict))
    return True

for rep in ['rep0', 'rep1','rep2','rep3','rep4','rep5,','rep6','rep7','rep8','rep9']:
        tree = 'simulated/20cell_halfevents/' + rep + '/MCC.nex'
        out_file = 'simulated/20cell_halfevents/' + rep + '/MCCtree_scaled.nw'
        dict_file = 'simulated/20cell_halfevents/' + rep + '/bd_mcctree.csv.dict'
        scale_branch(tree, out_file,dict_file)
