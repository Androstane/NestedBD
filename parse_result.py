# -*- coding: utf-8 -*-
"""
Created on Mon Mar 29 15:14:10 2021

@author: eveli
"""
import ete3
from ete3 import Tree
def get_map(path):
    #read file into list
    with open(path) as f:
        lines = f.readlines()
    L1_index = lines.index('Begin trees;\n')
    L2_index = [i for i,d in enumerate(lines) if d==';\n']
    m = lines[L1_index+2:L2_index[0]]
    f.close()
    dic = {}
    for element in m:
        element = element.replace('\t', '').replace('\n', '').replace(',', '')
        l = element.split(' ')
        dic[l[-2]+':'] =l[-1] 
    return dic

def get_tree(path):
    with open(path) as f:
        lines = f.readlines()
    return lines[-2].split('=')[1]

def replace(tree, map_dict):
    for key in map_dict:
        tree = tree.replace(','+ key, ','+ map_dict[key]+':').replace('('+ key, '('+ map_dict[key]+':')
    tree = tree.replace('leaf', '')
    return tree


def parse_result(path):
    tree = get_tree(path)
    #print(tree)
    map_dict = get_map(path)
    #print(map_dict)
    tree = replace(tree, map_dict)
    return tree
    
 
T = parse_result('beast.1616826147469.trees')
ti = Tree(T)
print(ti)
tt = '((((((83,(123,124)84)45,(((177,178)167,168)95,96)46)23,((71,(((87,88)77,78)73,((131,((163,164)149,(157,158)150)132)75,(169,170)76)74)72)47,((105,(109,(129,130)110)106)65,66)48)24)3,((((133,(183,(197,198)184)134)91,(171,172)92)57,((153,154)85,(181,182)86)58)11,(25,((33,(59,(107,(113,114)108)60)34)27,((119,(173,(179,180)174)120)79,(93,94)80)28)26)12)4)1,((((((((161,162)121,(165,166)122)97,(103,(147,148)104)98)49,50)35,(137,138)36)29,(185,186)30)7,((((145,146)127,128)17,((67,68)19,(((195,196)187,188)81,82)20)18)13,((((117,118)55,56)37,((159,(193,194)160)99,(191,192)100)38)15,((53,(125,126)54)21,(39,((61,(63,(151,152)64)62)51,52)40)22)16)14)8)5,(((89,(101,(139,140)102)90)41,(69,70)42)9,((143,144)31,(((135,(141,142)136)111,(115,116)112)43,((189,190)155,(175,176)156)44)32)10)6)2)0);'
tt = Tree(tt)
print(ti.robinson_foulds(tt)[0])
