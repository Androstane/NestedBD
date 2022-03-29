import argparse
import pandas as pd


import argparse
import pandas as pd
import numpy as np
def parse_xml(path_to_samplexml):
    with open(path_to_samplexml) as f:
        lines = f.readlines()
    L1_index = lines.index('<data id="CNAprofile" dataType="integer">\n')
    L2_index = lines.index('</data>\n')
    L1 = lines[0:L1_index+1]
    L2 = lines[L2_index:]
    f.close()
    return L1,L2

def parse(path):
    cnp = pd.read_csv(path, sep = '\t')
    cnp.drop(cnp.columns[[0,1,2]], axis = 1, inplace = True)
    cnp.columns = cnp.columns.str.replace(' ', '')
    bins_array = np.arange(0, cnp.shape[0], step = 20)
    cnp = cnp.iloc[cnp.index[bins_array]]
    cnp = cnp.astype(int)
    print("number of bins:", cnp.shape[0])
    d = cnp.to_dict('list')
    return d

def cap(l):
    L = []
    for i in l:
        if i > 9:
            L.append(9)
        else:
            L.append(i)
    return L
def inplace_change(filename, old_string, new_string):
    # Safely read the input filename using 'with'
    with open(filename) as f:
        s = f.read()
        if old_string not in s:
            print('"{old_string}" not found in {filename}.'.format(**locals()))
            return
    # Safely write the changed content, if found in the file
    with open(filename, 'w') as f:
        print('Changing "{old_string}" to "{new_string}" in {filename}'.format(**locals()))
        s = s.replace(old_string, new_string)
        f.write(s)

def write_xml(path_to_outxml, data, list1, list2):
    f=open(path_to_outxml,'w')
    f.writelines(list1)
    #prof_length = len(next(iter(data.values())))
    #string = "<sequence taxon=\"diploid \">"+ "2"*prof_length + "</sequence>" + "\n"
    #f.write(string)
    for key in data:
        string = "<sequence taxon=\"" + str(key) + "\">"+ "\n" + str(','.join(map(str, cap(data[key]))))+ "\n" + "</sequence>" + "\n"
        f.write(string)
    f.writelines(list2)
    f.close()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-data","--path to file contain dataset",required=True, help="Path to data")
    ap.add_argument("-out","--path to output xml",required=True, help="Path to output xml ")
    args = vars(ap.parse_args())
    if args['path to file contain dataset']!=None:
        data = args['path to file contain dataset']
    if args['path to output xml']!=None:
        path = args['path to output xml']
    L1, L2 = parse_xml("BD.xml")
    data_dic = parse(data)
    write_xml(path, data_dic, L1, L2)
