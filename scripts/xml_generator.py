import argparse
import pandas as pd


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
    d = cnp.to_dict('list')   # exchange again
    return d

def cap(l):
    L = []
    for i in l:
        if i > 29:
            L.append(29)
        else:
            L.append(i)
    return L

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
    L1, L2 = parse_xml("test_data/BD.xml")
    data_dic = parse(data)
    write_xml(path, data_dic, L1, L2)
    