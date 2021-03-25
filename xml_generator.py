import argparse
import sys


def parse_xml(path_to_samplexml):
    with open('filename') as f:
        lines = f.readlines()
    L1_index = lines.index('<data id="CNAprofile" dataType="integer">')
    L2_index = lines.index('</data>')
    L1 = lines[0:L1_index]
    L2 = lines[L2_index:]
    f.close()
    return L1,L2

def parse_data(path_to_datafile):
    return dic
    



def write_xml(path_to_outxml, data, list1, list2):
    f=open(path_to_outxml,'w')
    f.writelines(list1)
    prof_length = next(iter(data.values()))
    string = "<sequence taxon=\"diploid \">"+ "2"*prof_length + "</sequence>" + "\n"
    f.write(string)
    for key, value in data:
        string = "<sequence taxon=\"" + str(key) + ">\"" + str(value) + "</sequence>" + "\n"
        f.write(string)
    f.writelines(list2)


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
    data_dic = parse_data(data)
    write_xml(path, data_dic, L1, L2)
    