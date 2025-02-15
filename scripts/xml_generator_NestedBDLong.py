import os
import argparse
import pandas as pd

def convert_dict_to_xml(dict_seqs, append_str, sample_rate=20):
    xml_output = ""
    for name, sequence in dict_seqs.items():
        full_name = name + append_str
        # Convert sequence elements to strings, cap at 9, and sample the sequence
        sequence = [str(min(x, 9)) for x in sequence]
        sampled_sequence = sequence[::sample_rate]
        sequence_formatted = ','.join(sampled_sequence)
        xml_output += f'<sequence id="{full_name}" taxon="{full_name}">\n'
        xml_output += f'    {sequence_formatted}\n'
        xml_output += '</sequence>\n'
    return xml_output

def generate_taxon_values(dict_seqs, taxon_value, append_str):
    result = ""
    for name, sequence in dict_seqs.items():
        # Append the taxon value to the taxon name
        full_name = name + append_str
        result += f'    {full_name}={taxon_value},\n'
    return result.rstrip(',\n')

def update_xml_file_SA(xml_path, new_data_xml, new_trait_xml, output_path=None):
    # Read the original XML content from the file
    with open(xml_path, 'r') as file:
        xml_content = file.read()

    # Replace the old data section
    start_tag = '<data id="bears" name="alignment" dataType="integer">'
    end_tag = "</data>"
    start_index = xml_content.find(start_tag)
    end_index = xml_content.find(end_tag, start_index) + len(end_tag)
    xml_content = (
        xml_content[:start_index] +
        start_tag + '\n' + new_data_xml + end_tag + '\n' +
        xml_content[end_index:]
    )

    # Replace the old trait section
    start_tag = '<trait id="dateTrait.t:bears" spec="beast.evolution.tree.TraitSet" traitname="date-backward">'
    end_tag = '<taxa id="TaxonSet.bears" spec="TaxonSet">'
    start_index = xml_content.find(start_tag)
    end_index = xml_content.find(end_tag, start_index) + len(end_tag)
    xml_content = (
        xml_content[:start_index] +
        start_tag + '\n' + new_trait_xml + end_tag + '\n' +
        xml_content[end_index:]
    )

    # Write the modified XML to file
    if output_path is None:
        output_path = 'modified_output.xml'
    with open(output_path, 'w') as file:
        file.write(xml_content)

def generate_data_time(output_path, ntimepoints=3):
    data, trait = '', ''
    if ntimepoints == 3:
        time_for_SA = [0.0, 1.0, 3.0]
        name_for_SA = ['post', 'mid', 'pre']
    else:
        time_for_SA = [0.0, 1.0, 2.0, 3.0]
        name_for_SA = ['post', 'mid', 'mid1', 'pre']
    for i in range(ntimepoints):
        file_path = os.path.join(output_path, 'sampled_profiles_' + name_for_SA[i] + '.cnp')
        sampled_node_df = pd.read_csv(file_path, sep='\t')
        # Assume the columns 'chrom', 'start', and 'end' are not part of the node list
        node_list = sampled_node_df.columns.difference(['chrom', 'start', 'end']).tolist()
        d = {node: sampled_node_df[node].tolist() for node in node_list}
        data += convert_dict_to_xml(d, '_' + name_for_SA[i])
        trait += generate_taxon_values(d, time_for_SA[i], '_' + name_for_SA[i]) + '\n'
    return data, trait

def main():
    parser = argparse.ArgumentParser(
        description="Update an XML file with sampled data and trait values based on timepoints."
    )
    parser.add_argument("--output_dir", required=True,
                        help="Output directory where the sampled profile files (.cnp) are located")
    parser.add_argument("--xml_template", required=True,
                        help="Path to the XML template file")
    parser.add_argument("--output_xml", required=True,
                        help="Path for the output (modified) XML file")
    parser.add_argument("--ntimepoints", type=int, default=3,
                        help="Number of timepoints in the dataset (e.g., 3 or 4)")
    args = parser.parse_args()

    # Ensure the output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Generate data and trait strings from the sampled profile files
    data, trait = generate_data_time(args.output_dir, ntimepoints=args.ntimepoints)

    # Update the XML template with the new data and trait values
    update_xml_file_SA(args.xml_template, data, trait, output_path=args.output_xml)
    print(f"XML file written to {args.output_xml}")

if __name__ == '__main__':
    main()
