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

def generate_data_time(output_dir, names_for_SA, times_for_SA):
    data, trait = '', ''
    ntimepoints = len(names_for_SA)
    for i in range(ntimepoints):
        file_path = os.path.join(output_dir, 'sampled_profiles_' + names_for_SA[i] + '.cnp')
        try:
            sampled_node_df = pd.read_csv(file_path, sep='\t')
        except Exception as e:
            raise IOError(f"Error reading file '{file_path}': {e}")
        # Assume the columns 'chrom', 'start', and 'end' are not part of the node list
        node_list = sampled_node_df.columns.difference(['chrom', 'start', 'end']).tolist()
        d = {node: sampled_node_df[node].tolist() for node in node_list}
        data += convert_dict_to_xml(d, '_' + names_for_SA[i])
        trait += generate_taxon_values(d, times_for_SA[i], '_' + names_for_SA[i]) + '\n'
    return data, trait

def main():
    parser = argparse.ArgumentParser(
        description="Update an XML file with sampled data and trait values based on SA timepoints."
    )
    parser.add_argument("--output_dir", required=True,
                        help="Output directory where the sampled profile files (.cnp) are located")
    parser.add_argument("--xml_template", required=True,
                        help="Path to the XML template file")
    parser.add_argument("--output_xml", required=True,
                        help="Path for the output (modified) XML file")
    parser.add_argument("--ntimepoints", type=int, default=3,
                        help="Number of timepoints in the dataset (e.g., 3 or 4)")
    parser.add_argument("--names_of_SA", default=None,
                        help="Comma-separated list of names for SA segments (e.g., 'post,mid,pre')")
    parser.add_argument("--times_of_SA", default=None,
                        help="Comma-separated list of times for SA segments (e.g., '0.0,1.0,3.0')")
    args = parser.parse_args()

    # Ensure the output directory exists
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Process names_of_SA and times_of_SA command line inputs
    if args.names_of_SA:
        names_for_SA = [name.strip() for name in args.names_of_SA.split(',')]
    else:
        if args.ntimepoints == 3:
            names_for_SA = ['post', 'mid', 'pre']
        elif args.ntimepoints == 4:
            names_for_SA = ['post', 'mid', 'mid1', 'pre']
        else:
            names_for_SA = [f"SA{i}" for i in range(args.ntimepoints)]

    if args.times_of_SA:
        try:
            times_for_SA = [float(t.strip()) for t in args.times_of_SA.split(',')]
        except ValueError as e:
            raise ValueError("Times of SA must be numbers separated by commas.") from e
    else:
        if args.ntimepoints == 3:
            times_for_SA = [0.0, 1.0, 3.0]
        elif args.ntimepoints == 4:
            times_for_SA = [0.0, 1.0, 2.0, 3.0]
        else:
            times_for_SA = [float(i) for i in range(args.ntimepoints)]

    if len(names_for_SA) != len(times_for_SA):
        raise ValueError("The number of SA names and SA times must be the same.")

    # Generate data and trait strings from the sampled profile files
    data, trait = generate_data_time(args.output_dir, names_for_SA, times_for_SA)

    # Update the XML template with the new data and trait values
    update_xml_file_SA(args.xml_template, data, trait, output_path=args.output_xml)
    print(f"XML file written to {args.output_xml}")

if __name__ == '__main__':
    main()
