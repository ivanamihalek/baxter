#! /usr/bin/env python
from pprint import pprint

import requests
from rdkit import Chem

def to_canonical(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), True)

def get_drug_code(pdb_id, query_smiles):
    # Define the GraphQL query
    query = f"""
    {{
          entry(entry_id: "{pdb_id}") {{
                nonpolymer_entities {{
                    nonpolymer_comp {{
                        chem_comp{{
                             name
                             id
                        }}
                        rcsb_chem_comp_descriptor{{
                            SMILES
                            
                        }}           
                    }}
                }}
          }}
    }}
    """

    print(query)

    # Set up the request headers and URL
    url = "https://data.rcsb.org/graphql"
    headers = {"Content-Type": "application/json"}

    # Send the request to the PDB GraphQL API
    response = requests.post(url, json={"query": query}, headers=headers)

    # Check if the request was successful
    if response.status_code == 200:
        data = response.json()
        if data['data'] is None:
            print("Errors:")
            pprint(data['errors'])
            exit()

        pprint(data['data'])

        nonpoly_entries = data['data']['entry']['nonpolymer_entities']

        pdb_naming = []
        for nonpoly_entry in nonpoly_entries:
            smiles =  nonpoly_entry['nonpolymer_comp']['rcsb_chem_comp_descriptor']['SMILES']
            entry_dict = nonpoly_entry['nonpolymer_comp']['chem_comp']
            drug_name = entry_dict['name']
            pdb_code = entry_dict['id']

            if to_canonical(smiles) == query_smiles:
                pdb_naming.append((drug_name, pdb_code))
    else:
        return f"Error: {response.status_code} - {response.text}"

    return pdb_naming

# Example usage
pdb_id = "3TZF"
inchi_key = "JLKIGFTWXXRPMT-UHFFFAOYSA-N"
smiles = "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"

(drug_name, pdb_code) = get_drug_code(pdb_id, smiles)[0]

print(f"The 3-character code for '{drug_name}' in PDB entry '{pdb_id}' is: {pdb_code}")
