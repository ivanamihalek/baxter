#! /usr/bin/env python
from pprint import pprint

import requests

def get_drug_code(pdb_id, query_drug_name):
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

        # pprint(data['data'])

        nonpoly_entries = data['data']['entry']['nonpolymer_entities']

        pdb_codes = []
        for nonpoly_entry in nonpoly_entries:
            entry_dict = nonpoly_entry['nonpolymer_comp']['chem_comp']
            drug_name = entry_dict['name']
            pdb_code = entry_dict['id']
            if drug_name.lower() == query_drug_name.lower():
                pdb_codes.append(pdb_code)
    else:
        return f"Error: {response.status_code} - {response.text}"

    return pdb_codes

# Example usage
pdb_id = "3TZF"
drug_name = "sulfamethoxazole"
drug_codes = get_drug_code(pdb_id, drug_name)

print(f"The 3-character code for '{drug_name}' in PDB entry '{pdb_id}' is: {drug_codes}")
