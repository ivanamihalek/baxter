#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
import pubchempy as pcp
import requests

from rcsbsearchapi.search import TextQuery, AttributeQuery, Attr

def find_pdb_by_inchikey(inchi_key):

    # Create a query for chemical attributes using InChIKey
    query = AttributeQuery(attribute="rcsb_chem_comp.id", operator="exact_match", value=inchi_key)

    # Execute the query and retrieve results
    results = query()

    # Extract and return PDB IDs from results
    pdb_ids = [entry['identifier'] for entry in results]
    return pdb_ids


def run():

    inchi_key = "PLEGMCYXNQPJNV-JIEXBSPMSA-N"  # this is Simocyclinone D8, and should be present in 2Y3P
    find_pdb_by_inchikey(inchi_key)
    exit()

    # List of antibiotic names
    antibiotic_names = ['Penicillin', 'Amoxicillin', 'Ciprofloxacin', 'Tetracycline']

    # Dictionary to hold antibiotic names and their InChIKeys
    inchi_keys = {}
    smiles = {}
    # Loop through each antibiotic name
    for name in antibiotic_names:
        # Search for the compound by name
        compounds = pcp.get_compounds(name, 'name')

        # If compounds were found, retrieve the InChIKey
        if compounds:
            inchi_keys[name] = compounds[0].inchikey  # Get the first compound's InChIKey
            smiles[name] =  compounds[0].canonical_smiles

    # Output the results
    for name, inchi_key in inchi_keys.items():
        print("************************")
        print(f"{name}: {inchi_key}  {smiles[name]}")
        find_pdb_by_inchikey(inchi_key)
        print()

#######################
if __name__ == "__main__":
    run()
