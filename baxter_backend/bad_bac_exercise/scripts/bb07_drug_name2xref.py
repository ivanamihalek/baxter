#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
import pubchempy as pcp
import requests


def run():


    # List of antibiotic names
    antibiotic_names = ['Penicillin', 'Amoxicillin', 'Ciprofloxacin', 'Tetracycline', 'Simocyclinone D8']

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

#######################
if __name__ == "__main__":
    run()
