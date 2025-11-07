#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import os
from sys import argv

import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')
django.setup()

from models.bad_bac_models import PDBStructure, Drug, Pdb2Drug
from pprint import pprint

import numpy as np
import requests

from rdkit import Chem

def to_canonical(smiles):
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), True)

# InchiKey is almost unique (gosh)


def get_drug_code_by_smiles (pdb_id, query_smiles) -> tuple[str, str]:
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
        if nonpoly_entries is None:
            print(f"\tnonpoly_entries is None  for {pdb_id}")
            return "", ""


        pdb_naming = []
        for nonpoly_entry in nonpoly_entries:
            if nonpoly_entry is None:
                print(f"\t nonpoly_entry is None for {pdb_id}")
                continue

            if 'nonpolymer_comp' not in nonpoly_entry:
                print(f"\t nonpolymer_comp  not found for {pdb_id}")
                continue

            if nonpoly_entry['nonpolymer_comp'] is None:
                print(f"nonpoly_entry['nonpolymer_comp'] is None for {pdb_id}")
                continue

            if 'rcsb_chem_comp_descriptor' not in nonpoly_entry['nonpolymer_comp']:
                print(f"\t  rcsb_chem_comp_descriptor not found for {pdb_id}")
                continue

            if nonpoly_entry['nonpolymer_comp']['rcsb_chem_comp_descriptor'] is None:
                print(f"nonpoly_entry['nonpolymer_comp']['rcsb_chem_comp_descriptor'] is None for {pdb_id}")
                continue

            if 'SMILES' not in  nonpoly_entry['nonpolymer_comp']['rcsb_chem_comp_descriptor']:
                print(f"\t SMILES not found for {pdb_id}  (expected: {query_smiles})")
                continue

            smiles =  to_canonical(nonpoly_entry['nonpolymer_comp']['rcsb_chem_comp_descriptor']['SMILES'])
            entry_dict = nonpoly_entry['nonpolymer_comp']['chem_comp']
            drug_name = entry_dict['name']
            pdb_code = entry_dict['id']
            if smiles == query_smiles:
                pdb_naming.append((drug_name, pdb_code))

    else:
        print(f"Error: {response.status_code} - {response.text}")
        exit()

    if len(pdb_naming) == 0:
        print(f"No 3-character identifiers found in {pdb_id} for {query_smiles}")
        return "", ""

    if len(pdb_naming) > 1:
        print(f"Multiple 3-character identifiers found in {pdb_id} for {query_smiles}:\n{pdb_naming}")
        exit()

    (drug_name, pdb_code) = pdb_naming[0]
    if len(pdb_code) > 3:
        print(f"Pdb res id longer than 3 characters?: {pdb_code}")

    return drug_name, pdb_code


def run():

    for pdb_2_drug_entry in Pdb2Drug.objects.all():
        # if pdb_2_drug_entry.drug_residue_name is not None and pdb_2_drug_entry.drug_name_in_pdb is not None: continue
        pdb_id = pdb_2_drug_entry.pdb_id
        drug_id = pdb_2_drug_entry.drug_id
        pdb_entry  = PDBStructure.objects.get(pk=pdb_id)
        drug_entry = Drug.objects.get(pk=drug_id)
        (drug_name, pdb_code) = get_drug_code_by_smiles(pdb_entry.pdb_id, drug_entry.canonical_smiles)
        print(pdb_id, drug_entry.name, drug_name, pdb_code)
        if drug_name and pdb_code:
            pdb_2_drug_entry.drug_name_in_pdb  = drug_name
            pdb_2_drug_entry.drug_residue_name = pdb_code
            pdb_2_drug_entry.save()


def test_run():
    pdb_id = "3TZF"
    inchi_key = "JLKIGFTWXXRPMT-UHFFFAOYSA-N"
    smiles = "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"
    (drug_name, pdb_code) = get_drug_code_by_smiles(pdb_id, smiles)
    print(f"The 3-character code for '{drug_name}' in PDB entry '{pdb_id}' is: {pdb_code}")



#######################
def main():
    if len(argv) < 2:
        exit("Tell me what to do - test or run?")
    if argv[1] == 'run':
        run()
    elif argv[1] == 'test':
        test_run()
    else:
        print(f"I don't know what is '{argv[1]}'. What should I do - test or run?")

#######################
if __name__ == "__main__":
    main()
