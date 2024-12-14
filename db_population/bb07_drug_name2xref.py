#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
from pprint import pprint

import pubchempy as pcp
import requests


def run():
    from bad_bac_exercise.models import Drug

    for drug in Drug.objects.all():
        if not drug.is_discrete_structure: continue
        if drug.canonical_smiles is not None: continue
        print(drug.name)
        # Search for the compound by name
        compounds = pcp.get_compounds(drug.name, 'name')

        if (number_of_comps_found := len(compounds)) == 0:
            print(f"No compounds in PubChem found corresponding to {drug.name}")
            # some names do not correspond to a discrete structure, e.g. polymyxin B
            # (In chemistry, the term "discrete structure" refers to distinct, well-defined arrangements
            # of atoms or molecules that can be clearly identified and classified.)

            continue
        elif number_of_comps_found > 1:
            # these are presumably isomers - make sure they all lead to the same canonical smiles
            if any([c.canonical_smiles != compounds[0].canonical_smiles for c in compounds[1:]]):
                print(f"Not all compounds have the same canonical string for {drug.name}")
                for c in compounds:
                    print()
                    print(c.inchikey)
                    print(c.canonical_smiles)
                    print("atom count:", c.heavy_atom_count)
                    print(c.iupac_name)
                    exit()
                continue
        drug.inchi_key = compounds[0].inchikey
        drug.canonical_smiles = compounds[0].canonical_smiles
        drug.save()
        print(compounds[0].inchikey, compounds[0].canonical_smiles)


def test_run():
    # List of antibiotic names
    antibiotic_names = ['Penicillin', 'Amoxicillin', 'Ciprofloxacin', 'Tetracycline']

    # # Dictionary to hold antibiotic names and their InChIKeys
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
    test_run()
