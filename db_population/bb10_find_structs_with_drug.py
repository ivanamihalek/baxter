#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb08_...
# in that case django will take care of the paths and also check for migrations and such
import os
from sys import argv

import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')
django.setup()

from models.bad_bac_models import Drug, PDBStructure, Pdb2Drug

# note: this is installed with
# pip install rcsb-api
from rcsbapi.search import ChemSimilarityQuery


def find_pdb_by_smiles(smiles) -> list[str]:

    # it looks like PDB cannot be searched by inchikey, for some reason:
    # https://www.rcsb.org/docs/search-and-browse/advanced-search/chemical-similarity-search

    query = ChemSimilarityQuery(query_type="descriptor", descriptor_type="SMILES",  match_type="graph-exact", value=smiles)

    # Execute the query and retrieve results
    try:
        results = query.exec()
    except Exception as e:
        print(f"Problem running chem similarity query for {smiles}: {e}")
        return []
    if not results: return []

    results_dict = results.to_dict()
    if not results_dict: return []

    # results is a dict, with the keys
    # ['query_id', 'result_type', 'total_count', 'result_set']
    pdb_ids = results_dict['result_set']  # the result will be a list of PDB identifiers

    return pdb_ids


def run():
    """
    Updates PDBStructure, adds row to pdb_entry_3_drug through table
    """

    for drug in Drug.objects.all():
        if not drug.is_discrete_structure: continue
        pdbids = find_pdb_by_smiles(drug.canonical_smiles)
        # todo  -check ir therei any peptide ther at all
        # pdb has some NMR stuff with drugs only, and such, see for example 1T5N
        if not pdbids: continue
        print(f"{drug.name}  {pdbids}")
        for pdb_id in pdbids:
            (pdb_entry, was_created) = PDBStructure.objects.update_or_create(pdb_id=pdb_id)
            pdb_entry.save()
            Pdb2Drug.objects.update_or_create(pdb_id=pdb_entry.id, drug_id=drug.id)


def test_run():

    # this is Simocyclinone D8
    # one of the results should be 2Y3P
    smiles = "C[C@@H]1[C@H]([C@@H](C[C@@H](O1)C2=C(C3=C(C=C2)C(=O)[C@]45[C@]([C@@H]3O)(O4)"\
             "[C@H](C[C@]6([C@@]5(C(=O)C=C(C6)C)O)O)O)O)OC(=O)/C=C/C=C/C=C/C=C/C(=O)"\
             "NC7=C(C8=C(C(=C(C=C8)O)Cl)OC7=O)O)OC(=O)C"
    pdbids = find_pdb_by_smiles(smiles)
    print(pdbids)


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
