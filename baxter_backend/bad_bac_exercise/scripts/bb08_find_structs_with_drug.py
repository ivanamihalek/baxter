#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb08_...
# in that case django will take care of the paths and also check for migrations and such

from rcsbsearchapi.search import ChemSimilarityQuery


def find_pdb_by_smiles(smiles) -> list[str]:

    # Create a query for chemical attributes using InChIKey
    query = ChemSimilarityQuery(query_type="descriptor", descriptor_type="SMILES",  match_type="graph-exact", value=smiles)

    # Execute the query and retrieve results
    results = query.exec().to_dict()
    # results is a dict, with the keys
    # ['query_id', 'result_type', 'total_count', 'result_set']
    pdb_ids = results['result_set']  # the result will be a list of PDB identifiers

    return pdb_ids


def run():
    # this is Simocyclinone D8
    # one of the results should be 2Y3P
    smiles = "C[C@@H]1[C@H]([C@@H](C[C@@H](O1)C2=C(C3=C(C=C2)C(=O)[C@]45[C@]([C@@H]3O)(O4)[C@H](C[C@]6([C@@]5(C(=O)C=C(C6)C)O)O)O)O)OC(=O)/C=C/C=C/C=C/C=C/C(=O)NC7=C(C8=C(C(=C(C=C8)O)Cl)OC7=O)O)OC(=O)C"  # this is Simocyclinone D8, and should be present in 2Y3P
    pdbids = find_pdb_by_smiles(smiles)
    print(pdbids)


#######################
if __name__ == "__main__":
    run()
