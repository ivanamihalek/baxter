#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb08_...
# in that case django will take care of the paths and also check for migrations and such

from rcsbsearchapi.search import ChemSimilarityQuery


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

    from bad_bac_exercise.models import Drug, PDBStructure

    for drug in Drug.objects.all():
        if not drug.is_discrete_structure: continue
        pdbids = find_pdb_by_smiles(drug.canonical_smiles)
        if not pdbids: continue
        print(f"{drug.name}  {pdbids}")
        for pdb_id in pdbids:
            (pdb_entry, was_created) = PDBStructure.objects.update_or_create(pdb_id=pdb_id)
            pdb_entry.drugs.add(drug)
            pdb_entry.save()


def run_test():
    # it should be possible to search the

    # this is Simocyclinone D8
    # one of the results should be 2Y3P
    smiles = "C[C@@H]1[C@H]([C@@H](C[C@@H](O1)C2=C(C3=C(C=C2)C(=O)[C@]45[C@]([C@@H]3O)(O4)[C@H](C[C@]6([C@@]5(C(=O)C=C(C6)C)O)O)O)O)OC(=O)/C=C/C=C/C=C/C=C/C(=O)NC7=C(C8=C(C(=C(C=C8)O)Cl)OC7=O)O)OC(=O)C"  # this is Simocyclinone D8, and should be present in 2Y3P
    pdbids = find_pdb_by_smiles(smiles)
    print(pdbids)


#######################
if __name__ == "__main__":
    run_test()
