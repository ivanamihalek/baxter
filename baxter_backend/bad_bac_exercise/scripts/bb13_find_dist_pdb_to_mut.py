#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
import json
from pprint import pprint

import requests
import pandas as pd

import subprocess
from Bio.Blast import NCBIXML

from .utils import is_nonempty_file
from bad_bac_exercise.models import CARDModel, Gene, AntibioticResMutation, Drug, DrugClass, PDBStructure



def run():

    scratch = "/home/ivana/scratch/baxter/blast"

    # we fil calculate the distance from the mutated residue to the nearest drug, but only for those
    # pdbs that map to one of our genes of interest
    for abr_mutation in AntibioticResMutation.objects.all():
        drugs_affected = list(abr_mutation.drugs_affected.all())
        if len(drugs_affected) == 0: continue

        pdb_structures = list(abr_mutation.gene.pdbstructure_set.all())
        if len(pdb_structures) == 0: continue

        # find drug<->pdb pairs
        mapped_pdb_drug_pairs = []
        for pdb_entry in pdb_structures:
            drugs_in_pdb = list(pdb_entry.drugs.all())
            for drug_entry in drugs_affected:
                if drug_entry in drugs_in_pdb:
                    mapped_pdb_drug_pairs.append((pdb_entry, drug_entry))
        if not mapped_pdb_drug_pairs: continue

        print()
        print(abr_mutation.mutation, abr_mutation.gene.name)
        for p, d in mapped_pdb_drug_pairs:
            print("\t", p.pdb_id, d.name)

#######################
if __name__ == "__main__":
    run()
