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

        print()
        print(abr_mutation.mutation, abr_mutation.gene.name, [d.name for d in  drugs_affected])

#######################
if __name__ == "__main__":
    run()
