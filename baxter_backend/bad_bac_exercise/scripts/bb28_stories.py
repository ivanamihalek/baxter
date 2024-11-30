#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import os

from bad_bac_exercise.models import AntibioticResMutation, Pdb2Gene, Pdb2Drug, Pdb2Mutation, PDBStructure
from bad_bac_exercise.models import UCSCAssembly, Gene2UCSCAssembly


def run():

    # "arm" = antibiotic resistance mutation

    # start with pdb, arm pairs where the mutation is close to ligand
    # this is supposed to save us time from investigating arms that

    visited = set()
    for pdb_2_arm_entry in Pdb2Mutation.objects.filter(dist_to_drug__lte=9.0):
        arm_entry = pdb_2_arm_entry.antibio_res_mutation
        dist = pdb_2_arm_entry.dist_to_drug
        if arm_entry.id in visited: continue
        visited.add(arm_entry.id)
        print("***************************************************")
        assmb_entry =  arm_entry.assemblies.first()
        print(pdb_2_arm_entry.pdb.pdb_id, arm_entry.gene.name, arm_entry.mutation, str(round(dist, 1)), assmb_entry.refseq_assembly_id)
        # select region from  assmb_entry.refseq_assembly_id, + 2 decoys  --> task: find problematic species using NCBI Blast
        # toy output from sequencer - the name xxx  --> task: find problematic mutation using Galaxy
        # ---> task: find genomic location from the previous step using UCSC Genome Browser and some calculation - gene and protein location
        # ---> task what is the impact of mutation - use codon translation table Wikipedia
        # ---> is this mutation known?  CARD database -- should get the drug name
        # ---> PDB search for species + gene + drug name

#######################
if __name__ == "__main__":
    run()
