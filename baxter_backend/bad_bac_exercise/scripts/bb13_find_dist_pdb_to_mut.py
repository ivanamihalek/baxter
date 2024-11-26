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
from bad_bac_exercise.models import AntibioticResMutation, Pdb2Gene


def run():

    scratch = "/home/ivana/scratch/baxter/blast"

    # we fil calculate the distance from the mutated residue to the nearest drug, but only for those
    # pdbs that map to one of our genes of interest
    distinct_genes = set()
    total_candidates = 0
    for abr_mutation in AntibioticResMutation.objects.all():
        mutation_pos = int(abr_mutation.mutation[1:-1])
        mutation_from = abr_mutation.mutation[0]
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

        outstr = f"\n{abr_mutation.mutation}, {abr_mutation.gene.name}\n"
        mutation_mapped = False
        for p, d in mapped_pdb_drug_pairs:
            for mapping_entry in Pdb2Gene.objects.filter(pdb=p, gene=abr_mutation.gene):
                if not (mapping_entry.gene_seq_start <= mutation_pos <= mapping_entry.gene_seq_end): continue

                # the numbering may not match for some legit reason, but
                # we are not going to deal with this now
                pos_on_gene_seq = mutation_pos - mapping_entry.gene_seq_start
                gene_aa = mapping_entry.gene_seq_aligned[pos_on_gene_seq]
                if mutation_from != gene_aa: continue

                pos_on_pdb_seq = mutation_pos - mapping_entry.pdb_seq_start
                pdb_aa = mapping_entry.pdb_seq_aligned[pos_on_pdb_seq]
                if mutation_from != pdb_aa: continue

                mutation_mapped = True
                chain_id =  f"{p.pdb_id}{mapping_entry.pdb_chain}"
                outstr += f"\t {chain_id}, {mapping_entry.gene_seq_start}, {mapping_entry.gene_seq_end}  "
                outstr += f"{mutation_from}  {gene_aa}  {pdb_aa}\n"
                distinct_genes.add(abr_mutation.gene.name)

        if mutation_mapped:
            total_candidates += 1
            print(outstr)
    print(distinct_genes)
    print(total_candidates)
#######################
if __name__ == "__main__":
    run()
