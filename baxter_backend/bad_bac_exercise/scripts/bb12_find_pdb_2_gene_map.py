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


def run_blast(query_file, output_file):
    # Define parameters
    db = "/storage/databases/pdb/blast/pdb_seqres.fasta"

    print(f"running blastp")
    # Run the blastp command
    subprocess.run(["blastp", "-db", db, "-query", query_file, "-out", output_file, "-outfmt", "5"])

    # todo check success


def parse_blast_results(gene_entry, blast_results_file):
    blast_res_handle = open(blast_results_file)
    blast_records = NCBIXML.parse(blast_res_handle)
    for blast_record in blast_records:
        print("************************")
        print(f"Query: {blast_record.query}")
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:

                if hsp.expect > 0: continue
                pdb_id = alignment.hit_id.split("|")[1].upper()
                try:
                    pdb_struct_entry = PDBStructure.objects.get(pdb_id=pdb_id)
                except:
                    # print(f"{pdb_id} not in db")
                    continue
                identitiy_fraction = hsp.identities/hsp.align_length
                print(pdb_id, f"{identitiy_fraction:.2f}")
                if identitiy_fraction > 0.9:
                    pdb_struct_entry.genes.add(gene_entry)
                    pdb_struct_entry.save()



def run():

    scratch = "/home/ivana/scratch/baxter/blast"

    for gene_entry in Gene.objects.all():
        blast_results_file = f"{scratch}/{gene_entry.name}_blastout.xml"
        if is_nonempty_file(blast_results_file):
            pass
        else:
            query_file  = f"{scratch}/{gene_entry.name}.fasta"
            with open(query_file, "w") as outf:
                print(f"> {gene_entry.name}", file=outf)
                print(gene_entry.protein_seq, file=outf)
            run_blast(query_file, blast_results_file)

        parse_blast_results(gene_entry, blast_results_file)


#######################
if __name__ == "__main__":
    run()
