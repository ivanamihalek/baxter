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
from bad_bac_exercise.models import CARDModel, Gene, AntibioticResMutation, Drug, DrugClass


def run_blast(query_file, output_file):
    # Define parameters
    db = "/storage/databases/pdb/blast/pdb_seqres.fasta"

    print(f"running blastp")
    # Run the blastp command
    subprocess.run(["blastp", "-db", db, "-query", query_file, "-out", output_file, "-outfmt", "5"])

    # todo check success


def fetch_small_molecule_info(pdb_id):
    # I could not get the graphQL doodoo to work (tried pypdb)
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    print(response.status_code)  # Print HTTP status code
    response_dict = response.json()
    # pprint(response_dict)
    print(response_dict['rcsb_entry_info']['nonpolymer_bound_components'])
    # get  InChIKey or some such from comoponents.cif or some such
    # store inchi key into drug table that I have from CARD database


def parse_blast_results(blast_results_file):
    blast_res_handle = open(blast_results_file)
    blast_records = NCBIXML.parse(blast_res_handle)
    for blast_record in blast_records:
        print(f"Query: {blast_record.query}")
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:

                if hsp.expect > 0: continue
                pdb_id = alignment.hit_id.split("|")[1].upper()
                print("************************")
                print(pdb_id)
                fetch_small_molecule_info(pdb_id)
                print(f"Hit: {alignment.hit_id}  {alignment.hit_def}")
                print(f"Match: {hsp.identities}, alignment length: {hsp.align_length},   E-value: {hsp.expect:.1e}")
                print(f"alignment length: {hsp.align_length}    gaps: {hsp.gaps}")
                # print(f"Alignment:\n{hsp.sbjct}\n")
        exit()


def run():

    fetch_small_molecule_info("2Y3P")
    exit()
    scratch = "/home/ivana/scratch/baxter/blast"

    for gene_entry in Gene.objects.all():
        print(gene_entry.name)
        blast_results_file = f"{scratch}/{gene_entry.name}_blastout.xml"
        if is_nonempty_file(blast_results_file):
            pass
        else:
            query_file  = f"{scratch}/{gene_entry.name}.fasta"
            with open(query_file, "w") as outf:
                print(f"> {gene_entry.name}", file=outf)
                print(gene_entry.protein_seq, file=outf)
            run_blast(query_file, blast_results_file)

        parse_blast_results(blast_results_file)

        exit()


#######################
if __name__ == "__main__":
    run()
