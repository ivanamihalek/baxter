#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
import requests
import pandas as pd

import subprocess
from Bio.Blast import NCBIXML
from pypdb.clients.data.data_types import DataType, DataFetcher

from .utils import is_nonempty_file
from bad_bac_exercise.models import CARDModel, Gene, AntibioticResMutation, Drug, DrugClass


def run_blast(query_file, output_file):
    # Define parameters
    db = "/storage/databases/pdb/blast/pdb_seqres.fasta"

    print(f"running blastp")
    # Run the blastp command
    subprocess.run(["blastp", "-db", db, "-query", query_file, "-out", output_file, "-outfmt", "5"])

    # todo check success


def fetch_pdb_info(pdb_id):
    # try this maybe? https://github.com/rcsb/py-rcsb-api
    # https://github.com/rcsb/py-rcsb-api
    # see https://data.rcsb.org/rest/v1/schema/entry
    entry = DataFetcher([pdb_id], DataType.ENTRY)
    property = {
        "exptl": ["method", "details"],
        "nonpolymer_comp": ["chem_comp"]
    }
    entry.add_property(property)
    entry.fetch_data()
    print(entry.response)
    #
    # # see https://data.rcsb.org/rest/v1/schema/nonpolymer_entity
    # entry = DataFetcher([pdb_id], DataType.NONPOLYMER_ENTITY_INSTANCE)
    # property = {
    #     "nonpolymer_entities": ["nonpolymer_comp"]
    # }
    # entry.add_property(property)
    # entry.fetch_data()
    # print(entry.response)


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
                fetch_pdb_info(pdb_id)
                print(f"Hit: {alignment.hit_id}  {alignment.hit_def}")
                print(f"Match: {hsp.identities}, alignment length: {hsp.align_length},   E-value: {hsp.expect:.1e}")
                print(f"alignment length: {hsp.align_length}    gaps: {hsp.gaps}")
                # print(f"Alignment:\n{hsp.sbjct}\n")
        exit()


def run():

    fetch_pdb_info("2Y3P")
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
