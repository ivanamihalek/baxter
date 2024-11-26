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
from bad_bac_exercise.models import CARDModel, Gene, AntibioticResMutation, Drug, DrugClass, PDBStructure, Pdb2Gene


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
                # print(alignment.title, hsp.query, hsp.match, hsp.sbjct)
                if hsp.expect > 0: continue
                hit_identifier_fields = alignment.hit_id.split("|")
                if len(hit_identifier_fields) < 2:
                    print(f"hullo, where's my pdb id? {alignment.hit_id}")
                    exit(1)
                pdb_id    = hit_identifier_fields[1]
                pdb_chain =  hit_identifier_fields[2] if len(hit_identifier_fields) > 2 else ""
                try:
                    pdb_struct_entry = PDBStructure.objects.get(pdb_id=pdb_id)
                except:
                    # print(f"{pdb_id} not in db")
                    continue
                pct_identity = int(round((hsp.identities/hsp.align_length*100), 0))
                if pct_identity > 90:
                    junction_table_entry = {
                        "pdb": pdb_struct_entry,
                        "gene": gene_entry,
                        "pdb_chain": pdb_chain,
                        "pct_identity": pct_identity,
                        "gene_seq_aligned": hsp.query,
                        "gene_seq_start": hsp.query_start,
                        "gene_seq_end": hsp.query_end,
                        # there could be all sorts of funny things going on here,
                        # e.g. 9ggqD that has two residues prepended to the protein seq
                        # but then the residue numbering starts from 0, for good measure
                        "pdb_seq_start": hsp.sbjct_start,
                        "pdb_seq_end": hsp.sbjct_end,
                        "pdb_seq_aligned": hsp.sbjct,
                    }
                    (pdb2gene_entry, was_created) = Pdb2Gene.objects.update_or_create(**junction_table_entry)
                    pdb2gene_entry.save()

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
