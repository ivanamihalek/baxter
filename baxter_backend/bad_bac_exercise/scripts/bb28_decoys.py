#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import random
from os import listdir

from Bio import Entrez

from bad_bac_exercise.models import AntibioticResMutation, Pdb2Gene, Pdb2Drug, Pdb2Mutation, PDBStructure, Decoy
from bad_bac_exercise.models import UCSCAssembly, Gene2UCSCAssembly

# todo
# decoy list
# identifier, source, species name, sequence


def create_phage_decoy(seqs_dir):
    random_fasta_file = random.choice(list(listdir(seqs_dir)))
    print(random_fasta_file)
    with open(f"{seqs_dir}/{random_fasta_file}") as inf:
       seq = "".join(inf.read().replace(" ", "").split("\n")[1:])

    sample_length = random.choice(range(800, 1800))
    sample_start  = random.choice(range(1, len(seq)-1800)) - 1
    if sample_start < 0:
        return "seq too short"
    sample = seq[sample_start:sample_start+sample_length]
    if "N" in sample:
        return "sequencing error present"
    # otherwise store
    print(sample_start, sample_length)
    print(sample)
    decoy_table_entry = {
        "source": Decoy.Source.translate("phages"),
        "identifier": random_fasta_file.replace(".fasta", ''),
        "dna_seq": sample,
    }
    Decoy(**decoy_table_entry).save()

    return "ok"


def run():
    phagi_fasta_dir  = "/storage/databases/phagescope/RefSeq"  # one sequence per data file
    human_microbiome = "/storage/databases/ncbi/human_microbiome"
    plant_genomes    = "/storage/databases/ncbi/plant_genomes"
    for _ in range(10):
        create_phage_decoy(phagi_fasta_dir)
