#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
import os
import random
from glob import glob
from os import listdir

from Bio import Entrez

from bad_bac_exercise.models import Decoy


def select_from_single_seq_fasta(filepath):

    with open(filepath) as inf:
        seq = "".join(inf.read().replace(" ", "").split("\n")[1:])
    if not seq:
        return "Error: something went wrong"
    if len(seq) < 5000:
        return "Error: seq too short"
    sample_length = random.choice(range(800, 1800))
    sample_start  = random.choice(range(1, len(seq)-1800)) - 1
    if sample_start < 0:
        return "Error: seq too short"
    sample = seq[sample_start:sample_start+sample_length]
    if "N" in sample:
        return "Error: sequencing error present"
    print(sample_start, sample_length)
    print(sample)
    return sample


def create_phage_decoy(seqs_dir):
    random_fasta_file = random.choice(list(listdir(seqs_dir)))
    print(random_fasta_file)
    sample =  select_from_single_seq_fasta(f"{seqs_dir}/{random_fasta_file}")
    if sample[:3] == "Err": return sample

    # otherwise store
    decoy_table_entry = {
        "source": Decoy.Source.translate("phages"),
        "identifier": random_fasta_file.replace(".fasta", ''),
        "dna_seq": sample,
    }
    Decoy(**decoy_table_entry).save()

    return "ok"


def create_human_microbiome_decoy(seqs_dir):
    i_came_from = os.getcwd()
    random_species =  random.choice(list(listdir(seqs_dir)))
    os.chdir(f"{seqs_dir}/{random_species}")
    print(random_species)
    random_fasta_file = random.choice(list(glob("*.fna")))
    print(random_fasta_file)
    sample = select_from_single_seq_fasta(random_fasta_file)
    if sample[:3] == "Err": return sample
    # otherwise store
    decoy_table_entry = {
        "source": Decoy.Source.translate("human_microbiome"),
        "identifier": random_fasta_file.replace(".fna", ''),
        "dna_seq": sample,
    }
    Decoy(**decoy_table_entry).save()
    os.chdir(i_came_from)


def run():
    # https://phagescope.deepomics.org/download#fasta
    phagi_fasta_dir  = "/storage/databases/phagescope/RefSeq"  # one sequence per data file
    # https://ftp.ncbi.nlm.nih.gov/genomes/HUMAN_MICROBIOM/
    human_microbiome = "/storage/databases/ncbi/human_microbiome"
    # https://ftp.ncbi.nlm.nih.gov/genomes/genbank/plant
    plant_genomes    = "/storage/databases/ncbi/plant_genomes"
    # for _ in range(10):
    #     create_phage_decoy(phagi_fasta_dir)
    # for _ in range(10):
    #     create_human_microbiome_decoy(human_microbiome)
