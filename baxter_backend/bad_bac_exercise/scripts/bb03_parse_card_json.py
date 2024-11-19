#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import json
from pprint import pprint


# reconstructing this, for example:
# https://card.mcmaster.ca/ontology/40192

# start building the database
# is the reference sequence available in genome browser
# https://hgdownload.soe.ucsc.edu/hubs/bacteria/index.html
# any of those usable? is there any way to get to this data except page scraping?
# https://genome.ucsc.edu/cgi-bin/hgGateway paste in the assembly
# if there aren't too many, can request here https://genome.ucsc.edu/assemblyRequest.html

# from tax id can get to assembly, for example
# https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=242231
# search for PDB entry bound to a small molecule, using the protein sequence
# check if the small molecule belongs to one of the classes that are supposed to be disrupted
# https://en.wikipedia.org/wiki/List_of_antibiotics
# download structures and check proximity of the residues to the small molecule

def parse_card_entry(card_entry: dict):
    # pprint(card_entry)
    print(card_entry["ARO_name"])
    print(card_entry["ARO_description"])
    for identifier, sub_entry in card_entry["model_sequences"]["sequence"].items():
        print()
        print("sequence entry:", identifier)
        pprint(sub_entry)
    print()
    for identifier, sub_entry in card_entry["ARO_category"].items():
        if "category_aro_class_name" not in sub_entry: continue
        if sub_entry["category_aro_class_name"] not in ["Drug Class", "Antibiotic"]: continue
        print()
        print("drug entry:", identifier)
        print("class name:", sub_entry["category_aro_class_name"])
        print("drug name:", sub_entry["category_aro_name"])  # this may belong to an additional ddrug class

    print("*************************\n")


def process_card_json(card_home: str, accession_numbers: dict):
    jsonfile = f"{card_home}/card.json"
    card_dict = json.load(open(jsonfile))

    acc_number2card_short = {v: k for k, v in accession_numbers.items()}

    for card_id, card_entry in card_dict.items():
        if not isinstance(card_entry, dict): continue  # card_dict['_version'] == '3.3.0'
        if 'ARO_accession' not in card_entry: continue  # card_dict['description']

        if card_entry['ARO_accession'] not in acc_number2card_short.keys(): continue
        print(card_id)
        parse_card_entry(card_entry)


def process_card_snps(card_home: str) -> tuple[dict, dict]:
    snp_file = f"{card_home}/snps.txt"
    inf = open(snp_file)
    ct = 0

    mutations = {}
    accession_number = {}
    for line in inf:
        fields = line.strip().split("\t")
        if fields[0] == 'Accession': continue  # this is header
        if fields[3] != 'single resistance variant': continue  # this is header
        if fields[2] != 'protein variant model': continue  # this is header
        if len(fields[1]) < 5: continue  # not sure what these are
        card_short_name = fields[-1]

        if card_short_name not in mutations: mutations[card_short_name] = []
        mutations[card_short_name].append(fields[-2])

        if card_short_name in accession_number:
            if accession_number[card_short_name] != fields[0]:
                print(f"unexpected multiple accession numbers for a short CARD name:")
                print(f"{card_short_name}:  {accession_number[card_short_name]}  {fields[0]}")
        else:
            accession_number[card_short_name] = fields[0]
        ct += 1

    # print(f"total protein mutations: {ct}")
    # print(f"total card entries: {len(mutations)}")
    # for card_short_name, mutation_list in mutations.items():
    #     print(f"{accession_numbers[card_short_name]} {card_short_name}  {mutation_list}")
    return accession_number, mutations


def run():
    card_home = "/storage/databases/CARD-data"
    (accession_numbers, mutations) = process_card_snps(card_home)

    process_card_json(card_home, accession_numbers)


#######################
if __name__ == "__main__":
    run()
