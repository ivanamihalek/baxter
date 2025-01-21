#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_snp_txt
# in that case django will take care of the paths and also check for migrations and such

import json
from pprint import pprint
from models.bad_bac_models import CARDModel, Gene, AntibioticResMutation


def process_card_snps(card_home: str) -> tuple[dict, dict]:
    snp_file = f"{card_home}/snps.txt"
    inf = open(snp_file)
    ct = 0

    mutations = {}
    card_short2card_description = {}
    for line in inf:
        fields = line.strip().split("\t")
        if fields[0] == 'Accession': continue  # this is header
        if fields[3] != 'single resistance variant': continue
        if fields[2] != 'protein variant model': continue
        if len(fields[1]) < 5: continue  # not sure what these are
        card_short_name = fields[-1]
        if len(card_short_name.split("_")) < 3:  continue  # not sure what these are either

        if card_short_name not in mutations: mutations[card_short_name] = []
        mutations[card_short_name].append(fields[-2])

        if card_short_name in card_short2card_description:
            if card_short2card_description[card_short_name] != fields[1]:
                print(f"unexpected multiple descriptions for a short CARD name:")
                print(f"{card_short_name}:  {card_short2card_description[card_short_name]}  {fields[1]}")
                exit()
        else:
            card_short2card_description[card_short_name] = fields[1]
        ct += 1

    return card_short2card_description, mutations


def run():
    card_home = "/storage/databases/CARD-data"
    (card_short2card_description, mutations) = process_card_snps(card_home)
    for card_short, descr in card_short2card_description.items():
        gene_name = card_short.split("_")[1]
        print(len(card_short), card_short, descr, gene_name, mutations[card_short])
        (card_entry, was_created) = CARDModel.objects.update_or_create(card_name=card_short)
        if was_created:
            card_entry.card_description = descr
            card_entry.save()

        (gene_entry, was_created) = Gene.objects.update_or_create(name=gene_name)
        for mutation in mutations[card_short]:
            fields = {'mutation': mutation, 'gene': gene_entry}
            (mutation_entry, was_created) = AntibioticResMutation.objects.update_or_create(**fields)
            mutation_entry.card_models.add(card_entry)


#######################
if __name__ == "__main__":
    run()
