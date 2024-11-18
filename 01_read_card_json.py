#! /usr/bin/env python

import json
from pprint import pprint


# reconstructing this, for example:
# https://card.mcmaster.ca/ontology/40192


def parse_card_entry(card_entry: dict):
    pprint(card_entry)
    exit()


def main():
    jsonfile = "/storage/databases/CARD-data/card.json"
    card_dict = json.load(open(jsonfile))

    for card_id, card_entry in card_dict.items():
        if not isinstance(card_entry, dict): continue  # card_dict['_version'] == '3.3.0'
        if 'ARO_accession' not in card_entry: continue # card_dict['description']

        if card_entry['ARO_accession'] == '3003582':
            print(card_id)
            parse_card_entry(card_entry)


#######################
if __name__ == "__main__":
    main()
