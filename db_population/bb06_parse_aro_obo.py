#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb06
# in that case django will take care of the paths and also check for migrations and such
import json

import os
import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')
django.setup()

from models.bad_bac_models import CARDModel, Drug, Publication

def parse_aro_obo():
    card_home = "/storage/databases/CARD"
    # it helps to keep in mind that 'card' here stands for the 'CARD' database name
    # The Comprehensive Antibiotic Resistance Database
    aro_file = f"{card_home}/aro.obo"
    inf = open(aro_file)

    aro_id2pubchem = {}
    card_name2pubmed_ids = {}
    reading = False
    [aro_id, pubchem, card_name_short, pubmed_ids] = [None, None, None, None]
    for line in inf:
        line = line.strip()
        if line[:len('[Term]')] == '[Term]':
            reading = True
            continue
        elif len(line) == 0:
            if aro_id is not None and pubchem is not None:
                aro_id2pubchem[aro_id] = pubchem
            if card_name_short is not None and pubmed_ids is not None:
                card_name2pubmed_ids[card_name_short] = pubmed_ids
            [aro_id, pubchem, card_name_short, pubmed_ids] = [None, None, None, None]
            reading = False
            continue
        elif reading:
            if line[:len(linestart := 'id: ARO:')] == linestart:
                aro_id = line.replace(linestart, '')
            elif line[:len(linestart := 'xref: pubchem.compound:')] == linestart:
                pubchem = line.replace(linestart, '')
            elif line[:len(linestart := 'xref: PubChem:')] == linestart:
                pubchem = line.replace(linestart, '')
            elif "EXACT CARD_Short_Name" in line:
                card_name_short = line.split()[1].replace('"', '')
            elif "PMID:" in (possible_pmid_list := line.split('[')[-1]):
                pmid_list = possible_pmid_list.replace(']', '').split(",")
                pubmed_ids = []
                for pmid in pmid_list:
                    try:
                        pubmed_id = int(pmid.replace("PMID:", ''))
                    except:
                        print(f"found {pmid} among PUBMED ids")
                        continue
                    pubmed_ids.append(pubmed_id)
    inf.close()
    return aro_id2pubchem, card_name2pubmed_ids


def run():
    print(f"Updates drugs with pubchem ids")
    print(f"Creates through table card_model_2_publication")
    (aro_id2pubchem, card_name2pubmed_ids) = parse_aro_obo()
    for drug in Drug.objects.all():
        if drug.aro_id not in aro_id2pubchem:
            print(f"pubchem id not found for {drug.name}, aro id: {drug.aro_id}")
            exit()
        drug.pubchem_id = int(aro_id2pubchem[drug.aro_id])
        drug.save()
    
    for card_model in CARDModel.objects.all():
        if card_model.card_name not in card_name2pubmed_ids:
            print(f"pubmed ids not found for {card_model.card_name}")
            continue
        for pmid in card_name2pubmed_ids[card_model.card_name]:
            (publication_entry, was_created) = Publication.objects.update_or_create(pubmed_id=pmid)
            card_model.publications.add(publication_entry)
        card_model.save()


#######################
if __name__ == "__main__":
    run()
