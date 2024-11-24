#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb06
# in that case django will take care of the paths and also check for migrations and such

from bad_bac_exercise.models import CARDModel, Gene, AntibioticResMutation, Drug, DrugClass

def parse_aro_obo():
    card_home = "/storage/databases/CARD-data"
    # it helps to keep in mind that 'card' here stands for the 'CARD' database name
    # The Comprehensive Antibiotic Resistance Database
    aro_file = f"{card_home}/aro.obo"
    inf = open(aro_file)

    aro_id2pubchem = {}
    reading = False
    [aro_id, pubchem] = [None, None]
    for line in inf:
        line = line.strip()
        if line[:len('[Term]')] == '[Term]':
            reading = True
            continue
        elif len(line) == 0:
            if aro_id is not  None and pubchem is not None:
                aro_id2pubchem[aro_id] = pubchem
                aro_id  = None
                pubchem = None
            reading = False
            continue
        elif reading:
            if line[:len(linestart := 'id: ARO:')] == linestart:
                aro_id = line.replace(linestart, '')
            elif line[:len(linestart := 'xref: pubchem.compound:')] == linestart:
                pubchem = line.replace(linestart, '')
    inf.close()
    return aro_id2pubchem


def run():

    aro_id2pubchem = parse_aro_obo()
    for drug in Drug.objects.all():
        if drug.aro_id not in aro_id2pubchem:
            print(f"pubchem id not found for {drug.name}")
            exit()
        drug.pubchem_id = int(aro_id2pubchem[drug.aro_id])
        drug.save()


#######################
if __name__ == "__main__":
    run()
