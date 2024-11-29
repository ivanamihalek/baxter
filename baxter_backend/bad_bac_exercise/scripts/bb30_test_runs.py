#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

from bad_bac_exercise.models import CARDModel, Gene, AntibioticResMutation, Drug, DrugClass


def run():
    for gene_entry in Gene.objects.all():
        print(gene_entry.name)
        for abr in AntibioticResMutation.objects.filter(gene=gene_entry):
            print(f"\t {abr.mutation}")
            for drug_entry in abr.drugs_affected.all():
                print(f"\t\t {drug_entry.name}")


#######################
if __name__ == "__main__":
    run()
