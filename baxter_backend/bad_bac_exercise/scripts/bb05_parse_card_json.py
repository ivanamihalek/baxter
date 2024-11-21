#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import json
from pprint import pprint
from Bio.Seq import Seq
from bad_bac_exercise.models import CARDModel, Gene, AntibioticResMutation, Drug, DrugClass


# reconstructing this, for example:
# https://card.mcmaster.ca/ontology/40192

# store the card info, including the protein and cdna seqs

# to map to available assemblies, download all uccs sdna seqs from ncbi, for ex
# https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/012/525/GCF_000012525.1_ASM1252v1/
# create searchable db (blat format)
# then find mapping to assembly by blasting

# any of ucsc assemblies usable? is there any way to get to this data except page scraping?
# https://genome.ucsc.edu/cgi-bin/hgGateway paste in the assembly
# if there aren't too many, can request here https://genome.ucsc.edu/assemblyRequest.html

# from tax id can get to assembly, for example
# https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=242231
# search for PDB entry bound to a small molecule, using the protein sequence
# check if the small molecule belongs to one of the classes that are supposed to be disrupted
# https://en.wikipedia.org/wiki/List_of_antibiotics
# download structures and check proximity of the residues to the small molecule

def check_translation_sanity(pseq, dseq, card_name):
    biopython_dseq = Seq(dseq)
    # the table number refers to NCBI tables
    # for bacterial, see here: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11
    biopython_translation =  biopython_dseq.translate(table=11)
    translation = str(biopython_translation)[:-1]  # the last one is the stop codon (*)
    # For some reason Biopython translates start codon as M,
    # even though the actual translation might be V or L, as it sometimes happens in e.g. E Coli.
    # Also, Mycoplasma read UGA as a tryptophan rather than a stop codon

    if translation[1:].replace("*", "W") != pseq[1:]:
        print(f"Sequence sanity check failed for {card_name}")
        print("protein sequence")
        print(pseq)
        print("dna translation")
        print(translation)
        print("differenc points")
        for i in range(len(pseq)):
            if pseq[i] == translation[i]: continue
            print("\t", i, pseq[i], translation[i])
        print()
        exit()


def store_gene_info(card_entry: CARDModel, card_dict: dict) -> Gene:
    card_name = card_entry.card_name
    number_of_seqs = len(card_dict["model_sequences"]["sequence"])
    if number_of_seqs > 1:
        print(f"A problem: {card_name} has {number_of_seqs} associated with the protein.")
        print(f"(At the time of this writing this problem did not exist.)")
        exit()

    sequence_dict = list(card_dict["model_sequences"]["sequence"].values()).pop()
    if sequence_dict['dna_sequence']['partial'] != '0':
        print(card_name, "the dna seq is partial")
        print(f"(At the time of this writing this problem did not exist.)")
        exit()

    pseq = sequence_dict['protein_sequence']['sequence']
    dseq = sequence_dict['dna_sequence']['sequence']
    check_translation_sanity(pseq, dseq, card_name)
    gene_name  = card_name.split("_")[1]
    gene_entry = Gene.objects.get(name=gene_name)
    gene_entry.protein_seq = pseq
    gene_entry.dna_seq = dseq
    gene_entry.save()

    return gene_entry


def store_drug_info(card_dict: dict):

    drug_entries = []
    drug_class_entries = []
    for identifier, sub_entry in card_dict["ARO_category"].items():
        if "category_aro_class_name" not in sub_entry: continue
        if sub_entry["category_aro_class_name"] not in ["Drug Class", "Antibiotic"]: continue
        if sub_entry["category_aro_class_name"] == "Drug Class":
            drug_class_name = sub_entry["category_aro_name"]
            # some ad hoc input cleanup
            if drug_class_name[-len(" antibiotic"):] == " antibiotic":
                drug_class_name = drug_class_name[:-len(" antibiotic")]

            print(f"#### {drug_class_name}")
            (drug_class_entry, was_created) = DrugClass.objects.update_or_create(name=drug_class_name)
            drug_class_entries.append(drug_class_entry)
        else:
            drug_name = sub_entry["category_aro_name"]
            print(f"***** {drug_name}")
            (drug_entry, was_created) = Drug.objects.update_or_create(name=drug_name)
            drug_entries.append(drug_entry)

    return drug_entries, drug_class_entries


def run():
    card_home = "/storage/databases/CARD-data"
    # it helps to keep in mind that 'card' here stands for the 'CARD' database name
    # The Comprehensive Antibiotic Resistance Database
    jsonfile = f"{card_home}/card.json"
    card_dict = json.load(open(jsonfile))
    for card_id, card_dict in card_dict.items():
        if not isinstance(card_dict, dict): continue  # card_dict['_version'] == '3.3.0'
        if 'ARO_accession' not in card_dict: continue  # card_dict['description']
        card_name  = card_dict['CARD_short_name']
        try:
            card_entry = CARDModel.objects.get(card_name=card_name)
        except CARDModel.DoesNotExist:
            continue
        print()
        print(card_name)
        gene_entry = store_gene_info(card_entry, card_dict)
        (drug_entries, drug_class_entries) = store_drug_info(card_dict)
        # mutation to drug and drug class mapping
        for abr in AntibioticResMutation.objects.filter(gene=gene_entry):
            print(abr.mutation)
            for drug_entry in drug_entries:
                print(drug_entry.name)
                abr.drugs_affected.add(drug_entry)
            for drug_class_entry in drug_class_entries:
                print(drug_class_entry.name)
                abr.drug_classes_affected.add(drug_class_entry)


#######################
if __name__ == "__main__":
    run()
