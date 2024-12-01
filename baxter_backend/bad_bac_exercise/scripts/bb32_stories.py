#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import os

from bad_bac_exercise.models import AntibioticResMutation, Pdb2Gene, Pdb2Drug, Pdb2Mutation, PDBStructure, Decoy
from bad_bac_exercise.models import UCSCAssembly, Gene2UCSCAssembly
import random

def story_to_pptx(arm_entry, pdb_entries, ecount):


    pdb_ids = [pe.pdb_id for pe in pdb_entries]

    aa_from = arm_entry.mutation[0]
    aa_to = arm_entry.mutation[-1]
    pos = int(arm_entry.mutation[1:-1])
    # dist = pdb_2_arm_entry.dist_to_drug
    print(f"\n******************** EXERCISE {ecount}  *******************************")
    assmb_entry = arm_entry.assemblies.first()
    print(arm_entry.gene.name, arm_entry.mutation, assmb_entry.refseq_assembly_id, pdb_ids)

    ########################################
    # FINDING THE SPECIES
    # select region from  assmb_entry.refseq_assembly_id, + 2 decoys  --> task: find problematic species using NCBI Blast
    fingerprint = random.choice(assmb_entry.fingerprint_set.all()).dna_seq
    decoy1 = random.choice(Decoy.objects.all()).dna_seq
    decoy2 = random.choice(Decoy.objects.all()).dna_seq
    q = f"""Q1:  In the field, using your Flongle, you detected the following three DNA snippets:
    {fingerprint[:10]}   {decoy1[:10]}   {decoy2[:10]}   
    Which one is problematic?  Which species does it belong to? """
    a = f"{assmb_entry.common_name}  "
    print(q, "\n", a)

    ########################################
    # DETECTING THE VARIANT
    # toy output from sequencer - the name xxx  --> task: find variants using Galaxy
    q = f"""Q2: You took the sample to the lab for the full sequencing. The lab returned the following sequencing file. 
    The lab also determined that the strain that this sample contains is {assmb_entry.refseq_assembly_id}
     Which variants does this strain cary? What are their genomic coords?"""
    # TODO genomic coords
    # TODO toy sequencer output
    a = f" genomic coords"
    print(q, "\n", a)

    ########################################
    # PLACING THE VARIANT IN GENOMIC CONTEXT
    # ---> task: find genomic location from the previous step using UCSC Genome Browser and some calculation - gene and protein location
    q = """Q3: Which genomic locations do these variants hit? Which gene, if any? 
    Which codon is hit on the protein level, if any?"""
    a = f"{arm_entry.gene.name}, {aa_from} at the position {pos}"
    print(q, "\n", a)

    ########################################
    # DETERMINING THE VARIANT CONSEQUENCE ON GENOMIC LEVEL
    # ---> task what is the impact of mutation - use codon translation table Wikipedia
    q = """Q4: What is the impact of the nucleotide change? """
    a = f"{aa_from}>{aa_to}"
    print(q, "\n", a)

    ########################################
    # SEARCHING FOR THE PREVIOUS KNOWLEDGE
    # ---> is this mutation known?  CARD database -- should get the drug name
    q = """Q5: Is this mutation already known in the literature? Which drug resistance it may cause?"""
    # TODO fetch pubmed IDs
    drugs = [drug.name for drug in arm_entry.drugs_affected.all()]
    a = f"{drugs} pubmed ids -- any should do"
    print(q, "\n", a)

    ########################################
    # UNDERSTANDING THE IMPACT ON THE SYSTEMC LEVEL
    # ---> PDB search for species + gene + drug name -- make illustration
    q = """Q6: Is the structure of the protein-drug complex known for this bacterial species?
    Use it to inspect the impact that the mutation will have on the drug binding."""
    a = f"any of the pdbs that I have here, {pdb_ids}"
    print(q, "\n", a)


def run():

    # "arm" = antibiotic resistance mutation

    # start with pdb, arm pairs where the mutation is close to ligand
    # this is supposed to save us time from investigating arms that

    ecount = 0
    # for pdb_2_arm_entry in Pdb2Mutation.objects.filter(dist_to_drug__lte=9.0):
    for arm_entry in AntibioticResMutation.objects.filter(assemblies__isnull=False).distinct():
        pdb_entries = PDBStructure.objects.filter(
            pdb2mutation__antibio_res_mutation=arm_entry,
            pdb2mutation__dist_to_drug__lt=9
        ).distinct()
        if len(pdb_entries) < 1: continue
        ecount += 1
        story_to_pptx(arm_entry, pdb_entries, ecount)

        # store to db?
        # question_and_answers
        # category, title, question, answers, specific

#######################
if __name__ == "__main__":
    run()
