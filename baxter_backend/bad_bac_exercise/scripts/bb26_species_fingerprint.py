#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import os
from glob import glob
import random
from pprint import pprint

from Bio import Entrez

from bad_bac_exercise.models import AntibioticResMutation, PDBStructure, Fingerprint


from bad_bac_exercise.scripts.utils import run_subprocess


def get_assembly_contigs(blastdb_home) -> dict:
    assembly_contigs = {}
    for txttfile in glob(f"{blastdb_home}/*.contents.txt"):
        refseq_assembly_id = txttfile.split("/")[-1].replace(".contents.txt", "")
        assembly_contigs[refseq_assembly_id] = []
        with open(txttfile) as inf:
            for contig in inf.read().replace(">", "").split("\n"):
                if not contig: continue
                assembly_contigs[refseq_assembly_id].append(contig)
    return assembly_contigs


def assemblies_of_interest() -> set:
    # at some point later this may be all species in the db
    # for pdb_2_arm_entry in Pdb2Mutation.objects.filter(dist_to_drug__lte=9.0):
    assembly_entries = set()
    for arm_entry in AntibioticResMutation.objects.filter(assemblies__isnull=False).distinct():
        pdb_entries = PDBStructure.objects.filter(
            pdb2mutation__antibio_res_mutation=arm_entry,
            pdb2mutation__dist_to_drug__lt=9
        ).distinct()
        if len(pdb_entries) < 1: continue
        assembly_entries.add(arm_entry.assemblies.first())
    return assembly_entries


def contig_length(blastdbcmd, db, contig) -> int:
    # use blastdbcmd to fin the contig length:
    cmd = f'{blastdbcmd} -db {db} -entry {contig} -outfmt "%l"'
    return int(run_subprocess(cmd))


def run():
    blastdbcmd = "/usr/third/blast-2.16.0/bin/blastdbcmd"
    blastdb_home = "/storage/databases/ucsc/bacterial_genomes"
    db = f"{blastdb_home}/ucsc_bac_genomes.fa"
    assembly_contigs = get_assembly_contigs(blastdb_home)

    for assembly_entry in assemblies_of_interest():
        refseq_assembly_id = assembly_entry.refseq_assembly_id
        # find the longest contig, to stay on the safe side
        contig_ids = assembly_contigs.get(refseq_assembly_id, None)
        if not contig_ids:
            print(f"no contigs found for {assembly_entry.refseq_assembly_id}")
            exit()
        length_of_contig = {contig_id: contig_length(blastdbcmd, db, contig_id) for contig_id in contig_ids}
        the_longest_contig_id = sorted(contig_ids, key=lambda c: length_of_contig[c])[-1]

        print()
        print(refseq_assembly_id)
        pprint(length_of_contig)
        print(the_longest_contig_id)

        # let's get a couple of fingerprints, just for the giggles
        for _ in range(5):
            sample_length = random.choice(range(800, 1800))
            sample_start  = random.choice(range(1, length_of_contig[the_longest_contig_id]-1800)) - 1
            if sample_start < 0:
                print(f"theo longest contig in {assembly_entry.refseq_assembly_id} too short (?)")
                exit()
            sample_end = sample_start + sample_length - 1
            cmd  = f"{blastdbcmd} -db {db} -entry {the_longest_contig_id} -range {sample_start}-{sample_end}"
            ret = run_subprocess(cmd)
            sample = "".join(ret.split("\n")[1:])
            if "N" in sample:
                continue
            fingerprint_table_entry = {
                "identifier": the_longest_contig_id,
                "dna_seq": sample,
                "assembly": assembly_entry
            }
            Fingerprint(**fingerprint_table_entry).save()
