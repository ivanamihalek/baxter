#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

from Bio.Seq import Seq

from bad_bac_exercise.models import Pdb2Mutation
from bad_bac_exercise.models import Gene2UCSCAssembly

from bad_bac_exercise.scripts.utils import run_subprocess


def assmb_does_not_match_ref_aa(assmb_entry, gene_entry, pos, aa_from) -> bool:
    blastdbcmd = "/usr/third/blast-2.16.0/bin/blastdbcmd"
    db = "/storage/databases/ucsc/bacterial_genomes/ucsc_bac_genomes.fa"

    try:
        gene_2_assm_entry = Gene2UCSCAssembly.objects.get(gene_id=gene_entry.id, assembly_id=assmb_entry.id)
    except:
        return True
    # print(f"\t{assmb_entry.refseq_assembly_id}  {gene_2_assm_entry.pct_identity} {gene_2_assm_entry.contig} ")
    strand = gene_2_assm_entry.get_strand_on_contig_display()
    start  = gene_2_assm_entry.start_on_contig
    end    = gene_2_assm_entry.end_on_contig
    # print(f"\t  {start} {end}  {strand} ")
    cmd  = f"{blastdbcmd}  -db {db} -entry {gene_2_assm_entry.contig}  -range {start}-{end} -strand {strand}"
    ret = run_subprocess(cmd)
    seq_on_assembly = "".join(ret.split("\n")[1:])
    # print(seq_on_assembly)
    biopython_dseq = Seq(seq_on_assembly)
    biopython_translation = biopython_dseq.translate(table=11)
    # print(aa_from, pos)
    # print(biopython_translation)
    # print(biopython_translation[pos-1])
    # print(arm_entry.gene.protein_seq)
    # print(arm_entry.gene.protein_seq[pos-1])
    return biopython_translation[pos-1] != gene_entry.protein_seq[pos-1]


def run():

    # "arm" = antibiotic resistance mutation

    # start with pdb, arm pairs where the mutation is close to ligand
    # this is supposed to save us time from investigating arms that
    # we do not have crystallized, or are far from the ligand
    # for those arms, for all genes  they possilbly map to, find the assemblies
    # if the assembly has the right reference mutation, store arm2assembly map

    visited = set()
    for pdb_2_arm_entry in Pdb2Mutation.objects.filter(dist_to_drug__lte=9.0):
        arm_entry = pdb_2_arm_entry.antibio_res_mutation
        if arm_entry.id in visited: continue
        visited.add(arm_entry.id)
        print("***************************************************")
        aa_from = arm_entry.mutation[0]
        aa_tp = arm_entry.mutation[-1]
        pos = int(arm_entry.mutation[1:-1])
        print(pdb_2_arm_entry.dist_to_drug)
        print(pdb_2_arm_entry.pdb.pdb_id, arm_entry.gene.name, arm_entry.mutation)
        # a gene can map to an assembly, however, the assembly might belong
        # to a strain that does not have the expected reference aa ath the mutation position
        # now we are going to check for that
        for assmb_entry in arm_entry.gene.ucsc_assemblies.all():
            if assmb_does_not_match_ref_aa(assmb_entry, arm_entry.gene, pos, aa_from): continue
            # ok, the aa_from matches
            arm_entry.assemblies.add(assmb_entry)


#######################
if __name__ == "__main__":
    run()
