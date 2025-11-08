#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import os

import numpy as np
from Bio.PDB import PDBParser, PDBList, Select

from utils import is_nonempty_file

import os
from sys import argv

import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')
django.setup()
from models.bad_bac_models import AntibioticResMutation, Pdb2Gene, Pdb2Drug, Pdb2Mutation

def download(pdbdir, pdb_id):

    pdbfile = f"{pdbdir}/{pdb_id.lower()}.pdb"
    if is_nonempty_file(pdbfile):
        # print(f"found file: {pdbfile}")
        pass
    else:
        pdb_list = PDBList()
        try:
            pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir=pdbdir, file_format="pdb")
        except Exception as e:
            # print(e)
            return ""
        if not is_nonempty_file(pdbfile):
            return ""
        # print(f"Downloaded file: {pdb_filename}")
        os.rename(pdb_filename, pdbfile)
        # print(f"renamed to: {pdbfile}")

    return pdbfile


def ligand_distance(pdb_id, ligand_res_name, query_chain_id, protein_res_id) -> float:

    parser = PDBParser(QUIET=True)
    pdbfile = download("/storage/databases/pdb/structures", pdb_id)
    if not pdbfile:
        return -1
    if not is_nonempty_file(pdbfile):
        print(f"{pdbfile} seems to be empty")
        return -1
    structure = parser.get_structure(pdb_id, pdbfile)

    # for the ligand, the chain label may correspond to one of the peptide labels, or can be on its own
    ligand_residues = []
    for chain in structure.get_chains():
        for residue in chain:
            if residue.get_resname() != ligand_res_name: continue
            ligand_residues.append(residue)

    if not ligand_residues:
        print(f"ligand {ligand_res_name} not found in {pdb_id}")
        return -1

    # Find the target residue in the specified chain
    target_residue = None
    chain = structure[0][query_chain_id]  # this is not NMR
    for residue in chain:
        res_id = residue.get_full_id()[3][1]
        if res_id == protein_res_id:
            target_residue = residue
            break

    if target_residue is None:
        print(f"residue {protein_res_id} not found in {pdb_id}, chain {query_chain_id}")
        return -1

    # Calculate minimum distance if both are found
    min_distance = float('inf')
    for ligand_residue in ligand_residues:
        for ligand_atom in ligand_residue.get_atoms():
            for target_atom in target_residue.get_atoms():
                distance = np.linalg.norm(ligand_atom.coord - target_atom.coord)
                if distance < min_distance:
                    min_distance = distance

    return min_distance


def run():

    # we will calculate the distance from the mutated residue to the nearest drug, but only for those
    # pdbs that map to one of our genes of interest
    distinct_genes = set()
    total_candidates = 0
    for abr_mutation in AntibioticResMutation.objects.all():
        # not sure what is the purpose of the "Var" thing, but cannot deal with it now
        if abr_mutation.mutation[-3:] == "Var": continue
        mutation_pos = int(abr_mutation.mutation[1:-1])
        mutation_from = abr_mutation.mutation[0]
        mutation_to = abr_mutation.mutation[-1]
        gene_name = abr_mutation.gene.name
        print(f"{gene_name}  {mutation_pos}  {mutation_from} -> {mutation_to}")

        drugs_affected = list(abr_mutation.drugs_affected.all())
        if len(drugs_affected) == 0:
            print(f"\t no drugs affected")
            continue
        else:
            drug_names = [d.name for d in drugs_affected]
            print(f"\t drugs affected {drug_names}")

        pdb_structures = list(abr_mutation.gene.pdbstructure_set.all())
        if len(pdb_structures) == 0:
            print(f"\t no pdb structures")
            continue
        else:
            pdb_ids = [p.pdb_id for p in pdb_structures]
            print(f"\t pdb structures {pdb_ids}")

        # find drug<->pdb pairs
        mapped_pdb_drug_pairs = []
        for pdb_entry in pdb_structures:
            drugs_in_pdb = list(pdb_entry.drugs.all())
            for drug_entry in drugs_affected:
                if drug_entry in drugs_in_pdb:
                    mapped_pdb_drug_pairs.append((pdb_entry, drug_entry))
        if not mapped_pdb_drug_pairs: continue

        outstr = f"\n{abr_mutation.mutation}, {abr_mutation.gene.name}\n"
        mutation_mapped = False
        for p, d in mapped_pdb_drug_pairs:
            pdb_2_drug_entry =  Pdb2Drug.objects.filter(pdb_id=p, drug_id=d)[0]
            ligand_res_name = pdb_2_drug_entry.drug_residue_name
            if ligand_res_name is None: continue

            for mapping_entry in Pdb2Gene.objects.filter(pdb_id=p, gene=abr_mutation.gene):

                # is the residue outside this piece of structure?
                if not (mapping_entry.gene_seq_start <= mutation_pos <= mapping_entry.gene_seq_end):
                    # print(f"the model residue outside of the petide range")
                    continue

                # the numbering may not match for some legit reason, but
                # we are not going to deal with this now
                pos_on_gene_seq = mutation_pos - mapping_entry.gene_seq_start
                gene_aa = mapping_entry.gene_seq_aligned[pos_on_gene_seq]
                if mutation_from != gene_aa:
                    # print(f"aa mismatch between mutation position and the protein translation for the {gene_name} gene")
                    continue

                pos_on_pdb_seq = mutation_pos - mapping_entry.pdb_seq_start
                pdb_aa = mapping_entry.pdb_seq_aligned[pos_on_pdb_seq]
                if mutation_from != pdb_aa:
                    # print(f"aa mismatch between mutation position and the pdb sequence for the {gene_name} gene")
                    continue

                chain_id = f"{mapping_entry.pdb_chain}"
                distance = ligand_distance(p.pdb_id, ligand_res_name, chain_id,  mutation_pos)
                if distance > 11: continue  # this is too far to make difference, for the ligand et least
                if distance < 0: continue  # something went wrong

                (pdb_2_abr, was_created) = Pdb2Mutation.objects.update_or_create(pdb_id=p.id, antibio_res_mutation_id=abr_mutation.id)
                pdb_2_abr.dist_to_drug = round(distance, 1)
                pdb_2_abr.save()

                outstr += f"\t {d.name}, {chain_id}, {mapping_entry.gene_seq_start}, {mapping_entry.gene_seq_end}  "
                outstr += f"{mutation_from}  {gene_aa}  {pdb_aa}\n"
                outstr += f"{ligand_res_name}  {distance} \n"

                distinct_genes.add(abr_mutation.gene.name)
                mutation_mapped = True

        if mutation_mapped:
            total_candidates += 1
            print(outstr)
    print(distinct_genes)
    print(total_candidates)


def test_run():
    distance = ligand_distance(pdb_id="1c14", ligand_res_name="TCL", query_chain_id="A", protein_res_id=93)
    print(distance)

#######################
def main():
    if len(argv) < 2:
        exit("Tell me what to do - test or run?")
    if argv[1] == 'run':
        run()
    elif argv[1] == 'test':
        test_run()
    else:
        print(f"I don't know what is '{argv[1]}'. What should I do - test or run?")

#######################
if __name__ == "__main__":
    main()

