#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import os

import subprocess
from glob import glob

from Bio.Blast import NCBIXML
from bad_bac_exercise.models import AntibioticResMutation, Pdb2Gene, Pdb2Drug, Pdb2Mutation, PDBStructure
from bad_bac_exercise.models import UCSCAssembly, Gene2UCSCAssembly


from .utils import is_nonempty_file


def map_contigs_to_assembly(blastdb_home):
    contigs2assembly = {}
    for txttfile in glob(f"{blastdb_home}/*.contents.txt"):
        refseq_assembly_id = txttfile.split("/")[-1].replace(".contents.txt", "")
        with open(txttfile) as inf:
            for contig in inf.read().replace(">", "").split("\n"):
                if not contig: continue
                contigs2assembly[contig] = refseq_assembly_id
    return contigs2assembly


def run_blast(blastdb_home, query_file, output_file):
    # Define parameters
    db = f"{blastdb_home}/uscs_bac_genomes.fa"

    print(f"running blastp")
    # Run the blastp command
    subprocess.run(["blastn", "-db", db, "-query", query_file, "-out", output_file, "-outfmt", "5"])


def parse_blast_results(gene_entry, blast_results_file, contigs2assembly):
    blast_res_handle = open(blast_results_file)
    blast_records = NCBIXML.parse(blast_res_handle)
    for blast_record in blast_records:
        print("\n************************")
        print(f"Query: {blast_record.query} length {len(gene_entry.dna_seq)}")

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.strand[0] == "Minus":
                    print(f"The query was interpreted as being on the neg strand - how could that happen?")
                    exit()
                # print(alignment.title, hsp.query, hsp.match, hsp.sbjct)
                if hsp.expect > 0: continue

                pct_identity = round((hsp.identities/hsp.align_length*100), 2)
                if pct_identity < 98: continue

                assembly = contigs2assembly.get(alignment.hit_def, "unk")
                if assembly == "unk": continue
                qrylen = len(gene_entry.dna_seq)
                qry_aligned   = hsp.query_end - hsp.query_start + 1
                plus_strand =  hsp.strand[1] == "Plus"
                if plus_strand:
                    sbjct_aligned = hsp.sbjct_end - hsp.sbjct_start + 1
                else:
                    sbjct_aligned = hsp.sbjct_start - hsp.sbjct_end + 1
                if qry_aligned != qrylen or sbjct_aligned != qrylen: continue

                # strand          Tuple of (query, target) strand e.g. ('Plus', 'Minus')
                # not sure how or why would the fist one be minus
                print(alignment.hit_def, assembly, pct_identity, hsp.query_start, hsp.query_end,
                      hsp.sbjct_start, hsp.sbjct_end, hsp.strand)
                if hsp.strand[0] == "Minus":
                    print(f"The query was interpreted as being on the neg strand - how could that happen?")
                    exit()
                assembly_entry = UCSCAssembly.objects.filter(refseq_assembly_id=assembly)[0]
                junction_table_entry = {
                    "gene": gene_entry,
                    "assembly": assembly_entry,
                    "pct_identity": pct_identity,
                    "contig": alignment.hit_def,
                    "start_on_contig": hsp.sbjct_start if plus_strand else hsp.sbjct_end,
                    "end_on_contig": hsp.sbjct_end if plus_strand else hsp.sbjct_start,
                    "strand_on_contig": Gene2UCSCAssembly.Strand.translate(hsp.strand[1])
                }
                (gene2assm_entry, was_created) = Gene2UCSCAssembly.objects.update_or_create(**junction_table_entry)
                gene2assm_entry.save()


def run():
    blastdb_home = "/storage/databases/ucsc/bacterial_genomes"
    scratch = "/home/ivana/scratch/baxter/blast_against_genomes"

    genes = set()
    for pdb2mut in Pdb2Mutation.objects.filter(dist_to_drug__lt=10):
        distance = pdb2mut.dist_to_drug
        pdb_entry = PDBStructure.objects.get(pk=pdb2mut.pdb_id)
        abr_entry = AntibioticResMutation.objects.get(pk=pdb2mut.antibio_res_mutation_id)
        # print(f"{pdb_entry.pdb_id}   {abr_entry.gene.name}   {abr_entry.mutation}   {distance:.1f}")
        genes.add(abr_entry.gene)

    contigs2assembly = map_contigs_to_assembly(blastdb_home)

    for gene_entry in genes:
        blast_results_file = f"{scratch}/{gene_entry.name}_blastout.xml"
        if is_nonempty_file(blast_results_file):
            pass
        else:
            query_file  = f"{scratch}/{gene_entry.name}.fasta"
            with open(query_file, "w") as outf:
                print(f"> {gene_entry.name}", file=outf)
                print(gene_entry.dna_seq, file=outf)
            run_blast(blastdb_home, query_file, blast_results_file)

        parse_blast_results(gene_entry, blast_results_file, contigs2assembly)

#######################
if __name__ == "__main__":
    run()
