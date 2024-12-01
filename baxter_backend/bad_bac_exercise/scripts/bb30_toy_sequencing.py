#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import gffutils
import os.path
import portion as P
import random
import subprocess

import wget
from Bio.Data import CodonTable

from Bio.Seq import Seq

from bad_bac_exercise.scripts.utils import run_subprocess

from bad_bac_exercise.models import AntibioticResMutation, PDBStructure, Decoy, Gene2UCSCAssembly


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


def download_annotation_files(storage_dir, verbose=False) -> list:
    # for each assembly of interest
    # check that I have the gff annotation file
    # example
    # https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2_genomic.gff.gz
    base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all"

    i_came_from = os.getcwd()
    os.chdir(storage_dir)

    annotation_files = []
    for assembly_entry in assemblies_of_interest():
        refseq = assembly_entry.refseq_assembly_id
        ncbi_acc = assembly_entry.ncbi_accession_number
        fnm_unzipped = f"{refseq}_{ncbi_acc}_genomic.gff"
        if os.path.exists(fnm_unzipped):
            if verbose: print(f"{fnm_unzipped} found in {storage_dir}")
        else:
            fnm = fnm_unzipped + ".gz"
            url = f"{base_url}/{refseq[:3]}/{refseq[4:7]}/{refseq[7:10]}/{refseq[10:13]}/{refseq}_{ncbi_acc}/{fnm}"
            if verbose: print(url)
            wget.download(url, bar=wget.bar_thermometer)
            subprocess.run(["gunzip", fnm])
        annotation_files.append(fnm_unzipped)
    os.chdir(i_came_from)
    return annotation_files


def genes_of_interest() -> set:
    gene_entries = set()
    for arm_entry in AntibioticResMutation.objects.filter(assemblies__isnull=False).distinct():
        pdb_entries = PDBStructure.objects.filter(
            pdb2mutation__antibio_res_mutation=arm_entry,
            pdb2mutation__dist_to_drug__lt=9
        ).distinct()
        if len(pdb_entries) < 1: continue
        gene_entries.add(arm_entry.gene)

    return gene_entries


def ann_files_to_db(storage_dir, scratch_dir, annotation_files, verbose=False):
    for annf in annotation_files:
        infnm = f"{storage_dir}/{annf}"
        dbfnm = f"{scratch_dir}/{annf}.db"
        if os.path.exists(dbfnm):
            if verbose: print(f"found {dbfnm} in {scratch_dir}")
        else:
            gffutils.create_db(infnm, dbfn=dbfnm, force=True, keep_order=True,
                               merge_strategy='merge', sort_attribute_values=True)


def annotation_sanity_check(scratch_dir, verbose=False):

    # for each gene of interest, check that the start/end coords match gff
    # https://daler.github.io/gffutils/
    for gene_entry in genes_of_interest():
        assembly = gene_entry.ucsc_assemblies.first()
        refseq = assembly.refseq_assembly_id
        ncbi_acc = assembly.ncbi_accession_number
        fnm = f"{refseq}_{ncbi_acc}_genomic.gff"
        gene_2_assmb_entry : Gene2UCSCAssembly = Gene2UCSCAssembly.objects.get(gene_id=gene_entry.id,assembly_id=assembly.id)
        if verbose:
            print()
            print(gene_entry.name, assembly.refseq_assembly_id, fnm)
            print(gene_2_assmb_entry.contig, gene_2_assmb_entry.start_on_contig, gene_2_assmb_entry.end_on_contig)
        dbfnm = f"{scratch_dir}/{fnm}.db"
        db = gffutils.FeatureDB(dbfnm, keep_order=True)
        s = gene_2_assmb_entry.start_on_contig
        e = gene_2_assmb_entry.end_on_contig

        coding_seqs = [feature for feature in db.region(start=s-1, end=e+1, completely_within=True) if feature.attributes.get('gbkey', ["none"])[0]=="CDS"]
        if len(coding_seqs) == 0:
            print(f"Error: no CDS found in the specified interval")
            exit(1)
        if len(coding_seqs) > 1:
            print(f"Error: multiple CDS found in the specified interval")
            for feature in coding_seqs:
                print()
                print(feature.start, feature.end)
                for k, v in feature.attributes.items():
                    print(k, v)
            exit(1)
        cds = coding_seqs[0]
        if verbose: print(cds.start, cds.end)
        if cds.start != s or cds.end != e:
            print(f"Error: annotation mismatch. Gene coords in gff: {cds.start}, {cds.end}")
            exit(1)


def arms_of_interest() -> set:
    arm_entries = set()
    for arm_entry in AntibioticResMutation.objects.filter(assemblies__isnull=False).distinct():
        pdb_entries = PDBStructure.objects.filter(
            pdb2mutation__antibio_res_mutation=arm_entry,
            pdb2mutation__dist_to_drug__lt=9
        ).distinct()
        if len(pdb_entries) < 1: continue
        arm_entries.add(arm_entry)
    return arm_entries


def mutation_spec(from_codon, aa_to):
    # according to transl table 11, how do I mutate my codon do get aa_to?
    codon_table = CodonTable.unambiguous_dna_by_id[11]
    to_codons = [codon for codon, amino_acid in codon_table.forward_table.items() if amino_acid == aa_to]

    mutation_codon_pos  = None
    mutation_nucleotide = None
    for to_codon in to_codons:
        diff_indices = [i for i in range(3) if from_codon[i] != to_codon[i]]
        # for now. we'll go for single nucl substitution
        if len(diff_indices) == 0:
            print(f"bleep")
            exit(1)
        elif len(diff_indices) > 1:
            continue
        else:
            i = diff_indices[0]
            mutation_codon_pos = i
            mutation_nucleotide = to_codon[i]

    if mutation_codon_pos is None or mutation_nucleotide is None:
        print(f"bleep 2")
        exit(1)

    return mutation_codon_pos, mutation_nucleotide


def create_sample_genome(scratch_dir, arm_entry):
    blastdbcmd = "/usr/third/blast-2.16.0/bin/blastdbcmd"
    db = "/storage/databases/ucsc/bacterial_genomes/ucsc_bac_genomes.fa"

    # extend the region around the gene to create "sample genome" to use in inSilicoSeq
    # make sure that the region contains some unnannotatd intervals place some decoy mutations there - single nucl, or small deletion
    # in the same gene, away from the mutation of interest put a silent mutation there

    print(arm_entry.mutation)
    aa_from = arm_entry.mutation[0]
    aa_to = arm_entry.mutation[-1]
    pos = int(arm_entry.mutation[1:-1])

    gene_entry = arm_entry.gene
    assmb_entry = gene_entry.ucsc_assemblies.first()
    refseq = assmb_entry.refseq_assembly_id
    ncbi_acc = assmb_entry.ncbi_accession_number
    fnm = f"{refseq}_{ncbi_acc}_genomic.gff"
    gene_2_assmb_entry: Gene2UCSCAssembly = Gene2UCSCAssembly.objects.get(gene_id=gene_entry.id, assembly_id=assmb_entry.id)

    # create mutation that I need
    # what are the gnomic coords of my codon
    # is the coding strand plus or minus
    strand = gene_2_assmb_entry.get_strand_on_contig_display()
    if strand=="plus": return

    actual_gene_start  = gene_2_assmb_entry.start_on_contig
    actual_gene_end    = gene_2_assmb_entry.end_on_contig
    gene_length = actual_gene_end - actual_gene_start + 1
    extension = int(2.e4)

    cmd   = f"{blastdbcmd}  -db {db} -entry {gene_2_assmb_entry.contig} "
    cmd  += f"-range {actual_gene_start-extension}-{actual_gene_end+extension} "
    ret = run_subprocess(cmd)
    toy_genome = "".join(ret.split("\n")[1:])

    # print(seq_on_assembly)
    # according to transl table 11, is the aa for the codon correct?
    toy_gene_start = extension + 1
    toy_gene_end   = toy_gene_start + gene_length - 1

    print("**********************************")
    biopython_dseq = Seq(toy_genome[toy_gene_start-1:toy_gene_end])
    if strand == "minus":
        biopython_dseq = biopython_dseq.reverse_complement()
    biopython_translation = biopython_dseq.translate(table=11)
    print(aa_from, biopython_translation[pos-1], aa_to, strand)
    codon_start = (pos-1)*3
    from_codon = biopython_dseq[codon_start:codon_start+3]
    print(from_codon, from_codon.translate(table=11))
    print(biopython_translation)

    mutation_codon_pos, mutation_nucleotide = mutation_spec(from_codon, aa_to)

    mutation_cdna_pos = codon_start + mutation_codon_pos
    # this is offset 1 coordinate, liked by blast
    mutation_genomic_position = toy_gene_start + mutation_cdna_pos

    # TODO - I am here - recallute genomic position for the neg strand
    mutated_toy_genome = toy_genome[:mutation_genomic_position-1] + mutation_nucleotide + toy_genome[mutation_genomic_position:]
    print("------------------------------------")
    biopython_dseq = Seq(mutated_toy_genome[toy_gene_start-1:toy_gene_end])
    if strand == "minus":
        biopython_dseq = biopython_dseq.reverse_complement()
    biopython_translation = biopython_dseq.translate(table=11)
    print(aa_from, biopython_translation[pos-1], aa_to, strand)
    codon_start = (pos-1)*3
    from_codon = biopython_dseq[codon_start:codon_start+3]
    print(from_codon, from_codon.translate(table=11))
    print(biopython_translation)
    print()

    dbfnm = f"{scratch_dir}/{fnm}.db"
    db = gffutils.FeatureDB(dbfnm, keep_order=True)
    s = gene_2_assmb_entry.start_on_contig - 20000
    e = gene_2_assmb_entry.end_on_contig + 20000

    # one decoy: SNV in an un-annotated region
    # complement = P.closed(s, e)
    # for feature in db.region(start=s, end=e):
    #     # feature.attributes is a dict
    #     if 'genome' in feature.attributes: continue  # this is the top-top-level interval
    #     if 'Parent' in feature.attributes: continue  # we are looking for parent-less intervals
    #     complement = complement - P.closed(feature.start, feature.end)
    # for interval in complement:  # put a SNV here as a decoy
    #     print(interval)
    #
    # another decoy: silent mutation in one of the CDS
    # coding_seqs = [feature for feature in db.region(start=s-1, end=e+1, completely_within=True) if feature.attributes.get('gbkey', ["none"])[0]=="CDS"]
    # # pick one at random
    # random_coding_seq = coding_seqs[0]


def run():

    storage_dir = "/storage/databases/ncbi/genome_annotation/bacterial"
    scratch_dir = "/home/ivana/scratch/baxter"

    annotation_files = download_annotation_files(storage_dir)
    ann_files_to_db(storage_dir, scratch_dir, annotation_files)

    annotation_sanity_check(scratch_dir)

    for arm_entry in arms_of_interest():
        create_sample_genome(scratch_dir, arm_entry)
    exit(1)
    # iss generate --genomes sample_genome.fa --model miseq --output sample_reads -n 0.5M -p 8
    # I probably don't need 0.5M - maybe 10K (?) will have to experimant some with that

    # store location of the fastq files in the database
