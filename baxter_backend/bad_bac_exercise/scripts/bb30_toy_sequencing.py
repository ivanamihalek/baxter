#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such

import gffutils
import os.path
import random
import subprocess

import wget
from os import listdir

from Bio import Entrez

from bad_bac_exercise.models import AntibioticResMutation, Pdb2Gene, Pdb2Drug, Pdb2Mutation, PDBStructure, Decoy, \
    Gene2UCSCAssembly


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


def download_annotation_files(storage_dir) -> list:
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
            print(f"{fnm_unzipped} found in {storage_dir}")
        else:
            fnm = fnm_unzipped + ".gz"
            url = f"{base_url}/{refseq[:3]}/{refseq[4:7]}/{refseq[7:10]}/{refseq[10:13]}/{refseq}_{ncbi_acc}/{fnm}"
            print(url)
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


def ann_files_to_db(storage_dir, scratch_dir, annotation_files):
    for annf in annotation_files:
        infnm = f"{storage_dir}/{annf}"
        dbfnm = f"{scratch_dir}/{annf}.db"
        if os.path.exists(dbfnm):
            print(f"found {dbfnm} in {scratch_dir}")
        else:
            gffutils.create_db(infnm, dbfn=dbfnm, force=True, keep_order=True,
                               merge_strategy='merge', sort_attribute_values=True)


def annotation_sanity_check(scratch_dir):

    # for each gene of interest, check that the start/end coords match gff
    # https://daler.github.io/gffutils/
    for gene_entry in genes_of_interest():
        assembly = gene_entry.ucsc_assemblies.first()
        refseq = assembly.refseq_assembly_id
        ncbi_acc = assembly.ncbi_accession_number
        fnm = f"{refseq}_{ncbi_acc}_genomic.gff"
        gene_2_assmb_entry : Gene2UCSCAssembly = Gene2UCSCAssembly.objects.get(gene_id=gene_entry.id,
                                                                               assembly_id=assembly.id)
        print()
        print(gene_entry.name, assembly.refseq_assembly_id, fnm)
        print(gene_2_assmb_entry.contig, gene_2_assmb_entry.start_on_contig, gene_2_assmb_entry.end_on_contig)

        dbfnm = f"{scratch_dir}/{fnm}.db"
        db = gffutils.FeatureDB(dbfnm, keep_order=True)
        gene_annotation = db[gene_entry.name]


def run():

    storage_dir = "/storage/databases/ncbi/genome_annotation/bacterial"
    scratch_dir = "/home/ivana/scratch/baxter"

    annotation_files = download_annotation_files(storage_dir)
    ann_files_to_db(storage_dir, scratch_dir, annotation_files)

    annotation_sanity_check(scratch_dir)

    # extend the region around the gene to create "sample genome" to use in inSilicoSeq
    # make sure that the region contains some unnannotatd intervals place some decoy mutations there - single nucl, or small deletion
    # in the same gene, away from the mutation of interest put a silent mutation there

    # iss generate --genomes sample_genome.fa --model miseq --output sample_reads -n 0.5M -p 8
    # I probably don't need 0.5M - maybe 10K (?) will have to experimant some with that

    # store location of the fastq files in the database

    pass
