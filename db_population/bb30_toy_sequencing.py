#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
from glob import glob

import gffutils
import os.path
import portion as P
import random
import subprocess

import wget
from Bio.Data import CodonTable

from Bio.Seq import Seq

from db_population.utils import run_subprocess

from bad_bac_exercise.models import AntibioticResMutation, PDBStructure, Gene2UCSCAssembly
from db_population.thrid_party_tools import bwa_mem2_alignment, gatk_haplotyper_variant_caller


class VariantDescription:
    assembly: str | None = None
    contig: str | None = None
    genomic_coordinate: int | None = None
    toy_coordinate: int | None = None
    strand: str | None = None
    codon_number: int | None = None
    codon_coordinate: int | None = None
    nt_from: str | None = None
    nt_to: str | None = None
    codon_from: str | None = None
    codon_to: str | None = None
    aa_from: str | None = None
    aa_to: str | None = None

    def __init__(self, assembly=None, contig=None, genomic_coordinate=None):
        self.assembly = assembly
        self.contig   = contig
        self.genomic_coordinate = genomic_coordinate

    def __str__(self):
        return f"""address {self.assembly} {self.contig} {self.genomic_coordinate}"
            strand: {self.strand}
            nt: {self.nt_from} > {self.nt_to}
            codon: {self.codon_from} > {self.codon_to}
            aa:  {self.aa_from} > {self.aa_to} """

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
    for arm_entry in AntibioticResMutation.objects.filter(assemblies__isnull=False).filter(flagged=False).distinct():
        pdb_entries = PDBStructure.objects.filter(
            pdb2mutation__antibio_res_mutation=arm_entry,
            pdb2mutation__dist_to_drug__lt=9
        ).distinct()
        if len(pdb_entries) < 1: continue
        arm_entries.add(arm_entry)
    return arm_entries


def mutation_breakdown(mutation):
    aa_from = mutation[0]
    aa_to = mutation[-1]
    pos = int(mutation[1:-1])
    return aa_from, pos, aa_to


def print_string_100(string: str, file=None):
    while string:
        print(string[:100], file=file)
        string = string[100:]


def sequence_from_blastdb(contig_name, start, end) -> str:
    blastdbcmd = "/usr/third/blast-2.16.0/bin/blastdbcmd"
    db = "/storage/databases/ucsc/bacterial_genomes/ucsc_bac_genomes.fa"
    cmd = f"{blastdbcmd} -db {db} -entry {contig_name} -range {start}-{end} "
    ret = run_subprocess(cmd)
    return "".join(ret.split("\n")[1:])


def mutation_spec(from_codon, aa_to):
    # according to transl table 11, how do I mutate my codon do get aa_to?
    codon_table = CodonTable.unambiguous_dna_by_id[11]
    to_codons = [codon for codon, amino_acid in codon_table.forward_table.items() if amino_acid == aa_to]

    mutation_codon_pos  = None
    mutation_nucleotide = None
    to_codon = None
    for to_codon in to_codons:
        diff_indices = [i for i in range(3) if from_codon[i] != to_codon[i]]
        # for now. we'll go for single nucl substitution
        if len(diff_indices) == 0:
            # the codon {from_codon} already codes for {aa_to}
            # this might be a desired effect if we are looking for a silent mutation
            continue
        elif len(diff_indices) > 1:
            continue
        else:
            i = random.choice(diff_indices)
            mutation_nucleotide = to_codon[i]
            mutation_codon_pos = i  # note this is 0-offset
            break
            # print(from_codon, to_codon)

    if mutation_codon_pos is None or mutation_nucleotide is None:
        print(f"bleep - I did not find a single nucleotide change to go from {from_codon} to {aa_to}")
        print(f"to codons: {to_codons}")
        return None, None

    return mutation_codon_pos, mutation_nucleotide, to_codon


def insert_toy_mutation(genome_seq, gene_start, gene_end, strand, mutation, verbose=False) -> (str, VariantDescription):
    # print(seq_on_assembly)
    # according to transl table 11, is the aa for the codon correct?
    (aa_from, pos, aa_to) = mutation_breakdown(mutation)

    if verbose: print("**********************************")
    biopython_mutated_coding_dseq = Seq(genome_seq[gene_start - 1:gene_end])
    if strand == "minus":
        biopython_mutated_coding_dseq = biopython_mutated_coding_dseq.reverse_complement()
    biopython_translation = biopython_mutated_coding_dseq.translate(table=11)
    if verbose: print(aa_from, biopython_translation[pos - 1], aa_to, strand)
    if aa_from != biopython_translation[pos - 1]:
        errmsg  = f"Error the mutation is from {aa_from}. "
        errmsg += f"However, the sequence at {pos} is {biopython_translation[pos - 1]}."
        print(errmsg)
        exit()

    codon_start = (pos - 1) * 3  # 0-offset
    condon_on_the_orig_dseq = biopython_mutated_coding_dseq[codon_start:codon_start + 3]
    if verbose:
        print(condon_on_the_orig_dseq, condon_on_the_orig_dseq.translate(table=11))
        print_string_100(biopython_translation)

    (mutation_codon_pos, mutation_nucleotide, to_codon) = mutation_spec(condon_on_the_orig_dseq, aa_to)
    if mutation_codon_pos is None: exit(1)

    # both  codon_start and mutation_codon_pos are 0-offset, so we can slice the sequence string
    mutation_cdna_pos = codon_start + mutation_codon_pos
    if strand == 'plus':
        # this is offset 1 coordinate, liked by blast
        mutation_genomic_position = gene_start + mutation_cdna_pos
    else:
        # we need to translate back to the plus strand
        mutation_nucleotide = str(Seq(mutation_nucleotide).complement())
        mutation_genomic_position = gene_end - mutation_cdna_pos

    mutated_toy_genome = (genome_seq[:mutation_genomic_position - 1]
                          + mutation_nucleotide + genome_seq[mutation_genomic_position:])

    if verbose: print("------------------------------------")
    biopython_mutated_coding_dseq = Seq(mutated_toy_genome[gene_start - 1:gene_end])
    if strand == "minus":
        biopython_mutated_coding_dseq = biopython_mutated_coding_dseq.reverse_complement()
    biopython_translation = biopython_mutated_coding_dseq.translate(table=11)
    if verbose: print(aa_from, biopython_translation[pos - 1], aa_to, strand)
    if aa_to != biopython_translation[pos - 1]:
        errmsg  = f"Error the mutation should be to {aa_to}. "
        errmsg += f"However, the mutated sequence at {pos} is {biopython_translation[pos - 1]}."
        print(errmsg)
        exit()

    codon_start = (pos - 1) * 3
    condon_on_the_mutated_dseq = biopython_mutated_coding_dseq[codon_start:codon_start + 3]
    if verbose:
        print(condon_on_the_mutated_dseq, condon_on_the_mutated_dseq.translate(table=11))
        print_string_100(biopython_translation)
        print()

    var_descr = VariantDescription()
    var_descr.strand = strand
    var_descr.toy_coordinate = mutation_genomic_position
    var_descr.nt_from = genome_seq[mutation_genomic_position-1]
    var_descr.nt_to = mutation_nucleotide
    var_descr.codon_from = condon_on_the_orig_dseq
    var_descr.codon_to = condon_on_the_mutated_dseq
    var_descr.aa_from = aa_from
    var_descr.aa_to = aa_to
    return mutated_toy_genome, var_descr


def insert_irrelevant_decoy(toy_genome, db, toy_genome_start, toy_genome_end):
    # db refers to the actual genome
    # toy_genome_start and toy_genome_end are the coordinates of the toy genome
    # on the actual genome

    # one decoy: SNV in an un-annotated region
    complement = P.closed(toy_genome_start, toy_genome_end)
    for feature in db.region(start=toy_genome_start, end=toy_genome_end):
        # feature.attributes is a dict
        if 'genome' in feature.attributes: continue  # this is the top-top-level interval
        if 'Parent' in feature.attributes: continue  # we are looking for parent-less intervals
        complement = complement - P.closed(feature.start, feature.end)

    reasonable_complement = [atomic for atomic in complement if (atomic.upper - atomic.lower) >= 10]
    if len(reasonable_complement) == 0:
        print(f"Error = there seem to be fewer un-annotated regions than I anticipated.")
        exit()

    # blank - un-annotated, without features
    random_blank_interval = random.choice(reasonable_complement)
    random_blank_range_clipped = range(random_blank_interval.lower+3, random_blank_interval.upper-3)
    random_position = random.choice(random_blank_range_clipped)
    position_on_toy_genome = random_position - toy_genome_start - 1
    nt = toy_genome[position_on_toy_genome]
    replacement = random.choice(list({'A', 'C', 'T', 'G'}.difference({nt})))
    # print(random_blank_interval, random_position, position_on_toy_genome, nt, replacement)
    return toy_genome[:position_on_toy_genome] + replacement + toy_genome[position_on_toy_genome+1:]


def insert_random_silent_mutation(toy_genome, db, contig_name, toy_genome_start, toy_genome_end, verbose=False):
    # another decoy: silent mutation in one of the CDS
    coding_seqs = [feature for feature in db.region(start=toy_genome_start-1, end=toy_genome_end+1, completely_within=True)
                   if feature.attributes.get('gbkey', ["none"])[0]=="CDS"]
    # pick one at random
    random_coding_seq = coding_seqs[0]

    if verbose:
        print(random_coding_seq.start, random_coding_seq.end, random_coding_seq.strand)
        print(random_coding_seq.attributes)
    # translate to make sure it's coding
    gene_seq = sequence_from_blastdb(contig_name, random_coding_seq.start, random_coding_seq.end)
    if verbose: print("**********************************")
    biopython_dseq = Seq(gene_seq)
    if random_coding_seq.strand == "-":
        biopython_dseq = biopython_dseq.reverse_complement()
    biopython_translation = biopython_dseq.translate(table=11)

    # pick codon at random
    (mutation_codon_pos, mutation_nucleotide) = (None, None)
    # pick a substitution that leads to the same amino acid
    # sometimes a codon is unique (M in table 11)
    panic_ctr = 0
    random_codon_start = 1
    while mutation_codon_pos is None:
        if (panic_ctr := panic_ctr + 1) > 10:
            print(f"something went wrong with my mutation scheme")
            exit(1)
        random_codon_start = random.choice(range(len(biopython_translation)))
        random_codon = biopython_dseq[3*random_codon_start:3*random_codon_start+3]
        random_aa = biopython_translation[random_codon_start]
        if verbose:
            print(biopython_translation)
            print(random_codon_start, random_aa, random_codon)
        # look for the silent mutation; yes, tee same one; we are looking for a silent mutation
        (mutation_codon_pos, mutation_nucleotide) = mutation_spec(random_codon, random_aa)
        if verbose: print(mutation_codon_pos, mutation_nucleotide)

    # both  codon_start and mutation_codon_pos are 0-offset, so we can slice the sequence string

    mutation_cdna_pos = random_codon_start + mutation_codon_pos
    if  random_coding_seq.strand == '+':
        # this is offset 1 coordinate, liked by blast
        mutation_genomic_position = random_coding_seq.start + mutation_cdna_pos
    else:
        mutation_nucleotide = str(Seq(mutation_nucleotide).complement())
        mutation_genomic_position = random_coding_seq.end - mutation_cdna_pos

    mutation_toy_genomic_position = mutation_genomic_position - toy_genome_start + 1
    mutated_toy_genome = (toy_genome[:mutation_toy_genomic_position - 1]
                          + mutation_nucleotide + toy_genome[mutation_toy_genomic_position:])

    if verbose: print("------------------------------------")
    muated_biopython_dseq = Seq(mutated_toy_genome[random_coding_seq.start - 1:random_coding_seq.end])
    if random_coding_seq.strand == "-":
        muated_biopython_dseq = biopython_dseq.reverse_complement()
    muated_biopython_translation = biopython_dseq.translate(table=11)
    if verbose:
        print(muated_biopython_dseq == biopython_dseq)
        print(muated_biopython_translation == biopython_translation)

    # recalculate the coordinate from cds to genome, and from genome to toy genome
    return mutated_toy_genome


def create_sample_genome(scratch_dir, arm_entry):

    # extend the region around the gene to create "sample genome" to use in inSilicoSeq
    # make sure that the region contains some unnannotatd intervals place some decoy mutations there - single nucl, or small deletion
    # in the same gene, away from the mutation of interest put a silent mutation there

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
    # if strand == 'plus': return

    actual_gene_start  = gene_2_assmb_entry.start_on_contig
    actual_gene_end    = gene_2_assmb_entry.end_on_contig
    gene_length = actual_gene_end - actual_gene_start + 1
    extension = int(2.e4)
    toy_genome_start = actual_gene_start - extension
    toy_genome_end   = actual_gene_end + extension

    toy_genome = sequence_from_blastdb(gene_2_assmb_entry.contig, toy_genome_start, toy_genome_end)
    gene_start_on_toy_genome = extension + 1
    gene_end_on_toy_genome   = gene_start_on_toy_genome + gene_length - 1

    var_descr: VariantDescription
    (toy_genome, var_descr) = insert_toy_mutation(toy_genome, gene_start_on_toy_genome, gene_end_on_toy_genome,
                                     strand, arm_entry.mutation, verbose=False)
    var_descr.assembly = refseq
    var_descr.contig = gene_2_assmb_entry.contig
    var_descr.genomic_coordinate = toy_genome_start + var_descr.toy_coordinate - 1
    print(f"inserted toy variant to result in {arm_entry.mutation} in {arm_entry.gene.name}")
    print(var_descr)
    print()
    #
    # dbfnm = f"{scratch_dir}/{fnm}.db"
    # db = gffutils.FeatureDB(dbfnm, keep_order=True)
    # toy_genome = insert_irrelevant_decoy(toy_genome, db, toy_genome_start, toy_genome_end)
    #
    # toy_genome = insert_random_silent_mutation(toy_genome, db, gene_2_assmb_entry.contig, toy_genome_start, toy_genome_end)
    return toy_genome


def toy_sequencing(outdir):
    i_came_from = os.getcwd()
    os.chdir(outdir)
    # not that inSilicoSequencing, of which iss is part, must be installed in the local env
    sample_reads_name_root = "sample_reads"
    cmd = f"iss generate --genomes sample_genome.fa --model miseq --output {sample_reads_name_root} -n 10k -p 8"
    ret = run_subprocess(cmd)
    # ISS creates some junk
    os.remove(f"{sample_reads_name_root}_abundance.txt")
    for fnm in glob(f"{sample_reads_name_root}*tmp*.vcf"):
        os.remove(fnm)
    os.chdir(i_came_from)


def check_variants(assmb_entry, outdir):

    print("checking the synthetic variants")

    # align the sequences
    print(f"\treads alignment")
    reference_genome = f"/storage/databases/ucsc/bacterial_genomes/{assmb_entry.refseq_assembly_id}.fa"
    sample_reads = [f"sample_reads_R{i}.fastq" for i in range(1, 3)]
    bam_file = f"toy_alignment.bam"
    bwa_mem2_alignment(reference_genome, outdir, sample_reads,  bam_file)

    # call the variants
    print(f"\tvariants calling")
    out_vcf = "toy_vars_called.vcf"
    gatk_haplotyper_variant_caller(reference_genome, outdir, bam_file, out_vcf)

    # annotate # TODO I am here
    # after a lot of trouble I go the dockerized version of VeP to run, but it got
    # a single simple annotation wrong
    # snpEff - the docu so bad I could not figure out what the input should be

    # compare with the intended
    # vcfanno seems the most reasonable, but it needs the gff reformatted https://brentp.github.io/vcfanno/
    # it cannot annotate the effect on the protein level
    # use vcfanno to at least check if the toy variants are inside the gene or in an intergenic region


def run():

    storage_dir = "/storage/databases/ncbi/genome_annotation/bacterial"
    scratch_dir = "/home/ivana/scratch/baxter"

    annotation_files = download_annotation_files(storage_dir)
    ann_files_to_db(storage_dir, scratch_dir, annotation_files)

    annotation_sanity_check(scratch_dir)

    for arm_entry in arms_of_interest():
        gene_entry  = arm_entry.gene
        assmb_entry = gene_entry.ucsc_assemblies.first()
        # if arm_entry.mutation != "G406D": continue
        # if gene_entry.id != 40: continue
        # if assmb_entry.id != 66: continue
        print(gene_entry.id, assmb_entry.id, arm_entry.mutation)
        print(assmb_entry.refseq_assembly_id, assmb_entry.common_name)

        gene_2_assmb_entry: Gene2UCSCAssembly = Gene2UCSCAssembly.objects.get(gene_id=gene_entry.id, assembly_id=assmb_entry.id)
        # print(arm_entry.gene.name, arm_entry.mutation, assmb_entry.ncbi_accession_number, assmb_entry.refseq_assembly_id)
        # print(gene_2_assmb_entry.start_on_contig, gene_2_assmb_entry.end_on_contig)
        # print(assmb_entry)
        toy_genome = create_sample_genome(scratch_dir, arm_entry)

        outdir = f"{scratch_dir}/toy_genomes/{arm_entry.gene.name}_{arm_entry.mutation}"
        if not os.path.exists(outdir): os.makedirs(outdir)

        outfnm = f"sample_genome.fa"
        with open(f"{outdir}/{outfnm}", "w") as outf:
            print(f">sample", file=outf)
            print_string_100(toy_genome, file=outf)
            # store location of the fastq files in the database

        toy_sequencing(outdir)

        check_variants(assmb_entry, outdir)

        exit()
