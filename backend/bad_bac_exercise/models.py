from Bio import PDB
from django.db import models

"""
MAIN TABLES
    ucsc_assemblies
    antibio_res_mutations
    genes
    drugs
    drug_classes
    pdb_structures

AUX TABLES
    card_model_descriptions
    taxonomy_names

JUNCTION TABLES
# https://zerotobyte.com/django-many-to-many-relationship-explained/
# many to many with additional fields
# https://stackoverflow.com/questions/4443190/djangos-manytomany-relationship-with-additional-fields
# also
# https://stackoverflow.com/questions/37650362/understanding-manytomany-fields-in-django-with-a-through-model
    mutation2drug
    mutation2drug_class
    gene2assembly
    structure2gene
    structure2drug
    structure2drug_class
"""


class Publication(models.Model):
    pubmed_id = models.IntegerField(null=False, unique=True, default=None)

    class Meta:
        db_table = 'publications'


# 'model' here is the pathogenicity model; nothing to do with Django
class CARDModel(models.Model):
    card_name = models.CharField(max_length=20, blank=False, null=False, unique=True)
    card_description = models.TextField()
    publications = models.ManyToManyField(Publication, db_table="card_model_2_publication")

    class Meta:
        db_table = 'card_models'
        db_table_comment = 'CARD refers to Comprehensive Antibiotic Resistance Database'


class TaxonomyName(models.Model):
    tax_id = models.IntegerField(null=False, unique=True, default=None)
    tax_name = models.CharField(max_length=100, blank=False, null=False, default=None)

    class Meta:
        db_table = 'taxonomy_names'


class UCSCAssembly(models.Model):
    ncbi_accession_number = models.CharField(max_length=50, blank=False, null=False, unique=True)
    refseq_assembly_id = models.CharField(max_length=50, blank=False, null=False, unique=True)
    common_name = models.CharField(max_length=200, blank=False, null=False)
    scientific_name = models.CharField(max_length=100, blank=False, null=False)
    biosample = models.CharField(max_length=50, blank=False, null=False)
    bioproject = models.CharField(max_length=50, blank=False, null=False)
    assembly_date = models.DateField()

    class Meta:
        db_table = 'ucsc_assemblies'


class Drug(models.Model):
    name = models.CharField(max_length=30, blank=False, null=False)
    aro_id = models.CharField(max_length=7, blank=False, null=True)
    pubchem_id = models.IntegerField(blank=False, null=True)
    is_discrete_structure = models.BooleanField(blank=False, null=True)
    # InchiKey spec says its 27 characters
    inchi_key = models.CharField(max_length=27, blank=False, null=True)
    # there can be isomers, with different smiles
    canonical_smiles = models.TextField(blank=False, null=True)

    class Meta:
        db_table = 'drugs'


class DrugClass(models.Model):
    name = models.CharField(max_length=50, blank=False, null=False)

    class Meta:
        db_table = 'drug_classes'


class Decoy(models.Model):
    # decoys are some random sequences from innocuous species
    class Source(models.IntegerChoices):
        PHAGES = 1, 'phages'
        HUMAN_MICROBIOME = 2, 'human_microbiome'
        PLANTS = 3, 'plants'

        @classmethod
        def translate(self, instr):
            if str(instr).lower() == 'phages':
                return self.PHAGES
            elif str(instr).lower() == 'human_microbiome':
                return self.HUMAN_MICROBIOME
            elif str(instr).lower() == 'plants':
                return self.PLANTS

            raise Exception(f"Unrecognized source for dexoy seqs: {instr}")

    source = models.IntegerField(choices=Source.choices,
                                 blank=False, null=False,
                                 db_comment='The meaning of the source fields: 1=pages, 2=human_microbiome, 3=plants')
    identifier   = models.CharField(max_length=30, blank=False, null=True)
    species_name = models.CharField(max_length=200, blank=False, null=True)
    dna_seq = models.TextField(blank=False, null=False)

    class Meta:
        db_table = 'decoys'


class Fingerprint(models.Model):
    identifier = models.CharField(max_length=30, blank=False, null=False)
    dna_seq  = models.TextField(blank=False, null=False)
    assembly = models.ForeignKey(UCSCAssembly, blank=False, null=False, on_delete=models.CASCADE)

    class Meta:
        db_table = 'fingerprints'


class Gene(models.Model):
    name = models.CharField(max_length=20, blank=False, null=False, unique=True)
    protein_seq = models.TextField(blank=False, null=False)
    dna_seq = models.TextField(blank=False, null=False)
    ucsc_assemblies = models.ManyToManyField(UCSCAssembly, through='Gene2UCSCAssembly')

    class Meta:
        db_table = 'genes'


class AntibioticResMutation(models.Model):
    mutation = models.CharField(max_length=20, blank=False, null=False)
    gene = models.ForeignKey(Gene, on_delete=models.PROTECT)
    # mapping to card model is here only for sanity checking
    # and reconstructing how we got to here
    card_models = models.ManyToManyField(CARDModel, db_table="abrm_2_card_model")
    drugs_affected = models.ManyToManyField(Drug, db_table="abrm_2_drug")
    drug_classes_affected = models.ManyToManyField(DrugClass, db_table="abrm_2_drug_class")
    # note that becasue the class PDBStructure is defined below,
    # we use a string reference
    pdbs_available = models.ManyToManyField('PDBStructure',  through="Pdb2Mutation")
    # assembly that has the reference amino acid right
    assemblies = models.ManyToManyField(UCSCAssembly, db_table="abrm_2_assembly")
    flagged = models.BooleanField(blank=False, null=False, default=False)
    comment = models.CharField(max_length=100, blank=False, null=True)

    class Meta:
        db_table = 'antibiotic_res_mutations'
        db_table_comment = 'Mutations conferring antibiotic resistance'


class PDBStructure(models.Model):
    pdb_id = models.CharField(max_length=4, blank=False, null=False, unique=True)
    drugs = models.ManyToManyField(Drug, through="Pdb2Drug")
    drug_classes = models.ManyToManyField(DrugClass, db_table="pdb_2_drug_class")
    abr_mutations = models.ManyToManyField(AntibioticResMutation, through="Pdb2Mutation")
    genes = models.ManyToManyField(Gene, through="Pdb2Gene")

    class Meta:
        db_table = 'pdb_structures'


class Gene2UCSCAssembly(models.Model):
    class Strand(models.IntegerChoices):
        PLUS =  1, 'plus'
        MINUS = 2, 'minus'

        @classmethod
        def translate(self, instr):
            if str(instr).lower() in ['plus', '+', '1']:
                return self.PLUS
            elif str(instr).lower() in ['minus', '-', '-1']:
                return self.MINUS
            raise Exception(f"Unrecognized strand symbol: {instr}")

    gene = models.ForeignKey(Gene, on_delete=models.CASCADE)
    assembly     = models.ForeignKey(UCSCAssembly, on_delete=models.CASCADE)
    pct_identity = models.IntegerField(blank=False, null=True)
    contig           = models.CharField(max_length=50, blank=False, null=False, unique=True)
    start_on_contig  = models.IntegerField(blank=False, null=True)
    end_on_contig    = models.IntegerField(blank=False, null=True)
    strand_on_contig = models.IntegerField(choices=Strand.choices, blank=False, null=True)

    class Meta:
        db_table = 'gene_2_ucsc_assembly'
        # adding some constraints to mke sure that the start is always amller than end
        # I am already storing the infor about the strandedness, so I don't need to indicate it by making start > end
        constraints = [
            models.CheckConstraint(
                check=models.Q(start_on_contig__lte=models.F('end_on_contig')),
                name='start_less_than_end'
            ),
            models.CheckConstraint(
                check=models.Q(start_on_contig__gte=0) & models.Q(end_on_contig__gte=0),
                name='both_columns_non_negative'
            ),
        ]


class Pdb2Drug(models.Model):
    pdb = models.ForeignKey(PDBStructure, on_delete=models.CASCADE)
    drug = models.ForeignKey(Drug, on_delete=models.CASCADE)
    # we are looking for something like this
    # HETNAM     MRC MUPIROCIN
    # except that if we are not lucky, the pdb entry might have
    # a different name for the same drug --> use canonical smiles
    drug_name_in_pdb  = models.CharField(max_length=500, blank=True, null=True)
    drug_residue_name = models.CharField(max_length=3, blank=True, null=True)

    class Meta:
        db_table = 'pdb_2_drug'


class Pdb2Gene(models.Model):
    # django will turn this  "pdb_id"
    pdb = models.ForeignKey(PDBStructure, on_delete=models.CASCADE)
    # django will turn this name int "gene_id"
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE)
    pdb_chain        = models.CharField(max_length=4, blank=True, null=True)
    pct_identity     = models.IntegerField(db_comment="Alignment identity")
    gene_seq_aligned = models.TextField()
    gene_seq_start   = models.IntegerField()
    gene_seq_end     = models.IntegerField()
    pdb_seq_aligned  = models.TextField()
    pdb_seq_start    = models.IntegerField()
    pdb_seq_end      = models.IntegerField()

    class Meta:
        db_table = 'pdb_2_gene'


class Pdb2Mutation(models.Model):
    # django will turn this name int "pdb_id"
    pdb = models.ForeignKey(PDBStructure, on_delete=models.CASCADE)
    antibio_res_mutation = models.ForeignKey(AntibioticResMutation, on_delete=models.CASCADE)
    dist_to_drug = models.FloatField(db_comment="Angstroms", null=True, blank=False)

    class Meta:
        db_table = 'pdb_2_mutation'
