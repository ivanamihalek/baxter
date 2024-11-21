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


#######################################################################
#  AUX TABLES

# 'model' here is the pathogenicity model; nothing to do with Django
class CARDModel(models.Model):
    card_name = models.CharField(max_length=20, blank=False, null=False, unique=True)
    card_description = models.TextField()

    class Meta:
        db_table = 'card_models'
        db_table_comment = 'CARD refers to Comprehensive Antibiotic Resistance Database'


class TaxonomyName(models.Model):
    tax_id   = models.IntegerField(null=False, unique=True, default=None)
    tax_name = models.CharField(max_length=100, blank=False, null=False, default=None)

    class Meta:
        db_table = 'taxonomy_names'


#######################################################################
#  MAIN TABLES
class UCSCAssembly(models.Model):
    ncbi_accession_number = models.CharField(max_length=50, blank=False, null=False, unique=True)
    refseq_assembly_id   = models.CharField(max_length=50, blank=False, null=False)
    common_name          = models.CharField(max_length=200, blank=False, null=False)
    scientific_name      = models.CharField(max_length=100, blank=False, null=False)
    biosample            = models.CharField(max_length=50, blank=False, null=False)
    bioproject           = models.CharField(max_length=50, blank=False, null=False)
    assembly_date        = models.DateField()

    class Meta:
        db_table = 'ucsc_assemblies'


class Gene(models.Model):
    name        = models.CharField(max_length=20, blank=False, null=False, unique=True)
    protein_seq = models.TextField(blank=False, null=False)
    dna_seq     = models.TextField(blank=False, null=False)
    ucsc_assemblies = models.ManyToManyField(UCSCAssembly, through='Gene2UCSCAssembly')

    class Meta:
        db_table = 'genes'


class Gene2UCSCAssembly(models.Model):
    gene = models.ForeignKey(Gene, on_delete=models.CASCADE)
    assembly = models.ForeignKey(UCSCAssembly, on_delete=models.CASCADE)
    pct_identity = models.IntegerField()

    class Meta:
        db_table = 'gene_2_ucsc_assembly'


class Drug(models.Model):
    name  = models.CharField(max_length=30, blank=False, null=False)

    class Meta:
        db_table = 'drugs'


class DrugClass(models.Model):
    name = models.CharField(max_length=50, blank=False, null=False)

    class Meta:
        db_table = 'drug_classes'


class AntibioticResMutation(models.Model):
    mutation       = models.CharField(max_length=20, blank=False, null=False)
    gene           = models.ForeignKey(Gene, on_delete=models.PROTECT)
    # mapping to card model is here only for sanity checking
    # and reconstructing how we got to here
    card_models    = models.ManyToManyField(CARDModel, db_table="abrm_2_card_model")
    drugs_affected = models.ManyToManyField(Drug, db_table="abrm_2_drug")
    drug_classes_affected = models.ManyToManyField(DrugClass, db_table="abrm_2_drug_class")

    class Meta:
        db_table = 'antibiotic_res_mutations'
        db_table_comment = 'Mutations conferring antibiotic resistance'


class PDBStructure(models.Model):
    pdb_id        = models.CharField(max_length=20, blank=False, null=False, unique=True)
    drugs         = models.ManyToManyField(Drug, db_table="pdb_2_drug")
    drug_classes  = models.ManyToManyField(DrugClass, db_table="pdb_2_drug_class")
    abr_mutations = models.ManyToManyField(AntibioticResMutation, through="Pdb2Mutation")

    class Meta:
        db_table = 'pdb_structures'


class Pdb2Mutation(models.Model):
    pdb_id               = models.ForeignKey(PDBStructure, on_delete=models.CASCADE)
    antibio_res_mutation = models.ForeignKey(AntibioticResMutation, on_delete=models.CASCADE)
    dist_to_drug         = models.IntegerField(db_comment="Angstroms")

    class Meta:
        db_table = 'pdb_2_mutation'

