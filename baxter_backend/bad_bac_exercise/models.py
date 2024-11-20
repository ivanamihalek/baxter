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
    mutation2drug
    mutation2drug_class
    mutation2gene
    gene2assembly
    structure2gene
    structure2drug
    structure2drug_class
"""

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


class AntibioResMutation(models.Model):

    class Meta:
        db_table = 'ucsc_assemblies'
        db_table_comment = 'Mutations conferring antibiotic resistance'


class Gene(models.Model):

    class Meta:
        db_table = 'genes'


class Drug(models.Model):

    class Meta:
        db_table = 'drugs'


class DrugClass(models.Model):

    class Meta:
        db_table = 'drug_classes'


class PDBStructure(models.Model):

    class Meta:
        db_table = 'pdb_structures'


#######################################################################
#  AUX TABLES
class CARDModelDescription(models.Model):

    class Meta:
        db_table = 'card_model_description'
        db_table_comment = 'CARD refers to Comprehensive Antibiotic Resistance Database'


class TaxonomyName(models.Model):

    class Meta:
        db_table = 'taxonomy_names'

#######################################################################
# JUNCTION TABLES
# https://zerotobyte.com/django-many-to-many-relationship-explained/
# many to many with additional fields
# https://stackoverflow.com/questions/4443190/djangos-manytomany-relationship-with-additional-fields
#     mutation2drug
#     mutation2drug_class
#     mutation2gene
#     gene2assembly
#     structure2gene
#     structure2drug
#     structure2drug_class
