from django.db import models


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
