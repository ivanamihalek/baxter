from django.db import models

# Create your models here.
#
# class Variant(models.Model):
#     gdna_start             = models.IntegerField(blank=True, null=True)
#     gdna_end               = models.IntegerField(blank=True, null=True)
#     mod_type               = models.CharField(max_length=6, blank=True, null=True)
#     nt_ref                 = models.CharField(max_length=50, blank=False, null=False)
#     nt_alt                 = models.CharField(max_length=50, blank=False, null=False)
#     cdna                   = models.CharField(max_length=100, blank=False, null=False)
#     protein                = models.CharField(max_length=100, blank=True, null=True)
#     splicing               = models.CharField(max_length=100, blank=True, null=True)
#     gnomad_freqs_id        = models.ForeignKey(GnomadFreqs, on_delete=models.SET_NULL, null=True)
#     protein_domain         = models.CharField(max_length=6, blank=True, null=True)
#     sec_structure_element  = models.CharField(max_length=6, blank=True, null=True)
#     present_in_ortho_verts = models.IntegerField(blank=True, null=True)
#     present_in_para_verts  = models.IntegerField(blank=True, null=True)
#     systems_effect         = models.CharField(max_length=9, blank=True, null=True)
#
#     class Meta:
#         db_table = 'variants'
