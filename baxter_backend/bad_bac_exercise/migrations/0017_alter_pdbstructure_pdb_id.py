# Generated by Django 5.1.3 on 2024-11-24 21:20

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bad_bac_exercise', '0016_drug_aro_id_drug_pubchem_id_and_more'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pdbstructure',
            name='pdb_id',
            field=models.CharField(max_length=4, unique=True),
        ),
    ]
