# Generated by Django 5.1.3 on 2024-11-24 15:12

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('bad_bac_exercise', '0013_drug_inchi_key_drug_smiles'),
    ]

    operations = [
        migrations.RenameField(
            model_name='drug',
            old_name='smiles',
            new_name='canonical_smiles',
        ),
    ]
