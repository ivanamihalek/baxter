# Generated by Django 5.1.3 on 2024-11-30 18:26

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('bad_bac_exercise', '0019_alter_decoy_species_name'),
    ]

    operations = [
        migrations.AlterModelTableComment(
            name='decoy',
            table_comment='The meaning of source fields: 1=pages, 2=human_microbiome, 3=plants',
        ),
    ]