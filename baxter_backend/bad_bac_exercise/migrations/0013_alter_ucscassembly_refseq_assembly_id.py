# Generated by Django 5.1.3 on 2024-11-29 16:45

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bad_bac_exercise', '0012_alter_pdb2mutation_dist_to_drug'),
    ]

    operations = [
        migrations.AlterField(
            model_name='ucscassembly',
            name='refseq_assembly_id',
            field=models.CharField(max_length=50, unique=True),
        ),
    ]