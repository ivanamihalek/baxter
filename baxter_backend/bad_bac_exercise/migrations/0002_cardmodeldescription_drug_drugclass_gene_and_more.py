# Generated by Django 5.1.3 on 2024-11-21 09:34

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bad_bac_exercise', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='CARDModelDescription',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'db_table': 'card_model_description',
                'db_table_comment': 'CARD refers to Comprehensive Antibiotic Resistance Database',
            },
        ),
        migrations.CreateModel(
            name='Drug',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=20)),
            ],
            options={
                'db_table': 'drugs',
            },
        ),
        migrations.CreateModel(
            name='DrugClass',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=20)),
            ],
            options={
                'db_table': 'drug_classes',
            },
        ),
        migrations.CreateModel(
            name='Gene',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('gene_name', models.CharField(max_length=20)),
                ('protein_seq', models.TextField()),
                ('dna_seq', models.TextField()),
                ('taxonomy_id', models.IntegerField()),
            ],
            options={
                'db_table': 'genes',
            },
        ),
        migrations.CreateModel(
            name='TaxonomyName',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'db_table': 'taxonomy_names',
            },
        ),
        migrations.CreateModel(
            name='AntibioticResMutation',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('protein', models.CharField(max_length=20)),
                ('dna', models.CharField(blank=True, max_length=20)),
                ('drugs_affected', models.ManyToManyField(to='bad_bac_exercise.drug')),
                ('drug_classes_affected', models.ManyToManyField(to='bad_bac_exercise.drugclass')),
                ('gene', models.ManyToManyField(to='bad_bac_exercise.gene')),
            ],
            options={
                'db_table': 'antibiotic_res_mutations',
                'db_table_comment': 'Mutations conferring antibiotic resistance',
            },
        ),
        migrations.CreateModel(
            name='PDBStructure',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('pdb_id', models.CharField(max_length=20, unique=True)),
                ('drug', models.ManyToManyField(to='bad_bac_exercise.drug')),
            ],
            options={
                'db_table': 'pdb_structures',
            },
        ),
    ]