# Generated by Django 5.1.3 on 2024-11-21 11:59

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bad_bac_exercise', '0002_cardmodeldescription_drug_drugclass_gene_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='pdbstructure',
            name='drug_class',
            field=models.ManyToManyField(to='bad_bac_exercise.drugclass'),
        ),
        migrations.CreateModel(
            name='Gene2UCSCAssembly',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('pct_identity', models.IntegerField()),
                ('assembly', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='bad_bac_exercise.ucscassembly')),
                ('gene', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='bad_bac_exercise.gene')),
            ],
            options={
                'db_table': 'gene_2_ucsc_assembly',
            },
        ),
        migrations.AddField(
            model_name='gene',
            name='ucsc_assemblies',
            field=models.ManyToManyField(through='bad_bac_exercise.Gene2UCSCAssembly', to='bad_bac_exercise.ucscassembly'),
        ),
        migrations.CreateModel(
            name='Pdb2Mutation',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('dist_to_drug', models.IntegerField(db_comment='Angstroms')),
                ('antibio_res_mutation', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='bad_bac_exercise.antibioticresmutation')),
                ('pdb_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='bad_bac_exercise.pdbstructure')),
            ],
            options={
                'db_table': 'pdb_2_mutation',
            },
        ),
        migrations.AddField(
            model_name='pdbstructure',
            name='antibio_res_mutation',
            field=models.ManyToManyField(through='bad_bac_exercise.Pdb2Mutation', to='bad_bac_exercise.antibioticresmutation'),
        ),
    ]
