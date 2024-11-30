# Generated by Django 5.1.3 on 2024-11-30 20:11

import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bad_bac_exercise', '0021_alter_decoy_table_comment_alter_decoy_source'),
    ]

    operations = [
        migrations.CreateModel(
            name='Fingerprint',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('identifier', models.CharField(max_length=30)),
                ('dna_seq', models.TextField()),
                ('assembly', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='bad_bac_exercise.ucscassembly')),
            ],
            options={
                'db_table': 'fingerprints',
            },
        ),
    ]
