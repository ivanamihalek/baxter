# Generated by Django 5.1.3 on 2024-11-28 19:24

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bad_bac_exercise', '0010_alter_pdb2drug_drug_name_in_pdb'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pdb2mutation',
            name='dist_to_drug',
            field=models.FloatField(db_comment='Angstroms'),
        ),
    ]
