# Generated by Django 5.1.3 on 2024-11-28 14:30

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('bad_bac_exercise', '0009_alter_pdb2drug_drug_name_in_pdb'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pdb2drug',
            name='drug_name_in_pdb',
            field=models.CharField(blank=True, max_length=500, null=True),
        ),
    ]
