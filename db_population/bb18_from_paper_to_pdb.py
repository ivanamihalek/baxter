#! /usr/bin/env python
from __future__ import annotations
from dotenv import load_dotenv
load_dotenv()

import os
import django
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'settings')
django.setup()


from models.bad_bac_models import EnvironmentPublication, EnvPublication2Mutation


def main():

    results = EnvironmentPublication.objects.prefetch_related(
        'mutation__antibiotic_resistant_mutations__pdbs_available__pdb_structures'
    )

    for a in results:
        print(a)
        # for a2b in EnvironmentPublication.mutation.all():
        #     for b2c in a2b.b.b2c_set.all():
        #         print(a.name, a2b.b.title, b2c.c.description)


if __name__ == "__main__":
    main()
