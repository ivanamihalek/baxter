#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
import json
from pprint import pprint

import numpy as np
from Bio.Blast import NCBIXML
from Bio.PDB import PDBParser, Selection

from .utils import is_nonempty_file
from bad_bac_exercise.models import AntibioticResMutation, Pdb2Gene


def run():
    # grep -i mupirocin 1qu3.pdb | grep HETNAM -- what if it fails?
    pass


#######################
if __name__ == "__main__":
    run()
