#! /usr/bin/env python
# this is meant to be run with
# ./manage.py runscript bb03_parse_card_json
# in that case django will take care of the paths and also check for migrations and such
import os
import wget
import subprocess
from bad_bac_exercise.models import UCSCAssembly


def download_all():
    for assembly_entry in UCSCAssembly.objects.all():
        refseq = assembly_entry.refseq_assembly_id
        fnm_unzipped = f"{refseq}.fa"
        if os.path.exists(fnm_unzipped): continue

        print(assembly_entry.ncbi_accession_number, assembly_entry.refseq_assembly_id)

        fnm = f"{fnm_unzipped}.gz"
        if not os.path.exists(fnm):
            # https://hgdownload.soe.ucsc.edu/hubs/GCA/902/728/005/GCA_902728005.1/GCA_902728005.1.fa.gz
            url = f"https://hgdownload.soe.ucsc.edu/hubs/{refseq[:3]}/{refseq[4:7]}/{refseq[7:10]}/{refseq[10:13]}/{refseq}/{fnm}"
            print(url)
            wget.download(url, bar=wget.bar_thermometer)

        subprocess.run(["gunzip", fnm])
        print()


def create_contents_files():

    for assembly_entry in UCSCAssembly.objects.all():
        refseq = assembly_entry.refseq_assembly_id
        fnm_unzipped = f"{refseq}.fa"
        if os.path.exists(fnm_unzipped):
            subprocess.run(f"grep '>' {fnm_unzipped}  > {refseq}.contents.txt", shell=True)
        else:
            print(f"{fnm_unzipped} not found - no content created")

def run():
    homedir = "/storage/databases/ucsc/bacterial_genomes"
    os.chdir(homedir)
    all_fa = "uscs_bac_genomes.fa"
    if not os.path.exists(all_fa):
        download_all()
        subprocess.run(f"cat G*.fa > {all_fa}", shell=True)

    create_contents_files()
    # subprocess.run(f"makeblastdb -in {all_fa} -dbtype nucl ", shell=True)

    pass


#######################
if __name__ == "__main__":
    run()
